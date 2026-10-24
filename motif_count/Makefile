# CPU, Metal, and CUDA binaries live separately. Apple Clang needs no OpenMP
# runtime for Metal builds. Run make from the motif_count directory.
PROGRAMS = signed_xor_motif_count drosophila_signed_xor_motif celegans_signed_xor_motif
METAL_C_OBJECTS = $(addprefix build/metal/,$(addsuffix .o,$(PROGRAMS)))
CUDA_C_OBJECTS = $(addprefix build/cuda/,$(addsuffix .o,$(PROGRAMS)))
REFERENCE_C_OBJECTS = $(addprefix build/kernel-reference/,$(addsuffix .o,$(PROGRAMS)))
INCLUDES = -Iinclude -Ithird_party -Isrc/gpu
COUNTER_HEADERS = third_party/uthash.h include/motif_gpu.h include/motif_omp.h src/gpu/motif_gpu_variants.h
CC = cc
CXX = c++
METAL_CC = /usr/bin/clang
CUDA_CC ?= gcc
CUDA_CXX ?= g++
NVCC ?= nvcc
# Build native machine code for the visible GPU by default. This avoids relying
# on driver-side PTX JIT compatibility when nvcc is newer than the GPU driver.
# Override with (for example) CUDA_ARCH=sm_86 when cross-compiling.
CUDA_ARCH ?= native
NVCCFLAGS ?= -O3 -std=c++14
CFLAGS ?= -O3 -Wall -Wextra
OPENMP ?= -fopenmp
CUDA_OPENMP ?= -fopenmp
PYTHON ?= python3

.PHONY: all cpu metal cuda cuda-toolchain test test-gpu test-cuda test-kernel clean
.SECONDARY:
all: cpu
cpu: $(addprefix build/cpu/,$(PROGRAMS))
metal: $(addprefix build/metal/,$(PROGRAMS))
cuda: $(addprefix build/cuda/,$(PROGRAMS))

cuda-toolchain:
	@command -v "$(NVCC)" >/dev/null 2>&1 || { echo "CUDA compiler not found: $(NVCC). Install the NVIDIA CUDA Toolkit on Linux; use 'make metal' on Apple Silicon." >&2; exit 1; }

build/cpu/%: src/%.c $(COUNTER_HEADERS)
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) $(OPENMP) $< -o $@ -lm

build/metal/%: build/metal/%.o build/metal/motif_gpu.o build/metal/motif_gpu_metal.o
	$(METAL_CC) $^ -o $@ -framework Metal -framework Foundation -lm

$(METAL_C_OBJECTS): build/metal/%.o: src/%.c $(COUNTER_HEADERS)
	@mkdir -p $(@D)
	$(METAL_CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -Wno-unknown-pragmas -DMOTIF_USE_METAL -c $< -o $@

build/metal/motif_gpu.o: src/gpu/motif_gpu.c include/motif_gpu.h
	@mkdir -p $(@D)
	$(METAL_CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -DMOTIF_USE_METAL -c $< -o $@

build/metal/motif_gpu_metal.o: src/gpu/motif_gpu_metal.m include/motif_gpu.h build/metal/motif_gpu_source.h
	$(METAL_CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -fobjc-arc -DMOTIF_USE_METAL -Ibuild/metal -c $< -o $@

build/metal/motif_gpu_source.h: src/gpu/motif_gpu_kernel.metal src/gpu/motif_gpu_search.h
	@mkdir -p $(@D)
	$(PYTHON) -c 'import json,pathlib; p=pathlib.Path("$<"); s=p.read_text().replace("#include \"motif_gpu_search.h\"", p.with_name("motif_gpu_search.h").read_text()); pathlib.Path("$@").write_text("static const char motif_gpu_source[] =\n" + "\n".join(json.dumps(line) for line in s.splitlines(keepends=True)) + ";\n")'

build/cuda/%: build/cuda/%.o build/cuda/motif_gpu.o build/cuda/motif_gpu_cuda.o
	$(NVCC) -ccbin "$(CUDA_CXX)" $^ -o $@ $(if $(strip $(CUDA_OPENMP)),-Xcompiler "$(CUDA_OPENMP)") -lm

$(CUDA_C_OBJECTS): build/cuda/%.o: src/%.c $(COUNTER_HEADERS) | cuda-toolchain
	@mkdir -p $(@D)
	$(CUDA_CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) $(CUDA_OPENMP) -Wno-unknown-pragmas -DMOTIF_USE_CUDA -c $< -o $@

build/cuda/motif_gpu.o: src/gpu/motif_gpu.c include/motif_gpu.h | cuda-toolchain
	@mkdir -p $(@D)
	$(CUDA_CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) $(CUDA_OPENMP) -DMOTIF_USE_CUDA -c $< -o $@

build/cuda/motif_gpu_cuda.o: src/gpu/motif_gpu_cuda.cu include/motif_gpu.h src/gpu/motif_gpu_search.h Makefile | cuda-toolchain
	@mkdir -p $(@D)
	$(NVCC) $(CPPFLAGS) $(INCLUDES) $(NVCCFLAGS) -arch=$(CUDA_ARCH) -ccbin "$(CUDA_CXX)" -DMOTIF_USE_CUDA -c $< -o $@

test: cpu
	$(PYTHON) tests/test_gpu.py --bin-dir build/cpu

test-gpu: metal
	$(PYTHON) tests/test_gpu.py --bin-dir build/metal --gpu

test-cuda: cuda
	$(PYTHON) tests/test_gpu.py --bin-dir build/cuda --gpu

# Hardware-independent checks of the exact shared Metal/CUDA kernel logic.
$(REFERENCE_C_OBJECTS): build/kernel-reference/%.o: src/%.c $(COUNTER_HEADERS)
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -Wno-unknown-pragmas -DMOTIF_USE_CUDA -c $< -o $@

build/kernel-reference/motif_gpu.o: src/gpu/motif_gpu.c include/motif_gpu.h
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -DMOTIF_USE_CUDA -c $< -o $@

build/kernel-reference/motif_gpu_reference.o: tests/motif_gpu_reference.cpp include/motif_gpu.h src/gpu/motif_gpu_search.h
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -std=c++14 -DMOTIF_USE_CUDA -c $< -o $@

build/kernel-reference/%: build/kernel-reference/%.o build/kernel-reference/motif_gpu.o build/kernel-reference/motif_gpu_reference.o
	$(CXX) $^ -o $@ -lm

test-kernel: $(addprefix build/kernel-reference/,$(PROGRAMS))
	$(PYTHON) tests/test_gpu.py --bin-dir build/kernel-reference --gpu

clean:
	rm -rf build/cpu build/metal build/cuda build/kernel-reference
