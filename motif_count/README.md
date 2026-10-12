# motif_count — signed-XOR motif finders

C tools that enumerate the 8-node, 12-edge signed-XOR motif in a connectome
edge list. One finder per dataset format; the motif definition is identical
across them.

| file                                                                     | input format                                 |
| ------------------------------------------------------------------------ | -------------------------------------------- |
| [`src/signed_xor_motif_count.c`](src/signed_xor_motif_count.c)           | Allen V1 / MICrONS proofread CSV (5 columns) |
| [`src/drosophila_signed_xor_motif.c`](src/drosophila_signed_xor_motif.c) | FlyWire / FAFB edge list                     |
| [`src/celegans_signed_xor_motif.c`](src/celegans_signed_xor_motif.c)     | C. elegans (Cook et al.) connectome          |

## Folder layout

```text
motif_count/
├── Makefile
├── README.md
├── src/                # The three dataset counters
│   └── gpu/            # Shared GPU search, graph packing, Metal and CUDA backends
├── include/            # Shared GPU API and OpenMP compatibility headers
├── third_party/        # Bundled uthash.h and its license notice
├── scripts/            # Allen V1 and C. elegans input preparation
├── tests/              # Regression suite and shared-kernel CPU interpreter
└── build/              # Generated executables and intermediate files (Git-ignored)
    ├── cpu/
    ├── metal/
    ├── cuda/
    ├── kernel-reference/
    └── legacy/         # Preserved executables from before the reorganization
```

Run the build and usage commands below from `motif_count/`. Build directories
are created as needed. `make clean` removes generated backend builds and
retains `build/legacy/`. Data preparation commands using `scripts/` are in
[`../data/README.md`](../data/README.md).

## Build

The CPU builds use a C compiler with OpenMP and the bundled `third_party/uthash.h`.
Optional GPU builds use Apple's Metal framework on Apple Silicon Macs or
NVIDIA's CUDA Toolkit on Linux (x86-64 and ARM64, including DGX Spark).
Windows users can use the Linux builds through WSL 2, as described below.

### Linux

```bash
make cpu CC=gcc
```

### macOS

Apple's `/usr/bin/gcc` is Apple Clang and does not support `-fopenmp` by
default. Install Homebrew LLVM, which includes the OpenMP runtime:

```bash
brew install llvm
```

Then build all three counters with the Homebrew compiler. Using `brew --prefix`
works with both Apple Silicon and Intel Homebrew installations:

```bash
make cpu CC="$(brew --prefix llvm)/bin/clang"
```

CPU executables are written to `build/cpu/`. To add compiler options, set
`CFLAGS`, for example `CFLAGS="-O3 -march=native -Wall -Wextra"`.

### Windows (Ubuntu on WSL 2)

Use Ubuntu on Windows Subsystem for Linux 2 (WSL 2) to build and run these
programs on Windows. The current Makefile uses a Unix shell and GCC; the
commands below run inside Ubuntu and produce Linux executables for WSL.

#### Install WSL and build the CPU programs

On Windows 11 or Windows 10 version 2004/build 19041 or later, open
**PowerShell as Administrator** and install Ubuntu:

```powershell
wsl --install -d Ubuntu
```

Restart Windows if prompted, open Ubuntu, and create its Linux username and
password. From PowerShell, update WSL, check that Ubuntu's `VERSION` is `2`,
and enter its shell:

```powershell
wsl --update
wsl --list --verbose
wsl -d Ubuntu
```

If Ubuntu is listed as version `1`, run `wsl --set-version Ubuntu 2` in
PowerShell before entering Ubuntu. See Microsoft's
[WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install).

Run the remaining commands in the **Ubuntu terminal**:

```bash
sudo apt update
sudo apt install -y build-essential python3

# Example: the repository is on your Windows desktop.
# Replace YOUR_USERNAME and the folder path with your actual location.
cd "/mnt/c/Users/YOUR_USERNAME/Desktop/signed-xor-motifs/motif_count"

# Rebuild for WSL if this checkout contains binaries copied from another OS.
make clean
make cpu CC=gcc
make test CC=gcc
```

Windows drives are mounted under `/mnt/`: for example,
`C:\Users\Alice\Desktop\signed-xor-motifs` becomes
`/mnt/c/Users/Alice/Desktop/signed-xor-motifs`. Quote paths containing spaces.
For better file I/O performance, keep the project and large datasets in the
WSL filesystem, such as `~/signed-xor-motifs`. See Microsoft's
[filesystem guide](https://learn.microsoft.com/en-us/windows/wsl/filesystems).

From `motif_count/` in Ubuntu, run whichever counter matches your input:

```bash
export OMP_NUM_THREADS="$(nproc)"

./build/cpu/signed_xor_motif_count --min-weight 0.1 \
  ../data/allen_v1/allen_v1_edges.csv

./build/cpu/drosophila_signed_xor_motif \
  ../data/drosophila/flywire_fafb_v783_connections.csv

./build/cpu/celegans_signed_xor_motif \
  ../data/celegans/cook2019_hermaphrodite_chemical_nt_approx.csv
```

Supply your own datasets at these paths, or replace the input paths with
their actual locations. These CPU builds use OpenMP and require no GPU.

#### Optional NVIDIA CUDA on Windows through WSL 2

For a compatible NVIDIA GPU, install an up-to-date NVIDIA **Windows** driver.
Then install the CUDA Toolkit **inside Ubuntu**, following NVIDIA's
[CUDA on WSL guide](https://docs.nvidia.com/cuda/wsl-user-guide/index.html) and
the [CUDA download instructions](https://developer.nvidia.com/cuda-downloads)
for **Linux → x86_64 → WSL-Ubuntu**. Install the toolkit-only package; NVIDIA
instructs WSL users to keep the Windows driver and avoid installing a Linux
NVIDIA driver or the `cuda`/`cuda-drivers` packages inside WSL.

Ensure the toolkit's `bin` directory is on Ubuntu's `PATH` as described in
NVIDIA's [environment setup instructions](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#environment-setup).
From `motif_count/` in the **Ubuntu terminal**, verify the GPU and compiler,
then build and run:

```bash
nvidia-smi
nvcc --version
make cuda
make test-cuda

./build/cuda/drosophila_signed_xor_motif --gpu \
  ../data/drosophila/flywire_fafb_v783_connections.csv
```

If `nvidia-smi` is missing from `PATH`, try `/usr/lib/wsl/lib/nvidia-smi`.
The Allen and C. elegans CUDA commands are in the NVIDIA GPU section below.
Metal builds require macOS. Project execution on Windows/WSL 2 has not yet
been validated; run the regression commands above on your Windows machine.

## Run

Each tool takes an edge-list CSV and prints per-type motif counts. See the
header comment in each `src/*.c` file for its command-line options (e.g.
`-n <max_per_type>`, `--weight-stats <prefix>` in `signed_xor_motif_count.c`).

```bash
OMP_NUM_THREADS=$(nproc) ./build/cpu/signed_xor_motif_count <edges.csv>
```

On macOS, use `sysctl` instead of `nproc`. From the `motif_count` directory,
the downloaded datasets can be processed with:

```bash
OMP_NUM_THREADS="$(sysctl -n hw.logicalcpu)" \
  ./build/cpu/signed_xor_motif_count --min-weight 0.1 \
  ../data/allen_v1/allen_v1_edges.csv

OMP_NUM_THREADS="$(sysctl -n hw.logicalcpu)" \
  ./build/cpu/drosophila_signed_xor_motif \
  ../data/drosophila/flywire_fafb_v783_connections.csv

OMP_NUM_THREADS="$(sysctl -n hw.logicalcpu)" \
  ./build/cpu/celegans_signed_xor_motif \
  ../data/celegans/cook2019_hermaphrodite_chemical_nt_approx.csv
```

For a quick Allen V1 smoke test before processing all 70 million rows, add
`-t 1000000`:

```bash
OMP_NUM_THREADS="$(sysctl -n hw.logicalcpu)" \
  ./build/cpu/signed_xor_motif_count -t 1000000 --min-weight 0.1 \
  ../data/allen_v1/allen_v1_edges.csv
```

## GPU search on Apple Silicon (Metal)

All three programs support `--gpu` on Apple Silicon, including the M4 Pro.
Build from this directory using Apple's Command Line Tools and Python 3:

```bash
make metal
```

The executables are written to `build/metal/`. The Metal shader is embedded in the
executable and compiled by the Metal runtime; the separate Xcode Metal
Toolchain download is unnecessary. No Homebrew OpenMP installation is needed
for this build.

```bash
./build/metal/signed_xor_motif_count --gpu --min-weight 0.1 \
  ../data/allen_v1/allen_v1_edges.csv

./build/metal/drosophila_signed_xor_motif --gpu \
  ../data/drosophila/flywire_fafb_v783_connections.csv

./build/metal/celegans_signed_xor_motif --gpu \
  ../data/celegans/cook2019_hermaphrodite_chemical_nt_approx.csv
```

## GPU search on NVIDIA GPUs (CUDA)

On a Linux machine with an NVIDIA GPU, a compatible NVIDIA driver, and the
CUDA Toolkit (`nvcc` on `PATH`), build all three CUDA executables:

```bash
make cuda

./build/cuda/signed_xor_motif_count --gpu --min-weight 0.1 \
  ../data/allen_v1/allen_v1_edges.csv

./build/cuda/drosophila_signed_xor_motif --gpu \
  ../data/drosophila/flywire_fafb_v783_connections.csv

./build/cuda/celegans_signed_xor_motif --gpu \
  ../data/celegans/cook2019_hermaphrodite_chemical_nt_approx.csv
```

The CUDA build uses GCC for the C readers and OpenMP CPU path, and `nvcc`
with G++ for the CUDA backend and final linking. Configure these with
`CUDA_CC`, `CUDA_CXX`, and `NVCC` when necessary; the selected host compiler
must be supported by your CUDA Toolkit. Set `CUDA_OPENMP=` for a serial
CPU reference if OpenMP is unavailable.

The default `CUDA_ARCH=sm_75` includes Turing machine code and PTX that the
driver can compile for newer GPUs. Select a native architecture for the
machine where you build, or set one explicitly:

```bash
make cuda CUDA_ARCH=native

# DGX Spark / GB10, using a CUDA Toolkit that supports sm_121
make cuda CUDA_ARCH=sm_121
```

These are alternative build commands. When changing architecture or compiler
flags after a build, remove `build/cuda/` and rebuild. NVIDIA lists GB10 / DGX
Spark as compute capability 12.1 in its [GPU table](https://developer.nvidia.com/cuda/gpus);
architecture flags and PTX behavior are described in the
[nvcc documentation](https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-architecture-arch).

CUDA uses the first visible GPU by default. To select a specific GPU, use
`CUDA_VISIBLE_DEVICES`, for example:

```bash
CUDA_VISIBLE_DEVICES=1 ./build/cuda/drosophila_signed_xor_motif --gpu <edges.csv>
```

CUDA binaries live in `build/cuda/` and Metal binaries in `build/metal/`.
CUDA cannot run on the Apple GPU in this Mac; use the CUDA build on your
NVIDIA machine.

## GPU behavior and validation

CPU search remains the default. Pass `--gpu` to select the backend compiled
into the executable or `--cpu` to select the CPU explicitly. A GPU request
fails with an error if the binary lacks GPU support, the GPU/driver is
unavailable, allocation fails, or a GPU command fails. On macOS, sandboxed
processes may need GPU access from their host application.

Existing options remain available: `--induced` and `--out` for the FlyWire
and C. elegans counters, and `-t`, `-n`, `--pivot`, `--min-weight`,
`--weight-stats`, and `--sample-size` for the Allen counter. The GPU preserves
the current CPU algorithms, including A/B/AB variants, sign and identity
checks, duplicate-edge multiplicity, and Allen's original inhibitory target
pair order. CSV row order and reservoir samples can differ between backends;
counts and motif mappings should agree.

The CPU parses the CSV, prepares the graph, and generates bounded batches of
search jobs. The GPU searches for XOR and feedback nodes and validates induced
subgraphs using sorted adjacency lists. Counts use 64-bit integers. Counting
without `--induced` or exports also counts eligible feedback pairs
by inclusion/exclusion for FlyWire and C. elegans, avoiding a quadratic pair
loop while counting overlapping A/B matches once. Motif exports and Allen
weight statistics stream through bounded output pages;
dense output batches may be searched again to produce subsequent pages.
CSV parsing, host graph preparation, and writing results still take CPU time.
The graph must fit in host/GPU memory, each Metal buffer must fit the device's
buffer limit, and the packed adjacency storage has a 32-bit entry limit.
The CUDA backend uploads the packed graph once and allocates export staging
buffers only for runs that request motif CSV output or weight statistics.

GPU acceleration depends on graph size and degree distribution. Small graphs
such as the supplied C. elegans dataset can run faster on the CPU because
Metal startup and graph preparation cost more than the search itself. Start
with a count-only run when measuring performance; exporting every motif or
recording its weights adds processing and I/O.

The standard Metal build contains a serial CPU reference. For a parallel CPU
comparison, use the OpenMP build described above or:

```bash
# Linux
make cpu CC=gcc

# macOS with Homebrew LLVM
make cpu CC="$(brew --prefix llvm)/bin/clang"
```

These executables are written to `build/cpu/`.

### Validation

```bash
# CPU regression tests (use CC as above on macOS)
make test CC="$(brew --prefix llvm)/bin/clang"

# CPU/GPU comparison on an Apple Silicon Mac with Metal access
make test-gpu

# CPU/CUDA comparison on a Linux NVIDIA machine
make test-cuda

# Shared-kernel logic on a CPU; no GPU or CUDA Toolkit required
make test-kernel
```

## Data

Downloaded source files, provenance, validation counts, and preparation
commands are documented in [`../data/README.md`](../data/README.md). The large
data artifacts are intentionally ignored by Git and must be obtained from the
original providers under their respective terms.

## Third-party

`third_party/uthash.h` © Troy D. Hanson, distributed under the revised BSD license
(notice retained in the file header). Not covered by this repository's
CC BY-NC 4.0 license.
