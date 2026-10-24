#include <cuda_runtime.h>
#include <stdlib.h>
#include "motif_gpu.h"

using Node = motif_gpu_node;
using Task = motif_gpu_task;
using Match = motif_gpu_match;
using Params = motif_gpu_params;
using gpu_uint = uint32_t;
using gpu_count = uint64_t;
enum { EXC = GPU_EXC, INH = GPU_INH, IN_EXC = GPU_IN_EXC, OTHER = GPU_OTHER };
#define GPU_INLINE __device__ __forceinline__
#define GPU_DEVICE
#define GPU_THREAD
#define GPU_PARAM
#include "motif_gpu_search.h"

static_assert(sizeof(Node) == 44, "CUDA Node layout");
static_assert(sizeof(Task) == 28, "CUDA Task layout");
static_assert(sizeof(Match) == 36, "CUDA Match layout");
static_assert(sizeof(Params) == 32, "CUDA Params layout");

__global__ void motif_search(const Node *nodes, const int32_t *edges,
                            const Task *tasks, uint64_t *counts,
                            const uint64_t *offsets, Match *matches, Params p) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    motif_search_task(nodes, edges, tasks, counts, offsets, matches, p, tid);
}

struct MotifCuda {
    Node *nodes;
    int32_t *edges;
    Task *tasks;
    uint64_t *counts, *offsets, *host_offsets;
    Match *matches, *host_matches;
};

static int check_cuda(cudaError_t status, const char *operation, char *error, size_t size) {
    if (status == cudaSuccess) return 1;
    snprintf(error, size, "CUDA %s: %s", operation, cudaGetErrorString(status));
    return 0;
}

void motif_backend_destroy(void *backend) {
    MotifCuda *m = static_cast<MotifCuda *>(backend);
    if (!m) return;
    if (m->nodes) cudaFree(m->nodes);
    if (m->edges) cudaFree(m->edges);
    if (m->tasks) cudaFree(m->tasks);
    if (m->counts) cudaFree(m->counts);
    if (m->offsets) cudaFree(m->offsets);
    if (m->matches) cudaFree(m->matches);
    free(m->host_offsets);
    free(m->host_matches);
    free(m);
}

static int upload(void **target, const void *source, size_t size,
                  char *error, size_t error_size) {
    if (!check_cuda(cudaMalloc(target, size ? size : 4), "graph allocation", error, error_size)) return 0;
    return !size || check_cuda(cudaMemcpy(*target, source, size, cudaMemcpyHostToDevice),
                               "graph upload", error, error_size);
}

void *motif_backend_create(const motif_gpu_graph *g, int quiet, char *error, size_t error_size) {
    int count = 0, device = 0;
    if (!check_cuda(cudaGetDeviceCount(&count), "device discovery", error, error_size)) return NULL;
    if (!count) {
        snprintf(error, error_size, "No NVIDIA CUDA GPU is visible (check the driver and CUDA_VISIBLE_DEVICES)");
        return NULL;
    }
    if (!check_cuda(cudaGetDevice(&device), "device selection", error, error_size)) return NULL;
    cudaDeviceProp properties;
    if (!check_cuda(cudaGetDeviceProperties(&properties, device), "device properties", error, error_size)) return NULL;
    MotifCuda *m = static_cast<MotifCuda *>(calloc(1, sizeof(*m)));
    if (!m) {
        snprintf(error, error_size, "Cannot allocate CUDA backend state");
        return NULL;
    }
    if (!upload(reinterpret_cast<void **>(&m->nodes), g->nodes, g->node_count * sizeof(*g->nodes), error, error_size) ||
        !upload(reinterpret_cast<void **>(&m->edges), g->edges, g->edge_count * sizeof(*g->edges), error, error_size) ||
        !check_cuda(cudaMalloc(reinterpret_cast<void **>(&m->tasks), MOTIF_GPU_TASK_BATCH * sizeof(Task)), "task allocation", error, error_size) ||
        !check_cuda(cudaMalloc(reinterpret_cast<void **>(&m->counts), MOTIF_GPU_TASK_BATCH * sizeof(uint64_t)), "count allocation", error, error_size) ||
        !check_cuda(cudaMalloc(reinterpret_cast<void **>(&m->offsets), MOTIF_GPU_TASK_BATCH * sizeof(uint64_t)), "offset allocation", error, error_size)) {
        motif_backend_destroy(m);
        return NULL;
    }
    if (!quiet) fprintf(stderr, "GPU: %s (CUDA, device %d, compute %d.%d); graph upload %.1f MiB\n",
                        properties.name, device, properties.major, properties.minor,
                        (g->node_count * sizeof(*g->nodes) + g->edge_count * sizeof(*g->edges)) / 1048576.0);
    return m;
}

static int dispatch(MotifCuda *m, Params params, char *error, size_t error_size) {
    const uint32_t threads = 128;
    uint32_t blocks = (params.task_count + threads - 1) / threads;
    motif_search<<<blocks, threads>>>(m->nodes, m->edges, m->tasks, m->counts,
                                      m->offsets, m->matches, params);
    if (!check_cuda(cudaGetLastError(), "kernel launch", error, error_size)) return 0;
    return check_cuda(cudaDeviceSynchronize(), "kernel execution", error, error_size);
}

int motif_backend_run(void *backend, const Task *tasks, size_t n, int mode, int induced,
                      uint64_t *counts, motif_gpu_emit emit, void *context,
                      char *error, size_t error_size) {
    MotifCuda *m = static_cast<MotifCuda *>(backend);
    if (!n || n > MOTIF_GPU_TASK_BATCH) {
        snprintf(error, error_size, "Invalid CUDA task batch size");
        return 0;
    }
    if (!check_cuda(cudaMemcpy(m->tasks, tasks, n * sizeof(*tasks), cudaMemcpyHostToDevice),
                    "task upload", error, error_size)) return 0;
    Params params = {static_cast<uint32_t>(n), static_cast<uint32_t>(mode),
                     static_cast<uint32_t>(induced), 0, 0, 0};
    if (!dispatch(m, params, error, error_size) ||
        !check_cuda(cudaMemcpy(counts, m->counts, n * sizeof(*counts), cudaMemcpyDeviceToHost),
                    "count download", error, error_size)) return 0;
    if (!emit) return 1;
    /* Allocate export buffers only when requested, leaving count-only runs
     * with just the graph and fixed-size task/count buffers in device memory. */
    if (!m->host_offsets) {
        m->host_offsets = static_cast<uint64_t *>(malloc(MOTIF_GPU_TASK_BATCH * sizeof(uint64_t)));
        m->host_matches = static_cast<Match *>(malloc(MOTIF_GPU_OUTPUT_PAGE * sizeof(Match)));
        if (!m->host_offsets || !m->host_matches) {
            snprintf(error, error_size, "Cannot allocate CUDA export staging buffers");
            return 0;
        }
        if (!check_cuda(cudaMalloc(reinterpret_cast<void **>(&m->matches), MOTIF_GPU_OUTPUT_PAGE * sizeof(Match)),
                        "export allocation", error, error_size)) return 0;
    }
    uint64_t total = 0;
    for (size_t i = 0; i < n; ++i) {
        m->host_offsets[i] = total;
        if (counts[i] > INT64_MAX - total) {
            snprintf(error, error_size, "Motif count exceeds signed 64-bit output range");
            return 0;
        }
        total += counts[i];
    }
    if (!check_cuda(cudaMemcpy(m->offsets, m->host_offsets, n * sizeof(uint64_t), cudaMemcpyHostToDevice),
                    "export offsets upload", error, error_size)) return 0;
    params.emit = 1;
    for (uint64_t start = 0; start < total;) {
        uint64_t remaining = total - start;
        uint64_t length = remaining < MOTIF_GPU_OUTPUT_PAGE ? remaining : MOTIF_GPU_OUTPUT_PAGE;
        params.output_begin = start;
        params.output_end = start + length;
        if (!dispatch(m, params, error, error_size) ||
            !check_cuda(cudaMemcpy(m->host_matches, m->matches, length * sizeof(Match), cudaMemcpyDeviceToHost),
                        "motif download", error, error_size)) return 0;
        emit(m->host_matches, static_cast<size_t>(length), context);
        start += length;
    }
    return 1;
}
