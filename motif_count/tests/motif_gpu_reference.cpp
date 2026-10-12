/* Test-only host interpreter for the shared kernel. It exercises the CUDA
 * graph adapter and C/C++ ABI without claiming to execute CUDA hardware. */
#include "motif_gpu.h"
#include <stdlib.h>
#include <vector>

using Node = motif_gpu_node;
using Task = motif_gpu_task;
using Match = motif_gpu_match;
using Params = motif_gpu_params;
using gpu_uint = uint32_t;
using gpu_count = uint64_t;
enum { EXC = GPU_EXC, INH = GPU_INH, IN_EXC = GPU_IN_EXC, OTHER = GPU_OTHER };
#define GPU_INLINE inline
#define GPU_DEVICE
#define GPU_THREAD
#define GPU_PARAM
#include "motif_gpu_search.h"

void *motif_backend_create(const motif_gpu_graph *g, int quiet, char *, size_t) {
    if (!quiet) fprintf(stderr, "Kernel reference: CPU execution, no GPU hardware\n");
    return const_cast<motif_gpu_graph *>(g);
}

void motif_backend_destroy(void *) {}

int motif_backend_run(void *backend, const Task *tasks, size_t n, int mode, int induced,
                      uint64_t *counts, motif_gpu_emit emit, void *context,
                      char *error, size_t error_size) {
    const auto *g = static_cast<const motif_gpu_graph *>(backend);
    Params p = {static_cast<uint32_t>(n), static_cast<uint32_t>(mode),
                static_cast<uint32_t>(induced), 0, 0, 0};
    for (size_t i = 0; i < n; ++i)
        motif_search_task(g->nodes, g->edges, tasks, counts, nullptr, nullptr, p, i);
    if (!emit) return 1;
    std::vector<uint64_t> offsets(n), replay_counts(n);
    uint64_t total = 0;
    for (size_t i = 0; i < n; ++i) {
        offsets[i] = total;
        if (counts[i] > INT64_MAX - total) {
            snprintf(error, error_size, "Motif count exceeds signed 64-bit output range");
            return 0;
        }
        total += counts[i];
    }
    std::vector<Match> matches(MOTIF_GPU_OUTPUT_PAGE);
    p.emit = 1;
    for (uint64_t start = 0; start < total;) {
        uint64_t remaining = total - start;
        size_t length = remaining < MOTIF_GPU_OUTPUT_PAGE ? remaining : MOTIF_GPU_OUTPUT_PAGE;
        p.output_begin = start; p.output_end = start + length;
        for (size_t i = 0; i < n; ++i)
            motif_search_task(g->nodes, g->edges, tasks, replay_counts.data(),
                              offsets.data(), matches.data(), p, i);
        emit(matches.data(), length, context);
        start += length;
    }
    return 1;
}
