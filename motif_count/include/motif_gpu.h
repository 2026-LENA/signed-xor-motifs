#ifndef MOTIF_GPU_H
#define MOTIF_GPU_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#if defined(MOTIF_USE_METAL) && defined(MOTIF_USE_CUDA)
#error "Select one GPU backend per build: Metal or CUDA"
#endif
#if defined(MOTIF_USE_METAL) || defined(MOTIF_USE_CUDA)
#define MOTIF_USE_GPU 1
#endif

/* The ordinary one-file C builds still work. Explicit GPU requests must never
 * silently run on the CPU. Link either backend to opt in. */
static inline int motif_gpu_check_build(int requested) {
#ifndef MOTIF_USE_GPU
    if (requested) {
        fprintf(stderr, "GPU support is not compiled in. Run 'make metal' (Apple Silicon)\n"
                        "or 'make cuda' (NVIDIA) in motif_count/ and use the executable\n"
                        "in build/metal/ or build/cuda/ with --gpu.\n");
        return 0;
    }
#else
    (void)requested;
#endif
    return 1;
}

#ifdef MOTIF_USE_GPU
#define MOTIF_GPU_TASK_BATCH 32768
#define MOTIF_GPU_OUTPUT_PAGE 65536
enum { GPU_EXC, GPU_INH, GPU_IN_EXC, GPU_OTHER, GPU_INH_ORDER, GPU_LISTS };
enum { GPU_ALLEN, GPU_VARIANTS };

/* Fixed-width layouts are shared with the GPU kernels. CSR offsets address
 * entries, not bytes. Adjacencies retain duplicates to match the CPU counters. */
typedef struct {
    uint32_t offset[GPU_LISTS], count[GPU_LISTS];
    uint32_t identity;  /* 1: excitatory, 2: inhibitory, 0: ambiguous / unknown */
} motif_gpu_node;

typedef struct {
    motif_gpu_node *nodes;
    int32_t *edges;
    size_t node_count, edge_count;
} motif_gpu_graph;

typedef struct {
    const int *lists[4];
    size_t counts[4];
    uint32_t identity;
} motif_gpu_node_view;

typedef struct {
    int32_t anchor, n2, n4, n1;
    uint32_t begin, end, group;
} motif_gpu_task;

typedef struct {
    int32_t ids[8]; /* n1,n2,n3,n4,n5,n6,F+,F- */
    uint32_t variants; /* 1 = A, 2 = B, 3 = AB; Allen uses B */
} motif_gpu_match;

typedef struct {
    uint32_t task_count, mode, induced, emit;
    uint64_t output_begin, output_end;
} motif_gpu_params;

typedef void (*motif_gpu_emit)(const motif_gpu_match *, size_t, void *);

#ifdef __cplusplus
extern "C" {
#endif
int motif_gpu_pack(motif_gpu_graph *, const motif_gpu_node_view *, size_t,
                   int preserve_inh_order, char *, size_t);
void motif_gpu_free(motif_gpu_graph *);
int motif_gpu_count(const motif_gpu_graph *, const int *anchors, size_t nanchors,
                    int mode, int induced, uint64_t *counts,
                    motif_gpu_emit emit, void *context, int quiet,
                    char *error, size_t error_size);

/* Internal backend interface. Every batch completes before host reads results. */
void *motif_backend_create(const motif_gpu_graph *, int quiet, char *, size_t);
int motif_backend_run(void *, const motif_gpu_task *, size_t, int mode, int induced,
                    uint64_t *, motif_gpu_emit, void *, char *, size_t);
void motif_backend_destroy(void *);
#ifdef __cplusplus
}
#endif
#endif
#endif
