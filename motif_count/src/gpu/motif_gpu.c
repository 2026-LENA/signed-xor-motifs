#include "motif_gpu.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define TASK_BATCH MOTIF_GPU_TASK_BATCH
#define N3_CHUNK 8

static int compare_int(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

int motif_gpu_pack(motif_gpu_graph *g, const motif_gpu_node_view *v, size_t n,
                   int preserve_inh_order, char *error, size_t error_size) {
    memset(g, 0, sizeof(*g));
    if (n > INT_MAX) goto too_large;
    size_t total = 0;
    for (size_t i = 0; i < n; ++i) {
        for (int k = 0; k < 4; ++k) {
            if (v[i].counts[k] > UINT32_MAX - total) goto too_large;
            total += v[i].counts[k];
        }
        if (preserve_inh_order) {
            if (v[i].counts[GPU_INH] > UINT32_MAX - total) goto too_large;
            total += v[i].counts[GPU_INH];
        }
    }
    g->nodes = calloc(n ? n : 1, sizeof(*g->nodes));
    g->edges = malloc((total ? total : 1) * sizeof(*g->edges));
    if (!g->nodes || !g->edges) {
        snprintf(error, error_size, "Cannot allocate the GPU graph");
        motif_gpu_free(g);
        return 0;
    }
    g->node_count = n;
    g->edge_count = total;
    size_t off = 0;
    for (size_t i = 0; i < n; ++i) {
        motif_gpu_node *node = &g->nodes[i];
        node->identity = v[i].identity;
        for (int k = 0; k < GPU_LISTS; ++k) {
            if (k == GPU_INH_ORDER && !preserve_inh_order) {
                node->offset[k] = node->offset[GPU_INH];
                node->count[k] = node->count[GPU_INH];
                continue;
            }
            int source = k == GPU_INH_ORDER ? GPU_INH : k;
            size_t count = v[i].counts[source];
            node->offset[k] = (uint32_t)off;
            node->count[k] = (uint32_t)count;
            if (count) memcpy(g->edges + off, v[i].lists[source], count * sizeof(int));
            for (size_t j = 0; j < count; ++j) {
                if (g->edges[off + j] < 0 || (size_t)g->edges[off + j] >= n) {
                    snprintf(error, error_size, "Invalid compact node ID in GPU graph");
                    motif_gpu_free(g);
                    return 0;
                }
            }
            if (k != GPU_INH_ORDER && count > 1)
                qsort(g->edges + off, count, sizeof(int), compare_int);
            off += count;
        }
    }
    return 1;
too_large:
    snprintf(error, error_size, "Graph exceeds the GPU backend's 32-bit CSR limits");
    return 0;
}

void motif_gpu_free(motif_gpu_graph *g) {
    free(g->nodes);
    free(g->edges);
    memset(g, 0, sizeof(*g));
}

static int has_edge(const motif_gpu_graph *g, int s, int t, int kind) {
    const motif_gpu_node *n = &g->nodes[s];
    const int32_t *a = g->edges + n->offset[kind];
    size_t lo = 0, hi = n->count[kind];
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (a[mid] == t) return 1;
        if (a[mid] < t) lo = mid + 1; else hi = mid;
    }
    return 0;
}

static int flush_tasks(void *backend, motif_gpu_task *tasks, size_t n,
                       int mode, int induced, uint64_t *totals,
                       motif_gpu_emit emit, void *context,
                       char *error, size_t error_size) {
    uint64_t found[TASK_BATCH];
    if (!motif_backend_run(backend, tasks, n, mode, induced, found,
                        emit, context, error, error_size)) return 0;
    for (size_t i = 0; i < n; ++i) {
        uint64_t *total = &totals[tasks[i].group];
        if (found[i] > INT64_MAX - *total) {
            snprintf(error, error_size, "Motif count exceeds signed 64-bit output range");
            return 0;
        }
        *total += found[i];
    }
    return 1;
}

int motif_gpu_count(const motif_gpu_graph *g, const int *anchors, size_t nanchors,
                    int mode, int induced, uint64_t *counts,
                    motif_gpu_emit emit, void *context, int quiet,
                    char *error, size_t error_size) {
    if (nanchors > UINT32_MAX) {
        snprintf(error, error_size, "Too many GPU anchors");
        return 0;
    }
    memset(counts, 0, nanchors * sizeof(*counts));
    void *backend = motif_backend_create(g, quiet, error, error_size);
    if (!backend) return 0;
    motif_gpu_task *tasks = malloc(TASK_BATCH * sizeof(*tasks));
    if (!tasks) {
        snprintf(error, error_size, "Cannot allocate GPU work queue");
        motif_backend_destroy(backend);
        return 0;
    }
    size_t pending = 0;
    int ok = 1;
    uint64_t jobs = 0;
    struct timespec stamp;
    clock_gettime(CLOCK_MONOTONIC, &stamp);
    double last_progress = stamp.tv_sec + stamp.tv_nsec * 1e-9;
    /* Split each inhibitory anchor into target pairs, input neurons, and small
     * slices of the second input list. This exposes parallelism even for a
     * single --pivot, without materializing all possible eight-node tuples. */
    for (size_t r = 0; r < nanchors && ok; ++r) {
        int anchor = anchors[r];
        const motif_gpu_node *n6 = &g->nodes[anchor];
        const int32_t *out = g->edges + n6->offset[GPU_INH_ORDER];
        for (uint32_t i = 0; i < n6->count[GPU_INH_ORDER] && ok; ++i) {
            int n2 = out[i];
            if (mode == GPU_ALLEN && g->nodes[n2].identity != 1) continue;
            const motif_gpu_node *two = &g->nodes[n2];
            for (uint32_t j = mode == GPU_ALLEN ? i + 1 : 0;
                 j < n6->count[GPU_INH_ORDER] && ok; ++j) {
                int n4 = out[j];
                if (mode == GPU_VARIANTS && (i == j || n2 == n4)) continue;
                if (mode == GPU_ALLEN && g->nodes[n4].identity != 1) continue;
                const motif_gpu_node *four = &g->nodes[n4];
                if (!two->count[GPU_EXC] || !four->count[GPU_EXC]) continue;
                for (uint32_t a = 0; a < two->count[GPU_IN_EXC] && ok; ++a) {
                    int n1 = g->edges[two->offset[GPU_IN_EXC] + a];
                    if (mode == GPU_ALLEN) {
                        if (g->nodes[n1].identity != 1 ||
                            !has_edge(g, n1, anchor, GPU_EXC)) continue;
                    } else if (n1 == n4 || n1 == anchor ||
                               !has_edge(g, n1, anchor, GPU_EXC)) continue;
                    for (uint32_t b = 0; b < four->count[GPU_IN_EXC] && ok;) {
                        uint32_t end = b + (four->count[GPU_IN_EXC] - b < N3_CHUNK ?
                                           four->count[GPU_IN_EXC] - b : N3_CHUNK);
                        tasks[pending++] = (motif_gpu_task){anchor, n2, n4, n1,
                                                           b, end, (uint32_t)r};
                        b = end;
                        if (pending == TASK_BATCH) {
                            ok = flush_tasks(backend, tasks, pending, mode, induced,
                                             counts, emit, context, error, error_size);
                            pending = 0;
                            jobs += TASK_BATCH;
                            if (!quiet) {
                                clock_gettime(CLOCK_MONOTONIC, &stamp);
                                double now = stamp.tv_sec + stamp.tv_nsec * 1e-9;
                                if (now - last_progress >= 5.0) {
                                    fprintf(stderr, "GPU search: anchor %zu/%zu, %llu jobs completed\n",
                                            r + 1, nanchors, (unsigned long long)jobs);
                                    last_progress = now;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (ok && pending)
        ok = flush_tasks(backend, tasks, pending, mode, induced,
                         counts, emit, context, error, error_size);
    free(tasks);
    motif_backend_destroy(backend);
    return ok;
}
