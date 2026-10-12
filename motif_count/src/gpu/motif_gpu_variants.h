/* Adapter shared by the FlyWire and C. elegans counters. Include after their
 * graph and string-ID tables have been declared. */
#ifdef MOTIF_USE_GPU
#include <limits.h>

static void emit_gpu_variants(const motif_gpu_match *matches, size_t n, void *context) {
    FILE *fp = context;
    for (size_t i = 0; i < n; ++i) {
        const motif_gpu_match *m = &matches[i];
        fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
                id_to_name[m->ids[0]], id_to_name[m->ids[1]],
                id_to_name[m->ids[2]], id_to_name[m->ids[3]],
                id_to_name[m->ids[4]], id_to_name[m->ids[5]],
                id_to_name[m->ids[6]], id_to_name[m->ids[7]],
                m->variants == 3 ? "AB" : m->variants == 1 ? "A" : "B");
    }
}

static int count_gpu_variants(const int *anchors, size_t nanchors, int induced,
                              FILE *out, int quiet, long long *raw) {
    char error[2048] = "";
    motif_gpu_node_view *views = calloc(N_nodes ? (size_t)N_nodes : 1, sizeof(*views));
    uint64_t *counts = calloc(nanchors ? nanchors : 1, sizeof(*counts));
    if (!views || !counts) {
        fprintf(stderr, "Cannot allocate GPU graph views\n");
        free(views); free(counts);
        return 0;
    }
    for (int i = 0; i < N_nodes; ++i) {
        views[i] = (motif_gpu_node_view){
            {graph[i].out_exc, graph[i].out_inh, graph[i].in_exc, graph[i].out_other},
            {graph[i].out_exc_n, graph[i].out_inh_n, graph[i].in_exc_n, graph[i].out_other_n},
            graph[i].is_exc && !graph[i].is_inh ? 1u :
            graph[i].is_inh && !graph[i].is_exc ? 2u : 0u
        };
    }
    motif_gpu_graph packed;
    int ok = motif_gpu_pack(&packed, views, N_nodes, 0, error, sizeof(error));
    free(views);
    if (ok) {
        ok = motif_gpu_count(&packed, anchors, nanchors, GPU_VARIANTS, induced,
                             counts, out ? emit_gpu_variants : NULL, out, quiet,
                             error, sizeof(error));
        motif_gpu_free(&packed);
    }
    uint64_t sum = 0;
    if (ok) {
        for (size_t i = 0; i < nanchors; ++i) {
            if (counts[i] > LLONG_MAX - sum) {
                snprintf(error, sizeof(error), "Motif count exceeds signed 64-bit output range");
                ok = 0;
                break;
            }
            sum += counts[i];
        }
    }
    free(counts);
    if (!ok) fprintf(stderr, "GPU error: %s\n", error);
    else *raw = (long long)sum;
    return ok;
}
#endif
