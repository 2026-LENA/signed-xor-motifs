/* Shared Metal/CUDA search. The wrapper supplies GPU address-space and
 * inline macros, fixed-width scalar aliases, and Node/Task/Match/Params types.
 * This header also supports host execution for kernel regression tests. */
#ifndef MOTIF_GPU_SEARCH_H
#define MOTIF_GPU_SEARCH_H

GPU_INLINE bool has(GPU_DEVICE const Node *nodes, GPU_DEVICE const int *edges,
                int s, int t, gpu_uint kind) {
    gpu_uint lo = 0, hi = nodes[s].count[kind], off = nodes[s].offset[kind];
    while (lo < hi) {
        gpu_uint mid = lo + (hi - lo) / 2;
        int v = edges[off + mid];
        if (v == t) return true;
        if (v < t) lo = mid + 1; else hi = mid;
    }
    return false;
}

/* Roles here are n1,n2,n3,n4,n5,n6,F+,F-. Allen's XOR is n5,
 * and its inhibitory pivot is n6, matching the other two counters. */
GPU_INLINE gpu_uint expected(gpu_uint s, gpu_uint t, gpu_uint sign5, gpu_uint variants) {
    if ((s == 0 && (t == 1 || t == 5)) ||
        (s == 2 && (t == 3 || t == 5)) ||
        ((s == 1 || s == 3) && t == 4)) return 1;
    if (s == 5 && (t == 1 || t == 3)) return 2;
    if (s == 4 && (t == 6 || t == 7)) return sign5;
    if ((variants & 1) && ((s == 0 && t == 6) || (s == 2 && t == 7))) return 1;
    if ((variants & 2) && ((s == 1 && t == 6) || (s == 3 && t == 7))) return 1;
    return 0;
}

GPU_INLINE bool induced_ok(GPU_DEVICE const Node *nodes, GPU_DEVICE const int *edges,
                       GPU_THREAD const int *ids, gpu_uint n, gpu_uint sign5, gpu_uint variants) {
    for (gpu_uint i = 0; i < n; ++i) {
        for (gpu_uint j = 0; j < n; ++j) {
            if (i == j) continue;
            /* The six-node core was already checked before feedback search. */
            if (n == 8 && i < 6 && j < 6) continue;
            gpu_uint want = expected(i, j, sign5, variants);
            bool e = has(nodes, edges, ids[i], ids[j], EXC);
            bool h = has(nodes, edges, ids[i], ids[j], INH);
            bool o = has(nodes, edges, ids[i], ids[j], OTHER);
            if (want == 1) { if (!e || h || o) return false; }
            else if (want == 2) { if (!h || e || o) return false; }
            else if (e || h || o) return false;
        }
    }
    return true;
}

GPU_INLINE void motif_search_task(GPU_DEVICE const Node *nodes,
                         GPU_DEVICE const int *edges,
                         GPU_DEVICE const Task *tasks,
                         GPU_DEVICE gpu_count *counts,
                         GPU_DEVICE const gpu_count *offsets,
                         GPU_DEVICE Match *matches,
                         GPU_PARAM const Params &p,
                         gpu_uint tid) {
    if (tid >= p.task_count) return;
    Task task = tasks[tid];
    int n6 = task.anchor, n2 = task.n2, n4 = task.n4, n1 = task.n1;
    bool allen = p.mode == 0;
    gpu_uint c_source = n2, c_other = n4;
    /* Preserve the CPU's duplicate-edge multiplicity: variants iterate the
     * shorter outgoing list; Allen always iterates E2's outgoing list. */
    if (!allen && nodes[n2].count[EXC] > nodes[n4].count[EXC]) {
        c_source = n4; c_other = n2;
    }
    gpu_count found = 0;
    for (gpu_uint b = task.begin; b < task.end; ++b) {
        int n3 = edges[nodes[n4].offset[IN_EXC] + b];
        if (n3 == n1) continue;
        if (allen) {
            if (nodes[n3].identity != 1) continue;
        } else if (n3 == n2 || n3 == n6) continue;
        if (!has(nodes, edges, n3, n6, EXC)) continue;
        for (gpu_uint c = 0; c < nodes[c_source].count[EXC]; ++c) {
            int n5 = edges[nodes[c_source].offset[EXC] + c];
            if (!has(nodes, edges, c_other, n5, EXC)) continue;
            if (n5 == n1 || n5 == n2 || n5 == n3 || n5 == n4 || n5 == n6) continue;
            gpu_uint sign5 = nodes[n5].identity;
            if (allen ? sign5 != 1 : sign5 == 0) continue;
            int ids[8] = {n1,n2,n3,n4,n5,n6,0,0};
            /* Reject extra edges in the six-node core before the costly
             * feedback search. A/B only adds edges to the two new nodes. */
            if ((allen || p.induced) && !induced_ok(nodes, edges, ids, 6, sign5, 0)) continue;
            gpu_uint kind = sign5 == 1 ? EXC : INH;
            gpu_uint off = nodes[n5].offset[kind], len = nodes[n5].count[kind];
            if (!allen && !p.induced && !p.emit) {
                /* For count-only, non-induced A/B motifs, F+ and F- have
                 * opposite pure identities and cannot coincide. Count the
                 * eligible F- lists once, then apply inclusion/exclusion for
                 * each F+. This replaces a quadratic feedback-pair loop with
                 * two linear scans, preserving AB and duplicate multiplicity. */
                gpu_count fm_a = 0, fm_b = 0, fm_ab = 0;
                for (gpu_uint mi = 0; mi < len; ++mi) {
                    int fm = edges[off + mi];
                    if (fm == n1 || fm == n2 || fm == n3 || fm == n4 || fm == n5 || fm == n6) continue;
                    if (nodes[fm].identity != 2) continue;
                    bool a = has(nodes, edges, n3, fm, EXC);
                    bool b_ok = has(nodes, edges, n4, fm, EXC);
                    fm_a += gpu_count(a); fm_b += gpu_count(b_ok); fm_ab += gpu_count(a && b_ok);
                }
                if (!fm_a && !fm_b) continue;
                for (gpu_uint fi = 0; fi < len; ++fi) {
                    int fp = edges[off + fi];
                    if (fp == n1 || fp == n2 || fp == n3 || fp == n4 || fp == n5 || fp == n6) continue;
                    if (nodes[fp].identity != 1) continue;
                    bool a = has(nodes, edges, n1, fp, EXC);
                    bool b_ok = has(nodes, edges, n2, fp, EXC);
                    found += (a ? fm_a : 0) + (b_ok ? fm_b : 0) - (a && b_ok ? fm_ab : 0);
                }
                continue;
            }
            for (gpu_uint fi = 0; fi < len; ++fi) {
                int fp = edges[off + fi];
                if (fp == n1 || fp == n2 || fp == n3 || fp == n4 || fp == n5 || fp == n6) continue;
                if (nodes[fp].identity != 1) continue;
                bool a_fp = !allen && has(nodes, edges, n1, fp, EXC);
                bool b_fp = has(nodes, edges, n2, fp, EXC);
                if (!a_fp && !b_fp) continue;
                for (gpu_uint mi = 0; mi < len; ++mi) {
                    int fm = edges[off + mi];
                    if (fm == fp || fm == n1 || fm == n2 || fm == n3 ||
                        fm == n4 || fm == n5 || fm == n6) continue;
                    if (nodes[fm].identity != 2) continue;
                    bool a_ok = a_fp && has(nodes, edges, n3, fm, EXC);
                    bool b_ok = b_fp && has(nodes, edges, n4, fm, EXC);
                    if (!a_ok && !b_ok) continue;
                    gpu_uint variants = gpu_uint(a_ok) | (gpu_uint(b_ok) << 1);
                    ids[6] = fp; ids[7] = fm;
                    if ((allen || p.induced) && !induced_ok(nodes, edges, ids, 8, sign5, variants)) continue;
                    if (p.emit) {
                        gpu_count index = offsets[tid] + found;
                        if (index >= p.output_begin && index < p.output_end) {
                            for (gpu_uint k = 0; k < 8; ++k)
                                matches[index - p.output_begin].ids[k] = ids[k];
                            matches[index - p.output_begin].variants = variants;
                        }
                    }
                    ++found;
                }
            }
        }
    }
    counts[tid] = found;
}

#endif
