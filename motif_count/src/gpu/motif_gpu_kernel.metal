#include <metal_stdlib>
using namespace metal;

enum { EXC, INH, IN_EXC, OTHER, INH_ORDER, LISTS };
struct Node { uint offset[LISTS], count[LISTS], identity; };
struct Task { int anchor, n2, n4, n1; uint begin, end, group; };
struct Match { int ids[8]; uint variants; };
struct Params {
    uint task_count, mode, induced, emit;
    ulong output_begin, output_end;
};

using gpu_uint = uint;
using gpu_count = ulong;
#define GPU_INLINE inline
#define GPU_DEVICE device
#define GPU_THREAD thread
#define GPU_PARAM constant
#include "motif_gpu_search.h"

kernel void motif_search(device const Node *nodes [[buffer(0)]],
                         device const int *edges [[buffer(1)]],
                         device const Task *tasks [[buffer(2)]],
                         device ulong *counts [[buffer(3)]],
                         device const ulong *offsets [[buffer(4)]],
                         device Match *matches [[buffer(5)]],
                         constant Params &p [[buffer(6)]],
                         uint tid [[thread_position_in_grid]]) {
    motif_search_task(nodes, edges, tasks, counts, offsets, matches, p, tid);
}
