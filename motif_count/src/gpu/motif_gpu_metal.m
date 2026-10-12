#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "motif_gpu.h"
#include "motif_gpu_source.h"
#include <limits.h>

#define OUTPUT_PAGE MOTIF_GPU_OUTPUT_PAGE
#define TASK_BATCH MOTIF_GPU_TASK_BATCH

_Static_assert(sizeof(motif_gpu_node) == 44, "Metal Node layout");
_Static_assert(sizeof(motif_gpu_task) == 28, "Metal Task layout");
_Static_assert(sizeof(motif_gpu_match) == 36, "Metal Match layout");
_Static_assert(sizeof(motif_gpu_params) == 32, "Metal Params layout");

@interface MotifMetal : NSObject
@property(nonatomic, strong) id<MTLDevice> device;
@property(nonatomic, strong) id<MTLCommandQueue> queue;
@property(nonatomic, strong) id<MTLComputePipelineState> pipeline;
@property(nonatomic, strong) id<MTLBuffer> nodes, edges, tasks, counts, offsets, matches;
@end
@implementation MotifMetal
@end

static void describe_error(char *error, size_t size, NSString *message) {
    snprintf(error, size, "%s", message.UTF8String ?: "Unknown Metal error");
}

static id<MTLBuffer> upload(id<MTLDevice> device, const void *bytes, size_t size) {
    if (!size) return [device newBufferWithLength:4 options:MTLResourceStorageModeShared];
    if (size > device.maxBufferLength) return nil;
    return [device newBufferWithBytes:bytes length:size options:MTLResourceStorageModeShared];
}

void *motif_backend_create(const motif_gpu_graph *g, int quiet, char *error, size_t error_size) {
    @autoreleasepool {
        MotifMetal *m = [MotifMetal new];
        m.device = MTLCreateSystemDefaultDevice();
        if (!m.device) {
            describe_error(error, error_size, @"No Metal GPU is accessible to this process (check sandbox/session permissions)");
            return NULL;
        }
        if (![m.device supportsFamily:MTLGPUFamilyApple3]) {
            describe_error(error, error_size, @"The Metal backend requires an Apple Silicon GPU");
            return NULL;
        }
        NSError *err = nil;
        MTLCompileOptions *options = [MTLCompileOptions new];
        options.languageVersion = MTLLanguageVersion2_3;
        id<MTLLibrary> library = [m.device newLibraryWithSource:@(motif_gpu_source)
                                                       options:options error:&err];
        if (!library) {
            describe_error(error, error_size, err.localizedDescription);
            return NULL;
        }
        m.pipeline = [m.device newComputePipelineStateWithFunction:[library newFunctionWithName:@"motif_search"]
                                                            error:&err];
        m.queue = [m.device newCommandQueue];
        m.nodes = upload(m.device, g->nodes, g->node_count * sizeof(*g->nodes));
        m.edges = upload(m.device, g->edges, g->edge_count * sizeof(*g->edges));
        m.tasks = [m.device newBufferWithLength:TASK_BATCH * sizeof(motif_gpu_task) options:MTLResourceStorageModeShared];
        m.counts = [m.device newBufferWithLength:TASK_BATCH * sizeof(uint64_t) options:MTLResourceStorageModeShared];
        m.offsets = [m.device newBufferWithLength:TASK_BATCH * sizeof(uint64_t) options:MTLResourceStorageModeShared];
        m.matches = [m.device newBufferWithLength:OUTPUT_PAGE * sizeof(motif_gpu_match) options:MTLResourceStorageModeShared];
        if (!m.pipeline || !m.queue || !m.nodes || !m.edges || !m.tasks ||
            !m.counts || !m.offsets || !m.matches) {
            describe_error(error, error_size, err ? err.localizedDescription : @"Cannot allocate Metal buffers for this graph");
            return NULL;
        }
        if (!quiet) fprintf(stderr, "GPU: %s (Metal); graph upload %.1f MiB\n", m.device.name.UTF8String,
                            (g->node_count * sizeof(*g->nodes) + g->edge_count * sizeof(*g->edges)) / 1048576.0);
        return (__bridge_retained void *)m;
    }
}

static int dispatch(MotifMetal *m, motif_gpu_params params, char *error, size_t error_size) {
    id<MTLCommandBuffer> command = [m.queue commandBuffer];
    id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
    if (!command || !encoder) {
        describe_error(error, error_size, @"Cannot create Metal command buffer");
        return 0;
    }
    [encoder setComputePipelineState:m.pipeline];
    [encoder setBuffer:m.nodes offset:0 atIndex:0];
    [encoder setBuffer:m.edges offset:0 atIndex:1];
    [encoder setBuffer:m.tasks offset:0 atIndex:2];
    [encoder setBuffer:m.counts offset:0 atIndex:3];
    [encoder setBuffer:m.offsets offset:0 atIndex:4];
    [encoder setBuffer:m.matches offset:0 atIndex:5];
    [encoder setBytes:&params length:sizeof(params) atIndex:6];
    NSUInteger width = MIN((NSUInteger)64, m.pipeline.maxTotalThreadsPerThreadgroup);
    [encoder dispatchThreads:MTLSizeMake(params.task_count, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
    [encoder endEncoding];
    [command commit];
    [command waitUntilCompleted];
    if (command.status != MTLCommandBufferStatusCompleted) {
        describe_error(error, error_size, command.error.localizedDescription ?: @"Metal command did not complete");
        return 0;
    }
    return 1;
}

int motif_backend_run(void *backend, const motif_gpu_task *tasks, size_t n,
                    int mode, int induced, uint64_t *counts,
                    motif_gpu_emit emit, void *context, char *error, size_t error_size) {
    @autoreleasepool {
        MotifMetal *m = (__bridge MotifMetal *)backend;
        if (!n || n > TASK_BATCH) {
            describe_error(error, error_size, @"Invalid Metal task batch size");
            return 0;
        }
        memcpy(m.tasks.contents, tasks, n * sizeof(*tasks));
        motif_gpu_params p = {(uint32_t)n, (uint32_t)mode, (uint32_t)induced, 0, 0, 0};
        if (!dispatch(m, p, error, error_size)) return 0;
        memcpy(counts, m.counts.contents, n * sizeof(*counts));
        if (!emit) return 1;
        /* A prefix sum gives each thread a disjoint output range, avoiding
         * atomic contention and 32-bit counter overflow. Page large results
         * by replaying this batch instead of ever truncating exported motifs. */
        uint64_t total = 0, *offsets = m.offsets.contents;
        for (size_t i = 0; i < n; ++i) {
            offsets[i] = total;
            if (counts[i] > INT64_MAX - total) {
                describe_error(error, error_size, @"Motif count exceeds signed 64-bit output range");
                return 0;
            }
            total += counts[i];
        }
        p.emit = 1;
        for (uint64_t start = 0; start < total;) {
            uint64_t length = MIN((uint64_t)OUTPUT_PAGE, total - start);
            p.output_begin = start;
            p.output_end = start + length;
            if (!dispatch(m, p, error, error_size)) return 0;
            emit(m.matches.contents, (size_t)length, context);
            start += length;
        }
        return 1;
    }
}

void motif_backend_destroy(void *backend) {
    @autoreleasepool {
        MotifMetal *m = (__bridge_transfer MotifMetal *)backend;
        (void)m;
    }
}
