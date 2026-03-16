#import "MetalUtils.h"
#include <iostream>
#include <stdexcept>

// ============================================================================
// MetalContext implementation
// ============================================================================

MetalContext& MetalContext::get() {
    static MetalContext ctx;
    return ctx;
}

void MetalContext::init(const std::string& metallib_path) {
    device = MTLCreateSystemDefaultDevice();
    if (!device)
        throw std::runtime_error("MetalContext::init: no Metal device found");

    queue = [device newCommandQueue];
    if (!queue)
        throw std::runtime_error("MetalContext::init: failed to create command queue");

    // Load precompiled .metallib
    NSString* path = metallib_path.empty()
        ? @METALLIB_PATH
        : [NSString stringWithUTF8String:metallib_path.c_str()];

    NSError* err = nil;
    NSURL* url = [NSURL fileURLWithPath:path];
    library = [device newLibraryWithURL:url error:&err];
    if (!library) {
        std::string msg = "MetalContext::init: failed to load metallib at '";
        msg += [path UTF8String];
        msg += "': ";
        if (err) msg += [[err localizedDescription] UTF8String];
        throw std::runtime_error(msg);
    }

    std::cerr << "[metal] device: " << [[device name] UTF8String] << std::endl;
    std::cerr << "[metal] metallib: " << [path UTF8String] << std::endl;
}

id<MTLComputePipelineState> MetalContext::getPipeline(const std::string& name) {
    auto it = pipelines.find(name);
    if (it != pipelines.end()) return it->second;

    NSString* ns_name = [NSString stringWithUTF8String:name.c_str()];
    id<MTLFunction> fn = [library newFunctionWithName:ns_name];
    if (!fn)
        throw std::runtime_error("MetalContext::getPipeline: kernel not found: " + name);

    NSError* err = nil;
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:fn error:&err];
    if (!pso) {
        std::string msg = "MetalContext::getPipeline: PSO creation failed for " + name + ": ";
        if (err) msg += [[err localizedDescription] UTF8String];
        throw std::runtime_error(msg);
    }

    pipelines[name] = pso;
    return pso;
}

id<MTLCommandBuffer> MetalContext::createCommandBuffer() {
    id<MTLCommandBuffer> cmd = [queue commandBuffer];
    if (!cmd)
        throw std::runtime_error("MetalContext::createCommandBuffer: failed to create command buffer");
    return cmd;
}

void MetalContext::commitAndWait(id<MTLCommandBuffer> cmd) {
    [cmd commit];
    [cmd waitUntilCompleted];

    if (cmd.status == MTLCommandBufferStatusError) {
        std::string msg = "MetalContext::commitAndWait: command buffer error";
        if (cmd.error)
            msg += ": " + std::string([[cmd.error localizedDescription] UTF8String]);
        throw std::runtime_error(msg);
    }
}

void MetalContext::encodeDispatch(id<MTLCommandBuffer> cmd,
                                  const std::string& kernel,
                                  NSUInteger count,
                                  const std::vector<id<MTLBuffer>>& buffers,
                                  const std::vector<std::vector<uint8_t>>& constants)
{
    encodeDispatch(cmd, kernel, count, buffers, constants, 0);
}

void MetalContext::encodeDispatch(id<MTLCommandBuffer> cmd,
                                  const std::string& kernel,
                                  NSUInteger count,
                                  const std::vector<id<MTLBuffer>>& buffers,
                                  const std::vector<std::vector<uint8_t>>& constants,
                                  NSUInteger forced_threads_per_group)
{
    if (count == 0) return;

    id<MTLComputePipelineState> pso = getPipeline(kernel);
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    [enc setComputePipelineState:pso];

    for (NSUInteger i = 0; i < (NSUInteger)buffers.size(); ++i)
        [enc setBuffer:buffers[i] offset:0 atIndex:i];

    NSUInteger next_slot = (NSUInteger)buffers.size();
    for (const auto& c : constants) {
        [enc setBytes:c.data() length:c.size() atIndex:next_slot++];
    }

    // Use a threadgroup size aligned to the SIMD execution width rather than
    // blindly taking the maximum. Kernels like process_reads_kernel have large
    // private state, so oversized groups can reduce occupancy and increase
    // register spilling on Apple GPUs.
    NSUInteger tg_size = forced_threads_per_group;
    if (tg_size == 0) {
        NSUInteger exec_width = std::max((NSUInteger)pso.threadExecutionWidth, (NSUInteger)1);
        tg_size = exec_width * 4;
        tg_size = std::min(tg_size, (NSUInteger)pso.maxTotalThreadsPerThreadgroup);
        tg_size = std::min(tg_size, count);
        tg_size = std::max(tg_size, exec_width);
    } else {
        tg_size = std::min(tg_size, (NSUInteger)pso.maxTotalThreadsPerThreadgroup);
    }

    MTLSize threads_per_grid = MTLSizeMake(count, 1, 1);
    MTLSize threadgroup_size = MTLSizeMake(tg_size, 1, 1);

    [enc dispatchThreads:threads_per_grid threadsPerThreadgroup:threadgroup_size];
    [enc endEncoding];
}

void MetalContext::dispatch(const std::string& kernel,
                            NSUInteger count,
                            const std::vector<id<MTLBuffer>>& buffers,
                            const std::vector<std::vector<uint8_t>>& constants)
{
    id<MTLCommandBuffer> cmd = createCommandBuffer();
    encodeDispatch(cmd, kernel, count, buffers, constants);
    commitAndWait(cmd);
}

void MetalContext::dispatch(const std::string& kernel,
                            NSUInteger count,
                            const std::vector<id<MTLBuffer>>& buffers,
                            const std::vector<std::vector<uint8_t>>& constants,
                            NSUInteger threads_per_group)
{
    id<MTLCommandBuffer> cmd = createCommandBuffer();
    encodeDispatch(cmd, kernel, count, buffers, constants, threads_per_group);
    commitAndWait(cmd);
}
