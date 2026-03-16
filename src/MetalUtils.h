#pragma once

#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <algorithm>
#include <cassert>

// ============================================================================
// MetalContext — singleton holding device, queue, library, and pipeline cache
// ============================================================================

struct MetalContext {
    id<MTLDevice>       device  = nil;
    id<MTLCommandQueue> queue   = nil;
    id<MTLLibrary>      library = nil;
    std::unordered_map<std::string, id<MTLComputePipelineState>> pipelines;

    static MetalContext& get();
    void init(const std::string& metallib_path = "");

    /// Returns (or creates) a compiled compute pipeline for the named kernel.
    id<MTLComputePipelineState> getPipeline(const std::string& name);

    id<MTLCommandBuffer> createCommandBuffer();
    void commitAndWait(id<MTLCommandBuffer> cmd);

    /// Encode + commit + wait for a compute dispatch.
    /// @param kernel     Name of the MSL kernel function.
    /// @param count      Total number of threads to dispatch (1-D grid).
    /// @param buffers    Array of MTLBuffer to bind at slots 0, 1, 2, …
    /// @param constants  Byte arrays for setBytes at slots after buffers.
    void dispatch(const std::string& kernel,
                  NSUInteger count,
                  const std::vector<id<MTLBuffer>>& buffers,
                  const std::vector<std::vector<uint8_t>>& constants = {});

    void dispatch(const std::string& kernel,
                  NSUInteger count,
                  const std::vector<id<MTLBuffer>>& buffers,
                  const std::vector<std::vector<uint8_t>>& constants,
                  NSUInteger threads_per_group);

    /// Encode a compute dispatch onto an existing command buffer.
    void encodeDispatch(id<MTLCommandBuffer> cmd,
                        const std::string& kernel,
                        NSUInteger count,
                        const std::vector<id<MTLBuffer>>& buffers,
                        const std::vector<std::vector<uint8_t>>& constants = {});

    void encodeDispatch(id<MTLCommandBuffer> cmd,
                        const std::string& kernel,
                        NSUInteger count,
                        const std::vector<id<MTLBuffer>>& buffers,
                        const std::vector<std::vector<uint8_t>>& constants,
                        NSUInteger threads_per_group);
};

// ============================================================================
// MetalBuffer<T> — typed wrapper around MTLBuffer (shared storage mode)
//
// On Apple Silicon the buffer is directly addressable from CPU and GPU via
// the same physical memory — no copies needed.
// ============================================================================

template<typename T>
struct MetalBuffer {
    id<MTLBuffer> buf = nil;
    size_t        _size = 0;  // element count (not byte count)

    MetalBuffer() = default;
    explicit MetalBuffer(size_t n, T init_val = T{}) { resize(n, init_val); }

    // Prevent accidental copies (buffers own GPU memory)
    MetalBuffer(const MetalBuffer&) = delete;
    MetalBuffer& operator=(const MetalBuffer&) = delete;

    MetalBuffer(MetalBuffer&& o) noexcept
        : buf(o.buf), _size(o._size) { o.buf = nil; o._size = 0; }
    MetalBuffer& operator=(MetalBuffer&& o) noexcept {
        buf = o.buf; _size = o._size; o.buf = nil; o._size = 0; return *this;
    }

    ~MetalBuffer() { buf = nil; }

    // -----------------------------------------------------------------------

    size_t size() const { return _size; }
    bool   empty() const { return _size == 0; }

    /// CPU-visible pointer (valid on Apple Silicon unified memory)
    T*       data()       { return buf ? static_cast<T*>(buf.contents) : nullptr; }
    const T* data() const { return buf ? static_cast<const T*>(buf.contents) : nullptr; }

    T& operator[](size_t i)       { return data()[i]; }
    const T& operator[](size_t i) const { return data()[i]; }

    id<MTLBuffer> metalBuffer() const { return buf; }

    // -----------------------------------------------------------------------
    // Resize — preserves existing data up to min(old, new) elements
    void resize(size_t n, T init_val = T{}) {
        if (n == _size) return;

        if (n == 0) { buf = nil; _size = 0; return; }

        id<MTLDevice> dev = MetalContext::get().device;
        NSUInteger bytes = static_cast<NSUInteger>(n) * sizeof(T);

        id<MTLBuffer> new_buf = [dev newBufferWithLength:bytes
                                                 options:MTLResourceStorageModeShared];
        if (!new_buf)
            throw std::runtime_error("MetalBuffer::resize: MTLBuffer allocation failed ("
                                     + std::to_string(bytes) + " bytes)");

        T* dst = static_cast<T*>(new_buf.contents);

        if (n > _size) {
            // Copy old data, fill new region with init_val
            size_t old_bytes = _size * sizeof(T);
            if (buf && _size > 0)
                std::memcpy(dst, buf.contents, old_bytes);
            std::fill(dst + _size, dst + n, init_val);
        } else {
            // Shrink: copy only n elements
            if (buf && _size > 0)
                std::memcpy(dst, buf.contents, n * sizeof(T));
        }

        buf   = new_buf;
        _size = n;
    }

    // Assign from a host range [begin, end)
    template<typename Iter>
    void assign(Iter begin, Iter end) {
        size_t n = static_cast<size_t>(std::distance(begin, end));
        resize_raw(n);
        T* dst = data();
        size_t i = 0;
        for (auto it = begin; it != end; ++it, ++i)
            dst[i] = static_cast<T>(*it);
    }

    // Fill all elements with value
    void fill(T val) {
        T* p = data();
        std::fill(p, p + _size, val);
    }

    // Copy to host vector
    std::vector<T> to_host() const {
        std::vector<T> out(_size);
        if (_size > 0 && buf)
            std::memcpy(out.data(), buf.contents, _size * sizeof(T));
        return out;
    }

    // Copy from host vector
    void from_host(const std::vector<T>& h) {
        resize_raw(h.size());
        if (!h.empty())
            std::memcpy(data(), h.data(), h.size() * sizeof(T));
    }

private:
    // Resize without initialization (just allocate)
    void resize_raw(size_t n) {
        if (n == _size) return;
        if (n == 0) { buf = nil; _size = 0; return; }

        id<MTLDevice> dev = MetalContext::get().device;
        NSUInteger bytes = static_cast<NSUInteger>(n) * sizeof(T);
        buf = [dev newBufferWithLength:bytes options:MTLResourceStorageModeShared];
        if (!buf)
            throw std::runtime_error("MetalBuffer::resize_raw: allocation failed");
        _size = n;
    }
};

// ============================================================================
// Convenience: push a pod value as bytes into the constants vector
// ============================================================================

template<typename T>
inline std::vector<uint8_t> as_constant(const T& val) {
    std::vector<uint8_t> bytes(sizeof(T));
    std::memcpy(bytes.data(), &val, sizeof(T));
    return bytes;
}
