#pragma once

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <chrono>

// CUDA error checking macro
#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = call;                                                \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,   \
                    cudaGetErrorString(err));                                  \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0)

// Timer class for benchmarking
class Timer {
public:
    Timer() : start_time_(std::chrono::high_resolution_clock::now()) {}
    
    void reset() {
        start_time_ = std::chrono::high_resolution_clock::now();
    }
    
    double elapsed_ms() const {
        auto end_time = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end_time - start_time_).count();
    }
    
    double elapsed_s() const {
        return elapsed_ms() / 1000.0;
    }
    
    // Calculate throughput in GB/s
    static double throughput_gbps(size_t bytes, double seconds) {
        return (bytes / (1024.0 * 1024.0 * 1024.0)) / seconds;
    }
    
private:
    std::chrono::high_resolution_clock::time_point start_time_;
};

// GPU Timer using CUDA events
class GpuTimer {
public:
    GpuTimer() {
        CUDA_CHECK(cudaEventCreate(&start_));
        CUDA_CHECK(cudaEventCreate(&stop_));
    }
    
    ~GpuTimer() {
        cudaEventDestroy(start_);
        cudaEventDestroy(stop_);
    }
    
    void start(cudaStream_t stream = 0) {
        CUDA_CHECK(cudaEventRecord(start_, stream));
    }
    
    void stop(cudaStream_t stream = 0) {
        CUDA_CHECK(cudaEventRecord(stop_, stream));
    }
    
    float elapsed_ms() {
        CUDA_CHECK(cudaEventSynchronize(stop_));
        float ms = 0;
        CUDA_CHECK(cudaEventElapsedTime(&ms, start_, stop_));
        return ms;
    }
    
private:
    cudaEvent_t start_, stop_;
};

// Sequence descriptor: offset and length in the buffer
struct SeqDescriptor {
    uint32_t offset;
    uint32_t length;
};

// Print memory usage
inline void print_gpu_memory() {
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    printf("GPU Memory: %.2f GB free / %.2f GB total\n",
           free_mem / (1024.0 * 1024.0 * 1024.0),
           total_mem / (1024.0 * 1024.0 * 1024.0));
}

// Default batch size (256 MB)
constexpr size_t DEFAULT_BATCH_SIZE = 256 * 1024 * 1024;


