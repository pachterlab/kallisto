#ifndef GPU_EM_CUH
#define GPU_EM_CUH

#include "GPUIndex.cuh"
#include "BenchmarkStats.h"
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/copy.h>
#include <cuda_runtime.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

static const double EM_TOLERANCE = std::numeric_limits<double>::denorm_min();

// --- Original approach: one thread per EC, atomicAdd scatter ---

__global__ void em_iteration_kernel(
    const double* __restrict__ alpha,
    const double* __restrict__ eff_lens,
    const int* __restrict__ transcripts,
    const uint64_t* __restrict__ offsets,
    const int* __restrict__ ec_counts,
    double* __restrict__ alpha_next,
    size_t num_ecs,
    double tolerance)
{
  size_t e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= num_ecs) return;

  int count = ec_counts[e];
  if (count == 0) return;

  uint64_t start = offsets[e];
  uint64_t end   = offsets[e + 1];
  if (start >= end) return;

  if (end - start == 1) {
    int t = transcripts[start];
    atomicAdd(&alpha_next[t], static_cast<double>(count));
    return;
  }

  double denom = 0.0;
  for (uint64_t i = start; i < end; ++i) {
    int t = transcripts[i];
    denom += alpha[t] / eff_lens[t];
  }

  if (denom < tolerance) return;

  double count_dbl = static_cast<double>(count);
  for (uint64_t i = start; i < end; ++i) {
    int t = transcripts[i];
    atomicAdd(&alpha_next[t], count_dbl * alpha[t] / eff_lens[t] / denom);
  }
}

// --- Transpose approach: denom per EC, then gather per transcript ---

__global__ void compute_denom_kernel(
    const double* __restrict__ alpha,
    const double* __restrict__ eff_lens,
    const int* __restrict__ transcripts,
    const uint64_t* __restrict__ offsets,
    const int* __restrict__ ec_counts,
    double* __restrict__ denom,
    size_t num_ecs)
{
  size_t e = blockIdx.x * blockDim.x + threadIdx.x;
  if (e >= num_ecs) return;

  if (ec_counts[e] == 0) {
    denom[e] = 0.0;
    return;
  }

  uint64_t start = offsets[e];
  uint64_t end   = offsets[e + 1];

  double d = 0.0;
  for (uint64_t i = start; i < end; ++i) {
    int t = transcripts[i];
    d += alpha[t] / eff_lens[t];
  }
  denom[e] = d;
}

// One thread per transcript: gather contributions from all ECs containing it.
// No atomicAdd -- each thread writes its own alpha_next[t] directly.
__global__ void em_gather_kernel(
    const double* __restrict__ alpha,
    const double* __restrict__ eff_lens,
    const int* __restrict__ trans_ecs,
    const int* __restrict__ trans_ec_offsets,
    const int* __restrict__ ec_counts,
    const double* __restrict__ denom,
    double* __restrict__ alpha_next,
    size_t num_trans,
    double tolerance)
{
  size_t t = blockIdx.x * blockDim.x + threadIdx.x;
  if (t >= num_trans) return;

  int start = trans_ec_offsets[t];
  int end   = trans_ec_offsets[t + 1];

  if (start >= end) {
    alpha_next[t] = 0.0;
    return;
  }

  double a_over_el = alpha[t] / eff_lens[t];
  double sum = 0.0;

  for (int i = start; i < end; ++i) {
    int e = trans_ecs[i];
    int count = ec_counts[e];
    if (count == 0) continue;
    double d = denom[e];
    if (d < tolerance) continue;
    sum += static_cast<double>(count) * a_over_el / d;
  }

  alpha_next[t] = sum;
}


struct GPUEM {
  thrust::device_vector<double> d_alpha;
  thrust::device_vector<double> d_alpha_next;
  thrust::device_vector<double> d_eff_lens;
  size_t num_trans;
  size_t num_ecs;

  // Transpose CSR: for each transcript, list of ECs containing it
  thrust::device_vector<int> d_trans_ecs;
  thrust::device_vector<int> d_trans_ec_offsets;
  thrust::device_vector<double> d_denom;

  GPUEM(size_t num_trans_, size_t num_ecs_,
        const std::vector<double>& eff_lens)
    : num_trans(num_trans_), num_ecs(num_ecs_),
      d_alpha(num_trans_, 1.0 / static_cast<double>(num_trans_)),
      d_alpha_next(num_trans_, 0.0),
      d_eff_lens(eff_lens.begin(), eff_lens.end())
  {}

  void build_transpose(const GPUECMap& ecmap) {
    std::vector<int> h_transcripts(ecmap.transcripts.size());
    thrust::copy(ecmap.transcripts.begin(), ecmap.transcripts.end(), h_transcripts.begin());
    std::vector<uint64_t> h_offsets(ecmap.offsets.size());
    thrust::copy(ecmap.offsets.begin(), ecmap.offsets.end(), h_offsets.begin());

    std::vector<std::vector<int>> t2e(num_trans);
    for (size_t e = 0; e < num_ecs; ++e) {
      uint64_t start = h_offsets[e];
      uint64_t end   = h_offsets[e + 1];
      for (uint64_t i = start; i < end; ++i) {
        int t = h_transcripts[i];
        if (t >= 0 && t < (int)num_trans) {
          t2e[t].push_back(static_cast<int>(e));
        }
      }
    }

    std::vector<int> flat_ecs;
    std::vector<int> flat_offsets(num_trans + 1, 0);
    for (size_t t = 0; t < num_trans; ++t) {
      flat_offsets[t + 1] = flat_offsets[t] + static_cast<int>(t2e[t].size());
      for (int e : t2e[t]) {
        flat_ecs.push_back(e);
      }
    }

    d_trans_ecs.assign(flat_ecs.begin(), flat_ecs.end());
    d_trans_ec_offsets.assign(flat_offsets.begin(), flat_offsets.end());
    d_denom.resize(num_ecs);

    std::cerr << "[   em] built transpose map: " << flat_ecs.size()
              << " (transcript,EC) pairs" << std::endl;
  }

  // Original scatter-based EM (one thread per EC, atomicAdd)
  int run(const GPUECMap& ecmap,
          const thrust::device_vector<int>& ec_counts,
          int max_iter = 10000,
          int min_rounds = 50)
  {
    const double alpha_change_limit = 1e-2;
    const double alpha_change = 1e-2;
    const double alpha_limit = 1e-7;
    const int block_size = 256;

    const int* d_transcripts = ecmap.transcripts.data().get();
    const uint64_t* d_offsets = ecmap.offsets.data().get();
    const int* d_ec_counts = ec_counts.data().get();
    const double* d_el = d_eff_lens.data().get();

    size_t actual_num_ecs = ec_counts.size();
    if (actual_num_ecs > num_ecs) {
      num_ecs = actual_num_ecs;
    }

    int grid_ecs = (num_ecs + block_size - 1) / block_size;

    std::cerr << "[   em] quantifying the abundances (scatter) ..."; std::cerr.flush();

    const int BATCH_SIZE = 10;
    bool finalRound = false;
    int i = 0;

    while (i < max_iter) {
      int batch_end = std::min(i + BATCH_SIZE, max_iter);

      for (int j = i; j < batch_end; ++j) {
        thrust::fill(d_alpha_next.begin(), d_alpha_next.end(), 0.0);
        em_iteration_kernel<<<grid_ecs, block_size>>>(
          d_alpha.data().get(), d_el, d_transcripts, d_offsets,
          d_ec_counts, d_alpha_next.data().get(), num_ecs, EM_TOLERANCE);
        d_alpha.swap(d_alpha_next);
      }

      i = batch_end;
      cudaDeviceSynchronize();

      int chcount = 0;
      {
        const double* a_ptr = d_alpha_next.data().get();
        const double* an_ptr = d_alpha.data().get();
        size_t nt = num_trans;
        double acl = alpha_change_limit;
        double ac = alpha_change;
        chcount = thrust::count_if(
          thrust::device,
          thrust::make_counting_iterator<size_t>(0),
          thrust::make_counting_iterator<size_t>(nt),
          [a_ptr, an_ptr, acl, ac] __device__ (size_t t) {
            double next = an_ptr[t];
            if (next <= acl) return false;
            double prev = a_ptr[t];
            double rel = (next > prev) ? (next - prev) : (prev - next);
            return (rel / next) > ac;
          }
        );
      }

      bool stopEM = (chcount == 0 && i > min_rounds);

      if (finalRound) break;

      if (stopEM) {
        finalRound = true;
        double lim = alpha_limit / 10.0;
        thrust::transform(d_alpha.begin(), d_alpha.end(),
                          d_alpha.begin(),
                          [lim] __device__ (double x) { return (x < lim) ? 0.0 : x; });
      }
    }

    std::cerr << " done" << std::endl;
    std::cerr << "[   em] the Expectation-Maximization algorithm ran for "
              << i << " rounds" << std::endl;
    return i;
  }

  // Transpose gather-based EM (denom per EC, then one thread per transcript gathers)
  int run_transpose(const GPUECMap& ecmap,
                    const thrust::device_vector<int>& ec_counts,
                    int max_iter = 10000,
                    int min_rounds = 50)
  {
    const double alpha_change_limit = 1e-2;
    const double alpha_change = 1e-2;
    const double alpha_limit = 1e-7;
    const int block_size = 256;

    const int* d_transcripts = ecmap.transcripts.data().get();
    const uint64_t* d_offsets = ecmap.offsets.data().get();
    const int* d_ec_counts = ec_counts.data().get();
    const double* d_el = d_eff_lens.data().get();

    size_t actual_num_ecs = ec_counts.size();
    if (actual_num_ecs > num_ecs) {
      num_ecs = actual_num_ecs;
      d_denom.resize(num_ecs);
    }

    int grid_ecs = (num_ecs + block_size - 1) / block_size;
    int grid_trans = (num_trans + block_size - 1) / block_size;

    std::cerr << "[   em] quantifying the abundances (gather) ..."; std::cerr.flush();

    const int BATCH_SIZE = 10;
    bool finalRound = false;
    int i = 0;

    while (i < max_iter) {
      int batch_end = std::min(i + BATCH_SIZE, max_iter);

      for (int j = i; j < batch_end; ++j) {
        compute_denom_kernel<<<grid_ecs, block_size>>>(
          d_alpha.data().get(), d_el, d_transcripts, d_offsets,
          d_ec_counts, d_denom.data().get(), num_ecs);

        em_gather_kernel<<<grid_trans, block_size>>>(
          d_alpha.data().get(), d_el,
          d_trans_ecs.data().get(), d_trans_ec_offsets.data().get(),
          d_ec_counts, d_denom.data().get(),
          d_alpha_next.data().get(), num_trans, EM_TOLERANCE);

        d_alpha.swap(d_alpha_next);
      }

      i = batch_end;
      cudaDeviceSynchronize();

      int chcount = 0;
      {
        const double* a_ptr = d_alpha_next.data().get();
        const double* an_ptr = d_alpha.data().get();
        size_t nt = num_trans;
        double acl = alpha_change_limit;
        double ac = alpha_change;
        chcount = thrust::count_if(
          thrust::device,
          thrust::make_counting_iterator<size_t>(0),
          thrust::make_counting_iterator<size_t>(nt),
          [a_ptr, an_ptr, acl, ac] __device__ (size_t t) {
            double next = an_ptr[t];
            if (next <= acl) return false;
            double prev = a_ptr[t];
            double rel = (next > prev) ? (next - prev) : (prev - next);
            return (rel / next) > ac;
          }
        );
      }

      bool stopEM = (chcount == 0 && i > min_rounds);

      if (finalRound) break;

      if (stopEM) {
        finalRound = true;
        double lim = alpha_limit / 10.0;
        thrust::transform(d_alpha.begin(), d_alpha.end(),
                          d_alpha.begin(),
                          [lim] __device__ (double x) { return (x < lim) ? 0.0 : x; });
      }
    }

    std::cerr << " done" << std::endl;
    std::cerr << "[   em] the Expectation-Maximization algorithm ran for "
              << i << " rounds" << std::endl;
    return i;
  }
};

#endif // GPU_EM_CUH
