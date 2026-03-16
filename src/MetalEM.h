#pragma once

#include "MetalUtils.h"
#include "MetalIndex.h"
#include "BenchmarkStats.h"

#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

// Metal GPU shaders only support float (not double); use float for GPU EM buffers.
// Convergence check on CPU casts back to double for the tolerance test.
static const float METAL_EM_TOLERANCE = std::numeric_limits<float>::min();

// ============================================================================
// MetalEM — gather-variant EM algorithm
//
// compute_denom_kernel: one thread per EC → denom[e]
// em_gather_kernel:    one thread per transcript → alpha_next[t] (no atomics)
//
// Convergence check is CPU-side reading shared buffers directly.
// ============================================================================

struct MetalEM {
    MetalBuffer<float> d_alpha;
    MetalBuffer<float> d_alpha_next;
    MetalBuffer<float> d_eff_lens;
    MetalBuffer<float> d_denom;
    MetalBuffer<int>    d_trans_ecs;
    MetalBuffer<int>    d_trans_ec_offsets;

    size_t num_trans = 0;
    size_t num_ecs   = 0;

    MetalEM(size_t num_trans_,
            size_t num_ecs_,
            const std::vector<double>& eff_lens);

    void build_transpose(const MetalECMap& ecmap);

    int run_transpose(const MetalECMap& ecmap,
                      const MetalBuffer<int>& ec_counts,
                      int max_iter = 10000,
                      int min_rounds = 50);
};
