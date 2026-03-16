#import "MetalEM.h"

#include <algorithm>
#include <cstring>

MetalEM::MetalEM(size_t num_trans_,
                 size_t num_ecs_,
                 const std::vector<double>& eff_lens)
    : num_trans(num_trans_), num_ecs(num_ecs_)
{
    float init_alpha = 1.0f / static_cast<float>(num_trans_);

    d_alpha.resize(num_trans_);
    d_alpha_next.resize(num_trans_);
    d_eff_lens.resize(num_trans_);

    float* a  = d_alpha.data();
    float* an = d_alpha_next.data();
    float* el = d_eff_lens.data();

    for (size_t i = 0; i < num_trans_; ++i) {
        a[i]  = init_alpha;
        an[i] = 0.0f;
        el[i] = static_cast<float>(eff_lens[i]);
    }
}

void MetalEM::build_transpose(const MetalECMap& ecmap) {
    // Build transpose CSR on CPU (ecmap buffers are shared memory, readable directly)
    const int*      tx  = ecmap.transcripts.data();
    const uint64_t* off = ecmap.offsets.data();
    size_t n_ecs = ecmap.num_ecs;

    std::vector<std::vector<int>> t2e(num_trans);
    for (size_t e = 0; e < n_ecs; ++e) {
        uint64_t start = off[e];
        uint64_t end   = off[e + 1];
        for (uint64_t i = start; i < end; ++i) {
            int t = tx[i];
            if (t >= 0 && (size_t)t < num_trans)
                t2e[t].push_back((int)e);
        }
    }

    std::vector<int> flat_ecs;
    std::vector<int> flat_offsets(num_trans + 1, 0);
    for (size_t t = 0; t < num_trans; ++t) {
        flat_offsets[t + 1] = flat_offsets[t] + (int)t2e[t].size();
        for (int e : t2e[t]) flat_ecs.push_back(e);
    }

    d_trans_ecs.from_host(flat_ecs);
    d_trans_ec_offsets.from_host(flat_offsets);
    d_denom.resize(num_ecs);

    std::cerr << "[   em] built transpose map: " << flat_ecs.size()
              << " (transcript,EC) pairs" << std::endl;
}

int MetalEM::run_transpose(const MetalECMap& ecmap,
                            const MetalBuffer<int>& ec_counts,
                            int max_iter,
                            int min_rounds)
{
    const double alpha_change_limit = 1e-2;
    const double alpha_change       = 1e-2;
    const double alpha_limit        = 1e-7;

    size_t actual_num_ecs = ec_counts.size();
    if (actual_num_ecs > num_ecs) {
        num_ecs = actual_num_ecs;
        d_denom.resize(num_ecs);
    }

    std::cerr << "[   em] quantifying the abundances (gather) ..."; std::cerr.flush();

    const int BATCH_SIZE = 10;
    bool finalRound = false;
    int i = 0;

    float tolerance = METAL_EM_TOLERANCE;

    while (i < max_iter) {
        int batch_end = std::min(i + BATCH_SIZE, max_iter);

        for (int j = i; j < batch_end; ++j) {
            // Zero alpha_next
            float* an = d_alpha_next.data();
            for (size_t t = 0; t < num_trans; ++t) an[t] = 0.0f;

            // compute_denom_kernel
            MetalContext::get().dispatch(
                "compute_denom_kernel",
                (NSUInteger)num_ecs,
                {
                    d_alpha.metalBuffer(),
                    d_eff_lens.metalBuffer(),
                    ecmap.transcripts.metalBuffer(),
                    ecmap.offsets.metalBuffer(),
                    ec_counts.metalBuffer(),
                    d_denom.metalBuffer()
                },
                { as_constant((uint64_t)num_ecs) }
            );

            // em_gather_kernel
            MetalContext::get().dispatch(
                "em_gather_kernel",
                (NSUInteger)num_trans,
                {
                    d_alpha.metalBuffer(),
                    d_eff_lens.metalBuffer(),
                    d_trans_ecs.metalBuffer(),
                    d_trans_ec_offsets.metalBuffer(),
                    ec_counts.metalBuffer(),
                    d_denom.metalBuffer(),
                    d_alpha_next.metalBuffer()
                },
                {
                    as_constant((uint64_t)num_trans),
                    as_constant(tolerance)
                }
            );

            // Swap alpha <-> alpha_next (both are CPU-visible shared buffers)
            std::swap(d_alpha, d_alpha_next);
        }

        i = batch_end;

        // Convergence check (CPU reads shared buffers directly)
        const float* a  = d_alpha_next.data();  // previous (was alpha before swap)
        const float* an = d_alpha.data();        // current
        int chcount = 0;
        for (size_t t = 0; t < num_trans; ++t) {
            float next = an[t];
            if (next <= (float)alpha_change_limit) continue;
            float prev = a[t];
            float rel  = std::fabs(next - prev);
            if (rel / next > (float)alpha_change) ++chcount;
        }

        bool stopEM = (chcount == 0 && i > min_rounds);
        if (finalRound) break;
        if (stopEM) {
            finalRound = true;
            float lim = (float)(alpha_limit / 10.0);
            float* alpha_ptr = d_alpha.data();
            for (size_t t = 0; t < num_trans; ++t)
                if (alpha_ptr[t] < lim) alpha_ptr[t] = 0.0f;
        }
    }

    std::cerr << " done" << std::endl;
    std::cerr << "[   em] the Expectation-Maximization algorithm ran for "
              << i << " rounds" << std::endl;
    return i;
}
