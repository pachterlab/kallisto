#ifndef KALLISTO_BOOTSTRAP_H
#define KALLISTO_BOOTSTRAP_H

#include "KmerIndex.h"
#include "weights.h"
#include "EMAlgorithm.h"
#include "Multinomial.hpp"

class Bootstrap {
    // needs:
    // - "true" counts
    // - ecmap
    // - target_names
    // - eff_lens
    public:
        Bootstrap(const std::vector<int>& true_counts,
                const EcMap& ecmap,
                const std::vector<std::string>& target_names,
                const std::vector<double>& eff_lens,
                size_t seed) :
            ecmap_(ecmap),
            target_names_(target_names),
            eff_lens_(eff_lens),
            seed_(seed),
            mult_(true_counts, seed_)
        {}

        // EM Algorithm generates a sample from the Multinomial, then returns
        // an "EMAlgorithm" that has already run the EM as well as compute the
        // rho values
        EMAlgorithm run_em(const EMAlgorithm& em_start);

    private:
        const EcMap& ecmap_;
        const std::vector<std::string>& target_names_;
        const std::vector<double>& eff_lens_;
        size_t seed_;
        Multinomial mult_;
};

#endif
