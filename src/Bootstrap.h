#ifndef KALLISTO_BOOTSTRAP_H
#define KALLISTO_BOOTSTRAP_H

#include "KmerIndex.h"
#include "MinCollector.h"
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
            const KmerIndex& index,
            const MinCollector& tc,
            const std::vector<double>& eff_lens,
            double mean,
            size_t seed) :
    index_(index),
    tc_(tc),
    eff_lens_(eff_lens),
    mean_fl(mean),
    seed_(seed),
    mult_(true_counts, seed_)
    {}
  
  // EM Algorithm generates a sample from the Multinomial, then returns
  // an "EMAlgorithm" that has already run the EM as well as compute the
        // rho values
  EMAlgorithm run_em(const EMAlgorithm& em_start);
  
private:
  const KmerIndex& index_;
  const MinCollector& tc_;
  const std::vector<double>& eff_lens_;
  double mean_fl;
  size_t seed_;
  Multinomial mult_;
};

#endif
