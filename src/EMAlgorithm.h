#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "KmerIndex.h"
#include "weights.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>


// smallest weight we expect is ~10^-4
// on most machines, TOLERANCE should be 2.22045e-15
const double TOLERANCE = std::numeric_limits<double>::epsilon() * 10;

struct EMAlgorithm {
  // ecmap is the ecmap from KmerIndex
  // counts is vector from collector, with indices corresponding to ec ids
  // target_names is the target_names_ from collector
  // TODO: initialize alpha a bit more intelligently
  EMAlgorithm(const EcMap& ecmap,
              const std::vector<int>& counts,
              const std::vector<std::string>& target_names,
              const std::vector<double>& eff_lens,
              const WeightMap& wm) :
    //idx_(idx),
    num_trans_(target_names.size()),
    ecmap_(ecmap),
    counts_(counts),
    target_names_(target_names),
    eff_lens_(eff_lens),
    weight_map_(wm),
    alpha_(num_trans_, 1.0/num_trans_), // uniform distribution over transcripts
    rho_(num_trans_, 0.0),
    rho_set_(false)
  {
      assert(target_names_.size() == eff_lens.size());
  }

  ~EMAlgorithm() {}

  void run(size_t n_iter = 500) {
    std::vector<double> next_alpha(alpha_.size(), 0.0);

    assert(weight_map_.size() <= counts_.size());

    double denom;

    std::cout << "[em]\tfishing for the right mixture (. = 50 rounds)" <<
              std::endl;

    for (auto i = 0; i < n_iter; ++i) {
      if (i % 50 == 0) {
        std::cout << ".";
        std::cout.flush();
        if (i % 500 == 0 && i > 0) {
          std::cout << std::endl;
        }
      }

      for (auto& ec_kv : ecmap_ ) {
        denom = 0.0;

        // first, compute the denominator: a normalizer
        // iterate over transcripts in EC map
        auto w_search = weight_map_.find(ec_kv.first);

        // everything in ecmap should be in weight_map
        assert( w_search != weight_map_.end() );
        assert( w_search->second.size() == ec_kv.second.size() );

        for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
          denom += alpha_[ec_kv.second[t_it]] * w_search->second[t_it];
        }

        if (denom < TOLERANCE) {
          continue;
        }

        // compute the update step
        for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
          next_alpha[ec_kv.second[t_it]] += counts_[ec_kv.first] *
                                            ((w_search->second[t_it] * alpha_[ec_kv.second[t_it]]) / denom);
        }
      }

      // TODO: check for relative difference for convergence in EM

      // reassign alpha_ to next_alpha
      std::copy(next_alpha.begin(), next_alpha.end(), alpha_.begin());

      // clear all next_alpha values 0 for next iteration
      std::fill(next_alpha.begin(), next_alpha.end(), 0.0);
    }

    std::cout << std::endl;
    std::cout.flush();
  }

  void compute_rho() {
    if (rho_set_) {
      // rho has already been set, let's clear it
      std::fill(rho_.begin(), rho_.end(), 0.0);
    }

    double total {0.0};
    for (auto i = 0; i < alpha_.size(); ++i) {
      if (eff_lens_[i] < TOLERANCE) {
        std::cerr << "Should actually never really get here... tid: "  << i <<
            std::endl;
        continue;
      }
      rho_[i] = alpha_[i] / eff_lens_[i];
      total += rho_[i];
    }

    for (auto& r : rho_) {
      r /= total;
    }

    rho_set_ = true;
  }

  void write(const std::string& out_fname) const {
    std::ofstream out;
    out.open(out_fname, std::ios::out);

    if (!out.is_open()) {
      std::cerr << "Error opening '" << out_fname << "'" <<
          std::endl;
      exit(1);
    }

    out.precision(15);

    out <<
        "target_id" << "\t" <<
        "kallisto_id" << "\t" <<
        "rho" << "\t" <<
        "tpm" << "\t" <<
        "est_counts" <<
        std::endl;

    const double MILLION = 1e6;

    for (auto i = 0; i < rho_.size(); ++i) {
      out <<
          target_names_[i] << "\t" <<
          i << "\t" <<
          rho_[i] << "\t" <<
          rho_[i] * MILLION << "\t" <<
          alpha_[i] <<
          std::endl;
    }

    out.flush();
    out.close();
  }

  int num_trans_;
  const EcMap& ecmap_;
  const std::vector<int>& counts_;
  const std::vector<std::string>& target_names_;
  const std::vector<double>& eff_lens_;
  const WeightMap& weight_map_;
  std::vector<double> alpha_;
  std::vector<double> rho_;
  bool rho_set_;
};

#endif // KALLISTO_EMALGORITHM_H
