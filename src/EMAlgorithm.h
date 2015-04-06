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
//const double TOLERANCE = std::numeric_limits<double>::epsilon() * 10;
//const double TOLERANCE = 1e-100;
const double TOLERANCE = std::numeric_limits<double>::denorm_min();

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

  void run(size_t n_iter = 1000, size_t min_rounds=50) {
    std::vector<double> next_alpha(alpha_.size(), 0.0);

    assert(weight_map_.size() <= counts_.size());

    double denom;
    const double alpha_limit = 1e-7;
    const double alpha_change = 1e-3;
    
    std::cerr << "[em]\tfishing for the right mixture (. = 50 rounds)" <<
              std::endl;
    int i;
    for (i = 0; i < n_iter; ++i) {
      if (i % 50 == 0) {
        std::cerr << ".";
        std::cerr.flush();
        if (i % 500 == 0 && i > 0) {
          std::cerr << std::endl;
        }
      }
      
      //for (auto& ec_kv : ecmap_ ) {
      for (int ec = 0; ec < num_trans_; ec++) {
        next_alpha[ec] = counts_[ec];
      }
      
      for (int ec = num_trans_; ec < ecmap_.size();  ec++) {
        denom = 0.0;

        // first, compute the denominator: a normalizer
        // iterate over transcripts in EC map
        auto& wv = weight_map_[ec];

        // everything in ecmap should be in weight_map
        //assert( w_search != weight_map_.end() );
        //assert( w_search->second.size() == ec_kv.second.size() );

        // wv is weights vector
        // v is ec vector


        auto& v = ecmap_[ec]; //ecmap_.find(ec)->second;
        auto numEC = v.size();
          
        for (auto t_it = 0; t_it < numEC; ++t_it) {
          denom += alpha_[v[t_it]] * wv[t_it];
        }

        if (denom < TOLERANCE) {
          continue;
        }

        // compute the update step
        auto countNorm = counts_[ec] / denom;
        for (auto t_it = 0; t_it < numEC; ++t_it) {
          next_alpha[v[t_it]] +=  (wv[t_it] * alpha_[v[t_it]]) * countNorm;
        }

      }

      // TODO: check for relative difference for convergence in EM

      bool stopEM = (i >= min_rounds); // false initially
      //double maxChange = 0.0;

      for (int ec = 0; ec < num_trans_; ec++) {
          
        if (stopEM && next_alpha[ec] >= alpha_limit && abs(next_alpha[ec]-alpha_[ec]) / next_alpha[ec] >= alpha_change) {
          stopEM = false;
        }

        /*
        if (next_alpha[ec] > alpha_limit) {
          maxChange = std::max(maxChange,abs(next_alpha[ec]-alpha_[ec]) / next_alpha[ec]);
        }
        */

        
        // reassign alpha_ to next_alpha
        alpha_[ec] = next_alpha[ec];
        
        // clear all next_alpha values 0 for next iteration
        next_alpha[ec] = 0;
        
      }
      
      // std::cout << maxChange << std::endl;
      if (stopEM) {
        break;
      }
      


      //std::copy(next_alpha.begin(), next_alpha.end(), alpha_.begin());


      //std::fill(next_alpha.begin(), next_alpha.end(), 0.0);
    }


    std::cerr << std::endl << "Ran for " << i << " rounds of EM";
    
    std::cerr << std::endl;
    std::cerr.flush();
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
