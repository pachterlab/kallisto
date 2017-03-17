#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "KmerIndex.h"
#include "MinCollector.h"
#include "weights.h"

#include <algorithm>
#include <numeric>
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
  EMAlgorithm(const std::vector<int>& counts,
              const KmerIndex& index,
              const MinCollector& tc,
              const std::vector<double>& all_means,
              const ProgramOptions& opt) :
    index_(index),
    tc_(tc),
    num_trans_(index.target_names_.size()),
    ecmap_(index.ecmap),
    counts_(counts),
    target_names_(index.target_names_),
    post_bias_(4096,1.0),
    alpha_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    rho_(num_trans_, 0.0),
    rho_set_(false),
    all_fl_means(all_means),
    opt(opt)
  {
    assert(all_fl_means.size() == index_.target_lens_.size());
    eff_lens_ = calc_eff_lens(index_.target_lens_, all_fl_means);
    weight_map_ = calc_weights (tc_.counts, ecmap_, eff_lens_);    
    assert(target_names_.size() == eff_lens_.size());
  }

  ~EMAlgorithm() {}

  void run(size_t n_iter = 10000, size_t min_rounds=50, bool verbose = true, bool recomputeEffLen = true) {
    std::vector<double> next_alpha(alpha_.size(), 0.0);

    assert(weight_map_.size() <= counts_.size());

    double denom;
    const double alpha_limit = 1e-7;
    const double alpha_change_limit = 1e-2;
    const double alpha_change = 1e-2;
    bool finalRound = false;

    if (verbose) {
      std::cerr << "[   em] quantifying the abundances ..."; std::cerr.flush();
    }

    int i;
    for (i = 0; i < n_iter; ++i) {
      if (recomputeEffLen && (i == min_rounds || i == min_rounds + 500)) {
        eff_lens_ = update_eff_lens(all_fl_means, tc_, index_, alpha_, eff_lens_, post_bias_, opt);
        weight_map_ = calc_weights (tc_.counts, ecmap_, eff_lens_);
      }


      //for (auto& ec_kv : ecmap_ ) {
      for (int ec = 0; ec < num_trans_; ec++) {
        next_alpha[ec] = counts_[ec];
      }


      for (int ec = num_trans_; ec < ecmap_.size();  ec++) {
        denom = 0.0;

        if (counts_[ec] == 0) {
          continue;
        }

        // first, compute the denominator: a normalizer
        // iterate over targets in EC map
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

      bool stopEM = false; //!finalRound && (i >= min_rounds); // false initially
      //double maxChange = 0.0;
      int chcount = 0;
      for (int ec = 0; ec < num_trans_; ec++) {
        if (next_alpha[ec] > alpha_change_limit && (std::fabs(next_alpha[ec] - alpha_[ec]) / next_alpha[ec]) > alpha_change) {
          chcount++;
        }

        //if (stopEM && next_alpha[ec] >= alpha_limit) {

          /* double reldiff = abs(next_alpha[ec]-alpha_[ec]) / next_alpha[ec];
          if (reldiff >= alpha_change) {
            stopEM = false;
            }*/
        //}

        /*
        if (next_alpha[ec] > alpha_limit) {
          maxChange = std::max(maxChange,std::fabs(next_alpha[ec]-alpha_[ec]) / next_alpha[ec]);
        }
        */
        // reassign alpha_ to next_alpha
        alpha_[ec] = next_alpha[ec];

        // clear all next_alpha values 0 for next iteration
        next_alpha[ec] = 0.0;
      }

      //std::cout << chcount << std::endl;
      if (chcount == 0 && i > min_rounds) {

        stopEM=true;
      }

      if (finalRound) {
        break;
      }

      // std::cout << maxChange << std::endl;
      if (stopEM) {
        finalRound = true;
        alpha_before_zeroes_.resize( alpha_.size() );
        for (int ec = 0; ec < num_trans_; ec++) {
          alpha_before_zeroes_[ec] = alpha_[ec];
          if (alpha_[ec] < alpha_limit/10.0) {
            alpha_[ec] = 0.0;
          }
        }
      }

    }

    // ran for the maximum number of iterations
    if (n_iter == i) {
      alpha_before_zeroes_.resize( alpha_.size() );
      for (int ec = 0; ec < num_trans_; ec++) {
        alpha_before_zeroes_[ec] = alpha_[ec];
      }
    }

    if (verbose) {
      std::cerr << " done" << std::endl;
      std::cerr << "[   em] the Expectation-Maximization algorithm ran for "
        << pretty_num(i) << " rounds";
      std::cerr << std::endl;
      std::cerr.flush();
    }

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

  // DEPRECATED:
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

  void set_start(const EMAlgorithm& em_start) {
    assert(em_start.alpha_before_zeroes_.size() == alpha_.size());
    double big = 1.0;
    double sum_counts = std::accumulate(counts_.begin(), counts_.end(), 0.0);
    double sum_big = 0.0;
    int count_big = 0;
    for (auto x : em_start.alpha_before_zeroes_) {
      if (x >= big) {
        sum_big += x;
        count_big++;
      }
    }
    int n = alpha_.size();
    for (auto i = 0; i < n; i++) {
      if (em_start.alpha_before_zeroes_[i] >= big) {
        alpha_[i] = em_start.alpha_before_zeroes_[i];
      } else {
        alpha_[i] = sum_counts/(n - count_big);
      }
    }

    //std::cout << sum_big << " " << count_big << " " << n << std::endl;

    std::copy(em_start.alpha_before_zeroes_.begin(), em_start.alpha_before_zeroes_.end(),
        alpha_.begin());
  }


  int num_trans_;
  const KmerIndex& index_;
  const MinCollector& tc_;
  const EcMap& ecmap_;
  const std::vector<int>& counts_;
  const std::vector<std::string>& target_names_;
  const std::vector<double>& all_fl_means;
  std::vector<double> eff_lens_;
  std::vector<double> post_bias_;
  WeightMap weight_map_;
  std::vector<double> alpha_;
  std::vector<double> alpha_before_zeroes_;
  std::vector<double> rho_;
  bool rho_set_;
  const ProgramOptions& opt;
};


#endif // KALLISTO_EMALGORITHM_H
