#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "weights.h"

#include <algorithm>
#include <iostream>
#include <vector>


const double TOLERANCE = 1e-5;;

template <typename Index>
struct EMAlgorithm {

    // counts is vector from collector, with indices corresponding to ec ids
    // TODO: initialize alpha a bit more intelligently
    // TODO: refactor to remove dependence on Index
	EMAlgorithm(const ProgramOptions& opt, const Index& idx,
	        const std::vector<int>& counts,
	        const std::vector<double>& eff_lens,
	        const WeightMap& wm) :
	    idx_(idx),
	    counts_(counts),
	    num_trans_(idx.num_trans),
	    eff_lens_(eff_lens),
	    weight_map_(wm),
	    alpha_(idx.num_trans, 1.0/idx.num_trans), // uniform distribution over transcripts
	    rho_(idx.num_trans, 0.0),
	    rho_set_(false)
	{}

	~EMAlgorithm() {}

	void run(size_t n_iter = 500) {
        std::vector<double> next_alpha(alpha_.size(), 0.0);

        assert(weight_map_.size() == counts_.size());

        double denom;
	    for (auto i = 0; i < n_iter; ++i) {

            std::cout << "it: " << i << std::endl;
            for (auto& ec_kv : idx_.ecmap ) {

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

                /* std::cout << std::endl; */

                if (denom < TOLERANCE) {
                    continue;
                }

                /* std::cout << "denom: " << denom << std::endl; */

                // compute the update step
                for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
                    next_alpha[ec_kv.second[t_it]] += counts_[ec_kv.first] *
                        (w_search->second[t_it] * alpha_[ec_kv.second[t_it]] / denom);
                }
            }

            // TODO: check for relative difference for convergence in EM

            // reassign alpha_ to next_alpha
            std::copy(next_alpha.begin(), next_alpha.end(), alpha_.begin());

            // clear all next_alpha values 0 for next iteration
            std::fill(next_alpha.begin(), next_alpha.end(), 0.0);
	    }
	}

    void compute_rho() {

        if (rho_set_) {
            // rho has already been set, let's clear it
            std::fill(rho_.begin(), rho_.end(), 0.0);
        }

        for (auto& ec_kv : idx_.ecmap ) {
            double denom {0.0};

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

            for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
                rho_[ec_kv.second[t_it]] += alpha_[ec_kv.second[t_it]] *
                    w_search->second[t_it] / denom;
            }
        }

        rho_set_ = true;
    }

    void write(const std::string& rho_out) const {
        std::ofstream out;
        out.open(rho_out, std::ios::out);

        out.precision(15);
        std::cout.precision(15);

        std::cout << "alphas" << std::endl;
        for (auto i = 0; i < rho_.size(); ++i) {
            out << i << "\t" << std::fixed << rho_[i] << std::endl;
            std::cout << std::fixed << alpha_[i] << std::endl;
        }

        out.close();
    }

	int num_trans_;
	const Index &idx_;
	const std::vector<int>& counts_;
	const std::vector<double>& eff_lens_;
	const WeightMap& weight_map_;
    std::vector<double> alpha_;
    std::vector<double> rho_;
    bool rho_set_;
};

#endif // KALLISTO_EMALGORITHM_H
