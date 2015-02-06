#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "weights.h"

#include <algorithm>
#include <iostream>
#include <vector>


const double TOLERANCE = 1e-5;

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
	    num_trans_(idx.num_trans),
	    counts_(counts),
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

        std::cout << "[em]\tfishing for the right mixture (. = 50 rounds)" <<
            std::endl;

	    for (auto i = 0; i < n_iter; ++i) {
            if (i % 50 == 0) {
                std::cout << ".";
                if (i % 500 == 0 && i > 0) {
                    std::cout << std::endl;
                }
            }

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
            // TODO: consider what the right tolerance is
            if (eff_lens_[i] < TOLERANCE) {
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

    void write(const std::string& dir_out) const {
        const std::string out_fname = "/expression.txt";

        std::ofstream out;
        out.open(dir_out + out_fname, std::ios::out);

        if (!out.is_open()) {
            std::cerr << "Error opening '" << dir_out + out_fname << "'" <<
                std::endl;
            exit(1);
        }

        out.precision(15);

        out <<
            "target_id" << "\t" <<
            "kallisto_id" << "\t" <<
            "rho" << "\t" <<
            "tpm" << "\t" <<
            "expected_counts" <<
            std::endl;

        const double MILLION = 1e6;

        for (auto i = 0; i < rho_.size(); ++i) {
            out <<
                idx_.target_names_[i] << "\t" <<
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
	const Index &idx_;
	const std::vector<int>& counts_;
	const std::vector<double>& eff_lens_;
	const WeightMap& weight_map_;
    std::vector<double> alpha_;
    std::vector<double> rho_;
    bool rho_set_;
};

#endif // KALLISTO_EMALGORITHM_H
