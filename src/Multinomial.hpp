#ifndef KALLISTO_MULTINOMIAL_H
#define KALLISTO_MULTINOMIAL_H

#include <random>

class Multinomial {
    public:
        Multinomial(const std::vector<int>& counts, size_t seed = 42) :
            counts_(counts),
            gen_(seed),
            n_(0),
            dd_(counts_.begin(), counts_.end())
        {
            for (auto c : counts_) {
                n_ += c;
            }
        }

        std::vector<int> sample() {
            std::vector<int> samp(counts_.size(), 0.0);
            for (auto i = 0; i < n_; ++i) {
                ++samp[dd_(gen_)];
            }

            return samp;
        }

        int n() { return n_; }

    private:
        const std::vector<int>& counts_;
        std::default_random_engine gen_;
        std::discrete_distribution<int> dd_;
        int n_;
};

#endif
