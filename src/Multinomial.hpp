#ifndef KALLISTO_MULTINOMIAL_H
#define KALLISTO_MULTINOMIAL_H

#include <stdexcept>
#include <random>

class Multinomial {
    public:
        Multinomial(const std::vector<int>& counts, size_t seed = 42) :
            counts_(counts),
            gen_(seed),
            dd_(counts_.begin(), counts_.end()),
            n_(0)
        {
            for (auto c : counts_) {
                n_ += c;
            }
        }

        /**
         * Generate a sample from the multinomial distribution
         *
         * Note that calling this function with nsamp != n() will not result in
         * a "proper" multinomial, though, it might be useful as long as it is
         * understood that only samples that have the same nsamp are
         * comparable. Call sample() for a standard multinomial.
         *
         * @param nsamp the number of samples. default == -1, which means it
         * will default to n_
         * @return a vector of counts
         *
         */
        std::vector<int> sample(int nsamp) {
            if (nsamp < 1) {
                throw std::domain_error("nsamp must be -1 or >=1");
            }

            std::vector<int> samp(counts_.size(), 0.0);
            for (auto i = 0; i < nsamp; ++i) {
                ++samp[dd_(gen_)];
            }

            return samp;
        }

        /**
         * @return a vector of ints with nsamp == n()
         */
        std::vector<int> sample() {
            return sample(n_);
        }

        int n() const { return n_; }
        const std::vector<int>& counts() { return counts_; }

    private:
        const std::vector<int>& counts_;
        std::default_random_engine gen_;
        std::discrete_distribution<int> dd_;
        int n_;
};

#endif
