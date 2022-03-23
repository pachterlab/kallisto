#ifndef BIFROST_KMER_ITERATOR_HPP
#define BIFROST_KMER_ITERATOR_HPP

#include <iostream>
#include <iterator>
#include "Kmer.hpp"


/* Short description:
 *  - Easily iterate through kmers in a read
 *  - If the read contains any N, then the N is skipped and checked whether
 *    there is a kmer to the right of the N
 *  - iter->first gives the kmer, iter->second gives the position within the reads
 * */
class KmerIterator : public std::iterator<std::input_iterator_tag, std::pair<Kmer, int>, int> {

    public:

        KmerIterator() : str(nullptr), invalid(true), pos_s(0), pos_e(0) {

            p.second = 0;
        }

        KmerIterator(const char* s) : str(s), invalid(false), pos_s(0), pos_e(0) {

            p.second = 0;

           operator++();
        }

        KmerIterator& operator++();
        KmerIterator& operator+=(const int len);

        BFG_INLINE KmerIterator operator++(int) {

            const KmerIterator tmp(*this);

            operator++();

            return tmp;
        }

        BFG_INLINE bool operator==(const KmerIterator& o) const {

            if (invalid || o.invalid) return invalid && o.invalid;

            return (str == o.str) && (p == o.p);
        }

        BFG_INLINE bool operator!=(const KmerIterator& o) const {

            return !operator==(o);
        }

        BFG_INLINE const pair<Kmer, int>& operator*() const {

            return p;
        }

        BFG_INLINE const pair<Kmer, int>* operator->() const {

            return &p;
        }

    private:

        const char* str;
        bool invalid;
        pair<Kmer, int> p;
        int pos_s, pos_e;
};

template<class HF>
class KmerHashIterator {

    public:

        KmerHashIterator(const char* _s, const int _length, const int _k) : s(_s), n(_length), k(_k), hf(HF(_k)), p_(0,-1), invalid(true) {

            if ((_s != nullptr) || (n >= k)) {

                invalid = false;
                p_.second = 0;
                init();
            }
        }

        KmerHashIterator() : s(nullptr), n(0), k(0), hf(HF(0)), p_(0,-1), invalid(true) {}

        KmerHashIterator(const KmerHashIterator& o) : s(o.s), n(o.n), k(o.k), hf(o.hf), p_(o.p_), invalid(o.invalid) {}

        inline bool operator==(const KmerHashIterator& o) {

            if (invalid || o.invalid) return invalid && o.invalid;
            return (s == o.s) && (n == o.n) && (k == o.k) && (p_ == o.p_);
        }

        bool operator!=(const KmerHashIterator& o) { return !this->operator==(o); }

        inline void init() {

            invalid = (p_.second >= n - k + 1);

            if (!invalid){

                int j = p_.second + k - 1;

                while (j >= p_.second) {

                    const char c = s[j] & 0xDF; // mask lowercase bit

                    if (isDNA(c)) --j;
                    else {

                        p_.second = j + 1;
                        j += k;

                        if (p_.second >= n - k + 1) { // out of bounds

                            invalid = true;
                            p_ = {0,-1};
                            return;
                        }
                    }
                }

                hf.init(&s[p_.second]);
                p_.first = hf.hash();
            }
            else p_ = {0,-1};
        }

        KmerHashIterator& operator++() {

            if (invalid) return *this;

            ++(p_.second); // advance to next k-mer

            if (p_.second >= n - k + 1) { // out of bounds

                invalid = true;
                p_ = {0, -1};
            }
            else {

                const int j = p_.second + k - 1;
                const char c = s[j] & 0xDF; // mask lowercase bit

                if (isDNA(c)){

                    hf.update(s[j-k], s[j]);
                    p_.first = hf.hash();
                }
                else {

                    p_.second = j + 1;
                    init();
                }
            }

            return *this;
        }

        KmerHashIterator operator++(int) {

            KmerHashIterator tmp(*this);
            operator++();

            return tmp;
        }

        //Move iterator to next VALID position >= p_.second + length
        KmerHashIterator& operator+=(const int length){

            const size_t next_pos = p_.second + length;

            while (!invalid && p_.second < next_pos) operator++();

            return *this;
        }

        inline std::pair<uint64_t, int>& operator*() { return p_; }

        inline std::pair<uint64_t, int>* operator->() { return &p_; }

        const char *s; // K-mers are from a sequence s
        int n; // Length of sequence s
        int k; // Length of k-mers
        HF hf; // Rolling hash function for k-mers of s
        std::pair<uint64_t, int> p_; // <hash, position> current k-mer
        bool invalid; // If sequence is invalid (iterating on k-mers out of bounds, etc.)
};

#endif // BFG_KMER_ITERATOR_HPP
