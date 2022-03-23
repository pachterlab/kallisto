#ifndef BIFROST_BLOCKEDBLOOMFILTER_HPP
#define BIFROST_BLOCKEDBLOOMFILTER_HPP

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <unordered_set>

#include "libdivide.h"
#include "libpopcnt.h"
#include "wyhash.h"

#define NB_BITS_BLOCK (0x800ULL)
#define MASK_BITS_BLOCK (0x7ffULL)
#define NB_ELEM_BLOCK (32)

/* Short description:
 *  - Extended BloomFilter which hashes into 64-bit blocks
 *    that can be accessed very fast from the CPU cache
 * */

/*#if defined(__AVX2__)
#include <x86intrin.h>
#endif*/

class BlockedBloomFilter {

    public:

        BlockedBloomFilter();
        BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem);
        BlockedBloomFilter(const BlockedBloomFilter& o);
        BlockedBloomFilter(BlockedBloomFilter&& o);

        ~BlockedBloomFilter();

        BlockedBloomFilter& operator=(const BlockedBloomFilter& o);
        BlockedBloomFilter& operator=(BlockedBloomFilter&& o);

        int contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit) const;
        bool contains(const uint64_t kmh, const uint64_t minh) const;
        int contains_bids(const uint64_t kmh, const uint64_t minh) const;

        inline bool insert(const uint64_t kmh, const uint64_t minh, const bool multi_threaded = false) {

            return (multi_threaded ? insert_par(kmh, minh) : insert_unpar(kmh, minh));
        }

        bool WriteBloomFilter(FILE *fp) const;
        bool ReadBloomFilter(FILE *fp);

        void clear();

        inline uint64_t getNbBlocks() const { return blocks_; }

        inline double getOccupancy(const uint64_t block_id) const {

            if (block_id >= blocks_) return 0.0;

            return (static_cast<double>(table_[block_id].bits_occupancy) / static_cast<double>(NB_BITS_BLOCK));
        }

        inline void printOccupancy() const {

            size_t nb_overloaded = 0;
            size_t nb_underloaded = 0;

            for (uint64_t i = 0; i != blocks_; ++i) {

                const double occupancy = getOccupancy(i);

                nb_overloaded += (occupancy > 0.66);
                nb_underloaded += (occupancy < 0.35);

                std::cout << "[" << i << "] = " << (occupancy * 100.0) << "%" << std::endl;
            }

            std::cout << (static_cast<double>(nb_overloaded) / static_cast<double>(blocks_)) * 100 << " % blocks are overloaded." <<std::endl;
            std::cout << (static_cast<double>(nb_underloaded) / static_cast<double>(blocks_)) * 100 << " % blocks are underloaded." <<std::endl;
        }

    private:

        struct BBF_Block {

            BBF_Block() {

                clear();
            }

            inline void clear() {

                bits_occupancy = 0;

                lck.clear();
                memset(block, 0, NB_ELEM_BLOCK * sizeof(uint64_t));
            }

            inline void lock() {

                while (lck.test_and_set(std::memory_order_acquire));
            }

            inline void unlock() {

                lck.clear(std::memory_order_release);
            }

            uint64_t block[NB_ELEM_BLOCK];
            uint64_t bits_occupancy;
            std::atomic_flag lck = ATOMIC_FLAG_INIT;
        };

        BBF_Block* table_; //Bit array

        uint64_t blocks_; //Nb blocks

        int k_; //Nb hash functions

        libdivide::divider<uint64_t> fast_div_; // fast division

        uint64_t seed1, seed2; // Random seeds for hash functions

        std::unordered_set<uint64_t> ush;

        std::atomic_flag lck_ush = ATOMIC_FLAG_INIT;

        /*#if defined(__AVX2__)

        // The following functions are from the libpopcnt library 
        // (https://github.com/kimwalisch/libpopcnt) by Kim Walisch.
        // The algorithm is based on the paper "Faster Population Counts
        // using AVX2 Instructions" by Daniel Lemire, Nathan Kurz and
        // Wojciech Mula (23 Nov 2016). https://arxiv.org/abs/1611.07612

        static inline void CSA256(__m256i* h, __m256i* l, __m256i a, __m256i b, __m256i c) {

            __m256i u = a ^ b;
            *h = (a & b) | (u & c);
            *l = u ^ c;
        }

        static inline __m256i popcnt256(__m256i v) {

            __m256i lookup1 = _mm256_setr_epi8(
              4, 5, 5, 6, 5, 6, 6, 7,
              5, 6, 6, 7, 6, 7, 7, 8,
              4, 5, 5, 6, 5, 6, 6, 7,
              5, 6, 6, 7, 6, 7, 7, 8
            );

            __m256i lookup2 = _mm256_setr_epi8(
              4, 3, 3, 2, 3, 2, 2, 1,
              3, 2, 2, 1, 2, 1, 1, 0,
              4, 3, 3, 2, 3, 2, 2, 1,
              3, 2, 2, 1, 2, 1, 1, 0
            );

            __m256i low_mask = _mm256_set1_epi8(0x0f);
            __m256i lo = v & low_mask;
            __m256i hi = _mm256_srli_epi16(v, 4) & low_mask;
            __m256i popcnt1 = _mm256_shuffle_epi8(lookup1, lo);
            __m256i popcnt2 = _mm256_shuffle_epi8(lookup2, hi);

            return _mm256_sad_epu8(popcnt1, popcnt2);
        }

        static inline uint64_t popcnt_avx2(const __m256i* data, uint64_t size) {

            __m256i cnt = _mm256_setzero_si256();
            __m256i ones = _mm256_setzero_si256();
            __m256i twos = _mm256_setzero_si256();
            __m256i fours = _mm256_setzero_si256();
            __m256i eights = _mm256_setzero_si256();
            __m256i sixteens = _mm256_setzero_si256();
            __m256i twosA, twosB, foursA, foursB, eightsA, eightsB;

            uint64_t i = 0;
            uint64_t limit = size - size % 16;
            uint64_t* cnt64;

            for (; i < limit; i += 16) {

                CSA256(&twosA, &ones, ones, data[i+0], data[i+1]);
                CSA256(&twosB, &ones, ones, data[i+2], data[i+3]);
                CSA256(&foursA, &twos, twos, twosA, twosB);
                CSA256(&twosA, &ones, ones, data[i+4], data[i+5]);
                CSA256(&twosB, &ones, ones, data[i+6], data[i+7]);
                CSA256(&foursB, &twos, twos, twosA, twosB);
                CSA256(&eightsA, &fours, fours, foursA, foursB);
                CSA256(&twosA, &ones, ones, data[i+8], data[i+9]);
                CSA256(&twosB, &ones, ones, data[i+10], data[i+11]);
                CSA256(&foursA, &twos, twos, twosA, twosB);
                CSA256(&twosA, &ones, ones, data[i+12], data[i+13]);
                CSA256(&twosB, &ones, ones, data[i+14], data[i+15]);
                CSA256(&foursB, &twos, twos, twosA, twosB);
                CSA256(&eightsB, &fours, fours, foursA, foursB);
                CSA256(&sixteens, &eights, eights, eightsA, eightsB);

                cnt = _mm256_add_epi64(cnt, popcnt256(sixteens));
            }

            cnt = _mm256_slli_epi64(cnt, 4);
            cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(popcnt256(eights), 3));
            cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(popcnt256(fours), 2));
            cnt = _mm256_add_epi64(cnt, _mm256_slli_epi64(popcnt256(twos), 1));
            cnt = _mm256_add_epi64(cnt, popcnt256(ones));

            for(; i < size; i++) cnt = _mm256_add_epi64(cnt, popcnt256(data[i]));

            cnt64 = (uint64_t*) &cnt;

            return (cnt64[0] + cnt64[1] + cnt64[2] + cnt64[3]);
        }

        uint64_t hashes_mask[4];

        static const __m256i mask_and_div; // All MASK_BITS_BLOCK LSB of each 16 bits word
        static const __m256i mask_and_mod; // All 4 LSB of each 16 bits word
        static const __m256i one2shift_lsb; // All 1 LSB of each 32 bits word
        static const __m256i one2shift_msb; // Set the 17th LSB bit of each 32 bits word
        static const __m256i mask_lsb; // All 16 LSB bits of each 32 bits word

        uint16_t mask_ksup4, mask_ksup8, mask_ksup12;

        #endif*/

        void init_table();

        inline double fpp(size_t bits, int k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }

        bool insert_par(const uint64_t kmer_hash, const uint64_t min_hash);
        bool insert_unpar(const uint64_t kmer_hash, const uint64_t min_hash);
};

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
