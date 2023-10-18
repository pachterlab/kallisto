#ifndef BIFROST_BLOCKEDBLOOMFILTER_HPP
#define BIFROST_BLOCKEDBLOOMFILTER_HPP

#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <unordered_set>

#include "fastmod.h"
#include "libpopcnt.h"
#include "wyhash.h"

#define BBF_NB_BITS_BLOCK (0x800ULL)
#define BBF_MASK_BITS_BLOCK (0x7ffULL)
#define BBF_NB_ELEM_BLOCK (32)

#define BBF_OVERLOAD_RATIO (0.55)

class DualBlockedBloomFilter; // Pre-declaration for friend-class declaration

class BlockedBloomFilter {

    friend class DualBlockedBloomFilter;

    public:

        BlockedBloomFilter();
        BlockedBloomFilter(const size_t nb_elem, const size_t bits_per_elem);
        BlockedBloomFilter(const BlockedBloomFilter& o);
        BlockedBloomFilter(BlockedBloomFilter&& o);

        ~BlockedBloomFilter();

        BlockedBloomFilter& operator=(const BlockedBloomFilter& o);
        BlockedBloomFilter& operator=(BlockedBloomFilter&& o);

        void initialize(const size_t nb_elem, const size_t bits_per_elem);

        int contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit) const;
        std::array<int64_t, 4> contains_bids(const uint64_t (&kmh)[4], const uint64_t minh, const int limit) const;

        int64_t contains_bids(const uint64_t kmh, const uint64_t minh) const;

        bool write(FILE *fp) const;
        bool read(FILE *fp);

        void clear();
        void reset();

        DualBlockedBloomFilter transferToDBBF(const uint64_t idx_bbf);

        inline bool contains(const uint64_t kmh, const uint64_t minh) const {

            return (contains_bids(kmh, minh) != -1);
        }

        inline bool insert(const uint64_t kmh, const uint64_t minh, const bool multi_threaded = false) {

            return static_cast<bool>((multi_threaded ? insert_par(kmh, minh) : insert_unpar(kmh, minh)) & 0x1ULL);
        }

        inline uint64_t getNbBlocks() const {

            return (blocks_ + 1); // Number of blocks +  the one bucket of overflowing k-mers
        }

        inline double getOccupancy(const uint64_t block_id) const {

            if (block_id >= blocks_) return 0.0;

            return (static_cast<double>(table_[block_id].bits_occupancy) / static_cast<double>(BBF_NB_BITS_BLOCK));
        }

        inline void printOccupancy() const {

            size_t nb_overloaded = 0;
            size_t nb_underloaded = 0;

            for (uint64_t i = 0; i != blocks_; ++i) {

                const double occupancy = getOccupancy(i);

                nb_overloaded += (occupancy > BBF_OVERLOAD_RATIO);
                nb_underloaded += (occupancy < 0.35);

                std::cout << "[" << i << "] = " << (occupancy * 100.0) << "%" << std::endl;
            }

            std::cout << (static_cast<double>(nb_overloaded) / static_cast<double>(blocks_)) * 100 << " % blocks are overloaded." <<std::endl;
            std::cout << (static_cast<double>(nb_underloaded) / static_cast<double>(blocks_)) * 100 << " % blocks are underloaded." <<std::endl;
        }

        inline double getFPrate() const {

            if ((table_ == nullptr) || (k_ == 0) || (nb_bits_per_elem == 0)) return 1.0;

            return fpp(nb_bits_per_elem, k_);
        }

    protected:

        struct BBF_Block {

            BBF_Block() {

                clear();
            }

            inline void clear() {

                bits_occupancy = 0;

                lck.clear();
                memset(block, 0, BBF_NB_ELEM_BLOCK * sizeof(uint64_t));
            }

            inline void lock() {

                while (lck.test_and_set(std::memory_order_acquire));
            }

            inline void unlock() {

                lck.clear(std::memory_order_release);
            }

            uint64_t block[BBF_NB_ELEM_BLOCK];
            uint64_t bits_occupancy;
            std::atomic_flag lck = ATOMIC_FLAG_INIT;
        };

        BBF_Block* table_; //Array of Bloom filter blocks

        uint64_t blocks_; //Nb blocks

        uint64_t nb_bits_per_elem;

        int k_; //Nb hash functions

        __uint128_t M_u64; // for fast division/modulo

        uint64_t seed1, seed2; // Random seeds for hash functions

        std::unordered_set<uint64_t> ush;

        std::atomic_flag lck_ush = ATOMIC_FLAG_INIT;

        void init_arrays();

        inline double fpp(const size_t bits, const int k) const {

            return pow(1 - exp((-static_cast<double>(k))/static_cast<double>(bits)), static_cast<double>(k));
        }

        uint64_t insert_par(const uint64_t kmer_hash, const uint64_t min_hash);
        uint64_t insert_unpar(const uint64_t kmer_hash, const uint64_t min_hash);
};

class DualBlockedBloomFilter {

    friend class BlockedBloomFilter;

    public:

        DualBlockedBloomFilter();
        DualBlockedBloomFilter(const size_t nb_elem, const size_t bits_per_elem);
        DualBlockedBloomFilter(const DualBlockedBloomFilter& o);
        DualBlockedBloomFilter(DualBlockedBloomFilter&& o);

        ~DualBlockedBloomFilter();

        DualBlockedBloomFilter& operator=(const DualBlockedBloomFilter& o);
        DualBlockedBloomFilter& operator=(DualBlockedBloomFilter&& o);

        void initialize(const size_t nb_elem, const size_t bits_per_elem);

        int contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit, const uint64_t idx_bbf) const;
        std::array<int64_t, 4> contains_bids(const uint64_t (&kmh)[4], const uint64_t minh, const int limit, const uint64_t idx_bbf) const;

        int64_t contains_bids(const uint64_t kmh, const uint64_t minh, const uint64_t idx_bbf) const;

        BlockedBloomFilter transferToBBF(const uint64_t idx_bbf);

        bool write(FILE *fp) const;
        bool writeAsBBF(FILE *fp, const uint64_t idx_bbf) const;
        bool read(FILE *fp);
        bool readFromBBF(FILE *fp, const uint64_t idx_bbf);

        void clear();
        void reset();

        inline bool contains(const uint64_t kmh, const uint64_t minh, const uint64_t idx_bbf) const {

            return (contains_bids(kmh, minh, idx_bbf) != -1);
        }

        inline bool insert(const uint64_t kmh, const uint64_t minh, const uint64_t idx_bbf, const bool multi_threaded = false) {

            return static_cast<bool>((multi_threaded ? insert_par(kmh, minh, idx_bbf) : insert_unpar(kmh, minh, idx_bbf)) & 0x1ULL);
        }

        inline uint64_t getNbBlocks() const {

            return (blocks_ + 1); // Number of blocks +  the one bucket of overflowing k-mers
        }

        inline double getOccupancy(const uint64_t block_id, const uint64_t idx_bbf) const {

            if (block_id >= blocks_) return 0.0;

            const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

            return (static_cast<double>(table_[(block_id << 1) + idx_bbf_norm].bits_occupancy) / static_cast<double>(BBF_NB_BITS_BLOCK));
        }

        inline void printOccupancy(const uint64_t idx_bbf) const {

            const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

            size_t nb_overloaded = 0;
            size_t nb_underloaded = 0;

            for (uint64_t i = 0; i != blocks_; ++i) {

                const double occupancy = getOccupancy(i, idx_bbf);

                nb_overloaded += (occupancy > BBF_OVERLOAD_RATIO);
                nb_underloaded += (occupancy < 0.35);

                std::cout << "[" << idx_bbf_norm << "][" << i << "] = " << (occupancy * 100.0) << "%" << std::endl;
            }

            std::cout << (static_cast<double>(nb_overloaded) / static_cast<double>(blocks_)) * 100 << " % blocks are overloaded." <<std::endl;
            std::cout << (static_cast<double>(nb_underloaded) / static_cast<double>(blocks_)) * 100 << " % blocks are underloaded." <<std::endl;
        }

        inline double getFPrate() const {

            if ((table_ == nullptr) || (k_ == 0) || (nb_bits_per_elem == 0)) return 1.0;

            return fpp(nb_bits_per_elem, k_);
        }

    protected:

        struct BBF_Block {

            BBF_Block() {

                clear();
            }

            inline void clear() {

                bits_occupancy = 0;

                lck.clear();
                memset(block, 0, BBF_NB_ELEM_BLOCK * sizeof(uint64_t));
            }

            inline void lock() {

                while (lck.test_and_set(std::memory_order_acquire));
            }

            inline void unlock() {

                lck.clear(std::memory_order_release);
            }

            uint64_t block[BBF_NB_ELEM_BLOCK];
            uint64_t bits_occupancy;
            std::atomic_flag lck = ATOMIC_FLAG_INIT;
        };

        BBF_Block* table_; //Array of Bloom filter blocks

        uint64_t blocks_; //Nb blocks

        uint64_t nb_bits_per_elem;

        int k_; //Nb hash functions

        __uint128_t M_u64; // for fast division/modulo

        uint64_t seed1, seed2; // Random seeds for hash functions

        std::unordered_set<uint64_t> ush[2];

        std::atomic_flag lck_ush[2] = {ATOMIC_FLAG_INIT, ATOMIC_FLAG_INIT};

        void init_arrays();

        inline double fpp(size_t bits, int k) const {

            return pow(1-exp(-((double)k)/((double)bits)),(double)k);
        }

        uint64_t insert_par(const uint64_t kmer_hash, const uint64_t min_hash, const uint64_t idx_bbf);
        uint64_t insert_unpar(const uint64_t kmer_hash, const uint64_t min_hash, const uint64_t idx_bbf);
};

/*class CountingBlockedBloomFilter : public BlockedBloomFilter {

    public:

        CountingBlockedBloomFilter();
        CountingBlockedBloomFilter(size_t nb_elem, size_t bits_per_elem);
        CountingBlockedBloomFilter(const CountingBlockedBloomFilter& o);
        CountingBlockedBloomFilter(CountingBlockedBloomFilter&& o);

        ~CountingBlockedBloomFilter();

        CountingBlockedBloomFilter& operator=(const CountingBlockedBloomFilter& o);
        CountingBlockedBloomFilter& operator=(CountingBlockedBloomFilter&& o);

        void clear();

        bool write(FILE *fp) const;
        bool read(FILE *fp);

        int64_t contains_bids(const uint64_t kmh, const uint64_t minh, const uint64_t min_count) const;

        int contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit, const uint64_t min_count) const;

        inline bool contains(const uint64_t kmh, const uint64_t minh, const uint64_t min_count) const {

            return (contains_bids(kmh, minh, min_count) != -1);
        }

        inline bool insert(const uint64_t kmh, const uint64_t minh, const bool multi_threaded = false) {

            return (((multi_threaded ? insert_par(kmh, minh) : insert_unpar(kmh, minh)) & 0x1ULL) == 1);
        }

    private:

        uint64_t insert_par(const uint64_t kmer_hash, const uint64_t min_hash);
        uint64_t insert_unpar(const uint64_t kmer_hash, const uint64_t min_hash);

        void init_arrays();

        uint64_t elems_per_block; //Nb elements

        uint64_t* hashbit; //Bitmap of used hashes for approximate counts (BBHash-style)
        uint8_t* counts; //Counts
};*/

#endif // BFG_BLOCKEDBLOOMFILTER_HPP
