#ifndef BIFROST_STREAMCOUNTER_HPP
#define BIFROST_STREAMCOUNTER_HPP

#include <sstream>
#include <stdint.h>
#include <string>
#include <cmath>
#include <fstream>
#include <assert.h>

#include "Common.hpp"

class StreamCounter {

    public:

        void clear() {

            e = 0.0;
            seed = 0;
            sz = 0;
            mask = 0;

            if (blocks != nullptr) {

                delete[] blocks;

                blocks = nullptr;
            }
        }

        StreamCounter() :   blocks(nullptr), shift_div_idx(__builtin_ffsll(block_sz) - 1), mask_mod_idx(block_sz - 1) {

            clear();
        }

        StreamCounter(const double e_, const int seed_ = 0) : blocks(nullptr), shift_div_idx(__builtin_ffsll(block_sz) - 1), mask_mod_idx(block_sz - 1) {

           initialize(e_, seed_);
        }

        StreamCounter(const StreamCounter& o) : e(o.e), seed(o.seed), mask(o.mask), sz(o.sz),
                                                shift_div_idx(o.shift_div_idx), mask_mod_idx(o.mask_mod_idx) {

            const size_t nb_blocks = ((sz * MAX_TABLE) + block_sz - 1) / block_sz;

            blocks = new Block[nb_blocks];

            for (size_t i = 0; i < nb_blocks; ++i) blocks[i] = o.blocks[i];
        }

        StreamCounter& operator=(const StreamCounter& o) {

            if (this != &o) {

                clear();

                e = o.e;
                seed = o.seed;
                mask = o.mask;
                sz = o.sz;

                const size_t nb_blocks = ((sz * MAX_TABLE) + block_sz - 1) / block_sz;

                blocks = new Block[nb_blocks];

                for (size_t i = 0; i < nb_blocks; ++i) blocks[i] = o.blocks[i];
            }

            return *this;
        };

        StreamCounter& operator=(StreamCounter&& o) {

            if (this != &o) {

                clear();

                e = o.e;
                seed = o.seed;
                mask = o.mask;
                sz = o.sz;
                blocks = o.blocks;

                o.blocks = nullptr;

                o.clear();
            }

            return *this;
        };

        ~StreamCounter() {

            clear();
        }

        BFG_INLINE int getSeed() const {

            return seed;
        }

        void initialize(const double e_, const int seed_ = 0) {

            clear();

            e = e_;
            seed = seed_;

            const size_t numcounts = max(static_cast<size_t>(48.0/(e*e) + 1), static_cast<size_t>(8192));

            sz = rndup((numcounts + countsPerLong - 1) / countsPerLong);
            mask = (sz * countsPerLong) - 1;

            blocks = new Block[((sz * MAX_TABLE) + block_sz - 1) / block_sz];
        }

        void update(const size_t h) {

            // hashval is XXX .. XXX1000.. (w times) ..00 0
            const size_t w = std::min(static_cast<size_t>(__builtin_ctzll(h)), MAX_TABLE - 1);

            // hashval is now XXXX...XX, random
            const size_t hval = h >> (w + 1); // shift away pattern of XXX10000
            const size_t index = hval & mask;
            const size_t wordindex = w * sz + (index / countsPerLong);
            const size_t pos_block = wordindex & mask_mod_idx;
            const size_t bitshift = countWidth * (index & (countsPerLong - 1));
            const size_t bitmask = maxVal << bitshift;

            Block& block = blocks[wordindex >> shift_div_idx];

            const size_t val = (block.table[pos_block] & bitmask) >> bitshift;

            if (val != maxVal) block.table[pos_block] = ((((val + 1) & maxVal) << bitshift) & bitmask) | (block.table[pos_block] & ~bitmask);
        }

        void update_p(const size_t h) {

            // hashval is XXX .. XXX1000.. (w times) ..00 0
            const size_t w = std::min(static_cast<size_t>(__builtin_ctzll(h)), MAX_TABLE - 1);

            // hashval is now XXXX...XX, random
            const size_t hval = h >> (w + 1); // shift away pattern of XXX10000
            const size_t index = hval & mask;
            const size_t wordindex = w * sz + (index / countsPerLong);
            const size_t pos_block = wordindex & mask_mod_idx;
            const size_t bitshift = countWidth * (index & (countsPerLong - 1));
            const size_t bitmask = maxVal << bitshift;
            const size_t rev_bitmask = ~bitmask;

            Block& block = blocks[wordindex >> shift_div_idx];

            block.lock.acquire();

            const size_t val = (block.table[pos_block] & bitmask) >> bitshift;

            if (val != maxVal) block.table[pos_block] = ((((val + 1) & maxVal) << bitshift) & bitmask) | (block.table[pos_block] & rev_bitmask);

            block.lock.release();
        }

        /*void update_p(const vector<size_t>& vh) {

            vector<h_field_t> vh_block;

            size_t prev_block_id = 0xffffffffffffffffULL;

            auto comp_hf = [](const h_field_t& p1, const h_field_t& p2) {

                return (p1.wordindex < p2.wordindex);
            };

            vh_block.reserve(vh.size());

            for (const auto h : vh) {

                h_field_t hf;

                hf.w = std::min(static_cast<size_t>(__builtin_ctzll(h)), MAX_TABLE - 1);
                hf.index = (h >> (hf.w + 1)) & mask;
                hf.wordindex = hf.w * sz + (hf.index / countsPerLong);
                hf.pos_block = hf.wordindex & mask_mod_idx;
                hf.block_id = hf.wordindex >> shift_div_idx;
                hf.bitshift = countWidth * (hf.index & (countsPerLong - 1));
                hf.bitmask = maxVal << hf.bitshift;
                hf.rev_bitmask = ~(hf.bitmask);

                vh_block.push_back(move(hf));
            }

            sort(vh_block.begin(), vh_block.end(), comp_hf);

            for (const auto& hf : vh_block) {

                Block& block = blocks[hf.block_id];

                if (hf.block_id != prev_block_id) {

                    if (prev_block_id != 0xffffffffffffffffULL) blocks[prev_block_id].lock.release();

                    prev_block_id = hf.block_id;

                    block.lock.acquire();
                }

                const size_t val = (block.table[hf.pos_block] & hf.bitmask) >> hf.bitshift;

                if (val != maxVal) block.table[hf.pos_block] = ((((val + 1) & maxVal) << hf.bitshift) & hf.bitmask) | (block.table[hf.pos_block] & hf.rev_bitmask);
            }

            if (prev_block_id != 0xffffffffffffffffULL) blocks[prev_block_id].lock.release();
        }*/

        bool join(const StreamCounter& o) {

            if ((sz != o.sz) || (seed != o.seed)) return false;

            const size_t R = sz * countsPerLong;

            for (size_t i = 0; i < MAX_TABLE; ++i) {

                for (size_t j = 0; j < R; ++j) setVal(j, i, getVal(j,i) + o.getVal(j,i));

            }

            return true;
        }

        size_t F0() const {

            int n = 0;

            double sum = 0.0;
            double limit = 0.2;

            const size_t R = sz * countsPerLong;

            while (n == 0 && limit > 1e-8) {

                for (size_t i = 0; i < MAX_TABLE; ++i) {

                    size_t ts = 0;

                    for (size_t j = 0; j < R; ++j) ts += static_cast<size_t>(getVal(j,i) > 0);

                    if (ts == 0) { // nothing in this level, hack to break out of loop

                        limit = 0;

                        break;
                    }

                    if ((ts <= (1 - limit) * R) && ((ts >= limit * R) || (i == 0))) {

                        sum += (log(1.0 - ts/static_cast<double>(R)) / log(1.0-1.0 / R)) * pow(2.0, i+1);
                        ++n;
                        break;
                    }
                }

                limit /= 1.5;
            }

            if (n == 0) return 0;

            return (size_t)(sum/n);
        }

        size_t f1() const {

            int n = 0;

            double sum = 0;
            double limit = 0.2;

            const size_t R = sz * countsPerLong;

            while (n == 0 && limit > 1e-8) {

                for (size_t i = 0; i < MAX_TABLE; ++i) {

                    size_t r1 = 0, r0 = 0;

                    for (size_t j = 0; j < R; ++j) {

                        const uint64_t val = getVal(j,i);

                        r0 += static_cast<size_t>(val == 0);
                        r1 += static_cast<size_t>(val == 1);
                    }

                    if (r0 == R) { // empty level

                        limit = 0;
                        break;
                    }

                    if (((r0 <= ((1 - limit) * R)) || (i == 0)) && (r0 >= (limit * R))) {
                        // i==0 takes care of small first levels where R is too large
                        sum += (R-1) * (r1/static_cast<double>(r0)) * pow(2.0, i+1);
                        ++n;
                        break;
                    }
                }

                limit /= 1.5;
            }

            if (n == 0) return 0;

            return (size_t) (sum/n);
        }

    private:

        size_t getVal(const size_t index, const size_t w) const {

            const size_t wordindex = w * sz + (index / countsPerLong);
            const size_t bitindex = index & (countsPerLong - 1) ;
            const size_t bitmask = maxVal << (countWidth * bitindex);

            return (blocks[wordindex >> shift_div_idx].table[wordindex & mask_mod_idx] & bitmask) >> (countWidth * bitindex);
        }

        void setVal(const size_t index, const size_t w, size_t val) {

            if (val > maxVal) val = maxVal;

            const size_t wordindex = w * sz + (index / countsPerLong);
            const size_t pos_block = wordindex & mask_mod_idx;
            const size_t bitindex = index & (countsPerLong - 1);
            const size_t bitmask = maxVal << (countWidth * bitindex);

            Block& block = blocks[wordindex >> shift_div_idx];

            block.table[pos_block] = (((val & maxVal) << (countWidth * bitindex)) & bitmask) | (block.table[pos_block] & ~bitmask);
        }

        static const size_t block_sz = 256; // Always a power of 2
        static const size_t MAX_TABLE = 32;
        static const size_t maxVal = 3; // has to be a power of 2-1
        static const size_t countWidth = 2;  // number of bits per count, even number
        static const size_t countsPerLong = 32; // ugly

        /*struct h_field_t {

            size_t w; // hashval is XXX .. XXX1000.. (w times) ..00 0
            size_t index;
            size_t wordindex;
            size_t pos_block;
            size_t bitshift;
            size_t bitmask;
            size_t rev_bitmask;
            size_t block_id;
        };*/

        struct Block {

            SpinLock lock;

            size_t table[block_sz];

            Block() {

                clear();
            }

            Block(const Block& o) {

                std::memcpy(table, o.table, block_sz * sizeof(size_t));
            }

            Block& operator=(const Block& o) {

                if (this != &o) {

                    std::memcpy(table, o.table, block_sz * sizeof(size_t));
                }

                return *this;
            }

            void clear() {

                std::memset(table, 0, block_sz * sizeof(size_t));
            }
        };

        int seed;

        double e;

        size_t sz;
        size_t mask;

        size_t nb_blocks;

        const size_t shift_div_idx;
        const size_t mask_mod_idx;

        Block* blocks;
};

#endif
