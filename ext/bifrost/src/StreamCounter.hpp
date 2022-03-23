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

        StreamCounter(const double e_, const int seed_ = 0) :   e(e_), seed(seed_), sumCount(0),
                                                                shift_div_idx(__builtin_ffsll(block_sz) - 1), mask_mod_idx(block_sz - 1) {

            const size_t numcounts = max(static_cast<size_t>(48.0/(e*e) + 1), static_cast<size_t>(8192));

            sz = rndup((numcounts + countsPerLong - 1) / countsPerLong);
            mask = (sz * countsPerLong) - 1;
            maxCount = sz * countsPerLong * maxVal;

            std::memset(M, 0, MAX_TABLE * sizeof(size_t));

            blocks = new Block[((sz * MAX_TABLE) + block_sz - 1) / block_sz];

            for (size_t i = 0; i < MAX_TABLE; ++i) atomic_M[i].val = 0;

        }

        StreamCounter(const StreamCounter& o) : e(o.e), seed(o.seed), mask(o.mask), maxCount(o.maxCount), sumCount(o.sumCount),
                                                shift_div_idx(o.shift_div_idx), mask_mod_idx(o.mask_mod_idx), sz(o.sz) {

            std::memcpy(M, o.M, MAX_TABLE * sizeof(size_t));

            const size_t nb_blocks = ((sz * MAX_TABLE) + block_sz - 1) / block_sz;

            blocks = new Block[nb_blocks];

            for (size_t i = 0; i < nb_blocks; ++i) blocks[i] = o.blocks[i];
            for (size_t i = 0; i < MAX_TABLE; ++i) atomic_M[i].val = o.atomic_M[i].val.load();
        }

        ~StreamCounter() {

            delete[] blocks;
        }

        BFG_INLINE int getSeed() const {

            return seed;
        }

        void init_threads() {

            const size_t nb_blocks = ((sz * MAX_TABLE) + block_sz - 1) / block_sz;

            for (size_t i = 0; i < nb_blocks; ++i) blocks[i].clear();
            for (size_t i = 0; i < MAX_TABLE; ++i) atomic_M[i].val = 0;
        }

        void release_threads() {

            const size_t nb_blocks = ((sz * MAX_TABLE) + block_sz - 1) / block_sz;

            sumCount = 0;

            for (size_t i = 0; i < nb_blocks; ++i) sumCount += blocks[i].sumCount;
            for (size_t i = 0; i < MAX_TABLE; ++i) M[i] = atomic_M[i].val.load();
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

            ++sumCount;

            if (M[w] < maxCount){

                const size_t val = (block.table[pos_block] & bitmask) >> bitshift;

                if (val != maxVal) {

                    block.table[pos_block] = ((((val + 1) & maxVal) << bitshift) & bitmask) | (block.table[pos_block] & ~bitmask);

                    ++M[w];
                }
            }
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

            ++block.sumCount;

            if (block.M[w] < maxCount){

                const size_t val = (block.table[pos_block] & bitmask) >> bitshift;

                if (val != maxVal) {

                    block.table[pos_block] = ((((val + 1) & maxVal) << bitshift) & bitmask) | (block.table[pos_block] & rev_bitmask);
                    block.M[w] = ++atomic_M[w].val;
                }
            }

            block.lock.release();
        }

        bool join(const StreamCounter& o) {

            if ((sz != o.sz) || (seed != o.seed)) return false;

            const size_t R = sz * countsPerLong;

            for (size_t i = 0; i < MAX_TABLE; ++i) {

                M[i] = 0;

                for (size_t j = 0; j < R; ++j) {

                    setVal(j, i, getVal(j,i) + o.getVal(j,i));

                    M[i] += getVal(j,i);
                }
            }

            sumCount += o.sumCount;

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

        BFG_INLINE size_t F1() const {

            return sumCount;
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

        struct Block {

            SpinLock lock;

            size_t sumCount;

            size_t table[block_sz];
            size_t M[MAX_TABLE];

            Block() {

                clear();
            }

            Block(const Block& o) : sumCount(o.sumCount) {

                std::memcpy(table, o.table, block_sz * sizeof(size_t));
                std::memcpy(M, o.M, MAX_TABLE * sizeof(size_t));
            }

            Block& operator=(const Block& o) {

                if (this != &o) {

                    sumCount = o.sumCount;

                    std::memcpy(table, o.table, block_sz * sizeof(size_t));
                    std::memcpy(M, o.M, MAX_TABLE * sizeof(size_t));
                }

                return *this;
            }

            void clear() {

                sumCount = 0;

                std::memset(table, 0, block_sz * sizeof(size_t));
                std::memset(M, 0, MAX_TABLE * sizeof(size_t));
            }
        };

        struct atomic_uint64_t {

            std::atomic<size_t> val;
            const int padding[14];

            atomic_uint64_t() : padding{0} {

                val = 0;
            }
        };

        int seed;

        double e;

        size_t sz;
        size_t maxCount;
        size_t sumCount;
        size_t mask;

        size_t nb_blocks;

        const size_t shift_div_idx;
        const size_t mask_mod_idx;

        size_t M[MAX_TABLE];
        atomic_uint64_t atomic_M[MAX_TABLE];

        Block* blocks;
};

/*class StreamCounter {

    public:

        StreamCounter(double e_, int seed_ = 0) : MAX_TABLE(32), maxVal(3ULL), countWidth(2), countsPerLong(32), e(e_), seed(seed_), sumCount(0) {

            const size_t numcounts = max(static_cast<size_t>(48.0/(e*e) + 1), static_cast<size_t>(8192)); // approx 3 std-dev, true with 0.001 prob.

            size = (numcounts + countsPerLong - 1) / countsPerLong; // size is number of uint64_t's use
            size = roundUpPowerOfTwo(size);
            mask = (size * countsPerLong) - 1;
            maxCount = size * countsPerLong * maxVal;

            M = new size_t[MAX_TABLE]();
            table = new uint64_t[size * MAX_TABLE]();

        }

        StreamCounter(const StreamCounter& o) : seed(o.seed), e(o.e), size(o.size), maxCount(o.maxCount), mask(o.mask),
                                                MAX_TABLE(o.MAX_TABLE), countWidth(o.countWidth), countsPerLong(o.countsPerLong),
                                                maxVal(o.maxVal), sumCount(o.sumCount) {
            // copy constructor, creates object of same size with same seed, but empty data
            M = new size_t[MAX_TABLE]();
            table = new uint64_t[size * MAX_TABLE]();

            memcpy(M, o.M, MAX_TABLE * sizeof(size_t));
            memcpy(table, o.table, size * MAX_TABLE * sizeof(uint64_t));
        }

        ~StreamCounter() {

            delete[] table;
            delete[] M;
        }

        inline int getSeed() const { return seed; }

        void operator()(const uint64_t hashval) {

            ++sumCount;

            // hashval is XXX .. XXX1000.. (w times) ..00 0
            size_t w = __builtin_ctzll(hashval);

            if (w >= MAX_TABLE) w = MAX_TABLE - 1;
            if (M[w] == maxCount) return;

            // hashval is now XXXX...XX, random
            const uint64_t hval = hashval >> (w + 1); // shift away pattern of XXX10000
            const uint64_t index = hval & mask;
            const uint64_t val = getVal(index,w);

            if (val != maxVal) { // max count

                setVal(index, w, val+1);
                ++M[w];
            }
        }

        bool join(const StreamCounter& o) {

            if (size != o.size || seed != o.seed) return false;

            for (size_t i = 0; i < MAX_TABLE; ++i) {

                M[i] = 0;

                for (size_t j = 0; j < size * countsPerLong; ++j) {

                    setVal(j, i, getVal(j,i) + o.getVal(j,i));

                    M[i] += getVal(j,i);
                }
            }

            sumCount += o.sumCount;

            return true;
        }

        size_t F0() const {

            int n = 0;
            double sum = 0;
            double limit = 0.2;

            const size_t R = size * countsPerLong;

            while (n == 0 && limit > 1e-8) {

                for (size_t i = 0; i < MAX_TABLE; ++i) {

                    size_t ts = 0;

                    for (size_t j = 0; j < R; ++j) {

                        if (getVal(j,i) > 0) ++ts;
                    }

                    if (ts == 0) { // nothing in this level, hack to break out of loop

                        limit = 0;
                        break;
                    }

                    if ((ts <= (1 - limit) * R) && ((ts >= limit * R) || (i == 0))) {

                        sum += (log(1.0 - ts/((double) R)) / log(1.0-1.0 / R)) * pow(2.0, i+1);
                        ++n;
                        break;
                    }
                }

                limit = limit / 1.5;
            }

            if (n == 0) return 0;
            return (size_t)(sum/n);
        }

        inline size_t F1() const { return sumCount; }

        size_t f1() const {

            int n = 0;

            double sum = 0;
            double limit = 0.2;

            const size_t R = size * countsPerLong;

            while (n == 0 && limit > 1e-8) {

                for (size_t i = 0; i < MAX_TABLE; ++i) {

                    size_t r1 = 0, r0 = 0;

                    for (size_t j = 0; j < R; ++j) {

                        const uint64_t val = getVal(j,i);

                        if (val == 0) ++r0;
                        if (val == 1) ++r1;

                    }

                    if (r0 == R) { // empty level

                        limit = 0;
                        break;
                    }

                    if (((r0 <= ((1 - limit) * R)) || (i == 0)) && (r0 >= (limit * R))) {
                        // i==0 takes care of small first levels where R is too large
                        sum += (R-1) * (r1/((double) r0)) * pow(2.0,i+1);
                        ++n;
                        break;
                    }
                }

                limit = limit/1.5;
            }

            if (n == 0) return 0;
            return (size_t) (sum/n);
        }

        std::string humanReport() const {

            std::stringstream s;

            const size_t eF0 = F0();
            const size_t ef1 = f1();
            const size_t eF1 = sumCount;

            s << readable(eF0-ef1) << " repeated, " << readable(eF0) << " distinct, " <<
            readable(ef1) << " singletons, " << readable(eF1) << " total k-mers processed";

            return s.str();
        }

        std::string readable(size_t x) const {

            std::stringstream s;

            if (x < (1ULL<<10)) s << x;
            else if (x < (1ULL<<20)) s << (x/1024) << "K";
            else if (x < (1ULL<<30)) s << (x/(1ULL<<20)) << "M";
            else s << (x/(1ULL<<30)) << "G";

            return s.str();

        }

        std::string report(bool useTSV = false) const {

            std::stringstream s;

            if (useTSV) s << F0() << "\t" << f1() << "\t" << F1() << std::endl;
            else {

                s << "F0 = " << F0() << std::endl;
                s << "f1 = " << f1() << std::endl;
                s << "F1 = " << sumCount << std::endl;
            }

            return s.str();
        }

        bool writeBinary(const std::string& fn) {

            std::ofstream out;

            out.open(fn.c_str(), std::ios::out | std::ios::binary);

            if (!out.is_open()) {

                std::cerr << "Error: could not write to file " << fn << std::endl;
                return false;
            }

            // 1. write out seed
            out.write((char*)&seed, sizeof(seed));
            // 2. write out size
            out.write((char*)&size, sizeof(size));
            // 3. write out e
            out.write((char*)&e, sizeof(e));
            // 3.5 write out sumCount
            out.write((char*)&sumCount, sizeof(sumCount));
            // 4. write out MAX_TABLE
            out.write((char*)&MAX_TABLE, sizeof(MAX_TABLE));
            // 5. write out M
            out.write((char*)M, MAX_TABLE*sizeof(M[0]));
            // 6. write out table
            out.write((char*)table, size*MAX_TABLE*sizeof(table[0]));

            out.flush();
            out.close();

            return true;
        }

        bool loadBinary(const std::string& fn) {

            std::ifstream in;
            const size_t oldsize = size;

            in.open(fn.c_str(), std::ios::in | std::ios::binary);

            // 1. read seed
            in.read((char*)&seed,sizeof(seed));
            // 2. read size
            in.read((char*)&size,sizeof(size));
            // 3. read e
            in.read((char*)&e, sizeof(e));
            // 3.5 read sumCount
            in.read((char*)&sumCount, sizeof(sumCount));

            size_t max_table;
            // 4. read MAX_TABLE
            in.read((char*)&max_table, sizeof(max_table));

            if(MAX_TABLE != max_table) {

                std::cerr <<"Error: Max table size doesn't match" << std::endl;
                std::cerr << "MAX_TABLE = " << MAX_TABLE << std::endl << "max_table = " << max_table << std::endl;
                exit(1);
            }

            // fill in other variables
            mask = (size * countsPerLong) - 1;
            maxCount = size * countsPerLong * maxVal;

            // allocate space
            if (oldsize != size) {

                if (M != 0) {

                    delete[] M;
                    M = new size_t[MAX_TABLE];
                }

                if (table != 0) {

                    delete[] table;
                    table = new uint64_t[size * MAX_TABLE];
                }
            }

            // 5. read M
            in.read((char*)M, MAX_TABLE * sizeof(M[0]));

            // 6. read T
            in.read((char*)table, size * MAX_TABLE * sizeof(table[0]));
            in.close();

            return true;
        }

    private:

        static size_t roundUpPowerOfTwo(size_t size) {

            --size;

            size |= size >> 1;
            size |= size >> 2;
            size |= size >> 4;
            size |= size >> 8;
            size |= size >> 16;
            size |= size >> 32;

            ++size;

            return size;
        }

        uint64_t getVal(const size_t index, const size_t w) const {

            const size_t wordindex = w * size + (index / countsPerLong);
            const size_t bitindex = index & (countsPerLong - 1) ;
            const uint64_t bitmask = maxVal << (countWidth * bitindex);

            return (table[wordindex] & bitmask) >> (countWidth * bitindex);
        }

        void setVal(const size_t index, const size_t w, uint64_t val) {

            if (val > maxVal) val = maxVal;

            const size_t wordindex = w * size + (index / countsPerLong);
            const size_t bitindex = index & (countsPerLong - 1);
            const uint64_t bitmask = maxVal << (countWidth * bitindex);

            table[wordindex] = (((val & maxVal) << (countWidth * bitindex)) & bitmask) | (table[wordindex] & ~bitmask);
        }

        int seed;
        double e;
        uint64_t* table;
        size_t sumCount;
        size_t *M;
        size_t size;
        size_t maxCount;
        uint64_t mask;
        const size_t MAX_TABLE;
        const size_t countWidth; // number of bits per count, even number
        const size_t countsPerLong;  // fix this
        const uint64_t maxVal; // has to be a power of 2-1
};*/

#endif
