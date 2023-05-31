#ifndef BIFROST_COMPRESSED_COVERAGE_HPP
#define BIFROST_COMPRESSED_COVERAGE_HPP

#include <cstring>
#include <stdint.h>
#include <vector>
#include <iostream>

#include "BitContainer.hpp"
#include "Common.hpp"

/* Short description:
 *  - Tagged pointer union that is either
 *    - a pointer to a char array that stores 2-bit integers
 *    - a local 2-bit integer array
 *  - The bits are stored either as
 *    pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|ppppppF0
 *    dddddddd|dddddddd|dddddddd|dddddddd|dddddddd|dddddddd|dddddddd|ssssssF1
 *
 *    - Last bit is 0 for pointer and 1 for local array
 *    - Second last bit is 1 for a full coverage and 0 for non-full
 *    - For local array last 8 except last 2 bits store size of the array
 *    - For the pointer version first 62 bits encode the pointer (last two bits
 *      are zeroed out before dereferencing)
 *    - The pointer points to an array of bytes, where the first 8 encode the
 *      size of the array used in uint32_t and the number of full positions.
 *    - The remainder of bytes are 2-bit encoded integers.
 *    - If the full bit is set then the pointer must be 0 and the memory released
 *
 */
class CompressedCoverage {

    public:

        CompressedCoverage(size_t sz=0, bool full=false);
        CompressedCoverage(const CompressedCoverage& o); // Copy constructors
        CompressedCoverage(CompressedCoverage&& o); // move constructors
        ~CompressedCoverage();

        CompressedCoverage& operator=(const CompressedCoverage& o); // Copy assignment
        CompressedCoverage& operator=(CompressedCoverage&& o); // Move assignment

        BFG_INLINE void clear() { releasePointer(); };

        void initialize(const size_t sz, const bool full);
        void initialize(const size_t sz, const size_t init_cov);

        void cover(size_t start, size_t end);
        void uncover(size_t start, size_t end);

        uint8_t covAt(const size_t index) const;

        bool isFull() const;
        void setFull();

        vector<pair<int, int>> splittingVector() const;
        pair<size_t, size_t> lowCoverageInfo() const;

        size_t size() const;
        std::string toString() const; // for debugging

        static void setFullCoverage(size_t cov_max) {

            if ((cov_max == 1) || (cov_max == 2)){

                cov_full = cov_max;
                localCoverageMask = 0;

                for (size_t i = 0; i < size_limit; ++i) localCoverageMask = (localCoverageMask << 2) | cov_full;
            }
        }

        BFG_INLINE static size_t getFullCoverage() { return cov_full; }

    private:

        static const size_t size_limit = 28; // 56 bit array, 28 2-bit integers

        static const uintptr_t tagMask = 1; // local array bit
        static const uintptr_t fullMask = 2; // full bit
        static const uintptr_t sizeMask = 0xFC; // 0b11111100
        static const uintptr_t pointerMask = ~(tagMask | fullMask); // rest of bits

        static size_t cov_full;

        static uintptr_t localCoverageMask;

        BFG_INLINE size_t round_to_bytes(const size_t len) const { return (len + 3) / 4; }

        BFG_INLINE uint8_t* get8Pointer() const { return reinterpret_cast<uint8_t*>(asBits & pointerMask); }
        BFG_INLINE uint32_t* get32Pointer() const { return reinterpret_cast<uint32_t*>(asBits & pointerMask); }
        BFG_INLINE const uint32_t* getConst32Pointer() const { return reinterpret_cast<const uint32_t*>(asBits & pointerMask); }

        void releasePointer();

        union {

            uint8_t* asPointer;
            uintptr_t asBits;
        };
};

/*class CompressedCoverage {

    public:

        CompressedCoverage(size_t sz_ = 0, bool full_ = false);

        CompressedCoverage(const CompressedCoverage& o); // Copy constructors
        CompressedCoverage(CompressedCoverage&& o); // move constructors

        CompressedCoverage& operator=(const CompressedCoverage& o); // Copy assignment
        CompressedCoverage& operator=(CompressedCoverage&& o); // Move assignment

        BFG_INLINE void clear() {

            sz = 0;
            bc.clear();
        };

        BFG_INLINE size_t size() const {

            return sz;
        };

        void cover(size_t start, size_t end);
        void uncover(size_t start, size_t end);

        uint8_t covAt(const size_t idx) const;
        size_t covSum() const;

        vector<pair<int, int>> splittingVector() const;

        BFG_INLINE static void setFullCoverage(const size_t cov_max) {

            cov_full = min(cov_max, static_cast<size_t>(2));
            cov_full = max(cov_full, static_cast<size_t>(1));
        }

        BFG_INLINE static size_t getFullCoverage() {

            return cov_full;
        }

        BFG_INLINE bool isFull() const {

            return bc.contains(0);
        }

        BFG_INLINE void setFull() {

            bc.clear();
            bc.add(0);
        }

    private:

        static size_t cov_full;

        size_t sz;

        BitContainer bc;
};*/

template<typename T> struct CompressedCoverage_t {

    CompressedCoverage_t(size_t sz = 0, bool full = false) : ccov(sz, full) {}
    CompressedCoverage_t(const CompressedCoverage& c) : ccov(c) {}
    CompressedCoverage_t(CompressedCoverage&& c) : ccov(move(c)) {}

    BFG_INLINE const T* getData() const { return &data; }
    BFG_INLINE T* getData() { return &data; }

    CompressedCoverage ccov;
    T data;
};

template<> struct CompressedCoverage_t<void> {

    CompressedCoverage_t(size_t sz = 0, bool full = false) : ccov(sz, full) {}
    CompressedCoverage_t(const CompressedCoverage& c) : ccov(c) {}
    CompressedCoverage_t(CompressedCoverage&& c) : ccov(move(c)) {}

    BFG_INLINE const void* getData() const { return nullptr; }
    BFG_INLINE void* getData() { return nullptr; }

    CompressedCoverage ccov;
};

#endif // BFG_COMPRESSED_COVERAGE_HPP
