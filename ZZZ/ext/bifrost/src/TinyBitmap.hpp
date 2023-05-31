#ifndef BIFROST_TINYBITMAP_HPP
#define BIFROST_TINYBITMAP_HPP

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <vector>

/* TinyBitmap is a compressed bitmap that mimics the behavior of a CRoaring container.
*  Its main purpose is to store a tiny set of unsigned integers, up to 65488 uint.
*  Main differences with a CRoaring bitmap are:
*  + For one value inserted, CRoaring allocates >100 bytes, TinyBitmap allocates 24.
*  + CRoaring containers have a variable size except in bitmap mode (fixed 8kB). TinyBitmap
*    container have variable size in all modes (bitmap, list, RLE).
*  + Accessing a value in a CRoaring container is >3 cache-miss. TinyBitmap is 1 cache-miss.
*  - CRoaring can store >65488 values.
*  - CRoaring is SIMD optimized.
*  - CRoaring has a lot more functions (set intersection, union, etc.).
*/

using namespace std;

class TinyBitmap {

    class TinyBitmapIterator : public std::iterator<std::input_iterator_tag, uint32_t> {

        friend class TinyBitmap;

        public:

            TinyBitmapIterator();

            TinyBitmapIterator& operator++();
            TinyBitmapIterator operator++(int);

            inline bool operator==(const TinyBitmapIterator& o) const {

                if (invalid || o.invalid) return invalid && o.invalid;

                return  (tiny_bmp == o.tiny_bmp) && (sz == o.sz) && (mode == o.mode) && (card == o.card) &&
                        (i == o.i) && (j == o.j) && (e == o.e) && (offset == o.offset) && (val == o.val);
            }

            inline bool operator!=(const TinyBitmapIterator& o) const { return !operator==(o); }

            inline const uint32_t& operator*() const { return val; }
            inline const uint32_t* operator->() const { return &val; }

        private:

            TinyBitmapIterator(const TinyBitmap& t_bmp_, const bool start);

            uint16_t sz;
            uint16_t mode;
            uint16_t card;

            uint16_t i;
            uint16_t j;
            uint16_t e;

            uint32_t offset;
            uint32_t val;

            bool invalid;

            const uint16_t* tiny_bmp;
    };

    friend class TinyBitmapIterator;

    public:

        typedef TinyBitmapIterator const_iterator;

        TinyBitmap();
        TinyBitmap(const TinyBitmap& o);
        TinyBitmap(TinyBitmap&& o);
        TinyBitmap(uint16_t** o_ptr);

        ~TinyBitmap();

        TinyBitmap& operator=(const TinyBitmap& o);
        TinyBitmap& operator=(TinyBitmap&& o);
        TinyBitmap& operator=(uint16_t** o_ptr);

        void clear();

        bool add(const uint32_t val);
        bool remove(const uint32_t val);

        bool contains(const uint32_t val) const;
        bool containsRange(const uint32_t val_start, const uint32_t val_end) const;

        uint32_t maximum() const;

        bool write(ostream& stream_out) const;
        bool read(istream& stream_in);

        inline void toArray(uint32_t* values) const {

            for (const auto v : *this) *(values++) = v;
        }

        size_t getSizeInBytes() const;
        size_t size() const;
        size_t size(uint32_t start_value, const uint32_t end_value) const;

        size_t runOptimize();
        size_t shrinkSize();

        const_iterator begin() const;
        const_iterator end() const;

        static bool test(const bool verbose = true);

        inline uint16_t* detach() {

            uint16_t* ret = tiny_bmp;
            tiny_bmp = nullptr;

            return ret;
        }

    private:

        void print() const;

        bool change_sz(const uint16_t sz_min);
        bool switch_mode(const uint16_t sz_min, const uint16_t new_mode);

        inline uint16_t getSize() const { return (tiny_bmp[0] & sz_mask) >> 3; }
        inline uint16_t getMode() const { return (tiny_bmp[0] & mode_mask); }
        inline uint16_t getBits() const { return (tiny_bmp[0] & bits_mask); }

        inline uint16_t getCardinality() const { return tiny_bmp[1]; }
        inline uint16_t getOffset() const { return tiny_bmp[2]; }

        static inline uint16_t getNextSize(const uint16_t sz) {

            uint16_t idx = 0;

            while (sizes[idx] < sz) ++idx;

            return sizes[idx];
        }

        static inline uint16_t rndup(uint16_t v) {

            v--;
            v |= v >> 1;
            v |= v >> 2;
            v |= v >> 4;
            v |= v >> 8;
            v++;

            return v;
        }

        static const uint16_t sz_mask;
        static const uint16_t mode_mask;
        static const uint16_t bits_mask;

        static const uint16_t bmp_mode;
        static const uint16_t list_mode;
        static const uint16_t rle_list_mode;

        static const uint16_t bits_16;
        static const uint16_t bits_32;

        static const uint16_t sizes[];
        static const uint16_t nb_sizes;

        uint16_t* tiny_bmp;
};

#endif
