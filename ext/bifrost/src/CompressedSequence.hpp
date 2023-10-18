#ifndef BIFROST_COMPRESSED_SEQUENCE_HPP
#define BIFROST_COMPRESSED_SEQUENCE_HPP

#include <cstring>
#include <string>
#include <stdint.h>

#include "Kmer.hpp"
#include "wyhash.h"

/* Short description:
 *  - Compress a DNA string by using 2 bits per base instead of 8
 *  - Easily get the DNA string back from the compressed format
 *  - Create a sequence from a kmer
 *  - Get kmers from a sequence
 *  - Get length of a sequence
 *  - Easily get length of matching substring from a given string
 * */
class CompressedSequence {

    public:

        CompressedSequence();
        ~CompressedSequence();

        CompressedSequence(const char *s);
        CompressedSequence(const string& s);
        CompressedSequence(const Kmer& km);

        CompressedSequence(const CompressedSequence& o); // Copy constructor
        CompressedSequence(CompressedSequence&& o); // Move constructor

        CompressedSequence& operator=(const CompressedSequence& o); // Copy assignment
        CompressedSequence& operator=(CompressedSequence&& o); // Move assignment

        void clear();

        Kmer getKmer(size_t offset) const;
        Minimizer getMinimizer(size_t offset) const;

        //bool compareKmer(const size_t offset, const Kmer& km) const;
        bool compareKmer(const size_t offset, const size_t length, const Kmer& km) const;

        void setSequence(const CompressedSequence& o, const size_t offset_o, const size_t length_o, const size_t offset = 0, const bool reversed = false);
        void setSequence(const char *s, const size_t offset, const size_t length, const bool reversed = false);

        void toString(char *s, const size_t offset, const size_t length) const;
        string toString(const size_t offset, const size_t length) const;

        BFG_INLINE string toString() const {

            return toString(0, size());
        }

        BFG_INLINE void toString(char *s) const {

            toString(s, 0, size());
        }


        BFG_INLINE void setSequence(const string& s, const size_t offset, const size_t length, const bool reversed = false) {

            setSequence(s.c_str(), offset, length, reversed);
        }


        BFG_INLINE void setSequence(const Kmer& km, const size_t offset, const size_t length, const bool reversed = false) {

            char s[Kmer::MAX_K + 1];

            km.toString(s);
            setSequence(s, offset, length, reversed);
        }

        BFG_INLINE CompressedSequence rev() const {

            CompressedSequence r;

            r.setSequence(*this, 0, size(), 0, true);

            return r;
        }

        size_t jump(const char *s, const size_t i, int pos, const bool reversed) const;
        //size_t bw_jump(const char *s, const size_t i, int pos, const bool reversed) const;

        int64_t findKmer(const Kmer& km) const;

        bool write(std::ostream& stream_out) const;
        bool read(std::istream& stream_in);

        BFG_INLINE void reserveLength(const size_t new_length) {

           _resize_and_copy(round_to_bytes(new_length), size());
        }

        BFG_INLINE char operator[](const size_t offset) const {

            return alpha[(getPointer()[offset >> 2] >> ((offset & 0x3) << 1)) & 0x03];
        }

        BFG_INLINE char getChar(const size_t offset) const {

            return alpha[(getPointer()[offset >> 2] >> ((offset & 0x3) << 1)) & 0x03];
        }

        BFG_INLINE bool isShort() const {

            return ((asBits._size & 0x01) == 1);
        }

        BFG_INLINE size_t size() const {

            if (isShort()) return (asBits._size >> 1);

            return (asPointer._length >> 1);

            //const bool is_short = isShort();
            //return ((static_cast<size_t>(asBits._size) & (static_cast<size_t>(!is_short)-1)) + (static_cast<size_t>(asPointer._length) & (static_cast<size_t>(is_short)-1))) >> 1;
        }

        BFG_INLINE uint64_t hash(const uint64_t seed = 0) const {

            return wyhash(getPointer(), round_to_bytes(size()), seed, _wyp);
        }

    private:

        void _resize_and_copy(const size_t new_cap, const size_t copy_limit);

        BFG_INLINE size_t round_to_bytes(const size_t len) const {

            return (len+3)/4;
        }

        BFG_INLINE void initShort() {

            asBits._size = 1; // short and size 0

            memset(asBits._arr, 0, 15); // clear other bits
        }

        BFG_INLINE void setSize(const size_t size) {

            if (isShort()) asBits._size = ((0x7f & size) << 1) | 0x01; // set short flag
            else asPointer._length = (0x7fffffffffffffff & size) << 1;
        }

        BFG_INLINE const unsigned char* getPointer() const {

            if (isShort()) return &(asBits._arr[0]);

            return asPointer._data;
        }

        static const uint8_t revBits[256];

        union {

            struct {

                uint64_t _length; // size of sequence
                unsigned char *_data; // 0-based 2bit compressed dna string

            } asPointer;

            struct {

                uint8_t _size; // 7 bits can index up to 128
                unsigned char _arr[15]; // can store 124 nucleotides

            } asBits;
        };
};

#endif // BFG_COMPRESSED_SEQUENCE_HPP
