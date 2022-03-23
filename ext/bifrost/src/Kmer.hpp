#ifndef BIFROST_KMER_HPP
#define BIFROST_KMER_HPP

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 32
#endif

#ifndef MAX_GMER_SIZE
#define MAX_GMER_SIZE MAX_KMER_SIZE
#endif

#include <stdint.h>
#include <stdio.h>

#include <bitset>
#include <cassert>
#include <cstring>
#include <iostream>
#include <string>

#include "Common.hpp"

/** @file src/Kmer.hpp
* Interface for the class Kmer:
* - Store k-mer strings by using 2 bits per base
* - Easily return reverse-complements of k-mers, e.g. TTGG -> CCAA
* - Easily compare k-mers
* - Provide hash of k-mers
* - Get last and next kmer, e.g. ACGT -> CGTT or ACGT -> AACGT
*/

class CompressedSequence;

/** @class Kmer
* @brief Interface to store and manipulate k-mers.
* The default maximum supported k-mer size is 31. To work with larger k, you must replace the macro defined
* MAX_KMER_SIZE with a larger number which is a multiple of 32. For such a number K, the available  maximum
* supported k-mer size will be K-1. For example, if MAX_KMER_SIZE is defined as 64, the maximum k allowed is 63.
* Keep in mind that increasing MAX_KMER_SIZE increases memory usage (k=31 uses 8 bytes of memory per k-mer
* while k=63 uses 16 bytes of memory per k-mer).
*/
class Kmer {

    friend class CompressedSequence;

    public:

        /** Constructor (initialize a k-mer with 'A' k times).
        */
        Kmer();

        /** Copy constructor (copy a k-mer). After the call to this function, the same k-mer exists twice in
        * memory.
        * @param o is a constant reference to the k-mer to copy
        */
        Kmer(const Kmer& o);

        /** Constructor. Initialize a k-mer from a string of characters.
        * @param s is a constant pointer to an array of characters. The array of characters must contain valid DNA
        * characters ('A', 'C', 'G', 'T') and be of length at least k.
        */
        explicit Kmer(const char* s);

        /** Copy assignment operator (copy a k-mer). After the call to this function, the same k-mer exists
        * twice in memory.
        * @param o is a constant reference to the k-mer to copy
        */
        Kmer& operator=(const Kmer& o);

        /** Set a k-mer as "deleted".
        */
        BFG_INLINE void set_deleted() {

            for (size_t i = 0; i < MAX_K/32; ++i) longs[i] = 0xffffffffffffffffULL;
        }

        /** Set a k-mer as "empty".
        */
        BFG_INLINE void set_empty() {

            for (size_t i = 0; i < MAX_K/32; ++i) longs[i] = 0xfffffffffffffffeULL;
        }

        /** Check whether a k-mer is "deleted"
        * @return a boolean indicating if the k-mer is "deleted" (true) or not (false).
        */
        BFG_INLINE bool isDeleted() const {

            return (longs[(MAX_K/32)-1] == 0xffffffffffffffffULL);
        }

        /** Check whether a k-mer is "empty"
        * @return a boolean indicating if the k-mer is "empty" (true) or not (false).
        */
        BFG_INLINE bool isEmpty() const {

            return (longs[(MAX_K/32)-1] == 0xfffffffffffffffeULL);
        }

        bool operator<(const Kmer& o) const;

        /** Equality comparison operator.
        * @param o is a constant reference to the k-mer to compare to.
        * @return a boolean indicating if two k-mers are the same (true), including their "empty" and "deleted"
        * flags, or not (false).
        */
        bool operator==(const Kmer& o) const;

        /** Inequality comparison operator.
        * @param o is a constant reference to the k-mer to compare to.
        * @return a boolean indicating if two k-mers are the different (true), including their "empty" and
        * "deleted" flags, or not (false).
        */
        bool operator!=(const Kmer& o) const;

        /** Get the hash of a k-mer.
        * @param seed is a seed number for the hash function (0 by default).
        * @return a hash of the k-mer.
        */
        BFG_INLINE uint64_t hash(const uint64_t seed = 0) const {

            return wyhash(bytes, MAX_K/4, seed, _wyp);
        }

        /** Get the reverse-complement of a k-mer.
        * @return a new k-mer which is the reverse-complement.
        */
        Kmer twin() const;

        /** Get the canonical k-mer (lexicographically smallest between a k-mer and its reverse-complement).
        * @return a new k-mer which is the canonical k-mer.
        */
        Kmer rep() const;

        /** Get a new k-mer which is the shift of the current k-mer of one base on the left with one
        * new character on the right. For a k-mer km of length k, km.forwardBase(b) = km[1..k-1] + b.
        * @param b is a new character to add on the right of the shifted k-mer (as a last character).
        * It must be either 'A', 'C', 'G' or 'T'.
        * @return a new k-mer which is the shifted k-mer.
        */
        Kmer forwardBase(const char b) const;

        /** Get a new k-mer which is the shift of the current k-mer of one base on the right with one
        * new character on the left. For a k-mer km of length k, km.backwardBase(b) = b + km[0..k-2].
        * @param b is a new character to add on the left of the shifted k-mer (as a first character).
        * It must be either 'A', 'C', 'G' or 'T'.
        * @return a new k-mer which is the shifted k-mer.
        */
        Kmer backwardBase(const char b) const;

        /** Shift the current k-mer of one base on the left with one new character on the right.
        * @param b is a new character to add on the right (as a last character) after shifting the current
        * k-mer. It must be either 'A', 'C', 'G' or 'T'.
        */
        void selfForwardBase(const char b);

        /** Get the character at a given position in a k-mer.
        * @param offset is the position of the character to get in the k-mer.
        * @return the character to get in the k-mer
        */
        char getChar(const size_t offset) const;

       /** Set a character at a given position in a k-mer.
        * @param offset is the position of the character to set in the k-mer.
        * @param b is the character to set. It must be either 'A', 'C', 'G' or 'T'.
        * @return a boolean indicating if the character was succesfully set. The boolean will be equal to
        * false if offset > k-mer length or if the character to set is not 'A', 'C', 'G' or 'T'.
        */
        bool setChar(const size_t offset, const char b);

        /** Get the string of a k-mer.
        * @param s is a pointer to an array of characters that will be set with the string of a k-mer.
        */
        void toString(char *s) const;

        /** Get the string of a k-mer.
        * @return the string of a k-mer.
        */
        std::string toString() const;

        /** Write a k-mer (binary) to a stream. The stream must be opened prior to this function call
        * and it won't be closed by this function.
        * @param stream_out is an out stream to which the k-mer must be written to. The stream must
        * be opened prior to this function call and it won't be closed by this function.
        * @return a boolean indicating if the writing was successful (true) or not.
        */
        bool write(std::ostream& stream_out) const;

        /** Read a k-mer (binary) from a stream. The stream must be opened prior to this function call
        * and it won't be closed by this function.
        * @param stream_in is an in stream from which the k-mer must be read. The stream must
        * be opened prior to this function call and it won't be closed by this function.
        * @return a boolean indicating if the reading was successful (true) or not.
        */
        bool read(std::istream& stream_in);

        /** Set the length of k-mers. After this function call, all Kmer functions considers
        * _k to be the k-mer length.
        * @param _k is the k-mer length to set.
        */
        static void set_k(const unsigned int _k);

        static unsigned int k;

    private:

        static const unsigned int MAX_K = MAX_KMER_SIZE;

        union {

            uint8_t bytes[MAX_K/4];
            uint64_t longs[MAX_K/32];
        };

        void set_kmer(const char *s);

        Kmer getLink(const size_t index) const;

        std::string getBinary() const;
};


struct KmerHash {

    inline size_t operator()(const Kmer& km) const {

        return km.hash();
    }
};

///@cond NO_DOC
class Minimizer {

    public:

        Minimizer();
        Minimizer(const Minimizer& o);
        explicit Minimizer(const char *s);

        Minimizer& operator=(const Minimizer& o);

        bool operator<(const Minimizer& o) const;
        bool operator==(const Minimizer& o) const;
        bool operator!=(const Minimizer& o) const;

        void set_minimizer(const char *s);

        BFG_INLINE uint64_t hash(const uint64_t seed = 0) const {

            return wyhash(bytes, MAX_G/4, seed, _wyp);
        }

        BFG_INLINE void set_deleted() {

            for (size_t i = 0; i < MAX_G/32; ++i) longs[i] = 0xffffffffffffffffULL;
        }

        BFG_INLINE void set_empty() {

            for (size_t i = 0; i < MAX_G/32; ++i) longs[i] = 0xfffffffffffffffeULL;
        }

        BFG_INLINE bool isDeleted() const {

            return (longs[(MAX_G/32)-1] == 0xffffffffffffffffULL);
        }

        BFG_INLINE bool isEmpty() const {

            return (longs[(MAX_G/32)-1] == 0xfffffffffffffffeULL);
        }

        Minimizer twin() const;
        Minimizer rep() const;

        Minimizer getLink(const size_t index) const;

        Minimizer forwardBase(const char b) const;
        Minimizer backwardBase(const char b) const;

        std::string getBinary() const;

        void toString(char *s) const;
        std::string toString() const;

        // static functions
        static void set_g(unsigned int _g);

        static unsigned int g;

    private:

        static const unsigned int MAX_G = MAX_GMER_SIZE;

        // data fields
        union {

            uint8_t bytes[MAX_G/4];
            uint64_t longs[MAX_G/32];
        };

        // By default MAX_K == 64 so the union uses 16 bytes
        // However sizeof(Kmer) == 24
        // Are the 8 extra bytes alignment?
};


struct MinimizerHash {

    inline size_t operator()(const Minimizer& minz) const {

        return minz.hash();
    }
};
///@endcond

#endif // BFG_KMER_HPP
