#ifndef BIFROST_UNITIG_HPP
#define BIFROST_UNITIG_HPP

#include "Common.hpp"
#include "Kmer.hpp"
#include "CompressedSequence.hpp"
#include "CompressedCoverage.hpp"

/** @file src/Unitig.hpp
* The Unitig interface.
* Code snippets using these interface are provided in snippets/test.cpp.
*/

 /** @class Unitig
* @brief Represent a unitig which is a vertex of the Compacted de Bruijn graph.
* The first template argument T is the type of data associated with the unitigs
* of the graph.
* @var Unitig::data
* Data associated with the unitigs of the graph.
*/
template<typename T = void>
class Unitig {

    public:

        Unitig() {}
        Unitig(const CompressedSequence& s, const CompressedCoverage& c) : seq(s), cov(c) {}
        Unitig(CompressedSequence&& s, CompressedCoverage&& c) : seq(move(s)), cov(move(c)) {}
        Unitig(const char* s, bool full = false) : seq(s), cov(seq.size() - Kmer::k + 1, full) {}

        BFG_INLINE size_t numKmers() const {

            return seq.size( ) - Kmer::k + 1;
        }

        BFG_INLINE size_t length() const {

            return seq.size();
        }

        BFG_INLINE const CompressedSequence& getSeq() const {

            return seq;
        }

        BFG_INLINE CompressedSequence& getSeq() {

            return seq;
        }

        BFG_INLINE const CompressedCoverage& getCov() const {

            return cov;
        }

        BFG_INLINE CompressedCoverage& getCov() {

            return cov;
        }

        BFG_INLINE const T* getData() const {

            return &data;
        }

        BFG_INLINE T* getData() {

            return &data;
        }

    private:

        CompressedSequence seq;
        CompressedCoverage cov;

        T data;
};

///@cond NO_DOC
template<>
class Unitig<void> {

    public:

        Unitig() {}
        Unitig(const CompressedSequence& s, const CompressedCoverage& c) : seq(s), cov(c) {}
        Unitig(CompressedSequence&& s, CompressedCoverage&& c) : seq(move(s)), cov(move(c)) {}
        Unitig(const char* s, bool full = false) : seq(s), cov(seq.size() - Kmer::k + 1, full) {}

        BFG_INLINE size_t numKmers() const {

            return seq.size( ) - Kmer::k + 1;
        }

        BFG_INLINE size_t length() const {

            return seq.size();
        }

        BFG_INLINE const CompressedSequence& getSeq() const {

            return seq;
        }

        BFG_INLINE CompressedSequence& getSeq() {

            return seq;
        }

        BFG_INLINE const CompressedCoverage& getCov() const {

            return cov;
        }

        BFG_INLINE CompressedCoverage& getCov() {

            return cov;
        }

        BFG_INLINE const void* getData() const {

            return nullptr;
        }

        BFG_INLINE void* getData() {

            return nullptr;
        }

    private:

        CompressedSequence seq;
        CompressedCoverage cov;
};
///@endcond

#endif // BFG_CONTIG_HPP
