#ifndef BIFROST_UNITIGMAP_HPP
#define BIFROST_UNITIGMAP_HPP

#include <string>
#include "Common.hpp"
#include "Kmer.hpp"

/** @file src/UnitigMap.hpp
* UnitigMap type interface.
* Code snippets using this interface are provided in snippets/test.cpp.
*/

template<typename U> class Unitig;
template<typename U, typename G> class CompactedDBG;
template<typename U, typename G, bool is_const> class BackwardCDBG;
template<typename U, typename G, bool is_const> class ForwardCDBG;
template<typename U, typename G, bool is_const> class neighborIterator;

/** @class UnitigMapBase
* @brief Structure containing the basic information of a unitig mapping. This structure is independent
* from the graph. It is the base class for the class UnitigMap.
* @var UnitigMapBase::isEmpty
* True if there is no mapping.
* @var UnitigMapBase::dist
* Start position of the mapping (0-based distance) from the start of the reference unitig.
* @var UnitigMapBase::len
* Length of the mapping on the reference unitig, in k-mers.
* @var UnitigMapBase::size
* Length of the reference unitig.
* @var UnitigMapBase::strand
* True if the mapped k-mer or sequence matches the forward strand, false if it matches its reverse-complement.
*/
struct UnitigMapBase {

    /** UnitigMapBase constructor (isEmpty = true).
    * @param length is the length of the mapping in k-mers (default is 1 k-mer).
    * @return an empty UnitigMapBase.
    */
    UnitigMapBase(const size_t length = 1);

    /** UnitigMapBase constructor (isEmpty = false).
    * @param start is the start position of the mapping (0-based distance) from the start of the reference unitig.
    * @param length is the length of the mapping in k-mers.
    * @param unitig_sz is the length of the reference unitig used for the mapping.
    * @param strand indicates if the mapped k-mer or sequence matches the forward strand (true) or the reverse-complement (false).
    * @return a UnitigMapBase.
    */
    UnitigMapBase(const size_t start, const size_t length, const size_t unitig_sz, const bool strand);

    /** Equality operator: check if two UnitigMapBase are the same.
    * @return a boolean indicating if two UnitigMapBase are the same.
    */
    bool operator==(const UnitigMapBase& o) const;

    /** Inequality operator: check if two UnitigMapBase are different.
    * @return a boolean indicating if the two UnitigMapBase are different.
    */
    bool operator!=(const UnitigMapBase& o) const;

    size_t dist;
    size_t len;
    size_t size;

    bool strand;
    bool isEmpty;
};

/** @class UnitigMap
* @brief Contain all the information for the mapping of a k-mer or a sequence to a unitig
* of a Compacted de Bruijn graph. A UnitigMap object has 3 template parameters: the type of data
* associated with the unitigs of the graph, the type of data associated with the graph and a boolean
* indicating if this is a constant UnitigMap (const_UnitigMap) or not. A const_UnitigMap can be
* modified but you can't modify the CompactedDBG you can access using UnitigMap::getCompactedDBG.
* The unitig data and graph data types should be the same as the ones used for the CompactedDBG.
* \code{.cpp}
* UnitigMap<> um_1; // No unitig data, no graph data, NOT constant and its content IS NOT constant
* UnitigMap<void, void, false> um_2; // Equivalent to previous notation
* const UnitigMap<> um_3; // No unitig data, no graph data, constant BUT its content IS NOT constant
* UnitigMap<void, void, true> um_4; // No unitig data, no graph data, NOT constant BUT its content IS constant
* const UnitigMap<void, void, true> um_5; // No unitig data, no graph data, constant AND its content IS constant
* UnitigMap<myUnitigData, myGraphData> um_6; // Unitig data of type myUnitigData for each unitig, graph data of type myGraphData, not constant

* CompactedDBG<>* cdbg_ptr_1 = um_1.getGraph(); // Associated CompactedDBG can be modified from the UnitigMap
* CompactedDBG<>* cdbg_ptr_2 = um_2.getGraph(); // Associated CompactedDBG can be modified from the UnitigMap
* CompactedDBG<>* cdbg_ptr_3 = um_3.getGraph(); // Associated CompactedDBG can be modified from the UnitigMap
* const CompactedDBG<>* cdbg_ptr_4 = um_4.getGraph(); // Associated CompactedDBG cannot be modified from the UnitigMap
* const CompactedDBG<>* cdbg_ptr_5 = um_5.getGraph(); // Associated CompactedDBG cannot be modified from the UnitigMap
* CompactedDBG<myUnitigData, myGraphData>* cdbg_ptr_6 = um_6.getGraph(); // Associated CompactedDBG can be modified from the UnitigMap
* \endcode
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class UnitigMap : public UnitigMapBase {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    template<typename U, typename G> friend class CompactedDBG;
    template<typename U, typename G, bool C> friend class BackwardCDBG;
    template<typename U, typename G, bool C> friend class ForwardCDBG;
    template<typename U, typename G, bool C> friend class unitigIterator;
    template<typename U, typename G, bool C> friend class UnitigMap;

    typedef typename std::conditional<is_const, const CompactedDBG<U, G>*, CompactedDBG<U, G>*>::type CompactedDBG_ptr_t;
    typedef typename std::conditional<is_const, const U*, U*>::type Unitig_data_ptr_t;

    public:

        typedef BackwardCDBG<U, G, is_const> UnitigMap_BW;
        typedef ForwardCDBG<U, G, is_const> UnitigMap_FW;

        /** UnitigMap constructor.
        * @param length is the length of the mapping in k-mers (default is 1 k-mer).
        * @param cdbg_ is a pointer to the CompactedDBG containing the reference unitig used for the mapping (default is nullptr).
        * @return an empty UnitigMap.
        */
        UnitigMap(size_t length = 1, CompactedDBG_ptr_t cdbg_ = nullptr);

        /** UnitigMap constructor.
        * This constructor is used to create temporary mappings and must not be used to extract information from the graph.
        * @param start is the start position of the mapping (0-based distance) from the start of the reference unitig.
        * @param length is the length of the mapping in k-mers.
        * @param unitig_sz is the length of the reference unitig used for the mapping.
        * @param strand indicates if the mapped k-mer or sequence matches the forward strand (true) or the reverse-complement (false).
        * @return a UnitigMap.
        */
        UnitigMap(const size_t start, const size_t length, const size_t unitig_sz, const bool strand);

        /** Equality operator: check if two UnitigMap are the same.
        * @return a boolean indicating if two UnitigMap are the same.
        */
        bool operator==(const UnitigMap& o) const;

        /** Inequality operator: check if two UnitigMap are different.
        * @return a boolean indicating if the two UnitigMap are different.
        */
        bool operator!=(const UnitigMap& o) const;

        /** check if two UnitigMap have the same unitig as reference.
        * @return a boolean indicating if two UnitigMap have the same unitig as reference.
        */
        bool isSameReferenceUnitig(const UnitigMap& o) const;

        /** Create a string containing the sequence corresponding to the mapping.
        * @return a string containing the sequence corresponding to the mapping or
        * an empty string if there is no mapping (UnitigMap::isEmpty = true).
        */
        std::string mappedSequenceToString() const;

        /** Create a string containing the sequence of the reference unitig used the mapping.
        * @return a string containing the sequence of the reference unitig used the mapping or
        * an empty string if there is no mapping (UnitigMap::isEmpty = true).
        */
        std::string referenceUnitigToString() const;

        /** Compute the length of the longest common prefix between a given sequence and
        * the reference unitig used in the mapping.
        * @param s is a pointer to an array of characters representing the sequence from
        * which the length of the longest common prefix must be computed.
        * @param pos_s is the start position in s from which the longest common prefix must
        * be computed.
        * @param pos_um_seq is the start position in the reference unitig of the mapping from
        * which the longest common prefix must be computed.
        * @param um_reversed indicates if the longest common prefix must be computed from
        * the reverse-complement of the reference unitig used in the mapping (true) or not (false).
        * @return the length of the longest common prefix
        */
        size_t lcp(const char* s, const size_t pos_s = 0, const size_t pos_um_seq = 0, const bool um_reversed = false) const;

        /** Get the head k-mer of the reference unitig used for the mapping.
        * @return a Kmer object which is either the head k-mer of the mapped unitig or
        * an empty k-mer if there is no mapping (UnitigMap::isEmpty = true).
        */
        Kmer getUnitigHead() const;

        /** Get the tail k-mer of the reference unitig used for the mapping.
        * @return a Kmer object which is either the tail k-mer of the mapped unitig or
        * an empty k-mer if there is no mapping (UnitigMap::isEmpty = true).
        */
        Kmer getUnitigTail() const;

        /** Get the k-mer starting at position pos in the reference unitig used for the mapping.
        * @param pos is the start position of the k-mer to extract.
        * @return a Kmer object which is either the k-mer starting at position pos in the mapped unitig or
        * an empty k-mer if there is either no mapping (UnitigMap::isEmpty = true) or the position is
        * greater than length(unitig)-k.
        */
        Kmer getUnitigKmer(const size_t pos) const;

        /** Get the head k-mer of the mapped sequence.
        * @return a Kmer object which is either the head k-mer of the mapped sequence or
        * an empty k-mer if there is no mapping (UnitigMap::isEmpty = true).
        */
        Kmer getMappedHead() const;

        /** Get the tail k-mer of the mapped sequence.
        * @return a Kmer object which is either the tail k-mer of the mapped sequence or
        * an empty k-mer if there is no mapping (UnitigMap::isEmpty = true).
        */
        Kmer getMappedTail() const;

        /** Get the k-mer starting at position pos in the mapped sequence.
        * @param pos is the start position of the k-mer to extract within the mapped sequence.
        * @return a Kmer object which is either the k-mer starting at position pos in the mapped sequence or
        * an empty k-mer if there is either no mapping (UnitigMap::isEmpty = true) or the position is
        * greater than length(unitig)-k.
        */
        Kmer getMappedKmer(const size_t pos) const;

        /** Create a new UnitigMap object which is the mapping of a k-mer on a reference unitig
        * @param pos is the start position of the k-mer to map in the reference unitig used for the current mapping.
        * @return a UnitigMap object which is the mapping.
        */
        UnitigMap<U, G, is_const> getKmerMapping(const size_t pos) const;

        /** Get a pointer to the data associated with the reference unitig used in the mapping.
        * @return a pointer to the data associated with the reference unitig used in the mapping
        * or a nullptr pointer if:
        * - there is no mapping (UnitigMap::isEmpty = true)
        * - there is no data associated with the unitigs (U = void).
        * The returned pointer is constant if the UnitigMap is constant (UnitigMap<U, G, true>).
        */
        Unitig_data_ptr_t getData() const;

        /** Create a UnitigMap_BW object that can create iterators (through UnitigMap_BW::begin() and
        * UnitigMap_BW::end()) over the predecessors of the reference unitig used in the mapping.
        * @return a UnitigMap_BW object that can create iterators (through UnitigMap_BW::begin() and
        * UnitigMap_BW::end()) over the predecessors of the reference unitig used in the mapping.
        */
        UnitigMap_BW getPredecessors() const;

        /** Create a UnitigMap_FW object that can create iterators (through UnitigMap_FW::begin() and
        * UnitigMap_FW::end()) over the successors of the reference unitig used in the mapping.
        * @return a UnitigMap_FW object that can create iterators (through UnitigMap_FW::begin() and
        * UnitigMap_FW::end()) over the successors of the reference unitig used in the mapping.
        */
        UnitigMap_FW getSuccessors() const;

        /** Get a pointer to the CompactedDBG containing the reference unitig used in the mapping.
        * @return a pointer to the CompactedDBG containing the reference unitig used in the mapping.
        * If the mapping is empty, a nullptr pointer is returned. The pointer is a constant pointer
        * if the UnitigMap is constant (UnitigMap<U, G, true>).
        */
        inline CompactedDBG_ptr_t getGraph() const { return cdbg; }

        operator UnitigMap<U, G, true>() const {

            UnitigMap<U, G, true> um(pos_unitig, dist, len, size, isShort, isAbundant, strand, cdbg);

            um.isEmpty = isEmpty;

            return um;
        }

        void setFullCoverage() const;
        void increaseCoverage() const;
        void decreaseCoverage() const;

        bool isCoverageFull() const;
        size_t getCoverage(const size_t pos) const;

    private:

        UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd, CompactedDBG_ptr_t cdbg_);

        neighborIterator<U, G, is_const> bw_begin() const;
        neighborIterator<U, G, is_const> bw_end() const;

        neighborIterator<U, G, is_const> fw_begin() const;
        neighborIterator<U, G, is_const> fw_end() const;

        template<bool is_void> typename std::enable_if<!is_void, Unitig<U>>::type splitData_(const bool last_split) const;
        template<bool is_void> typename std::enable_if<is_void, Unitig<U>>::type splitData_(const bool last_split) const;

        Unitig<U> splitData(const bool last_split) const;

        template<bool is_void> typename std::enable_if<!is_void, Unitig_data_ptr_t>::type getData_() const;
        template<bool is_void> typename std::enable_if<is_void, Unitig_data_ptr_t>::type getData_() const;

        void partialCopy(const UnitigMap<U, G, is_const>& um);

        size_t pos_unitig; // unitig pos. in v_unitigs or km_unitigs or h_kmers

        bool isShort; // true if the unitig has length k
        bool isAbundant; // true if the unitig has length k and has an abundant minimizer

        CompactedDBG_ptr_t cdbg;
};

template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
struct UnitigMapHash {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    size_t operator()(const UnitigMap<U, G, is_const>& um) const {

        struct UnitigMapTMP {

            size_t pos_unitig; // unitig pos. in v_unitigs or km_unitigs or h_kmers
            size_t dist;
            size_t len;
            size_t size;

            bool strand;
            bool isEmpty;

            bool isShort; // true if the unitig has length k
            bool isAbundant; // true if the unitig has length k and has an abundant minimizer

            const void* cdbg;

            UnitigMapTMP(const UnitigMap<U, G, is_const>& um) : pos_unitig(um.pos_unitig), dist(um.dist), len(um.len), size(um.size),
                                                                strand(um.strand), isEmpty(um.isEmpty), isShort(um.isShort),
                                                                isAbundant(um.isAbundant), cdbg(static_cast<const void*>(um.cdbg)) {};
        };

        UnitigMapTMP tmp(um);

        //return static_cast<size_t>(XXH64(static_cast<const void*>(&tmp), sizeof(UnitigMapTMP), 0));
        return static_cast<size_t>(wyhash(&tmp, sizeof(UnitigMapTMP), 0, _wyp));
    }
};

#include "UnitigMap.tcc"

#endif
