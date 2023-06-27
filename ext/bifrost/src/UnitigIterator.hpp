#ifndef BIFROST_UNITIG_ITERATOR_HPP
#define BIFROST_UNITIG_ITERATOR_HPP

#include "UnitigMap.hpp"
#include "KmerHashTable.hpp"
#include "CompressedCoverage.hpp"

/** @file src/UnitigIterator.hpp
* The unitigIterator type interface.
* Code snippets using this interface are provided in snippets/test.cpp.
*/

template<typename U, typename G> class CompactedDBG;

/** @class unitigIterator
* @brief Iterator for the unitigs of a Compacted de Bruijn graph.
* A unitigIterator object has 3 template parameters: the type of data associated with the unitigs
* of the graph, the type of data associated with the graph and a boolean indicating if this is a
* constant iterator or not. Note that you are supposed to use this class as the iterator of the
* class CompactedDBG (CompactedDBG::iterator and CompactedDBG::const_iterator) so you shouldn't
* have to instantiate an object unitigIterator and its template parameters yourself.
* The unitig data and graph data types should be the same as the ones used for the CompactedDBG
* the iterator is from. No specific order (such as a lexicographic one) is assumed during iteration.
* \code{.cpp}
* CompactedDBG<> cdbg;
* ... // Some more code, cdbg construction
* for (const auto& unitig : cdbg){
*   cout << unitig.toString() << endl; // unitig is of type UnitigMap
* }
* for (CompactedDBG<>::const_iterator it = cdbg.begin(); it != cdbg.end(); ++it){
*   cout << it.toString() << endl;
* }
* \endcode
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class unitigIterator : public std::iterator<std::input_iterator_tag, UnitigMap<Unitig_data_t, Graph_data_t, is_const>, int> {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        typedef typename std::conditional<is_const, const CompactedDBG<U, G>*, CompactedDBG<U, G>*>::type CompactedDBG_ptr_t;

        /** Constructor.
        * @return an empty unitigIterator.
        */
        unitigIterator();

        /** Constructor.
        * @param cdbg_ is a pointer to a Compacted de Bruijn graph. The unitigIterator will iterate over the unitigs of this graph.
        * If the unitigIterator is a constant, the CompactedDBG pointer is a constant pointer and the accessed unitigs can't be modified.
        * @return a unitigIterator.
        */
        unitigIterator(CompactedDBG_ptr_t cdbg_);

        /** Copy constructor.
        * @return a copy of a unitigIterator.
        */
        unitigIterator(const unitigIterator& o);

        /** Prefix increment, iterate over the next unitig.
        * Note that no specific order (such as a lexicographic one) is assumed during iteration.
        */
        unitigIterator& operator++();

        /** Postfix increment, iterate over the next unitig.
        * Note that no specific order (such as a lexicographic one) is assumed during iteration.
        */
        unitigIterator operator++(int);

        /** Equality operator: check if two unitigIterator are the same.
        * @param o is another unitigIterator.
        * @return a boolean indicating whether the two unitigIterator are the same.
        */
        bool operator==(const unitigIterator& o) const;

        /** Inequality operator: check if two unitigIterator are different.
        * @param o is another unitigIterator.
        * @return a boolean indicating whether the two unitigIterator are different.
        */
        bool operator!=(const unitigIterator& o) const;

        /** Indirection operator.
        * @return a UnitigMap object reference which contains information about the mapped unitig.
        */
        const UnitigMap<U, G, is_const>& operator*() const;

        /** Dereference operator.
        * @return a UnitigMap object pointer which contains information about the mapped unitig.
        */
        const UnitigMap<U, G, is_const>* operator->() const;

    private:

        size_t i;

        size_t v_unitigs_sz;
        size_t v_kmers_sz;
        size_t h_kmers_ccov_sz;
        size_t sz;

        bool invalid;

        typename KmerHashTable<CompressedCoverage_t<U>>::const_iterator it_h_kmers_ccov;

        UnitigMap<U, G, is_const> um;

        CompactedDBG_ptr_t cdbg;
};

#include "UnitigIterator.tcc"

#endif
