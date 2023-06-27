#ifndef BIFROST_NEIGHBOR_ITERATOR_HPP
#define BIFROST_NEIGHBOR_ITERATOR_HPP

#include "Kmer.hpp"

/** @file src/NeighborIterator.hpp
* The neighborIterator, BackwardCDBG and ForwardCDBG type interfaces.
* Code snippets using these interfaces are provided in snippets/test.cpp.
*/

template<typename U, typename G> class CompactedDBG;
template<typename U, typename G, bool is_const> class UnitigMap;

/** @class neighborIterator
* @brief Iterator for the neighbors (predecessors or successors) of a reference unitig used in a UnitigMap object.
* A neighborIterator object has 3 template parameters: the type of data associated with the unitigs of the graph,
* the type of data associated with the graph and a boolean indicating if this is a constant iterator or not.
* Note that you are supposed to use this class as the iterator of a BackwardCDBG or ForwardCDBG object, which can
* be obtained respectively from UnitigMap::getPredecessors() and UnitigMap::getSuccessors(), so you shouldn't
* have to instantiate an object neighborIterator and its template parameters yourself. The unitig data and graph
* data types should be the same as the ones used for the CompactedDBG the iterator is from. No specific order
* (such as a lexicographic one) is assumed during iteration.
* \code{.cpp}
* CompactedDBG<> cdbg;
* ... // Some more code, cdbg construction
* for (const auto& unitig : cdbg){
*   cout << unitig.toString() << endl; // unitig is of type UnitigMap
*   for (const auto& pred : unitig.getPredecessors()) cout << pred.toString() << endl;
*   for (const auto& succ : unitig.getSuccessors()) cout << succ.toString() << endl;
* }
* \endcode
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class neighborIterator : public std::iterator<std::input_iterator_tag, UnitigMap<Unitig_data_t, Graph_data_t, is_const>, int> {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        typedef typename std::conditional<is_const, const CompactedDBG<U, G>*, CompactedDBG<U, G>*>::type CompactedDBG_ptr_t;

        /** Constructor.
        * @return an empty neighborIterator.
        */
        neighborIterator();

        /** Constructor.
        * @param um_ is a UnitigMap object: the constructed neighborIterator will iterate over the neighbors of the UnitigMap reference unitig
        * @param is_forward_ indicates if the iterator must iterate over the successors (true) or the predecessors (false)
        * @return a neighborIterator.
        */
        neighborIterator(const UnitigMap<U, G, is_const>& um_, const bool is_forward_);

        /** Copy constructor.
        * @return a copy of a neighborIterator.
        */
        neighborIterator(const neighborIterator& o);

        /** Prefix increment, iterate over the next neighbor (predecessor or successor).
        * Note that no specific order (such as a lexicographic one) is assumed during iteration.
        */
        neighborIterator& operator++();

        /** Postfix increment, iterate over the next neighbor (predecessor or successor).
        * Note that no specific order (such as a lexicographic one) is assumed during iteration.
        */
        neighborIterator operator++(int);

        /** Equality operator: check if two neighborIterator are the same.
        * @param o is another neighborIterator.
        * @return a boolean indicating whether the two neighborIterator are the same.
        */
        bool operator==(const neighborIterator& o) const;

        /** Inequality operator: check if two neighborIterator are different.
        * @param o is another neighborIterator.
        * @return a boolean indicating whether the two neighborIterator are different.
        */
        bool operator!=(const neighborIterator& o) const;

        /** Indirection operator
        * @return a UnitigMap reference which contains information about the mapped neighbor.
        */
        const UnitigMap<U, G, is_const>& operator*() const;

        /** Dereference operator.
        * @return a UnitigMap pointer which contains information about the mapped neighbor.
        */
        const UnitigMap<U, G, is_const>* operator->() const;

    private:

        int i;

        bool is_fw;

        Kmer km_head;
        Kmer km_tail;

        UnitigMap<U, G, is_const> um;

        CompactedDBG_ptr_t cdbg;
};

/** @class BackwardCDBG
* @brief Wrapper for class neighborIterator to iterate over the predecessors of a reference unitig used in a UnitigMap object.
* A BackwardCDBG object has 3 template parameters: the type of data associated with the unitigs of the graph,
* the type of data associated with the graph and a boolean indicating if this is a constant neighborIterator or not.
* Note that you are supposed to obtain an instance of this class from UnitigMap::getPredecessors() so you shouldn't
* have to instantiate an object BackwardCDBG and its template parameters yourself. The unitig data and graph
* data types should be the same as the ones used for the CompactedDBG the iterator is from. No specific order
* (such as a lexicographic one) is assumed during iteration.
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class BackwardCDBG {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        /** Constructor.
        * @param um_ is a UnitigMap object: the wrapper will provide an iterator over the predecessors of the UnitigMap reference unitig
        */
        explicit BackwardCDBG(const UnitigMap<U, G, is_const>& um_);

        /** Check if the unitig has at least one predecessor.
        * @return a boolean indicating if the unitig has at least one predecessor (true) or not (false).
        */
        bool hasPredecessors() const;

        size_t cardinality() const;

        /** Return an iterator over the predecessors of a reference unitig. The returned iterator is initialized
        * over the first such predecessor, if there is one.
        * @return an iterator over the predecessors of a reference unitig. It is initialized over the first such
        * predecessor, if there is one.
        */
        neighborIterator<U, G, is_const> begin() const;

        /** Return an iterator over the past-the-last predecessor of a reference unitig.
        * @return an iterator over the past-the-last predecessor of a reference unitig.
        */
        neighborIterator<U, G, is_const> end() const;

    private:

        UnitigMap<U, G, is_const> um;
};

/** @class ForwardCDBG
* @brief Wrapper for class neighborIterator to iterate over the predecessors of a reference unitig used in a UnitigMap object.
* A ForwardCDBG object has 3 template parameters: the type of data associated with the unitigs of the graph,
* the type of data associated with the graph and a boolean indicating if this is a constant neighborIterator or not.
* Note that you are supposed to obtain an instance of this class from UnitigMap::getSuccessors() so you shouldn't
* have to instantiate an object ForwardCDBG and its template parameters yourself. The unitig data and graph
* data types should be the same as the ones used for the CompactedDBG the iterator is from. No specific order
* (such as a lexicographic one) is assumed during iteration.
*/
template<typename Unitig_data_t = void, typename Graph_data_t = void, bool is_const = false>
class ForwardCDBG {

    typedef Unitig_data_t U;
    typedef Graph_data_t G;

    public:

        /** Constructor.
        * @param um_ is a UnitigMap object: the wrapper will provide an iterator over the successors of the UnitigMap reference unitig
        */
        explicit ForwardCDBG(const UnitigMap<U, G, is_const>& um_);

        /** Check if the unitig has at least one successor.
        * @return a boolean indicating if the unitig has at least one successor (true) or not (false).
        */
        bool hasSuccessors() const;

        size_t cardinality() const;

        /** Return an iterator over the successors of a reference unitig. The returned iterator is initialized
        * over the first such successor, if there is one.
        * @return an iterator over the successors of a reference unitig. It is initialized over the first such
        * successor, if there is one.
        */
        neighborIterator<U, G, is_const> begin() const;

        /** Return an iterator over the past-the-last successor of a reference unitig.
        * @return an iterator over the past-the-last successor of a reference unitig.
        */
        neighborIterator<U, G, is_const> end() const;

    private:

        UnitigMap<U, G, is_const> um;
};

#include "NeighborIterator.tcc"

#endif
