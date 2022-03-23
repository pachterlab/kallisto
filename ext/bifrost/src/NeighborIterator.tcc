#ifndef BIFROST_NEIGHBOR_ITERATOR_TCC
#define BIFROST_NEIGHBOR_ITERATOR_TCC

#include "CompactedDBG.hpp"
#include "UnitigMap.hpp"

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const>::neighborIterator() : i(4), is_fw(true), cdbg(nullptr) {}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const>::neighborIterator(const UnitigMap<U, G, is_const>& um_, const bool is_forward_) : i(-1), is_fw(is_forward_), cdbg(um_.getGraph()) {

    if (um_.isEmpty || (cdbg == nullptr) || cdbg->invalid) i = 4;
    else {

        km_head = um_.strand ? um_.getUnitigHead() : um_.getUnitigTail().twin();
        km_tail = um_.strand ? um_.getUnitigTail() : um_.getUnitigHead().twin();
    }
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const>::neighborIterator(const neighborIterator& o) : i(o.i), is_fw(o.is_fw), um(o.um), km_head(o.km_head), km_tail(o.km_tail), cdbg(o.cdbg) {}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const>& neighborIterator<U, G, is_const>::operator++() {

    if ((cdbg == NULL) || cdbg->invalid || (i >= 4)) return *this;

    ++i;

    while (i < 4){

        um = cdbg->find(is_fw ? km_tail.forwardBase(alpha[i]) : km_head.backwardBase(alpha[i]), true);

        if (!um.isEmpty){

            um.dist = 0;
            um.len = um.size - cdbg->getK() + 1;

            break;
        }

        ++i;
    }

    return *this;
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> neighborIterator<U, G, is_const>::operator++(int) {

    neighborIterator tmp(*this);
    operator++();

    return tmp;
}

template<typename U, typename G, bool is_const>
bool neighborIterator<U, G, is_const>::operator==(const neighborIterator& o) const {

    if ((i >= 4) || (o.i >= 4)) return (i >= 4) && (o.i >= 4);
    return (is_fw == o.is_fw) && (km_head == o.km_head) && (km_tail == o.km_tail) && (cdbg == o.cdbg) && (um == o.um);
}

template<typename U, typename G, bool is_const>
bool neighborIterator<U, G, is_const>::operator!=(const neighborIterator& o) const { return !operator==(o); }

template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>& neighborIterator<U, G, is_const>::operator*() const { return um; }

template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>* neighborIterator<U, G, is_const>::operator->() const { return &um; }




template<typename U, typename G, bool is_const>
BackwardCDBG<U, G, is_const>::BackwardCDBG(const UnitigMap<U, G, is_const>& um_) : um(um_) {}

template<typename U, typename G, bool is_const>
bool BackwardCDBG<U, G, is_const>::hasPredecessors() const { return (um.bw_begin() != um.bw_end()); }

template<typename U, typename G, bool is_const>
size_t BackwardCDBG<U, G, is_const>::cardinality() const {

    neighborIterator<U, G, is_const> it_start = begin();
    neighborIterator<U, G, is_const> it_end = end();

    size_t card = 0;

    while (it_start != it_end){

        ++it_start;
        ++card;
    }

    return card;
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> BackwardCDBG<U, G, is_const>::begin() const { return um.bw_begin(); }

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> BackwardCDBG<U, G, is_const>::end() const { return um.bw_end(); }




template<typename U, typename G, bool is_const>
ForwardCDBG<U, G, is_const>::ForwardCDBG(const UnitigMap<U, G, is_const>& um_) : um(um_) {}

template<typename U, typename G, bool is_const>
bool ForwardCDBG<U, G, is_const>::hasSuccessors() const { return (um.fw_begin() != um.fw_end()); }

template<typename U, typename G, bool is_const>
size_t ForwardCDBG<U, G, is_const>::cardinality() const {

    neighborIterator<U, G, is_const> it_start = begin();
    neighborIterator<U, G, is_const> it_end = end();

    size_t card = 0;

    while (it_start != it_end){

        ++it_start;
        ++card;
    }

    return card;
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> ForwardCDBG<U, G, is_const>::begin() const { return um.fw_begin(); }

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> ForwardCDBG<U, G, is_const>::end() const { return um.fw_end(); }

#endif
