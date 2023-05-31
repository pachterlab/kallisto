#ifndef BIFROST_UNITIG_ITERATOR_TCC
#define BIFROST_UNITIG_ITERATOR_TCC

#include "CompactedDBG.hpp"

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator() :  i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(nullptr) {}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator(CompactedDBG_ptr_t cdbg_) :
                i(0), v_unitigs_sz(0), v_kmers_sz(0), h_kmers_ccov_sz(0), sz(0), invalid(true), cdbg(cdbg_),
                it_h_kmers_ccov((cdbg_ == nullptr) || cdbg_->invalid ? typename KmerHashTable<CompressedCoverage_t<U>>::const_iterator() : cdbg_->h_kmers_ccov.begin()){

    if ((cdbg != nullptr) && !cdbg->invalid && (cdbg->size() != 0)){

        invalid = false;

        v_unitigs_sz = cdbg->v_unitigs.size();
        v_kmers_sz = cdbg->km_unitigs.size();
        h_kmers_ccov_sz = cdbg->h_kmers_ccov.size();

        sz = v_unitigs_sz + v_kmers_sz + h_kmers_ccov_sz;
    }
}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>::unitigIterator(const unitigIterator& o) :   i(o.i), v_unitigs_sz(o.v_unitigs_sz), v_kmers_sz(o.v_kmers_sz),
                                                                            it_h_kmers_ccov(o.it_h_kmers_ccov), h_kmers_ccov_sz(o.h_kmers_ccov_sz),
                                                                            sz(o.sz), invalid(o.invalid), um(o.um), cdbg(o.cdbg) {}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const>& unitigIterator<U, G, is_const>::operator++() {

    if (invalid) return *this;

    if ((cdbg == nullptr) || cdbg->invalid || (i >= sz)){

        invalid = true;
        return *this;
    }

    if (i < v_unitigs_sz){

        um = UnitigMap<U, G, is_const>(i, 0, cdbg->v_unitigs[i]->getSeq().size() - cdbg->getK() + 1,
                                       cdbg->v_unitigs[i]->getSeq().size(), false, false, true, cdbg);
    }
    else if (i < (v_unitigs_sz + v_kmers_sz)){

        um = UnitigMap<U, G, is_const>(i - v_unitigs_sz, 0, 1, cdbg->getK(), true, false, true, cdbg);
    }
    else {

        um = UnitigMap<U, G, is_const>(it_h_kmers_ccov.getHash(), 0, 1, cdbg->getK(), false, true, true, cdbg);

        ++it_h_kmers_ccov;
    }

    ++i;

    return *this;
}

template<typename U, typename G, bool is_const>
unitigIterator<U, G, is_const> unitigIterator<U, G, is_const>::operator++(int) {

    unitigIterator<U, G, is_const> tmp(*this);
    operator++();

    return tmp;
}

template<typename U, typename G, bool is_const>
bool unitigIterator<U, G, is_const>::operator==(const unitigIterator& o) const {

    if (invalid || o.invalid) return invalid && o.invalid;
    return  (i == o.i) && (v_unitigs_sz == o.v_unitigs_sz) && (v_kmers_sz == o.v_kmers_sz) &&
            (h_kmers_ccov_sz == o.h_kmers_ccov_sz) && (sz == o.sz) && (it_h_kmers_ccov == o.it_h_kmers_ccov) &&
            (cdbg == o.cdbg) && (um == o.um);
}

template<typename U, typename G, bool is_const>
bool unitigIterator<U, G, is_const>::operator!=(const unitigIterator& o) const { return !operator==(o); }

template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>& unitigIterator<U, G, is_const>::operator*() const { return um; }

template<typename U, typename G, bool is_const>
const UnitigMap<U, G, is_const>* unitigIterator<U, G, is_const>::operator->() const { return &um; }

#endif
