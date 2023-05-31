#ifndef BIFROST_UNITIGMAP_TCC
#define BIFROST_UNITIGMAP_TCC

#include "CompactedDBG.hpp"
#include "NeighborIterator.hpp"

template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const>::UnitigMap(const size_t length, CompactedDBG_ptr_t cdbg_) :   UnitigMapBase(length), pos_unitig(0), isShort(false),
                                                                                        isAbundant(false), cdbg(cdbg_) {}

template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const>::UnitigMap(const size_t start, const size_t length, const size_t unitig_sz, const bool strand) :
                                    UnitigMapBase(start, length, unitig_sz, strand), pos_unitig(0),
                                    isShort(false), isAbundant(false), cdbg(nullptr) {}

template<typename U, typename G, bool is_const>
bool UnitigMap<U, G, is_const>::operator==(const UnitigMap& o) const {

    return  UnitigMapBase::operator==(o) && (pos_unitig == o.pos_unitig) &&
            (isShort == o.isShort) && (isAbundant == o.isAbundant) && (cdbg == o.cdbg);
}

template<typename U, typename G, bool is_const>
bool UnitigMap<U, G, is_const>::operator!=(const UnitigMap& o) const {

    return  UnitigMapBase::operator!=(o) || (pos_unitig != o.pos_unitig) ||
            (isShort != o.isShort) || (isAbundant != o.isAbundant) || (cdbg != o.cdbg);
}

template<typename U, typename G, bool is_const>
bool UnitigMap<U, G, is_const>::isSameReferenceUnitig(const UnitigMap& o) const {

    return  (pos_unitig == o.pos_unitig) && (isShort == o.isShort) && (isAbundant == o.isAbundant) && (cdbg == o.cdbg);
}

template<typename U, typename G, bool is_const>
string UnitigMap<U, G, is_const>::mappedSequenceToString() const {

    if (isEmpty) return string();

    if (strand){

        if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig).toString();
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().toString();

        return cdbg->v_unitigs[pos_unitig]->getSeq().toString(dist, len + cdbg->k_ - 1);
    }
    else {

        if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig).twin().toString();
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin().toString();

        return reverse_complement(cdbg->v_unitigs[pos_unitig]->getSeq().toString(dist, len + cdbg->k_ - 1));
    }
}

template<typename U, typename G, bool is_const>
string UnitigMap<U, G, is_const>::referenceUnitigToString() const {

    if (isEmpty) return string();
    if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig).toString();
    if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().toString();

    return cdbg->v_unitigs[pos_unitig]->getSeq().toString();
}

template<typename U, typename G, bool is_const>
size_t UnitigMap<U, G, is_const>::lcp(const char* s, const size_t pos_s, const size_t pos_um_seq, const bool um_reversed) const {

    if (isEmpty || (pos_s >= strlen(s))) return 0;

    if (isShort || isAbundant){

        if (pos_um_seq >= Kmer::k) return 0;

        char km_str[MAX_KMER_SIZE];

        const Kmer km = isShort ? cdbg->km_unitigs.getKmer(pos_unitig) : cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        um_reversed ? km.twin().toString(km_str) : km.toString(km_str);

        return cstrMatch(&s[pos_s], &km_str[pos_um_seq]);
    }

    if (pos_um_seq >= cdbg->v_unitigs[pos_unitig]->length()) return 0;

    return cdbg->v_unitigs[pos_unitig]->getSeq().jump(s, pos_s, pos_um_seq, um_reversed);
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getUnitigHead() const {

    if (!isEmpty){

        if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig);
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(0);
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getUnitigTail() const {

    if (!isEmpty){

        if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig);
        if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(cdbg->v_unitigs[pos_unitig]->numKmers() - 1);
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getUnitigKmer(const size_t pos) const {

    if (!isEmpty){

        if (isShort && (pos == 0)) return cdbg->km_unitigs.getKmer(pos_unitig);
        if (isAbundant && (pos == 0)) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

        if (!isShort && !isAbundant && (pos < cdbg->v_unitigs[pos_unitig]->numKmers())) {

            return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(pos);
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getMappedHead() const {

    if (!isEmpty){

        if (strand){

            if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig);
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

            return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(dist);
        }
        else {

            if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig).twin();
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin();

            return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(dist + len - 1).twin();
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getMappedTail() const {

    if (!isEmpty){

        if (strand){

            if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig);
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

            return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(dist + len - 1);
        }
        else {

            if (isShort) return cdbg->km_unitigs.getKmer(pos_unitig).twin();
            if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin();

            return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(dist).twin();
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
Kmer UnitigMap<U, G, is_const>::getMappedKmer(const size_t pos) const {

    if (!isEmpty && (pos < len)){

        if (strand){

            if (isShort && (pos + dist == 0)) return cdbg->km_unitigs.getKmer(pos_unitig);
            if (isAbundant && (pos + dist == 0)) return cdbg->h_kmers_ccov.find(pos_unitig).getKey();

            if (!isShort && !isAbundant && (pos < cdbg->v_unitigs[pos_unitig]->numKmers())) {

                return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(pos);
            }
        }
        else {

            if (isShort && (pos + dist == 0)) return cdbg->km_unitigs.getKmer(pos_unitig).twin();
            if (isAbundant && (pos + dist == 0)) return cdbg->h_kmers_ccov.find(pos_unitig).getKey().twin();

            if (!isShort && !isAbundant && (pos < cdbg->v_unitigs[pos_unitig]->numKmers())) {

                return cdbg->v_unitigs[pos_unitig]->getSeq().getKmer(pos).twin();
            }
        }
    }

    Kmer km;

    km.set_empty();

    return km;
}

template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const> UnitigMap<U, G, is_const>::getKmerMapping(const size_t pos) const {

    UnitigMap<U, G, is_const> um(*this);

    um.dist = pos;
    um.len = 1;

    if (pos > um.size - cdbg->getK()) um.isEmpty = true;

    return um;
}

template<typename U, typename G, bool is_const>
typename UnitigMap<U, G, is_const>::Unitig_data_ptr_t UnitigMap<U, G, is_const>::getData() const {

    return getData_<is_void<U>::value>();
}

template<typename U, typename G, bool is_const>
typename UnitigMap<U, G, is_const>::UnitigMap_BW UnitigMap<U, G, is_const>::getPredecessors() const {

    return BackwardCDBG<U, G, is_const>(*this);
}

template<typename U, typename G, bool is_const>
typename UnitigMap<U, G, is_const>::UnitigMap_FW UnitigMap<U, G, is_const>::getSuccessors() const {

    return ForwardCDBG<U, G, is_const>(*this);
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<!is_void, typename UnitigMap<U, G, is_const>::Unitig_data_ptr_t>::type UnitigMap<U, G, is_const>::getData_() const {

    if (isEmpty) return nullptr;
    if (isShort) return cdbg->km_unitigs.getData(pos_unitig);
    if (isAbundant) return cdbg->h_kmers_ccov.find(pos_unitig)->getData();

    return cdbg->v_unitigs[pos_unitig]->getData();
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<is_void, typename UnitigMap<U, G, is_const>::Unitig_data_ptr_t>::type UnitigMap<U, G, is_const>::getData_() const {

    return nullptr;
}


template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<!is_void, Unitig<U>>::type UnitigMap<U, G, is_const>::splitData_(const bool last_split) const {

    Unitig<U> unitig;

    unitig.getData()->extract(*this, last_split);

    return unitig;
}

template<typename U, typename G, bool is_const>
template<bool is_void>
typename std::enable_if<is_void, Unitig<U>>::type UnitigMap<U, G, is_const>::splitData_(const bool last_split) const {

    return Unitig<U>();
}

template<typename U, typename G, bool is_const>
Unitig<U> UnitigMap<U, G, is_const>::splitData(const bool last_split) const {

    return splitData_<is_void<U>::value>(last_split);
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::bw_begin() const {

    neighborIterator<U, G, is_const> it(*this, false);
    it++;
    return it;
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::bw_end() const { return neighborIterator<U, G, is_const>(); }

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::fw_begin() const {

    neighborIterator<U, G, is_const> it(*this, true);
    it++;
    return it;
}

template<typename U, typename G, bool is_const>
neighborIterator<U, G, is_const> UnitigMap<U, G, is_const>::fw_end() const { return neighborIterator<U, G, is_const>(); }

template<typename U, typename G, bool is_const>
void UnitigMap<U, G, is_const>::partialCopy(const UnitigMap<U, G, is_const>& o) {

    pos_unitig = o.pos_unitig;
    dist = o.dist;
    len = o.len;
    size = o.size;

    strand = o.strand;
    isEmpty = o.isEmpty;

    isShort = o.isShort;
    isAbundant = o.isAbundant;
}

template<typename U, typename G, bool is_const>
void UnitigMap<U, G, is_const>::setFullCoverage() const {

    if (!isEmpty){

        if (isShort) cdbg->km_unitigs.setFull(pos_unitig);
        else if (isAbundant) cdbg->h_kmers_ccov.find(pos_unitig)->ccov.setFull();
        else cdbg->v_unitigs[pos_unitig]->getCov().setFull();
    }
}

template<typename U, typename G, bool is_const>
void UnitigMap<U, G, is_const>::increaseCoverage() const {

    if (isEmpty) return; // nothing maps, move on

    if (isShort) cdbg->km_unitigs.cover(pos_unitig, dist, dist + len - 1);
    else if (isAbundant) cdbg->h_kmers_ccov.find(pos_unitig)->ccov.cover(dist, dist + len - 1);
    else cdbg->v_unitigs[pos_unitig]->getCov().cover(dist, dist + len - 1);
}

template<typename U, typename G, bool is_const>
void UnitigMap<U, G, is_const>::decreaseCoverage() const {

    if (isEmpty) return; // nothing maps, move on

    if (isShort) cdbg->km_unitigs.uncover(pos_unitig, dist, dist + len - 1);
    else if (isAbundant) cdbg->h_kmers_ccov.find(pos_unitig)->ccov.uncover(dist, dist + len - 1);
    else cdbg->v_unitigs[pos_unitig]->getCov().uncover(dist, dist + len - 1);
}

template<typename U, typename G, bool is_const>
bool UnitigMap<U, G, is_const>::isCoverageFull() const {

    if (isEmpty) return false; // nothing maps, move on

    if (isShort) cdbg->km_unitigs.isFull(pos_unitig);
    else if (isAbundant) cdbg->h_kmers_ccov.find(pos_unitig)->ccov.isFull();
    else cdbg->v_unitigs[pos_unitig]->getCov().isFull();
}

template<typename U, typename G, bool is_const>
size_t UnitigMap<U, G, is_const>::getCoverage(const size_t pos) const {

    if (isEmpty || (pos > size - cdbg->getK())) return 0; // nothing maps, move on

    if (isShort) cdbg->km_unitigs.covAt(pos_unitig);
    else if (isAbundant) cdbg->h_kmers_ccov.find(pos_unitig)->ccov.covAt(0);
    else cdbg->v_unitigs[pos_unitig]->getCov().covAt(pos);
}


template<typename U, typename G, bool is_const>
UnitigMap<U, G, is_const>::UnitigMap(size_t p_unitig, size_t i, size_t l, size_t sz, bool short_, bool abundance, bool strd,
                                     CompactedDBG_ptr_t cdbg_) : UnitigMapBase(i, l, sz, strd), pos_unitig(p_unitig),
                                     isShort(short_), isAbundant(abundance), cdbg(cdbg_) {}

#endif
