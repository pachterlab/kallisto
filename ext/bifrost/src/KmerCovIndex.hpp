#ifndef BIFROST_KMER_COV_IDX_HPP
#define BIFROST_KMER_COV_IDX_HPP

#include <algorithm>

#include "Common.hpp"
#include "BitContainer.hpp"
#include "Kmer.hpp"
#include "rw_spin_lock.h"

template<typename T = void>
class KmerCovIndex {

    template<typename X> friend class KmerCovIndex;

    public:

        KmerCovIndex();
        ~KmerCovIndex();

        KmerCovIndex(const KmerCovIndex& o); // Copy constructor
        KmerCovIndex(KmerCovIndex&& o); // Move constructor

        KmerCovIndex<T>& operator=(const KmerCovIndex<T>& o); // Copy assignment
        KmerCovIndex<T>& operator=(KmerCovIndex<T>&& o); // Move assignment

        void clear();

        BFG_INLINE size_t size() const {

            return sz;
        }

        void push_back(const Kmer& km);
        bool set(const size_t idx, const Kmer& km);
        bool set(const size_t idx, const Kmer& km, const size_t cov);
        bool swap(const size_t idx1, const size_t idx2);
        void remove(const size_t idx);
        void resize(const size_t new_sz);

        BFG_INLINE static void setFullCoverage(const size_t cov_max) {

            cov_full = min(cov_max, static_cast<size_t>(2));
            cov_full = max(cov_full, static_cast<size_t>(1));
        }

        BFG_INLINE static size_t getFullCoverage() {

            return cov_full;
        }

        BFG_INLINE bool isFull(const size_t idx) const {

            if (idx < sz) return v_blocks[idx >> shift_div]->bc_cov.contains((idx & mask_mod) * cov_full + cov_full - 1);

            return false;
        }

        void setFull(const size_t idx);
        int covAt(const size_t idx) const;

        void cover(const size_t idx);
        void uncover(const size_t idx);

        Kmer getKmer(const size_t idx) const;

        const T* getData(const size_t idx) const;
        T* getData(const size_t idx);

        BFG_INLINE void lock(const size_t idx) {

            if (idx < sz) v_blocks[idx >> shift_div]->lck.acquire();
        }

        BFG_INLINE void unlock(const size_t idx) {

            if (idx < sz) v_blocks[idx >> shift_div]->lck.release();
        }

        BFG_INLINE void cover_thread_safe(const size_t idx) {

            if (idx < sz){

                v_blocks[idx >> shift_div]->lck.acquire();

                cover(idx);

                v_blocks[idx >> shift_div]->lck.release();
            }
        }

        BFG_INLINE void uncover_thread_safe(const size_t idx) {

            if (idx < sz){

                v_blocks[idx >> shift_div]->lck.acquire();

                uncover(idx);

                v_blocks[idx >> shift_div]->lck.release();
            }
        }

        KmerCovIndex<T>& toData(KmerCovIndex<void>&& o, const size_t nb_threads = 1);

    private:

        static const size_t block_sz = 1024; // Always a power of 2

        template<typename U>
        struct Block {

            Kmer km_block[block_sz];
            U data_block[block_sz];

            SpinLock lck;
            BitContainer bc_cov;
        };

        static size_t cov_full;

        size_t shift_div;
        size_t mask_mod;
        size_t sz;

        vector<Block<T>*> v_blocks;
};

// Declare template specializations for type void
// I think it's ugly but without it, compiler can't shut up about explicit specialization after instanciation
// I get the point but still... it's ugly

template<>
template<>
struct KmerCovIndex<void>::Block<void> {

    Kmer km_block[block_sz];

    SpinLock lck;
    BitContainer bc_cov;
};

template<typename T> size_t KmerCovIndex<T>::cov_full = 2;

template<> inline KmerCovIndex<void>::KmerCovIndex(const KmerCovIndex& o);
template<> inline KmerCovIndex<void>& KmerCovIndex<void>::operator=(const KmerCovIndex<void>& o);
template<> inline bool KmerCovIndex<void>::swap(const size_t idx1, const size_t idx2);
template<> inline void KmerCovIndex<void>::remove(const size_t idx1);
template<> inline void KmerCovIndex<void>::resize(const size_t new_sz);
template<> inline const void* KmerCovIndex<void>::getData(const size_t idx) const;
template<> inline void* KmerCovIndex<void>::getData(const size_t idx);

#include "KmerCovIndex.tcc"

#endif
