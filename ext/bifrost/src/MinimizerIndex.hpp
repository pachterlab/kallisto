#ifndef BIFROST_MINIMIZER_IDX_HPP
#define BIFROST_MINIMIZER_IDX_HPP

#include <utility>
#include <string>
#include <iterator>
#include <algorithm>

#include "Kmer.hpp"
#include "Lock.hpp"
#include "TinyVector.hpp"

class MinimizerIndex {

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, packed_tiny_vector> {

        public:

            typedef typename std::conditional<is_const, const MinimizerIndex*, MinimizerIndex*>::type MI_ptr_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector&, packed_tiny_vector&>::type MI_tinyv_ref_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector*, packed_tiny_vector*>::type MI_tinyv_ptr_t;
            typedef typename std::conditional<is_const, const uint8_t&, uint8_t&>::type MI_tinyv_sz_ref_t;
            typedef typename std::conditional<is_const, const uint8_t*, uint8_t*>::type MI_tinyv_sz_ptr_t;

            iterator_() : ht(nullptr), h(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_) : ht(ht_), h(ht_->size_) {}
            iterator_(MI_ptr_t ht_, size_t h_) :  ht(ht_), h(h_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
            iterator_& operator=(const iterator_& o) { ht=o.ht; h=o.h; return *this; }

            BFG_INLINE Minimizer getKey() const {

                return ht->table_keys[h];
            }

            BFG_INLINE size_t getHash() const {

                return h;
            }

            BFG_INLINE MI_tinyv_sz_ref_t getVectorSize() const {

                return ht->table_tinyv_sz[h];
            }

            BFG_INLINE MI_tinyv_ref_t getVector() const {

                return ht->table_tinyv[h];
            }

            MI_tinyv_ref_t operator*() const {

                return ht->table_tinyv[h];
            }

            MI_tinyv_ptr_t operator->() const {

                return &(ht->table_tinyv[h]);
            }

            iterator_ operator++(int) {

                const iterator_ tmp(*this);
                operator++();
                return tmp;
            }

            iterator_& operator++() {

                if (h == ht->size_) return *this;

                ++h;

                for (; h < ht->size_; ++h) {

                    if (!ht->table_keys[h].isEmpty() && !ht->table_keys[h].isDeleted()) break;
                }

                return *this;
            }

            BFG_INLINE bool operator==(const iterator_ &o) const {

                return (ht == o.ht) && (h == o.h);
            }

            BFG_INLINE bool operator!=(const iterator_ &o) const {

                return (ht != o.ht) || (h != o.h);
            }

            friend class iterator_<true>;

        //private:

            MI_ptr_t ht;
            size_t h;
    };

    public:

        template<bool is_const> friend class iterator_;

        typedef iterator_<true> const_iterator;
        typedef iterator_<false> iterator;

        MinimizerIndex();
        MinimizerIndex(const size_t sz);

        MinimizerIndex(const MinimizerIndex& o);
        MinimizerIndex(MinimizerIndex&& o);

        MinimizerIndex& operator=(const MinimizerIndex& o);
        MinimizerIndex& operator=(MinimizerIndex&& o);

        ~MinimizerIndex();

        BFG_INLINE size_t size() const {

            return pop;
        }

        BFG_INLINE bool empty() const {

            return (pop == 0);
        }

        void clear();

        iterator find(const Minimizer& key);
        const_iterator find(const Minimizer& key) const;

        iterator find(const size_t h);
        const_iterator find(const size_t h) const;

        iterator erase(const_iterator it);
        size_t erase(const Minimizer& minz);

        pair<iterator, bool> insert(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag);

        void init_threads();
        void release_threads();

        iterator find_p(const Minimizer& key);
        const_iterator find_p(const Minimizer& key) const;

        iterator find_p(const size_t h);
        const_iterator find_p(const size_t h) const;

        void release_p(const_iterator it) const;
        void release_p(iterator it);

        size_t erase_p(const Minimizer& minz);

        pair<iterator, bool> insert_p(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag);

        iterator begin();
        const_iterator begin() const;

        iterator end();
        const_iterator end() const;

    private:

        void clear_tables();
        void init_tables(const size_t sz);
        void reserve(const size_t sz);

        size_t size_, pop, num_empty;

        Minimizer* table_keys;
        packed_tiny_vector* table_tinyv;
        uint8_t* table_tinyv_sz;

        mutable vector<SpinLock> lck_min;
        mutable SpinLockRW lck_edit_table;

        atomic<size_t> pop_p, num_empty_p;

        // For future myself: lck_block_sz must be a poswer of 2. If you change it, change lck_block_div_shift accordingly.
        // I could automate this in a much better looking implementation but I'm way too tired today.
        static const size_t lck_block_sz;
        static const size_t lck_block_div_shift;
};

#endif
