#ifndef BIFROST_MINIMIZER_IDX_HPP
#define BIFROST_MINIMIZER_IDX_HPP

#include <algorithm>
#include <atomic>
#include <iterator>
#include <mutex>
#include <thread>
#include <string>
#include <utility>

#include "fastmod.h"

#include "Kmer.hpp"
#include "Lock.hpp"
#include "TinyVector.hpp"
#include "BooPHF.h"

#define BIFROST_MI_MAX_OCCUPANCY 0.95
#define BIFROST_MI_INIT_SZ 1024

typedef boomphf::mphf<Minimizer, MinimizerHash> boophf_t;

class MinimizerIndex {

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, packed_tiny_vector> {

        public:

            typedef typename std::conditional<is_const, const MinimizerIndex*, MinimizerIndex*>::type MI_ptr_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector&, packed_tiny_vector&>::type MI_tinyv_ref_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector*, packed_tiny_vector*>::type MI_tinyv_ptr_t;
            typedef typename std::conditional<is_const, const uint8_t&, uint8_t&>::type MI_tinyv_sz_ref_t;
            typedef typename std::conditional<is_const, const uint8_t*, uint8_t*>::type MI_tinyv_sz_ptr_t;

            iterator_() : ht(nullptr), h(0xffffffffffffffffULL), psl(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_) : ht(ht_), h(0xffffffffffffffffULL), psl(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_, size_t h_) :  ht(ht_), h(h_), psl(0xffffffffffffffffULL) {}
            iterator_(MI_ptr_t ht_, size_t h_, size_t psl_) :  ht(ht_), h(h_), psl(psl_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h), psl(o.psl) {}

            iterator_& operator=(const iterator_& o) {
  
              if (this != &o) {
                
                ht=o.ht;
                h=o.h;
                psl=o.psl;
              }
              
              return *this;
            }
  
            BFG_INLINE Minimizer getKey() const {

                return ht->table_keys[h];
            }

            BFG_INLINE size_t getHash() const {

                return h;
            }
            
            BFG_INLINE size_t getPSL() const {
              
              return psl;
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

                h += static_cast<size_t>(h < ht->size_);

                while ((h < ht->size_) && ht->table_keys[h].isEmpty()) ++h;

                h |= static_cast<size_t>(h < ht->size_) - 1;
                psl = 0xffffffffffffffffULL;

                return *this;
            }

            BFG_INLINE bool operator==(const iterator_ &o) const {

                return (ht == o.ht) && (h == o.h);
            }

            BFG_INLINE bool operator!=(const iterator_ &o) const {

                return (ht != o.ht) || (h != o.h);
            }
            
            void get_to_first() {
              
              h = 0xffffffffffffffffULL;
              psl = 0xffffffffffffffffULL;
              
              if ((ht != nullptr) && (ht->size_ != 0)) {
                
                h = 0;
                
                while ((h < ht->size_) && ht->table_keys[h].isEmpty()) ++h;
                
                h |= static_cast<size_t>(h < ht->size_) - 1;
              }
            }

            friend class iterator_<true>;

        //private:

            MI_ptr_t ht;
            size_t h;
            size_t psl;
    };

    public:

        template<bool is_const> friend class iterator_;

        typedef iterator_<true> const_iterator;
        typedef iterator_<false> iterator;

        MinimizerIndex();
        MinimizerIndex(const size_t sz, const double max_ratio_occupancy = BIFROST_MI_MAX_OCCUPANCY);

        MinimizerIndex(const MinimizerIndex& o);
        MinimizerIndex(MinimizerIndex&& o);

        MinimizerIndex& operator=(const MinimizerIndex& o);
        MinimizerIndex& operator=(MinimizerIndex&& o);

        ~MinimizerIndex();

        BFG_INLINE size_t size() const {

            return pop;
        }
        
        BFG_INLINE size_t capacity() const {
          
          return size_;
        }

        BFG_INLINE bool empty() const {

            return (pop == 0);
        }

        // Generates an MPHF for all the minimizers
        // gamma=1. yields lowest bit/elem ratio. Higher values yield faster
        // construction and query times. gamma=2. is a good trade-off
        void generate_mphf(std::vector<Minimizer>& minimizers, uint32_t threads=1, float gamma=2.0);
        void register_mphf(boophf_t* mphf_);

        void to_static(uint32_t threads=1, float gamma=1.0);
        void drop_table_keys();

        void clear();
        void clearPTV();

        iterator find(const Minimizer& key);
        const_iterator find(const Minimizer& key) const;

        iterator find(const size_t h);
        const_iterator find(const size_t h) const;

        size_t erase(const_iterator it);
        
        BFG_INLINE size_t erase(const Minimizer& minz) {
          
          const const_iterator it = find(minz);
          
          return erase(it);
        }
        
        pair<iterator, bool> insert(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag);

        pair<iterator, bool> add_unitig_p(const Minimizer& key, const size_t pos_id_unitig); // only if static

        iterator begin();
        const_iterator begin() const;

        iterator end();
        const_iterator end() const;
        
        void recomputeMaxPSL(const size_t nb_threads = 1);

        BFG_INLINE size_t get_mean_psl() const {

          return sum_psl / (pop + 1); // Slightly biased computation but avoids to check for (psl == 0). Fine since we just need an approximate result.
        }
        
        BFG_INLINE size_t get_max_psl() const {

          return max_psl;
        }

        bool is_static;

    private:

        void clear_tables();
        void init_tables(const size_t sz);
        void reserve(const size_t sz);
        void swap(const size_t i, const size_t j);

        double max_ratio_occupancy;

        __uint128_t M_u64;

        size_t size_, pop, max_psl, sum_psl;

        Minimizer* table_keys;
        packed_tiny_vector* table_tinyv;
        uint8_t* table_tinyv_sz;

        boophf_t* mphf;
};

/*class CompactedMinimizerIndex {

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, packed_tiny_vector> {

        public:

            typedef typename std::conditional<is_const, const CompactedMinimizerIndex*, CompactedMinimizerIndex*>::type CMI_ptr_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector&, packed_tiny_vector&>::type CMI_tinyv_ref_t;
            typedef typename std::conditional<is_const, const packed_tiny_vector*, packed_tiny_vector*>::type CMI_tinyv_ptr_t;
            typedef typename std::conditional<is_const, const uint8_t&, uint8_t&>::type CMI_tinyv_sz_ref_t;
            typedef typename std::conditional<is_const, const uint8_t*, uint8_t*>::type CMI_tinyv_sz_ptr_t;

            iterator_() : cmi(nullptr), pos(0xffffffffffffffffULL) {}
            iterator_(CMI_ptr_t cmi_) : cmi(cmi_), pos(cmi_->size_), mi_it(cmi_->mi_overflow.end()) {}
            iterator_(CMI_ptr_t cmi_, size_t pos_) :  cmi(cmi_), pos(pos_), mi_it(cmi_->mi_overflow.begin()) {}
            iterator_(CMI_ptr_t cmi_, size_t pos_, MinimizerIndex::iterator_<is_const> mi_it_) :  cmi(cmi_), pos(pos_), mi_it(mi_it_) {}
            iterator_(const iterator_<false>& o) : cmi(o.cmi), pos(o.pos), mi_it(o.mi_it) {}

            iterator_& operator=(const iterator_& o) {

                if (this != &o) {

                    cmi = o.cmi;
                    pos = o.pos;
                    mi_it = o.mi_it;
                }

                return *this;
            }

            BFG_INLINE Minimizer getKey() const {

                if (pos < cmi->size) return cmi->table_keys[pos];
                
                return mi_it->getKey();
            }

            BFG_INLINE CMI_tinyv_sz_ref_t getVectorSize() const {

                if (pos < cmi->size) return cmi->table_tinyv_sz[pos];

                return mi_it->getVectorSize();
            }

            BFG_INLINE CMI_tinyv_ref_t getVector() const {

                if (pos < cmi->size) return cmi->table_tinyv[pos];

                return mi_it->getVector();
            }

            CMI_tinyv_ref_t operator*() const {

                if (pos < cmi->size) return cmi->table_tinyv[pos];

                return mi_it.operator*();
            }

            CMI_tinyv_ptr_t operator->() const {

                if (pos < cmi->size) return &(cmi->table_tinyv[pos]);

                return mi_it.operator->();
            }

            iterator_ operator++(int) {

                const iterator_ tmp(*this);
                operator++();
                return tmp;
            }

            iterator_& operator++() {

                if (pos == cmi->size) {

                    if (mi_it != cmi->mi_overflow.end()) mi_it.operator++();
                }
                else {

                    for (++pos; pos < cmi->size; ++pos) {

                        if (!cmi->table_keys[pos].isEmpty() && !cmi->table_keys[pos].isDeleted()) break;
                    }

                    // At this point, if (pos == cmi->size), mi_it should already be initialized to mi_overflow.begin() by the constructor
                }

                return *this;
            }

            BFG_INLINE bool operator==(const iterator_ &o) const {

                return (cmi == o.cmi) && (pos == o.pos) && (mi_it == o.mi_it);
            }

            BFG_INLINE bool operator!=(const iterator_ &o) const {

                return (cmi != o.cmi) || (pos != o.pos) || (mi_it != o.mi_it);
            }

            friend class iterator_<true>;

        //private:

            MinimizerIndex::iterator_<is_const> mi_it;
            CMI_ptr_t cmi;
            size_t pos;
    };

    public:

        CompactedMinimizerIndex();
        CompactedMinimizerIndex(const MinimizerIndex& mi, const size_t nb_hashes, const size_t nb_threads = 1);
        CompactedMinimizerIndex(MinimizerIndex&& mi, const size_t nb_hashes, const size_t nb_threads = 1);

        CompactedMinimizerIndex& operator=(const CompactedMinimizerIndex& cmi);
        CompactedMinimizerIndex& operator=(CompactedMinimizerIndex&& cmi);

        void clear();

        iterator begin();
        const_iterator begin() const;

        iterator end();
        const_iterator end() const;

        iterator find(const Minimizer& key);
        const_iterator find(const Minimizer& key) const;

        iterator find(const size_t h);
        const_iterator find(const size_t h) const;

    private:

        size_t size, nb_h, hbits_per_minz;

        Minimizer* table_keys;
        packed_tiny_vector* table_tinyv;
        uint8_t* table_tinyv_sz;
        uint64_t* table_hash_bits;

        MinimizerIndex mi_overflow;
};*/

#endif
