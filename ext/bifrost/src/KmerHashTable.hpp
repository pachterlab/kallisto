#ifndef BIFROST_KMER_HASHTABLE_HPP
#define BIFROST_KMER_HASHTABLE_HPP

#include <utility>
#include <string>
#include <iterator>
#include <algorithm>

#include "Kmer.hpp"

#define BIFROST_KHT_MAX_OCCUPANCY 0.95
#define BIFROST_KHT_INIT_SZ 1024

template<typename T>
struct KmerHashTable {

    double max_ratio_occupancy;

    __uint128_t M_u64;

    size_t size_, pop, max_psl, sum_psl;

    Kmer* table_keys;
    T* table_values;

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, T> {

        public:

            typedef typename std::conditional<is_const, const KmerHashTable*, KmerHashTable*>::type KHT_ptr_t;
            typedef typename std::conditional<is_const, const T&, T&>::type KHT_val_ref_t;
            typedef typename std::conditional<is_const, const T*, T*>::type KHT_val_ptr_t;

            iterator_() : ht(nullptr), h(0xffffffffffffffffULL), psl(0xffffffffffffffffULL) {}
            iterator_(KHT_ptr_t ht_) : ht(ht_), h(0xffffffffffffffffULL), psl(0xffffffffffffffffULL) {}
            iterator_(KHT_ptr_t ht_, size_t h_) :  ht(ht_), h(h_), psl(0xffffffffffffffffULL) {}
            iterator_(KHT_ptr_t ht_, size_t h_, size_t psl_) :  ht(ht_), h(h_), psl(psl_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h), psl(o.psl) {}

            iterator_& operator=(const iterator_& o) {

                if (this != &o) {

                    ht=o.ht;
                    h=o.h;
                    psl=o.psl;
                }

                return *this;
            }

            BFG_INLINE Kmer getKey() const {

                return ht->table_keys[h];
            }

            BFG_INLINE size_t getHash() const {

                return h;
            }

            BFG_INLINE size_t getPSL() const {

                return psl;
            }

            BFG_INLINE KHT_val_ref_t operator*() const {

                return ht->table_values[h];
            }

            BFG_INLINE KHT_val_ptr_t operator->() const {

                return &(ht->table_values[h]);
            }

            iterator_ operator++(int) {

                const iterator_ tmp(*this);
                operator++();
                return tmp;
            }

            iterator_& operator++() {

                h += static_cast<size_t>(h < ht->size_);
                psl = 0xffffffffffffffffULL;

                while ((h < ht->size_) && ht->table_keys[h].isEmpty()) ++h;

                h |= static_cast<size_t>(h < ht->size_) - 1;

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

            KHT_ptr_t ht;

            size_t h;
            size_t psl;
    };

    typedef iterator_<true> const_iterator;
    typedef iterator_<false> iterator;

    KmerHashTable() :   table_keys(nullptr), table_values(nullptr) {

        clear();
    }

    KmerHashTable(const size_t sz, const double ratio_occupancy = BIFROST_KHT_MAX_OCCUPANCY) : table_keys(nullptr), table_values(nullptr) {

        clear();

        max_ratio_occupancy = ratio_occupancy;

        if (sz != 0) {

            const size_t sz_with_empty = static_cast<size_t>((1.0 + (1.0 - ratio_occupancy)) * sz);

            init_tables(sz_with_empty);
        }
    }

    KmerHashTable(const KmerHashTable& o) : table_keys(nullptr), table_values(nullptr), max_ratio_occupancy(o.max_ratio_occupancy),
                                            size_(o.size_), pop(o.pop), sum_psl(o.sum_psl), max_psl(o.max_psl), M_u64(o.M_u64) {

        if (size_ != 0) {

            table_keys = new Kmer[size_];
            table_values = new T[size_];

            std::copy(o.table_keys, o.table_keys + size_, table_keys);
            std::copy(o.table_values, o.table_values + size_, table_values);
        }
    }

    KmerHashTable(KmerHashTable&& o) :  table_keys(o.table_keys), table_values(o.table_values),
                                        size_(o.size_), pop(o.pop), sum_psl(o.sum_psl), max_psl(o.max_psl),
                                        M_u64(o.M_u64), max_ratio_occupancy(o.max_ratio_occupancy) {

        o.table_keys = nullptr;
        o.table_values = nullptr;

        o.clear();
    }

    KmerHashTable& operator=(const KmerHashTable& o) {

        if (this != &o) {

            clear();

            size_ = o.size_;
            pop = o.pop;
            sum_psl = o.sum_psl;
            max_psl = o.max_psl;
            max_ratio_occupancy = o.max_ratio_occupancy;
            M_u64 = o.M_u64;

            if (size_ != 0) {

                table_keys = new Kmer[size_];
                table_values = new T[size_];

                std::copy(o.table_keys, o.table_keys + size_, table_keys);
                std::copy(o.table_values, o.table_values + size_, table_values);
            }
        }

        return *this;
    }

    KmerHashTable& operator=(KmerHashTable&& o){

        if (this != &o) {

            clear();

            size_ = o.size_;
            pop = o.pop;
            sum_psl = o.sum_psl;
            max_psl = o.max_psl;
            max_ratio_occupancy = o.max_ratio_occupancy;
            M_u64 = o.M_u64;

            table_keys = o.table_keys;
            table_values = o.table_values;

            o.table_keys = nullptr;
            o.table_values = nullptr;

            o.clear();
        }

        return *this;
    }

    ~KmerHashTable() {

        clear();
    }

    void clear() {

        if (table_keys != nullptr) {

            delete[] table_keys;
            table_keys = nullptr;
        }

        if (table_values != nullptr) {

            delete[] table_values;
            table_values = nullptr;
        }

        size_ = 0;
        pop  = 0;
        sum_psl = 0;
        max_psl = 1;
        M_u64 = 0;
        max_ratio_occupancy = BIFROST_KHT_MAX_OCCUPANCY;
    }

    KmerHashTable::iterator find(const Kmer& key) {

        if ((pop != 0) && (size_ != 0)) {

            const size_t end_table = size_-1;
            const size_t mean_psl = get_mean_psl();

            size_t psl = 0;
            size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);

            if (mean_psl <= 2) {

                while ((psl != max_psl) && !table_keys[h].isEmpty() && (table_keys[h] != key)) {

                    h = (h+1) & (static_cast<size_t>(h == end_table) - 1);
                    ++psl;
                }

                if ((psl != max_psl) && (table_keys[h] == key)) return iterator(this, h, psl);
            }
            else {

                size_t h_mean = fastmod::fastmod_u64(h + mean_psl, M_u64, size_);
                size_t h_inc = h_mean, h_dec = h_mean;

                bool has_empty_key = false;

                size_t i = 0;

                // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
                for (; !has_empty_key && (i <= mean_psl); ++i) {

                    if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);

                    has_empty_key = table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty();

                    if (!has_empty_key && (table_keys[h_inc] == key)) return iterator(this, h_inc, mean_psl + i);

                    h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
                    h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
                }

                // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
                for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

                    if (table_keys[h_inc] == key) return iterator(this, h_inc, mean_psl + i);

                    has_empty_key = table_keys[h_inc].isEmpty();
                    h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
                }

                if (has_empty_key) { // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered

                    for (; (i <= mean_psl); ++i) {

                        if (table_keys[h_dec] == key) return iterator(this, h_dec, mean_psl - i);

                        h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
                    }
                }
            }
        }

        return iterator(this);
    }

    KmerHashTable::const_iterator find(const Kmer& key) const {

         if ((pop != 0) && (size_ != 0)) {

            const size_t end_table = size_-1;
            const size_t mean_psl = get_mean_psl();

            size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);

            if (mean_psl <= 2) {

                size_t psl = 0;

                while ((psl != max_psl) && !table_keys[h].isEmpty() && (table_keys[h] != key)) {

                    h = (h+1) & (static_cast<size_t>(h==end_table)-1);
                    ++psl;
                }

                if ((psl != max_psl) && (table_keys[h] == key)) return const_iterator(this, h, psl);
            }
            else {

                size_t h_mean = fastmod::fastmod_u64(h + mean_psl, M_u64, size_);
                size_t h_inc = h_mean, h_dec = h_mean;

                bool has_empty_key = false;

                size_t i = 0;

                // Check all elements located at positions mean + i and mean -i. Stop if empty key encountered. Stop if minimum key is encountered.
                for (; !has_empty_key && (i <= mean_psl); ++i) {

                    if (table_keys[h_dec] == key) return const_iterator(this, h_dec, mean_psl - i);

                    has_empty_key = table_keys[h_dec].isEmpty() || table_keys[h_inc].isEmpty();

                    if (!has_empty_key && (table_keys[h_inc] == key)) return const_iterator(this, h_inc, mean_psl + i);

                    h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
                    h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
                }

                // Only check remaining acending positions if neither the empty key nor the minimum key were encountered
                for (; !has_empty_key && (mean_psl + i <= max_psl); ++i) {

                    if (table_keys[h_inc] == key) return const_iterator(this, h_inc, mean_psl + i);

                    has_empty_key = table_keys[h_inc].isEmpty();
                    h_inc = (h_inc + 1) & (static_cast<size_t>(h_inc == end_table) - 1);
                }

                // Only check remaining descending positions if empty key was previously encountered but minimum key was not encountered
                if (has_empty_key) {

                    for (; (i <= mean_psl); ++i) {

                        if (table_keys[h_dec] == key) return const_iterator(this, h_dec, mean_psl - i);

                        h_dec = (h_dec - 1) + ((static_cast<size_t>(h_dec != 0) - 1) & size_);
                    }
                }
            }
        }

        return const_iterator(this);
    }

    KmerHashTable::iterator find(const size_t h) {

        if ((h < size_) && !table_keys[h].isEmpty()) return iterator(this, h);

        return iterator(this);
    }

    KmerHashTable::const_iterator find(const size_t h) const {

        if ((h < size_) && !table_keys[h].isEmpty()) return const_iterator(this, h);

        return const_iterator(this);
    }

    size_t erase(const_iterator it) {

        if ((size_ == 0) || (it == end())) return 0;

        const size_t end_table = size_-1;
        
        // Remove <key, value>
        {
            if (it.psl != 0xffffffffffffffffULL) sum_psl -= it.psl;
            else {

                const size_t h = fastmod::fastmod_u64(table_keys[it.h].hash(), M_u64, size_);

                sum_psl -= (size_ - h + it.h) & (static_cast<size_t>(it.h >= h) - 1);
                sum_psl -= (it.h - h) & (static_cast<size_t>(it.h < h) - 1);
            }

            table_keys[it.h].set_empty();
            table_values[it.h] = T();

            --pop;
        }

        // Robin-hood hash: push the tombstone further away if subsequent keys can be closer to where they are supposed to be
        {
            size_t i = 0;
            size_t j1 = it.h;
            size_t j2 = (it.h + 1) & (static_cast<size_t>(it.h == end_table) - 1);

            while ((i != size_) && !table_keys[j2].isEmpty() && (fastmod::fastmod_u64(table_keys[j2].hash(), M_u64, size_) != j2)) {

                KmerHashTable::swap(j1, j2);

                j1 = j2;
                j2 = (j1 + 1) & (static_cast<size_t>(j1 == end_table) - 1);

                --sum_psl; 
                ++i;
            }
        }

        return 1;
    }

    BFG_INLINE size_t erase(const Kmer& key) {

        const const_iterator it = find(key);

        return erase(it);
    }

    // Insert with Robin Hood hashing
    pair<KmerHashTable::iterator, bool> insert(const Kmer& key, const T& val) {

        if (size_ == 0) init_tables(BIFROST_KHT_INIT_SZ);
        else if (pop >= static_cast<size_t>(size_ * max_ratio_occupancy)) {

            size_t resize = 1.2 * size_;

            while (pop >= static_cast<size_t>(resize * max_ratio_occupancy)) resize *= 1.2;

            reserve(resize);
        }

        const size_t end_table = size_-1;

        bool has_rich_psl = false, cascade_ins = false;

        size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);
        size_t h_rich_psl_ins = 0;
        size_t psl_ins_key = 0, psl_rich_key = 0, psl_curr_key = 0;

        pair<KmerHashTable::iterator, bool> it_ret;

        Kmer l_key = key;
        T l_val = val;

        ++pop; // Pre-emptively increase population

        while (true) {

            if (table_keys[h].isEmpty() || (has_rich_psl && (cascade_ins || (psl_ins_key >= max_psl)))) {

                if (has_rich_psl) {

                    h = h_rich_psl_ins;

                    // Swap keys
                    {
                        Kmer key_tmp(std::move(table_keys[h]));

                        table_keys[h] = std::move(l_key); 
                        l_key = std::move(key_tmp);
                    }

                    // Swap values
                    {
                        T val_tmp(std::move(table_values[h]));

                        table_values[h] = std::move(l_val); 
                        l_val = std::move(val_tmp);
                    }

                    if (!cascade_ins) it_ret = {iterator(this, h, psl_rich_key), true};

                    max_psl = max(max_psl, psl_rich_key + 1);
                    sum_psl -= psl_curr_key;
                    sum_psl += psl_rich_key;

                    psl_ins_key = psl_curr_key;
                    has_rich_psl = false;
                    cascade_ins = true;
                }
                else {

                    table_keys[h] = std::move(l_key);
                    table_values[h] = std::move(l_val);

                    max_psl = max(max_psl, psl_ins_key + 1);
                    sum_psl += psl_ins_key;

                    if (!cascade_ins) it_ret = {iterator(this, h, psl_ins_key), true};

                    return it_ret;
                }
            }
            else if (table_keys[h] == l_key) {

                --pop; // Key already in there, pre-emptive population increase was not necessary

                return {iterator(this, h, psl_ins_key), false}; // Can only happen when inserting the input key
            }
            else if (!has_rich_psl) {

                const size_t h_curr = fastmod::fastmod_u64(table_keys[h].hash(), M_u64, size_);
                
                psl_curr_key = ((size_ - h_curr + h) & (static_cast<size_t>(h >= h_curr) - 1)) + ((h - h_curr) & (static_cast<size_t>(h < h_curr) - 1));

                if (psl_ins_key > psl_curr_key) {

                    h_rich_psl_ins = h;
                    psl_rich_key = psl_ins_key;
                    has_rich_psl = true;
                }
            }

            h = (h + 1) & (static_cast<size_t>(h == end_table) - 1);
            ++psl_ins_key;
        }
    }

    // Insert with Robin Hood hashing
    pair<KmerHashTable::iterator, bool> insert(Kmer&& key, T&& val) {

        if (size_ == 0) init_tables(BIFROST_KHT_INIT_SZ);
        else if (pop >= static_cast<size_t>(size_ * max_ratio_occupancy)) {

            size_t resize = 1.2 * size_;

            while (pop >= static_cast<size_t>(resize * max_ratio_occupancy)) resize *= 1.2;

            reserve(resize);
        }

        const size_t end_table = size_-1;

        bool has_rich_psl = false, cascade_ins = false;

        size_t h = fastmod::fastmod_u64(key.hash(), M_u64, size_);
        size_t h_rich_psl_ins = 0;
        size_t psl_ins_key = 0, psl_rich_key = 0, psl_curr_key = 0;

        pair<KmerHashTable::iterator, bool> it_ret;

        Kmer l_key = std::move(key);
        T l_val = std::move(val);

        ++pop; // Pre-emptively increase population

        while (true) {

            if (table_keys[h].isEmpty() || (has_rich_psl && (cascade_ins || (psl_ins_key >= max_psl)))) {

                if (has_rich_psl) {

                    h = h_rich_psl_ins;

                    // Swap keys
                    {
                        Kmer key_tmp(std::move(table_keys[h]));

                        table_keys[h] = std::move(l_key); 
                        l_key = std::move(key_tmp);
                    }

                    // Swap values
                    {
                        T val_tmp(std::move(table_values[h]));

                        table_values[h] = std::move(l_val); 
                        l_val = std::move(val_tmp);
                    }

                    if (!cascade_ins) it_ret = {iterator(this, h, psl_rich_key), true};

                    max_psl = max(max_psl, psl_rich_key + 1);
                    sum_psl -= psl_curr_key;
                    sum_psl += psl_rich_key;

                    psl_ins_key = psl_curr_key;
                    has_rich_psl = false;
                    cascade_ins = true;
                }
                else {

                    table_keys[h] = std::move(l_key);
                    table_values[h] = std::move(l_val);

                    max_psl = max(max_psl, psl_ins_key + 1);
                    sum_psl += psl_ins_key;

                    if (!cascade_ins) it_ret = {iterator(this, h, psl_ins_key), true};

                    return it_ret;
                }
            }
            else if (table_keys[h] == l_key) {

                --pop; // Key already in there, pre-emptive population increase was not necessary

                return {iterator(this, h, psl_ins_key), false}; // Can only happen when inserting the input key
            }
            else if (!has_rich_psl) {

                const size_t h_curr = fastmod::fastmod_u64(table_keys[h].hash(), M_u64, size_);
                
                psl_curr_key = ((size_ - h_curr + h) & (static_cast<size_t>(h >= h_curr) - 1)) + ((h - h_curr) & (static_cast<size_t>(h < h_curr) - 1));

                if (psl_ins_key > psl_curr_key) {

                    h_rich_psl_ins = h;
                    psl_rich_key = psl_ins_key;
                    has_rich_psl = true;
                }
            }

            h = (h + 1) & (static_cast<size_t>(h == end_table) - 1);
            ++psl_ins_key;
        }
    }

    void recomputeMaxPSL(const size_t nb_threads = 1) {

        max_psl = 1;

        if ((pop != 0) && (size_ != 0)) {

            if (nb_threads <= 1){

                for (size_t i = 0; i != size_; ++i) {

                    if (!table_keys[i].isEmpty()) {

                        const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                        const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                        max_psl = max(max_psl, psl + 1);
                    }
                }
            }
            else {

                const size_t chunk_per_thread = (size_ + nb_threads - 1) / nb_threads;

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mtx_max_psl;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&, t]{

                            const size_t chunk_start = t * chunk_per_thread;
                            const size_t chunk_end = min(((t+1) * chunk_per_thread), size_);

                            size_t l_max_psl = 1;

                            for (size_t i = chunk_start; i < chunk_end; ++i) {

                                if (!table_keys[i].isEmpty()) {

                                    const size_t h = fastmod::fastmod_u64(table_keys[i].hash(), M_u64, size_);
                                    const size_t psl = ((size_ - h + i) & (static_cast<size_t>(i >= h) - 1)) + ((i - h) & (static_cast<size_t>(i < h) - 1));

                                    l_max_psl = max(l_max_psl, psl + 1);
                                }
                            }

                            {
                                unique_lock<mutex> lock(mtx_max_psl);

                                max_psl = max(max_psl, l_max_psl);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }
        }
    }

    BFG_INLINE size_t size() const {

        return pop;
    }

    BFG_INLINE size_t capacity() const {

        return size_;
    }

    BFG_INLINE bool empty() const {

        return (pop == 0);
    }

    KmerHashTable::iterator begin() {

        iterator it(this);

        it.get_to_first();

        return it;
    }

    KmerHashTable::const_iterator begin() const {

        const_iterator it(this);

        it.get_to_first();

        return it;
    }

    KmerHashTable::iterator end() {

        return iterator(this);
    }

    KmerHashTable::const_iterator end() const {

        return const_iterator(this);
    }

    BFG_INLINE size_t get_mean_psl() const {

        return sum_psl / (pop + 1); // Slightly biased computation but avoids to check for (psl == 0). Fine since we just need an approximate result.
    }

    BFG_INLINE size_t get_max_psl() const {

        return max_psl;
    }

    void init_tables(const size_t sz) {

        clear();

        Kmer empty_key;

        pop = 0;
        size_ = sz;
        M_u64 = fastmod::computeM_u64(size_);

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        empty_key.set_empty();

        std::fill(table_keys, table_keys + size_, empty_key);
    }

    void reserve(const size_t sz) {

        if (sz <= size_) return;

        if (size_ == 0) init_tables(sz);
        else {

            const size_t old_size_ = size_;

            Kmer empty_key;

            Kmer* old_table_keys = table_keys;
            T* old_table_values = table_values;

            size_ = sz;
            pop = 0;
            sum_psl = 0;
            max_psl = 1;
            M_u64 = fastmod::computeM_u64(size_);

            table_keys = new Kmer[size_];
            table_values = new T[size_];

            empty_key.set_empty();

            std::fill(table_keys, table_keys + size_, empty_key);

            for (size_t i = 0; i < old_size_; ++i) {

                if (!old_table_keys[i].isEmpty()) insert(std::move(old_table_keys[i]), std::move(old_table_values[i]));
            }

            delete[] old_table_keys;
            delete[] old_table_values;
        }
    }

    void swap(const size_t i, const size_t j) {

       // Swap keys
        {
            Kmer key_tmp(std::move(table_keys[i]));

            table_keys[i] = std::move(table_keys[j]); 
            table_keys[j] = std::move(key_tmp);
        }

        // Swap values
        {
            T val_tmp(std::move(table_values[i]));

            table_values[i] = std::move(table_values[j]); 
            table_values[j] = std::move(val_tmp);
        }
    }
};

/*template<typename T>
struct KmerHashTable {

    size_t size_, pop, num_empty;

    Kmer* table_keys;
    T* table_values;

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, T> {

        public:

            typedef typename std::conditional<is_const, const KmerHashTable *, KmerHashTable *>::type MHT_ptr_t;
            typedef typename std::conditional<is_const, const T&, T&>::type MHT_val_ref_t;
            typedef typename std::conditional<is_const, const T*, T*>::type MHT_val_ptr_t;

            MHT_ptr_t ht;
            size_t h;

            iterator_() : ht(nullptr), h(0) {}
            iterator_(MHT_ptr_t ht_) : ht(ht_), h(ht_->size_) {}
            iterator_(MHT_ptr_t ht_, size_t h_) :  ht(ht_), h(h_) {}
            iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
            iterator_& operator=(const iterator_& o) { ht=o.ht; h=o.h; return *this; }

            MHT_val_ref_t operator*() const { return ht->table_values[h]; }
            MHT_val_ptr_t operator->() const { return &(ht->table_values[h]); }

            const Kmer& getKey() const { return ht->table_keys[h]; }

            size_t getHash() const { return h; }

            void find_first() {

                h = 0;

                if ((ht != nullptr) && (ht->size_ > 0) && (ht->table_keys[h].isEmpty() || ht->table_keys[h].isDeleted())) operator++();
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

            bool operator==(const iterator_ &o) const { return (ht == o.ht) && (h == o.h); }
            bool operator!=(const iterator_ &o) const { return (ht != o.ht) || (h != o.h); }

            friend class iterator_<true>;
        };

    typedef iterator_<true> const_iterator;
    typedef iterator_<false> iterator;

    // --- hash table
    KmerHashTable() : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        init_tables(1024);
    }

    KmerHashTable(const size_t sz) : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        if (sz < 2) init_tables(2);
        else {

            const size_t sz_with_empty = static_cast<size_t>(1.2 * sz);

            size_t rdnup_sz = rndup(sz);

            while (rdnup_sz < sz_with_empty) rdnup_sz <<= 1;

            init_tables(max(rdnup_sz, static_cast<size_t>(2)));
        }
    }

    KmerHashTable(const KmerHashTable& o) : size_(o.size_), pop(o.pop), num_empty(o.num_empty) {

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);
    }

    KmerHashTable(KmerHashTable&& o){

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = o.table_keys;
        table_values = o.table_values;

        o.table_keys = nullptr;
        o.table_values = nullptr;

        o.clear_tables();
    }

    KmerHashTable& operator=(const KmerHashTable& o) {

        clear_tables();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);

        return *this;
    }

    KmerHashTable& operator=(KmerHashTable&& o){

        if (this != &o) {

            clear_tables();

            size_ = o.size_;
            pop = o.pop;
            num_empty = o.num_empty;

            table_keys = o.table_keys;
            table_values = o.table_values;

            o.table_keys = nullptr;
            o.table_values = nullptr;

            o.clear_tables();
        }

        return *this;
    }

    ~KmerHashTable() {

        clear_tables();
    }

    // DUMMY TO DELETE
    void recomputeMaxPSL(const size_t nb_threads = 1) {
    }

    inline size_t size() const {

        return pop;
    }

    inline bool empty() const {

        return pop == 0;
    }

    void clear() {

        Kmer empty_key;

        empty_key.set_empty();

        std::fill(table_keys, table_keys + size_, empty_key);

        pop = 0;
        num_empty = size_;
    }

    void clear_tables() {

        if (table_keys != nullptr) {

            delete[] table_keys;
            table_keys = nullptr;
        }

        if (table_values != nullptr) {

            delete[] table_values;
            table_values = nullptr;
        }

        size_ = 0;
        pop  = 0;
        num_empty = 0;
    }

    void init_tables(const size_t sz) {

        clear_tables();

        size_ = rndup(sz);

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        clear();
    }

    void reserve(const size_t sz) {

        if (sz <= size_) return;

        const size_t old_size_ = size_;

        Kmer empty_key;

        Kmer* old_table_keys = table_keys;
        T* old_table_values = table_values;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table_keys = new Kmer[size_];
        table_values = new T[size_];

        empty_key.set_empty();

        std::fill(table_keys, table_keys + size_, empty_key);

        for (size_t i = 0; i < old_size_; ++i) {

            if (!old_table_keys[i].isEmpty() && !old_table_keys[i].isDeleted()){

                insert(std::move(old_table_keys[i]), std::move(old_table_values[i]));
            }
        }

        delete[] old_table_keys;
        delete[] old_table_values;
    }

    iterator find(const Kmer& key) {

        const size_t end_table = size_-1;

        size_t h = key.hash() & end_table;

        const size_t end_h = (h-1) & end_table;

        for (; h != end_h; h = (h+1) & end_table) {

            if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;
        }

        if ((h != end_h) && (table_keys[h] == key)) return iterator(this, h);

        return iterator(this);
    }

    const_iterator find(const Kmer& key) const {

        const size_t end_table = size_-1;

        size_t h = key.hash() & end_table;

        const size_t end_h = (h-1) & end_table;

        for (; h != end_h; h = (h+1) & end_table) {

            if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;
        }

        if ((h != end_h) && (table_keys[h] == key)) return const_iterator(this, h);

        return const_iterator(this);
    }

    iterator find(const size_t h) {

        if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return iterator(this, h);
        return iterator(this);
    }

    const_iterator find(const size_t h) const {

        if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return const_iterator(this, h);
        return const_iterator(this);
    }

    iterator erase(const_iterator pos) {

        if (pos == end()) return end();

        table_keys[pos.h].set_deleted();
        --pop;

        return ++iterator(this, pos.h); // return pointer to next element
    }

    size_t erase(const Kmer& minz) {

        const_iterator pos = find(minz);

        size_t oldpop = pop;

        if (pos != end()) erase(pos);

        return oldpop - pop;
    }

    std::pair<iterator, bool> insert(const Kmer& key, const T& value) {

        if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

        bool is_deleted = false;

        const size_t end_table = size_-1;

        for (size_t h = key.hash() & end_table, h_tmp;; h = (h+1) & end_table) {

            if (table_keys[h].isEmpty()) {

                is_deleted ? h = h_tmp : --num_empty;

                table_keys[h] = key;
                table_values[h] = value;

                ++pop;

                return {iterator(this, h), true};
            }
            else if (table_keys[h] == key) return {iterator(this, h), false};
            else if (!is_deleted && table_keys[h].isDeleted()) {

                is_deleted = true;
                h_tmp = h;
            }
        }
    }

    std::pair<iterator, bool> insert(Kmer&& key, T&& value) {

        if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

        bool is_deleted = false;

        const size_t end_table = size_-1;

        for (size_t h = key.hash() & end_table, h_tmp;; h = (h+1) & end_table) {

            if (table_keys[h].isEmpty()) {

                is_deleted ? h = h_tmp : --num_empty;

                table_keys[h] = key;
                table_values[h] = value;

                ++pop;

                return {iterator(this, h), true};
            }
            else if (table_keys[h] == key) return {iterator(this, h), false};
            else if (!is_deleted && table_keys[h].isDeleted()) {

                is_deleted = true;
                h_tmp = h;
            }
        }
    }

    iterator begin() {

        iterator it(this);
        it.find_first();
        return it;
    }

    const_iterator begin() const {

        const_iterator it(this);
        it.find_first();
        return it;
    }

    iterator end() {

        return iterator(this);
    }

    const_iterator end() const {

        return const_iterator(this);
    }
};*/

#endif
