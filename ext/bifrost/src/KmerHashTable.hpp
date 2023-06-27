#ifndef BIFROST_KMER_HASHTABLE_HPP
#define BIFROST_KMER_HASHTABLE_HPP

#include <utility>
#include <string>
#include <iterator>
#include <algorithm>

#include "Kmer.hpp"

template<typename T>
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
};

template<typename T>
struct MinimizerHashTable {

    size_t size_, pop, num_empty;

    Minimizer* table_keys;
    T* table_values;

    template<bool is_const = true>
    class iterator_ : public std::iterator<std::forward_iterator_tag, T> {

        public:

            typedef typename std::conditional<is_const, const MinimizerHashTable *, MinimizerHashTable *>::type MHT_ptr_t;
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

            const Minimizer& getKey() const { return ht->table_keys[h]; }

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

    MinimizerHashTable() : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        init_tables(1024);
    }

    MinimizerHashTable(const size_t sz) : table_keys(nullptr), table_values(nullptr), size_(0), pop(0), num_empty(0) {

        if (sz < 2) init_tables(2);
        else {

            const size_t sz_with_empty = static_cast<size_t>(1.2 * sz);

            size_t rdnup_sz = rndup(sz);

            while (rdnup_sz < sz_with_empty) rdnup_sz <<= 1;

            init_tables(max(rdnup_sz, static_cast<size_t>(2)));
        }
    }

    MinimizerHashTable(const MinimizerHashTable& o) :   size_(o.size_), pop(o.pop), num_empty(o.num_empty) {

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);
    }

    MinimizerHashTable(MinimizerHashTable&& o){

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = o.table_keys;
        table_values = o.table_values;

        o.table_keys = nullptr;
        o.table_values = nullptr;

        o.clear_tables();
    }

    MinimizerHashTable& operator=(const MinimizerHashTable& o) {

        clear_tables();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        std::copy(o.table_keys, o.table_keys + size_, table_keys);
        std::copy(o.table_values, o.table_values + size_, table_values);

        return *this;
    }

    MinimizerHashTable& operator=(MinimizerHashTable&& o){

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

    ~MinimizerHashTable() {

        clear_tables();
    }

    inline size_t size() const {

        return pop;
    }

    inline bool empty() const {

        return pop == 0;
    }

    void clear() {

        Minimizer empty_key;

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

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        clear();
    }

    void reserve(const size_t sz) {

        if (sz <= size_) return;

        const size_t old_size_ = size_;

        Minimizer empty_key;

        Minimizer* old_table_keys = table_keys;
        T* old_table_values = table_values;

        size_ = rndup(sz);
        pop = 0;
        num_empty = size_;

        table_keys = new Minimizer[size_];
        table_values = new T[size_];

        empty_key.set_empty();

        std::fill(table_keys, table_keys + size_, empty_key);

        for (size_t i = 0; i < old_size_; i++) {

            if (!old_table_keys[i].isEmpty() && !old_table_keys[i].isDeleted()){

                insert(std::move(old_table_keys[i]), std::move(old_table_values[i]));
            }
        }

        delete[] old_table_keys;
        delete[] old_table_values;
    }

    iterator find(const Minimizer& key) {

        const size_t end_table = size_-1;

        size_t h = key.hash() & end_table;

        const size_t end_h = (h-1) & end_table;

        for (; h != end_h; h = (h+1) & end_table) {

            if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;
        }

        if ((h != end_h) && (table_keys[h] == key)) return iterator(this, h);

        return iterator(this);
    }

    const_iterator find(const Minimizer& key) const {

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

    size_t erase(const Minimizer& minz) {

        const_iterator pos = find(minz);

        size_t oldpop = pop;

        if (pos != end()) erase(pos);

        return oldpop - pop;
    }

    std::pair<iterator, bool> insert(const Minimizer& key, const T& value) {

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

    std::pair<iterator, bool> insert(Minimizer&& key, T&& value) {

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
};

#endif
