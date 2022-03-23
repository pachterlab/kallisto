#include "MinimizerIndex.hpp"

MinimizerIndex::MinimizerIndex() :  table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                    size_(0), pop(0), num_empty(0)  {

    init_tables(max(static_cast<size_t>(1024), lck_block_sz));
}

MinimizerIndex::MinimizerIndex(const size_t sz) :   table_keys(nullptr), table_tinyv(nullptr), table_tinyv_sz(nullptr),
                                                    size_(0), pop(0), num_empty(0) {

    if (sz == 0) init_tables(lck_block_sz);
    else {

        const size_t sz_with_empty = static_cast<size_t>(1.2 * sz);

        size_t rdnup_sz = rndup(sz);

        while (rdnup_sz < sz_with_empty) rdnup_sz <<= 1;

        init_tables(max(rdnup_sz, lck_block_sz));
    }
}

MinimizerIndex::MinimizerIndex(const MinimizerIndex& o) :   size_(o.size_), pop(o.pop), num_empty(o.num_empty) {

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    lck_min = vector<SpinLock>(o.lck_min.size());

    std::copy(o.table_keys, o.table_keys + size_, table_keys);

    for (size_t i = 0; i < size_; ++i){

        table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
        table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
    }
}

MinimizerIndex::MinimizerIndex(MinimizerIndex&& o){

    size_ = o.size_;
    pop = o.pop;
    num_empty = o.num_empty;

    table_keys = o.table_keys;
    table_tinyv = o.table_tinyv;
    table_tinyv_sz = o.table_tinyv_sz;

    lck_min = vector<SpinLock>(o.lck_min.size());

    o.table_keys = nullptr;
    o.table_tinyv = nullptr;
    o.table_tinyv_sz = nullptr;

    o.clear();
}

MinimizerIndex& MinimizerIndex::operator=(const MinimizerIndex& o) {

    if (this != &o) {

        clear();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = new Minimizer[size_];
        table_tinyv = new packed_tiny_vector[size_];
        table_tinyv_sz = new uint8_t[size_];

        lck_min = vector<SpinLock>(o.lck_min.size());

        std::copy(o.table_keys, o.table_keys + size_, table_keys);

        for (size_t i = 0; i < size_; ++i){

            table_tinyv_sz[i] = packed_tiny_vector::FLAG_EMPTY;
            table_tinyv[i].copy(table_tinyv_sz[i], o.table_tinyv[i], o.table_tinyv_sz[i]);
        }
    }

    return *this;
}

MinimizerIndex& MinimizerIndex::operator=(MinimizerIndex&& o){

    if (this != &o) {

        clear();

        size_ = o.size_;
        pop = o.pop;
        num_empty = o.num_empty;

        table_keys = o.table_keys;
        table_tinyv = o.table_tinyv;
        table_tinyv_sz = o.table_tinyv_sz;

        lck_min = vector<SpinLock>(o.lck_min.size());

        o.table_keys = nullptr;
        o.table_tinyv = nullptr;
        o.table_tinyv_sz = nullptr;

        o.clear();
    }

    return *this;
}

MinimizerIndex::~MinimizerIndex() {

    clear();
}

void MinimizerIndex::clear() {

    if (table_tinyv != nullptr){

        for (size_t i = 0; i < size_; ++i) table_tinyv[i].destruct(table_tinyv_sz[i]);
    }

    clear_tables();

    lck_min.clear();
    lck_edit_table.release_all();
}

MinimizerIndex::iterator MinimizerIndex::find(const Minimizer& key) {

    const size_t end_table = size_-1;

    size_t h = key.hash() & end_table;
    size_t i = 0;

    while (i != size_) {

        if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == key)) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const Minimizer& key) const {

    const size_t end_table = size_-1;

    size_t h = key.hash() & end_table;
    size_t i = 0;

    while (i != size_) {

        if (table_keys[h].isEmpty() || (table_keys[h] == key)) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == key)) return const_iterator(this, h);

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::find(const size_t h) {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return iterator(this, h);

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find(const size_t h) const {

    if ((h < size_) && !table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return const_iterator(this, h);

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::erase(const_iterator it) {

    if (it == end()) return end();

    table_keys[it.h].set_deleted();
    table_tinyv[it.h].destruct(table_tinyv_sz[it.h]);
    table_tinyv_sz[it.h] = packed_tiny_vector::FLAG_EMPTY;

    --pop;

    return ++iterator(this, it.h); // return pointer to next element
}

size_t MinimizerIndex::erase(const Minimizer& minz) {

    const size_t end_table = size_-1;
    const size_t oldpop = pop;

    size_t h = minz.hash() & end_table;
    size_t i = 0;

    while (i != size_) {

        if (table_keys[h].isEmpty() || (table_keys[h] == minz)) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == minz)){

        table_keys[h].set_deleted();
        table_tinyv[h].destruct(table_tinyv_sz[h]);
        table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

        --pop;
    }

    return oldpop - pop;
}

// Insert with Robin Hood hashing
/*
pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

    const size_t end_table = size_-1;

    bool is_deleted = false, has_rich_psl = false;

    size_t h = key.hash() & end_table;
    size_t h_del = 0, h_rich_psl = 0;
    size_t psl_ins_key = 0, psl_rich_key = 0;

    pair<MinimizerIndex::iterator, bool> it_ret;

    Minimizer l_key = key;

    packed_tiny_vector l_ptv(ptv, flag);

    uint8_t l_flag = flag;

    while (true) {

        if (table_keys[h].isEmpty()) {

            if (has_rich_psl) {

                h = h_rich_psl;

                // Swap keys
                swap(table_keys[h], l_key);

                // Swap values
                packed_tiny_vector l_ptv_swap;
                uint8_t l_flag_swap;

                l_ptv_swap.move(l_flag_swap, move(table_tinyv[h]), move(table_tinyv_sz[h]));
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));
                l_ptv.move(l_flag, move(l_ptv_swap), move(l_flag_swap));

                psl_ins_key = psl_rich_key;

                is_deleted = false;
                has_rich_psl = false;

                if (table_keys[h] == key) it_ret = {iterator(this, h), true};
            }
            else {

                h = ((static_cast<size_t>(is_deleted) - 1) & h) + ((static_cast<size_t>(!is_deleted) - 1) & h_del);
                num_empty -= static_cast<size_t>(!is_deleted);

                table_keys[h] = l_key;
                table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;
                table_tinyv[h].move(table_tinyv_sz[h], move(l_ptv), move(l_flag));

                if (table_keys[h] == key) {

                    ++pop;
                    it_ret = {iterator(this, h), true};
                }

                return it_ret;
            }
        }
        else if (table_keys[h] == l_key) return {iterator(this, h), false}; // Can only happen when inserting the input key
        else if (table_keys[h].isDeleted()) {

            h_del = ((static_cast<size_t>(!is_deleted) - 1) & h_del) + ((static_cast<size_t>(is_deleted) - 1) & h);
            is_deleted = true;
        }
        else if (!is_deleted && !has_rich_psl) {

            const size_t h_curr = table_keys[h].hash() & end_table;
            const size_t psl_curr_key = (h < h_curr) ? (size_ - h_curr + h) : (h - h_curr);

            if (psl_ins_key > psl_curr_key) {

                has_rich_psl = true;
                h_rich_psl = h;
                psl_rich_key = psl_curr_key;
            }
        }

        h = (h+1) & end_table;
        ++psl_ins_key;
    }
}*/

pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert(const Minimizer& key, const packed_tiny_vector& ptv, const uint8_t& flag) {

    if ((5 * num_empty) < size_) reserve(2 * size_); // if more than 80% full, resize

    const size_t end_table = size_-1;
    
    size_t h = key.hash() & end_table, h_del;

    bool is_deleted = false;

    while (true) {

        if (table_keys[h].isEmpty()) {

            h = ((static_cast<size_t>(is_deleted) - 1) & h) + ((static_cast<size_t>(!is_deleted) - 1) & h_del);
            num_empty -= static_cast<size_t>(!is_deleted);

            table_keys[h] = key;
            table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

            table_tinyv[h].copy(table_tinyv_sz[h], ptv, flag);

            ++pop;

            return {iterator(this, h), true};
        }
        else if (table_keys[h] == key) return {iterator(this, h), false};
        else if (table_keys[h].isDeleted()) {

            h_del = ((static_cast<size_t>(!is_deleted) - 1) & h_del) + ((static_cast<size_t>(is_deleted) - 1) & h);
            is_deleted = true;
        }

        h = (h+1) & end_table;
    }
}

void MinimizerIndex::init_threads() {

    lck_min = vector<SpinLock>((size_ + lck_block_sz - 1) / lck_block_sz);

    pop_p = pop;
    num_empty_p = num_empty;
}

void MinimizerIndex::release_threads() {

    pop = pop_p;
    num_empty = num_empty_p;

    lck_min.clear();
    lck_edit_table.release_all();
}

MinimizerIndex::iterator MinimizerIndex::find_p(const Minimizer& key) {

    lck_edit_table.acquire_reader();

    const size_t end_table = size_-1;

    size_t i = 0;
    size_t h = key.hash() & end_table;
    size_t id_block = h >> lck_block_div_shift;

    lck_min[id_block].acquire();

    while (i != size_) {

        if ((h >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h].isEmpty()){

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            iterator(this);
        }
        else if (table_keys[h] == key) return iterator(this, h);

        h = (h+1) & end_table;
        ++i;
    }

    lck_min[id_block].release();
    lck_edit_table.release_reader();

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find_p(const Minimizer& key) const {

    lck_edit_table.acquire_reader();

    const size_t end_table = size_-1;

    size_t i = 0;
    size_t h = key.hash() & end_table;
    size_t id_block = h >> lck_block_div_shift;

    lck_min[id_block].acquire();

    while (i != size_) {

        if ((h >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h].isEmpty()){

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            const_iterator(this);
        }
        else if (table_keys[h] == key) return const_iterator(this, h);

        h = (h+1) & end_table;
        ++i;
    }

    lck_min[id_block].release();
    lck_edit_table.release_reader();

    return const_iterator(this);
}

MinimizerIndex::iterator MinimizerIndex::find_p(const size_t h) {

    lck_edit_table.acquire_reader();

    if (h < size_){

        const size_t id_block = h >> lck_block_div_shift;

        lck_min[id_block].acquire();

        if (!table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return iterator(this, h);

        lck_min[id_block].release();
    }

    lck_edit_table.release_reader();

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::find_p(const size_t h) const {

    lck_edit_table.acquire_reader();

    if (h < size_){

        const size_t id_block = h >> lck_block_div_shift;

        lck_min[id_block].acquire();

        if (!table_keys[h].isEmpty() && !table_keys[h].isDeleted()) return const_iterator(this, h);

        lck_min[id_block].release();
    }

    lck_edit_table.release_reader();

    return const_iterator(this);
}

void MinimizerIndex::release_p(const_iterator it) const {

    if (it != end()){

        lck_min[it.h >> lck_block_div_shift].release();
        lck_edit_table.release_reader();
    }
}

void MinimizerIndex::release_p(iterator it) {

    if (it != end()){

        lck_min[it.h >> lck_block_div_shift].release();
        lck_edit_table.release_reader();
    }
}

size_t MinimizerIndex::erase_p(const Minimizer& minz) {

    lck_edit_table.acquire_reader();

    const size_t end_table = size_ - 1;

    size_t i = 0;
    size_t h = minz.hash() & end_table;
    size_t id_block = h >> lck_block_div_shift;
    size_t l_pop = pop;

    lck_min[id_block].acquire();

    while (i != size_) {

        if ((h >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h].isEmpty()){

            lck_min[id_block].release();
            lck_edit_table.release_reader();

            return 0;
        }
        else if (table_keys[h] == minz) break;

        h = (h+1) & end_table;
        ++i;
    }

    if ((i != size_) && (table_keys[h] == minz)){

        table_keys[h].set_deleted();
        table_tinyv[h].destruct(table_tinyv_sz[h]);
        table_tinyv_sz[h] = packed_tiny_vector::FLAG_EMPTY;

        lck_min[id_block].release();

        --pop;
    }

    l_pop -= pop;

    lck_edit_table.release_reader();

    return l_pop;
}

pair<MinimizerIndex::iterator, bool> MinimizerIndex::insert_p(const Minimizer& key, const packed_tiny_vector& v, const uint8_t& flag) {

    bool is_deleted = false;

    lck_edit_table.acquire_reader();

    if ((5 * num_empty_p) < size_){

        lck_edit_table.release_reader();
        lck_edit_table.acquire_writer();

        reserve(2 * size_); // if more than 80% full, resize

        pop_p = pop;
        num_empty_p = num_empty;

        lck_edit_table.release_writer_acquire_reader();
    }

    const size_t end_table = size_-1;
    const size_t h = key.hash() & end_table;

    size_t id_block = h >> lck_block_div_shift;

    size_t h1 = h;
    size_t h2;

    lck_min[id_block].acquire();

    while (true) {

        if ((h1 >> lck_block_div_shift) != id_block){

            lck_min[id_block].release();
            id_block = h1 >> lck_block_div_shift;
            lck_min[id_block].acquire();
        }

        if (table_keys[h1].isEmpty()) {

            if (is_deleted){

                const size_t id_block2 = h2 >> lck_block_div_shift;

                if (id_block2 == id_block) h1 = h2;
                else {

                    lck_min[id_block2].acquire();

                    if (table_keys[h2].isDeleted()){

                        lck_min[id_block].release();
                        h1 = h2;
                    }
                    else {

                        lck_min[id_block2].release();
                        --num_empty_p;
                    }
                }
            }
            else --num_empty_p;

            table_keys[h1] = key;
            table_tinyv_sz[h1] = packed_tiny_vector::FLAG_EMPTY;

            table_tinyv[h1].copy(table_tinyv_sz[h1], v, flag);

            ++pop_p;

            return {iterator(this, h1), true};
        }
        else if (table_keys[h1] == key){

            return {iterator(this, h1), false};
        }
        else if (!is_deleted && table_keys[h1].isDeleted()) {

            is_deleted = true;
            h2 = h1;
        }

        h1 = (h1+1) & end_table;
    }

    lck_min[id_block].release(); // Just for safety
    lck_edit_table.release_reader(); // Just for safety
}

MinimizerIndex::iterator MinimizerIndex::begin() {

    iterator it(this);
    it.operator++();
    return it;
}

MinimizerIndex::const_iterator MinimizerIndex::begin() const {

    const_iterator it(this);
    it.operator++();
    return it;
}

MinimizerIndex::iterator MinimizerIndex::end() {

    return iterator(this);
}

MinimizerIndex::const_iterator MinimizerIndex::end() const {

    return const_iterator(this);
}

void MinimizerIndex::clear_tables() {

    if (table_keys != nullptr) {

        delete[] table_keys;
        table_keys = nullptr;
    }

    if (table_tinyv != nullptr) {

        delete[] table_tinyv;
        table_tinyv = nullptr;
    }

    if (table_tinyv_sz != nullptr) {

        delete[] table_tinyv_sz;
        table_tinyv_sz = nullptr;
    }

    size_ = 0;
    pop  = 0;
    num_empty = 0;
}

void MinimizerIndex::init_tables(const size_t sz) {

    clear_tables();

    Minimizer empty_key;

    pop = 0;
    size_ = rndup(sz);
    num_empty = size_;

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    empty_key.set_empty();

    std::fill(table_keys, table_keys + size_, empty_key);

    memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));
}

void MinimizerIndex::reserve(const size_t sz) {

    if (sz <= size_) return;

    const size_t old_size_ = size_;

    Minimizer empty_key;

    Minimizer* old_table_keys = table_keys;
    packed_tiny_vector* old_table_tinyv = table_tinyv;
    uint8_t* old_table_tinyv_sz = table_tinyv_sz;

    size_ = rndup(sz);
    pop = 0;
    num_empty = size_;

    table_keys = new Minimizer[size_];
    table_tinyv = new packed_tiny_vector[size_];
    table_tinyv_sz = new uint8_t[size_];

    empty_key.set_empty();

    if (!lck_min.empty()) lck_min = vector<SpinLock>((size_ + lck_block_sz - 1) / lck_block_sz);

    std::fill(table_keys, table_keys + size_, empty_key);

    memset(table_tinyv_sz, packed_tiny_vector::FLAG_EMPTY, size_ * sizeof(uint8_t));

    for (size_t i = 0; i < old_size_; ++i) {

        if (!old_table_keys[i].isEmpty() && !old_table_keys[i].isDeleted()){

            insert(old_table_keys[i], old_table_tinyv[i], old_table_tinyv_sz[i]);

            old_table_tinyv[i].destruct(old_table_tinyv_sz[i]);
        }
    }

    delete[] old_table_keys;
    delete[] old_table_tinyv;
    delete[] old_table_tinyv_sz;
}

const size_t MinimizerIndex::lck_block_sz = 64;
const size_t MinimizerIndex::lck_block_div_shift = 6;
