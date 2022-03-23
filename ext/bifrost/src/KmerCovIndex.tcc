template<typename T>
KmerCovIndex<T>::KmerCovIndex() : sz(0), shift_div(__builtin_ffsll(block_sz) - 1), mask_mod(block_sz - 1) {}

template<typename T>
KmerCovIndex<T>::KmerCovIndex(const KmerCovIndex& o) : sz(o.sz), shift_div(o.shift_div), mask_mod(o.mask_mod) {

    v_blocks = vector<Block<T>*>(o.v_blocks.size());

    for (size_t i = 0; i < v_blocks.size(); ++i) {

        v_blocks[i] = new Block<T>;
        v_blocks[i]->bc_cov = o.v_blocks[i]->bc_cov;

        std::copy(o.v_blocks[i]->km_block, o.v_blocks[i]->km_block + block_sz, v_blocks[i]->km_block);
        std::copy(o.v_blocks[i]->data_block, o.v_blocks[i]->data_block + block_sz, v_blocks[i]->data_block);
    }
}

template<>
inline KmerCovIndex<void>::KmerCovIndex(const KmerCovIndex& o) : sz(o.sz), shift_div(o.shift_div), mask_mod(o.mask_mod) {

    v_blocks = vector<Block<void>*>(o.v_blocks.size());

    for (size_t i = 0; i < v_blocks.size(); ++i) {

        v_blocks[i] = new Block<void>;
        v_blocks[i]->bc_cov = o.v_blocks[i]->bc_cov;

        std::copy(o.v_blocks[i]->km_block, o.v_blocks[i]->km_block + block_sz, v_blocks[i]->km_block);
    }
}

template<typename T>
KmerCovIndex<T>::KmerCovIndex(KmerCovIndex&& o) :  sz(o.sz), shift_div(o.shift_div), mask_mod(o.mask_mod), v_blocks(move(o.v_blocks)) {

    o.clear();
}

template<typename T>
KmerCovIndex<T>::~KmerCovIndex() {

    clear();
}

template<typename T>
KmerCovIndex<T>& KmerCovIndex<T>::operator=(const KmerCovIndex<T>& o) {

    if (this != &o) {

        clear();

        sz = o.sz;
        shift_div = o.shift_div;
        mask_mod = o.mask_mod;

        v_blocks = vector<Block<T>*>(o.v_blocks.size());

        for (size_t i = 0; i < v_blocks.size(); ++i) {

            v_blocks[i] = new Block<T>;
            v_blocks[i]->bc_cov = o.v_blocks[i]->bc_cov;

            std::copy(o.v_blocks[i]->km_block, o.v_blocks[i]->km_block + block_sz, v_blocks[i]->km_block);
            std::copy(o.v_blocks[i]->data_block, o.v_blocks[i]->data_block + block_sz, v_blocks[i]->data_block);
        }
    }

    return *this;
}

template<typename T>
KmerCovIndex<T>& KmerCovIndex<T>::toData(KmerCovIndex<void>&& o, const size_t nb_threads) {

    sz = o.sz;
    shift_div = o.shift_div;
    mask_mod = o.mask_mod;

    v_blocks = vector<Block<T>*>(o.v_blocks.size(), nullptr);

    auto copyBlock = [&](const size_t start, const size_t end){

        for (size_t i = start; i < end; ++i) {

            v_blocks[i] = new Block<T>;
            v_blocks[i]->bc_cov = move(o.v_blocks[i]->bc_cov);

            std::copy(o.v_blocks[i]->km_block, o.v_blocks[i]->km_block + block_sz, v_blocks[i]->km_block);

            delete o.v_blocks[i];

            o.v_blocks[i] = nullptr;
        }
    };

    if ((nb_threads == 1) || (v_blocks.size() < nb_threads)) copyBlock(0, v_blocks.size());
    else {

        vector<thread> workers;

        const size_t slice = (v_blocks.size() / nb_threads) + 1;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    const size_t start = t * slice;
                    const size_t end = min(start + slice, v_blocks.size());

                    if (start < v_blocks.size()) copyBlock(start, end);
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    o.clear();

    return *this;
}

template<>
inline KmerCovIndex<void>& KmerCovIndex<void>::operator=(const KmerCovIndex<void>& o) {

    if (this != &o) {

        clear();

        sz = o.sz;
        shift_div = o.shift_div;
        mask_mod = o.mask_mod;

        v_blocks = vector<Block<void>*>(o.v_blocks.size());

        for (size_t i = 0; i < v_blocks.size(); ++i) {

            v_blocks[i] = new Block<void>;
            v_blocks[i]->bc_cov = o.v_blocks[i]->bc_cov;

            std::copy(o.v_blocks[i]->km_block, o.v_blocks[i]->km_block + block_sz, v_blocks[i]->km_block);
        }
    }

    return *this;
}

template<typename T>
KmerCovIndex<T>& KmerCovIndex<T>::operator=(KmerCovIndex<T>&& o) {

    if (this != &o) {

        clear();

        sz = o.sz;
        shift_div = o.shift_div;
        mask_mod = o.mask_mod;

        v_blocks = std::move(o.v_blocks);

        o.clear();
    }

    return *this;
}

template<typename T>
void KmerCovIndex<T>::clear() {

    sz = 0;

    for (auto block : v_blocks) {

        if (block != nullptr) delete block;
    }

    v_blocks.clear();
}

template<typename T>
void KmerCovIndex<T>::push_back(const Kmer& km) {

    const size_t mod = sz & mask_mod;

    if (mod == 0) { // Current block is full

        v_blocks.push_back(nullptr);
        v_blocks.back() = new Block<T>;
    }

    v_blocks[sz >> shift_div]->km_block[mod] = km;

    ++sz;
}

template<typename T>
bool KmerCovIndex<T>::set(const size_t idx, const Kmer& km) {

    if (idx >= sz) return false;

    const int cov_ = covAt(idx);

    const size_t idx_mod = idx & mask_mod;

    Block<T>* block = v_blocks[idx >> shift_div];

    block->km_block[idx_mod] = km;

    if (cov_ != 0) {

        block->bc_cov.remove(idx_mod * cov_full + cov_ - 1);
        block->bc_cov.runOptimize();
    }

    return true;
}

template<typename T>
bool KmerCovIndex<T>::set(const size_t idx, const Kmer& km, const size_t cov) {

    if (idx >= sz) return false;

    const int cov_ = covAt(idx);

    const size_t idx_mod = idx & mask_mod;

    Block<T>* block = v_blocks[idx >> shift_div];

    block->km_block[idx_mod] = km;

    if (cov_ != cov) {

        if (cov_ != 0) block->bc_cov.remove(idx_mod * cov_full + cov_ - 1);
        if (cov != 0) block->bc_cov.add(idx_mod * cov_full + min(cov, static_cast<size_t>(2)) - 1);

        block->bc_cov.runOptimize();
    }

    return true;
}

template<typename T>
void KmerCovIndex<T>::setFull(const size_t idx) {

    if (idx < sz){

        Block<T>* block = v_blocks[idx >> shift_div];

        const size_t idx_mod = idx & mask_mod;
        const size_t pos = idx_mod * cov_full;
        const size_t pos_end = pos + cov_full;

        for (size_t i = pos; i < pos_end; ++i) block->bc_cov.remove(i);

        block->bc_cov.add(pos_end - 1);
        block->bc_cov.runOptimize();
    }
}

template<typename T>
int KmerCovIndex<T>::covAt(const size_t idx) const {

    if (idx < sz){

        Block<T>* block = v_blocks[idx >> shift_div];

        const size_t idx_mod = idx & mask_mod;
        const size_t pos = idx_mod * cov_full;
        const size_t pos_end = pos + cov_full;

        for (size_t i = pos; i < pos_end; ++i) {

            if (block->bc_cov.contains(i)) return i - pos + 1;
        }

        return 0;
    }

    return -1;
}

template<typename T>
void KmerCovIndex<T>::cover(const size_t idx) {

    if (idx < sz){

        const int cov = covAt(idx);

        if (cov != cov_full){

            Block<T>* block = v_blocks[idx >> shift_div];

            const size_t idx_mod = idx & mask_mod;

            if (cov != 0) block->bc_cov.remove(idx_mod * cov_full + cov - 1);

            block->bc_cov.add(idx_mod * cov_full + cov);
            block->bc_cov.runOptimize();
        }
    }
}

template<typename T>
void KmerCovIndex<T>::uncover(const size_t idx) {

    if (idx < sz){

        const int cov = covAt(idx);

        if (cov != 0){

            Block<T>* block = v_blocks[idx >> shift_div];

            const size_t idx_mod = idx & mask_mod;

            block->bc_cov.remove(idx_mod * cov_full + cov - 1);

            if (cov != 1) block->bc_cov.add(idx_mod * cov_full + cov - 2);

            block->bc_cov.runOptimize();
        }
    }
}

template<typename T>
bool KmerCovIndex<T>::swap(const size_t idx1, const size_t idx2) {

    if ((idx1 < sz) && (idx2 < sz)) {

        if (idx1 != idx2){

            const int cov1 = covAt(idx1), cov2 = covAt(idx2);

            const size_t idx1_mod = idx1 & mask_mod;
            const size_t idx2_mod = idx2 & mask_mod;

            Block<T>* block1 = v_blocks[idx1 >> shift_div];
            Block<T>* block2 = v_blocks[idx2 >> shift_div];

            std::swap(block1->km_block[idx1_mod], block2->km_block[idx2_mod]);
            std::swap(block1->data_block[idx1_mod], block2->data_block[idx2_mod]);

            if (cov1 != cov2){

                if (cov1 != 0) block1->bc_cov.remove(idx1_mod * cov_full + cov1 - 1);
                if (cov2 != 0) block2->bc_cov.remove(idx2_mod * cov_full + cov2 - 1);

                if (cov1 != 0) block2->bc_cov.add(idx2_mod * cov_full + cov1 - 1);
                if (cov2 != 0) block1->bc_cov.add(idx1_mod * cov_full + cov2 - 1);

                block1->bc_cov.runOptimize();
                block2->bc_cov.runOptimize();
            }
        }

        return true;
    }

    return false;
}

template<>
inline bool KmerCovIndex<void>::swap(const size_t idx1, const size_t idx2) {

    if ((idx1 < sz) && (idx2 < sz)) {

        if (idx1 != idx2){

            const int cov1 = covAt(idx1), cov2 = covAt(idx2);

            const size_t idx1_mod = idx1 & mask_mod;
            const size_t idx2_mod = idx2 & mask_mod;

            Block<void>* block1 = v_blocks[idx1 >> shift_div];
            Block<void>* block2 = v_blocks[idx2 >> shift_div];

            std::swap(block1->km_block[idx1_mod], block2->km_block[idx2_mod]);

            if (cov1 != cov2){

                if (cov1 != 0) block1->bc_cov.remove(idx1_mod * cov_full + cov1 - 1);
                if (cov2 != 0) block2->bc_cov.remove(idx2_mod * cov_full + cov2 - 1);

                if (cov1 != 0) block2->bc_cov.add(idx2_mod * cov_full + cov1 - 1);
                if (cov2 != 0) block1->bc_cov.add(idx1_mod * cov_full + cov2 - 1);

                block1->bc_cov.runOptimize();
                block2->bc_cov.runOptimize();
            }
        }

        return true;
    }

    return false;
}

template<typename T>
void KmerCovIndex<T>::resize(const size_t new_sz) {

    if (new_sz == 0) clear();
    else if (new_sz < sz){ // resize down

        Kmer km_empty;

        const size_t new_v_block_sz = (new_sz >> shift_div) + ((new_sz & mask_mod) != 0);
        const size_t rounded_sz = min(new_v_block_sz << shift_div, sz);
        const size_t nb_last_block = new_sz & mask_mod;

        for (size_t i = new_v_block_sz; i < v_blocks.size(); ++i) {

            if (v_blocks[i] != nullptr) delete v_blocks[i];
        }

        v_blocks.resize(new_v_block_sz);

        Block<T>* block = v_blocks.back();

        if (nb_last_block != 0){

            for (size_t i = nb_last_block; i < block_sz; ++i) block->data_block[i] = T();
        }

        for (size_t i = new_sz; i < rounded_sz; ++i) {

            const int cov = covAt(i);

            if (cov > 0) block->bc_cov.remove((i & mask_mod) * cov_full + cov - 1);
        }

        block->bc_cov.runOptimize();

        sz = new_sz;
    }
    else if (new_sz > sz){

        Kmer km_empty;

        const size_t old_v_km_sz = v_blocks.size();
        const size_t new_v_block_sz = (new_sz >> shift_div) + ((new_sz & mask_mod) != 0);
        const size_t nb_last_block = sz & mask_mod;

        km_empty.set_empty();

        if (nb_last_block != 0) {

            Block<T>* block = v_blocks.back();

            std::fill(block->km_block + nb_last_block, block->km_block + block_sz, km_empty);

            for (size_t i = nb_last_block; i < block_sz; ++i) block->data_block[i] = T();
        }

        v_blocks.resize(new_v_block_sz);

        for (size_t i = old_v_km_sz; i < v_blocks.size(); ++i) {

            v_blocks[i] = new Block<T>;

            std::fill(v_blocks[i]->km_block, v_blocks[i]->km_block + block_sz, km_empty);
        }

        sz = new_sz;
    }
}

template<>
inline void KmerCovIndex<void>::resize(const size_t new_sz) {

    if (new_sz == 0) clear();
    else if (new_sz < sz){ // resize down

        Kmer km_empty;

        const size_t new_v_block_sz = (new_sz >> shift_div) + ((new_sz & mask_mod) != 0);
        const size_t rounded_sz = min(new_v_block_sz << shift_div, sz);

        for (size_t i = new_v_block_sz; i < v_blocks.size(); ++i) {

            if (v_blocks[i] != nullptr) delete v_blocks[i];
        }

        v_blocks.resize(new_v_block_sz);

        Block<void>* block = v_blocks.back();

        for (size_t i = new_sz; i < rounded_sz; ++i) {

            const int cov = covAt(i);

            if (cov > 0) block->bc_cov.remove((i & mask_mod) * cov_full + cov - 1);
        }

        block->bc_cov.runOptimize();

        sz = new_sz;
    }
    else if (new_sz > sz){

        Kmer km_empty;

        const size_t old_v_km_sz = v_blocks.size();
        const size_t new_v_block_sz = (new_sz >> shift_div) + ((new_sz & mask_mod) != 0);
        const size_t nb_last_block = sz & mask_mod;

        km_empty.set_empty();

        if (nb_last_block != 0) std::fill(v_blocks.back()->km_block + nb_last_block, v_blocks.back()->km_block + block_sz, km_empty);

        v_blocks.resize(new_v_block_sz);

        for (size_t i = old_v_km_sz; i < v_blocks.size(); ++i) {

            v_blocks[i] = new Block<void>;

            std::fill(v_blocks[i]->km_block, v_blocks[i]->km_block + block_sz, km_empty);
        }

        sz = new_sz;
    }
}

template<typename T>
const T* KmerCovIndex<T>::getData(const size_t idx) const {

    if (idx < sz) return &(v_blocks[idx >> shift_div]->data_block[idx & mask_mod]);

    return nullptr;
}

template<>
inline const void* KmerCovIndex<void>::getData(const size_t idx) const {

    return nullptr;
}

template<typename T>
T* KmerCovIndex<T>::getData(const size_t idx) {

    if (idx < sz) return &(v_blocks[idx >> shift_div]->data_block[idx & mask_mod]);

    return nullptr;
}

template<>
inline void* KmerCovIndex<void>::getData(const size_t idx) {

    return nullptr;
}

template<typename T>
Kmer KmerCovIndex<T>::getKmer(const size_t idx) const {

    if (idx < sz) return v_blocks[idx >> shift_div]->km_block[idx & mask_mod];

    Kmer empty_km;

    empty_km.set_empty();

    return empty_km;
}

template<typename T>
void KmerCovIndex<T>::remove(const size_t idx) {

    if (idx < sz){

        const size_t idx_mod = idx & mask_mod;

        Block<T>* block = v_blocks[idx >> shift_div];

        block->km_block[idx_mod].set_deleted();
        block->data_block[idx_mod] = T();

        const int cov = covAt(idx);

        if (cov != 0) {

            block->bc_cov.remove(idx_mod * cov_full + cov - 1);
            block->bc_cov.runOptimize();
        }
    }
}

template<>
inline void KmerCovIndex<void>::remove(const size_t idx) {

    if (idx < sz){

        const size_t idx_mod = idx & mask_mod;

        Block<void>* block = v_blocks[idx >> shift_div];

        block->km_block[idx_mod].set_deleted();

        const int cov = covAt(idx);

        if (cov != 0) {

            block->bc_cov.remove(idx_mod * cov_full + cov - 1);
            block->bc_cov.runOptimize();
        }
    }
}
