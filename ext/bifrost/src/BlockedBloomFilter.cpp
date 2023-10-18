#include "BlockedBloomFilter.hpp"

BlockedBloomFilter::BlockedBloomFilter() : table_(nullptr) {

    clear();
}

BlockedBloomFilter::BlockedBloomFilter(const size_t nb_elem, const size_t bits_per_elem) : table_(nullptr) {

    initialize(nb_elem, bits_per_elem);
}

BlockedBloomFilter::BlockedBloomFilter(const BlockedBloomFilter& o) :   table_(nullptr), blocks_(o.blocks_), nb_bits_per_elem(o.nb_bits_per_elem),
                                                                        k_(o.k_), M_u64(o.M_u64), seed1(o.seed1), seed2(o.seed2), ush(o.ush) {

    if (blocks_ != 0) {

        init_arrays();

        for (uint64_t i = 0; i != blocks_; ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), BBF_NB_ELEM_BLOCK * sizeof(uint64_t));

            table_[i].bits_occupancy = o.table_[i].bits_occupancy;
        }
    }
}

BlockedBloomFilter::BlockedBloomFilter(BlockedBloomFilter&& o) :    table_(o.table_), blocks_(o.blocks_), nb_bits_per_elem(o.nb_bits_per_elem),
                                                                    k_(o.k_), M_u64(o.M_u64), seed1(o.seed1), seed2(o.seed2), ush(move(o.ush)) {

    o.table_ = nullptr;

    o.clear();
}

BlockedBloomFilter::~BlockedBloomFilter() {

    clear();
}

BlockedBloomFilter& BlockedBloomFilter::operator=(const BlockedBloomFilter& o) {

    clear();

    blocks_ = o.blocks_;
    k_ = o.k_;
    nb_bits_per_elem = o.nb_bits_per_elem;
    M_u64 = o.M_u64;
    seed1 = o.seed1;
    seed2 = o.seed2;
    ush = o.ush;

    if (blocks_ != 0) {

        init_arrays();

        for (uint64_t i = 0; i != blocks_; ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), BBF_NB_ELEM_BLOCK * sizeof(uint64_t));

            table_[i].bits_occupancy = o.table_[i].bits_occupancy;
        }
    }

    return *this;
}

BlockedBloomFilter& BlockedBloomFilter::operator=(BlockedBloomFilter&& o) {

    if (this != &o) {

        clear();

        table_ = o.table_;
        blocks_ = o.blocks_;
        k_ = o.k_;
        nb_bits_per_elem = o.nb_bits_per_elem;
        M_u64 = o.M_u64;
        seed1 = o.seed1;
        seed2 = o.seed2;

        ush = move(o.ush);

        o.table_ = nullptr;

        o.clear();
    }

    return *this;
}

void BlockedBloomFilter::initialize(const size_t nb_elem, const size_t bits_per_elem) {

    clear();

    if ((nb_elem != 0) && (bits_per_elem != 0)){

        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<unsigned long long> distribution;

        blocks_ = (bits_per_elem * nb_elem + BBF_MASK_BITS_BLOCK) / BBF_NB_BITS_BLOCK;
        k_ = static_cast<int>(bits_per_elem * log(2));
        nb_bits_per_elem = bits_per_elem;

        if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) ++k_;

        seed1 = distribution(gen);
        seed2 = distribution(gen);

        init_arrays();
    }    
}

void BlockedBloomFilter::clear() {

    if (table_ != nullptr){

        delete[] table_;
        table_ = nullptr;
    }

    blocks_ = 0;
    k_ = 0;
    seed1 = 0;
    seed2 = 0;
    nb_bits_per_elem = 0;
    M_u64 = 0;

    ush.clear();

    lck_ush.clear(std::memory_order_release);
}

void BlockedBloomFilter::reset() {

    for (size_t i = 0; (table_ != nullptr) && (i < blocks_); ++i) table_[i].clear();

    ush.clear();
    lck_ush.clear(std::memory_order_release);
}

void BlockedBloomFilter::init_arrays() {

    M_u64 = fastmod::computeM_u64(blocks_);
    table_ = new BBF_Block[blocks_];
}

DualBlockedBloomFilter BlockedBloomFilter::transferToDBBF(const uint64_t idx_bbf) {

    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    DualBlockedBloomFilter dbbf;

    dbbf.blocks_ = blocks_;
    dbbf.k_ = k_;
    dbbf.nb_bits_per_elem = nb_bits_per_elem;
    dbbf.M_u64 = M_u64;
    dbbf.seed1 = seed1;
    dbbf.seed2 = seed2;

    dbbf.ush[idx_bbf_norm] = move(ush);
    dbbf.table_ = new DualBlockedBloomFilter::BBF_Block[blocks_ << 1];

    for (size_t i = 0; i < blocks_; ++i) {

            memcpy(&(dbbf.table_[(i << 1) + idx_bbf_norm].block), &(table_[i].block), BBF_NB_ELEM_BLOCK * sizeof(uint64_t));

            dbbf.table_[(i << 1) + idx_bbf_norm].bits_occupancy = table_[i].bits_occupancy;
    }

    clear();

    return dbbf;
}

int BlockedBloomFilter::contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit) const {

    const std::array<int64_t, 4> bids = contains_bids(kmh, minh, limit);

    pres[0] = (bids[0] != -1);
    pres[1] = (bids[1] != -1);
    pres[2] = (bids[2] != -1);
    pres[3] = (bids[3] != -1);

    return (static_cast<int>(pres[0]) + static_cast<int>(pres[1]) + static_cast<int>(pres[2]) + static_cast<int>(pres[3]));
}

std::array<int64_t, 4> BlockedBloomFilter::contains_bids(const uint64_t (&kmh)[4], const uint64_t minh, const int limit) const {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    /*const uint64_t kmh_s1[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed1, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed1, _wyp)};*/
    const uint64_t kmh_s1[4] = {kmh[0], kmh[1], kmh[2], kmh[3]};

    const uint64_t kmh_s2[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed2, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed2, _wyp)};

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    std::array<int64_t, 4> bids = {-1, -1, -1, -1};

    int cpt = 0;

    while (nb_overflow < 8) {

        uint64_t bid1 = fastmod::fastmod_u64(minh1, M_u64, blocks_);

        const uint64_t bits_occupancy_1 = table_[bid1].bits_occupancy;

        for (uint64_t i = 0, j = 0; (j != 4) && (cpt != limit); ++j){

            if (bids[j] == -1){

                uint64_t kmh1 = kmh_s1[j];

                for (i = 0; i != k_; ++i, kmh1 += kmh_s2[j]) {

                    if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
                }

                bids[j] = static_cast<int64_t>((bid1 + 2) & (static_cast<uint64_t>(i != k_) - 1)) - 1;
                cpt += (i == k_);
            }
        }

        if (cpt != limit){

            minh1 += minh_s2;
            bid1 = fastmod::fastmod_u64(minh1, M_u64, blocks_);

            const uint64_t bits_occupancy_2 = table_[bid1].bits_occupancy;

            for (size_t i = 0, j = 0; (j != 4) && (cpt != limit); ++j){

                if (bids[j] == -1){

                    uint64_t kmh1 = kmh_s1[j];

                    for (i = 0; i != k_; ++i, kmh1 += kmh_s2[j]) {

                        if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
                    }

                    bids[j] = static_cast<int64_t>((bid1 + 2) & (static_cast<uint64_t>(i != k_) - 1)) - 1;
                    cpt += (i == k_);
                }
            }

            if ((bits_occupancy_1 < min_overload_bits) || (bits_occupancy_2 < min_overload_bits) || (cpt == limit)) break;

            minh1 += minh_s2;
            ++nb_overflow;
        }
        else break;
    }

    if (nb_overflow >= 8) {

        for (size_t j = 0; (j != 4) && (cpt != limit); ++j) {

            if ((bids[j] == -1) && (ush.find(kmh[j]) != ush.end())) {

                bids[j] = 0;
                ++cpt;
            }
        }
    }

    return bids;
}

int64_t BlockedBloomFilter::contains_bids(const uint64_t kmh, const uint64_t minh) const {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = /*wyhash(&kmh, sizeof(uint64_t), seed1, _wyp)*/ kmh;
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;
    uint64_t kmh1, bid1;

    int i = 0;

    while (nb_overflow < 8) {

        kmh1 = kmh_s1;

        uint64_t bid1 = fastmod::fastmod_u64(minh1, M_u64, blocks_);

        const uint64_t bits_occupancy_1 = table_[bid1].bits_occupancy;

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            minh1 += minh_s2;
            bid1 = fastmod::fastmod_u64(minh1, M_u64, blocks_);
            kmh1 = kmh_s1;

            const uint64_t bits_occupancy_2 = table_[bid1].bits_occupancy;

            const bool overloaded = (bits_occupancy_1 >= min_overload_bits) && (bits_occupancy_2 >= min_overload_bits);

            for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

                if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) {

                    if (overloaded) break;
                    
                    return -1;
                }
            }

            if (i == k_) return static_cast<int64_t>(bid1 + 1);

            minh1 += minh_s2;
            ++nb_overflow;
        }
        else return static_cast<int64_t>(bid1 + 1);
    }

    if ((nb_overflow >= 8) && (ush.find(kmh) != ush.end())) return 0;

    return -1;
}

bool BlockedBloomFilter::write(FILE *fp) const {

    if (fwrite(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&nb_bits_per_elem, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&k_, sizeof(int), 1, fp) != 1) return false;

    const uint64_t ush_sz = ush.size();

    if (fwrite(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (const uint64_t h : ush) {

        if (fwrite(&h, sizeof(uint64_t), 1, fp) != 1) return false;
    }

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fwrite(&(table_[i].block), sizeof(uint64_t), BBF_NB_ELEM_BLOCK, fp) != BBF_NB_ELEM_BLOCK) return false;
        if (fwrite(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

bool BlockedBloomFilter::read(FILE *fp) {

    uint64_t ush_sz = 0;

    clear();

    if (fread(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&nb_bits_per_elem, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&k_, sizeof(int), 1, fp) != 1) return false;

    if (fread(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (uint64_t i = 0, h = 0; i != ush_sz; ++i){

        if (fread(&h, sizeof(uint64_t), 1, fp) != 1) return false;

        ush.insert(h);
    }

    init_arrays();

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fread(&(table_[i].block), sizeof(uint64_t), BBF_NB_ELEM_BLOCK, fp) != BBF_NB_ELEM_BLOCK) return false;
        if (fread(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

uint64_t BlockedBloomFilter::insert_par(const uint64_t kmh, const uint64_t minh) {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = /*wyhash(&kmh, sizeof(uint64_t), seed1, _wyp)*/ kmh;
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int i = 0, j = 0;

    while (i != k_) {

        uint64_t kmh1 = kmh_s1, kmh2 = kmh_s1;
        uint64_t minh2 = minh1 + minh_s2;

        uint64_t bid1 = fastmod::fastmod_u64(minh1, M_u64, blocks_);
        uint64_t bid2 = fastmod::fastmod_u64(minh2, M_u64, blocks_);

        if (bid2 < bid1) std::swap(bid1, bid2);

        table_[bid1].lock();

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            if (bid2 == bid1){

                j = i;
                kmh2 = kmh1;
            }
            else {

                table_[bid2].lock();

                for (j = 0; j != k_; ++j, kmh2 += kmh_s2) {

                    if ((table_[bid2].block[(kmh2 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh2 & 0x3fULL))) == 0) break;
                }
            }

            if (j != k_){

                if ((table_[bid1].bits_occupancy < min_overload_bits) || (table_[bid2].bits_occupancy < min_overload_bits)) {

                    uint64_t nb_inserted_bits = 0;

                    if (table_[bid2].bits_occupancy < table_[bid1].bits_occupancy){

                        i = j;
                        kmh1 = kmh2;

                        std::swap(bid1, bid2);
                    }

                    for (; i != k_; ++i, kmh1 += kmh_s2) {

                        const uint64_t div = (kmh1 & BBF_MASK_BITS_BLOCK) >> 6;
                        const uint64_t mod = 1ULL << (kmh1 & 0x3fULL);

                        nb_inserted_bits += static_cast<uint64_t>((table_[bid1].block[div] & mod) == 0);
                        table_[bid1].block[div] |= mod;
                    }

                    table_[bid1].bits_occupancy += nb_inserted_bits;

                    if (bid2 != bid1) table_[bid2].unlock();

                    table_[bid1].unlock();

                    return (((bid1 + 1) << 1) | 0x1ULL);
                }
                else if (nb_overflow == 7) {

                    bool inserted = false;

                    if (bid2 != bid1) table_[bid2].unlock();

                    table_[bid1].unlock();

                    while (lck_ush.test_and_set(std::memory_order_acquire));

                    inserted = ush.insert(kmh).second;

                    lck_ush.clear(std::memory_order_release);

                    return static_cast<uint64_t>(inserted);
                }
            }
            else i = j;

            if (bid2 != bid1) table_[bid2].unlock();
        }

        table_[bid1].unlock();

        minh1 += minh_s2 + minh_s2;
        ++nb_overflow;
    }

    return 0;
}

uint64_t BlockedBloomFilter::insert_unpar(const uint64_t kmh, const uint64_t minh) {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = /*wyhash(&kmh, sizeof(uint64_t), seed1, _wyp)*/ kmh;
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int i = 0, j = 0;

    while (i != k_) {

        uint64_t kmh1 = kmh_s1, kmh2 = kmh_s1;
        uint64_t minh2 = minh1 + minh_s2;

        uint64_t bid1 = fastmod::fastmod_u64(minh1, M_u64, blocks_);
        uint64_t bid2 = fastmod::fastmod_u64(minh2, M_u64, blocks_);

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            if (bid2 == bid1){

                j = i;
                kmh2 = kmh1;
            }
            else {

    		    for (j = 0; j != k_; ++j, kmh2 += kmh_s2) {

    		        if ((table_[bid2].block[(kmh2 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh2 & 0x3fULL))) == 0) break;
    		    }
	       }

            if (j != k_){

                if ((table_[bid1].bits_occupancy < min_overload_bits) || (table_[bid2].bits_occupancy < min_overload_bits)) {

                    uint64_t nb_inserted_bits = 0;

                    if (table_[bid2].bits_occupancy < table_[bid1].bits_occupancy){

                        i = j;
                        kmh1 = kmh2;

                        std::swap(bid1, bid2);
                    }

                    for (; i != k_; ++i, kmh1 += kmh_s2) {

                        const uint64_t div = (kmh1 & BBF_MASK_BITS_BLOCK) >> 6;
                        const uint64_t mod = 1ULL << (kmh1 & 0x3fULL);

                        nb_inserted_bits += static_cast<uint64_t>((table_[bid1].block[div] & mod) == 0);
                        table_[bid1].block[div] |= mod;
                    }

                    table_[bid1].bits_occupancy += nb_inserted_bits;

                    return (((bid1 + 1) << 1) | 0x1ULL);
                }
                else if (nb_overflow == 7) {

                    const bool inserted = ush.insert(kmh).second;

                    return static_cast<uint64_t>(inserted);
                }
            }
            else i = j;
        }

        minh1 += minh_s2 + minh_s2;

        ++nb_overflow;
    }

    return 0;
}

DualBlockedBloomFilter::DualBlockedBloomFilter() : table_(nullptr) {

    clear();
}

DualBlockedBloomFilter::DualBlockedBloomFilter(const size_t nb_elem, const size_t bits_per_elem) : table_(nullptr) {

    initialize(nb_elem, bits_per_elem);
}

DualBlockedBloomFilter::DualBlockedBloomFilter(const DualBlockedBloomFilter& o) :   table_(nullptr), blocks_(o.blocks_), nb_bits_per_elem(o.nb_bits_per_elem),
                                                                                    k_(o.k_), M_u64(o.M_u64), seed1(o.seed1), seed2(o.seed2), ush{o.ush[0], o.ush[1]} {

    if (blocks_ != 0) {

        init_arrays();

        for (uint64_t i = 0; i != (blocks_ << 1); ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), BBF_NB_ELEM_BLOCK * sizeof(uint64_t));

            table_[i].bits_occupancy = o.table_[i].bits_occupancy;
        }
    }
}

DualBlockedBloomFilter::DualBlockedBloomFilter(DualBlockedBloomFilter&& o) :    table_(o.table_), blocks_(o.blocks_), nb_bits_per_elem(o.nb_bits_per_elem),
                                                                                k_(o.k_), M_u64(o.M_u64), seed1(o.seed1), seed2(o.seed2), ush{move(o.ush[0]), move(o.ush[1])} {

    o.table_ = nullptr;

    o.clear();
}

DualBlockedBloomFilter::~DualBlockedBloomFilter() {

    clear();
}

DualBlockedBloomFilter& DualBlockedBloomFilter::operator=(const DualBlockedBloomFilter& o) {

    clear();

    blocks_ = o.blocks_;
    nb_bits_per_elem = o.nb_bits_per_elem;
    k_ = o.k_;
    M_u64 = o.M_u64;
    seed1 = o.seed1;
    seed2 = o.seed2;

    ush[0] = o.ush[0];
    ush[1] = o.ush[1];

    if (blocks_ != 0) {

        init_arrays();

        for (uint64_t i = 0; i != (blocks_ << 1); ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), BBF_NB_ELEM_BLOCK * sizeof(uint64_t));

            table_[i].bits_occupancy = o.table_[i].bits_occupancy;
        }
    }

    return *this;
}

DualBlockedBloomFilter& DualBlockedBloomFilter::operator=(DualBlockedBloomFilter&& o) {

    if (this != &o) {

        clear();

        table_ = o.table_;
        blocks_ = o.blocks_;
        nb_bits_per_elem = o.nb_bits_per_elem;
        k_ = o.k_;
        M_u64 = o.M_u64;
        seed1 = o.seed1;
        seed2 = o.seed2;

        ush[0] = move(o.ush[0]);
        ush[1] = move(o.ush[1]);

        o.table_ = nullptr;

        o.clear();
    }

    return *this;
}

void DualBlockedBloomFilter::initialize(const size_t nb_elem, const size_t bits_per_elem) {

    clear();

    if ((nb_elem != 0) && (bits_per_elem != 0)){

        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<unsigned long long> distribution;

        blocks_ = (bits_per_elem * nb_elem + BBF_MASK_BITS_BLOCK) / BBF_NB_BITS_BLOCK;
        nb_bits_per_elem = bits_per_elem;
        k_ = static_cast<int>(bits_per_elem * log(2));

        if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) ++k_;

        seed1 = distribution(gen);
        seed2 = distribution(gen);

        init_arrays();
    }    
}

void DualBlockedBloomFilter::clear() {

    if (table_ != nullptr){

        delete[] table_;
        table_ = nullptr;
    }

    blocks_ = 0;
    nb_bits_per_elem = 0;
    k_ = 0;
    seed1 = 0;
    seed2 = 0;
    M_u64 = 0;

    ush[0].clear();
    ush[1].clear();

    lck_ush[0].clear(std::memory_order_release);
    lck_ush[1].clear(std::memory_order_release);
}

void DualBlockedBloomFilter::reset() {

    for (size_t i = 0; (table_ != nullptr) && (i < (blocks_ << 1)); ++i) table_[i].clear();

    ush[0].clear();
    ush[1].clear();

    lck_ush[0].clear(std::memory_order_release);
    lck_ush[1].clear(std::memory_order_release);
}

void DualBlockedBloomFilter::init_arrays() {

    M_u64 = fastmod::computeM_u64(blocks_);
    table_ = new BBF_Block[blocks_ << 1];
}

int DualBlockedBloomFilter::contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit, const uint64_t idx_bbf) const {

    const std::array<int64_t, 4> bids = contains_bids(kmh, minh, limit, idx_bbf);

    pres[0] = (bids[0] != -1);
    pres[1] = (bids[1] != -1);
    pres[2] = (bids[2] != -1);
    pres[3] = (bids[3] != -1);

    return (static_cast<int>(pres[0]) + static_cast<int>(pres[1]) + static_cast<int>(pres[2]) + static_cast<int>(pres[3]));
}

std::array<int64_t, 4> DualBlockedBloomFilter::contains_bids(const uint64_t (&kmh)[4], const uint64_t minh, const int limit, const uint64_t idx_bbf) const {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    /*const uint64_t kmh_s1[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed1, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed1, _wyp)};*/
    const uint64_t kmh_s1[4] = {kmh[0], kmh[1], kmh[2], kmh[3]};

    const uint64_t kmh_s2[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed2, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed2, _wyp)};

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;
    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    std::array<int64_t, 4> bids = {-1, -1, -1, -1};

    int cpt = 0;

    while (nb_overflow < 8) {

        uint64_t bid1 = (fastmod::fastmod_u64(minh1, M_u64, blocks_) << 1) + idx_bbf_norm;

        const uint64_t bits_occupancy_1 = table_[bid1].bits_occupancy;

        for (uint64_t i = 0, j = 0; (j != 4) && (cpt != limit); ++j){

            if (bids[j] == -1){

                uint64_t kmh1 = kmh_s1[j];

                for (i = 0; i != k_; ++i, kmh1 += kmh_s2[j]) {

                    if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
                }

                bids[j] = static_cast<int64_t>((bid1 + 2) & (static_cast<uint64_t>(i != k_) - 1)) - 1;
                cpt += (i == k_);
            }
        }

        if (cpt != limit){

            minh1 += minh_s2;
            bid1 = (fastmod::fastmod_u64(minh1, M_u64, blocks_) << 1) + idx_bbf_norm;

            const uint64_t bits_occupancy_2 = table_[bid1].bits_occupancy;

            for (size_t i = 0, j = 0; (j != 4) && (cpt != limit); ++j){

                if (bids[j] == -1){

                    uint64_t kmh1 = kmh_s1[j];

                    for (i = 0; i != k_; ++i, kmh1 += kmh_s2[j]) {

                        if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
                    }

                    bids[j] = static_cast<int64_t>((bid1 + 2) & (static_cast<uint64_t>(i != k_) - 1)) - 1;
                    cpt += (i == k_);
                }
            }

            if ((bits_occupancy_1 < min_overload_bits) || (bits_occupancy_2 < min_overload_bits) || (cpt == limit)) break;

            minh1 += minh_s2;
            ++nb_overflow;
        }
        else break;
    }

    if (nb_overflow >= 8) {

        for (size_t j = 0; (j != 4) && (cpt != limit); ++j) {

            if ((bids[j] == -1) && (ush[idx_bbf_norm].find(kmh[j]) != ush[idx_bbf_norm].end())) {

                bids[j] = 0;
                ++cpt;
            }
        }
    }

    return bids;
}

int64_t DualBlockedBloomFilter::contains_bids(const uint64_t kmh, const uint64_t minh, const uint64_t idx_bbf) const {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = /*wyhash(&kmh, sizeof(uint64_t), seed1, _wyp)*/ kmh;
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;
    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;
    uint64_t kmh1, bid1;

    int i = 0;

    while (nb_overflow < 8) {

        kmh1 = kmh_s1;

        uint64_t bid1 = (fastmod::fastmod_u64(minh1, M_u64, blocks_) << 1) + idx_bbf_norm;

        const uint64_t bits_occupancy_1 = table_[bid1].bits_occupancy;

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            minh1 += minh_s2;
            bid1 = (fastmod::fastmod_u64(minh1, M_u64, blocks_) << 1) + idx_bbf_norm;
            kmh1 = kmh_s1;

            const uint64_t bits_occupancy_2 = table_[bid1].bits_occupancy;

            const bool overloaded = (bits_occupancy_1 >= min_overload_bits) && (bits_occupancy_2 >= min_overload_bits);

            for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

                if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) {

                    if (overloaded) break;
                    
                    return -1;
                }
            }

            if (i == k_) return static_cast<int64_t>(bid1 + 1);

            minh1 += minh_s2;
            ++nb_overflow;
        }
        else return static_cast<int64_t>(bid1 + 1);
    }

    if ((nb_overflow >= 8) && (ush[idx_bbf_norm].find(kmh) != ush[idx_bbf_norm].end())) return 0;

    return -1;
}

BlockedBloomFilter DualBlockedBloomFilter::transferToBBF(const uint64_t idx_bbf) {

    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    BlockedBloomFilter bbf;

    bbf.blocks_ = blocks_;
    bbf.k_ = k_;
    bbf.seed1 = seed1;
    bbf.seed2 = seed2;
    bbf.M_u64 = M_u64;

    bbf.ush = move(ush[idx_bbf_norm]);
    bbf.table_ = new BlockedBloomFilter::BBF_Block[blocks_];

    for (size_t i = 0; i < blocks_; ++i) {

            memcpy(&(bbf.table_[i].block), &(table_[(i << 1) + idx_bbf_norm].block), BBF_NB_ELEM_BLOCK * sizeof(uint64_t));

            bbf.table_[i].bits_occupancy = table_[(i << 1) + idx_bbf_norm].bits_occupancy;
    }

    clear();

    return bbf;
}

bool DualBlockedBloomFilter::write(FILE *fp) const {

    if (fwrite(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&nb_bits_per_elem, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&k_, sizeof(int), 1, fp) != 1) return false;

    for (size_t i = 0; i < 2; ++i) {

        const uint64_t ush_sz = ush[i].size();

        if (fwrite(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

        for (const uint64_t h : ush[i]) {

            if (fwrite(&h, sizeof(uint64_t), 1, fp) != 1) return false;
        }
    }

    for (uint64_t i = 0; i != (blocks_ << 1); ++i){

        if (fwrite(&(table_[i].block), sizeof(uint64_t), BBF_NB_ELEM_BLOCK, fp) != BBF_NB_ELEM_BLOCK) return false;
        if (fwrite(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

bool DualBlockedBloomFilter::writeAsBBF(FILE *fp, const uint64_t idx_bbf) const {

    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    if (fwrite(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&nb_bits_per_elem, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&k_, sizeof(int), 1, fp) != 1) return false;

    const uint64_t ush_sz = ush[idx_bbf_norm].size();

    if (fwrite(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (const uint64_t h : ush[idx_bbf_norm]) {

        if (fwrite(&h, sizeof(uint64_t), 1, fp) != 1) return false;
    }

    for (size_t i = 0; i < blocks_; ++i) {

        if (fwrite(&(table_[(i << 1) + idx_bbf_norm].block), sizeof(uint64_t), BBF_NB_ELEM_BLOCK, fp) != BBF_NB_ELEM_BLOCK) return false;
        if (fwrite(&(table_[(i << 1) + idx_bbf_norm].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

bool DualBlockedBloomFilter::read(FILE *fp) {

    uint64_t ush_sz = 0;

    clear();

    if (fread(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&nb_bits_per_elem, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&k_, sizeof(int), 1, fp) != 1) return false;

    for (size_t i = 0; i < 2; ++i) {

        if (fread(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

        for (uint64_t j = 0, h = 0; j != ush_sz; ++j){

            if (fread(&h, sizeof(uint64_t), 1, fp) != 1) return false;

            ush[i].insert(h);
        }
    }

    init_arrays();

    for (uint64_t i = 0; i != (blocks_ << 1); ++i){

        if (fread(&(table_[i].block), sizeof(uint64_t), BBF_NB_ELEM_BLOCK, fp) != BBF_NB_ELEM_BLOCK) return false;
        if (fread(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

bool DualBlockedBloomFilter::readFromBBF(FILE *fp, const uint64_t idx_bbf) {

    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    uint64_t ush_sz = 0;

    clear();

    if (fread(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&nb_bits_per_elem, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&k_, sizeof(int), 1, fp) != 1) return false;

    if (fread(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (uint64_t j = 0, h = 0; j != ush_sz; ++j){

        if (fread(&h, sizeof(uint64_t), 1, fp) != 1) return false;

        ush[idx_bbf_norm].insert(h);
    }

    init_arrays();

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fread(&(table_[(i << 1) + idx_bbf_norm].block), sizeof(uint64_t), BBF_NB_ELEM_BLOCK, fp) != BBF_NB_ELEM_BLOCK) return false;
        if (fread(&(table_[(i << 1) + idx_bbf_norm].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

uint64_t DualBlockedBloomFilter::insert_par(const uint64_t kmh, const uint64_t minh, const uint64_t idx_bbf) {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = /*wyhash(&kmh, sizeof(uint64_t), seed1, _wyp)*/ kmh;
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;
    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int i = 0, j = 0;

    while (i != k_) {

        uint64_t kmh1 = kmh_s1, kmh2 = kmh_s1;
        uint64_t minh2 = minh1 + minh_s2;

        uint64_t bid1 = (fastmod::fastmod_u64(minh1, M_u64, blocks_) << 1) + idx_bbf_norm;
        uint64_t bid2 = (fastmod::fastmod_u64(minh2, M_u64, blocks_) << 1) + idx_bbf_norm;

        if (bid2 < bid1) std::swap(bid1, bid2);

        table_[bid1].lock();

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            if (bid2 == bid1){

                j = i;
                kmh2 = kmh1;
            }
            else {

                table_[bid2].lock();

                for (j = 0; j != k_; ++j, kmh2 += kmh_s2) {

                    if ((table_[bid2].block[(kmh2 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh2 & 0x3fULL))) == 0) break;
                }
            }

            if (j != k_){

                if ((table_[bid1].bits_occupancy < min_overload_bits) || (table_[bid2].bits_occupancy < min_overload_bits)) {

                    uint64_t nb_inserted_bits = 0;

                    if (table_[bid2].bits_occupancy < table_[bid1].bits_occupancy){

                        i = j;
                        kmh1 = kmh2;

                        std::swap(bid1, bid2);
                    }

                    for (; i != k_; ++i, kmh1 += kmh_s2) {

                        const uint64_t div = (kmh1 & BBF_MASK_BITS_BLOCK) >> 6;
                        const uint64_t mod = 1ULL << (kmh1 & 0x3fULL);

                        nb_inserted_bits += static_cast<uint64_t>((table_[bid1].block[div] & mod) == 0);
                        table_[bid1].block[div] |= mod;
                    }

                    table_[bid1].bits_occupancy += nb_inserted_bits;

                    if (bid2 != bid1) table_[bid2].unlock();

                    table_[bid1].unlock();

                    return (((bid1 + 1) << 1) | 0x1ULL);
                }
                else if (nb_overflow == 7) {

                    bool inserted = false;

                    if (bid2 != bid1) table_[bid2].unlock();

                    table_[bid1].unlock();

                    while (lck_ush[idx_bbf_norm].test_and_set(std::memory_order_acquire));

                    inserted = ush[idx_bbf_norm].insert(kmh).second;

                    lck_ush[idx_bbf_norm].clear(std::memory_order_release);

                    return static_cast<uint64_t>(inserted);
                }
            }
            else i = j;

            if (bid2 != bid1) table_[bid2].unlock();
        }

        table_[bid1].unlock();

        minh1 += minh_s2 + minh_s2;
        ++nb_overflow;
    }

    return 0;
}

uint64_t DualBlockedBloomFilter::insert_unpar(const uint64_t kmh, const uint64_t minh, const uint64_t idx_bbf) {

    const uint64_t minh_s1 = /*wyhash(&minh, sizeof(uint64_t), seed1, _wyp)*/ minh;
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = /*wyhash(&kmh, sizeof(uint64_t), seed1, _wyp)*/ kmh;
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = BBF_NB_BITS_BLOCK * BBF_OVERLOAD_RATIO;
    const uint64_t idx_bbf_norm = static_cast<uint64_t>(idx_bbf != 0);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int i = 0, j = 0;

    while (i != k_) {

        uint64_t kmh1 = kmh_s1, kmh2 = kmh_s1;
        uint64_t minh2 = minh1 + minh_s2;

        uint64_t bid1 = (fastmod::fastmod_u64(minh1, M_u64, blocks_) << 1) + idx_bbf_norm;
        uint64_t bid2 = (fastmod::fastmod_u64(minh2, M_u64, blocks_) << 1) + idx_bbf_norm;

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            if (bid2 == bid1){

                j = i;
                kmh2 = kmh1;
            }
            else {

                for (j = 0; j != k_; ++j, kmh2 += kmh_s2) {

                    if ((table_[bid2].block[(kmh2 & BBF_MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh2 & 0x3fULL))) == 0) break;
                }
           }

            if (j != k_){

                if ((table_[bid1].bits_occupancy < min_overload_bits) || (table_[bid2].bits_occupancy < min_overload_bits)) {

                    uint64_t nb_inserted_bits = 0;

                    if (table_[bid2].bits_occupancy < table_[bid1].bits_occupancy){

                        i = j;
                        kmh1 = kmh2;

                        std::swap(bid1, bid2);
                    }

                    for (; i != k_; ++i, kmh1 += kmh_s2) {

                        const uint64_t div = (kmh1 & BBF_MASK_BITS_BLOCK) >> 6;
                        const uint64_t mod = 1ULL << (kmh1 & 0x3fULL);

                        nb_inserted_bits += static_cast<uint64_t>((table_[bid1].block[div] & mod) == 0);
                        table_[bid1].block[div] |= mod;
                    }

                    table_[bid1].bits_occupancy += nb_inserted_bits;

                    return (((bid1 + 1) << 1) | 0x1ULL);
                }
                else if (nb_overflow == 7) {

                    const bool inserted = ush[idx_bbf_norm].insert(kmh).second;

                    return static_cast<uint64_t>(inserted);
                }
            }
            else i = j;
        }

        minh1 += minh_s2 + minh_s2;

        ++nb_overflow;
    }

    return 0;
}

/*CountingBlockedBloomFilter::CountingBlockedBloomFilter() : BlockedBloomFilter::BlockedBloomFilter(), hashbit(nullptr), counts(nullptr) {

    clear();
}

CountingBlockedBloomFilter::CountingBlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) :  BlockedBloomFilter::BlockedBloomFilter(nb_elem, bits_per_elem),
                                                                                                hashbit(nullptr), counts(nullptr), elems_per_block(0) {

    if ((BlockedBloomFilter::blocks_ != 0) && (nb_elem != 0)) {

        elems_per_block = std::max(static_cast<uint64_t>(nb_elem / BlockedBloomFilter::blocks_), static_cast<uint64_t>(1));
        elems_per_block = ((elems_per_block + 15) / 16) * 16; // Round it to largest multiple of 16

        init_arrays();
    }
}

CountingBlockedBloomFilter::CountingBlockedBloomFilter(const CountingBlockedBloomFilter& o) :   BlockedBloomFilter::BlockedBloomFilter(o),
                                                                                                hashbit(nullptr), counts(nullptr), elems_per_block(o.elems_per_block) {

    if ((BlockedBloomFilter::blocks_ != 0) && (elems_per_block != 0)) {

        init_arrays();

        memcpy(hashbit, o.hashbit, BlockedBloomFilter::blocks_ * ((elems_per_block + 15) / 16) * sizeof(uint64_t)); // 4 bits per elem -> 16 elem in a uint64_t
        memcpy(counts, o.counts, BlockedBloomFilter::blocks_ * elems_per_block * sizeof(uint8_t)); 
    }
}

CountingBlockedBloomFilter::CountingBlockedBloomFilter(CountingBlockedBloomFilter&& o) :    BlockedBloomFilter::BlockedBloomFilter(std::move(o)),
                                                                                            hashbit(o.hashbit), counts(o.counts), elems_per_block(o.elems_per_block) {
    o.hashbit = nullptr;
    o.counts = nullptr;

    o.clear();
}

CountingBlockedBloomFilter::~CountingBlockedBloomFilter() {

    clear();
}

CountingBlockedBloomFilter& CountingBlockedBloomFilter::operator=(const CountingBlockedBloomFilter& o) {

    if (this != &o) {

        clear();

        BlockedBloomFilter::operator=(o);

        elems_per_block = o.elems_per_block;

        if ((BlockedBloomFilter::blocks_ != 0) && (elems_per_block != 0)) {

            init_arrays();

            memcpy(hashbit, o.hashbit, BlockedBloomFilter::blocks_ * ((elems_per_block + 15) / 16) * sizeof(uint64_t)); // 4 bits per elem -> 16 elem in a uint64_t
            memcpy(counts, o.counts, BlockedBloomFilter::blocks_ * elems_per_block * sizeof(uint8_t)); 
        }
    }

    return *this;
}


CountingBlockedBloomFilter& CountingBlockedBloomFilter::operator=(CountingBlockedBloomFilter&& o) {

    if (this != &o) {

        clear();

        BlockedBloomFilter::operator=(std::move(o));

        hashbit = o.hashbit;
        counts = o.counts;

        elems_per_block = o.elems_per_block;

        o.hashbit = nullptr;
        o.counts = nullptr;

        o.clear();
    }

    return *this;
}

void CountingBlockedBloomFilter::clear() {

    if (hashbit != nullptr){

        delete[] hashbit;
        hashbit = nullptr;
    }

    if (counts != nullptr){

        delete[] counts;
        counts = nullptr;
    }

    elems_per_block = 0;

    BlockedBloomFilter::clear();
}

void CountingBlockedBloomFilter::init_arrays() {

    hashbit = new uint64_t[blocks_ * ((elems_per_block + 15) / 16)]();
    counts = new uint8_t[blocks_ * elems_per_block]();
}

uint64_t CountingBlockedBloomFilter::insert_par(const uint64_t kmh, const uint64_t minh) {

    const uint64_t bid_ins = BlockedBloomFilter::insert_par(kmh, minh); // Insert in the bloom filter
    const uint64_t bid = bid_ins >> 1; // Block ID of insertion 

    const bool ins = static_cast<bool>(bid_ins & 0x1ULL); // Boolean indicating whether insertion took place

    if (bid != 0) { // bid == 0 means all buckets to insert were overflowing

        const uint64_t pos_start_block = (bid - 1) * elems_per_block;

        const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
        const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

        uint64_t kmh1 = kmh_s1;
        uint64_t count_pos = 0xFFFFFFFFFFFFFFFFULL;

        table_[bid].lock();

        if (ins) { // Element was just inserted so it was not present in the graph before

            for (size_t i = 0; i < 4; ++i, kmh1 += kmh_s2) {

                const uint64_t pos_within_block = kmh1 % elems_per_block;
                const uint64_t pos_abs = pos_start_block + pos_within_block;

                const uint64_t pos_div = pos_abs >> 4; // Div by 16
                const uint64_t pos_mod = ((pos_abs & 0xFULL) << 2) + i;
                const uint64_t mask_ins = 0x1ULL << pos_mod;

                if ((hashbit[pos_div] & mask_ins) == 0){

                    hashbit[pos_div] |= mask_ins;
                    count_pos = pos_abs;

                    break;
                }
            }
        }
        
        if (!ins || (count_pos == 0xFFFFFFFFFFFFFFFFULL)) { 

            uint8_t count = 0xFF;

            for (size_t i = 0; i < 4; ++i, kmh1 += kmh_s2) {

                const uint64_t pos_within_block = kmh1 % elems_per_block;
                const uint64_t pos_abs = pos_start_block + pos_within_block;

                const uint64_t pos_div = pos_abs >> 4; // Div by 16
                const uint64_t pos_mod = ((pos_abs & 0xFULL) << 2) + i;
                const uint64_t mask_ins = 0x1ULL << pos_mod;

                if ((hashbit[pos_div] & mask_ins) == 0) break;
                else if (counts[pos_abs] <= count) {

                    count = counts[pos_abs];
                    count_pos = pos_abs;
                }
            }
        }

        counts[count_pos] += static_cast<uint8_t>(counts[count_pos] != 0xFF);

        table_[bid].unlock();
    }

    return bid_ins;
}

uint64_t CountingBlockedBloomFilter::insert_unpar(const uint64_t kmh, const uint64_t minh) {

    const uint64_t bid_ins = BlockedBloomFilter::insert_unpar(kmh, minh); // Insert in the bloom filter
    const uint64_t bid = bid_ins >> 1; // Block ID of insertion 

    const bool ins = static_cast<bool>(bid_ins & 0x1ULL); // Boolean indicating whether insertion took place

    if (bid != 0) { // bid == 0 means all buckets to insert were overflowing

        const uint64_t pos_start_block = (bid - 1) * elems_per_block;

        const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
        const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

        uint64_t kmh1 = kmh_s1;
        uint64_t count_pos = 0xFFFFFFFFFFFFFFFFULL;

        if (ins) { // Element was just inserted so it was not present in the graph before

            for (size_t i = 0; i < 4; ++i, kmh1 += kmh_s2) {

                const uint64_t pos_within_block = kmh1 % elems_per_block;
                const uint64_t pos_abs = pos_start_block + pos_within_block;

                const uint64_t pos_div = pos_abs >> 4; // Div by 16
                const uint64_t pos_mod = ((pos_abs & 0xFULL) << 2) + i;
                const uint64_t mask_ins = 0x1ULL << pos_mod;

                if ((hashbit[pos_div] & mask_ins) == 0){

                    hashbit[pos_div] |= mask_ins;
                    count_pos = pos_abs;

                    break;
                }
            }
        }
        
        if (!ins || (count_pos == 0xFFFFFFFFFFFFFFFFULL)) { 

            uint8_t count = 0xFF;

            for (size_t i = 0; i < 4; ++i, kmh1 += kmh_s2) {

                const uint64_t pos_within_block = kmh1 % elems_per_block;
                const uint64_t pos_abs = pos_start_block + pos_within_block;

                const uint64_t pos_div = pos_abs >> 4; // Div by 16
                const uint64_t pos_mod = ((pos_abs & 0xFULL) << 2) + i;
                const uint64_t mask_ins = 0x1ULL << pos_mod;

                if ((hashbit[pos_div] & mask_ins) == 0) break;
                else if (counts[pos_abs] <= count) {

                    count = counts[pos_abs];
                    count_pos = pos_abs;
                }
            }
        }

        counts[count_pos] += static_cast<uint8_t>(counts[count_pos] != 0xFF);
    }

    return bid_ins;
}

int64_t CountingBlockedBloomFilter::contains_bids(const uint64_t kmh, const uint64_t minh, const uint64_t min_count) const {

    const int64_t bid = BlockedBloomFilter::contains_bids(kmh, minh);

    if (bid != -1) {

        if (bid == 0) return 0;
        else {

            const uint64_t pos_start_block = (bid - 1) * elems_per_block;

            const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
            const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

            uint64_t kmh1 = kmh_s1;

            uint8_t count = 0;

            for (size_t i = 0; (i < 4) && (count < min_count); ++i, kmh1 += kmh_s2) {

                const uint64_t pos_within_block = kmh1 % elems_per_block;
                const uint64_t pos_abs = pos_start_block + pos_within_block;

                const uint64_t pos_div = pos_abs >> 4; // Div by 16
                const uint64_t pos_mod = ((pos_abs & 0xFULL) << 2) + i;
                const uint64_t mask_ins = 0x1ULL << pos_mod;

                if ((hashbit[pos_div] & mask_ins) != 0) count += counts[pos_abs];
                else break;
            }

            if (count >= min_count) return bid;
        }
    }

    return -1;
}

int CountingBlockedBloomFilter::contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit, const uint64_t min_count) const {

    const std::array<int64_t, 4> bids = BlockedBloomFilter::contains_bids(kmh, minh, 4);

    int cpt = 0;

    pres[0] = false;
    pres[1] = false;
    pres[2] = false;
    pres[3] = false;

    for (size_t i = 0; (i < 4) && (cpt < limit); ++i) {

        if (bids[i] != -1) {

            if (bids[i] == 0) {

                pres[i] = true;
                ++cpt;
            }
            else {

                const uint64_t pos_start_block = (bids[i] - 1) * elems_per_block;

                const uint64_t kmh_s1 = wyhash(&kmh[i], sizeof(uint64_t), seed1, _wyp);
                const uint64_t kmh_s2 = wyhash(&kmh[i], sizeof(uint64_t), seed2, _wyp);

                uint64_t kmh1 = kmh_s1;

                uint8_t count = 0;

                for (size_t ii = 0; (ii < 4) && (count < min_count); ++ii, kmh1 += kmh_s2) {

                    const uint64_t pos_within_block = kmh1 % elems_per_block;
                    const uint64_t pos_abs = pos_start_block + pos_within_block;

                    const uint64_t pos_div = pos_abs >> 4; // Div by 16
                    const uint64_t pos_mod = ((pos_abs & 0xFULL) << 2) + ii;
                    const uint64_t mask_ins = 0x1ULL << pos_mod;

                    if ((hashbit[pos_div] & mask_ins) != 0) count += counts[pos_abs];
                    else break;
                }

                if (count >= min_count) {

                    pres[i] = true;
                    ++cpt;
                }
            }
        }
    }

    return cpt;
}

bool CountingBlockedBloomFilter::write(FILE *fp) const {

    if (!BlockedBloomFilter::write(fp)) return false;

    if (fwrite(&elems_per_block, sizeof(uint64_t), 1, fp) != 1) return false;

    if (fwrite(hashbit, sizeof(uint64_t), blocks_ * ((elems_per_block + 15) / 16), fp) != (blocks_ * ((elems_per_block + 15) / 16))) return false;
    if (fwrite(counts, sizeof(uint8_t), blocks_ * elems_per_block, fp) != (blocks_ * elems_per_block)) return false;

    return true;
}

bool CountingBlockedBloomFilter::read(FILE *fp) {

    clear();

    if (!BlockedBloomFilter::read(fp)) return false;

    if (fread(&elems_per_block, sizeof(uint64_t), 1, fp) != 1) return false;

    init_arrays();

    if (fread(hashbit, sizeof(uint64_t), BlockedBloomFilter::blocks_ * ((elems_per_block + 15) / 16), fp) != (blocks_ * ((elems_per_block + 15) / 16))) return false;
    if (fread(counts, sizeof(uint8_t), BlockedBloomFilter::blocks_ * elems_per_block, fp) != (blocks_ * elems_per_block)) return false;

    return true;
}*/

//#endif
