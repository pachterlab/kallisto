#include "BlockedBloomFilter.hpp"

/*#if defined(__AVX2__)

BlockedBloomFilter::BlockedBloomFilter() : table_(nullptr) {

    clear();
}

BlockedBloomFilter::BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(nullptr) {

    clear();

    if ((nb_elem != 0) && (bits_per_elem != 0)){

        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<unsigned long long> distribution;

        blocks_ = (bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK;
        k_ = (int) (bits_per_elem * log(2));

        if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) ++k_;

        seed1 = distribution(gen);
        seed2 = distribution(gen);

        mask_ksup4 = static_cast<uint16_t>(k_ <= 4) - 1;
        mask_ksup8 = static_cast<uint16_t>(k_ <= 8) - 1;
        mask_ksup12 = static_cast<uint16_t>(k_ <= 12) - 1;

        if (k_ > 16){

            std::cerr << "BlockedBloomFilter(): The AVX2 Blocked Bloom filter does not support more than 16 hash functions." << std::endl;
            std::cerr << "Either use less bits per element to insert or recompile code of Bifrost with AVX2 deactivated." << std::endl;

            clear();
        }
        else {

            for (int k = 0; k != k_; ++k) hashes_mask[k/4] = (hashes_mask[k/4] << 16) | 0xffff;
        }

        init_table();
    }
}

BlockedBloomFilter::BlockedBloomFilter(const BlockedBloomFilter& o) :   table_(nullptr), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_),
                                                                        seed1(o.seed1), seed2(o.seed2), mask_ksup4(o.mask_ksup4), mask_ksup8(o.mask_ksup8),
                                                                        mask_ksup12(o.mask_ksup12), ush(o.ush) {

    std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

    if (blocks_ != 0){

        init_table();

        for (uint64_t i = 0; i != blocks_; ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));

            table_[i].bits_occupancy = o.table_[i].bits_occupancy;
        }
    }
}

BlockedBloomFilter::BlockedBloomFilter(BlockedBloomFilter&& o) :    table_(o.table_), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_),
                                                                    seed1(o.seed1), seed2(o.seed2), mask_ksup4(o.mask_ksup4), mask_ksup8(o.mask_ksup8),
                                                                    mask_ksup12(o.mask_ksup12), ush(move(o.ush)) {

    std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

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
    fast_div_ = o.fast_div_;
    seed1 = o.seed1;
    seed2 = o.seed2;
    mask_ksup4 = o.mask_ksup4;
    mask_ksup8 = o.mask_ksup8;
    mask_ksup12 = o.mask_ksup12;
    ush = o.ush;

    std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));

        table_[i].bits_occupancy = o.table_[i].bits_occupancy;
    }

    return *this;
}

BlockedBloomFilter& BlockedBloomFilter::operator=(BlockedBloomFilter&& o) {

    if (this != &o) {

        clear();

        table_ = o.table_;
        blocks_ = o.blocks_;
        k_ = o.k_;
        fast_div_ = o.fast_div_;
        seed1 = o.seed1;
        seed2 = o.seed2;
        mask_ksup4 = o.mask_ksup4;
        mask_ksup8 = o.mask_ksup8;
        mask_ksup12 = o.mask_ksup12;

        ush = move(o.ush);

        std::memcpy(hashes_mask, o.hashes_mask, 4 * sizeof(uint64_t));

        o.table_ = nullptr;

        o.clear();
    }

    return *this;
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

    mask_ksup4 = 0;
    mask_ksup8 = 0;
    mask_ksup12 = 0;

    memset(hashes_mask, 0, 4 * sizeof(uint64_t));

    ush.clear();

    lck_ush.clear(std::memory_order_release);
}

void BlockedBloomFilter::init_table(){

    fast_div_ = libdivide::divider<uint64_t>(blocks_);
    table_ = new BBF_Block[blocks_];
}

int BlockedBloomFilter::contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit) const {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed1, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed1, _wyp)};

    const uint64_t kmh_s2[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed2, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed2, _wyp)};

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int cpt = 0;

    __m256i table_gather;

    while (nb_overflow < 8) {

        uint64_t bid = minh1 - (minh1 / fast_div_) * blocks_;

        const uint64_t bits_occupancy_1 = table_[bid].bits_occupancy;

        //Gather and compare
        const uint16_t* table = reinterpret_cast<const uint16_t*>(table_[bid].block);

        for (uint8_t j = 0; (j != 4) && (cpt != limit); ++j){

            if (!pres[j]){

                const uint64_t h_div[4] __attribute__((aligned(32))) {  kmh_s1[j],
                                                                        kmh_s1[j] + kmh_s2[j],
                                                                        kmh_s1[j] + kmh_s2[j] + kmh_s2[j],
                                                                        kmh_s1[j] + kmh_s2[j] + kmh_s2[j] + kmh_s2[j]};

                const __m256i km_hashes = _mm256_and_si256(_mm256_set_epi64x(h_div[0], h_div[1], h_div[2], h_div[3]), mask_and_div);
                const __m256i h_shift = _mm256_and_si256(km_hashes, mask_and_mod);

                const __m256i h_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
                const __m256i h_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

                const __m256i h_gather = _mm256_and_si256(mask_h, _mm256_or_si256(h_gather_lsb, h_gather_msb));

                _mm256_stream_si256((__m256i*)h_div, _mm256_srli_epi16(km_hashes, 4));

                const uint16_t* h_div16 = reinterpret_cast<const uint16_t*>(h_div);

                table_gather = _mm256_set_epi16(
                        table[h_div16[15]] & mask_ksup12, table[h_div16[14]] & mask_ksup12, table[h_div16[13]] & mask_ksup12, table[h_div16[12]] & mask_ksup12,
                        table[h_div16[11]] & mask_ksup8, table[h_div16[10]] & mask_ksup8, table[h_div16[9]] & mask_ksup8, table[h_div16[8]] & mask_ksup8,
                        table[h_div16[7]] & mask_ksup4, table[h_div16[6]] & mask_ksup4, table[h_div16[5]] & mask_ksup4, table[h_div16[4]] & mask_ksup4,
                        table[h_div16[3]], table[h_div16[2]], table[h_div16[1]], table[h_div16[0]]
                        );

                const __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, h_gather), h_gather);

                pres[j] = (_mm256_testz_si256(xor_m256i, xor_m256i) != 0);
                cpt += pres[j];
            }
        }

        if (cpt != limit){

            minh1 += minh_s2;
            bid = minh1 - (minh1 / fast_div_) * blocks_;
            table = reinterpret_cast<const uint16_t*>(table_[bid].block);

            const uint64_t bits_occupancy_2 = table_[bid].bits_occupancy;

            for (uint8_t j = 0; (j != 4) && (cpt != limit); ++j){

                if (!pres[j]){

                    const uint64_t h_div[4] __attribute__((aligned(32))) {  kmh_s1[j],
                                                                            kmh_s1[j] + kmh_s2[j],
                                                                            kmh_s1[j] + kmh_s2[j] + kmh_s2[j],
                                                                            kmh_s1[j] + kmh_s2[j] + kmh_s2[j] + kmh_s2[j]};

                    const __m256i km_h = _mm256_and_si256(_mm256_set_epi64x(h_div[0], h_div[1], h_div[2], h_div[3]), mask_and_div);
                    const __m256i h_shift = _mm256_and_si256(km_h, mask_and_mod);

                    const __m256i h_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
                    const __m256i h_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

                    const __m256i h_gather = _mm256_and_si256(mask_h, _mm256_or_si256(h_gather_lsb, h_gather_msb));

                    _mm256_stream_si256((__m256i*)h_div, _mm256_srli_epi16(km_h, 4));

                    const uint16_t* h_div16 = reinterpret_cast<const uint16_t*>(h_div);

                    table_gather = _mm256_set_epi16(
                            table[h_div16[15]] & mask_ksup12, table[h_div16[14]] & mask_ksup12, table[h_div16[13]] & mask_ksup12, table[h_div16[12]] & mask_ksup12,
                            table[h_div16[11]] & mask_ksup8, table[h_div16[10]] & mask_ksup8, table[h_div16[9]] & mask_ksup8, table[h_div16[8]] & mask_ksup8,
                            table[h_div16[7]] & mask_ksup4, table[h_div16[6]] & mask_ksup4, table[h_div16[5]] & mask_ksup4, table[h_div16[4]] & mask_ksup4,
                            table[h_div16[3]], table[h_div16[2]], table[h_div16[1]], table[h_div16[0]]
                            );

                    const __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, h_gather), h_gather);

                    pres[j] = (_mm256_testz_si256(xor_m256i, xor_m256i) != 0);
                    cpt += pres[j];
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

            if (!pres[j] && (ush.find(kmh[j]) != ush.end())) {

                pres[j] = true;
                ++cpt;
            }
        }
    }

    return cpt;
}

bool BlockedBloomFilter::contains(const uint64_t kmh, const uint64_t minh) const {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    const uint64_t h_div[4] __attribute__((aligned(32))) {  kmh_s1,
                                                            kmh_s1 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2 + kmh_s2};

    const __m256i km_h = _mm256_and_si256(_mm256_set_epi64x(h_div[0], h_div[1], h_div[2], h_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_h, mask_and_mod);

    const __m256i h_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i h_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i h_gather = _mm256_and_si256(mask_h, _mm256_or_si256(h_gather_lsb, h_gather_msb));

    _mm256_stream_si256((__m256i*)h_div, _mm256_srli_epi16(km_h, 4));

    const uint16_t* h_div16 = reinterpret_cast<const uint16_t*>(h_div);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    __m256i table_gather;

    while (nb_overflow < 8) {

        uint64_t bid = minh1 - (minh1 / fast_div_) * blocks_;

        const uint64_t bits_occupancy_1 = table_[bid].bits_occupancy;

        //Gather and compare
        const uint16_t* table = reinterpret_cast<const uint16_t*>(table_[bid].block);

        table_gather = _mm256_set_epi16(
                            table[h_div16[15]] & mask_ksup12, table[h_div16[14]] & mask_ksup12, table[h_div16[13]] & mask_ksup12, table[h_div16[12]] & mask_ksup12,
                            table[h_div16[11]] & mask_ksup8, table[h_div16[10]] & mask_ksup8, table[h_div16[9]] & mask_ksup8, table[h_div16[8]] & mask_ksup8,
                            table[h_div16[7]] & mask_ksup4, table[h_div16[6]] & mask_ksup4, table[h_div16[5]] & mask_ksup4, table[h_div16[4]] & mask_ksup4,
                            table[h_div16[3]], table[h_div16[2]], table[h_div16[1]], table[h_div16[0]]
                            );

        __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, h_gather), h_gather);

        if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return true;

        minh1 += minh_s2;
        bid = minh1 - (minh1 / fast_div_) * blocks_;
        table = reinterpret_cast<const uint16_t*>(table_[bid].block);

        const uint64_t bits_occupancy_2 = table_[bid].bits_occupancy;

        const bool overloaded = (bits_occupancy_1 >= min_overload_bits) && (bits_occupancy_2 >= min_overload_bits);

        table_gather = _mm256_set_epi16(
                            table[h_div16[15]] & mask_ksup12, table[h_div16[14]] & mask_ksup12, table[h_div16[13]] & mask_ksup12, table[h_div16[12]] & mask_ksup12,
                            table[h_div16[11]] & mask_ksup8, table[h_div16[10]] & mask_ksup8, table[h_div16[9]] & mask_ksup8, table[h_div16[8]] & mask_ksup8,
                            table[h_div16[7]] & mask_ksup4, table[h_div16[6]] & mask_ksup4, table[h_div16[5]] & mask_ksup4, table[h_div16[4]] & mask_ksup4,
                            table[h_div16[3]], table[h_div16[2]], table[h_div16[1]], table[h_div16[0]]
                            );

        xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, h_gather), h_gather);

        if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return true;
        else if (!overloaded) return false;

        minh1 += minh_s2;
        ++nb_overflow;
    }

    if (nb_overflow >= 8) return (ush.find(kmh) != ush.end());

    return false;
}

int BlockedBloomFilter::contains_bids(const uint64_t kmh, const uint64_t minh) const {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    const uint64_t h_div[4] __attribute__((aligned(32))) {  kmh_s1,
                                                            kmh_s1 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2 + kmh_s2};

    const __m256i km_h = _mm256_and_si256(_mm256_set_epi64x(h_div[0], h_div[1], h_div[2], h_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_h, mask_and_mod);

    const __m256i h_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i h_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i h_gather = _mm256_and_si256(mask_h, _mm256_or_si256(h_gather_lsb, h_gather_msb));

    _mm256_stream_si256((__m256i*)h_div, _mm256_srli_epi16(km_h, 4));

    const uint16_t* h_div16 = reinterpret_cast<const uint16_t*>(h_div);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;
    uint64_t bid;

    while (nb_overflow < 8) {

        bid = minh1 - (minh1 / fast_div_) * blocks_;

        const uint64_t bits_occupancy_1 = table_[bid].bits_occupancy;

        const uint16_t* table = reinterpret_cast<const uint16_t*>(table_[bid].block);

        __m256i table_gather;

        table_gather = _mm256_set_epi16(
                            table[h_div16[15]] & mask_ksup12, table[h_div16[14]] & mask_ksup12, table[h_div16[13]] & mask_ksup12, table[h_div16[12]] & mask_ksup12,
                            table[h_div16[11]] & mask_ksup8, table[h_div16[10]] & mask_ksup8, table[h_div16[9]] & mask_ksup8, table[h_div16[8]] & mask_ksup8,
                            table[h_div16[7]] & mask_ksup4, table[h_div16[6]] & mask_ksup4, table[h_div16[5]] & mask_ksup4, table[h_div16[4]] & mask_ksup4,
                            table[h_div16[3]], table[h_div16[2]], table[h_div16[1]], table[h_div16[0]]
                            );

        __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, h_gather), h_gather);

        if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return bid;

        minh1 += minh_s2;
        bid = minh1 - (minh1 / fast_div_) * blocks_;
        table = reinterpret_cast<const uint16_t*>(table_[bid].block);

        const uint64_t bits_occupancy_2 = table_[bid].bits_occupancy;

        const bool overloaded = (bits_occupancy_1 >= min_overload_bits) && (bits_occupancy_2 >= min_overload_bits);

        table_gather = _mm256_set_epi16(
                            table[h_div16[15]] & mask_ksup12, table[h_div16[14]] & mask_ksup12, table[h_div16[13]] & mask_ksup12, table[h_div16[12]] & mask_ksup12,
                            table[h_div16[11]] & mask_ksup8, table[h_div16[10]] & mask_ksup8, table[h_div16[9]] & mask_ksup8, table[h_div16[8]] & mask_ksup8,
                            table[h_div16[7]] & mask_ksup4, table[h_div16[6]] & mask_ksup4, table[h_div16[5]] & mask_ksup4, table[h_div16[4]] & mask_ksup4,
                            table[h_div16[3]], table[h_div16[2]], table[h_div16[1]], table[h_div16[0]]
                            );

        xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather, h_gather), h_gather);

        if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return bid;
        if (!overloaded) return -1;

        minh1 += minh_s2;
        ++nb_overflow;
    }

    if ((nb_overflow >= 8) && (ush.find(kmh) != ush.end())) return bid;

    return -1;
}

bool BlockedBloomFilter::ReadBloomFilter(FILE *fp) {

    uint64_t ush_sz = 0;

    clear();

    if (fread(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
    if (fread(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&k_, sizeof(k_), 1, fp) != 1) return false;

    if (fread(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (uint64_t i = 0, h = 0; i != ush_sz; ++i){

        if (fread(&h, sizeof(uint64_t), 1, fp) != 1) return false;

        ush.insert(h);
    }

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fread(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
        if (fread(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    for (int k = 0; k != k_; ++k) hashes_mask[k/4] = (hashes_mask[k/4] << 16) | 0xffff;

    mask_ksup4 = static_cast<uint16_t>(k_ <= 4) - 1;
    mask_ksup8 = static_cast<uint16_t>(k_ <= 8) - 1;
    mask_ksup12 = static_cast<uint16_t>(k_ <= 12) - 1;

    return true;
}

bool BlockedBloomFilter::WriteBloomFilter(FILE *fp) const {

    if (fwrite(&blocks_, sizeof(blocks_), 1, fp) != 1) return false;
    if (fwrite(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&k_, sizeof(k_), 1, fp) != 1) return false;

    const uint64_t ush_sz = ush.size();

    if (fwrite(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (const uint64_t h : ush) {

        if (fwrite(&h, sizeof(uint64_t), 1, fp) != 1) return false;
    }

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fwrite(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
        if (fwrite(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

bool BlockedBloomFilter::insert_par(const uint64_t kmh, const uint64_t minh) {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    const uint64_t h_div[4] __attribute__((aligned(32))) {  kmh_s1,
                                                            kmh_s1 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2 + kmh_s2};

    const __m256i km_h = _mm256_and_si256(_mm256_set_epi64x(h_div[0], h_div[1], h_div[2], h_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_h, mask_and_mod);

    const __m256i h_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i h_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i h_gather = _mm256_and_si256(mask_h, _mm256_or_si256(h_gather_lsb, h_gather_msb));

    _mm256_stream_si256((__m256i*)h_div, _mm256_srli_epi16(km_h, 4));

    const uint16_t* h_div16 = reinterpret_cast<const uint16_t*>(h_div);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    while (nb_overflow < 8) {

        const uint64_t minh2 = minh1 + minh_s2;

        uint64_t bid1 = minh1 - (minh1 / fast_div_) * blocks_;
        uint64_t bid2 = minh2 - (minh2 / fast_div_) * blocks_;

        if (bid2 < bid1) std::swap(bid1, bid2);

        uint16_t* table1 = reinterpret_cast<uint16_t*>(table_[bid1].block);
        uint16_t* table2 = reinterpret_cast<uint16_t*>(table_[bid2].block);

         __m256i table_gather1, table_gather2;

        table_[bid1].lock();

        table_gather1 = _mm256_set_epi16(
                            table1[h_div16[15]] & mask_ksup12, table1[h_div16[14]] & mask_ksup12, table1[h_div16[13]] & mask_ksup12, table1[h_div16[12]] & mask_ksup12,
                            table1[h_div16[11]] & mask_ksup8, table1[h_div16[10]] & mask_ksup8, table1[h_div16[9]] & mask_ksup8, table1[h_div16[8]] & mask_ksup8,
                            table1[h_div16[7]] & mask_ksup4, table1[h_div16[6]] & mask_ksup4, table1[h_div16[5]] & mask_ksup4, table1[h_div16[4]] & mask_ksup4,
                            table1[h_div16[3]], table1[h_div16[2]], table1[h_div16[1]], table1[h_div16[0]]
                            );

        __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather1, h_gather), h_gather);

        if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0){

            table_[bid1].unlock();

            return false;
        }

        if (bid1 != bid2){

            table_[bid2].lock();

            table_gather2 = _mm256_set_epi16(
                                table2[h_div16[15]] & mask_ksup12, table2[h_div16[14]] & mask_ksup12, table2[h_div16[13]] & mask_ksup12, table2[h_div16[12]] & mask_ksup12,
                                table2[h_div16[11]] & mask_ksup8, table2[h_div16[10]] & mask_ksup8, table2[h_div16[9]] & mask_ksup8, table2[h_div16[8]] & mask_ksup8,
                                table2[h_div16[7]] & mask_ksup4, table2[h_div16[6]] & mask_ksup4, table2[h_div16[5]] & mask_ksup4, table2[h_div16[4]] & mask_ksup4,
                                table2[h_div16[3]], table2[h_div16[2]], table2[h_div16[1]], table2[h_div16[0]]
                                );

            xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather2, h_gather), h_gather);

            if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0){

                table_[bid1].unlock();
                table_[bid2].unlock();

                return false;
            }
        }

        if ((table_[bid1].bits_occupancy < min_overload_bits) || (table_[bid2].bits_occupancy < min_overload_bits)) {

            uint16_t h_mod[16] __attribute__((aligned(32)));

            _mm256_stream_si256(reinterpret_cast<__m256i*>(h_mod), h_gather);

            if (table_[bid2].bits_occupancy < table_[bid1].bits_occupancy) {

                table1 = table2;
                table_gather1 = table_gather2;

                std::swap(bid1, bid2);
            }

            const uint64_t popcnt_before = popcnt_avx2(&table_gather1, 1);

            switch(k_){
                case 16: table1[h_div16[15]] |= h_mod[15];
                case 15: table1[h_div16[14]] |= h_mod[14];
                case 14: table1[h_div16[13]] |= h_mod[13];
                case 13: table1[h_div16[12]] |= h_mod[12];
                case 12: table1[h_div16[11]] |= h_mod[11];
                case 11: table1[h_div16[10]] |= h_mod[10];
                case 10: table1[h_div16[9]] |= h_mod[9];
                case 9: table1[h_div16[8]] |= h_mod[8];
                case 8: table1[h_div16[7]] |= h_mod[7];
                case 7: table1[h_div16[6]] |= h_mod[6];
                case 6: table1[h_div16[5]] |= h_mod[5];
                case 5: table1[h_div16[4]] |= h_mod[4];
                case 4: table1[h_div16[3]] |= h_mod[3];
                case 3: table1[h_div16[2]] |= h_mod[2];
                case 2: table1[h_div16[1]] |= h_mod[1];
                case 1: table1[h_div16[0]] |= h_mod[0];
            }

            table_gather1 = _mm256_set_epi16(
                            table1[h_div16[15]] & mask_ksup12, table1[h_div16[14]] & mask_ksup12, table1[h_div16[13]] & mask_ksup12, table1[h_div16[12]] & mask_ksup12,
                            table1[h_div16[11]] & mask_ksup8, table1[h_div16[10]] & mask_ksup8, table1[h_div16[9]] & mask_ksup8, table1[h_div16[8]] & mask_ksup8,
                            table1[h_div16[7]] & mask_ksup4, table1[h_div16[6]] & mask_ksup4, table1[h_div16[5]] & mask_ksup4, table1[h_div16[4]] & mask_ksup4,
                            table1[h_div16[3]], table1[h_div16[2]], table1[h_div16[1]], table1[h_div16[0]]
                            );

            table_[bid1].bits_occupancy += popcnt_avx2(&table_gather1, 1) - popcnt_before;

            table_[bid1].unlock();
            
            if (bid1 != bid2) table_[bid2].unlock();

            return true;
        }
        else if (nb_overflow == 7) {

            table_[bid1].unlock();
            
            if (bid1 != bid2) table_[bid2].unlock();

            bool inserted = false;

            while (lck_ush.test_and_set(std::memory_order_acquire));

            inserted = ush.insert(kmh).second;

            lck_ush.clear(std::memory_order_release);

            return inserted;
        }

        table_[bid1].unlock();
            
        if (bid1 != bid2) table_[bid2].unlock();

        minh1 += minh_s2 + minh_s2;
        ++nb_overflow;
    }

    return false;
}

bool BlockedBloomFilter::insert_unpar(const uint64_t kmh, const uint64_t minh) {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    const uint64_t h_div[4] __attribute__((aligned(32))) {  kmh_s1,
                                                            kmh_s1 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2,
                                                            kmh_s1 + kmh_s2 + kmh_s2 + kmh_s2};

    const __m256i km_h = _mm256_and_si256(_mm256_set_epi64x(h_div[0], h_div[1], h_div[2], h_div[3]), mask_and_div);
    const __m256i h_shift = _mm256_and_si256(km_h, mask_and_mod);

    const __m256i h_gather_lsb = _mm256_sllv_epi32(one2shift_lsb, _mm256_and_si256(h_shift, mask_lsb));
    const __m256i h_gather_msb = _mm256_sllv_epi32(one2shift_msb, _mm256_srli_epi32(h_shift, 16));

    const __m256i mask_h = _mm256_set_epi64x(hashes_mask[3], hashes_mask[2], hashes_mask[1], hashes_mask[0]);

    const __m256i h_gather = _mm256_and_si256(mask_h, _mm256_or_si256(h_gather_lsb, h_gather_msb));

    _mm256_stream_si256((__m256i*)h_div, _mm256_srli_epi16(km_h, 4));

    const uint16_t* h_div16 = reinterpret_cast<const uint16_t*>(h_div);

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    while (nb_overflow < 8) {

        const uint64_t minh2 = minh1 + minh_s2;

        uint64_t bid1 = minh1 - (minh1 / fast_div_) * blocks_;
        uint64_t bid2 = minh2 - (minh2 / fast_div_) * blocks_;

        if (bid2 < bid1) std::swap(bid1, bid2);

        uint16_t* table1 = reinterpret_cast<uint16_t*>(table_[bid1].block);
        uint16_t* table2 = reinterpret_cast<uint16_t*>(table_[bid2].block);

        __m256i table_gather1, table_gather2;

        table_gather1 = _mm256_set_epi16(
                            table1[h_div16[15]] & mask_ksup12, table1[h_div16[14]] & mask_ksup12, table1[h_div16[13]] & mask_ksup12, table1[h_div16[12]] & mask_ksup12,
                            table1[h_div16[11]] & mask_ksup8, table1[h_div16[10]] & mask_ksup8, table1[h_div16[9]] & mask_ksup8, table1[h_div16[8]] & mask_ksup8,
                            table1[h_div16[7]] & mask_ksup4, table1[h_div16[6]] & mask_ksup4, table1[h_div16[5]] & mask_ksup4, table1[h_div16[4]] & mask_ksup4,
                            table1[h_div16[3]], table1[h_div16[2]], table1[h_div16[1]], table1[h_div16[0]]
                            );

        __m256i xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather1, h_gather), h_gather);

        if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return false;

        if (bid1 != bid2){

            table_gather2 = _mm256_set_epi16(
                                table2[h_div16[15]] & mask_ksup12, table2[h_div16[14]] & mask_ksup12, table2[h_div16[13]] & mask_ksup12, table2[h_div16[12]] & mask_ksup12,
                                table2[h_div16[11]] & mask_ksup8, table2[h_div16[10]] & mask_ksup8, table2[h_div16[9]] & mask_ksup8, table2[h_div16[8]] & mask_ksup8,
                                table2[h_div16[7]] & mask_ksup4, table2[h_div16[6]] & mask_ksup4, table2[h_div16[5]] & mask_ksup4, table2[h_div16[4]] & mask_ksup4,
                                table2[h_div16[3]], table2[h_div16[2]], table2[h_div16[1]], table2[h_div16[0]]
                                );

            xor_m256i = _mm256_xor_si256(_mm256_and_si256(table_gather2, h_gather), h_gather);

            if (_mm256_testz_si256(xor_m256i, xor_m256i) != 0) return false;
        }

        if ((table_[bid1].bits_occupancy < min_overload_bits) || (table_[bid2].bits_occupancy < min_overload_bits)) {

            uint16_t h_mod[16] __attribute__((aligned(32)));

            _mm256_stream_si256(reinterpret_cast<__m256i*>(h_mod), h_gather);

            if (table_[bid2].bits_occupancy < table_[bid1].bits_occupancy) {

                table1 = table2;
                table_gather1 = table_gather2;

                std::swap(bid1, bid2);
            }

            const uint64_t popcnt_before = popcnt_avx2(&table_gather1, 1);

            switch(k_){
                case 16: table1[h_div16[15]] |= h_mod[15];
                case 15: table1[h_div16[14]] |= h_mod[14];
                case 14: table1[h_div16[13]] |= h_mod[13];
                case 13: table1[h_div16[12]] |= h_mod[12];
                case 12: table1[h_div16[11]] |= h_mod[11];
                case 11: table1[h_div16[10]] |= h_mod[10];
                case 10: table1[h_div16[9]] |= h_mod[9];
                case 9: table1[h_div16[8]] |= h_mod[8];
                case 8: table1[h_div16[7]] |= h_mod[7];
                case 7: table1[h_div16[6]] |= h_mod[6];
                case 6: table1[h_div16[5]] |= h_mod[5];
                case 5: table1[h_div16[4]] |= h_mod[4];
                case 4: table1[h_div16[3]] |= h_mod[3];
                case 3: table1[h_div16[2]] |= h_mod[2];
                case 2: table1[h_div16[1]] |= h_mod[1];
                case 1: table1[h_div16[0]] |= h_mod[0];
            }

            table_gather1 = _mm256_set_epi16(
                            table1[h_div16[15]] & mask_ksup12, table1[h_div16[14]] & mask_ksup12, table1[h_div16[13]] & mask_ksup12, table1[h_div16[12]] & mask_ksup12,
                            table1[h_div16[11]] & mask_ksup8, table1[h_div16[10]] & mask_ksup8, table1[h_div16[9]] & mask_ksup8, table1[h_div16[8]] & mask_ksup8,
                            table1[h_div16[7]] & mask_ksup4, table1[h_div16[6]] & mask_ksup4, table1[h_div16[5]] & mask_ksup4, table1[h_div16[4]] & mask_ksup4,
                            table1[h_div16[3]], table1[h_div16[2]], table1[h_div16[1]], table1[h_div16[0]]
                            );

            table_[bid1].bits_occupancy += popcnt_avx2(&table_gather1, 1) - popcnt_before;

            return true;
        }
        else if (nb_overflow == 7) return ush.insert(kmh).second;

        minh1 += minh_s2 + minh_s2;
        ++nb_overflow;
    }

    return false;
}

const __m256i BlockedBloomFilter::mask_and_div = _mm256_set1_epi16(MASK_BITS_BLOCK);
// All 4 LSB of each 16 bits word
const __m256i BlockedBloomFilter::mask_and_mod = _mm256_set1_epi16(0xf);
// All 1 LSB of each 32 bits word
const __m256i BlockedBloomFilter::one2shift_lsb = _mm256_set1_epi32(0x1);
// Set the 17th LSB bit of each 32 bits word
const __m256i BlockedBloomFilter::one2shift_msb = _mm256_set1_epi32(0x10000);
// All 16 LSB bits of each 32 bits word
const __m256i BlockedBloomFilter::mask_lsb = _mm256_set1_epi64x(0x0000ffff0000ffff);

#else*/

BlockedBloomFilter::BlockedBloomFilter() : table_(nullptr) {

    clear();
}

BlockedBloomFilter::BlockedBloomFilter(size_t nb_elem, size_t bits_per_elem) : table_(nullptr){

    clear();

    if ((nb_elem != 0) && (bits_per_elem != 0)){

        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<unsigned long long> distribution;

        blocks_ = (bits_per_elem * nb_elem + MASK_BITS_BLOCK) / NB_BITS_BLOCK;
        k_ = (int) (bits_per_elem * log(2));

        if (fpp(bits_per_elem, k_) >= fpp(bits_per_elem, k_+1)) ++k_;

        seed1 = distribution(gen);
        seed2 = distribution(gen);

        init_table();
    }
}

BlockedBloomFilter::BlockedBloomFilter(const BlockedBloomFilter& o) : table_(nullptr), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_), seed1(o.seed1), seed2(o.seed2), ush(o.ush) {

    if (blocks_ != 0){

        init_table();

        for (uint64_t i = 0; i != blocks_; ++i){

            memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));

            table_[i].bits_occupancy = o.table_[i].bits_occupancy;
        }
    }
}

BlockedBloomFilter::BlockedBloomFilter(BlockedBloomFilter&& o) : table_(o.table_), blocks_(o.blocks_), k_(o.k_), fast_div_(o.fast_div_), seed1(o.seed1), seed2(o.seed2), ush(move(o.ush)) {

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
    fast_div_ = o.fast_div_;
    seed1 = o.seed1;
    seed2 = o.seed2;
    ush = o.ush;

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        memcpy(&(table_[i].block), &(o.table_[i].block), NB_ELEM_BLOCK * sizeof(uint64_t));

        table_[i].bits_occupancy = o.table_[i].bits_occupancy;
    }

    return *this;
}

BlockedBloomFilter& BlockedBloomFilter::operator=(BlockedBloomFilter&& o) {

    if (this != &o) {

        clear();

        table_ = o.table_;
        blocks_ = o.blocks_;
        k_ = o.k_;
        fast_div_ = o.fast_div_;
        seed1 = o.seed1;
        seed2 = o.seed2;

        ush = move(o.ush);

        o.table_ = nullptr;

        o.clear();
    }

    return *this;
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

    ush.clear();

    lck_ush.clear(std::memory_order_release);
}

void BlockedBloomFilter::init_table(){

    fast_div_ = libdivide::divider<uint64_t>(blocks_);
    table_ = new BBF_Block[blocks_];
}

int BlockedBloomFilter::contains(const uint64_t (&kmh)[4], const uint64_t minh, bool (&pres)[4], const int limit) const {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed1, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed1, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed1, _wyp)};

    const uint64_t kmh_s2[4] = {wyhash(&kmh[0], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[1], sizeof(uint64_t), seed2, _wyp),
                                wyhash(&kmh[2], sizeof(uint64_t), seed2, _wyp), wyhash(&kmh[3], sizeof(uint64_t), seed2, _wyp)};

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int cpt = 0;

    while (nb_overflow < 8) {

        uint64_t bid1 = minh1 - (minh1 / fast_div_) * blocks_;

        const uint64_t bits_occupancy_1 = table_[bid1].bits_occupancy;

        for (uint64_t i = 0, j = 0; (j != 4) && (cpt != limit); ++j){

            if (!pres[j]){

                uint64_t kmh1 = kmh_s1[j];

                for (i = 0; i != k_; ++i, kmh1 += kmh_s2[j]) {

                    if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
                }

                pres[j] = (i == k_);
                cpt += (i == k_);
            }
        }

        if (cpt != limit){

            minh1 += minh_s2;
            bid1 = minh1 - (minh1 / fast_div_) * blocks_;

            const uint64_t bits_occupancy_2 = table_[bid1].bits_occupancy;

            for (size_t i = 0, j = 0; (j != 4) && (cpt != limit); ++j){

                if (!pres[j]){

                    uint64_t kmh1 = kmh_s1[j];

                    for (i = 0; i != k_; ++i, kmh1 += kmh_s2[j]) {

                        if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
                    }

                    pres[j] = (i == k_);
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

            if (!pres[j] && (ush.find(kmh[j]) != ush.end())) {

                pres[j] = true;
                ++cpt;
            }
        }
    }

    return cpt;
}

bool BlockedBloomFilter::contains(const uint64_t kmh, const uint64_t minh) const {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int i = 0;

    while (nb_overflow < 8) {

        uint64_t kmh1 = kmh_s1;
        uint64_t bid1 = minh1 - (minh1 / fast_div_) * blocks_;

        const uint64_t bits_occupancy_1 = table_[bid1].bits_occupancy;

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            minh1 += minh_s2;
            bid1 = minh1 - (minh1 / fast_div_) * blocks_;
            kmh1 = kmh_s1;

            const uint64_t bits_occupancy_2 = table_[bid1].bits_occupancy;

            const bool overloaded = (bits_occupancy_1 >= min_overload_bits) && (bits_occupancy_2 >= min_overload_bits);

            for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

                if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) {

                    if (overloaded) break;
                    
                    return false;
                }
            }

            if (i == k_) return true;

            minh1 += minh_s2;
            ++nb_overflow;
        }
        else return true;
    }

    if (nb_overflow >= 8) return (ush.find(kmh) != ush.end());

    return false;
}

int BlockedBloomFilter::contains_bids(const uint64_t kmh, const uint64_t minh) const {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;
    uint64_t kmh1, bid1;

    int i = 0;

    while (nb_overflow < 8) {

        kmh1 = kmh_s1;
        bid1 = minh1 - (minh1 / fast_div_) * blocks_;

        const uint64_t bits_occupancy_1 = table_[bid1].bits_occupancy;

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            minh1 += minh_s2;
            bid1 = minh1 - (minh1 / fast_div_) * blocks_;
            kmh1 = kmh_s1;

            const uint64_t bits_occupancy_2 = table_[bid1].bits_occupancy;

            const bool overloaded = (bits_occupancy_1 >= min_overload_bits) && (bits_occupancy_2 >= min_overload_bits);

            for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

                if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) {

                    if (overloaded) break;
                    
                    return -1;
                }
            }

            if (i == k_) return bid1;

            minh1 += minh_s2;
            ++nb_overflow;
        }
        else return bid1;
    }

    if ((nb_overflow >= 8) && (ush.find(kmh) != ush.end())) return bid1;

    return -1;
}

bool BlockedBloomFilter::WriteBloomFilter(FILE *fp) const {

    if (fwrite(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fwrite(&k_, sizeof(int), 1, fp) != 1) return false;

    const uint64_t ush_sz = ush.size();

    if (fwrite(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (const uint64_t h : ush) {

        if (fwrite(&h, sizeof(uint64_t), 1, fp) != 1) return false;
    }

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fwrite(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
        if (fwrite(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

bool BlockedBloomFilter::ReadBloomFilter(FILE *fp) {

    uint64_t ush_sz = 0;

    clear();

    if (fread(&blocks_, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed1, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&seed2, sizeof(uint64_t), 1, fp) != 1) return false;
    if (fread(&k_, sizeof(int), 1, fp) != 1) return false;

    if (fread(&ush_sz, sizeof(uint64_t), 1, fp) != 1) return false;

    for (uint64_t i = 0, h = 0; i != ush_sz; ++i){

        if (fread(&h, sizeof(uint64_t), 1, fp) != 1) return false;

        ush.insert(h);
    }

    init_table();

    for (uint64_t i = 0; i != blocks_; ++i){

        if (fread(&(table_[i].block), sizeof(uint64_t), NB_ELEM_BLOCK, fp) != NB_ELEM_BLOCK) return false;
        if (fread(&(table_[i].bits_occupancy), sizeof(uint64_t), 1, fp) != 1) return false;
    }

    return true;
}

bool BlockedBloomFilter::insert_par(const uint64_t kmh, const uint64_t minh) {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int i = 0, j = 0;

    while (i != k_) {

        uint64_t kmh1 = kmh_s1, kmh2 = kmh_s1;
        uint64_t minh2 = minh1 + minh_s2;

        uint64_t bid1 = minh1 - (minh1 / fast_div_) * blocks_;
        uint64_t bid2 = minh2 - (minh2 / fast_div_) * blocks_;

        if (bid2 < bid1) std::swap(bid1, bid2);

        table_[bid1].lock();

        for (i = 0; i != k_; ++i, kmh1 += kmh_s2) {

            if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;
        }

        if (i != k_){

            if (bid2 == bid1){

                j = i;
                kmh2 = kmh1;
            }
            else {

                table_[bid2].lock();

                for (j = 0; j != k_; ++j, kmh2 += kmh_s2) {

                    if ((table_[bid2].block[(kmh2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh2 & 0x3fULL))) == 0) break;
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

                        const uint64_t div = (kmh1 & MASK_BITS_BLOCK) >> 6;
                        const uint64_t mod = 1ULL << (kmh1 & 0x3fULL);

                        nb_inserted_bits += static_cast<uint64_t>((table_[bid1].block[div] & mod) == 0);
                        table_[bid1].block[div] |= mod;
                    }

                    table_[bid1].bits_occupancy += nb_inserted_bits;

                    if (bid2 != bid1) table_[bid2].unlock();

                    table_[bid1].unlock();

                    return true;
                }
                else if (nb_overflow == 7) {

                    if (bid2 != bid1) table_[bid2].unlock();

                    table_[bid1].unlock();

                    bool inserted = false;

                    while (lck_ush.test_and_set(std::memory_order_acquire));

                    inserted = ush.insert(kmh).second;

                    lck_ush.clear(std::memory_order_release);

                    return inserted;
                }
            }
            else i = j;

            if (bid2 != bid1) table_[bid2].unlock();
        }

        table_[bid1].unlock();

        minh1 += minh_s2 + minh_s2;
        ++nb_overflow;
    }

    return false;
}

bool BlockedBloomFilter::insert_unpar(const uint64_t kmh, const uint64_t minh) {

    const uint64_t minh_s1 = wyhash(&minh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t minh_s2 = wyhash(&minh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t kmh_s1 = wyhash(&kmh, sizeof(uint64_t), seed1, _wyp);
    const uint64_t kmh_s2 = wyhash(&kmh, sizeof(uint64_t), seed2, _wyp);

    const uint64_t min_overload_bits = NB_BITS_BLOCK * 0.65;

    uint64_t nb_overflow = 0;
    uint64_t minh1 = minh_s1;

    int i = 0, j = 0;

    while (i != k_) {

        uint64_t kmh1 = kmh_s1;
        uint64_t bid1 = (minh1 - (minh1 / fast_div_) * blocks_);

        for (i = 0; i != k_; ++i) {

            if ((table_[bid1].block[(kmh1 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh1 & 0x3fULL))) == 0) break;

            kmh1 += kmh_s2;
        }

        if (i != k_){

            uint64_t kmh2 = kmh_s1;
            uint64_t minh2 = minh1 + minh_s2;
            uint64_t bid2 = (minh2 - (minh2 / fast_div_) * blocks_);

            for (j = 0; j != k_; ++j) {

                if ((table_[bid2].block[(kmh2 & MASK_BITS_BLOCK) >> 6] & (1ULL << (kmh2 & 0x3fULL))) == 0) break;

                kmh2 += kmh_s2;
            }

            if (j != k_){

                if ((table_[bid1].bits_occupancy < min_overload_bits) || (table_[bid2].bits_occupancy < min_overload_bits)) {

                    uint64_t nb_inserted_bits = 0;

                    if (table_[bid2].bits_occupancy < table_[bid1].bits_occupancy){

                        i = j;
                        bid1 = bid2;
                        kmh1 = kmh2;
                    }

                    for (; i != k_; ++i) {

                        const uint64_t div = (kmh1 & MASK_BITS_BLOCK) >> 6;
                        const uint64_t mod = 1ULL << (kmh1 & 0x3fULL);

                        nb_inserted_bits += static_cast<uint64_t>((table_[bid1].block[div] & mod) == 0);
                        table_[bid1].block[div] |= mod;

                        kmh1 += kmh_s2;
                    }

                    table_[bid1].bits_occupancy += nb_inserted_bits;

                    return true;
                }
                else if (nb_overflow == 7) return ush.insert(kmh).second;
            }
            else i = j;
        }

        minh1 = minh_s2 + minh_s2;

        ++nb_overflow;
    }

    return false;
}

//#endif
