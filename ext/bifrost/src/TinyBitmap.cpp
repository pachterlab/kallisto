#include "TinyBitmap.hpp"

TinyBitmap::TinyBitmap() : tiny_bmp(nullptr) {}

TinyBitmap::TinyBitmap(const TinyBitmap& o) : tiny_bmp(nullptr) {

    if (o.tiny_bmp != nullptr){

        const uint16_t sz = o.getSize();

        const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp), 8, sz * sizeof(uint16_t));

        if (aligned_alloc != 0){

            cerr << "TinyBitmap::TinyBitmap(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
            exit(1);
        }

        std::copy(o.tiny_bmp, o.tiny_bmp + sz, tiny_bmp);
    }
}

TinyBitmap::TinyBitmap(TinyBitmap&& o) {

    tiny_bmp = o.tiny_bmp;
    o.tiny_bmp = nullptr;
}

TinyBitmap::TinyBitmap(uint16_t** o_ptr) : tiny_bmp(*o_ptr) {

    *o_ptr = nullptr;
}

TinyBitmap& TinyBitmap::operator=(const TinyBitmap& o) {

    if (this != &o){

        clear();

        if (o.tiny_bmp != nullptr){

            const uint16_t sz = o.getSize();

            const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp), 8, sz * sizeof(uint16_t));

            if (aligned_alloc != 0){

                cerr << "TinyBitmap::operator=(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
                exit(1);
            }

            std::copy(o.tiny_bmp, o.tiny_bmp + sz, tiny_bmp);
        }
    }

    return *this;
}

TinyBitmap& TinyBitmap::operator=(TinyBitmap&& o) {

    if (this != &o){

        clear();

        tiny_bmp = o.tiny_bmp;
        o.tiny_bmp = nullptr;
    }

    return *this;
}

TinyBitmap& TinyBitmap::operator=(uint16_t** o_ptr) {

    clear();

    tiny_bmp = *o_ptr;
    *o_ptr = nullptr;

    return *this;
}

TinyBitmap::~TinyBitmap() { clear(); }

void TinyBitmap::clear() {

    if (tiny_bmp != nullptr){

        free(tiny_bmp);
        tiny_bmp = nullptr;
    }
}

bool TinyBitmap::add(const uint32_t val){

    const uint16_t val_div = val >> 16;
    const uint16_t val_mod = val & 0xFFFF;

    if (tiny_bmp == nullptr){

        const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp), 8, sizes[0] * sizeof(uint16_t));

        if (aligned_alloc != 0){

            cerr << "TinyBitmap::add(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
            exit(1);
        }

        std::memset(tiny_bmp, 0, sizes[0] * sizeof(uint16_t));

        tiny_bmp[0] = (sizes[0] << 3) | bmp_mode | bits_16;
        tiny_bmp[2] = val_div;
    }

    if (getOffset() != val_div) return false;

    uint16_t sz = getSize();
    uint16_t mode = getMode();
    uint16_t cardinality = getCardinality();

    // Compute if inserting new value triggers an increase of the container size
    if (((mode == bmp_mode) && (val_mod >= ((sz - 3) << 4))) || ((mode != bmp_mode) && (cardinality >= (sz - 3 - (mode == rle_list_mode))))){

        // Means that in its current mode, container size must be increased to add a value
        // We need to compute if which mode has the smaller container size

        if ((mode != bmp_mode) && contains(val)) return true;

        runOptimize();

        sz = getSize();
        mode = getMode();
        cardinality = getCardinality();

        if ((mode != rle_list_mode) || (cardinality > (sz - 5))){

            const uint16_t nb_uint_bmp = getNextSize((std::max(val_mod, static_cast<const uint16_t>(maximum() & 0xFFFF)) >> 4) + 4);

            bool res;

            if (mode == rle_list_mode){

                const uint16_t nb_val = size();
                const uint16_t nb_val_rle_list = getNextSize(getNextSize(sz) + 1);
                const uint16_t nb_val_list = getNextSize(nb_val + 4);
                const uint16_t nb_val_min = (nb_val > (0xFFFF - 48)) ? 0xFFFF : std::min(nb_val_rle_list, std::min(nb_val_list, nb_uint_bmp));

                if (nb_val_min > sizes[nb_sizes - 1]) return false;

                res = (nb_val_rle_list == nb_val_min);
                res = res ? change_sz(nb_val_min) : switch_mode(nb_val_min, (nb_val_list <= nb_uint_bmp) ? list_mode : bmp_mode);

                cardinality = getCardinality();
            }
            else {

                const uint16_t nb_val_list = getNextSize(cardinality + 4);

                if (mode == bmp_mode) res = (nb_uint_bmp <= nb_val_list) ? change_sz(nb_uint_bmp) : switch_mode(nb_val_list, list_mode);
                else res = (nb_val_list <= nb_uint_bmp) ? change_sz(nb_val_list) : switch_mode(nb_uint_bmp, bmp_mode);
            }

            if (!res) return false;

            mode = getMode();
        }
    }

    if (mode == bmp_mode){ // Bitmap mode

        uint16_t& div = tiny_bmp[(val_mod >> 4) + 3];

        const uint16_t mod = 1U << (val_mod & 0xF);

        tiny_bmp[1] += ((div & mod) == 0); // += (1 << 8) if not already set (increase cardinality in header), 0 otherwise
        div |= mod; // Insert new value
    }
    else if (mode == list_mode) { // Binary search

        if (cardinality == 0){

            ++tiny_bmp[1];
            tiny_bmp[3] = val_mod;
        }
        else {

            uint16_t imid;

            uint16_t imin = 3;
            uint16_t imax = cardinality + 2;

            while (imin < imax){

                imid = (imin + imax) >> 1;

                if (tiny_bmp[imid] < val_mod) imin = imid + 1;
                else imax = imid;
            }

            if (tiny_bmp[imin] != val_mod){

                if (tiny_bmp[imin] < val_mod){

                    std::memmove(&tiny_bmp[imin + 2], &tiny_bmp[imin + 1], (cardinality + 2 - imin) * sizeof(uint16_t)); // Shift values
                    tiny_bmp[imin + 1] = val_mod; // Insert value
                }
                else {

                    std::memmove(&tiny_bmp[imin + 1], &tiny_bmp[imin], (cardinality + 3 - imin) * sizeof(uint16_t)); // Shift values
                    tiny_bmp[imin] = val_mod; // Insert value
                }

                ++tiny_bmp[1]; // Increase cardinality
            }
        }
    }
    else { // Binary search

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_mod) imin = imid + 2;
            else imax = imid;
        }

        if ((val_mod < tiny_bmp[imin]) || (val_mod > tiny_bmp[imin + 1])){

            if (val_mod > tiny_bmp[imin + 1]){

                const bool before_end = imin < (cardinality + 1);

                if (val_mod == (tiny_bmp[imin + 1] + 1)){ // The run can be extended after its end

                    if (before_end && (tiny_bmp[imin + 2] == val_mod + 1)){
                        // The inserted value merges the gap between the current and next run, both runs can be merged
                        std::memmove(&tiny_bmp[imin + 1], &tiny_bmp[imin + 3], (cardinality - imin) * sizeof(uint16_t)); // Shift values
                        tiny_bmp[1] -= 2;
                    }
                    else ++tiny_bmp[imin + 1]; // Just extend the end of the run
                }
                else if (before_end && (val_mod == (tiny_bmp[imin + 2] - 1))) --tiny_bmp[imin + 2];
                else {

                    std::memmove(&tiny_bmp[imin + 4], &tiny_bmp[imin + 2], (cardinality + 1 - imin) * sizeof(uint16_t)); // Shift values

                    tiny_bmp[imin + 2] = val_mod; // Insert start run
                    tiny_bmp[imin + 3] = val_mod; // Insert end run
                    tiny_bmp[1] += 2; // Increase cardinality
                }
            }
            // val_mod < tiny_bmp[imin]
            else if (val_mod == (tiny_bmp[imin] - 1)){

                if ((imin > 3) && (tiny_bmp[imin - 1] == val_mod - 1)){
                    // The inserted value merges the gap between the current and previous run, both runs can be merged
                    std::memmove(&tiny_bmp[imin - 1], &tiny_bmp[imin + 1], (cardinality + 2 - imin) * sizeof(uint16_t)); // Shift values
                    tiny_bmp[1] -= 2;
                }
                else --tiny_bmp[imin]; // Just extend the beginning of the run
            }
            else if ((imin > 3) && (val_mod == (tiny_bmp[imin - 1] + 1))) ++tiny_bmp[imin - 1];
            else {

                std::memmove(&tiny_bmp[imin + 2], &tiny_bmp[imin], (cardinality + 3 - imin) * sizeof(uint16_t)); // Shift values

                tiny_bmp[imin] = val_mod; // Insert start run
                tiny_bmp[imin + 1] = val_mod; // Insert end run
                tiny_bmp[1] += 2; // Increase cardinality
            }
        }
    }

    return true;
}

bool TinyBitmap::contains(const uint32_t val) const {

    // If not allocated or cardinality is 0, val is not present
    if ((tiny_bmp == nullptr) || (getCardinality() == 0) || ((val >> 16) != getOffset())) return false;

    const uint16_t mode = getMode();
    const uint16_t cardinality = getCardinality();
    const uint16_t val_mod = val & 0xFFFF;

    // Bitmap mode
    if (mode == bmp_mode){

        if (val_mod >= ((getSize() - 3) << 4)) return false;

        return ((tiny_bmp[(val_mod >> 4) + 3] & (1U << (val_mod & 0xF))) != 0);
    }
    else if (mode == list_mode) {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 2;

        while (imin < imax){

            imid = (imin + imax) >> 1;

            if (tiny_bmp[imid] < val_mod) imin = imid + 1;
            else imax = imid;
        }

        return (tiny_bmp[imin] == val_mod);
    }
    else {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_mod) imin = imid + 2;
            else imax = imid;
        }

        return ((val_mod >= tiny_bmp[imin]) && (val_mod <= tiny_bmp[imin + 1]));
    }

    return false;
}

/*uint32_t TinyBitmap::and_cardinality(const TinyBitmap& rhs) const {

    // If not allocated or cardinality is 0, val is not present
    if ((tiny_bmp == nullptr) || (getCardinality() == 0)) return false;
    if ((rhs.tiny_bmp == nullptr) || (rhs.getCardinality() == 0)) return false;
    if (getOffset() != rhs.getOffset()) return false;

    const uint16_t mode = getMode(), rhs_mode = rhs.getMode();
    const uint16_t cardinality = getCardinality(), rhs_cardinality = rhs.getCardinality();

    if (mode == rhs_mode) { // Both TinyBitmaps are encoded the same way


    }

    // Bitmap mode
    if (mode == bmp_mode){

        if (val_mod >= ((getSize() - 3) << 4)) return false;

        return ((tiny_bmp[(val_mod >> 4) + 3] & (1U << (val_mod & 0xF))) != 0);
    }
    else if (mode == list_mode) {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 2;

        while (imin < imax){

            imid = (imin + imax) >> 1;

            if (tiny_bmp[imid] < val_mod) imin = imid + 1;
            else imax = imid;
        }

        return (tiny_bmp[imin] == val_mod);
    }
    else {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_mod) imin = imid + 2;
            else imax = imid;
        }

        return ((val_mod >= tiny_bmp[imin]) && (val_mod <= tiny_bmp[imin + 1]));
    }

    return false;
}*/

bool TinyBitmap::containsRange(const uint32_t val_start, const uint32_t val_end) const {

    // If not allocated or cardinality is 0, val is not present
    if ((tiny_bmp == nullptr) || (val_end < val_start)) return false;
    if (val_start == val_end) return contains(val_start);

    const uint16_t offset = getOffset();
    const uint16_t cardinality = getCardinality();

    if ((cardinality == 0) || ((val_start >> 16) != offset) || ((val_end >> 16) != offset)) return false;

    const uint16_t mode = getMode();

    uint16_t val_start_mod = val_start & 0xFFFF;
    const uint16_t val_end_mod = val_end & 0xFFFF;

    // Bitmap mode
    if (mode == bmp_mode){

        if (val_end_mod < ((getSize() - 3) << 4)){ // If size of container is big enough to contain the values

            if ((val_start_mod >> 4) == (val_end_mod >> 4)){ // If start value is on same slot as end value

                const uint16_t mask = ((1 << (val_end_mod & 0xF)) - 1) - ((1 << (val_start_mod & 0xF)) - 1);

                return ((tiny_bmp[(val_start_mod >> 4) + 3] & mask) == mask);
            }
            else {

                const uint16_t mask_start = ~((1 << (val_start_mod & 0xF)) - 1);

                if ((tiny_bmp[(val_start_mod >> 4) + 3] & mask_start) != mask_start) return false;

                const uint16_t mask_end = ((1 << (val_end_mod & 0xF)) << 1) - 1;

                if ((tiny_bmp[(val_end_mod >> 4) + 3] & mask_end) != mask_end) return false;

                for (size_t i = (val_start_mod >> 4) + 4; i != (val_end_mod >> 4) + 3; ++i){

                    if (tiny_bmp[i] != 0xFFFF) return false;
                }

                return true;
            }
        }
    }
    else if (mode == list_mode) {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 2;

        while (imin < imax){

            imid = (imin + imax) >> 1;

            if (tiny_bmp[imid] < val_start_mod) imin = imid + 1;
            else imax = imid;
        }

        if ((cardinality + 3 - imin) >= (val_end_mod - val_start_mod + 1)) {

            while ((val_start_mod <= val_end_mod) && (imin < cardinality + 3) && (tiny_bmp[imin] == val_start_mod)){

                ++imin;
                ++val_start_mod;
            }
        }

        return (val_start_mod > val_end_mod);
    }
    else {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_start_mod) imin = imid + 2;
            else imax = imid;
        }

        return ((val_start_mod >= tiny_bmp[imin]) && (val_end_mod <= tiny_bmp[imin + 1]));
    }

    return false;
}

bool TinyBitmap::remove(const uint32_t val){

    // If not allocated or cardinality is 0, val is not present
    if ((tiny_bmp == nullptr) || (getCardinality() == 0) || ((val >> 16) != getOffset())) return true;

    const uint16_t mode = getMode();
    const uint16_t cardinality = getCardinality();
    const uint16_t val_mod = val & 0xFFFF;

    bool try_decrease_sz = false;

    // Bitmap mode
    if (mode == bmp_mode){

        if (val_mod >= ((getSize() - 3) << 4)) return true;

        uint16_t& div = tiny_bmp[(val_mod >> 4) + 3];
        const uint16_t mod = 1U << (val_mod & 0xF);

        if ((div & mod) != 0){

            div &= ~mod;
            --tiny_bmp[1];
            try_decrease_sz = true;
        }
    }
    else if (mode == list_mode) {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 2;

        while (imin < imax){

            imid = (imin + imax) >> 1;

            if (tiny_bmp[imid] < val_mod) imin = imid + 1;
            else imax = imid;
        }

        if (tiny_bmp[imin] == val_mod){

            std::memmove(&tiny_bmp[imin], &tiny_bmp[imin + 1], (cardinality + 2 - imin) * sizeof(uint16_t)); // Shift values
            --tiny_bmp[1];
            try_decrease_sz = true;
        }
    }
    else {

        uint16_t imid;

        uint16_t imin = 3;
        uint16_t imax = cardinality + 1;

        while (imin < imax){

            imid = (imin + imax) >> 1;
            imid -= ((imid & 0x1) == 0);

            if (tiny_bmp[imid + 1] < val_mod) imin = imid + 2;
            else imax = imid;
        }

        if ((val_mod >= tiny_bmp[imin]) && (val_mod <= tiny_bmp[imin + 1])){

            if ((val_mod == tiny_bmp[imin]) && (val_mod == tiny_bmp[imin + 1])){ // The run is the value to delete

                std::memmove(&tiny_bmp[imin], &tiny_bmp[imin + 2], (cardinality + 1 - imin) * sizeof(uint16_t)); // Shift values
                tiny_bmp[1] -= 2;
            }
            else if (val_mod == tiny_bmp[imin]) ++tiny_bmp[imin];
            else if (val_mod == tiny_bmp[imin + 1]) --tiny_bmp[imin + 1];
            else if ((cardinality + 5) <= getSize()){ // There is enough space to insert a new run

                std::memmove(&tiny_bmp[imin + 3], &tiny_bmp[imin + 1], (cardinality + 2 - imin) * sizeof(uint16_t)); // Shift values

                tiny_bmp[imin + 1] = val_mod - 1; // Insert start run
                tiny_bmp[imin + 2] = val_mod + 1; // Insert start run
                tiny_bmp[1] += 2; // Increase cardinality
            }
            else {

                bool ret = switch_mode(size() + 3, bmp_mode);

                if (ret){

                    ret = remove(val);
                    runOptimize();
                }

                return ret;
            }
        }
    }

    if (tiny_bmp[1] == 0) clear();
    else if (try_decrease_sz){ // Can only be bmp_mode || list_mode

        const uint16_t nb_uint_bmp = getNextSize((static_cast<const uint16_t>(maximum() & 0xFFFF) >> 4) + 4);
        const uint16_t nb_val_list = getNextSize(cardinality + 2);

        if (std::min(nb_uint_bmp, nb_val_list) < getSize()){

            if (mode == bmp_mode) (nb_uint_bmp <= nb_val_list) ? change_sz(nb_uint_bmp) : switch_mode(nb_val_list, list_mode);
            else (nb_val_list <= nb_uint_bmp) ? change_sz(nb_val_list) : switch_mode(nb_uint_bmp, bmp_mode);
        }
    }

    return true;
}

uint32_t TinyBitmap::maximum() const {

    uint16_t card;

    if ((tiny_bmp == nullptr) || ((card = getCardinality()) == 0)) return 0; // If not allocated or cardinality is 0, return 0

    const uint32_t offset = static_cast<uint32_t>(getOffset()) << 16;

    if (getMode() == bmp_mode){

        for (uint16_t i = getSize() - 1; i != 2; --i){

            for (uint16_t j = 15, e = tiny_bmp[i]; e != 0; --j, e <<= 1){

                if ((e & 0x8000) != 0) return offset | (((i - 3) << 4) + j);
            }
        }
    }

    return offset | tiny_bmp[card + 2]; // mode = list_mode
}

size_t TinyBitmap::getSizeInBytes() const {

    if (tiny_bmp == nullptr) return sizeof(uint16_t*);
    return sizeof(uint16_t*) + getSize() * sizeof(uint16_t);
}

size_t TinyBitmap::size() const {

    if (tiny_bmp == nullptr) return 0;

    const uint16_t cardinality = getCardinality();

    if (getMode() == rle_list_mode){

        size_t card = 0;

        for (size_t i = 3; i < cardinality + 3; i += 2) card += tiny_bmp[i + 1] - tiny_bmp[i];

        return card + (cardinality >> 1);
    }

    return cardinality;
}

size_t TinyBitmap::size(uint32_t start_value, const uint32_t end_value) const {

    if ((tiny_bmp == nullptr) || (end_value < start_value)) return 0;

    size_t cpt = 0;

    while (start_value != end_value) cpt += contains(start_value++);

    return cpt;
}

size_t TinyBitmap::runOptimize() {

    if (tiny_bmp != nullptr){

        const uint16_t mode = getMode();
        const uint16_t cardinality = getCardinality();

        if ((mode != rle_list_mode) && (cardinality != 0)){

            const uint16_t sz = getSize();

            if (mode == bmp_mode){

                uint16_t nb_run = 0;
                uint16_t prev_val_pres = 0xFFFE; // This value cannot exist in bitmap mode
                uint16_t cardinality_cpy = cardinality;

                for (uint16_t i = 3; (i != sz) && (cardinality_cpy != 0); ++i){

                    for (uint16_t j = (i - 3) << 4, e = tiny_bmp[i]; e != 0; e >>= 1, ++j){

                        if (e & 0x1){

                            nb_run += (j != (prev_val_pres + 1));
                            --cardinality_cpy;
                            prev_val_pres = j;
                        }
                    }
                }

                const uint16_t new_cardinality = nb_run << 1;
                const uint16_t new_sz = getNextSize(new_cardinality + 3);

                if ((new_cardinality < cardinality) && (new_sz <= sizes[nb_sizes - 1])){

                    cardinality_cpy = cardinality;
                    prev_val_pres = 0xFFFE;

                    uint16_t k = 3;

                    uint16_t* tiny_bmp_new = nullptr;

                    const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp_new), 8, new_sz * sizeof(uint16_t));

                    if (aligned_alloc != 0){

                        cerr << "TinyBitmap::runOptimize(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
                        exit(1);
                    }

                    std::memset(tiny_bmp_new, 0, new_sz * sizeof(uint16_t));

                    for (uint16_t i = 3; (i != sz) && (cardinality_cpy != 0); ++i){

                        for (uint16_t j = (i - 3) << 4, e = tiny_bmp[i]; e != 0; e >>= 1, ++j){

                            if (e & 0x1){

                                if (j != (prev_val_pres + 1)){

                                    if (prev_val_pres != 0xFFFE) tiny_bmp_new[k++] = prev_val_pres;
                                    tiny_bmp_new[k++] = j;
                                }

                                --cardinality_cpy;
                                prev_val_pres = j;
                            }
                        }
                    }

                    tiny_bmp_new[k] = prev_val_pres;
                    tiny_bmp_new[0] = (new_sz << 3) | rle_list_mode | bits_16;
                    tiny_bmp_new[1] = new_cardinality;
                    tiny_bmp_new[2] = getOffset();

                    free(tiny_bmp);

                    tiny_bmp = tiny_bmp_new;

                    return sz - new_sz;
                }
            }
            else {

                uint16_t nb_run = 1;

                for (size_t i = 4; i < cardinality + 3; ++i) nb_run += (tiny_bmp[i] != (tiny_bmp[i-1] + 1));

                const uint16_t new_cardinality = nb_run << 1;
                const uint16_t new_sz = getNextSize(new_cardinality + 3);

                if ((new_cardinality < cardinality) && (new_sz <= sizes[nb_sizes - 1])){

                    uint16_t k = 4;

                    uint16_t* tiny_bmp_new = nullptr;

                    const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp_new), 8, new_sz * sizeof(uint16_t));

                    if (aligned_alloc != 0){

                        cerr << "TinyBitmap::runOptimize(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
                        exit(1);
                    }

                    std::memset(tiny_bmp_new, 0, new_sz * sizeof(uint16_t));

                    tiny_bmp_new[0] = (new_sz << 3) | rle_list_mode | bits_16;
                    tiny_bmp_new[1] = new_cardinality;
                    tiny_bmp_new[2] = getOffset();
                    tiny_bmp_new[3] = tiny_bmp[3];

                    for (size_t i = 4; i < cardinality + 3; ++i){

                        if (tiny_bmp[i] != (tiny_bmp[i-1] + 1)){

                            tiny_bmp_new[k++] = tiny_bmp[i-1];
                            tiny_bmp_new[k++] = tiny_bmp[i];
                        }
                    }

                    tiny_bmp_new[k] = tiny_bmp[cardinality + 2];

                    free(tiny_bmp);

                    tiny_bmp = tiny_bmp_new;

                    return sz - new_sz;
                }
            }
        }
    }

    return 0;
}

size_t TinyBitmap::shrinkSize() {

    if (tiny_bmp == nullptr) return 0;

    const uint16_t sz = getSize();
    const uint16_t mode = getMode();
    const uint16_t new_sz = (mode == bmp_mode) ? (static_cast<const uint16_t>(maximum() & 0xFFFF) >> 4) + 4 : getCardinality() + 3;

    uint16_t* new_t_bmp = nullptr;

    const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&new_t_bmp), 8, new_sz * sizeof(uint16_t));

    if (aligned_alloc != 0){

        cerr << "TinyBitmap::shrinkSize(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
        exit(1);
    }

    std::copy(tiny_bmp, tiny_bmp + new_sz, new_t_bmp);

    free(tiny_bmp);
    tiny_bmp = new_t_bmp;

    tiny_bmp[0] = (tiny_bmp[0] & ~sz_mask) | (new_sz << 3);

    return sz - new_sz;
}

bool TinyBitmap::write(ostream& stream_out) const {

    uint16_t header;

    bool ret;

    if (tiny_bmp == nullptr){

        header = bmp_mode | bits_16; // Size is 0
        ret = stream_out.write(reinterpret_cast<const char*>(&header), sizeof(uint16_t)).fail();
    }
    else {

        const uint16_t new_sz = (getMode() == bmp_mode) ? (static_cast<const uint16_t>(maximum() & 0xFFFF) >> 4) + 4 : getCardinality() + 3;

        header = (tiny_bmp[0] & ~sz_mask) | (new_sz << 3);

        if ((ret = stream_out.write(reinterpret_cast<const char*>(&header), sizeof(uint16_t)).fail()) == false){

            ret = stream_out.write(reinterpret_cast<const char*>(&tiny_bmp[1]), (new_sz - 1) * sizeof(uint16_t)).fail();
        }
    }

    return !ret;
}

bool TinyBitmap::read(istream& stream_in) {

    clear();

    uint16_t header;

    bool ret = stream_in.read(reinterpret_cast<char*>(&header), sizeof(uint16_t)).fail();

    if (ret == false){

        const uint16_t sz = header >> 3;

        if (sz != 0){

            const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp), 8, sz * sizeof(uint16_t));

            if (aligned_alloc != 0){

                cerr << "TinyBitmap::read(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
                exit(1);
            }

            ret = stream_in.read(reinterpret_cast<char*>(&tiny_bmp[1]), (sz - 1) * sizeof(uint16_t)).fail();

            tiny_bmp[0] = header;
        }
    }

    return !ret;
}

void TinyBitmap::print() const {

    if (tiny_bmp == nullptr) return;

    const uint16_t sz = getSize();
    const uint16_t mode = getMode();
    const uint16_t cardinality = getCardinality();

    cout << "sz = " << sz << endl;
    cout << "mode = " << mode << endl;
    cout << "cardinality = " << cardinality << endl;

    if (mode == bmp_mode){

        const uint32_t max_value = maximum() & 0xFFFF;

        for (size_t i = 3; i < (max_value >> 4) + 4; ++i) cout << "tiny_bmp[" << i << "] = " << tiny_bmp[i] << endl;
    }
    else {

        for (size_t i = 3; i < cardinality + 3; ++i) cout << "tiny_bmp[" << i << "] = " << tiny_bmp[i] << endl;
    }
}

bool TinyBitmap::change_sz(const uint16_t sz_min) {

    if (sz_min > sizes[nb_sizes - 1]) return false;

    const bool is_allocated = (tiny_bmp != nullptr);

    const uint16_t sz = is_allocated ? getSize() : 0;
    const uint16_t new_sz = getNextSize(sz_min);

    if (is_allocated){

        uint16_t* tiny_bmp_new = nullptr;

        const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp_new), 8, new_sz * sizeof(uint16_t));

        if (aligned_alloc != 0){

            cerr << "TinyBitmap::change_sz(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
            exit(1);
        }

        std::memset(tiny_bmp_new, 0, new_sz * sizeof(uint16_t));
        std::copy(tiny_bmp, tiny_bmp + (new_sz >= sz ? sz : sz_min), tiny_bmp_new);

        free(tiny_bmp);

        tiny_bmp = tiny_bmp_new;
        tiny_bmp[0] = (tiny_bmp[0] & ~sz_mask) | (new_sz << 3);
    }
    else {

        const int aligned_alloc = posix_memalign(reinterpret_cast<void**>(&tiny_bmp), 8, new_sz * sizeof(uint16_t));

        if (aligned_alloc != 0){

            cerr << "TinyBitmap::change_sz(): Aligned memory could not be allocated with error " << aligned_alloc << endl;
            exit(1);
        }

        std::memset(tiny_bmp, 0, new_sz * sizeof(uint16_t));

        tiny_bmp[0] = bmp_mode | (new_sz << 3);
    }

    return true;
}

bool TinyBitmap::switch_mode(const uint16_t sz_min, const uint16_t new_mode) {

    if (tiny_bmp == nullptr) return true;

    const uint16_t sz = getSize();
    const uint16_t mode = getMode();
    const uint16_t offset = getOffset();
    const uint16_t cardinality = getCardinality();

    uint16_t* tiny_bmp_new = nullptr;

    if ((mode == bmp_mode) && (new_mode == list_mode)) { // Switch to list mode

        const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>(cardinality + 3));

        if (new_sz_min > sizes[nb_sizes - 1]) return false;

        std::swap(tiny_bmp, tiny_bmp_new);
        change_sz(new_sz_min);

        for (uint16_t i_bmp = 3, i_list = 3, card = cardinality; (i_bmp < sz) && (card != 0); ++i_bmp){

            for (uint16_t j = (i_bmp - 3) << 4; tiny_bmp_new[i_bmp] != 0; tiny_bmp_new[i_bmp] >>= 1, ++j){

                if (tiny_bmp_new[i_bmp] & 0x1){

                    tiny_bmp[i_list++] = j;
                    --card;
                }
            }
        }

        tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | list_mode | bits_16;
        tiny_bmp[1] = cardinality;
        tiny_bmp[2] = offset;
    }
    else if ((mode == list_mode) && (new_mode == bmp_mode)) { // mode is list

        const uint16_t max_val_mod = static_cast<const uint16_t>(maximum() & 0xFFFF);
        const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>((max_val_mod >> 4) + 4));

        if (new_sz_min > sizes[nb_sizes - 1]) return false;

        std::swap(tiny_bmp, tiny_bmp_new);
        change_sz(new_sz_min);

        for (size_t i = 3; i < cardinality + 3; ++i) tiny_bmp[(tiny_bmp_new[i] >> 4) + 3] |= (1U << (tiny_bmp_new[i] & 0xF));

        tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | bmp_mode | bits_16;
        tiny_bmp[1] = cardinality;
        tiny_bmp[2] = offset;
    }
    else if (mode == rle_list_mode){

        if (new_mode == bmp_mode){

            uint16_t new_card = 0;

            const uint16_t max_val_mod = static_cast<const uint16_t>(maximum() & 0xFFFF);
            const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>((max_val_mod >> 4) + 4));

            if (new_sz_min > sizes[nb_sizes - 1]) return false;

            std::swap(tiny_bmp, tiny_bmp_new);
            change_sz(new_sz_min);

            for (size_t i = 3; i < cardinality + 3; i += 2){

                new_card += tiny_bmp_new[i + 1] - tiny_bmp_new[i] + 1;

                for (uint16_t j = tiny_bmp_new[i]; j <= tiny_bmp_new[i + 1]; ++j) tiny_bmp[(j >> 4) + 3] |= (1U << (j & 0xF));
            }

            tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | bmp_mode | bits_16;
            tiny_bmp[1] = new_card;
            tiny_bmp[2] = offset;
        }
        else if (new_mode == list_mode){

            const uint16_t new_card = size();
            const uint16_t new_sz_min = std::max(sz_min, static_cast<const uint16_t>(new_card + 3));

            if (new_sz_min > sizes[nb_sizes - 1]) return false;

            std::swap(tiny_bmp, tiny_bmp_new);
            change_sz(new_sz_min);

            for (size_t i = 3, k = 3; i < cardinality + 3; i += 2){

                for (uint32_t j = tiny_bmp_new[i]; j <= tiny_bmp_new[i + 1]; ++j) tiny_bmp[k++] = j;
            }

            tiny_bmp[0] = (tiny_bmp[0] & sz_mask) | list_mode | bits_16;
            tiny_bmp[1] = new_card;
            tiny_bmp[2] = offset;
        }
    }

    if (tiny_bmp_new != nullptr) free(tiny_bmp_new);

    return true;
}

TinyBitmap::const_iterator TinyBitmap::begin() const {

    const_iterator it(*this, true);
    ++it;
    return it;
}

TinyBitmap::const_iterator TinyBitmap::end() const {

    return const_iterator(*this, false);
}

bool TinyBitmap::test(const bool verbose) {

    TinyBitmap t_bmp;

    auto check_cardinality = [&](){

        if ((t_bmp.tiny_bmp != nullptr) && (t_bmp.getMode() != bmp_mode) && (t_bmp.getCardinality() + 3 > t_bmp.getSize())){

            cout << "TinyBitmap::test(): cardinality (" << t_bmp.getCardinality() << ") + 3 > sz (" << t_bmp.getSize() << ")" << endl;
            t_bmp.print();
            exit(1);
        }
    };

    auto check_sorting = [&](){

        if ((t_bmp.tiny_bmp != nullptr) && (t_bmp.getMode() != bmp_mode)){

            for (size_t i = 4; i < t_bmp.getCardinality() + 3; ++i){

                if (t_bmp.tiny_bmp[i] < t_bmp.tiny_bmp[i-1]){

                    cout << "TinyBitmap::test(): Not sorted " << endl;

                    t_bmp.print();
                    exit(1);
                }
            }
        }
    };

    const size_t nb_rounds = 10;

    if (verbose) cout << "TinyBitmap::test(): Adding values in sequential order from 0 to 65536-49" << endl;

    for (uint32_t i = 0; i != 65536-48; ++i){

        if (t_bmp.add(i) != t_bmp.contains(i)){

            if (verbose) cerr << "TinyBitmap::test(): Error while adding values" << endl;
            return false;
        }

        check_cardinality();
        check_sorting();
    }

    if (verbose) cout << "TinyBitmap::test(): Iterating over the values" << endl;

    for (const auto val : t_bmp) {}

    if (verbose) cout << "TinyBitmap::test(): Deleting values in sequential order from 0 to 65536-49" << endl;

    for (uint32_t i = 0; i != 65536-48; ++i){

        t_bmp.remove(i);

        check_cardinality();
        check_sorting();

        if (t_bmp.contains(i)){

            if (verbose) cerr << "TinyBitmap::test(): Error while removing values" << endl;
            return false;
        }
    }

    t_bmp.clear();

    for (size_t j = 0; j < nb_rounds; ++j){

        if (verbose) cout << "TinyBitmap::test(): Adding values in random order from 0 to 65536-49 (round " << j << ")" << endl;

        vector<uint32_t> val_added;

        for (uint32_t i = 0; i != 65536 - 48; ++i){

            const uint32_t val = rand() % (65536 - 48);

            val_added.push_back(val);

            if (t_bmp.add(val) != t_bmp.contains(val)){

                if (verbose) cerr << "TinyBitmap::test(): Error while adding values" << endl;
                return false;
            }

            check_cardinality();
            check_sorting();
        }

        if (verbose) cout << "TinyBitmap::test(): Iterating over the values" << endl;

        for (const auto val : t_bmp) {}

        if (verbose) cout << "TinyBitmap::test(): Removing values in random order from 0 to 65536-49 (round " << j << ")" << endl;

        std::random_shuffle(val_added.begin(), val_added.end());

        for (const auto val : val_added){

            t_bmp.remove(val);

            check_cardinality();
            check_sorting();

            if (t_bmp.contains(val)){

                if (verbose) cerr << "TinyBitmap::test(): Error while removing values" << endl;
                return false;
            }
        }

        t_bmp.clear();
    }

    if (verbose) cout << "TinyBitmap::test(): Adding values in sequential order from 65536 to 65536 + 4096 - 3" << endl;

    for (uint32_t i = 65536; i != 65536 + 4096 - 3; ++i){

        if (t_bmp.add(i) != t_bmp.contains(i)){

            if (verbose) cerr << "TinyBitmap::test(): Error while adding values" << endl;
            return false;
        }

        check_cardinality();
        check_sorting();
    }

    if (verbose) cout << "TinyBitmap::test(): Iterating over the values" << endl;

    for (const auto val : t_bmp) {}

    if (verbose) cout << "TinyBitmap::test(): Removing values in sequential order from 65536 to 65536 + 4096 - 3" << endl;

    for (uint32_t i = 65536; i != 65536 + 4096 - 3; ++i){

        t_bmp.remove(i);

        check_cardinality();
        check_sorting();

        if (t_bmp.contains(i)){

            if (verbose) cerr << "TinyBitmap::test(): Error while removing values" << endl;
            return false;
        }
    }

    t_bmp.clear();

    if (verbose) cout << "TinyBitmap::test(): Adding values in reverse sequential order from 65536 + 4096 - 4 to 65536" << endl;

    for (uint32_t i = 65536 + 4096 - 4; i >= 65536; --i){

        if (t_bmp.add(i) != t_bmp.contains(i)){

            if (verbose) cerr << "TinyBitmap::test(): Error while adding values" << endl;
            return false;
        }

        check_cardinality();
        check_sorting();
    }

    if (verbose) cout << "TinyBitmap::test(): Iterating over the values" << endl;

    for (const auto val : t_bmp) {}

    if (verbose) cout << "TinyBitmap::test(): Removing values in reverse sequential order from 65536 + 4096 - 4 to 65536" << endl;

    for (uint32_t i = 65536 + 4096 - 4; i >= 65536; --i){

        t_bmp.remove(i);

        check_cardinality();
        check_sorting();

        if (t_bmp.contains(i)){

            if (verbose) cerr << "TinyBitmap::test(): Error while removing values" << endl;
            return false;
        }
    }

    t_bmp.clear();

    for (size_t j = 0; j < nb_rounds; ++j){

        if (verbose) cout << "TinyBitmap::test(): Adding values in random order (round " << j << ")" << endl;

        vector<uint32_t> val_added;

        for (uint32_t i = 0; i != 4093; ++i){

            const uint32_t val = rand() % 0xFFFFUL;

            val_added.push_back(val);

            if (t_bmp.add(val) != t_bmp.contains(val)){

                if (verbose) cerr << "TinyBitmap::test(): Error while adding values" << endl;
                return false;
            }

            check_cardinality();
            check_sorting();
        }

        if (verbose) cout << "TinyBitmap::test(): Iterating over the values" << endl;

        for (const auto val : t_bmp) {}

        if (verbose) cout << "TinyBitmap::test(): Removing values in random order (round " << j << ")" << endl;

        std::random_shuffle(val_added.begin(), val_added.end());

        for (const auto val : val_added){

            t_bmp.remove(val);

            check_cardinality();
            check_sorting();

            if (t_bmp.contains(val)){

                if (verbose) cerr << "TinyBitmap::test(): Error while removing values" << endl;
                return false;
            }
        }

        t_bmp.clear();
    }

    return true;
}

TinyBitmap::TinyBitmapIterator::TinyBitmapIterator() :  sz(0), mode(0), card(0), offset(0), i(0xFFFF), j(0xFFFF), e(0xFFFF),
                                                        val(0xFFFFFFFF), invalid(true), tiny_bmp(nullptr) {}

TinyBitmap::TinyBitmapIterator::TinyBitmapIterator(const TinyBitmap& t_bmp, const bool start) :    sz(0), mode(0), card(0), offset(0), i(0xFFFF), j(0xFFFF),
                                                                                                    e(0xFFFF), val(0xFFFFFFFF), invalid(true),
                                                                                                    tiny_bmp(start ? t_bmp.tiny_bmp : nullptr) {

    if (start){

        sz = t_bmp.getSize();
        mode = t_bmp.getMode();
        card = t_bmp.getCardinality();

        offset = static_cast<uint32_t>(t_bmp.getOffset()) << 16;

        if (card != 0){

            i = 2;
            invalid = false;

            if (mode == bmp_mode) e = 0;
            else if (mode == rle_list_mode){

                i = 3;
                j = 4;
                val = (offset | static_cast<uint32_t>(tiny_bmp[i])) - 1;
            }
        }
    }
}

TinyBitmap::TinyBitmapIterator& TinyBitmap::TinyBitmapIterator::operator++() {

    if (invalid) return *this; // Iterator has ended

    if (mode == bmp_mode){

        ++j;
        e >>= 1;

        if (e == 0){

            ++i;
            j = 0;
            e = (i == sz) ? 0 : tiny_bmp[i];
        }

        while ((i != sz) && (card != 0)){

            for (; e != 0; e >>= 1, ++j){

                if (e & 0x1){

                    val = offset | (((i - 3) << 4) + j);
                    --card;
                    return *this;
                }
            }

            ++i;
            e = (i == sz) ? 0 : tiny_bmp[i];
        }

        invalid = true;
    }
    else if (mode == list_mode){

        ++i;

        if (i >= card + 3) invalid = true;
        else val = offset | tiny_bmp[i];
    }
    else {

        ++val;

        if ((val & 0xFFFF0000) != offset) invalid = true;
        else if ((val & 0xFFFF) > tiny_bmp[j]) {

            i += 2;
            j += 2;

            if (i >= card + 3) invalid = true;
            else val = offset | tiny_bmp[i];
        }
    }

    return *this;
}

TinyBitmap::TinyBitmapIterator TinyBitmap::TinyBitmapIterator::operator++(int){

    TinyBitmapIterator tmp(*this);
    operator++();
    return tmp;
}

const uint16_t TinyBitmap::sz_mask = 0xFFF8;    // 1111111111111000
const uint16_t TinyBitmap::mode_mask = 0x0006;  // 0000000000000110
const uint16_t TinyBitmap::bits_mask = 0x0001;  // 0000000000000001

const uint16_t TinyBitmap::bmp_mode = 0x0000;
const uint16_t TinyBitmap::list_mode = 0x0002;
const uint16_t TinyBitmap::rle_list_mode = 0x0004;

const uint16_t TinyBitmap::bits_16 = 0x0000;
const uint16_t TinyBitmap::bits_32 = 0x0001;

const uint16_t TinyBitmap::sizes[] = {  8, 16, 32, 64, 96, 160, 224, 320,
                                        448, 768, 1024, 1280, 1792, 2048, 2560, 3072,
                                        4096, 0xFFFF };

const uint16_t TinyBitmap::nb_sizes = 17;

/*const uint16_t TinyBitmap::sizes[] = {  8, 16, 24, 32, 40, 48, 56, 64,
                                        80, 96, 112, 128, 160, 192, 224, 256,
                                        320, 384, 448, 512, 640, 768, 896, 1024,
                                        1280, 1536, 1792, 2048, 2560, 3072, 3584, 4096,
                                        0xFFFF};

const uint16_t TinyBitmap::nb_sizes = 32;*/
