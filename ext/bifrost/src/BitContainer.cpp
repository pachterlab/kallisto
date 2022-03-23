#include "BitContainer.hpp"

BitContainer::BitContainer() : setBits(localBitVector) {}

BitContainer::BitContainer(const BitContainer& o) {

    const uintptr_t flag = o.setBits & flagMask;

    if (flag == ptrBitmap){

        Bitmap* setPtrBmp = new Bitmap;

        setPtrBmp->r = o.getConstPtrBitmap()->r;

        setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = o.getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        TinyBitmap t_bmp_cpy(t_bmp);

        t_bmp.detach();

        setBits = (reinterpret_cast<uintptr_t>(t_bmp_cpy.detach()) & pointerMask) | localTinyBitmap;
    }
    else setBits = o.setBits;
}

BitContainer::BitContainer(BitContainer&& o) : setBits(o.setBits) {

    o.setBits = localBitVector;
}

BitContainer::~BitContainer() {

    releaseMemory();
}

BitContainer& BitContainer::operator=(const BitContainer& o){

    if (this != &o) {

        const uintptr_t flag = o.setBits & flagMask;
        const uintptr_t this_flag = setBits & flagMask;

        if (flag == ptrBitmap){

            Bitmap* setPtrBmp = nullptr;

            if (this_flag == ptrBitmap) setPtrBmp = getPtrBitmap();
            else {

                releaseMemory();
                setPtrBmp = new Bitmap;
            }

            setPtrBmp->r = o.getConstPtrBitmap()->r;

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
        }
        else if (flag == localTinyBitmap){

            releaseMemory();

            uint16_t* setPtrTinyBmp = o.getPtrTinyBitmap();

            TinyBitmap t_bmp(&setPtrTinyBmp);
            TinyBitmap t_bmp_cpy(t_bmp);

            t_bmp.detach();

            setBits = (reinterpret_cast<uintptr_t>(t_bmp_cpy.detach()) & pointerMask) | localTinyBitmap;
        }
        else {

            releaseMemory();

            setBits = o.setBits;
        }
    }

    return *this;
}

BitContainer& BitContainer::operator=(BitContainer&& o){

    if (this != &o) {

        releaseMemory();

        setBits = o.setBits;
        o.setBits = localBitVector;
    }

    return *this;
}

bool BitContainer::operator==(const BitContainer& o) const {

    if (size() != o.size()) return false;

    BitContainer::const_iterator it(begin()), o_it(o.begin());
    const BitContainer::const_iterator it_end(end()), o_it_end(o.end());

    for (; (it != it_end) && (o_it != o_it_end); ++it, ++o_it){

        if (*it != *o_it) return false;
    }

    return ((it == it_end) && (o_it == o_it_end));
}

void BitContainer::clear(){

    releaseMemory();
    setBits = localBitVector;
}

size_t BitContainer::getSizeInBytes() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.getSizeInBytes() + sizeof(Bitmap) + sizeof(BitContainer);

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.getSizeInBytes();

        t_bmp.detach();

        return ret;
    }

    return sizeof(BitContainer); // if (flag == ptrSharedBitContainer), we do not own the pointed UnitigCors so its size is not considered
}

BitContainer BitContainer::operator|(const BitContainer& rhs) const {

    BitContainer lhs(*this);

    lhs |= rhs;
    return lhs;
}

BitContainer& BitContainer::operator|=(const BitContainer& rhs) {

    if (&rhs != this){

        if (cardinality() == 0) *this = rhs;
        else if (rhs.cardinality() != 0){

            const uintptr_t flag = setBits & flagMask;
            const uintptr_t r_flag = rhs.setBits & flagMask;

            if ((flag == r_flag) && ((flag == ptrBitmap) || (flag == localBitVector))) {

                if (flag == ptrBitmap) getPtrBitmap()->r |= rhs.getConstPtrBitmap()->r;
                if (flag == localBitVector) setBits |= rhs.setBits;
            }
            else {

                vector<uint32_t> new_ids;

                const_iterator it = begin(), it_end = end();
                const_iterator r_it = rhs.begin(), r_it_end = rhs.end();

                while ((it != it_end) && (r_it != r_it_end)){

                    if (*it > *r_it){

                        new_ids.push_back(*r_it);
                        ++r_it;
                    }
                    else if (*it < *r_it) ++it;
                    else {

                        ++it;
                        ++r_it;
                    }
                }

                while (r_it != r_it_end){

                    new_ids.push_back(*r_it);
                    ++r_it;
                }

                addSortedVector(new_ids);
            }
        }
    }

    return *this;
}

BitContainer BitContainer::operator&(const BitContainer& rhs) const {

    BitContainer lhs(*this);

    lhs &= rhs;
    return lhs;
}

BitContainer& BitContainer::operator&=(const BitContainer& rhs) {

    if (&rhs != this){

        if ((cardinality() == 0) || (rhs.cardinality() == 0)) clear();
        else {

            vector<uint32_t> old_ids;

            const_iterator it = begin(), it_end = end();
            const_iterator r_it = rhs.begin(), r_it_end = rhs.end();

            while (r_it != r_it_end){

                while ((it != it_end) && (*it < *r_it)){

                    old_ids.push_back(*it);
                    ++it;
                }

                if ((it != it_end) && (*it == *r_it)) ++it;

                ++r_it;
            }

            while (it != it_end){

                old_ids.push_back(*it);
                ++it;
            }

            removeSortedVector(old_ids);
        }
    }

    return *this;
}

BitContainer BitContainer::operator-(const BitContainer& rhs) const {

    BitContainer lhs(*this);

    lhs -= rhs;
    return lhs;
}

BitContainer& BitContainer::operator-=(const BitContainer& rhs) {

    if ((&rhs != this) && (rhs.cardinality() != 0)) {

        const_iterator it = begin(), it_end = end();
        const_iterator r_it = rhs.begin(), r_it_end = rhs.end();

        vector<uint32_t> old_ids;

        old_ids.reserve(cardinality());

        while ((it != it_end) && (r_it != r_it_end)){

            if (*it > *r_it) ++r_it;
            else if (*it < *r_it) ++it;
            else {

                old_ids.push_back(*r_it);

                ++it;
                ++r_it;
            }
        }

        removeSortedVector(old_ids);
    }

    return *this;
}

Roaring BitContainer::toRoaring() const {

    uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getPtrBitmap()->r;
    else {

        Roaring b;

        const_iterator it = begin();
        const_iterator it_end = end();

        while (it != it_end){

            b.add(*it);
            ++it;
        }

        return b;
    }
}

void BitContainer::add(const size_t pair_id) {

    uintptr_t flag = setBits & flagMask;

    if (flag == localSingleInt){

        const uintptr_t setBits_tmp = setBits >> shiftMaskBits;

        if (setBits_tmp != pair_id){

            if ((setBits_tmp < maxBitVectorIDs) && (pair_id < maxBitVectorIDs)){

                setBits = (1ULL << (setBits_tmp + shiftMaskBits)) | localBitVector;
            }
            else {

                TinyBitmap t_bmp;

                if (t_bmp.add(setBits_tmp)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                else {

                    t_bmp.clear();

                    Bitmap* setPtrBmp = new Bitmap;

                    setPtrBmp->r.add(setBits_tmp);

                    setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
                }
            }


            flag = setBits & flagMask;
        }
    }

    if (flag == localBitVector){

        if (setBits == localBitVector) setBits = (pair_id << shiftMaskBits) | localSingleInt;
        else if (pair_id < maxBitVectorIDs) setBits |= 1ULL << (pair_id + shiftMaskBits);
        else {

            uintptr_t setBits_tmp_tb = setBits >> shiftMaskBits;
            uintptr_t setBits_tmp_cr = setBits >> shiftMaskBits;

            TinyBitmap t_bmp;

            bool add_ok = true;

            for (size_t i = 0; (setBits_tmp_tb != 0) && add_ok; ++i, setBits_tmp_tb >>= 1) {

                if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(i);
            }

            if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                t_bmp.clear();

                Bitmap* setPtrBmp = new Bitmap;

                for (size_t i = 0; setBits_tmp_cr != 0; ++i, setBits_tmp_cr >>= 1) {

                    if (setBits_tmp_cr & 0x1) setPtrBmp->r.add(i);
                }

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localTinyBitmap) {

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        if (t_bmp.add(pair_id)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        else {

            const size_t sz_t_bmp = t_bmp.size();

            size_t i = 0;

            uint32_t* values = new uint32_t[sz_t_bmp];

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[i++] = *it;

            t_bmp.clear();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) getPtrBitmap()->r.add(pair_id); // flag == ptrBitmap
}

void BitContainer::addSortedVector(const vector<uint32_t>& v) { // Private, assumes vector is sorted

    if (v.empty()) return;

    uintptr_t flag = setBits & flagMask;

    size_t i = 0;

    if (flag == localSingleInt){

        const uintptr_t setBits_tmp = setBits >> shiftMaskBits;

        if ((setBits_tmp < maxBitVectorIDs) && (v[v.size() - 1] < maxBitVectorIDs)){

            setBits = (1ULL << (setBits_tmp + shiftMaskBits)) | localBitVector;
        }
        else {

            TinyBitmap t_bmp;

            if (t_bmp.add(setBits_tmp)) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                t_bmp.clear();

                Bitmap* setPtrBmp = new Bitmap;

                setPtrBmp->r.add(setBits_tmp);

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localBitVector){

        if ((setBits == localBitVector) && (v.size() == 1)) setBits = (v[0] << shiftMaskBits) | localSingleInt;
        else if (v[v.size() - 1] < maxBitVectorIDs){

            for (const uint32_t id : v) setBits |= 1ULL << (id + shiftMaskBits);
        }
        else {

            uintptr_t setBits_tmp_tb = setBits >> shiftMaskBits;
            uintptr_t setBits_tmp_cr = setBits >> shiftMaskBits;

            TinyBitmap t_bmp;

            bool add_ok = true;

            for (size_t j = 0; (setBits_tmp_tb != 0) && add_ok; ++j, setBits_tmp_tb >>= 1) {

                if (setBits_tmp_tb & 0x1) add_ok = t_bmp.add(j);
            }

            if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
            else {

                Bitmap* setPtrBmp = new Bitmap;

                t_bmp.clear();

                for (size_t j = 0; setBits_tmp_cr != 0; ++j, setBits_tmp_cr >>= 1) {

                    if (setBits_tmp_cr & 0x1) setPtrBmp->r.add(j);
                }

                setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            }
        }

        flag = setBits & flagMask;
    }

    if (flag == localTinyBitmap) {

        bool add_ok = true;

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        for (i = 0; (i < v.size()) && add_ok; ++i) add_ok = t_bmp.add(v[i]);

        if (add_ok) setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t j = 0;

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[j++] = *it;

            t_bmp.clear();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        bitmap->r.addMany(v.size() - i, &v[i]);
        bitmap->r.runOptimize();
    }
}

void BitContainer::remove(const size_t pair_id) {

    uintptr_t flag = setBits & flagMask;

    if (flag == localBitVector){

        if (pair_id < maxBitVectorIDs) setBits &= ~(1ULL << (pair_id + shiftMaskBits));
    }
    else if (flag == localSingleInt){

        if (pair_id == (setBits >> shiftMaskBits)) setBits = localBitVector;
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        const bool rm_ok = t_bmp.remove(pair_id);

        if (rm_ok){

            const size_t card = t_bmp.size();

            if (card == 0){

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                clear();
            }
            else if (card == 1){

                const uint32_t l_pair_id = *(t_bmp.begin());

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                clear();
                add(l_pair_id);
            }
            else if (/*(card <= maxBitVectorIDs) && */(t_bmp.maximum() < maxBitVectorIDs)){

                BitContainer new_uc;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) new_uc.add(*it);

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                *this = move(new_uc);
            }
            else setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t i = 0;

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[i++] = *it;

            t_bmp.clear();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        bitmap->r.remove(pair_id);

        const size_t card = bitmap->r.cardinality();

        if (card == 0) clear();
        else if (card == 1){

            const uint32_t l_pair_id = bitmap->r.minimum();

            clear();
            add(l_pair_id);
        }
        else if (/*(card <= maxBitVectorIDs) &&*/ (bitmap->r.maximum() < maxBitVectorIDs)){

            BitContainer new_uc;

            const_iterator it = begin(), it_end = end();

            for (; it != it_end; ++it) new_uc.add(*it);

            *this = move(new_uc);
        }
        else if ((setBits & flagMask) == ptrBitmap) bitmap->r.runOptimize();
    }
}

void BitContainer::removeSortedVector(const vector<uint32_t>& v) {

    if (v.empty()) return;

    uintptr_t flag = setBits & flagMask;

    size_t i = 0;

    if (flag == localBitVector){

        uintptr_t mask = 0;

        for (size_t j = 0; j < v.size(); ++j){

            if (v[j] < maxBitVectorIDs) mask |= 1ULL << (v[j] + shiftMaskBits);
        }

        setBits &= ~mask;
    }
    else if (flag == localSingleInt){

        const uintptr_t id = setBits >> shiftMaskBits;

        if (id <= v[v.size() - 1]){

            for (size_t j = 0; j < v.size(); ++j){

                if (v[j] == id) setBits = localBitVector;
                if (v[j] >= id) break;
            }
        }
    }
    else if (flag == localTinyBitmap){

        bool rm_ok = true;

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();

        TinyBitmap t_bmp(&setPtrTinyBmp);

        for (; (i < v.size()) && rm_ok; ++i) rm_ok = t_bmp.remove(v[i]);

        if (rm_ok){

            const size_t card = t_bmp.size();

            if (card == 0){

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
                clear();
            }
            else if (card == 1){

                const uint32_t id = *(t_bmp.begin());

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                clear();
                add(id);
            }
            else if (/*(card <= maxBitVectorIDs) &&*/ (t_bmp.maximum() < maxBitVectorIDs)){

                BitContainer new_uc;

                for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) new_uc.add(*it);

                setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;

                *this = move(new_uc);
            }
            else setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }
        else {

            const size_t sz_t_bmp = t_bmp.size();

            uint32_t* values = new uint32_t[sz_t_bmp];

            size_t j = 0;

            Bitmap* setPtrBmp = new Bitmap;

            for (TinyBitmap::const_iterator it = t_bmp.begin(), it_end = t_bmp.end(); it != it_end; ++it) values[j++] = *it;

            t_bmp.clear();
            setPtrBmp->r.addMany(sz_t_bmp, values);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;
            flag = ptrBitmap;

            delete[] values;
        }
    }

    if (flag == ptrBitmap) {

        Bitmap* bitmap = getPtrBitmap();

        for (; i < v.size(); ++i) bitmap->r.remove(v[i]);

        const size_t card = bitmap->r.cardinality();

        if (card == 0) clear();
        else if (card == 1){

            const uint32_t id = bitmap->r.minimum();

            clear();
            add(id);
        }
        else if (/*(card <= maxBitVectorIDs) &&*/ (bitmap->r.maximum() < maxBitVectorIDs)){

            BitContainer new_uc;

            for (const_iterator it = begin(), it_end = end(); it != it_end; ++it) new_uc.add(*it);

            *this = move(new_uc);
        }
        else if ((setBits & flagMask) == ptrBitmap) bitmap->r.runOptimize();
    }
}

bool BitContainer::contains(const size_t pair_id) const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.contains(pair_id);

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const bool ret = t_bmp.contains(pair_id);

        t_bmp.detach();

        return ret;
    }

    if (flag == localSingleInt) return (pair_id == (setBits >> shiftMaskBits));

    if (pair_id < maxBitVectorIDs){

        const uintptr_t setBits_tmp = 0x1ULL << (pair_id + shiftMaskBits);

        return ((setBits & setBits_tmp) != 0);
    }

    return false;
}

size_t BitContainer::maximum() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.maximum();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.maximum();

        t_bmp.detach();

        return ret;
    }

    if (flag == localSingleInt) return (setBits >> shiftMaskBits);

    const int nb_lead_0 = __builtin_clzll(setBits | flagMask);

    return (maxBitVectorIDs - nb_lead_0 - (nb_lead_0 != maxBitVectorIDs));
}

size_t BitContainer::minimum() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.minimum();

    const_iterator it = begin(), it_end = end();

    if (it != it_end) return *it;
    return 0;
}

size_t BitContainer::size() const {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) return getConstPtrBitmap()->r.cardinality();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);
        const size_t ret = t_bmp.size();

        t_bmp.detach();

        return ret;
    }

    if (flag == localBitVector) return __builtin_popcountll(setBits & pointerMask);

    return 1;
}

bool BitContainer::write(ostream& stream_out) const {

    if (stream_out.good()){

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrBitmap){

            const uint32_t expected_sz = getConstPtrBitmap()->r.getSizeInBytes();

            const uintptr_t flag_expected_sz = (static_cast<uintptr_t>(expected_sz) << shiftMaskBits) | flag;

            char* serialized = new char[expected_sz];

            getConstPtrBitmap()->r.write(serialized);

            stream_out.write(reinterpret_cast<const char*>(&flag_expected_sz), sizeof(uintptr_t));
            stream_out.write(serialized, expected_sz);

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
            TinyBitmap t_bmp(&setPtrTinyBmp);

            stream_out.write(reinterpret_cast<const char*>(&flag), sizeof(uintptr_t));

            t_bmp.write(stream_out);
            t_bmp.detach();
        }
        else stream_out.write(reinterpret_cast<const char*>(&setBits), sizeof(uintptr_t));

        return true;
    }

    return false;
}

bool BitContainer::read(istream& stream_in) {

    if (stream_in.good()){

        clear();

        stream_in.read(reinterpret_cast<char*>(&setBits), sizeof(uintptr_t));

        const uintptr_t flag = setBits & flagMask;

        if (flag == ptrBitmap){

            Bitmap* setPtrBmp = new Bitmap;

            const uint32_t expected_sz = static_cast<uint32_t>(setBits >> shiftMaskBits);

            char* serialized = new char[expected_sz];

            stream_in.read(serialized, expected_sz);

            setPtrBmp->r = Roaring::read(serialized);

            setBits = (reinterpret_cast<uintptr_t>(setPtrBmp) & pointerMask) | ptrBitmap;

            delete[] serialized;
        }
        else if (flag == localTinyBitmap){

            TinyBitmap t_bmp;

            t_bmp.read(stream_in);

            setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
        }

        return true;
    }

    return false;
}

void BitContainer::runOptimize() {

    const uintptr_t flag = setBits & flagMask;

    if (flag == ptrBitmap) getPtrBitmap()->r.runOptimize();

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = getPtrTinyBitmap();
        TinyBitmap t_bmp(&setPtrTinyBmp);

        t_bmp.runOptimize();

        setBits = (reinterpret_cast<uintptr_t>(t_bmp.detach()) & pointerMask) | localTinyBitmap;
    }
}

BitContainer::const_iterator BitContainer::begin() const {

    const_iterator it(this, true);
    return ++it;
}

BitContainer::const_iterator BitContainer::end() const {

    return const_iterator(this, false);
}

BitContainer::BitContainer_const_iterator::BitContainer_const_iterator() :    cs(nullptr), flag(localBitVector), it_setBits(0), it_roar(empty_roar.end()),
                                                            ck_id(0xffffffffffffffff), cs_sz(0) {}

BitContainer::BitContainer_const_iterator::BitContainer_const_iterator( const BitContainer* cs_, const bool beg) :  cs(cs_), ck_id(0xffffffffffffffff), it_setBits(0xffffffffffffffff),
                                                                                            it_roar(empty_roar.end()) {

    flag = cs->setBits & flagMask;

    if (flag == ptrBitmap){

        it_roar = beg ? cs->getConstPtrBitmap()->r.begin() : cs->getConstPtrBitmap()->r.end();
        cs_sz = cs->size();
    }
    else if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();

        t_bmp = &setPtrTinyBmp;
        it_t_bmp = beg ? t_bmp.begin() : t_bmp.end();
        cs_sz = cs->size();
    }
    else cs_sz = (flag == localSingleInt) ? 1 : maxBitVectorIDs;

    if (!beg) it_setBits = cs_sz;
}

BitContainer::BitContainer_const_iterator::BitContainer_const_iterator(const BitContainer_const_iterator& o) :  cs(o.cs), flag(o.flag), it_setBits(o.it_setBits), cs_sz(o.cs_sz),
                                                                                        ck_id(o.ck_id), it_roar(o.it_roar), it_t_bmp(o.it_t_bmp) {

    if (flag == localTinyBitmap){

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();
        t_bmp = &setPtrTinyBmp;
    }
}

BitContainer::BitContainer_const_iterator::~BitContainer_const_iterator() {

    t_bmp.detach();
}

BitContainer::BitContainer_const_iterator& BitContainer::BitContainer_const_iterator::operator=(const BitContainer_const_iterator& o) {

    cs = o.cs;

    flag = o.flag;

    it_setBits = o.it_setBits;
    cs_sz = o.cs_sz;

    ck_id = o.ck_id;

    it_roar = o.it_roar;
    it_t_bmp = o.it_t_bmp;

    if (flag == localTinyBitmap){

        t_bmp.detach();

        uint16_t* setPtrTinyBmp = cs->getPtrTinyBitmap();
        t_bmp = &setPtrTinyBmp;
    }

    return *this;
}

BitContainer::BitContainer_const_iterator BitContainer::BitContainer_const_iterator::operator++(int) {

    BitContainer_const_iterator tmp(*this);
    operator++();
    return tmp;
}

BitContainer::BitContainer_const_iterator& BitContainer::BitContainer_const_iterator::operator++() {

    if (it_setBits != cs_sz){

        ++it_setBits;

        if (flag == ptrBitmap) {

            if (it_setBits != 0) ++it_roar;
            if (it_roar != cs->getConstPtrBitmap()->r.end()) ck_id = *it_roar;
        }
        else if (flag == localTinyBitmap) {

            if (it_setBits != 0) ++it_t_bmp;
            if (it_t_bmp != t_bmp.end()) ck_id = *it_t_bmp;
        }
        else if (flag == localBitVector){

            while (it_setBits < maxBitVectorIDs){

                if (((cs->setBits >> (it_setBits + shiftMaskBits)) & 0x1) != 0){

                    ck_id = it_setBits;
                    break;
                }

                ++it_setBits;
            }
        }
        else ck_id = cs->setBits >> shiftMaskBits;
    }

    return *this;
}

bool BitContainer::BitContainer_const_iterator::operator==(const BitContainer_const_iterator& o) const {

    if ((cs == o.cs) && (flag == o.flag) && (cs_sz == o.cs_sz)){

        if (flag == ptrBitmap) return (it_roar == o.it_roar);
        if (flag == localTinyBitmap) return (it_t_bmp == o.it_t_bmp);

        return (it_setBits == o.it_setBits);
    }

    return false;
}

bool BitContainer::BitContainer_const_iterator::operator!=(const BitContainer_const_iterator& o) const {

    return !operator==(o);
}

const size_t BitContainer::maxBitVectorIDs = 61; // 64 bits - 3 bits for the color set type

const uintptr_t BitContainer::localTinyBitmap = 0x0;
const uintptr_t BitContainer::localBitVector = 0x1;
const uintptr_t BitContainer::localSingleInt = 0x2;
const uintptr_t BitContainer::ptrBitmap = 0x3;

const size_t BitContainer::shiftMaskBits = 3;

const uintptr_t BitContainer::flagMask = 0x7;
const uintptr_t BitContainer::pointerMask = 0xFFFFFFFFFFFFFFF8;
