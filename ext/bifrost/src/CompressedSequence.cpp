#include <iostream>

#include "Common.hpp"
#include "CompressedSequence.hpp"
#include "Kmer.hpp"

CompressedSequence::CompressedSequence() {

    initShort();
}

CompressedSequence::~CompressedSequence() {

    clear();
}

// use:  _cs = CompressedSequence(cs);
// pre:
// post: the DNA string in _cs and is the same as in cs
CompressedSequence::CompressedSequence(const CompressedSequence& o) {

    if (o.isShort()) {

        asBits._size = o.asBits._size;
        memcpy(asBits._arr, o.asBits._arr, 31);
    }
    else setSequence(o, 0, o.size()); // copy sequence and pointers etc.
}

CompressedSequence::CompressedSequence(CompressedSequence&& o) {

    if (o.isShort()) {

        asBits._size = o.asBits._size;
        memcpy(asBits._arr, o.asBits._arr, 31); // plain vanilla copy
    }
    else {

        asPointer._length = o.asPointer._length;
        asPointer._capacity = o.asPointer._capacity;
        asPointer._data = o.asPointer._data;

        o.initShort();
    }
}

// use:  _cs = cs;
// pre:
// post: the DNA string in _cs is the same as in cs
CompressedSequence& CompressedSequence::operator=(const CompressedSequence& o) {

    if (this != &o){

        if (o.isShort()) {

            asBits._size = o.asBits._size;
            memcpy(asBits._arr, o.asBits._arr,31); // plain vanilla copy
        }
        else setSequence(o, 0, o.size()); // copy sequence and pointers etc.
    }

    return *this;
}

CompressedSequence& CompressedSequence::operator=(CompressedSequence&& o) {

    if (this != &o) {

        if (o.isShort()) {

            asBits._size = o.asBits._size;
            memcpy(asBits._arr, o.asBits._arr, 31); // plain vanilla copy
        }
        else {

            clear();

            asPointer._length = o.asPointer._length;
            asPointer._capacity = o.asPointer._capacity;
            asPointer._data = o.asPointer._data;

            o.initShort();
        }
    }

    return *this;
}

// use:  s
// pre:  s has only the characters 'A','C','G' and 'T' and can have any length
// post: the DNA string in cs is now the same as s
CompressedSequence::CompressedSequence(const char *s) {

    initShort();

    if (s != NULL) setSequence(s, strlen(s));
}


// same as with char *s but with string
CompressedSequence::CompressedSequence(const string& s) {

    initShort();

    setSequence(s.c_str(), s.length());
}


// use:  cs = CompressedSequence(km);
// pre:
// post: the DNA string in cs is now the same as the DNA string in km
CompressedSequence::CompressedSequence(const Kmer& km) {

    initShort();

    setSequence(km, Kmer::k);
}

// use:  a.setSequence(b, start, length, offset, reversed);
// pre:  start+length <= b._length, offset <= a._length
// post: copies compressed sequence from b to a (reverse complement if reversed == true)
//       the string copied from b is from [start,...,start+length-1]
//          (reverse complement of [o._length-1-start-length,...,o._length-1-start] if reversed == true)
//       the positions in a that are updated are [offset,...,offset+length-1]
//       capacity of a might be updated to fit the new string.
void CompressedSequence::setSequence(const CompressedSequence& o, const size_t start, const size_t length, const size_t offset, const bool reversed) {

    assert(length + start <= o.size());

    if (round_to_bytes(length+offset) > capacity()) _resize_and_copy(round_to_bytes(length+offset),size());

    unsigned char* data = const_cast<unsigned char*>(getPointer());
    const unsigned char *odata = o.getPointer();

    size_t w_index = offset;
    size_t wi, wj, r_index;

    if (reversed){

        r_index = o.size() - start - 1;

        for (size_t i = 0; i < length; ++i, ++w_index, --r_index) {

            wi = w_index >> 2;
            wj = (w_index & 0x3) << 1;

            data[wi] &= ~(0x3 << wj); // clear bits
            data[wi] |= (3 - ((odata[r_index >> 2] >> ((r_index & 0x3) << 1)) & 0x3)) << wj;
        }
    }
    else {

        r_index = start;

        for (size_t i = 0; i < length; ++i, ++w_index, ++r_index) {

            wi = w_index >> 2;
            wj = (w_index & 0x3) << 1;

            data[wi] &= ~(0x3 << wj); // clear bits
            data[wi] |= ((odata[r_index >> 2] >> ((r_index & 0x3) << 1)) & 0x3) << wj;
        }
    }

    // new length?
    if (offset + length > size()) setSize(offset+length);
}


// use:  cs._resize_and_copy(new_cap, copy_limit);
// pre:
// post: The DNA string in cs has space for at least new_length bases
//       the first copy_limit characters of cs are the same as before this method
void CompressedSequence::_resize_and_copy(const size_t new_cap, const size_t copy_limit) {

    if (new_cap <= capacity()) return;

    unsigned char* new_data = new unsigned char[new_cap]; // allocate new storage
    size_t bytes = round_to_bytes(copy_limit);

    memcpy(new_data, getPointer(), bytes); // copy old data

    if (isShort()) {

        size_t sz = size();

        asBits._size = 0; // this is now a long sequence.

        setSize(sz);

        asPointer._data = new_data;
        asPointer._capacity = new_cap;
    }
    else {

        delete[] asPointer._data;

        asPointer._data = new_data;
        asPointer._capacity = new_cap;
    }
}


// use:  a.setSequence(s, length, offset, reversed);
// pre:  length <= strlen(s), offset <= a._length
// post: copies substring s[0,...,length-1] to a (reverse complement if reversed == true)
//       the positions in a that are updated are [offset,...,offset+length-1]
//       capacity of a might be updated to fit the new string.
void CompressedSequence::setSequence(const char *s, const size_t length, const size_t offset, const bool reversed) {

    const size_t len = offset + length;

    if (round_to_bytes(len) > capacity()) _resize_and_copy(round_to_bytes(length + offset), size());

    unsigned char* data = const_cast<unsigned char*>(getPointer());

    if (reversed) {

        for (size_t index = offset; index < len; ++index) {

            const size_t i = index >> 2;
            const size_t j = (index & 0x3) << 1;
            const uint8_t c = bases[0x03-bits[(uint8_t)*(s+len-index-1)]];

            data[i] &= ~(0x03 << j); // set bits to 0, default
            data[i] |= (bits[c] << j);
        }
    }
    else {

        for (size_t index = offset; index < len; ++index) {

            const size_t i = index >> 2;
            const size_t j = (index & 0x3) << 1;
            const uint8_t c = *(s+index-offset);

            data[i] &= ~(0x03 << j); // set bits to 0, default
            data[i] |= (bits[c] << j);
        }
    }

    if (len > size()) setSize(len);
}


// use:  cs.setSequence(s, start, length, offset, reversed);
// pre:  0 <= start + length < o.size()
// post: If reversed is false then: cs[offset,...,offset+length-1] = s[0,...,start+length-1]
//       else: cs[offset,...,offset+length-1] is the reverse complement of s[0,...,start+length-1]
void CompressedSequence::setSequence(const string& s, const size_t length, const size_t offset, const bool reversed) {

    setSequence(s.c_str(),length,offset,reversed);
}


// use:  cs.setSequence(km, length, offset, reversed);
// pre:  0 <= length < cs._length,
//       length <= Kmer::k
// post: If reversed is false then: cs[offset,...,offset+length-1]
//         is the first length characters from the DNA string in km
//       else: cs[offset,...,offset+length-1] is the first length characters from the
//         reverse complement of the DNA string in km
void CompressedSequence::setSequence(const Kmer& km, const size_t length, const size_t offset, const bool reversed) {

    char s[Kmer::MAX_K + 1];

    km.toString(s);
    setSequence(s, length, offset, reversed);
}


// use:  s = cs.toString(offset, length);
// pre:  offset + length <= cs.size(),
// post: s is the DNA string from c[offset,...,offset+length-1]
string CompressedSequence::toString(const size_t offset, const size_t length) const {

    const unsigned char* data = getPointer();

    string s(length, 0);

    for (size_t index = offset; index < offset+length; ++index) {

        s[index-offset] = bases[(data[index >> 2] >> ((index & 0x3) << 1)) & 0x03];
    }

    return s;
}


// use:  s = cs.toString(offset, length);
// pre:  offset + length <= cs.size()
//       s has space for length characters
// post: s is the same as cs[offset,...,offset+length-1]
void CompressedSequence::toString(char *s, const size_t offset, const size_t length) const {

    const unsigned char* data = getPointer();

    for (size_t index = offset; index < offset+length; ++index) {

        s[index-offset] = bases[(data[index >> 2] >> ((index & 0x3) << 1)) & 0x03];
    }

    s[length] = 0; // 0-terminated string
}

char CompressedSequence::getChar(const size_t offset) const {

    return bases[(getPointer()[offset >> 2] >> ((offset & 0x3) << 1)) & 0x03];
}

Kmer CompressedSequence::getKmer(const size_t offset) const {

    Kmer km;

    const unsigned char* data = getPointer();

    const size_t len = offset + Kmer::k;

    if ((offset & 0x3) == 0){ // Extraction byte to byte is possible

        const size_t nbytes = (Kmer::k + 3) / 4;

        size_t i = offset >> 2, j = 0;

        for (; j < nbytes - 1; ++i, ++j) km.bytes[(7 - (j & 0x7)) + ((j >> 3) << 3)] = revBits[data[i]];

        unsigned char tmp_km = 0;
        unsigned char tmp_data = data[i];

        for (i <<= 2; i < len; ++i, tmp_data >>= 2) tmp_km = (tmp_km << 2) | (tmp_data & 0x3);

        tmp_km <<= ((4 - (Kmer::k & 0x3)) << 1) & 0x7;
        km.bytes[(7 - (j & 0x7)) + ((j >> 3) << 3)] = tmp_km;
    }
    else { // Extraction 2 bits per 2 bits

        const size_t nlongs = (Kmer::k + 31) / 32;

        size_t j = 0;

        for (size_t i = offset; j < nlongs; ++j){

            uint64_t tmp_km = 0;
            const size_t end = len < i + 32 ? len : i + 32;

            for (; i != end; ++i) tmp_km = (tmp_km << 2) | ((data[i >> 2] >> ((i & 0x3) << 1)) & 0x3);

            km.longs[j] = tmp_km;
        }

        km.longs[j-1] <<= (32 - (Kmer::k & 0x1f)) << 1;
    }

    return km;
}

bool CompressedSequence::compareKmer(const size_t offset, const size_t length, const Kmer& km) const {

    const unsigned char* data = getPointer();

    if ((length > Kmer::k) || ((offset + length) > size()) || km.isEmpty() || km.isDeleted()) return false;

    if ((offset & 0x3) == 0){ // Comparison byte to byte is possible

        const size_t nbytes = (length + 3) / 4;

        size_t i = offset >> 2, j = 0;

        while (j < nbytes - 1){ // Check full bytes (bytes containing 4 characters to compare)

            j += ((static_cast<size_t>(data[i] == revBits[km.bytes[(7 - (j & 0x7)) + ((j >> 3) << 3)]]) - 1) & nbytes) + 1;
            ++i;
        }

        if (j != nbytes - 1) return false;
        if ((length & 0x3) == 0) return (data[i] == revBits[km.bytes[(7 - (j & 0x7)) + ((j >> 3) << 3)]]); // Last byte is also a full byte

        const char mask = (0x1ULL << ((length & 0x3) << 1)) - 1;

        return (data[i] & mask) == (revBits[km.bytes[(7 - (j & 0x7)) + ((j >> 3) << 3)]] & mask);
    }
    else { //Comparison 2 bits per 2 bits

        size_t iu = offset, ik = 0;

        while (ik < length){

            const size_t cu = (data[iu >> 2] >> ((iu & 0x3) << 1)) & 0x3;
            const size_t ck = (km.longs[ik >> 5] >> (62 - ((ik & 0x1F) << 1))) & 0x3;

            ik += ((static_cast<size_t>(cu == ck) - 1) & length) + 1;
            ++iu;
        }

        return (ik == length);
    }
}

int64_t CompressedSequence::findKmer(const Kmer& km) const {

    const int k = Kmer::k;

    const size_t sz = size();

    if (sz >= k){

        Kmer km_cs = getKmer(0);
        if (km_cs == km) return 0;

        if (sz > k){

            size_t i = k;
            const unsigned char* data = getPointer();
            unsigned char tmp = data[i >> 2] >> ((i & 0x3) << 1);

            for (; i < sz; ++i, tmp >>= 2){

                if ((i & 0x3) == 0) tmp = data[i >> 2];

                km_cs.selfForwardBase(bases[tmp & 0x3]);

                if (km_cs == km) return i-k+1;
            }
        }
    }

    return -1;
}

// use:  _cs = cs.rev();
// pre:
// post: _cs is the reverse complement CompressedSequence with respect to cs,
//       i.e. if the DNA string in cs is 'GTCA'
//          then the DNA string in _cs is 'TGAC'
CompressedSequence CompressedSequence::rev() const {

    CompressedSequence r;

    r.setSequence(*this, 0, size(), 0, true);

    return r;
}


// use:  j = cs.jump(s,i,pos,reversed)
// pre:  0 <= i < s.length, -1 <= pos < cs._length if reversed true, 0 <= pos <= cs._length if reversed false
// post: if reversed == false
//         s[i...i+j-1] == cs._data[pos...pos+j-1], 0 <= j <= min(s.length-i, cs._length-pos)
//       else
//         reverse_complement(s[i...i+j-1]) == cs._data[pos-j+1...pos], 0 <= j <= min(s.length-i, pos+1)
/*size_t CompressedSequence::jump(const char *s, const size_t i, int pos, const bool reversed) const {

    assert(i >= 0);
    assert(i < strlen(s));
    assert(pos >= -1);
    assert(0 <= size() - pos); // this prevents -1 <= _length from giving false

    const unsigned char* data = getPointer();

    size_t i_cpy = i;

    if (reversed){

        if (pos == -1) return 0;

        int idx_div = pos >> 2;
        int idx_mod = (pos & 0x3) << 1;

        for (; (s[i_cpy] != '\0') && (pos != -1); --pos, ++i_cpy, idx_mod -= 2) {

            if (idx_mod == -2){
                --idx_div;
                idx_mod = 6;
            }

            if (s[i_cpy] != bases[3-((data[idx_div] >> idx_mod) & 0x03)]) break;
        }
    }
    else {

        size_t cs_size = size();

        if (pos == cs_size) return 0;

        unsigned char tmp = data[pos >> 2] >> ((pos & 0x3) << 1);

        for (; (s[i_cpy] != '\0') && (pos != cs_size); ++pos, ++i_cpy, tmp >>= 2) {

            if ((pos & 0x3) == 0) tmp = data[pos >> 2];
            if (s[i_cpy] != bases[tmp & 0x3]) break;
        }
    }

    return i_cpy - i;
}*/

size_t CompressedSequence::jump(const char *s, const size_t i, int pos, const bool reversed) const {

    const unsigned char* data = getPointer();
    const char* s_tmp = &s[i];

    if (reversed){

        for (; (*s_tmp != '\0') && (pos != -1); --pos, ++s_tmp) {

            if (*s_tmp != bases[3 - ((data[pos >> 2] >> ((pos & 0x3) << 1)) & 0x03)]) break;
        }
    }
    else {

        const size_t cs_size = size();

        for (; (*s_tmp != '\0') && (pos < cs_size); ++pos, ++s_tmp) {

            if (*s_tmp != bases[(data[pos >> 2] >> ((pos & 0x3) << 1)) & 0x3]) break;
        }
    }

    return s_tmp - &s[i];
}

/*size_t CompressedSequence::bw_jump(const char *s, const size_t i, int pos, const bool reversed) const {

    assert(i >= 0);
    assert(i < strlen(s));
    assert(pos >= -1);
    assert(0 <= size() - pos); // this prevents -1 <= _length from giving false

    const unsigned char* data = getPointer();

    size_t i_cpy = i;

    if (reversed){

        size_t cs_size = size();

        if (pos == cs_size) return 0;

        unsigned char tmp = data[pos/4] >> (2 * (pos % 4));

        for (; (i_cpy != -1) && (pos != cs_size); ++pos, --i_cpy, tmp >>= 2) {

            if (pos%4 == 0) tmp = data[pos/4];
            if (s[i_cpy] != bases[3 - (tmp & 0x3)]) break;
        }
    }
    else {

        if (pos == -1) return 0;

        int idx_div = pos / 4;
        int idx_mod = 2 * (pos % 4);

        for (; (i_cpy != -1) && (pos != -1); --pos, --i_cpy, idx_mod -= 2) {

            if (idx_mod == -2){
                idx_div--;
                idx_mod = 6;
            }

            if (s[i_cpy] != bases[(data[idx_div] >> idx_mod) & 0x03]) break;
        }
    }

    return i - i_cpy;
}*/

void CompressedSequence::clear() {

    if (!isShort() && (asPointer._capacity > 0) && (asPointer._data != NULL)) {

        delete[] asPointer._data;

        asPointer._data = NULL;
    }
}

const char CompressedSequence::bases[256] = {
    'A','C','G','T','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N',  'N','N','N','N','N','N','N','N'
};

const uint8_t CompressedSequence::bits[256] = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
};

const uint8_t CompressedSequence::revBits[256] =
{
    0x0, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0, 0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0,
    0x4, 0x44, 0x84, 0xc4, 0x14, 0x54, 0x94, 0xd4, 0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74, 0xb4, 0xf4,
    0x8, 0x48, 0x88, 0xc8, 0x18, 0x58, 0x98, 0xd8, 0x28, 0x68, 0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8,
    0xc, 0x4c, 0x8c, 0xcc, 0x1c, 0x5c, 0x9c, 0xdc, 0x2c, 0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc,
    0x1, 0x41, 0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1, 0x21, 0x61, 0xa1, 0xe1, 0x31, 0x71, 0xb1, 0xf1,
    0x5, 0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5, 0x25, 0x65, 0xa5, 0xe5, 0x35, 0x75, 0xb5, 0xf5,
    0x9, 0x49, 0x89, 0xc9, 0x19, 0x59, 0x99, 0xd9, 0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9,
    0xd, 0x4d, 0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd, 0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d, 0xbd, 0xfd,
    0x2, 0x42, 0x82, 0xc2, 0x12, 0x52, 0x92, 0xd2, 0x22, 0x62, 0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2,
    0x6, 0x46, 0x86, 0xc6, 0x16, 0x56, 0x96, 0xd6, 0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6,
    0xa, 0x4a, 0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda, 0x2a, 0x6a, 0xaa, 0xea, 0x3a, 0x7a, 0xba, 0xfa,
    0xe, 0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde, 0x2e, 0x6e, 0xae, 0xee, 0x3e, 0x7e, 0xbe, 0xfe,
    0x3, 0x43, 0x83, 0xc3, 0x13, 0x53, 0x93, 0xd3, 0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3,
    0x7, 0x47, 0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7, 0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77, 0xb7, 0xf7,
    0xb, 0x4b, 0x8b, 0xcb, 0x1b, 0x5b, 0x9b, 0xdb, 0x2b, 0x6b, 0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb,
    0xf, 0x4f, 0x8f, 0xcf, 0x1f, 0x5f, 0x9f, 0xdf, 0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff
};

