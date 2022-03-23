#include "Kmer.hpp"

using namespace std;

static const uint64_t twin_table[256] = {
  0xFF, 0xBF, 0x7F, 0x3F, 0xEF, 0xAF, 0x6F, 0x2F,
  0xDF, 0x9F, 0x5F, 0x1F, 0xCF, 0x8F, 0x4F, 0x0F,
  0xFB, 0xBB, 0x7B, 0x3B, 0xEB, 0xAB, 0x6B, 0x2B,
  0xDB, 0x9B, 0x5B, 0x1B, 0xCB, 0x8B, 0x4B, 0x0B,
  0xF7, 0xB7, 0x77, 0x37, 0xE7, 0xA7, 0x67, 0x27,
  0xD7, 0x97, 0x57, 0x17, 0xC7, 0x87, 0x47, 0x07,
  0xF3, 0xB3, 0x73, 0x33, 0xE3, 0xA3, 0x63, 0x23,
  0xD3, 0x93, 0x53, 0x13, 0xC3, 0x83, 0x43, 0x03,
  0xFE, 0xBE, 0x7E, 0x3E, 0xEE, 0xAE, 0x6E, 0x2E,
  0xDE, 0x9E, 0x5E, 0x1E, 0xCE, 0x8E, 0x4E, 0x0E,
  0xFA, 0xBA, 0x7A, 0x3A, 0xEA, 0xAA, 0x6A, 0x2A,
  0xDA, 0x9A, 0x5A, 0x1A, 0xCA, 0x8A, 0x4A, 0x0A,
  0xF6, 0xB6, 0x76, 0x36, 0xE6, 0xA6, 0x66, 0x26,
  0xD6, 0x96, 0x56, 0x16, 0xC6, 0x86, 0x46, 0x06,
  0xF2, 0xB2, 0x72, 0x32, 0xE2, 0xA2, 0x62, 0x22,
  0xD2, 0x92, 0x52, 0x12, 0xC2, 0x82, 0x42, 0x02,
  0xFD, 0xBD, 0x7D, 0x3D, 0xED, 0xAD, 0x6D, 0x2D,
  0xDD, 0x9D, 0x5D, 0x1D, 0xCD, 0x8D, 0x4D, 0x0D,
  0xF9, 0xB9, 0x79, 0x39, 0xE9, 0xA9, 0x69, 0x29,
  0xD9, 0x99, 0x59, 0x19, 0xC9, 0x89, 0x49, 0x09,
  0xF5, 0xB5, 0x75, 0x35, 0xE5, 0xA5, 0x65, 0x25,
  0xD5, 0x95, 0x55, 0x15, 0xC5, 0x85, 0x45, 0x05,
  0xF1, 0xB1, 0x71, 0x31, 0xE1, 0xA1, 0x61, 0x21,
  0xD1, 0x91, 0x51, 0x11, 0xC1, 0x81, 0x41, 0x01,
  0xFC, 0xBC, 0x7C, 0x3C, 0xEC, 0xAC, 0x6C, 0x2C,
  0xDC, 0x9C, 0x5C, 0x1C, 0xCC, 0x8C, 0x4C, 0x0C,
  0xF8, 0xB8, 0x78, 0x38, 0xE8, 0xA8, 0x68, 0x28,
  0xD8, 0x98, 0x58, 0x18, 0xC8, 0x88, 0x48, 0x08,
  0xF4, 0xB4, 0x74, 0x34, 0xE4, 0xA4, 0x64, 0x24,
  0xD4, 0x94, 0x54, 0x14, 0xC4, 0x84, 0x44, 0x04,
  0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20,
  0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00
};

Kmer::Kmer() {

    for (size_t i = 0; i < MAX_K/32; ++i) longs[i] = 0;
}

Kmer::Kmer(const Kmer& o) {

    for (size_t i = 0; i < MAX_K/32; ++i) longs[i] = o.longs[i];
}

Kmer::Kmer(const char *s) {

    set_kmer(s);
}

Kmer& Kmer::operator=(const Kmer& o) {

    if (this != &o) {

        for (size_t i = 0; i < MAX_K/32; ++i) longs[i] = o.longs[i];
    }

    return *this;
}

bool Kmer::operator<(const Kmer& o) const {

    const size_t end = MAX_K/32;

    size_t i = 0;

    while ((i < end) && (longs[i] == o.longs[i])) ++i;

    return ((i < end) && (longs[i] < o.longs[i]));
}

bool Kmer::operator==(const Kmer& o) const {

    const size_t end = MAX_K/32;

    size_t i = 0;

    while ((i < end) && (longs[i] == o.longs[i])) ++i;

    return (i == end);
}

bool Kmer::operator!=(const Kmer& o) const {

    return !(*this == o);
}

void Kmer::set_kmer(const char* s)  {

    for (size_t i = 0; i < MAX_K/32; ++i) longs[i] = 0;

    for (size_t i = 0, j, l; i < k; ++i) {

        j = 62 - ((i & 0x1F) << 1);
        l = i >> 5;

        const size_t x = ((*s) & 4) >> 1;

        longs[l] |= ((x + ((x ^ (*s & 2)) >> 1)) << j);

        ++s;
    }
}

Kmer Kmer::rep() const {

    const Kmer tw = twin();

    return (tw < *this) ? tw : *this;
}

Kmer Kmer::twin() const {

    Kmer km(*this);

    const size_t nlongs = (k+31)/32;

    for (size_t i = 0; i < nlongs; ++i) {

        const uint64_t v = longs[i];

        km.longs[nlongs-1-i] =
          (twin_table[v & 0xFF] << 56) |
          (twin_table[(v>>8) & 0xFF] << 48) |
          (twin_table[(v>>16) & 0xFF] << 40) |
          (twin_table[(v>>24) & 0xFF] << 32) |
          (twin_table[(v>>32) & 0xFF] << 24) |
          (twin_table[(v>>40) & 0xFF] << 16) |
          (twin_table[(v>>48) & 0xFF] << 8)  |
          (twin_table[(v>>56)]);
    }

    const size_t mod = (k & 0x1f) << 1; // (k % 32) * 2
    const size_t shift = (64 - mod) & 0x3f; // mod ? 64 - mod : 0
    const size_t mod2 = 64 - shift;
    const uint64_t shiftmask = (static_cast<uint64_t>(!mod) - 1) & (((1ULL << shift) - 1) << mod2);

    km.longs[0] <<= shift;

    for (size_t i = 1; i < nlongs; ++i) {

        km.longs[i-1] |= (km.longs[i] & shiftmask) >> mod2;
        km.longs[i] <<= shift;
    }

    return km;
}

Kmer Kmer::getLink(const size_t index) const {

    char c;

    switch (index % 4) {

        case 0: c = 'A'; break;
        case 1: c = 'C'; break;
        case 2: c = 'G'; break;
        case 3: c = 'T'; break;
    }

    return (index < 4) ? forwardBase(c) : backwardBase(c);
}

Kmer Kmer::forwardBase(const char b) const {

    Kmer km(*this);

    const size_t nlongs = (k+31)/32;

    km.longs[0] <<= 2;

    for (size_t i = 1; i < nlongs; ++i) {

        km.longs[i-1] |= km.longs[i] >> 62;
        km.longs[i] <<= 2;
    }

    const uint64_t x = (b & 4) >> 1;

    km.longs[nlongs-1] |= (x + ((x ^ (b & 2)) >> 1)) << ((31-((k-1) & 0x1f)) << 1);

    return km;
}

void Kmer::selfForwardBase(const char b) {

    const size_t nlongs = (k+31)/32;

    longs[0] <<= 2;

    for (size_t i = 1; i < nlongs; ++i) {

        longs[i-1] |= longs[i] >> 62;
        longs[i] <<= 2;
    }

    const uint64_t x = (b & 4) >>1;

    longs[nlongs-1] |= (x + ((x ^ (b & 2)) >>1 )) << ((31-((k-1) & 0x1f)) << 1);
}

Kmer Kmer::backwardBase(const char b) const {

    const size_t nlongs = (k+31)/32 - 1;

    Kmer km(*this);

    km.longs[nlongs] >>= 2;
    km.longs[nlongs] &= (k & 0x1f) ? (((1ULL << ((k & 0x1f) << 1)) - 1) << ((32-(k & 0x1f)) << 1)) : ~0ULL;

    for (size_t i = 1; i < nlongs + 1; ++i) {

        km.longs[nlongs-i+1] |= (km.longs[nlongs-i] & 3ULL) << 62;
        km.longs[nlongs-i] >>= 2;
    }

    const uint64_t x = (b & 4) >> 1;

    km.longs[0] |= (x + ((x ^ (b & 2)) >> 1)) << 62;

    return km;
}


// use:  km.printBinary();
// pre:
// post: The bits in the binary representation of the
//       DNA string for km has been printed to stdout
std::string Kmer::getBinary() const {

    const size_t nlongs = MAX_K/32;

    std::string r;

    r.reserve(64*nlongs);

    for (size_t i = 0; i < nlongs; i++){

        r.append(std::bitset<64>(longs[i]).to_string<char,std::char_traits<char>,std::allocator<char>>());
    }

    return r;
}


// use:  km.toString(s);
// pre:  s has space for k+1 elements
// post: s[0,...,k-1] is the DNA string for the Kmer km and s[k] = '\0'
void Kmer::toString(char *s) const {

    const size_t nlongs = (k + 31) / 32;

    for (size_t i = 0, j = 0; j < nlongs; ++j) {

        uint64_t tmp = longs[j];
        const size_t end = k < i + 32 ? k : i + 32;

        for (; i < end; ++i, ++s, tmp <<= 2){

            const char v = tmp >> 62;
            *s = 0x40 | (v + 1) | (0x1 << ((v - 1) << 1));
        }
    }

    *s = '\0';
}

char Kmer::getChar(const size_t offset) const {

    assert(offset < Kmer::k);

    const char v = (longs[offset >> 5] >> (62 - ((offset & 0x1F) << 1))) & 0x3;

    return (0x40 | (v + 1) | (0x1 << ((v - 1) << 1)));
}

bool Kmer::setChar(const size_t offset, const char b)  {

    if (offset >= Kmer::k) return false;

    const size_t pos_shift = 62 - ((offset & 0x1F) << 1);
    const size_t pos_longs = offset >> 5;

    const size_t x = (b & 4) >> 1;

    longs[pos_longs] &= 0xffffffffffffffffULL - (0x3ULL << pos_shift);
    longs[pos_longs] |= (x + ((x ^ (b & 2)) >> 1)) << pos_shift;

    return true;
}

std::string Kmer::toString() const {

    char buf[MAX_K];

    toString(buf);

    return std::string(buf);
}

// use:  set_k(k);
// pre:  this method has not been called before and 0 < k < MAX_K
// post: The Kmer size has been set to k
void Kmer::set_k(const unsigned int _k) {

    k = _k;
}

bool Kmer::write(ostream& stream_out) const {

    const size_t nlongs = MAX_K/32;

    for (size_t i = 0; (i < nlongs) && stream_out.good(); ++i){

        stream_out.write(reinterpret_cast<const char*>(&longs[i]), sizeof(uint64_t));
    }

    return stream_out.good();
}

bool Kmer::read(istream& stream_in) {

    const size_t nlongs = MAX_K/32;

    for (size_t i = 0; (i < nlongs) && stream_in.good(); ++i){

        stream_in.read(reinterpret_cast<char*>(&longs[i]), sizeof(uint64_t));
    }

    return stream_in.good();
}

unsigned int Kmer::k = 0;

Minimizer::Minimizer() {

    for (size_t i = 0; i < MAX_G/32; ++i) longs[i] = 0;
}

Minimizer::Minimizer(const Minimizer& o) {

    for (size_t i = 0; i < MAX_G/32; ++i) longs[i] = o.longs[i];
}

Minimizer::Minimizer(const char *s) {

    set_minimizer(s);
}

Minimizer& Minimizer::operator=(const Minimizer& o) {

    if (this != &o) {

        for (size_t i = 0; i < MAX_G/32; ++i) longs[i] = o.longs[i];
    }

    return *this;
}

bool Minimizer::operator<(const Minimizer& o) const {

    const size_t end = MAX_G/32;

    size_t i = 0;

    while ((i < end) && (longs[i] == o.longs[i])) ++i;

    return ((i < end) && (longs[i] < o.longs[i]));
}



// use:  b = (km1 == km2);
// pre:
// post: b is true <==> the DNA strings in km1 and km2 are equal
bool Minimizer::operator==(const Minimizer& o) const {

    const size_t end = MAX_G/32;

    size_t i = 0;

    while ((i < end) && (longs[i] == o.longs[i])) ++i;

    return (i == end);
}

bool Minimizer::operator!=(const Minimizer& o) const {

    return !(*this == o);
}

void Minimizer::set_minimizer(const char *s)  {

    for (size_t i = 0; i < MAX_G/32; ++i) longs[i] = 0;

    for (size_t i = 0, j, l; i < g; ++i) {

        j = 62 - ((i & 0x1f) << 1);
        l = i >> 5;

        const size_t x = ((*s) & 4) >> 1;

        longs[l] |= ((x + ((x ^ (*s & 2)) >> 1)) << j);

        s++;
    }
}

Minimizer Minimizer::rep() const {

    const Minimizer tw = twin();

    return (tw < *this) ? tw : *this;
}

Minimizer Minimizer::twin() const {

    const size_t nlongs = (g+31)/32;

    Minimizer minz(*this);

    for (size_t i = 0; i < nlongs; ++i) {

        const uint64_t v = longs[i];

        minz.longs[nlongs-1-i] = (twin_table[v & 0xFF] << 56) | (twin_table[(v>>8) & 0xFF] << 48) |
                                    (twin_table[(v>>16) & 0xFF] << 40) | (twin_table[(v>>24) & 0xFF] << 32) |
                                    (twin_table[(v>>32) & 0xFF] << 24) | (twin_table[(v>>40) & 0xFF] << 16) |
                                    (twin_table[(v>>48) & 0xFF] << 8)  | (twin_table[(v>>56)]);
    }

    const size_t mod = (g & 0x1f) << 1; // (g % 32) * 2
    const size_t shift = (64 - mod) & 0x3f; // mod ? 64 - mod : 0
    const size_t mod2 = 64 - shift;

    const uint64_t shiftmask = (static_cast<uint64_t>(!mod) - 1) & (((1ULL << shift) - 1) << mod2);

    minz.longs[0] <<= shift;

    for (size_t i = 1; i < nlongs; ++i) {

        minz.longs[i-1] |= (minz.longs[i] & shiftmask) >> mod2;
        minz.longs[i] <<= shift;
    }

    return minz;
}


Minimizer Minimizer::getLink(const size_t index) const {

    assert(index >= 0 && index < 8);
    char c;

    switch (index % 4) {
        case 0: c = 'A'; break;
        case 1: c = 'C'; break;
        case 2: c = 'G'; break;
        case 3: c = 'T'; break;
    }

    return (index < 4) ? forwardBase(c) : backwardBase(c);
}

Minimizer Minimizer::forwardBase(const char b) const {

    const size_t nlongs = (g+31)/32;

    Minimizer minz(*this);

    minz.longs[0] <<= 2;

    for (size_t i = 1; i < nlongs; ++i) {

        minz.longs[i-1] |= minz.longs[i] >> 62;
        minz.longs[i] <<= 2;
    }

    const uint64_t x = (b & 4) >> 1;

    minz.longs[nlongs-1] |= (x + ((x ^ (b & 2)) >> 1)) << ((31-((g-1) & 0x1f)) << 1);

    return minz;
}

Minimizer Minimizer::backwardBase(const char b) const {

    const size_t nlongs = (g+31)/32 - 1;

    Minimizer minz(*this);

    minz.longs[nlongs] >>= 2;
    minz.longs[nlongs] &= (g & 0x1f) ? (((1ULL << ((g & 0x1f) << 1)) - 1) << ((32 - (g & 0x1f)) << 1)) : ~0ULL;

    for (size_t i = 1; i < nlongs + 1; ++i) {

        minz.longs[nlongs-i+1] |= minz.longs[nlongs-i] << 62;
        minz.longs[nlongs-i] >>= 2;
    }

    const uint64_t x = (b & 4) >> 1;

    minz.longs[0] |= (x + ((x ^ (b & 2)) >> 1)) << 62;

    return minz;
}

std::string Minimizer::getBinary() const {

    size_t nlongs = MAX_G/32;

    std::string r;

    r.reserve(64*nlongs);

    for (size_t i = 0; i < nlongs; i++) {

        r.append(std::bitset<64>(longs[i]).to_string<char,std::char_traits<char>,std::allocator<char>>());
    }

    return r;
}

void Minimizer::toString(char *s) const {

    const size_t nlongs = (g + 31) / 32;

    for (size_t i = 0, j = 0; j < nlongs; ++j) {

        uint64_t tmp = longs[j];

        const size_t end = g < i + 32 ? g : i + 32;

        for (; i < end; ++i, ++s, tmp <<= 2){

            const char v = tmp >> 62;

            *s = 0x40 | (v + 1) | (0x1 << ((v - 1) << 1));
        }
    }

    *s = '\0';
}

std::string Minimizer::toString() const {

    char buf[MAX_G];

    toString(buf);

    return std::string(buf);
}

void Minimizer::set_g(unsigned int _g) {

    g = _g;
}

unsigned int Minimizer::g = 0;
