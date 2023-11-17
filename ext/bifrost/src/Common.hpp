#ifndef BIFROST_COMMON_HPP
#define BIFROST_COMMON_HPP

#include <cassert>
#include <algorithm>
#include <fstream>

#include <stdint.h>
#include <sys/stat.h>
#include <zlib.h>

#include "wyhash.h"

#if defined(__GNUC__)
#define BFG_LIKELY(x) (__builtin_expect((x), 1))
#define BFG_UNLIKELY(x) (__builtin_expect((x), 0))
#else
#define BFG_LIKELY(x) (x)
#define BFG_UNLIKELY(x) (x)
#endif

#ifdef _MSC_VER
#define BFG_INLINE __forceinline
#elif defined(__clang__) || defined(__GNUC__)
#define BFG_INLINE inline __attribute__((__always_inline__))
#else
#define BFG_INLINE inline
#endif

#define BUFFER_SIZE 1048576

#define BFG_VERSION "1.0.6.4"

#define BFG_METABIN_FORMAT_VERSION 1
#define BFG_GRAPHBIN_FORMAT_VERSION 1
#define BFG_COLOREDCDBG_FORMAT_VERSION 3

#define BFG_METABIN_FORMAT_HEADER 0x267c3d5d
#define BFG_GRAPHBIN_FORMAT_HEADER 0x7e215f3f


static const char alpha[4] = {'A','C','G','T'};

BFG_INLINE bool isDNA(const char c) {

    static const size_t DNAbits[4] = {0x0ULL, static_cast<size_t>(0x10008A0010008AULL), 0x0ULL, 0x0ULL};

    return static_cast<bool>((DNAbits[c >> 6] >> (c & 0x3F)) & 0x1ULL);
}

BFG_INLINE size_t cstrMatch(const char* a, const char* b) {

    const char* a_ = a;

    while ((*a != '\0') && (*a == *b)){ ++a; ++b; }

    return a - a_;
}

BFG_INLINE char reverse_complement(const char nuc){

    const char c = nuc & 0xDF;

    if (isDNA(c)){

        const size_t x = (c & 4) >> 1;

        return alpha[3 - (x + ((x ^ (c & 2)) >> 1))];
    }
    else {

        switch(c) {

            case 'N': return 'N';
            case 'R': return 'Y';
            case 'Y': return 'R';
            case 'S': return 'S';
            case 'W': return 'W';
            case 'K': return 'M';
            case 'M': return 'K';
            case 'B': return 'V';
            case 'D': return 'H';
            case 'H': return 'D';
            case 'V': return 'B';
            case '.': return '.';
            default: return c;
        }
    }
}

BFG_INLINE std::string reverse_complement(const std::string& s){

    std::string seq(s);

    reverse(seq.begin(), seq.end());

    for (size_t i = 0; i < seq.length(); ++i){

        const char c = seq[i] & 0xDF;

        if (isDNA(c)){

            const size_t x = (c & 4) >> 1;

            seq[i] = alpha[3 - (x + ((x ^ (c & 2)) >> 1))];
        }
        else {

            switch(c) {

                case 'N': seq[i] = 'N'; break;
                case 'R': seq[i] = 'Y'; break;
                case 'Y': seq[i] = 'R'; break;
                case 'S': seq[i] = 'S'; break;
                case 'W': seq[i] = 'W'; break;
                case 'K': seq[i] = 'M'; break;
                case 'M': seq[i] = 'K'; break;
                case 'B': seq[i] = 'V'; break;
                case 'D': seq[i] = 'H'; break;
                case 'H': seq[i] = 'D'; break;
                case 'V': seq[i] = 'B'; break;
                case '.': seq[i] = '.'; break;
                default: seq[i] = c;
            }
        }
    }

    return seq;
}

BFG_INLINE std::string reverse_complement(const char* s){

    std::string seq(s);

    reverse(seq.begin(), seq.end());

    for (size_t i = 0; i < seq.length(); ++i){

        const char c = seq[i] & 0xDF;

        if (isDNA(c)){

            const size_t x = (c & 4) >> 1;

            seq[i] = alpha[3 - (x + ((x ^ (c & 2)) >> 1))];
        }
        else {

            switch(c) {

                case 'N': seq[i] = 'N'; break;
                case 'R': seq[i] = 'Y'; break;
                case 'Y': seq[i] = 'R'; break;
                case 'S': seq[i] = 'S'; break;
                case 'W': seq[i] = 'W'; break;
                case 'K': seq[i] = 'M'; break;
                case 'M': seq[i] = 'K'; break;
                case 'B': seq[i] = 'V'; break;
                case 'D': seq[i] = 'H'; break;
                case 'H': seq[i] = 'D'; break;
                case 'V': seq[i] = 'B'; break;
                case '.': seq[i] = '.'; break;
                default: seq[i] = c;
            }
        }
    }

    return seq;
}

BFG_INLINE size_t rndup(size_t v) {

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;

    return v;
}

BFG_INLINE uint32_t rndup(uint32_t v) {

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;

    return v;
}

BFG_INLINE uint16_t rndup(uint16_t v) {

    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v++;

    return v;
}

BFG_INLINE bool check_file_exists(const std::string& filename) {

    struct stat stFileInfo;

    return (stat(filename.c_str(), &stFileInfo) == 0);
}

BFG_INLINE uint32_t crc32_checksum(std::istream& in) {

    unsigned char buffer[65536];

    uint32_t crc = 0, prev_crc = crc + 1;

    while (in.good() && (crc != prev_crc)) {

        in.read(reinterpret_cast<char*>(buffer), 65536 * sizeof(char));

        const size_t len_buffer_read = in.gcount();

        prev_crc = crc;
        crc = crc32(crc, buffer, len_buffer_read);
    }

    return crc;
}

BFG_INLINE uint32_t crc32_checksum(const std::string& fn) {

    if ((fn.length() == 0) || !check_file_exists(fn)) return 0;

    std::ifstream infile;
    std::istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return crc32_checksum(in);
}

template<typename T> class wrapperData {

    public:

        BFG_INLINE const T* getData() const { return &data; }
        BFG_INLINE T* getData() { return &data; }

    private:

        T data;
};

template<> class wrapperData<void> {

    public:

        BFG_INLINE const void* getData() const { return nullptr; }
        BFG_INLINE void* getData() { return nullptr; }
};

#endif // BFG_COMMON_HPP
