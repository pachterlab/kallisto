#ifndef BIFROST_REPHASH_HPP
#define BIFROST_REPHASH_HPP

#include <cassert>
#include <stdint.h>

#include "wyhash.h"

static const unsigned char twin[32] = {
    0, 20, 2, 7, 4, 5, 6, 3,
    8,  9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 1, 21, 22, 23,
    24, 25, 26, 27, 28, 29, 30, 31
};

#if MAX_KMER_SIZE <= 64

static const uint64_t hvals[4] = {
    2053695854357871005ULL, 5073395517033431291ULL,
    10060236952204337488ULL, 7783083932390163561ULL
};

class RepHash {

    public:

        RepHash(const size_t _k = 0) {

            setK(_k);
        }

        void init(const char* _s) {

            h = 0;
            ht = 0;

            const unsigned char* s = (const unsigned char*) _s;

            for (size_t i = 0; i < k; ++i) {

                fastleftshift1(h);
                fastleftshift1(ht);

                h ^= hvals[charmask(s[i])];
                ht ^= hvals[twinmask(s[k-1-i])];
            }
        }

        inline void updateFW(const unsigned char out, const unsigned char in) {

            uint64_t z(hvals[charmask(out)]);
            uint64_t zt(hvals[twinmask(in)]);

            fastleftshiftk(z);
            fastleftshiftk(zt);

            fastleftshift1(h);

            h ^= z;
            h ^= hvals[charmask(in)];

            ht ^= zt;
            ht ^= hvals[twinmask(out)];

            fastrightshift1(ht);
        }

        inline void updateBW(const unsigned char out, const unsigned char in) {

            uint64_t z(hvals[twinmask(out)]);
            uint64_t zt(hvals[charmask(in)]);

            fastleftshiftk(z);
            fastleftshiftk(zt);

            fastleftshift1(ht);

            ht ^= z;
            ht ^= hvals[twinmask(in)];

            h ^= zt;
            h ^= hvals[charmask(out)];

            fastrightshift1(h);
        }

        inline void update(const unsigned char out, const unsigned char in){

            updateFW(out, in);
        }

        inline uint64_t hash() const {

            uint64_t hashes[2] = {h, ht};

            if (hashes[1] < hashes[0]) swap(hashes[0], hashes[1]);

            return wyhash(hashes, sizeof(uint64_t) + sizeof(uint64_t), 0, _wyp);
            //return (h ^ ht);
        }

        inline void setK(const size_t _k) {

            k = _k;
            h = 0;
            ht = 0;
        }

    private:

        inline uint64_t charmask (const unsigned char x) const {

            return (x & 6) >> 1;
        }

        inline uint64_t twinmask (const unsigned char x) const {

            return ((x ^ 4) & 6) >> 1;
        }

        inline void fastleftshiftk(uint64_t& x) const {

            x = (x << k) | (x >> (64-k));
        }

        inline void fastrightshiftk(uint64_t& x) const {

            x = (x >> k) | (x << (64-k));
        }

        inline void fastleftshift1(uint64_t& x) const {

            x = (x << 1) | (x >> 63);
        }

        inline void fastrightshift1(uint64_t& x) const {

            x = (x >> 1) | (x << 63);
        }

        size_t k;
        uint64_t h, ht;

        string str;
};

#else

struct rep_state_t {

    uint64_t hi;
    uint64_t lo;

    rep_state_t& operator^=(const rep_state_t& o) {

        hi ^= o.hi;
        lo ^= o.lo;

        return *this;
    }

    rep_state_t() : hi(0), lo(0) {}
    rep_state_t(const uint64_t hi_, const uint64_t lo_) : hi(hi_), lo(lo_) {}
};

static const rep_state_t hvals[32] = {
    rep_state_t(0x498bf4da68e4a5d2ULL, 0xcd18b2ed2719ae49ULL), rep_state_t(0x87613b9d792b27bfULL, 0x5a9f31e9650d3ac6ULL),
    rep_state_t(0x4395a1a049fa76a2ULL, 0x10e4b5bf779adbc8ULL), rep_state_t(0x7544bc518bb54d90ULL, 0x41f3e8de1fe48d6eULL),
    rep_state_t(0x8ed640367372ac6cULL, 0xa291bee592c5e1e7ULL), rep_state_t(0x4f1f5965708ec160ULL, 0xe71b6d1444ba099dULL),
    rep_state_t(0xc378038bae2f493bULL, 0x70d4008dad7e3aa5ULL), rep_state_t(0x615bf67b70de0cdfULL, 0x208f3adcc8ca1aa4ULL),
    rep_state_t(0xff239058a8bedba9ULL, 0x4e06dacb0b57d0d2ULL), rep_state_t(0x6221ac90a2ed270eULL, 0x1485df91c0b2a0efULL),
    rep_state_t(0x6aaa09ef06685b06ULL, 0xbddd68547f296a28ULL), rep_state_t(0xba49899319971eb2ULL, 0x105f2a98dd90b3f8ULL),
    rep_state_t(0x312a4ac3ddfdc014ULL, 0x412038b51ee959e5ULL), rep_state_t(0xc0cb77a29256270dULL, 0xbe6b7e9d4c6766deULL),
    rep_state_t(0xaa74044d51fdf18fULL, 0xa911bddbf8ffa08eULL), rep_state_t(0xcbdc765b04fa21d3ULL, 0xb18ad8b272439eecULL),
    rep_state_t(0x82765cedb75bd1c4ULL, 0x9445fa7b930d0503ULL), rep_state_t(0x1179cbc60b28e89cULL, 0x92104ea02ed0534fULL),
    rep_state_t(0xaee876b54b4df8daULL, 0xb3eaf5b4a10a08ccULL), rep_state_t(0xc34976b3ef5ee0a4ULL, 0x8210c76de6308193ULL),
    rep_state_t(0xee023d391e54a57bULL, 0x96e73dc6786b1f9aULL), rep_state_t(0x9929b98b3762d17fULL, 0xa1594e9916c0a563ULL),
    rep_state_t(0x395ecbebcfae31b1ULL, 0x79af8d5843ac918eULL), rep_state_t(0xd50a220a58e47a01ULL, 0x84d8e77ecba6959bULL),
    rep_state_t(0x4d710459ed613706ULL, 0x7b2bed1fd86d045eULL), rep_state_t(0xd671765f3784fd94ULL, 0x52d47f000a07bf45ULL),
    rep_state_t(0xfbb4a41a674fa47bULL, 0xaf81a75844adead1ULL), rep_state_t(0x55e243ef1134258fULL, 0x5de7beb874b7a97dULL),
    rep_state_t(0x85a7787202b8571cULL, 0x50f5fc4dd85517a9ULL), rep_state_t(0xb31bf1fcb1d83e09ULL, 0xa9cbec09c5de31b2ULL),
    rep_state_t(0x2a753d679e777b97ULL, 0x29641395bed797e9ULL), rep_state_t(0x1817b8f128494dbcULL, 0xb7fed4eae0bf9afcULL)
};

class RepHash {

    public:

        RepHash(const size_t _k = 0) {

            setK(_k);
        }

        inline void setK(const size_t _k) {

            full_k = _k;
            k = _k % 64;

            firstkmask = (1ULL << _k) - 1;
            lastkmask = firstkmask << (64 - k);
        }

        uint64_t hash() const {

            uint64_t hashes[2] = {h.lo, ht.lo};

            if (hashes[1] < hashes[0]) swap(hashes[0], hashes[1]);

            return wyhash(hashes, sizeof(uint64_t) + sizeof(uint64_t), 0, _wyp);
            //return (h.lo ^ ht.lo);
        }

        void init(const char *_s) {

            h = rep_state_t();
            ht = rep_state_t();

            const unsigned char *s = (const unsigned char*) _s;

            for (size_t i = 0; i < full_k; ++i) {

                fastleftshift1(h);
                fastleftshift1(ht);

                h ^= hvals[s[i] & charmask];
                ht ^= hvals[twin[s[full_k-1-i] & charmask]];
            }
        }

        inline void updateFW(const unsigned char out, const unsigned char in) {

            rep_state_t z(hvals[out & charmask]);
            rep_state_t zt(hvals[twin[in & charmask]]);

            fastleftshift1(h);
            fastleftshiftk(z);
            fastleftshiftk(zt);

            h ^= z;
            h ^= hvals[in & charmask];

            ht ^= hvals[twin[out & charmask]];
            ht ^= zt;

            fastrightshift1(ht);
        }

        inline void updateBW(const unsigned char out, const unsigned char in) {

            rep_state_t z(hvals[twin[out & charmask]]);
            rep_state_t zt(hvals[in & charmask]);

            fastleftshift1(ht);
            fastleftshiftk(z);
            fastleftshiftk(zt);

            ht ^= z;
            ht ^= hvals[twin[in & charmask]];

            h ^= hvals[out & charmask] ;
            h ^= zt;

            fastrightshift1(h);
        }

        inline void update(const unsigned char out, const unsigned char in){

            updateFW(out, in);
        }

    private:

        inline void fastleftshiftk(rep_state_t& x) const {

            const uint64_t upper = x.hi & lastkmask; // upper k of 64 bits of hi

            x.hi = (x.hi << k) | ((x.lo & lastkmask) >> (64 - k));
            x.lo = (x.lo << k) | (upper >> (64 - k));

            if (full_k & 64) swap(x.hi, x.lo);
        }

        inline void fastrightshiftk(rep_state_t& x) const {

            const uint64_t lower = x.hi & firstkmask; // lower k bits

            x.hi = (x.hi >> k) | ((x.lo & firstkmask) << (64 - k));
            x.lo = (x.lo >> k) | (lower << (64 - k));

            if (full_k & 64) swap(x.hi, x.lo);
        }

        inline void fastleftshift1(rep_state_t& x) const {

            uint64_t last1 = (x.hi & last1mask); // last bit of hi

            x.hi = (x.hi << 1) | ((x.lo & last1mask) >> 63);
            x.lo = (x.lo << 1) | (last1 >> 63);
        }

        inline void fastrightshift1(rep_state_t& x) const {

            const uint64_t first1 = (x.hi & 1ULL);  // first bit of hi

            x.hi = (x.hi >> 1) | ((x.lo & 1ULL) << 63);
            x.lo = (x.lo >> 1) | (first1 << 63);
        }

        size_t k, full_k;

        uint64_t lastkmask, firstkmask;

        static const uint64_t last1mask = 1ULL << 63;

        static const unsigned char charmask = 31;

        rep_state_t h, ht;
};

#endif

#endif
