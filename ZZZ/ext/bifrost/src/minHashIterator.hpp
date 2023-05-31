#ifndef BIFROST_MINHASHITERATOR_HPP
#define BIFROST_MINHASHITERATOR_HPP

#include <stdint.h>
#include <deque>
#include <vector>

#include <iostream>

#include "Kmer.hpp"

using namespace std;

struct minHashResult {

	minHashResult() : hash((uint64_t) -1),pos(-1) {}
	minHashResult(const uint64_t h, const int p) : hash(h), pos(p) {}
	minHashResult(const minHashResult& o) : hash(o.hash), pos(o.pos) {}

	uint64_t hash;
	int pos;
};

template<class HF>
struct minHashResultIterator;

template<class HF>
class minHashIterator {

    public:

        minHashIterator(const int _k, const int _g, const HF _h) :  s(nullptr), n(0), k(_k), g(_g), hf(_h),
                                                                    v(k-g+1), p(-1), invalid(true), nh(false) {

            hf.setK(g);
        }

        minHashIterator(const char* _s, const int _length, const int _k, const int _g, const HF _h, const bool _nh) :
                        k(_k), g(_g), hf(_h), v(k-g+1), p(-1), invalid(true), nh(_nh) {

            hf.setK(g);
            initString(_s, _length);
        }

        minHashIterator() : s(nullptr), n(0), k(0), g(0), hf(HF(0)), invalid(true), nh(false) {}

        minHashIterator(const minHashIterator& o) : s(o.s), n(o.n), k(o.k), g(o.g), hf(o.hf), v(o.v), p(o.p), invalid(o.invalid), nh(o.nh) {}

        inline void seed(const HF &ohf) {

            hf = ohf;
            invalid = true;
        }

        inline bool operator==(const minHashIterator& o) {

            if (invalid || o.invalid) return invalid && o.invalid;
            return (s == o.s) && (n == o.n) && (g == o.g) && (k == o.k) && (nh == o.nh);
        }

        inline bool operator!=(const minHashIterator& o) { return !operator==(o); }

        minHashIterator& operator++() {

            if (invalid) return *this;

            ++p; // advance to next k-mer

            if (p >= n-k+1) {

                invalid = true;
                return *this;
            }

            const int shift = static_cast<int>(nh);

            if (p==0) {

                hf.init(s + shift);

                v.push_back(minHashResult(hf.hash(), shift));

                for (int j = shift; j < k-g-shift;) {

                    hf.update(s[j], s[j+g]);

                    const uint64_t h = hf.hash();

                    int t = static_cast<int>(v.size()) - 1;

                    while ((t >= 0) && (v[t].hash > h)) {

                        v.pop_back();
                        --t;
                    }

                    ++j;
                    v.push_back(minHashResult(h, j));
                }
            }
            else {

                if (v[0].pos < p + shift) v.pop_front(); // remove first element, fell outside of window

                hf.update(s[p+k-g-1-shift],s[p+k-1-shift]);

                const uint64_t h = hf.hash();

                int t = static_cast<int>(v.size()) - 1;

                while ((t >= 0) && (v[t].hash > h)) {

                    v.pop_back();
                    --t;
                }

                v.push_back(minHashResult(h, p+k-g-shift));
            }

            return *this;
        }

        inline minHashIterator operator++(int) {

            minHashIterator tmp(*this);
            operator++();

            return tmp;
        }

        inline minHashIterator& operator+=(int i) {

            if (i >= k){

                const size_t inc = (i/k)*k;

                reinit(p + inc);
                i -= inc;
            }

            for (; i > 0; --i) operator++();

            return *this;
        }

        minHashResult getNewMin(const minHashResult& mhr_discard) const {

            if (invalid) return minHashResult();

            const int shift = static_cast<int>(nh);
            const int end = p+k-g-shift;

            int j = p + shift;

            HF hf_tmp;

            hf_tmp.setK(g);
            hf_tmp.init(&s[j]);

            while ((hf_tmp.hash() <= mhr_discard.hash) && (j < end)){

                hf_tmp.update(s[j], s[j+g]);
                ++j;
            }

            if ((j == end) && (hf_tmp.hash() <= mhr_discard.hash)) return mhr_discard; //No other minimizer can be found

            minHashResult mhr = minHashResult(hf_tmp.hash(), j);

            while (j < end) {

                hf_tmp.update(s[j], s[j+g]);
                ++j;

                const uint64_t h = hf_tmp.hash();

                if ((h <= mhr.hash) && (h > mhr_discard.hash)){

                    if (((h == mhr.hash) && (Minimizer(s + j).rep() < Minimizer(s + mhr.pos).rep())) || (h != mhr.hash)){

                        mhr.hash = h;
                        mhr.pos = j;
                    }
                }
            }

            return mhr;
        }

        BFG_INLINE minHashResultIterator<HF> operator*() const { return minHashResultIterator<HF>(this); }

        BFG_INLINE uint64_t getHash() const {

            return (invalid ? 0 : v[0].hash);
        }

        BFG_INLINE int getPosition() const {

            return (invalid ? 0 : v[0].pos);
        }

        BFG_INLINE int getKmerPosition() const {

            return p;
        }

        const char *s; //Minimizers are from k-mers, k-mers are from a sequence s
        int n; // Length of sequence s
        int k; // Length of k-mers
        int g; // Length of minimizers
        HF hf; // Rolling hash function
        deque<minHashResult> v; //Hash values and positions of a same minimizer with k-mer at position p
        int p; // Position of current k-mer traversed in the sequence
        bool invalid; // If sequence is invalid (iterating on k-mers out of bounds, etc.)
        bool nh; // If true, minimizer of k-mers km cannot start at position 0 or k-g

    private:

        void initString(const char* _s, const int _length) {

            n = _length;
            s = _s;
            p = -1;// reset position
            invalid = false;
            v.clear();

            // advance to first position or set to invalid
            if (n < k || k < g) invalid = true;
            else operator++();
        }

        void reinit(const size_t pos){

            if (invalid) return;

            p = pos; // advance to next k-mer

            if (p >= n-k+1){

                invalid = true;
                return;
            }

            const int shift = static_cast<int>(nh);

            hf.init(&s[p + shift]);

            v.clear();
            v.push_back(minHashResult(hf.hash(), p + shift));

            for (int j = p+shift; j < p+k-g-shift;) {

                hf.update(s[j], s[j+g]);

                const uint64_t h = hf.hash();

                int t = static_cast<int>(v.size()) - 1;

                while ((t >= 0) && (v[t].hash > h)) {

                    v.pop_back();
                    --t;
                }

                ++j;
                v.push_back(minHashResult(h, j));
            }
        }
};

template <class HF>
struct minHashResultIterator {

    minHashResultIterator(const minHashIterator<HF> *p) : p(p), invalid(false), pos(0), p_pos(p->p), p_s(p->s) {}

    minHashResultIterator() : invalid(true), p_pos(-1), p_s(nullptr) {}

    inline bool operator==(const minHashResultIterator& o) {

        if (o.invalid || invalid) return o.invalid && invalid;
        return (p_pos == o.p_pos) && (p_s == o.p_s) && (pos == o.pos);
    }

    inline bool operator!=(const minHashResultIterator& o) { return !operator==(o); }

    minHashResultIterator& operator++() {

        if (!invalid){

            if ((p_pos != p->p || p_s != p->s) // check if parent iterator has moved
                || (pos>=p->v.size()-1) // or if we advance past the end
                || (p->v[pos+1].hash != p->v[pos].hash))    { // or if the next position doesn't share the hash value
                // invalidate the iterator
                invalid = true;
            }
            else ++pos; // advance to next equal hash value
        }

        return *this;
    }


    inline minHashResultIterator operator++(int) {

        minHashResultIterator tmp(*this);
        operator++();

        return tmp;
    }

    BFG_INLINE const minHashResult& operator*() const { return p->v[pos]; }
    BFG_INLINE const minHashResult* operator->() const { return &(p->v[pos]); }

    // pos points to a minHashIterator, all the values from p.v[0] to p.v[pos] have the
    // same (minimum) hash value or the this.invalid is true
    // at the time this was created p.s==p_s and p.p=p_pos
    const minHashIterator<HF>* p;
    bool invalid;
    int pos;
    int p_pos;
    const char* p_s;
};


template<class HF>
struct preAllocMinHashResultIterator;

template<class HF>
class preAllocMinHashIterator {

    public:

        preAllocMinHashIterator() : s(nullptr), n(0), k(0), g(0), hf(HF(0)), p(-1), v(0),
                                    p_cur_start(0), p_cur_end(0), invalid(true), nh(false) {}

        preAllocMinHashIterator(const char* _s, const int _n, const int _k, const int _g, const HF _h, const bool _nh) :
                                s(_s), n(_n), k(_k), g(_g), hf(_h), p(-1), p_cur_start(0), p_cur_end(0), invalid(true), nh(_nh) {

            if ((s != nullptr) && (n >= k) && (k >= g)){

                invalid = false;

                v = vector<minHashResult>(n - g + 1);
                hf.setK(g);

                operator++();
            }
        }

        inline bool operator==(const preAllocMinHashIterator& o) {

            if (invalid || o.invalid) return invalid && o.invalid;
            return (s == o.s) && (n == o.n) && (g == o.g) && (k == o.k) && (nh == o.nh);
        }

        inline bool operator!=(const preAllocMinHashIterator& o) { return !this->operator==(o); }

        // invariant:
        //   v[0..sz] contains only ascending elements in s from [p,p+k-g+1)
        //
        //   there exists no a<b s.t. v[a] > v[b].
        preAllocMinHashIterator& operator++() {

            if (invalid) return *this;

            ++p; // advance to next k-mer

            if (p >= n-k+1) {
                // out of bounds
                invalid = true;
                return *this;
            }

            const int shift = static_cast<int>(nh);

            if (p == 0) {

                hf.init(&s[shift]);

                v[p_cur_end] = minHashResult(hf.hash(), shift);
                ++p_cur_end;

                for (int j = shift; j < k-g-shift;) {

                    hf.update(s[j], s[j+g]);

                    const uint64_t h = hf.hash();

                    while ((p_cur_end > p_cur_start) && (v[p_cur_end-1].hash > h)) --p_cur_end;

                    ++j;

                    v[p_cur_end] = minHashResult(h,j);

                    ++p_cur_end;
                }
            }
            else {

                p_cur_start += static_cast<size_t>(v[p_cur_start].pos < (p + shift));

                hf.update(s[p+k-g-1-shift], s[p+k-1-shift]);

                const uint64_t h = hf.hash();

                while ((p_cur_end > p_cur_start) && (v[p_cur_end-1].hash > h)) --p_cur_end;

                v[p_cur_end] = minHashResult(h, p+k-g-shift);

                ++p_cur_end;
            }

            return *this;
        }

        inline preAllocMinHashIterator operator++(int) {

            preAllocMinHashIterator tmp(*this);
            operator++();

            return tmp;
        }

        inline preAllocMinHashIterator& operator+=(int i) {

            for (; i > 0; --i) operator++();
            return *this;
        }

        BFG_INLINE preAllocMinHashResultIterator<HF> operator*() const { return preAllocMinHashResultIterator<HF>(*this); }

        BFG_INLINE uint64_t getHash() const { return ((static_cast<uint64_t>(invalid) - 1) & v[p_cur_start].hash); }

        BFG_INLINE int getPosition() const { return invalid ? 0 : v[p_cur_start].pos; }

        BFG_INLINE int getNbMin() const { return p_cur_end - p_cur_start; }

        BFG_INLINE int getKmerPosition() const { return p; }

        minHashResult getNewMin(const minHashResult& mhr_discard) const {

            if (invalid) return minHashResult();

            const int shift = static_cast<int>(nh);
            const int end = p+k-g-shift;

            int j = p + shift;

            HF hf_tmp;

            hf_tmp.setK(g);
            hf_tmp.init(&s[j]);

            while ((hf_tmp.hash() <= mhr_discard.hash) && (j < end)){

                hf_tmp.update(s[j], s[j+g]);
                ++j;
            }

            if ((j == end) && (hf_tmp.hash() <= mhr_discard.hash)) return mhr_discard;

            minHashResult mhr = minHashResult(hf_tmp.hash(), j);

            while (j < end) {

                hf_tmp.update(s[j], s[j+g]);
                ++j;

                const uint64_t h = hf_tmp.hash();

                if ((h <= mhr.hash) && (h > mhr_discard.hash)){

                    if (((h == mhr.hash) && (Minimizer(s + j).rep() < Minimizer(s + mhr.pos).rep())) || (h != mhr.hash)) {

                        mhr.hash = h;
                        mhr.pos = j;
                    }
                }
            }

            return mhr;
        }

        const char* s; //Minimizers are from k-mers, k-mers are from a sequence s
        int n; // Length of sequence s
        int k; // Length of k-mers
        int g; // Length of minimizers
        HF hf; // Rolling hash function
        vector<minHashResult> v; //Hash values and positions of a same minimizer with k-mer at position p
        size_t p_cur_start;
        size_t p_cur_end;
        int p; // Position of current k-mer traversed in the sequence
        bool invalid; // If sequence is invalid (iterating on k-mers out of bounds, etc.)
        bool nh;

        // private copy constructor
        preAllocMinHashIterator(const preAllocMinHashIterator& o) : s(o.s), n(o.n), k(o.k), g(o.g), hf(o.hf), v(o.v), p(o.p),
                                                                    p_cur_start(o.p_cur_start), p_cur_end(o.p_cur_end), invalid(o.invalid), nh(o.nh) {}

        preAllocMinHashIterator(const preAllocMinHashIterator& o, int len) : s(o.s + o.p), n(len), k(o.k), g(o.g), hf(o.hf), p(0), p_cur_start(0),
                                                                                p_cur_end(o.p_cur_end - o.p_cur_start), invalid(o.invalid), nh(o.nh) {

            if (!invalid && (o.p + n <= o.n)){

                vector<minHashResult> v_tmp(o.v.begin() + o.p_cur_start, o.v.begin() + o.p_cur_end);

                v = move(v_tmp);

                for (auto& min_h : v) min_h.pos -= o.p;
            }
            else invalid = true;
        }
};

template <class HF>
struct preAllocMinHashResultIterator {

    preAllocMinHashResultIterator(const preAllocMinHashIterator<HF>& _p) :  p(&_p), p_pos(_p.p), p_it(_p.p_cur_start),
                                                                            p_it_end(_p.p_cur_end), p_s(_p.s), invalid(false) {}

    preAllocMinHashResultIterator() : invalid(true), p(nullptr), p_pos(-1), p_s(nullptr), p_it(0), p_it_end(0) {}

    inline bool operator==(const preAllocMinHashResultIterator& o) const {

        if (o.invalid || invalid) return (o.invalid && invalid);
        return (p_pos == o.p_pos) && (p_s == o.p_s )&& (p_it == o.p_it) && (p_it_end == o.p_it_end);
    }

    inline bool operator!=(const preAllocMinHashResultIterator& o) const { return !operator==(o); }

    preAllocMinHashResultIterator& operator++() {

        if (!invalid){

            if ((p_s != p->s || p_pos != p->p || p_it_end != p->p_cur_end) // check if parent iterator has moved
                || (p_it >= p_it_end - 1) // or if we advance past the end
                || (p->v[p_it + 1].hash != p->v[p_it].hash)) { // or if the next position doesn't share the hash value
                // invalidate the iterator
                invalid = true;
            }
            else ++p_it; // advance to next equal hash value
        }

        return *this;
    }

    inline preAllocMinHashResultIterator operator++(int) {

        preAllocMinHashResultIterator tmp(*this);
        operator++();

        return tmp;
    }

    preAllocMinHashResultIterator& operator=(const preAllocMinHashResultIterator &o){

        p_it = o.p_it;
        p_it_end = o.p_it_end;
        invalid = o.invalid;

        invalid = operator!=(o);

        return *this;
    }

    BFG_INLINE const minHashResult& operator*() const { return p->v[p_it]; }
    BFG_INLINE const minHashResult* operator->() const { return &(p->v[p_it]); }

    // pos points to a minHashIterator, all the values from p.v[0] to p.v[pos] have the
    // same (minimum) hash value or the this.invalid is true
    // at the time this was created p.s==p_s and p.p=p_pos

    const preAllocMinHashIterator<HF>* p;
    bool invalid;
    size_t p_it;
    size_t p_it_end;
    int p_pos;
    const char* p_s;
};

template <class HF>
struct minHashKmer {

    public:

        minHashKmer(const char* _s, const int _k, const int _g, const HF _h, const bool neighbor_hash) :
                    s(_s), k(_k), g(_g), hf(_h), h(0), i(0), p(0), invalid(true), nh(neighbor_hash) {

            if ((s != nullptr) && ((n = strlen(s)) >= k) && (k >= g) && (k <= MAX_KMER_SIZE)){

                invalid = false;

                compute_min();
            }
        }

        minHashKmer() : s(nullptr), n(0), k(0), g(0), hf(HF(0)), h(0), i(0), p(0), invalid(true), nh(false) {}

        minHashKmer& operator++() {

            ++i;
            invalid = invalid || (i >= p);

            return *this;
        }

        BFG_INLINE minHashKmer operator++(int) {

            minHashKmer tmp(*this);
            operator++();

            return tmp;
        }

        minHashKmer& operator=(const minHashKmer &o){

            s = o.s;
            n = o.n;
            k = o.k;
            g = o.g;
            h = o.h;
            i = o.i;
            p = o.p;
            hf = o.hf;
            nh = o.nh;
            invalid = o.invalid;

            memcpy(pos, o.pos, p * sizeof(uint16_t));

            return *this;
        }

        BFG_INLINE bool operator==(const minHashKmer& o) {

            if (invalid || o.invalid) return (invalid && o.invalid);

            return  (s == o.s) && (n == o.n) && (g == o.g) && (k == o.k) && (nh == o.nh) &&
                    (p == o.p) && (memcmp(pos, o.pos, p * sizeof(uint16_t)) == 0);
        }

        BFG_INLINE bool operator!=(const minHashKmer& o) { return !operator==(o); }

        BFG_INLINE uint64_t getHash() const { return h; }
        BFG_INLINE int getPosition() const { return pos[i]; }

        BFG_INLINE void getNewMin() {

            const uint64_t prev_h = h;

            h = 0;
            i = 0;
            p = -1;

            compute_min(prev_h);
        }

    private:

        void compute_min(){

            if (invalid) return;

            const int shift = static_cast<int>(nh);

            hf.setK(g);
            hf.init(&s[shift]);

            p = 1;
            i = 0;
            h = hf.hash();
            pos[0] = shift;

            for (int j = shift; j < k-g-shift; ++j) {

                hf.update(s[j], s[j+g]);

                const uint64_t h_v = hf.hash();

                if (h_v < h){

                    h = h_v;
                    p = 1;
                    pos[0] = j + 1;
                }
                else if (h_v == h){

                    pos[p] = j + 1;
                    ++p;
                }
            }
        }

        void compute_min(const uint64_t min_v){

            if (invalid) return;

            const int shift = static_cast<int>(nh);

            hf.setK(g);
            hf.init(&s[shift]);

            i = 0;
            p = 0;

            if (hf.hash() > min_v){

                h = hf.hash();
                p = 1;
                pos[0] = shift;
            }

            for (int j = shift; j < k-g-shift; ++j) {

                hf.update(s[j], s[j+g]);

                const uint64_t h_v = hf.hash();

                if (h_v > min_v){

                    if ((p == 0) || (h_v < h) || ((h_v == h) && (Minimizer(s + j + 1).rep() < Minimizer(s + pos[0]).rep()))){

                        h = h_v;
                        p = 1;
                        pos[0] = j + 1;
                    }
                }
            }

            invalid = (p == 0);
        }

        const char* s;
        HF hf;
        uint64_t h;
        int n;
        int k;
        int g;
        int p;
        int i;
        uint16_t pos[MAX_KMER_SIZE];
        bool invalid;
        bool nh;
};

#endif // MINHASHITERATOR_H
