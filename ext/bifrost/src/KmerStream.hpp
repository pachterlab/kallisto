#ifndef BIFROST_KMERSTREAM_HPP
#define BIFROST_KMERSTREAM_HPP

#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <getopt.h>
#include <sys/stat.h>
#include <fstream>
#include <stdint.h>
#include <sstream>
#include <bitset>

#include <thread>
#include <atomic>

#include "minHashIterator.hpp"
#include "wyhash.h"

#include "File_Parser.hpp"
#include "Lock.hpp"
#include "RepHash.hpp"
#include "StreamCounter.hpp"

#define NB_STREAMCOUNTER_PER_READHASHER 1

using namespace std;

struct KmerStream_Build_opt {

    vector<string> files;

    bool verbose;

    double e;

    size_t k;
    size_t g;
    size_t q;
    size_t q_base;
    size_t threads;
    size_t chunksize;

    KmerStream_Build_opt() : q_base(33), q(0), k(31), g(23), verbose(false), e(0.01), threads(1), chunksize(64) {}
};

class ReadHasherMinimizer;

class ReadQualityHasherMinimizer {

    friend class ReadHasherMinimizer;

    public:

        ReadQualityHasherMinimizer() :    k(0), g(0), q_cutoff(0), q_base(0) {}

        /*ReadQualityHasherMinimizer(const double e_, const size_t q_base_) :    k(0), g(0), q_cutoff(0), q_base(q_base_), sc_min(e_) {

            for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc_km[i].initialize(e_);
        }

        ReadQualityHasherMinimizer(const ReadQualityHasherMinimizer& o) :   k(o.k), g(o.g), q_cutoff(o.q_cutoff), q_base(o.q_base), sc_min(o.sc_min) {

            for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc_km[i] = o.sc_km[i];
        }*/

        ReadQualityHasherMinimizer(const double e_, const size_t q_base_) :    k(0), g(0), q_cutoff(0), q_base(q_base_), sc_min(e_), sc_km(e_) {}
        ReadQualityHasherMinimizer(const ReadQualityHasherMinimizer& o) :   k(o.k), g(o.g), q_cutoff(o.q_cutoff), q_base(o.q_base), sc_min(o.sc_min), sc_km(o.sc_km) {}

        BFG_INLINE void setK(const size_t _k) {

            k = _k;
        }

        BFG_INLINE void setG(const size_t _g) {

            g = _g;
        }

        BFG_INLINE void setQualityCutoff(const size_t q) {

            q_cutoff = q;
        }

        void clear() {

            k = 0;
            g = 0;
            q_cutoff = 0;
            q_base = 0;

            sc_min.clear();

            /*for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i)*/ sc_km/*[i]*/.clear();
        }

        void initialize(const double e_, const size_t q_cutoff_, const size_t q_base_, const size_t _k, const size_t _g) {

            clear();

            k = _k;
            g = _g;

            q_cutoff = q_cutoff_;
            q_base = q_base_;

            sc_min.initialize(e_);

            /*for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i)*/ sc_km/*[i]*/.initialize(e_);
        }

        // create hashes for all k-mers
        // operate on hashes
        pair<size_t, size_t> update(const char* s, const size_t l, const char* q, const size_t ql) {

            if (l < k) return {0, 0};

            size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;
            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            bool last_valid = false;

            const char q_base_cut = static_cast<char>(q_base + q_cutoff);

            minHashIterator<RepHash> min_it = minHashIterator<RepHash>(s, l, k, g, RepHash(), true);

            RepHash hf;

            hf.setK(k);

            while (j < l) {
                // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed
                if (isDNA(s[j]) && (q[j] >= q_base_cut)) {

                    if (last_valid) { // s[i..j-1] was a valid k-mer k-mer, update

                        hf.update(s[i],s[j]);

                        ++i;
                        ++min_it;
                    }
                    else if (i + k - 1 == j) {

                        hf.init(s+i); // start the k-mer at position i

                        last_valid = true;
                        min_it += i - min_it.getKmerPosition();
                    }

                    ++j;
                }
                else { // invalid character, restart

                    ++j;

                    i = j;
                    last_valid = false;
                }

                if (last_valid){

                    const size_t min_pos = min_it.getPosition();

                    sc_km/*[0]*/.update(hf.hash());

                    ++l_sumCount_km_F1;

                    if (min_pos != prev_pos_min){

                        sc_min.update(min_it.getHash());

                        ++l_sumCount_min_F1;

                        prev_pos_min = min_pos;
                    }
                }
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        pair<size_t, size_t> update(const char* seq_buf, const char* qual_buf, const size_t buf_sz) {

            const char q_base_cut = static_cast<char>(q_base + q_cutoff);

            const char* str = seq_buf;
            const char* str_end = seq_buf + buf_sz;
            const char* q_str = qual_buf;

            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            while (str < str_end) { // for each input

                const int sl = strlen(str);

                if (sl >= k){

                    size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;

                    bool last_valid = false;

                    minHashIterator<RepHash> min_it = minHashIterator<RepHash>(str, sl, k, g, RepHash(), true);

                    RepHash hf;

                    hf.setK(k);

                    while (j < sl) {
                        // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed
                        if (isDNA(str[j]) && (q_str[j] >= q_base_cut)) {

                            if (last_valid) {
                                // s[i..j-1] was a valid k-mer k-mer, update
                                hf.update(str[i], str[j]);

                                ++i;
                                ++min_it;
                            }
                            else if (i + k - 1 == j) {

                                hf.init(str + i); // start the k-mer at position i

                                last_valid = true;
                                min_it += i - min_it.getKmerPosition();
                            }

                            ++j;
                        }
                        else { // invalid character, restart

                            ++j;
                            i = j;
                            last_valid = false;
                        }

                        if (last_valid){

                            const size_t min_pos = min_it.getPosition();

                            sc_km/*[0]*/.update(hf.hash());

                            ++l_sumCount_km_F1;

                            if (min_pos != prev_pos_min){

                                sc_min.update(min_it.getHash());

                                ++l_sumCount_min_F1;

                                prev_pos_min = min_pos;
                            }
                        }
                    }
                }

                str += sl + 1;
                q_str += sl + 1;
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        pair<size_t, size_t> update_p(const char* s, const size_t l, const char* q, const size_t ql) {

            if (l < k) return {0,0};

            size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;
            size_t min_idx = 0;
            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            bool last_valid = false;

            const char q_base_cut = static_cast<char>(q_base + q_cutoff);

            minHashIterator<RepHash> min_it = minHashIterator<RepHash>(s, l, k, g, RepHash(), true);

            RepHash hf;

            hf.setK(k);

            while (j < l) {
                // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed
                if (isDNA(s[j]) && (q[j] >= q_base_cut)) {

                    if (last_valid) { // s[i..j-1] was a valid k-mer k-mer, update

                        hf.update(s[i],s[j]);

                        ++i;
                        ++min_it;
                    }
                    else if (i + k - 1 == j) {

                        hf.init(s+i); // start the k-mer at position i

                        last_valid = true;
                        min_it += i - min_it.getKmerPosition();
                    }

                    ++j;
                }
                else { // invalid character, restart

                    ++j;
                    i = j;
                    last_valid = false;
                }

                if (last_valid) {

                    const size_t min_pos = min_it.getPosition();

                    /*if (min_pos != prev_pos_min){

                        const size_t min_h = min_it.getHash();

                        prev_pos_min = min_pos;
                        min_idx = min_h % NB_STREAMCOUNTER_PER_READHASHER;

                        sc_min.update_p(min_h);

                        ++l_sumCount_min_F1;
                    }
                    
                    sc_km[min_idx].update_p(hf.hash());*/

                    sc_km.update_p(hf.hash());

                    if (min_pos != prev_pos_min){

                        sc_min.update_p(min_it.getHash());

                        prev_pos_min = min_pos;

                        ++l_sumCount_min_F1;
                    }

                    ++l_sumCount_km_F1;
                }
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        pair<size_t, size_t> update_p(const char* seq_buf, const char* qual_buf, const size_t buf_sz) {

            const char q_base_cut = static_cast<char>(q_base + q_cutoff);

            const char* str = seq_buf;
            const char* str_end = &seq_buf[buf_sz];
            const char* q_str = qual_buf;

            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            while (str < str_end) { // for each input

                const int sl = strlen(str);

                if (sl >= k){

                    size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;
                    size_t min_idx = 0;

                    bool last_valid = false;

                    minHashIterator<RepHash> min_it = minHashIterator<RepHash>(str, sl, k, g, RepHash(), true);

                    RepHash hf;

                    hf.setK(k);

                    while (j < sl) {
                        // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed
                        if (isDNA(str[j]) && (q_str[j] >= q_base_cut)) {

                            if (last_valid) {
                                // s[i..j-1] was a valid k-mer k-mer, update
                                hf.update(str[i], str[j]);

                                ++i;
                                ++min_it;
                            }
                            else if (i + k - 1 == j) {

                                hf.init(str + i); // start the k-mer at position i

                                last_valid = true;
                                min_it += i - min_it.getKmerPosition();
                            }

                            ++j;
                        }
                        else { // invalid character, restart

                            ++j;
                            i = j;
                            last_valid = false;
                        }

                        if (last_valid) {

                            const size_t min_pos = min_it.getPosition();

                            /*if (min_pos != prev_pos_min){

                                const size_t min_h = min_it.getHash();

                                prev_pos_min = min_pos;
                                min_idx = min_h % NB_STREAMCOUNTER_PER_READHASHER;

                                sc_min.update_p(min_h);

                                ++l_sumCount_min_F1;
                            }
                            
                            sc_km[min_idx].update_p(hf.hash());*/

                            sc_km.update_p(hf.hash());

                            if (min_pos != prev_pos_min){

                                sc_min.update_p(min_it.getHash());

                                prev_pos_min = min_pos;

                                ++l_sumCount_min_F1;
                            }



                            ++l_sumCount_km_F1;
                        }
                    }
                }

                str += sl + 1;
                q_str += sl + 1;
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        BFG_INLINE bool join(const ReadQualityHasherMinimizer& o) {

            bool join_rhm = sc_min.join(o.sc_min);

            /*for (size_t i = 0; join_rhm && (i < NB_STREAMCOUNTER_PER_READHASHER); ++i)*/ join_rhm = join_rhm && (sc_km/*[i]*/.join(o.sc_km/*[i]*/));

            return join_rhm;
        }

        BFG_INLINE bool join(const ReadHasherMinimizer& o);

        /*BFG_INLINE size_t KmerF0() const {

            StreamCounter sc(sc_km[0]);

            for (size_t i = 1; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc.join(sc_km[i]);

            return sc.F0();
        }

        BFG_INLINE size_t Kmerf1() const {

            StreamCounter sc(sc_km[0]);

            for (size_t i = 1; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc.join(sc_km[i]);

            return sc.f1();
        }*/

        BFG_INLINE size_t KmerF0() const { return sc_km.F0(); }

        BFG_INLINE size_t Kmerf1() const { return sc_km.f1(); }

        BFG_INLINE size_t MinimizerF0() const {

            return sc_min.F0();
        }

        BFG_INLINE size_t Minimizerf1() const {

            return sc_min.f1();
        }

    private:

        size_t q_cutoff;
        size_t q_base;

        size_t k;
        size_t g;

        StreamCounter sc_min;
        StreamCounter sc_km/*[NB_STREAMCOUNTER_PER_READHASHER]*/;
};

class ReadHasherMinimizer {

    friend class ReadQualityHasherMinimizer;

    public:

        ReadHasherMinimizer() :  k(0), g(0) {}

        /*ReadHasherMinimizer(const double e) : k(0), g(0), sc_min(e) {

            for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc_km[i].initialize(e);
        }

        ReadHasherMinimizer(const ReadHasherMinimizer& o) : k(o.k), g(o.g), sc_min(o.sc_min) {

            for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc_km[i] = o.sc_km[i];
        }

        ReadHasherMinimizer(ReadQualityHasherMinimizer&& o) : k(o.k), g(o.g), sc_min(move(o.sc_min)) {

            for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc_km[i] = move(o.sc_km[i]);
        }

        ReadHasherMinimizer(const ReadQualityHasherMinimizer& o) : k(o.k), g(o.g), sc_min(o.sc_min) {

            for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc_km[i] = o.sc_km[i];
        }*/

        ReadHasherMinimizer(const double e) : k(0), g(0), sc_min(e), sc_km(e) {}

        ReadHasherMinimizer(const ReadHasherMinimizer& o) : k(o.k), g(o.g), sc_min(o.sc_min), sc_km(o.sc_km) {}

        ReadHasherMinimizer(ReadQualityHasherMinimizer&& o) : k(o.k), g(o.g), sc_min(move(o.sc_min)), sc_km(move(o.sc_km)) {}

        ReadHasherMinimizer(const ReadQualityHasherMinimizer& o) : k(o.k), g(o.g), sc_min(o.sc_min), sc_km(o.sc_km) {}

        void clear() {

            k = 0;
            g = 0;

            sc_min.clear();

            /*for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i)*/ sc_km/*[i]*/.clear();
        }

        void initialize(const double e_, const size_t _k, const size_t _g) {

            clear();

            k = _k;
            g = _g;

            sc_min.initialize(e_);

            /*for (size_t i = 0; i < NB_STREAMCOUNTER_PER_READHASHER; ++i)*/ sc_km/*[i]*/.initialize(e_);
        }

        BFG_INLINE void setK(const size_t _k) {

            k = _k;
        }

        BFG_INLINE void setG(const size_t _g) {

            g = _g;
        }

        // create hashes for all k-mers
        // operate on hashes
        pair<size_t, size_t> update(const char *s, const size_t l, const char *q, const size_t ql) {

            if (l < k) return {0,0};

            size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;
            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            bool last_valid = false;

            minHashIterator<RepHash> min_it = minHashIterator<RepHash>(s, l, k, g, RepHash(), true);

            RepHash hf;

            hf.setK(k);

            while (j < l) { // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed

                if (isDNA(s[j])) {

                    if (last_valid) { // s[i..j-1] was a valid k-mer k-mer, update

                        hf.update(s[i],s[j]);

                        ++i;
                        ++min_it;
                    }
                    else if (i + k - 1 == j) {

                        hf.init(s+i); // start the k-mer at position i

                        last_valid = true;
                        min_it += i - min_it.getKmerPosition();
                    }

                    ++j;
                }
                else { // invalid character, restart

                    ++j;
                    i = j;

                    last_valid = false;
                }

                if (last_valid){

                    const size_t min_pos = min_it.getPosition();

                    sc_km/*[0]*/.update(hf.hash());

                    ++l_sumCount_km_F1;

                    if (min_pos != prev_pos_min){

                        sc_min.update(min_it.getHash());

                        ++l_sumCount_min_F1;

                        prev_pos_min = min_pos;
                    }
                }
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        pair<size_t, size_t> update(const char* seq_buf, const size_t seq_buf_sz) {

            const char* str = seq_buf;
            const char* str_end = seq_buf + seq_buf_sz;

            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            while (str < str_end) { // for each input

                const int sl = strlen(str);

                if (sl >= k){

                    size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;

                    bool last_valid = false;

                    minHashIterator<RepHash> min_it = minHashIterator<RepHash>(str, sl, k, g, RepHash(), true);

                    RepHash hf;

                    hf.setK(k);

                    while (j < sl) { // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed

                        if (isDNA(str[j])) {

                            if (last_valid) {
                                // s[i..j-1] was a valid k-mer k-mer, update
                                hf.update(str[i], str[j]);

                                ++i;
                                ++min_it;
                            }
                            else if (i + k -1 == j) {

                                hf.init(str + i); // start the k-mer at position i

                                last_valid = true;
                                min_it += i - min_it.getKmerPosition();
                            }

                            ++j;
                        }
                        else { // invalid character, restart

                            ++j;
                            i = j;
                            last_valid = false;
                        }

                        if (last_valid){

                            const size_t min_pos = min_it.getPosition();

                            sc_km/*[0]*/.update(hf.hash());

                            ++l_sumCount_km_F1;

                            if (min_pos != prev_pos_min){

                                sc_min.update(min_it.getHash());

                                ++l_sumCount_min_F1;

                                prev_pos_min = min_pos;
                            }
                        }
                    }
                }

                str += sl + 1;
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        pair<size_t, size_t> update_p(const char *s, const size_t l, const char *q, const size_t ql) {

            if (l < k) return {0,0};

            size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;
            size_t min_idx = 0;

            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            bool last_valid = false;

            minHashIterator<RepHash> min_it = minHashIterator<RepHash>(s, l, k, g, RepHash(), true);

            RepHash hf;

            hf.setK(k);

            while (j < l) { // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed

                if (isDNA(s[j])) {

                    if (last_valid) { // s[i..j-1] was a valid k-mer k-mer, update

                        hf.update(s[i],s[j]);

                        ++i;
                        ++min_it;
                    }
                    else if (i + k - 1 == j) {

                        hf.init(s+i); // start the k-mer at position i

                        last_valid = true;
                        min_it += i - min_it.getKmerPosition();
                    }

                    ++j;
                }
                else { // invalid character, restart

                    ++j;
                    i = j;

                    last_valid = false;
                }

                if (last_valid) {

                    /*const size_t min_pos = min_it.getPosition();

                    if (min_pos != prev_pos_min){

                        const size_t min_h = min_it.getHash();

                        prev_pos_min = min_pos;
                        min_idx = min_h % NB_STREAMCOUNTER_PER_READHASHER;

                        ++l_sumCount_min_F1;

                        sc_min.update_p(min_h);
                    }

                    ++l_sumCount_km_F1;
                    
                    sc_km[min_idx].update_p(hf.hash());*/

                    sc_km.update_p(hf.hash());

                    ++l_sumCount_km_F1;

                    if (min_it.getPosition() != prev_pos_min){

                        sc_min.update_p(min_it.getHash());

                        prev_pos_min = min_it.getPosition();

                        ++l_sumCount_min_F1;
                    }
                }
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        pair<size_t, size_t> update_p(const char* seq_buf, const size_t seq_buf_sz) {

            const char* str = seq_buf;
            const char* str_end = seq_buf + seq_buf_sz;

            size_t l_sumCount_min_F1 = 0;
            size_t l_sumCount_km_F1 = 0;

            while (str < str_end) { // for each input

                const int sl = strlen(str);

                if (sl >= k){

                    size_t i = 0, j = 0, prev_pos_min = 0xffffffffffffffffULL;
                    size_t min_pos = 0, min_h = 0, min_idx = 0;

                    bool last_valid = false;

                    minHashIterator<RepHash> min_it = minHashIterator<RepHash>(str, sl, k, g, RepHash(), true);

                    RepHash hf;

                    hf.setK(k);

                    while (j < sl) { // s[i...j-1] is a valid string, all k-mers in s[..j-1] have been processed

                        if (isDNA(str[j])) {

                            if (last_valid) {
                                // s[i..j-1] was a valid k-mer k-mer, update
                                hf.update(str[i], str[j]);

                                ++i;
                                ++min_it;
                            }
                            else if (i + k -1 == j) {

                                hf.init(str + i); // start the k-mer at position i

                                last_valid = true;
                                min_it += i - min_it.getKmerPosition();
                            }

                            ++j;
                        }
                        else { // invalid character, restart

                            ++j;

                            i = j;
                            last_valid = false;
                        }

                        if (last_valid) {

                            /*const size_t min_pos = min_it.getPosition();

                            if (min_pos != prev_pos_min){

                                const size_t min_h = min_it.getHash();

                                prev_pos_min = min_pos;
                                min_idx = min_h % NB_STREAMCOUNTER_PER_READHASHER;

                                ++l_sumCount_min_F1;

                                sc_min.update_p(min_h);
                            }

                            ++l_sumCount_km_F1;
                            
                            sc_km[min_idx].update_p(hf.hash());*/

                            sc_km.update_p(hf.hash());

                            ++l_sumCount_km_F1;

                            if (min_it.getPosition() != prev_pos_min){

                                sc_min.update_p(min_it.getHash());

                                prev_pos_min = min_it.getPosition();

                                ++l_sumCount_min_F1;
                            }
                        }
                    }
                }

                str += sl + 1;
            }

            return {l_sumCount_min_F1, l_sumCount_km_F1};
        }

        BFG_INLINE void setQualityCutoff(const size_t q) {}

        /*BFG_INLINE bool join(const ReadHasherMinimizer& o) {

            bool join_rhm = sc_min.join(o.sc_min);

            for (size_t i = 0; join_rhm && (i < NB_STREAMCOUNTER_PER_READHASHER); ++i) join_rhm = join_rhm && (sc_km[i].join(o.sc_km[i]));

            return join_rhm;
        }*/
        BFG_INLINE bool join(const ReadHasherMinimizer& o) {

            const bool join_km = sc_km.join(o.sc_km);
            const bool join_min = sc_min.join(o.sc_min);

            return (join_km && join_min);
        }

        BFG_INLINE bool join(const ReadQualityHasherMinimizer& o);

        /*BFG_INLINE size_t KmerF0() const {

            StreamCounter sc(sc_km[0]);

            for (size_t i = 1; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc.join(sc_km[i]);

            return sc.F0();
        }

        BFG_INLINE size_t Kmerf1() const {

            StreamCounter sc(sc_km[0]);

            for (size_t i = 1; i < NB_STREAMCOUNTER_PER_READHASHER; ++i) sc.join(sc_km[i]);

            return sc.f1();
        }*/
        BFG_INLINE size_t KmerF0() const { return sc_km.F0(); }
        BFG_INLINE size_t Kmerf1() const { return sc_km.f1(); }

        BFG_INLINE size_t MinimizerF0() const {

            return sc_min.F0();
        }

        BFG_INLINE size_t Minimizerf1() const {

            return sc_min.f1();
        }

    private:

        size_t k;
        size_t g;

        StreamCounter sc_min;
        StreamCounter sc_km/*[NB_STREAMCOUNTER_PER_READHASHER]*/;
};

/*inline bool ReadQualityHasherMinimizer::join(const ReadHasherMinimizer& o) {

    bool join_rhm = sc_min.join(o.sc_min);

    for (size_t i = 0; join_rhm && (i < NB_STREAMCOUNTER_PER_READHASHER); ++i) join_rhm = join_rhm && (sc_km[i].join(o.sc_km[i]));

    return join_rhm;
}

inline bool ReadHasherMinimizer::join(const ReadQualityHasherMinimizer& o) {

    bool join_rhm = sc_min.join(o.sc_min);

    for (size_t i = 0; join_rhm && (i < NB_STREAMCOUNTER_PER_READHASHER); ++i) join_rhm = join_rhm && (sc_km[i].join(o.sc_km[i]));

    return join_rhm;
}*/
BFG_INLINE bool ReadQualityHasherMinimizer::join(const ReadHasherMinimizer& o) {

    const bool join_km = sc_km.join(o.sc_km);
    const bool join_min = sc_min.join(o.sc_min);

    return (join_km && join_min);
}

BFG_INLINE bool ReadHasherMinimizer::join(const ReadQualityHasherMinimizer& o) {

    const bool join_km = sc_km.join(o.sc_km);
    const bool join_min = sc_min.join(o.sc_min);

    return (join_km && join_min);
}

class KmerStream {

    public:

        KmerStream(const KmerStream_Build_opt& opt) :   k(opt.k), g(opt.g), q(opt.q), q_base(opt.q_base), e(opt.e), km_F1(0), min_F1(0),
                                                        nb_threads(opt.threads), chunksize(opt.chunksize), invalid(false), verbose(opt.verbose) {

            const size_t max_threads = std::thread::hardware_concurrency();

            if (nb_threads <= 0){

                cerr << "KmerStream::KmerStream(): Number of threads cannot be less than or equal to 0" << endl;
                invalid = true;
            }

            if (nb_threads > max_threads){

                cerr << "KmerStream::KmerStream(): Number of threads cannot be greater than or equal to " << max_threads << endl;
                invalid = true;
            }

            if (k == 0){

                cerr << "KmerStream::KmerStream(): Length k of k-mers cannot be less than or equal to 0" << endl;
                invalid = true;
            }

            if (g <= 0){

                cerr << "KmerStream::KmerStream(): Length g of minimizers cannot be less than or equal to 0" << endl;
                invalid = true;
            }

            if (g > k){

                cerr << "KmerStream::KmerStream(): Length g of minimizers cannot be greater than length k of k-mers" << endl;
                invalid = true;
            }

            if (e <= 0){

                cerr << "KmerStream::KmerStream(): Guaranteed error rate cannot be less than or equal to 0" << endl;
                invalid = true;
            }

            if ((q_base != 33) && (q_base != 64)){

                cerr << "KmerStream::KmerStream(): Quality score can only be PHREAD+64 (q_base=64) or PHREAD+33 (q_base=33)" << endl;
                invalid = true;
            }

            if (opt.files.size() == 0) {

                cerr << "KmerStream::KmerStream(): Missing input files" << endl;
                invalid = true;
            }
            else {

                struct stat stFileInfo;
                int intStat;

                for (const auto& s : opt.files) {

                    intStat = stat(s.c_str(), &stFileInfo);

                    if (intStat != 0) {

                        cerr << "KmerStream::KmerStream(): File not found: " << s << endl;
                        invalid = true;
                    }
                    else {

                        const int format = FileParser::getFileFormat(s.c_str());

                        if (format == -1){

                            cerr << "KmerStream::KmerStream(): Input file " << s << " does not exist, is ill-formed or is not in FASTA/FASTQ/GFA format." << endl;

                            invalid = true;
                        }
                        else if (format == 0) files_no_quality.push_back(s); // FASTA
                        else if (format == 1) files_with_quality.push_back(s); // FASTQ
                        else if (format == 2) files_no_quality.push_back(s); // GFA
                    }
                }
            }

            if (invalid) exit(1);

            {
                if (verbose) cout << "KmerStream::KmerStream(): Start computing k-mer cardinality estimations (1/2)" << endl;

                rqh.initialize(e, q, q_base, k, g);

                if (files_with_quality.size() != 0) nb_threads > 1 ? RunThreadedQualityStream() : RunQualityStream();

                rsh = move(rqh);

                if (verbose) cout << "KmerStream::KmerStream(): Start computing k-mer cardinality estimations (1/2)" << endl;

                if (files_no_quality.size() != 0) nb_threads > 1 ? RunThreadedSequenceStream() : RunSequenceStream();

                if (verbose) cout << "KmerStream::KmerStream(): Finished" << endl;

                //rsh.join(rqh);
            }
        }

        BFG_INLINE size_t KmerF0() const { return rsh.KmerF0(); }

        BFG_INLINE size_t KmerF1() const { return km_F1; }

        BFG_INLINE size_t Kmerf1() const { return rsh.Kmerf1(); }

        BFG_INLINE size_t MinimizerF0() const { return rsh.MinimizerF0(); }

        BFG_INLINE size_t MinimizerF1() const { return min_F1; }

        BFG_INLINE size_t Minimizerf1() const { return rsh.Minimizerf1(); }

    private:

        void RunQualityStream() {

            FileParser fp(files_with_quality);

            size_t file_id = 0;

            string seq;

            while (fp.read(seq, file_id)){

                const char* qss = fp.getQualityScoreString();

                std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

                const pair<size_t, size_t> F1 = rqh.update(seq.c_str(), seq.length(), qss, strlen(qss));

                min_F1 += F1.first;
                km_F1 += F1.second; 
            }

            fp.close();
        }

        void RunThreadedQualityStream() {

            size_t pos_read = 0;
            size_t len_read = 0;

            const size_t max_len_seq = rndup(static_cast<size_t>(1024 + k - 1));
            const size_t thread_seq_buf_sz = BUFFER_SIZE;

            string seq, qual;

            FileParser fp(files_with_quality);

            auto reading_function = [&](char* seq_buf, char* qual_buf, size_t& buf_sz) {

                size_t file_id = 0;

                const size_t seq_buf_sz = thread_seq_buf_sz - k;

                const char* s_str = seq.c_str();
                const char* q_str = qual.c_str();

                buf_sz = 0;

                while (buf_sz < seq_buf_sz) {

                    const bool new_reading = (pos_read >= len_read);

                    if (!new_reading || fp.read(seq, file_id)) {

                        if (new_reading) qual = fp.getQualityScoreString();

                        pos_read &= static_cast<size_t>(new_reading) - 1;

                        len_read = seq.length();
                        s_str = seq.c_str();
                        q_str = qual.c_str();

                        if (len_read >= k){

                            if ((thread_seq_buf_sz - buf_sz - 1) < (len_read - pos_read)){

                                strncpy(&seq_buf[buf_sz], &s_str[pos_read], thread_seq_buf_sz - buf_sz - 1);
                                strncpy(&qual_buf[buf_sz], &q_str[pos_read], thread_seq_buf_sz - buf_sz - 1);

                                seq_buf[thread_seq_buf_sz - 1] = '\0';

                                pos_read += seq_buf_sz - buf_sz;
                                buf_sz = thread_seq_buf_sz;

                                break;
                            }
                            else {

                                strcpy(&seq_buf[buf_sz], &s_str[pos_read]);
                                strcpy(&qual_buf[buf_sz], &q_str[pos_read]);

                                buf_sz += (len_read - pos_read) + 1;
                                pos_read = len_read;
                            }
                        }
                        else pos_read = len_read;
                    }
                    else return true;
                }

                return false;
            };

            {
                vector<thread> workers; // need to keep track of threads so we can join them

                bool stop = false;

                mutex mutex_file;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            char* buffer_seq = new char[thread_seq_buf_sz];
                            char* buffer_qual = new char[thread_seq_buf_sz];

                            size_t buffer_sz = 0;
                            
                            size_t l_min_F1 = 0;
                            size_t l_km_F1 = 0;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    if (stop){

                                        min_F1 += l_min_F1;
                                        km_F1 += l_km_F1;

                                        delete[] buffer_seq;
                                        delete[] buffer_qual;

                                        return;
                                    }

                                    stop = reading_function(buffer_seq, buffer_qual, buffer_sz);
                                }

                                for (char* s = buffer_seq; s != (buffer_seq + thread_seq_buf_sz); ++s) *s &= 0xDF;

                                const pair<size_t, size_t> F1 = rqh.update_p(buffer_seq, buffer_qual, buffer_sz);

                                l_min_F1 += F1.first;
                                l_km_F1 += F1.second;
                            }

                            delete[] buffer_seq;
                            delete[] buffer_qual;
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            fp.close();
        }

        void RunSequenceStream() {

            size_t file_id = 0;

            string seq;

            FileParser fp(files_no_quality);

            while (fp.read(seq, file_id)){

                std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

                const pair<size_t, size_t> F1 = rsh.update(seq.c_str(), seq.length());

                min_F1 += F1.first;
                km_F1 += F1.second; 
            }

            fp.close();
        }

        void RunThreadedSequenceStream() {

            size_t pos_read = 0;
            size_t len_read = 0;

            const size_t thread_seq_buf_sz = BUFFER_SIZE;

            FileParser fp(files_no_quality);

            string s;

            auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

                size_t file_id = 0;

                const size_t sz_buf = thread_seq_buf_sz - k;

                const char* s_str = s.c_str();

                seq_buf_sz = 0;

                while (seq_buf_sz < sz_buf) {

                    const bool new_reading = (pos_read >= len_read);

                    if (!new_reading || fp.read(s, file_id)) {

                        pos_read &= static_cast<size_t>(new_reading) - 1;

                        len_read = s.length();
                        s_str = s.c_str();

                        if (len_read >= k){

                            if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                                strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                                seq_buf[thread_seq_buf_sz - 1] = '\0';

                                pos_read += sz_buf - seq_buf_sz;
                                seq_buf_sz = thread_seq_buf_sz;

                                break;
                            }
                            else {

                                strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                                seq_buf_sz += (len_read - pos_read) + 1;
                                pos_read = len_read;
                            }
                        }
                        else pos_read = len_read;
                    }
                    else return true;
                }

                return false;
            };

            {
                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_file;

                bool stop = false;

                for (size_t t = 0; t != nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            char* buffer_seq = new char[thread_seq_buf_sz];

                            size_t buffer_seq_sz = 0;
                            size_t l_min_F1 = 0;
                            size_t l_km_F1 = 0;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    if (stop) {

                                        min_F1 += l_min_F1;
                                        km_F1 += l_km_F1;

                                        delete[] buffer_seq;

                                        return;
                                    }

                                    stop = reading_function(buffer_seq, buffer_seq_sz);
                                }

                                for (char* s = buffer_seq; s != (buffer_seq + thread_seq_buf_sz); ++s) *s &= 0xDF;

                                const pair<size_t, size_t> F1 = rsh.update_p(buffer_seq, buffer_seq_sz);

                                l_min_F1 += F1.first;
                                l_km_F1 += F1.second;
                            }

                            delete[] buffer_seq;
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            fp.close();
        }

        size_t k;
        size_t g;
        size_t q;
        size_t q_base;

        size_t km_F1;
        size_t min_F1;

        double e;

        ReadQualityHasherMinimizer rqh;
        ReadHasherMinimizer rsh;

        vector<string> files_no_quality;
        vector<string> files_with_quality;

        bool verbose;
        bool invalid;

        size_t nb_threads;
        size_t chunksize;
};

#endif
