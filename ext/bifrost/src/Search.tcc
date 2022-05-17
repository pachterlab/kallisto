#ifndef BIFROST_SEARCH_DBG_TCC
#define BIFROST_SEARCH_DBG_TCC

template<typename U, typename G>
vector<pair<size_t, UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence( const string& s, const bool exact, const bool insertion,
                                                                                const bool deletion, const bool substitution,
                                                                                const bool or_exclusive_match) {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    Roaring rpos;

    vector<pair<size_t, UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, UnitigMap<U, G>>& p1, const pair<size_t, UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!or_exclusive_match || (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end())))) {

                            const UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                ki_s += um.len - 1;
            }
        }

        if (or_exclusive_match && (insertion || deletion || substitution)){

            for (const auto& pum : v_um) {

                us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
                rpos.add(pum.first);
            }
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_func(true, false, false, i);
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, true, false, i);
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, false, true, i);
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence(   const string& s, const bool exact, const bool insertion,
                                                                            const bool deletion, const bool substitution,
                                                                            const double ratio_kmers, const bool or_exclusive_match) {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (ratio_kmers < 0.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is less than 0.0" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (ratio_kmers > 1.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is greater than 1.0" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, UnitigMap<U, G>>>();
    }

    const size_t nb_km_min = static_cast<double>(s.length() - k_ + 1) * ratio_kmers;

    Roaring rpos;

    vector<pair<size_t, UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, UnitigMap<U, G>>& p1, const pair<size_t, UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end()))) {

                            const UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                if (rpos.cardinality() >= nb_km_min) return;

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                if (v_um.size() >= nb_km_min) return v_um;

                ki_s += um.len - 1;
            }
        }

        for (const auto& pum : v_um) {

            us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
            rpos.add(pum.first);
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_func(true, false, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, true, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, false, true, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, const_UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence( const string& s, const bool exact, const bool insertion,
                                                                                const bool deletion, const bool substitution,
                                                                                const bool or_exclusive_match) const {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    Roaring rpos;

    vector<pair<size_t, const_UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, const_UnitigMap<U, G>>& p1, const pair<size_t, const_UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const const_UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!or_exclusive_match || (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end())))) {

                            const const_UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const const_UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                ki_s += um.len - 1;
            }
        }

        if (or_exclusive_match && (insertion || deletion || substitution)){

            for (const auto& pum : v_um) {

                us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
                rpos.add(pum.first);
            }
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_func(true, false, false, i);
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, true, false, i);
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, false, true, i);
        }
    }

    return v_um;
}

template<typename U, typename G>
vector<pair<size_t, const_UnitigMap<U, G>>> CompactedDBG<U, G>::searchSequence( const string& s, const bool exact, const bool insertion,
                                                                                const bool deletion, const bool substitution,
                                                                                const double ratio_kmers, const bool or_exclusive_match) const {

    struct hash_pair {

        size_t operator()(const pair<size_t, Kmer>& p) const {

            return wyhash(&(p.first), sizeof(size_t), 0, _wyp) ^ p.second.hash();
        }
    };

    if (invalid){

        cerr << "CompactedDBG::searchSequence(): Graph is invalid and cannot be searched" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (ratio_kmers < 0.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is less than 0.0" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (ratio_kmers > 1.0){

        cerr << "CompactedDBG::searchSequence(): Ratio of k-mers is greater than 1.0" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    if (s.length() < k_){

        cerr << "CompactedDBG::searchSequence(): Query length is shorter than k-mer size" << endl;

        return vector<pair<size_t, const_UnitigMap<U, G>>>();
    }

    const size_t nb_km_min = static_cast<double>(s.length() - k_ + 1) * ratio_kmers;

    Roaring rpos;

    vector<pair<size_t, const_UnitigMap<U, G>>> v_um;

    string s_inexact;

    unordered_set<pair<size_t, Kmer>, hash_pair> us_pos_km;

    auto comp_pair = [](const pair<size_t, const_UnitigMap<U, G>>& p1, const pair<size_t, const_UnitigMap<U, G>>& p2) {

        return (p1.first < p2.first);
    };

    auto worker_func = [&](const bool subst, const bool ins, const bool del, const size_t shift){

        const size_t ins_mask = static_cast<size_t>(!ins) - 1;
        const size_t del_mask = static_cast<size_t>(!del) - 1;

        const size_t s_len = s.length();
        const char* s_str = s.c_str();

        const size_t s_inexact_len = s_inexact.length();
        const char* s_inexact_str = s_inexact.c_str();

        const size_t k_1 = k_-1;

        auto processUnitigMap = [&](const const_UnitigMap<U, G>& um, const size_t pos_s){

            if (um.strand){

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + j - um.dist;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
            else {

                for (size_t j = um.dist; j < um.dist + um.len; ++j){

                    size_t l_pos_seq = pos_s + um.dist + um.len - j - 1;

                    const size_t shift_pos_seq = (l_pos_seq / k_) + (l_pos_seq % k_ > shift);

                    l_pos_seq -= (ins_mask & shift_pos_seq);
                    l_pos_seq += (del_mask & shift_pos_seq);

                    if ((l_pos_seq + k_1 < s_len) && us_pos_km.insert({l_pos_seq, um.getMappedKmer(j)}).second) {

                        v_um.push_back({l_pos_seq, um.getKmerMapping(j)});
                        rpos.add(l_pos_seq);
                    }
                }
            }
        };

        for (size_t i = 0; i != ((subst || ins) ? 4 : 1); ++i){

            if (ins) {

                for (size_t j = shift; j < s_inexact_len; j += k_) s_inexact[j] = alpha[i];
            }
            else if (subst) {

                for (size_t j = shift; j < s_inexact_len; j += k_) {

                    if (!isDNA(s[j]) || (alpha[i] == s[j])) s_inexact[j] = 'N';
                    else s_inexact[j] = alpha[i];
                }
            } 

            KmerIterator ki_s(s_inexact_str), ki_e;
            minHashIterator<RepHash> mhi = minHashIterator<RepHash>(s_inexact_str, s_inexact_len, k_, g_, RepHash(), true);

            minHashResultIterator<RepHash> it_min, it_min_end;
            minHashResult mhr;

            Minimizer minz;

            pair<size_t, bool> minz_pres = {0xffffffffffffffffULL, true};

            while (ki_s != ki_e) {

                const size_t pos_s = ki_s->second;

                mhi += (pos_s - mhi.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                it_min = *mhi;
                mhr = *it_min;

                // If minimizers of new kmer are different from minimizers of previous kmer
                // or if minimizers are the same but they were present, search them again
                if (minz_pres.second || (mhr.pos != minz_pres.first)){

                    if (mhr.pos != minz_pres.first){

                        minz = Minimizer(s_inexact_str + mhr.pos).rep();
                        minz_pres = {mhr.pos, hmap_min_unitigs.find(minz) != hmap_min_unitigs.end()};

                        for (++it_min; !minz_pres.second && (it_min != it_min_end); ++it_min){

                            mhr = *it_min;
                            minz = Minimizer(s_inexact_str + mhr.pos).rep();
                            minz_pres.second = (hmap_min_unitigs.find(minz) != hmap_min_unitigs.end());
                        }
                    }

                    if (minz_pres.second) { // If the k-mer has already been searched in the past, discard

                        const size_t shift_pos_seq = (pos_s / k_) + (pos_s % k_ > shift);
                        const size_t l_pos_s = pos_s - (ins_mask & shift_pos_seq) + (del_mask & shift_pos_seq);

                        if ((l_pos_s + k_1 < s_len) && isDNA(s_str[l_pos_s]) && isDNA(s_str[l_pos_s + k_1]) && (!rpos.contains(l_pos_s) && (us_pos_km.find({l_pos_s, ki_s->first}) == us_pos_km.end()))) {

                            const const_UnitigMap<U, G> um = findUnitig(s_inexact_str, pos_s, s_inexact_len, mhi);

                            if (!um.isEmpty){

                                processUnitigMap(um, pos_s);

                                if (rpos.cardinality() >= nb_km_min) return;

                                ki_s += um.len - 1;
                            }
                        }
                    }
                }

                ++ki_s;
            }
        }
    };

    if (exact){

        for (KmerIterator ki_s(s.c_str()), ki_e; ki_s != ki_e; ++ki_s) {

            const size_t pos_s = ki_s->second;
            const const_UnitigMap<U, G> um = findUnitig(s.c_str(), pos_s, s.length());

            if (!um.isEmpty) { // Read maps to a Unitig

                if (um.strand){

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + j - um.dist, um.getKmerMapping(j)});
                }
                else {

                    for (size_t j = um.dist; j < um.dist + um.len; ++j) v_um.push_back({pos_s + um.dist + um.len - j - 1, um.getKmerMapping(j)});
                }

                if (v_um.size() >= nb_km_min) return v_um;

                ki_s += um.len - 1;
            }
        }

        for (const auto& pum : v_um) {

            us_pos_km.insert({pum.first, pum.second.getMappedKmer(pum.second.dist)});
            rpos.add(pum.first);
        }
    }

    if (substitution){

        for (size_t i = 0; i != k_; ++i){

            s_inexact = s;

            worker_func(true, false, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (insertion){

        for (size_t i = 0; i != k_; ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ - 1) == 0) ss << alpha[0];

                ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, true, false, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    if (deletion && (s.length() >= (k_ + 1))){

        for (size_t i = 0; i != (k_ + 1); ++i){

            std::stringstream ss;

            for (size_t j = 0; j < i; ++j) ss << s[j];

            for (size_t j = i, cpt = 0; j < s.length(); ++j, ++cpt) {

                if (cpt % (k_ + 1) != 0) ss << s[j];
            }

            s_inexact = ss.str();

            worker_func(false, false, true, i);

            if (rpos.cardinality() >= nb_km_min) return v_um;
        }
    }

    return v_um;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::search(const vector<string>& query_filenames, const string& out_filename_prefix,
                                const double ratio_kmers, const bool inexact_search, const size_t nb_threads,
                                const size_t verbose) const {

     if (invalid){

        cerr << "CompactedDBG::search(): Graph is invalid and cannot be searched" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "CompactedDBG::search(): Number of threads cannot be greater than or equal to " << std::thread::hardware_concurrency() << "." << endl;
        return false;
    }

    if (nb_threads <= 0){

        cerr << "CompactedDBG::search(): Number of threads cannot be less than or equal to 0." << endl;
        return false;
    }

    const string out_tmp = out_filename_prefix + ".tsv";

    FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

    if (fp_tmp == NULL) {

        cerr << "CompactedDBG::search(): Could not open file " << out_tmp << " for writing." << endl;
        return false;
    }
    else {

        fclose(fp_tmp);

        if (std::remove(out_tmp.c_str()) != 0) cerr << "CompactedDBG::search(): Could not remove temporary file " << out_tmp << endl;
    }

    if (verbose) cout << "CompactedDBG::search(): Querying graph." << endl;

    const CompactedDBG<U, G>& dbg = *this;

    string s;

    size_t file_id = 0;

    //const size_t max_len_seq = 1024;
    //const size_t thread_seq_buf_sz = 64 * max_len_seq;
    const size_t max_len_seq = rndup(static_cast<size_t>(1024 + k_ - 1));
    const size_t thread_seq_buf_sz = BUFFER_SIZE;

    FileParser fp(query_filenames);

    ofstream outfile;
    ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());
    //out.sync_with_stdio(false);

    const char query_pres[3] = {'\t', '1', '\n'};
    const char query_abs[3] = {'\t', '0', '\n'};

    const size_t l_query_res = 3;

    // Write header to TSV file
    out << "query_name\tpresence_query\n";

    if (nb_threads == 1){

        char* buffer_res = new char[thread_seq_buf_sz];

        size_t pos_buffer_out = 0;
        size_t nb_queries_found = 0;

        while (fp.read(s, file_id)){

            bool is_found = false;

            const size_t nb_km_min = static_cast<double>(s.length() - k_ + 1) * ratio_kmers;
            const char* query_name = fp.getNameString();
            const size_t l_query_name = strlen(query_name);

            for (auto& c : s) c &= 0xDF;

            const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   s, true, inexact_search, inexact_search,
                                                                                        inexact_search, ratio_kmers, true);

            if (inexact_search){

                Roaring r;

                for (const auto& p : v) r.add(p.first);

                is_found = (r.cardinality() >= nb_km_min);
            }
            else is_found = (v.size() >= nb_km_min);

            if (pos_buffer_out + l_query_name + l_query_res >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                out.write(buffer_res, pos_buffer_out); // Write result buffer
                pos_buffer_out = 0; // Reset position to 0;
            }

            // Copy new result to buffer
            std::memcpy(buffer_res + pos_buffer_out, query_name, l_query_name * sizeof(char));

            if (is_found){

                std::memcpy(buffer_res + pos_buffer_out + l_query_name, query_pres, l_query_res * sizeof(char));

                ++nb_queries_found;
            }
            else std::memcpy(buffer_res + pos_buffer_out + l_query_name, query_abs, l_query_res * sizeof(char));

            pos_buffer_out += l_query_name + l_query_res;
        }

        // Flush unresult written to final output
        if (pos_buffer_out > 0) out.write(buffer_res, pos_buffer_out);

        delete[] buffer_res;

        if (verbose) cout << "CompactedDBG::search(): Found " << nb_queries_found << " queries. " << endl;
    }
    else {

        {
            bool stop = false;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_files_in, mutex_file_out;

            std::atomic<size_t> nb_queries_found;

            nb_queries_found = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        char* buffer_res = new char[thread_seq_buf_sz];

                        vector<string> buffers_seq;
                        vector<string> buffers_name;

                        while (true) {

                            {
                                if (stop) {

                                    delete[] buffer_res;

                                    return;
                                }

                                size_t buffer_sz = 0;

                                unique_lock<mutex> lock(mutex_files_in);

                                stop = !fp.read(s, file_id);

                                while (!stop){

                                    buffer_sz += s.length();

                                    buffers_seq.push_back(std::move(s));
                                    buffers_name.push_back(string(fp.getNameString()));

                                    if (buffer_sz >= thread_seq_buf_sz) break;
                                    else stop = !fp.read(s, file_id);
                                }
                            }

                            size_t pos_buffer_out = 0;

                            const size_t buffers_seq_sz = buffers_seq.size();

                            for (size_t i = 0; i < buffers_seq_sz; ++i){

                                bool is_found = false;

                                const size_t nb_km_min = static_cast<double>(buffers_seq[i].length() - k_ + 1) * ratio_kmers;
                                const size_t l_name = buffers_name[i].length();

                                for (auto& c : buffers_seq[i]) c &= 0xDF;

                                const vector<pair<size_t, const_UnitigMap<U, G>>> v = dbg.searchSequence(   buffers_seq[i], true, inexact_search, inexact_search,
                                                                                                            inexact_search, ratio_kmers, true);

                                if (inexact_search){

                                    Roaring r;

                                    for (const auto& p : v) r.add(p.first);

                                    is_found = (r.cardinality() >= nb_km_min);
                                }
                                else is_found = (v.size() >= nb_km_min);

                                if (pos_buffer_out + l_name + l_query_res >= thread_seq_buf_sz){ // If next result cannot fit in the buffer

                                    unique_lock<mutex> lock(mutex_file_out); // Get the output lock

                                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                                    pos_buffer_out = 0; // Reset position to 0;
                                }

                                // Copy new result to buffer
                                std::memcpy(buffer_res + pos_buffer_out, buffers_name[i].c_str(), l_name * sizeof(char));

                                if (is_found){

                                    std::memcpy(buffer_res + pos_buffer_out + l_name, query_pres, l_query_res * sizeof(char));

                                    ++nb_queries_found;
                                }
                                else std::memcpy(buffer_res + pos_buffer_out + l_name, query_abs, l_query_res * sizeof(char));

                                pos_buffer_out += l_name + l_query_res;
                            }

                            if (pos_buffer_out > 0){ // Flush unresult written to final output

                                unique_lock<mutex> lock(mutex_file_out);

                                out.write(buffer_res, pos_buffer_out);
                            }

                            // Clear buffers for next round
                            buffers_seq.clear();
                            buffers_name.clear();
                        }

                        delete[] buffer_res;
                    }
                );
            }

            for (auto& t : workers) t.join();

            if (verbose) cout << "CompactedDBG::search(): Found " << nb_queries_found << " queries. " << endl;
        }
    }

    outfile.close();
    fp.close();

    return true;
}

#endif
