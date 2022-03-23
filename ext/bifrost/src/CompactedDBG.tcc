#ifndef BIFROST_COMPACTED_DBG_TCC
#define BIFROST_COMPACTED_DBG_TCC

static const uint8_t bits[256] = {
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

template<typename U, typename G>
CompactedDBG<U, G>::CompactedDBG(const int kmer_length, const int minimizer_length) : invalid(false) {

    setKmerGmerLength(kmer_length, minimizer_length);
}

template<typename U, typename G>
CompactedDBG<U, G>::CompactedDBG(const CompactedDBG<U, G>& o) : k_(o.k_), g_(o.g_), invalid(o.invalid),
                                                                bf(o.bf), km_unitigs(o.km_unitigs), v_unitigs(o.v_unitigs.size(), nullptr),
                                                                data(o.data), h_kmers_ccov(o.h_kmers_ccov),
                                                                hmap_min_unitigs(o.hmap_min_unitigs){

    for (size_t i = 0; i < o.v_unitigs.size(); ++i){

        v_unitigs[i] = new Unitig<U>;
        *(v_unitigs[i]) = *(o.v_unitigs[i]);
    }
}

template<typename U, typename G>
CompactedDBG<U, G>::CompactedDBG(CompactedDBG<U, G>&& o) :  k_(o.k_), g_(o.g_), invalid(o.invalid),
                                                            bf(std::move(o.bf)), km_unitigs(std::move(o.km_unitigs)), data(std::move(o.data)),
                                                            v_unitigs(std::move(o.v_unitigs)), h_kmers_ccov(std::move(o.h_kmers_ccov)),
                                                            hmap_min_unitigs(std::move(o.hmap_min_unitigs)){

    o.clear();
}

template<typename U, typename G>
CompactedDBG<U, G>::~CompactedDBG() {

    clear();
}

template<typename U, typename G>
CompactedDBG<U, G>& CompactedDBG<U, G>::operator=(const CompactedDBG<U, G>& o){

    clear();

    k_ = o.k_;
    g_ = o.g_;

    invalid = o.invalid;

    km_unitigs = o.km_unitigs;

    h_kmers_ccov = o.h_kmers_ccov;
    hmap_min_unitigs = o.hmap_min_unitigs;

    bf = o.bf;

    data = o.data;

    v_unitigs = vector<Unitig<U>*>(o.v_unitigs.size(), nullptr);

    for (size_t i = 0; i < o.v_unitigs.size(); ++i){

        v_unitigs[i] = new Unitig<U>;
        *(v_unitigs[i]) = *(o.v_unitigs[i]);
    }

    return *this;
}

template<typename U, typename G>
CompactedDBG<U, G>& CompactedDBG<U, G>::toDataGraph(CompactedDBG<void, void>&& o, const size_t nb_threads) {

    clear();

    k_ = o.k_;
    g_ = o.g_;
    invalid = o.invalid;

    km_unitigs.toData(std::move(o.km_unitigs), nb_threads);

    hmap_min_unitigs = std::move(o.hmap_min_unitigs);
    bf = std::move(o.bf);

    data = wrapperData<G>();

    v_unitigs = vector<Unitig<U>*>(o.v_unitigs.size(), nullptr);

    auto moveUnitigs = [&](const size_t start, const size_t end){

        for (size_t i = start; i < end; ++i){

            v_unitigs[i] = new Unitig<U>(std::move(o.v_unitigs[i]->getSeq()), std::move(o.v_unitigs[i]->getCov()));

            delete o.v_unitigs[i];
        }
    };

    if ((nb_threads == 1) || (v_unitigs.size() < 1024)) moveUnitigs(0, v_unitigs.size());
    else {

        vector<thread> workers;

        const size_t slice = (v_unitigs.size() / nb_threads) + 1;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    const size_t start = t * slice;
                    const size_t end = min(start + slice, v_unitigs.size());

                    if (start < v_unitigs.size()) moveUnitigs(start, end);
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    o.v_unitigs.clear();

    KmerHashTable<CompressedCoverage_t<void>>::const_iterator it_s = o.h_kmers_ccov.begin();
    KmerHashTable<CompressedCoverage_t<void>>::const_iterator it_e = o.h_kmers_ccov.end();

    h_kmers_ccov = KmerHashTable<CompressedCoverage_t<U>>(o.h_kmers_ccov.size());

    while (it_s != it_e) {

        Kmer km = it_s.getKey();

        h_kmers_ccov.insert(move(km), move(it_s->ccov));

        ++it_s;
    }

    o.h_kmers_ccov.clear();
    o.clear();

    return *this;
}

template<typename U, typename G>
CompactedDBG<U, G>& CompactedDBG<U, G>::operator=(CompactedDBG<U, G>&& o){

    if (this != &o) {

        clear();

        k_ = o.k_;
        g_ = o.g_;

        invalid = o.invalid;

        km_unitigs = std::move(o.km_unitigs);
        v_unitigs = std::move(o.v_unitigs);

        h_kmers_ccov = std::move(o.h_kmers_ccov);
        hmap_min_unitigs = std::move(o.hmap_min_unitigs);

        bf = std::move(o.bf);

        data = std::move(o.data);

        o.clear();
    }

    return *this;
}

template<typename U, typename G>
CompactedDBG<U, G>& CompactedDBG<U, G>::operator+=(const CompactedDBG<U, G>& o){

    merge(o, false);

    return *this;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::operator==(const CompactedDBG<U, G>& o) const {

    if ((!invalid && !o.invalid) && (k_ == o.k_) && (size() == o.size())){

        for (const auto& unitig : *this){

            const_UnitigMap<U, G> unitig_o(o.find(unitig.getUnitigHead(), true));

            if (unitig_o.isEmpty) return false;
            else {

                unitig_o.dist = 0;
                unitig_o.len = unitig_o.size - k_ + 1;

                const string unitig_o_str = unitig_o.strand ? unitig_o.referenceUnitigToString() : reverse_complement(unitig_o.referenceUnitigToString());

                if (unitig_o_str != unitig.referenceUnitigToString()) return false;
            }
        }

        return true;
    }

    return false;
}

template<typename U, typename G>
inline bool CompactedDBG<U, G>::operator!=(const CompactedDBG<U, G>& o) const {

    return !operator==(o);
}

template<typename U, typename G>
void CompactedDBG<U, G>::clear(){

    k_ = 0;
    g_ = 0;

    invalid = true;

    for (auto unitig : v_unitigs) delete unitig;

    v_unitigs.clear();
    km_unitigs.clear();
    hmap_min_unitigs.clear();
    h_kmers_ccov.clear();
    bf.clear();
}

/*template<typename U, typename G>
bool CompactedDBG<U, G>::build(CDBG_Build_opt& opt){

    size_t max_threads = std::thread::hardware_concurrency();

    bool construct_finished = true;

    const int k_cpy = k_;
    const int g_cpy = g_;
    const bool invalid_cpy = invalid;

    if (invalid){

        cerr << "CompactedDBG::build(): Graph is invalid and cannot be built" << endl;
        construct_finished = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "CompactedDBG::build(): Number of threads cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "CompactedDBG::build(): Number of threads cannot exceed " << max_threads << "threads" << endl;
        construct_finished = false;
    }

    if (opt.outFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
            construct_finished = false;
        }
        else {

            fclose(fp);

            if (std::remove(opt.outFilenameBBF.c_str()) != 0){

                cerr << "CompactedDBG::build(): Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
            }
        }
    }

    if (opt.inFilenameBBF.length() != 0){

        if (check_file_exists(opt.inFilenameBBF)){

            FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

            if (fp == NULL) {

                cerr << "CompactedDBG::build(): Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                construct_finished = false;
            }
            else fclose(fp);
        }
        else {

            cerr << "CompactedDBG::build(): Input Blocked Bloom filter file " << opt.inFilenameBBF << " does not exist." << endl;
            construct_finished = false;
        }
    }

    if (opt.filename_seq_in.size() + opt.filename_ref_in.size() == 0){

        cerr << "CompactedDBG::build(): Number of FASTA/FASTQ files in input cannot be 0." << endl;
        construct_finished = false;
    }
    else {

        for (const auto& file : opt.filename_seq_in){

            if (check_file_exists(file)){

                FILE* fp = fopen(file.c_str(), "r");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << file << endl;
                    construct_finished = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "CompactedDBG::build(): Input file " << opt.inFilenameBBF << " does not exist." << endl;
                construct_finished = false;
            }
        }

        for (const auto& file : opt.filename_ref_in){

            if (check_file_exists(file)){

                FILE* fp = fopen(file.c_str(), "r");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << file << endl;
                    construct_finished = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "CompactedDBG::build(): Input file " << opt.inFilenameBBF << " does not exist." << endl;
                construct_finished = false;
            }
        }
    }

    clear();

    k_ = k_cpy;
    g_ = g_cpy;
    invalid = invalid_cpy;

    if (construct_finished){

        if ((opt.filename_seq_in.size() != 0) && (opt.filename_ref_in.size() != 0)){

            CompactedDBG<U, G> graph_seq(k_, g_);
            CompactedDBG<U, G> graph_ref(k_, g_);

            CDBG_Build_opt opt_seq(opt);
            CDBG_Build_opt opt_ref(opt);

            opt_seq.filename_ref_in.clear();
            opt_ref.filename_seq_in.clear();

            construct_finished = graph_seq.build(opt_seq);

            if (construct_finished) construct_finished = graph_ref.build(opt_ref);

            if (construct_finished){

                if (graph_ref.length() < graph_seq.length()){

                    construct_finished = graph_seq.merge(std::move(graph_ref), opt.nb_threads, opt.verbose);

                    if (construct_finished) *this = std::move(graph_seq);
                }
                else {

                    construct_finished = graph_ref.merge(std::move(graph_seq), opt.nb_threads, opt.verbose);

                    if (construct_finished) *this = std::move(graph_ref);
                }
            }

            setFullCoverage(2);
        }
        else {

            const bool reference_mode = (opt.filename_ref_in.size() != 0);

            const vector<string>& v_files = reference_mode ? opt.filename_ref_in : opt.filename_seq_in;

            setFullCoverage(reference_mode ? 1 : 2);

            KmerStream_Build_opt kms_opt;

            kms_opt.threads = opt.nb_threads;
            kms_opt.verbose = opt.verbose;
            kms_opt.k = k_;
            kms_opt.g = g_;
            kms_opt.q = 0;

            for (const auto& s : v_files) kms_opt.files.push_back(s);

            KmerStream kms(kms_opt);

            const size_t nb_unique_kmers = max(1UL, kms.KmerF0());
            const size_t nb_unique_minimizers = max(1UL, kms.MinimizerF0());
            const size_t nb_non_unique_kmers = reference_mode ? 0 : max(1UL, nb_unique_kmers - min(nb_unique_kmers, kms.Kmerf1()));
            const size_t nb_non_unique_minimizers = reference_mode ? 0 : max(1UL, nb_unique_minimizers - min(nb_unique_minimizers, kms.Minimizerf1()));

            if (opt.verbose){

                cout << "CompactedDBG::build(): Estimated number of k-mers occurring at least once: " << nb_unique_kmers << endl;
                cout << "CompactedDBG::build(): Estimated number of minimizer occurring at least once: " << nb_unique_minimizers << endl;

                if (!reference_mode){

                    cout << "CompactedDBG::build(): Estimated number of k-mers occurring twice or more: " << nb_non_unique_kmers << endl;
                    cout << "CompactedDBG::build(): Estimated number of minimizers occurring twice or more: " << nb_non_unique_minimizers << endl;
                }
            }

            if (opt.inFilenameBBF.length() != 0){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    construct_finished = false;
                }
                else {

                    construct_finished = bf.ReadBloomFilter(fp);

                    fclose(fp);
                }
            }
            else construct_finished = filter(opt, nb_unique_kmers, nb_non_unique_kmers);

            if (construct_finished){

                if (opt.outFilenameBBF.length() != 0){

                    FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

                    if (fp == NULL) {

                        cerr << "CompactedDBG::build(): Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                        construct_finished = false;
                    }
                    else {

                        bf.WriteBloomFilter(fp);

                        fclose(fp);
                    }
                }

                if (construct_finished) construct_finished = construct(opt, nb_unique_minimizers, nb_non_unique_minimizers); // Construction step

                bf.clear();
            }
        }
    }

    return construct_finished;
}*/

template<typename U, typename G>
bool CompactedDBG<U, G>::build(CDBG_Build_opt& opt){

    size_t max_threads = std::thread::hardware_concurrency();

    bool construct_finished = true;

    const int k_cpy = k_;
    const int g_cpy = g_;
    const bool invalid_cpy = invalid;

    if (invalid){

        cerr << "CompactedDBG::build(): Graph is invalid and cannot be built" << endl;
        construct_finished = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "CompactedDBG::build(): Number of threads cannot be less than or equal to 0" << endl;
        construct_finished = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "CompactedDBG::build(): Number of threads cannot exceed " << max_threads << "threads" << endl;
        construct_finished = false;
    }

    if (opt.outFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "CompactedDBG::build(): Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
            construct_finished = false;
        }
        else {

            fclose(fp);

            if (std::remove(opt.outFilenameBBF.c_str()) != 0){

                cerr << "CompactedDBG::build(): Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
            }
        }
    }

    if (opt.inFilenameBBF.length() != 0){

        if (check_file_exists(opt.inFilenameBBF)){

            FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

            if (fp == NULL) {

                cerr << "CompactedDBG::build(): Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                construct_finished = false;
            }
            else fclose(fp);
        }
        else {

            cerr << "CompactedDBG::build(): Input Blocked Bloom filter file " << opt.inFilenameBBF << " does not exist." << endl;
            construct_finished = false;
        }
    }

    if (opt.filename_seq_in.size() + opt.filename_ref_in.size() == 0){

        cerr << "CompactedDBG::build(): Number of FASTA/FASTQ files in input cannot be 0." << endl;
        construct_finished = false;
    }
    else {

        for (const auto& file : opt.filename_seq_in){

            if (check_file_exists(file)){

                FILE* fp = fopen(file.c_str(), "r");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << file << endl;
                    construct_finished = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "CompactedDBG::build(): Input file " << opt.inFilenameBBF << " does not exist." << endl;
                construct_finished = false;
            }
        }

        for (const auto& file : opt.filename_ref_in){

            if (check_file_exists(file)){

                FILE* fp = fopen(file.c_str(), "r");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input FASTA/FASTQ file " << file << endl;
                    construct_finished = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "CompactedDBG::build(): Input file " << opt.inFilenameBBF << " does not exist." << endl;
                construct_finished = false;
            }
        }
    }

    clear();

    k_ = k_cpy;
    g_ = g_cpy;
    invalid = invalid_cpy;

    if (construct_finished){

        if ((opt.filename_seq_in.size() != 0) && (opt.filename_ref_in.size() != 0)){

            CompactedDBG<void, void> graph_seq(k_, g_);
            CompactedDBG<void, void> graph_ref(k_, g_);

            CDBG_Build_opt opt_seq(opt);
            CDBG_Build_opt opt_ref(opt);

            opt_seq.filename_ref_in.clear();
            opt_ref.filename_seq_in.clear();

            construct_finished = graph_seq.build(opt_seq);

            if (construct_finished) construct_finished = graph_ref.build(opt_ref);

            if (construct_finished){

                if (graph_ref.length() < graph_seq.length()){

                    construct_finished = graph_seq.merge(std::move(graph_ref), opt.nb_threads, opt.verbose);

                    if (construct_finished) toDataGraph(std::move(graph_seq), opt.nb_threads);
                }
                else {

                    construct_finished = graph_ref.merge(std::move(graph_seq), opt.nb_threads, opt.verbose);

                    if (construct_finished) toDataGraph(std::move(graph_ref), opt.nb_threads);
                }
            }

            setFullCoverage(2);
        }
        else if (!is_void<U>::value) {

            CompactedDBG<void, void> graph(k_, g_);

            construct_finished = graph.build(opt);

            if (construct_finished) toDataGraph(std::move(graph), opt.nb_threads);
        }
        else {

            const bool reference_mode = (opt.filename_ref_in.size() != 0);

            size_t nb_unique_kmers, nb_unique_minimizers;
            size_t nb_non_unique_kmers, nb_non_unique_minimizers;

            {
                const vector<string>& v_files = reference_mode ? opt.filename_ref_in : opt.filename_seq_in;

                KmerStream_Build_opt kms_opt;

                kms_opt.threads = opt.nb_threads;
                kms_opt.verbose = opt.verbose;
                kms_opt.k = k_;
                kms_opt.g = g_;
                kms_opt.q = 0;

                for (const auto& s : v_files) kms_opt.files.push_back(s);

                KmerStream kms(kms_opt);

                nb_unique_kmers = max(1UL, kms.KmerF0());
                nb_unique_minimizers = max(1UL, kms.MinimizerF0());
                nb_non_unique_kmers = reference_mode ? 0 : max(1UL, nb_unique_kmers - min(nb_unique_kmers, kms.Kmerf1()));
                nb_non_unique_minimizers = reference_mode ? 0 : max(1UL, nb_unique_minimizers - min(nb_unique_minimizers, kms.Minimizerf1()));

                if (opt.verbose){

                    cout << "CompactedDBG::build(): Estimated number of k-mers occurring at least once: " << nb_unique_kmers << endl;
                    cout << "CompactedDBG::build(): Estimated number of minimizer occurring at least once: " << nb_unique_minimizers << endl;

                    if (!reference_mode){

                        cout << "CompactedDBG::build(): Estimated number of k-mers occurring twice or more: " << nb_non_unique_kmers << endl;
                        cout << "CompactedDBG::build(): Estimated number of minimizers occurring twice or more: " << nb_non_unique_minimizers << endl;
                    }
                }
            }

            setFullCoverage(reference_mode ? 1 : 2);

            if (opt.inFilenameBBF.length() != 0){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "CompactedDBG::build(): Could not open input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    construct_finished = false;
                }
                else {

                    construct_finished = bf.ReadBloomFilter(fp);

                    fclose(fp);
                }
            }
            else construct_finished = filter(opt, nb_unique_kmers, nb_non_unique_kmers);

            if (construct_finished){

                if (opt.outFilenameBBF.length() != 0){

                    FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

                    if (fp == NULL) {

                        cerr << "CompactedDBG::build(): Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                        construct_finished = false;
                    }
                    else {

                        bf.WriteBloomFilter(fp);

                        fclose(fp);
                    }
                }

                if (construct_finished) construct_finished = construct(opt, nb_unique_minimizers, nb_non_unique_minimizers); // Construction step

                bf.clear();
            }
        }
    }

    return construct_finished;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::simplify(const bool delete_short_isolated_unitigs, const bool clip_short_tips, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::simplify(): Graph is invalid and cannot be simplified" << endl;
        return false;
    }

    if (delete_short_isolated_unitigs || clip_short_tips){

        if (verbose) cout << endl << "CompactedDBG::simplify(): Removing isolated unitigs and/or clipping tips" << endl;

        vector<Kmer> v_joins;
        size_t joined = 0;

        size_t removed = removeUnitigs(delete_short_isolated_unitigs, clip_short_tips, v_joins);

        if (clip_short_tips) joined = joinUnitigs_<is_void<U>::value>(&v_joins);

        v_joins.clear();

        if (verbose){

            cout << "CompactedDBG::simplify(): After: " << size() << " unitigs" << endl;
            cout << "CompactedDBG::simplify(): Removed " << removed << " unitigs" << endl;
            cout << "CompactedDBG::simplify(): Joined " << joined << " unitigs" << endl;
        }

        return true;
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::write(const string& output_filename, const size_t nb_threads, const bool GFA_output, const bool verbose) const {

    if (invalid){

        cerr << "CompactedDBG::write(): Graph is invalid and cannot be written to disk" << endl;
        return false;
    }

    if (nb_threads <= 0){

        cerr << "CompactedDBG::write(): Number of threads cannot be less than 0" << endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        cerr << "CompactedDBG::write(): Number of threads cannot exceed " << std::thread::hardware_concurrency() << "threads" << endl;
        return false;
    }

    if (verbose) cout << endl << "CompactedDBG::write(): Writing graph to disk" << endl;

    const string out = output_filename + (GFA_output ? ".gfa" : ".fasta");

    FILE* fp = fopen(out.c_str(), "w");

    if (fp == NULL) {

        cerr << "CompactedDBG::write(): Could not open file " << out << " for writing graph" << endl;
        return false;
    }
    else {

        fclose(fp);
        if (std::remove(out.c_str()) != 0) cerr << "CompactedDBG::write(): Could not remove temporary file " << out << endl;
    }

    GFA_output ? writeGFA(out, nb_threads) : writeFASTA(out);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::read(const string& input_filename, const size_t nb_threads, const bool verbose){

    if (verbose) cout << endl << "CompactedDBG::read(): Reading graph from disk" << endl;

    const int format = FileParser::getFileFormat(input_filename.c_str());

    if (format == -1){

        cerr << "CompactedDBG::read(): Input file " << input_filename << " does not exist, is ill-formed or is not in FASTA/FASTQ/GFA format." << endl;

        return false;
    }
    else if ((format != 0) && (format != 2)){

        cerr << "CompactedDBG::read(): Input file must be in FASTA or GFA format." << endl;

        return false;
    }

    if (format == 2){ // GFA format

        FILE* fp = fopen(input_filename.c_str(), "r");

        if (fp == NULL) {

            cerr << "CompactedDBG::read(): Could not open file " << input_filename << " for reading graph" << endl;
            return false;
        }

        fclose(fp);

        char buffer[4096];

        ifstream graphfile_in(input_filename);
        istream graph_in(graphfile_in.rdbuf());

        graph_in.getline(buffer, 4096); // Read and discard header
        graphfile_in.close();

        const string header(buffer);

        if (header[0] != 'H'){

            cerr << "CompactedDBG::read(): An error occurred while reading input GFA file." << endl;
            return false;
        }

        stringstream hs(&buffer[2]); // Skip the first 2 char. of the line "H\t"
        string sub;

        int k = k_;
        int g = g_;

        while (hs.good()){ // Split line based on tabulation

            getline(hs, sub, '\t');

            const string tag = sub.substr(0, 5);

            if (tag == "KL:Z:") k = atoi(sub.c_str() + 5);
            else if (tag == "ML:Z:") g = atoi(sub.c_str() + 5);
        }

        clear();

        {
            KmerStream_Build_opt kms_opt;

            kms_opt.threads = nb_threads;
            kms_opt.verbose = verbose;
            kms_opt.k = k;
            kms_opt.g = g;
            kms_opt.q = 0;

            kms_opt.files.push_back(input_filename);

            KmerStream kms(kms_opt);

            MinimizerIndex hmap_min_unitigs_tmp(max(1UL, kms.MinimizerF0()) * 1.05);
            hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
        }

        setKmerGmerLength(k, g);

        if (!invalid) readGFA(input_filename, nb_threads);

        if (verbose) cout << endl << "CompactedDBG::read(): Finished reading graph from disk" << endl;

        return !invalid;
    }
    else {

        const int k = k_;
        const int g = g_;

        clear();

        {
            KmerStream_Build_opt kms_opt;

            kms_opt.threads = nb_threads;
            kms_opt.verbose = verbose;
            kms_opt.k = k;
            kms_opt.g = g;
            kms_opt.q = 0;

            kms_opt.files.push_back(input_filename);

            KmerStream kms(kms_opt);

            MinimizerIndex hmap_min_unitigs_tmp(max(1UL, kms.MinimizerF0()) * 1.05);

            hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
        }

        setKmerGmerLength(k, g);

        readFASTA(input_filename, nb_threads);
    }

    setFullCoverage(1);

    for (auto& unitig : *this) unitig.setFullCoverage();

    invalid = false;

    if (verbose) cout << endl << "CompactedDBG::read(): Finished reading graph from disk" << endl;

    return true;
}

template<typename U, typename G>
const_UnitigMap<U, G> CompactedDBG<U, G>::find(const char* s, const size_t pos_km, const minHashIterator<RepHash>& it_min, const bool extremities_only) const {

    if (invalid){

        cerr << "CompactedDBG::find(): Graph is invalid and cannot be searched" << endl;
        return const_UnitigMap<U, G>();
    }

    const Kmer km(s + pos_km);
    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        const size_t min_h_res_pos_km = min_h_res.pos - pos_km;
        Minimizer minz(Minimizer(s + min_h_res.pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return const_UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(s + mhr.pos).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if (((min_h_res_pos_km == unitig_id_pos) || (min_h_res_pos_km == diff - unitig_id_pos)) && (km_unitigs.getKmer(unitig_id) == km_rep)){

                            return const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km == km_rep, this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res_pos_km;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res_pos_km;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res_pos_km;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return const_UnitigMap<U, G>();
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::find(const char* s, const size_t pos_km, const minHashIterator<RepHash>& it_min, const bool extremities_only) {

    if (invalid){

        cerr << "CompactedDBG::find(): Graph is invalid and cannot be searched" << endl;
        return UnitigMap<U, G>();
    }

    const Kmer km(s + pos_km);
    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        const size_t min_h_res_pos_km = min_h_res.pos - pos_km;
        Minimizer minz(Minimizer(s + min_h_res.pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(s + mhr.pos).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if (((min_h_res_pos_km == unitig_id_pos) || (min_h_res_pos_km == diff - unitig_id_pos)) && (km_unitigs.getKmer(unitig_id) == km_rep)){

                            return UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km == km_rep, this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res_pos_km;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res_pos_km;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res_pos_km;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return UnitigMap<U, G>();
}

template<typename U, typename G>
const_UnitigMap<U, G> CompactedDBG<U, G>::find(const Kmer& km, const bool extremities_only) const {

    if (invalid){

        cerr << "CompactedDBG::find(): Graph is invalid and cannot be searched" << endl;
        return const_UnitigMap<U, G>();
    }

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km.toString(km_tmp); // Set k-mer to look-up in string version

    minHashKmer<RepHash> it_min(km_tmp, k_, g_, RepHash(), true), it_min2, it_min_end;

    while (it_min != it_min_end){

        int mhr_pos = it_min.getPosition();
        Minimizer minz(Minimizer(km_tmp + mhr_pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        it_min2 = it_min;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return const_UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        it_min2.getNewMin();

                        if (it_min2 != it_min_end){

                            minz = Minimizer(km_tmp + it_min2.getPosition()).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if (((mhr_pos == unitig_id_pos) || (mhr_pos == diff - unitig_id_pos)) && (km_unitigs.getKmer(unitig_id) == km_rep)){

                            return const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km == km_rep, this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - mhr_pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + mhr_pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + mhr_pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        ++it_min;
    }

    return const_UnitigMap<U, G>();
}

/*template<typename U, typename G>
vector<const_UnitigMap<U, G>> CompactedDBG<U, G>::find(const Minimizer& minz) const {

    if (invalid){

        cerr << "CompactedDBG::find(): Graph is invalid and cannot be searched" << endl;
        return vector<const_UnitigMap<U, G>>();
    }

    vector<const_UnitigMap<U, G>> v_um;

    const Minimizer minz_rep = minz.rep();
    const MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz_rep);

    if (it != hmap_min_unitigs.end()){ // If the minimizer is found

        const packed_tiny_vector& v = it.getVector();
        const uint8_t flag_v = it.getVectorSize();
        const int v_sz = v.size(flag_v);

        for (size_t i = 0; i < v_sz; ++i){

            size_t unitig_id_pos = v(i, flag_v);
            size_t unitig_id = unitig_id_pos >> 32;

            if (unitig_id == RESERVED_ID){

                if ((unitig_id_pos & RESERVED_ID) != 0) v_um.push_back(const_UnitigMap<U, G>()); // This minimizer has abundant k-mers
            }
            else {

                unitig_id_pos &= MASK_CONTIG_POS;

                if ((unitig_id_pos & MASK_CONTIG_TYPE) != 0) v_um.push_back(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, true, this));
                else {

                    const size_t len_unitig = v_unitigs[unitig_id]->length();

                    const size_t pos_s = (unitig_id_pos + g_ < k_) ? 0 : (unitig_id_pos + g_ - k_);
                    const size_t pos_e = (unitig_id_pos + k_ > len_unitig) ? len_unitig : (unitig_id_pos + k_);

                    const string subseq = v_unitigs[unitig_id]->getSeq().toString(pos_s, pos_e - pos_s);

                    const char* subseq_str = subseq.c_str();

                    size_t l_pos_s = 0xffffffffffffffffULL;
                    size_t l_pos_e = 0;

                    minHashIterator<RepHash> mhi_s = minHashIterator<RepHash>(subseq_str, pos_e - pos_s, k_, g_, RepHash(), true), mhi_e;
                    minHashResultIterator<RepHash> mhrit_s, mhrit_e;

                    while (mhi_s != mhi_e) {

                        mhrit_s = *mhi_s;

                        while (mhrit_s != mhrit_e) {

                            if (Minimizer(subseq_str + mhrit_s->pos).rep() == minz_rep) {

                                l_pos_s = min(l_pos_s, static_cast<size_t>(mhi_s.getKmerPosition()));
                                l_pos_e = max(l_pos_e, static_cast<size_t>(mhi_s.getKmerPosition()));

                                break;
                            }

                            ++mhrit_s;
                        }

                        ++mhi_s;
                    }

                    if (l_pos_e >= l_pos_s) v_um.push_back(const_UnitigMap<U, G>(unitig_id, pos_s + l_pos_s, l_pos_e - l_pos_s + 1, len_unitig + k_, false, false, true, this));
                }
            }
        }
    }

    return v_um;
}*/

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::find(const Kmer& km, const bool extremities_only) {

    if (invalid){

        cerr << "CompactedDBG::find(): Graph is invalid and cannot be searched" << endl;
        return UnitigMap<U, G>();
    }

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km.toString(km_tmp); // Set k-mer to look-up in string version

    minHashKmer<RepHash> it_min(km_tmp, k_, g_, RepHash(), true), it_min2, it_min_end;

    while (it_min != it_min_end){

        int mhr_pos = it_min.getPosition();

        Minimizer minz(Minimizer(km_tmp + mhr_pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        it_min2 = it_min;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()) return UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        it_min2.getNewMin();

                        if (it_min2 != it_min_end){

                            minz = Minimizer(km_tmp + it_min2.getPosition()).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if (((mhr_pos == unitig_id_pos) || (mhr_pos == diff - unitig_id_pos)) && (km_unitigs.getKmer(unitig_id) == km_rep)){

                            return UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km == km_rep, this);
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - mhr_pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + mhr_pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + mhr_pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                                return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        ++it_min;
    }

    return UnitigMap<U, G>();
}

template<typename U, typename G>
vector<const_UnitigMap<U, G>> CompactedDBG<U, G>::findPredecessors(const Kmer& km, const bool extremities_only) const {

    const Kmer km_pred[4] = {km.backwardBase('A'), km.backwardBase('C'), km.backwardBase('G'), km.backwardBase('T')};
    const Kmer km_rep[4] = {km_pred[0].rep(), km_pred[1].rep(), km_pred[2].rep(), km_pred[3].rep()};

    const Kmer km_twin_a = km_pred[0].twin();

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km_pred[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    vector<const_UnitigMap<U, G>> v_um(4, const_UnitigMap<U, G>(1, this));

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(km_tmp + min_h_res.pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end()){

                                v_um[j].partialCopy(const_UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km_pred[j] == km_rep[j], this));
                            }
                        }
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(km_tmp + mhr.pos).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if ((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)){

                            const Kmer km_unitig = km_unitigs.getKmer(unitig_id);

                            uint8_t idx = bits[km_unitig.getChar(0)];

                            if (km_unitig == km_rep[idx]) {

                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_pred[idx] == km_rep[idx], this));
                            }
                            else {

                                idx = 0x3 - bits[km_unitig.getChar(k_ - 1)];

                                if (km_unitig == km_rep[idx]) {

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_pred[idx] == km_rep[idx], this));
                                }
                            }
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km)){

                                const uint8_t idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km)){

                                const uint8_t idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];
                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return v_um;
}

template<typename U, typename G>
vector<UnitigMap<U, G>> CompactedDBG<U, G>::findPredecessors(const Kmer& km, const bool extremities_only) {

    const Kmer km_pred[4] = {km.backwardBase('A'), km.backwardBase('C'), km.backwardBase('G'), km.backwardBase('T')};
    const Kmer km_rep[4] = {km_pred[0].rep(), km_pred[1].rep(), km_pred[2].rep(), km_pred[3].rep()};

    const Kmer km_twin_a = km_pred[0].twin();

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km_pred[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    vector<UnitigMap<U, G>> v_um(4, UnitigMap<U, G>(1, this));

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(km_tmp + min_h_res.pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end()){

                                v_um[j] = UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km_pred[j] == km_rep[j], this);
                            }
                        }
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(km_tmp + mhr.pos).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if ((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)){

                            const Kmer km_unitig = km_unitigs.getKmer(unitig_id);

                            uint8_t idx = bits[km_unitig.getChar(0)];

                            if (km_unitig == km_rep[idx]) {

                                v_um[idx] = UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_pred[idx] == km_rep[idx], this);
                            }
                            else {

                                idx = 0x3 - bits[km_unitig.getChar(k_ - 1)];

                                if (km_unitig == km_rep[idx]) {

                                    v_um[idx] = UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_pred[idx] == km_rep[idx], this);
                                }
                            }
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km)){

                                const uint8_t idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];
                                v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];
                                v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km)){

                                const uint8_t idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];
                                v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_twin_a)){

                                const uint8_t idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];
                                v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return v_um;
}

template<typename U, typename G>
vector<const_UnitigMap<U, G>> CompactedDBG<U, G>::findSuccessors(const Kmer& km, const size_t limit, const bool extremities_only) const {

    vector<const_UnitigMap<U, G>> v_um(4, const_UnitigMap<U, G>(1, this));

    if (limit == 0) return v_um;

    const Kmer km_succ[4] = {km.forwardBase('A'), km.forwardBase('C'), km.forwardBase('G'), km.forwardBase('T')};
    const Kmer km_rep[4] = {km_succ[0].rep(), km_succ[1].rep(), km_succ[2].rep(), km_succ[3].rep()};

    const Kmer km_twin_a = km_succ[0].twin().forwardBase('A');

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    int nb_found = 0;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km_succ[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(km_tmp + min_h_res.pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if (v_um[j].isEmpty && ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end())){

                                v_um[j].partialCopy(const_UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km_succ[j] == km_rep[j], this));

                                if (++nb_found == limit) return v_um;
                            }
                        }
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(km_tmp + mhr.pos).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if ((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)){

                            const Kmer km_unitig = km_unitigs.getKmer(unitig_id);

                            uint8_t idx = bits[km_unitig.getChar(k_ - 1)];

                            if (v_um[idx].isEmpty && (km_unitig == km_rep[idx])){

                                v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_succ[idx] == km_rep[idx], this));

                                if (++nb_found == limit) return v_um;
                            }
                            else {

                                idx = 0x3 - bits[km_unitig.getChar(0)];

                                if (v_um[idx].isEmpty && (km_unitig == km_rep[idx])){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_succ[idx] == km_rep[idx], this));

                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_succ[0])){

                                const int idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));

                                    if (++nb_found == limit) return v_um;
                                }
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));

                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_succ[0])){

                                const int idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this));

                                    if (++nb_found == limit) return v_um;
                                }

                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx].partialCopy(const_UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this));

                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return v_um;
}

template<typename U, typename G>
vector<UnitigMap<U, G>> CompactedDBG<U, G>::findSuccessors(const Kmer& km, const size_t limit, const bool extremities_only) {

    vector<UnitigMap<U, G>> v_um(4, UnitigMap<U, G>(1, this));

    if (limit == 0) return v_um;

    const Kmer km_succ[4] = {km.forwardBase('A'), km.forwardBase('C'), km.forwardBase('G'), km.forwardBase('T')};
    const Kmer km_rep[4] = {km_succ[0].rep(), km_succ[1].rep(), km_succ[2].rep(), km_succ[3].rep()};

    const Kmer km_twin_a = km_succ[0].twin().forwardBase('A');

    bool isShort;

    size_t unitig_id_pos, unitig_id, len;

    int64_t pos_match;

    int nb_found = 0;

    const int diff = k_ - g_;

    char km_tmp[MAX_KMER_SIZE];
    km_succ[0].toString(km_tmp); // Set k-mer to look-up in string version

    preAllocMinHashIterator<RepHash> it_min(km_tmp, k_, k_, g_, RepHash(), true);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz(Minimizer(km_tmp + min_h_res.pos).rep());
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id_pos = v(i, flag_v);
                unitig_id = unitig_id_pos >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((unitig_id_pos & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km;

                        for (size_t j = 0; j != 4; ++j){

                            if (v_um[j].isEmpty && ((it_km = h_kmers_ccov.find(km_rep[j])) != h_kmers_ccov.end())){

                                v_um[j] = UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km_succ[j] == km_rep[j], this);

                                if (++nb_found == limit) return v_um;
                            }
                        }
                    }

                    if ((unitig_id_pos & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(km_tmp + mhr.pos).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (unitig_id_pos & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos &= MASK_CONTIG_POS;

                    if (isShort){

                        if ((min_h_res.pos == unitig_id_pos) || (min_h_res.pos == diff - unitig_id_pos)){

                            const Kmer km_unitig = km_unitigs.getKmer(unitig_id);

                            uint8_t idx = bits[km_unitig.getChar(k_ - 1)];

                            if (v_um[idx].isEmpty && (km_unitig == km_rep[idx])){

                                v_um[idx] = UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_succ[idx] == km_rep[idx], this);

                                if (++nb_found == limit) return v_um;
                            }
                            else {

                                idx = 0x3 - bits[km_unitig.getChar(0)];

                                if (v_um[idx].isEmpty && (km_unitig == km_rep[idx])){

                                    v_um[idx] = UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, km_succ[idx] == km_rep[idx], this);

                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                    }
                    else {

                        len = v_unitigs[unitig_id]->length() - k_;
                        pos_match = unitig_id_pos - min_h_res.pos;

                        if (extremities_only){

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_succ[0])){

                                const int idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);

                                    if (++nb_found == limit) return v_um;
                                }
                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if (((pos_match == 0) || (pos_match == len)) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);

                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                        else{

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_ - 1, km_succ[0])){

                                const int idx = bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match + k_ - 1)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);

                                    if (++nb_found == limit) return v_um;
                                }

                            }

                            pos_match = unitig_id_pos - diff + min_h_res.pos;

                            if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match + 1, k_ - 1, km_twin_a)){

                                const int idx = 0x3 - bits[v_unitigs[unitig_id]->getSeq().getChar(pos_match)];

                                if (v_um[idx].isEmpty){

                                    v_um[idx] = UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);

                                    if (++nb_found == limit) return v_um;
                                }
                            }
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return v_um;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::add(const string& seq, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::add(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::add(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    const char* str_seq = seq.c_str();

    string no_dup_km_seq;

    unordered_set<Kmer, KmerHash> km_seen;

    size_t last_start_pos = 0;

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end; ++it_km) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        if (!km_seen.insert(p.first.rep()).second){

            no_dup_km_seq = seq.substr(last_start_pos, p.second - last_start_pos + k_ - 1);
            last_start_pos = p.second;

            std::transform(no_dup_km_seq.begin(), no_dup_km_seq.end(), no_dup_km_seq.begin(), ::toupper);

            mergeUnitig(no_dup_km_seq, verbose);
            km_seen.clear();
        }
    }

    no_dup_km_seq = seq.substr(last_start_pos, seq.length() - last_start_pos + k_ - 1);

    std::transform(no_dup_km_seq.begin(), no_dup_km_seq.end(), no_dup_km_seq.begin(), ::toupper);

    mergeUnitig(no_dup_km_seq, verbose);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::remove(const const_UnitigMap<U, G>& um, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::remove(): Graph is invalid and no unitig can be removed from it" << endl;
        return false;
    }

    if (um.isEmpty) return false;

    vector<Kmer> v_km;

    const Kmer head(um.getUnitigHead());
    const Kmer tail(um.getUnitigTail());

    for (size_t i = 0; i != 4; ++i){

        const Kmer bw(head.backwardBase(alpha[i]));

        if (!find(bw, true).isEmpty) v_km.push_back(bw);
    }

    for (size_t i = 0; i != 4; ++i){

        const Kmer fw(tail.forwardBase(alpha[i]));

        if (!find(fw, true).isEmpty) v_km.push_back(fw);
    }

    if (um.isAbundant){

        deleteUnitig_<is_void<U>::value>(um.isShort, um.isAbundant, um.pos_unitig);
    }
    else {

        const size_t swap_position = (um.isShort ? km_unitigs.size() : v_unitigs.size()) - 1;

        if (um.pos_unitig != swap_position) swapUnitigs(um.isShort, um.pos_unitig, swap_position);

        deleteUnitig_<is_void<U>::value>(um.isShort, um.isAbundant, swap_position);

        if (um.isShort) km_unitigs.resize(swap_position);
        else v_unitigs.resize(swap_position);
    }

    joinUnitigs_<is_void<U>::value>(&v_km);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::merge(const CompactedDBG& o, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "CompactedDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    if (o.invalid){

         if (verbose) cerr << "CompactedDBG::merge(): Graph to merge is invalid." << endl;
         ret = false;
    }

    if (k_ != o.k_){

         if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
         ret = false;
    }

    if (g_ != o.g_){

         if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same minimizer length." << endl;
         ret = false;
    }

    if (this == &o){

         if (verbose) cerr << "CompactedDBG::merge(): Cannot merge graph with itself." << endl;
         ret = false;
    }

    if (ret){

        const size_t sz_before = size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        if (annotateSplitUnitigs(o, nb_threads, verbose)){

            const size_t sz_after = size();
            const pair<size_t, size_t> p = splitAllUnitigs();
            const size_t joined = (p.second != 0) ? joinUnitigs_<is_void<U>::value>() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p.first << " unitigs into " << p.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << size() << " unitigs after merging." << endl;
            }

            if (!is_void<U>::value) mergeData(o, nb_threads, verbose);

            return true;
        }
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::merge(const vector<CompactedDBG>& v, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) cerr << "CompactedDBG::merge(): Current graph is invalid." << endl;
         ret = false;
    }

    for (const auto& cdbg : v){

        if (cdbg.invalid){

             if (verbose) cerr << "CompactedDBG::merge(): One of the graph to merge is invalid." << endl;
             ret = false;
        }

        if (k_ != cdbg.k_){

             if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same k-mer length." << endl;
             ret = false;
        }

        if (g_ != cdbg.g_){

             if (verbose) cerr << "CompactedDBG::merge(): The graphs to merge do not have the same minimizer length." << endl;
             ret = false;
        }

        if (this == &cdbg){

             if (verbose) cerr << "CompactedDBG::merge(): Cannot merge graph with itself." << endl;
             ret = false;
        }
    }

    if (ret){

        const size_t sz_before = size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        for (const auto& cdbg : v){

            ret = annotateSplitUnitigs(cdbg, nb_threads, verbose);

            if (!ret) break;
        }

        if (ret){

            const size_t sz_after = size();
            const pair<size_t, size_t> p = splitAllUnitigs();
            const size_t joined = (p.second != 0) ? joinUnitigs_<is_void<U>::value>() : 0;

            if (verbose){

                cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Split " << p.first << " unitigs into " << p.second << " new unitigs." << endl;
                cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << endl;
                cout << "CompactedDBG::merge(): " << size() << " unitigs after merging." << endl;
            }

            if (!is_void<U>::value){

                for (const auto& cdbg : v) mergeData(cdbg, nb_threads, verbose);
            }

            return true;
        }
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::annotateSplitUnitigs(const CompactedDBG<U, G>& o, const size_t nb_threads, const bool verbose){

    if ((this != &o) && !invalid && !o.invalid) { // TODO: Check if k_ and g_ are the same

        if (verbose){

            cout << "CompactedDBG::annotateSplitUnitigs(): Current graph has " << size() << " unitigs." << endl;
            cout << "CompactedDBG::annotateSplitUnitigs(): Graph to merge has " << o.size() << " unitigs." << endl;
            cout << "CompactedDBG::annotateSplitUnitigs(): Start unitigs merging." << endl;
        }

        if (nb_threads == 1){

            for (const auto& unitig : o) annotateSplitUnitig(unitig.referenceUnitigToString(), false);
        }
        else {

            const size_t chunk = 100;
            const size_t nb_locks = nb_threads * 1024;

            vector<thread> workers; // need to keep track of threads so we can join them

            typename CompactedDBG<U, G>::const_iterator g_a(o.begin());
            typename CompactedDBG<U, G>::const_iterator g_b(o.end());

            LockGraph lck_g(nb_locks);

            mutex mutex_o_unitig;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&/*, t*/]{

                        typename CompactedDBG<U, G>::const_iterator l_a, l_b;

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_o_unitig);

                                if (g_a == g_b) return;

                                l_a = g_a;
                                l_b = g_a;

                                for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                                g_a = l_b;
                            }

                            for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                                annotateSplitUnitig(it_unitig->referenceUnitigToString(), lck_g/*, t*/, false);
                            }
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        if (verbose) cout << "CompactedDBG::annotateSplitUnitigs(): Merging unitigs finished." << endl;

        return true;
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::mergeData(const CompactedDBG<U, G>& o, const size_t nb_threads, const bool verbose){

    if ((this != &o) && !invalid && !o.invalid){ // TODO: Check if k_ and g_ are the same

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data started." << endl;

        const size_t nb_locks = nb_threads * 1024;

        std::atomic_flag* locks_unitig = new std::atomic_flag[nb_locks];

        for (size_t i = 0; i < nb_locks; ++i) locks_unitig[i].clear();

        auto worker_function = [&](typename CompactedDBG<U, G>::const_iterator it_a, typename CompactedDBG<U, G>::const_iterator it_b) {

            while (it_a != it_b) {

                const string str(it_a->referenceUnitigToString());
                const char* str_seq = str.c_str();

                for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

                    const std::pair<Kmer, int>& p = *it_km;

                    const UnitigMap<U, G> um_dest = findUnitig(p.first, str_seq, p.second);

                    if (um_dest.isEmpty) ++it_km;
                    else {

                        const_UnitigMap<U, G> um_src(*it_a);

                        um_src.dist = p.second;
                        um_src.len = um_dest.len;

                        const uint64_t id_lock = um_dest.pos_unitig % nb_locks;

                        while (locks_unitig[id_lock].test_and_set(std::memory_order_acquire));

                        mergeData_<is_void<U>::value>(um_dest, um_src);

                        locks_unitig[id_lock].clear(std::memory_order_release);

                        it_km += um_dest.len;
                    }
                }

                ++it_a;
            }
        };

        const size_t chunk = 100;

        vector<thread> workers; // need to keep track of threads so we can join them

        typename CompactedDBG<U, G>::const_iterator g_a = o.begin();
        typename CompactedDBG<U, G>::const_iterator g_b = o.end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    typename CompactedDBG<U, G>::const_iterator l_a, l_b;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        worker_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data finished." << endl;

        delete[] locks_unitig;

        return true;
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::mergeData(CompactedDBG<U, G>&& o, const size_t nb_threads, const bool verbose){

    if ((this != &o) && !invalid && !o.invalid){ // TODO: Check if k_ and g_ are the same

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data started." << endl;

        const size_t nb_locks = nb_threads * 1024;

        std::atomic_flag* locks_unitig = new std::atomic_flag[nb_locks];

        for (size_t i = 0; i < nb_locks; ++i) locks_unitig[i].clear();

        auto worker_function = [&](typename CompactedDBG<U, G>::iterator it_a, typename CompactedDBG<U, G>::iterator it_b) {

            while (it_a != it_b) {

                const string str(it_a->referenceUnitigToString());
                const char* str_seq = str.c_str();

                for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

                    const std::pair<Kmer, int>& p = *it_km;

                    const UnitigMap<U, G> um_dest = findUnitig(p.first, str_seq, p.second);

                    if (um_dest.isEmpty) ++it_km;
                    else {

                        const_UnitigMap<U, G> um_src(*it_a);

                        um_src.dist = p.second;
                        um_src.len = um_dest.len;

                        const uint64_t id_lock = um_dest.pos_unitig % nb_locks;

                        while (locks_unitig[id_lock].test_and_set(std::memory_order_acquire));

                        mergeData_<is_void<U>::value>(um_dest, um_src);

                        locks_unitig[id_lock].clear(std::memory_order_release);

                        it_km += um_dest.len;
                    }
                }

                it_a->getData()->clear(*it_a);

                ++it_a;
            }
        };

        const size_t chunk = 100;

        vector<thread> workers; // need to keep track of threads so we can join them

        typename CompactedDBG<U, G>::iterator g_a = o.begin();
        typename CompactedDBG<U, G>::iterator g_b = o.end();

        mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename CompactedDBG<U, G>::iterator l_a, l_b;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        worker_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        if (verbose) cout << "CompactedDBG::mergeData(): Merging data finished." << endl;

        delete[] locks_unitig;

        return true;
    }

    return false;
}

template<typename U, typename G>
typename CompactedDBG<U, G>::iterator CompactedDBG<U, G>::begin() {

    if (invalid) return iterator();

    iterator it(this);
    ++it;

    return it;
}

template<typename U, typename G>
typename CompactedDBG<U, G>::const_iterator CompactedDBG<U, G>::begin() const {

    if (invalid) return const_iterator();

    const_iterator it(this);
    ++it;

    return it;
}

template<typename U, typename G>
typename CompactedDBG<U, G>::iterator CompactedDBG<U, G>::end() { return iterator(); }

template<typename U, typename G>
typename CompactedDBG<U, G>::const_iterator CompactedDBG<U, G>::end() const { return const_iterator(); }

template<typename U, typename G>
size_t CompactedDBG<U, G>::length() const {

    size_t len = 0;

    if (!invalid){

        len = (km_unitigs.size() + h_kmers_ccov.size()) * k_;

        for (const auto& unitig : v_unitigs) len += unitig->length();
    }

    return len;
}

template<typename U, typename G>
size_t CompactedDBG<U, G>::nbKmers() const {

    size_t nb = 0;

    if (!invalid){

        nb = (km_unitigs.size() + h_kmers_ccov.size());

        for (const auto& unitig : v_unitigs) nb += unitig->length() - k_ + 1;
    }

    return nb;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::filter(const CDBG_Build_opt& opt, const size_t nb_unique_kmers, const size_t nb_non_unique_kmers) {

    if (invalid){

        cerr << "CompactedDBG::filter(): Graph is invalid and it cannot be built" << endl;
        return false;
    }

    BlockedBloomFilter bf_tmp;

    const bool reference_mode = (opt.filename_ref_in.size() != 0);

    if (reference_mode){

        BlockedBloomFilter tmp(nb_unique_kmers, opt.nb_bits_unique_kmers_bf);
        bf = std::move(tmp);
    }
    else {

        {
            BlockedBloomFilter tmp(nb_unique_kmers, opt.nb_bits_unique_kmers_bf);
            bf_tmp = std::move(tmp);
        }

        {
            BlockedBloomFilter tmp(nb_non_unique_kmers, opt.nb_bits_non_unique_kmers_bf);
            bf = std::move(tmp);
        }
    }

    string s;

    size_t len_read = 0;
    size_t pos_read = 0;
    size_t nb_seq = 0;

    //const size_t max_len_seq = 1024;
    //const size_t thread_seq_buf_sz = 64 * max_len_seq;
    const size_t max_len_seq = rndup(static_cast<size_t>(1024 + k_ - 1));
    const size_t thread_seq_buf_sz = BUFFER_SIZE;

    const bool multi_threaded = (opt.nb_threads != 1);

    atomic<uint64_t> num_kmers(0), num_ins(0);

    const vector<string>& filename_in = reference_mode ? opt.filename_ref_in : opt.filename_seq_in;

    FileParser fp(filename_in);

    // Main worker thread
    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz) {

        uint64_t l_num_kmers = 0, l_num_ins = 0;

        char* str = seq_buf;
        const char* str_end = &seq_buf[seq_buf_sz];

        while (str < str_end) { // for each input

            const int len = strlen(str);

            for (char* s = str; s != str + len; ++s) *s &= 0xDF; // Put characters in upper case

            KmerHashIterator<RepHash> it_kmer_h(str, len, k_), it_kmer_h_end;
            minHashIterator<RepHash> it_min(str, len, k_, g_, RepHash(), true);

            if (reference_mode){

                for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++l_num_kmers) {

                    const pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                    it_min += (p_.second - it_min.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.
                    l_num_ins += bf.insert(p_.first, it_min.getHash(), multi_threaded);
                }
            }
            else {

                for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h, ++l_num_kmers) {

                    const pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                    it_min += (p_.second - it_min.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                    const uint64_t min_hr = it_min.getHash();

                    if (bf_tmp.insert(p_.first, min_hr, multi_threaded)) ++l_num_ins;
                    else bf.insert(p_.first, min_hr, multi_threaded);
                }
            }

            str += len + 1;
        }

        // atomic adds
        num_kmers += l_num_kmers;
        num_ins += l_num_ins;
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

        size_t file_id = 0;

        const size_t sz_buf = thread_seq_buf_sz - k_;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_buf) {

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                nb_seq += new_reading;

                //pos_read = (new_reading ? 0 : pos_read);
                pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k_){

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
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    char* buffer_seq = new char[thread_seq_buf_sz]();

                    size_t buffer_seq_sz = 0;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) {

                                delete[] buffer_seq;
                                return;
                            }

                            stop = reading_function(buffer_seq, buffer_seq_sz);
                        }

                        worker_function(buffer_seq, buffer_seq_sz);
                    }

                    delete[] buffer_seq;
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    fp.close();

    if (opt.verbose) {

        cout << "CompactedDBG::filter(): Closed all fasta/fastq files" << endl;
        cout << "CompactedDBG::filter(): Processed " << num_kmers << " k-mers in " << nb_seq  << " reads" << endl;
        cout << "CompactedDBG::filter(): Found " << num_ins << " unique k-mers" << endl;
        cout << "CompactedDBG::filter(): Number of blocks in Bloom filter is " << bf.getNbBlocks() << endl;
    }

    if (opt.useMercyKmers && !reference_mode){

        const string mbbf_uniq_filename(opt.prefixFilenameOut + "_uniq");

        FILE* f_mbbf = fopen(mbbf_uniq_filename.c_str(), "wb");

        if (f_mbbf == NULL){

            cerr << "CompactedDBG::filter(): Minimizer Blocked Bloom filter of unique k-mers cannot be written to disk" << endl;
            return false;
        }

        bf_tmp.WriteBloomFilter(f_mbbf);

        fclose(f_mbbf);
    }

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::construct(const CDBG_Build_opt& opt, const size_t nb_unique_minimizers, const size_t nb_non_unique_minimizers){

    if (invalid){

        cerr << "CompactedDBG::construct(): Graph is invalid and cannot be built" << endl;
        return false;
    }

    const bool reference_mode = (opt.filename_ref_in.size() != 0);

    const vector<string>& filename_in = reference_mode ? opt.filename_ref_in : opt.filename_seq_in;

    FileParser fp(filename_in);

    string s;

    size_t len_read = 0;
    size_t pos_read = 0;

    //const size_t max_len_seq = 1024;
    //const size_t thread_seq_buf_sz = 64 * max_len_seq;
    const size_t max_len_seq = rndup(static_cast<size_t>(1024 + k_ - 1));
    const size_t thread_seq_buf_sz = BUFFER_SIZE;

    tiny_vector<Kmer, 2>* fp_candidate = nullptr;

    KmerHashTable<bool> ignored_km_tips;

    const size_t nb_locks = opt.nb_threads * 1024;

    std::atomic_flag lock_ignored_km_tips = ATOMIC_FLAG_INIT;

    vector<SpinLock> locks_fp;

    LockGraph lck_g(nb_locks);

    if (!reference_mode){

        fp_candidate = new tiny_vector<Kmer, 2>[bf.getNbBlocks()];
        locks_fp = vector<SpinLock>(nb_locks);

        MinimizerIndex hmap_min_unitigs_tmp(nb_non_unique_minimizers * 1.05);

        hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
    }
    else {

        MinimizerIndex hmap_min_unitigs_tmp(nb_unique_minimizers * 1.05);

        hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
    }

    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz) {

        vector<Kmer> l_ignored_km_tips;

        Kmer km;
        RepHash rep;

        char* str = seq_buf;
        const char* str_end = &seq_buf[seq_buf_sz];

        while (str < str_end) { // for each input

            const int len = strlen(str);

            for (char* s = str; s != &str[len]; ++s) *s &= 0xDF;

            for (size_t i = 0; i < len - k_ + 1; i += max_len_seq - k_ + 1){

                const int curr_len = min(len - i, max_len_seq);
                const char saved_char = str[i + curr_len];
                const char* str_tmp = &str[i];

                str[i + curr_len] = '\0';

                KmerHashIterator<RepHash> it_kmer_h(str_tmp, curr_len, k_), it_kmer_h_end;
                minHashIterator<RepHash> it_min(str_tmp, curr_len, k_, g_, rep, true);

                for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

                    const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                    it_min += (p_.second - it_min.getKmerPosition());

                    const uint64_t it_min_h = it_min.getHash();
                    const int bid = bf.contains_bids(p_.first, it_min_h);

                    if (bid != -1){

                        km = Kmer(str_tmp + p_.second);

                        lck_g.acquire_reader();

                        const UnitigMap<U, G> um = findUnitig(km, str_tmp, p_.second);

                        if (um.isEmpty) { // kmer did not map, push into queue for next unitig generation round

                            lck_g.release_reader();

                            string newseq;

                            bool isIsolated = false;

                            const size_t pos_match = findUnitigSequenceBBF(km, newseq, isIsolated, l_ignored_km_tips); //Build unitig from Bloom filter

                            if (!reference_mode && isIsolated){ // According to the BF, k-mer is isolated in the graph and is a potential false positive

                                const uint64_t id_lock = bid % nb_locks;
                                const Kmer km_rep(km.rep());

                                tiny_vector<Kmer, 2>& v = fp_candidate[bid];

                                size_t i = 0;

                                locks_fp[id_lock].acquire();

                                for (; i < v.size(); ++i){ // Search list of fp candidate for k-mer

                                    if (v[i] == km_rep) break;
                                }

                                if (i >= v.size()){

                                    v.push_back(km_rep);

                                    locks_fp[id_lock].release();
                                }
                                else {

                                    v.remove(i);

                                    locks_fp[id_lock].release();

                                    addUnitigSequenceBBF(km, newseq, pos_match, 1, lck_g);
                                    addUnitigSequenceBBF(km, newseq, pos_match, 1, lck_g);
                                }
                            }
                            else {

                                const size_t len_match_km = 1 + cstrMatch(&str_tmp[p_.second + k_], &(newseq.c_str()[pos_match + k_]));

                                addUnitigSequenceBBF(km, newseq, pos_match, len_match_km, lck_g);

                                it_kmer_h += len_match_km - 1;
                            }
                        }
                        else {

                            mapRead(um, lck_g);

                            lck_g.release_reader();

                            it_kmer_h += um.len - 1;
                        }
                    }
                }

                str[i + curr_len] = saved_char;
            }

            str += len + 1;
        }

        while (lock_ignored_km_tips.test_and_set(std::memory_order_acquire));

        for (const auto& km_tip : l_ignored_km_tips) ignored_km_tips.insert(km_tip, false);

        lock_ignored_km_tips.clear(std::memory_order_release);
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

        size_t file_id = 0;

        const size_t sz_buf = thread_seq_buf_sz - k_;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_buf) {

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k_){

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
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        if (opt.verbose) cout << "CompactedDBG::construct(): Extract approximate unitigs" << endl;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    char* buffer_seq = new char[thread_seq_buf_sz];

                    size_t buffer_seq_sz = 0;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) {

                                delete[] buffer_seq;
                                return;
                            }

                            stop = reading_function(buffer_seq, buffer_seq_sz);
                        }

                        worker_function(buffer_seq, buffer_seq_sz);
                    }

                    delete[] buffer_seq;
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    fp.close();

    bf.clear();
    lck_g.clear();
    locks_fp.clear();

    if (fp_candidate != nullptr) delete[] fp_candidate;

    if (opt.verbose) cout << "CompactedDBG::construct(): Closed all input files" << endl;

    const size_t unitigsBefore = size();

    if (opt.verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (1/2)" << endl;

    pair<size_t, size_t> unitigSplit = extractAllUnitigs();

    const int unitigsAfter1 = size();

    if (opt.verbose) cout << endl << "CompactedDBG::construct(): Splitting unitigs (2/2)" << endl;

    check_fp_tips(ignored_km_tips);
    ignored_km_tips.clear_tables();

    const int unitigsAfter2 = size();

    if (opt.verbose) {

        cout << "CompactedDBG::construct(): Before split: " << unitigsBefore << " unitigs" << endl;
        cout << "CompactedDBG::construct(): After split (1/" << (reference_mode ? "1" : "2" ) << "): " << unitigsAfter1 << " unitigs" <<  endl;
        if (!reference_mode) cout << "CompactedDBG::construct(): After split (2/2): " << unitigsAfter2 << " unitigs" <<  endl;
        cout << "CompactedDBG::construct(): Unitigs split: " << unitigSplit.first << endl;
        cout << "CompactedDBG::construct(): Unitigs deleted: " << unitigSplit.second << endl;

        cout << endl << "CompactedDBG::construct(): Joining unitigs" << endl;
    }

    const size_t joined = joinUnitigs_<is_void<U>::value>(nullptr, opt.nb_threads);

    const int unitigsAfter3 = size();

    if (opt.verbose) {

        cout << "CompactedDBG::construct(): After join: " << unitigsAfter3 << " unitigs" << endl;
        cout << "CompactedDBG::construct(): Joined " << joined << " unitigs" << endl;
    }

    if (opt.useMercyKmers && !reference_mode){

        string filename_mbbf_uniq_km = opt.prefixFilenameOut + "_uniq";

        joinTips(filename_mbbf_uniq_km, opt.nb_threads, opt.verbose);

        if (opt.verbose) cout << "CompactedDBG::construct(): After join tips using mercy k-mers: " << size() << " unitigs" << endl;

        if (std::remove(filename_mbbf_uniq_km.c_str()) != 0) {

            cerr << "CompactedDBG::construct(): Minimizer Blocked Bloom filter file of unique k-mers cannot be removed from disk" << endl;
        }
    }

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::addUnitigSequenceBBF(const Kmer km, const string& seq, const size_t pos_match_km, const size_t len_match_km, LockGraph& lck_g) {

    bool isAbundant = false;

    const bool isShort = (seq.length() == k_);

    lck_g.acquire_writer();

    const size_t id_unitig = isShort ? km_unitigs.size() : v_unitigs.size();

    UnitigMap<U, G> um = find(km); // Look if unitig was already inserted

    if (um.isEmpty) isAbundant = addUnitig(seq, id_unitig); // If it wasn't already inserted, does it

    lck_g.release_writer();
    lck_g.acquire_reader();

    um = find(km);
    um.len = len_match_km; // Prepare read mapping

    if (!um.isShort && !um.isAbundant && !um.strand) um.dist -= um.len - 1;

    if (um.dist + um.len > um.size - k_){ // This is a self loop

        KmerIterator it(&(seq.c_str()[pos_match_km]));

        for (size_t i = 0; i != len_match_km; ++it, ++i){

            um = find(it->first);

            mapRead(um, lck_g);
        }
    }
    else mapRead(um, lck_g);

    lck_g.release_reader();

    return !um.isEmpty;
}

// use:  cm.findUnitigSequenceBBF(km, s, selfLoop)
// pre:  km is in the bloom filter
// post: s is the unitig containing the kmer km
//       and the first k-mer in s is smaller (wrt. < operator)
//       than the last kmer
//       selfLoop is true of the unitig is a loop or hairpin
template<typename U, typename G>
size_t CompactedDBG<U, G>::findUnitigSequenceBBF(Kmer km, string& s, bool& isIsolated, vector<Kmer>& l_ignored_km_tip) {

    string fw_s, bw_s;

    Kmer end(km);
    Kmer last(end);

    const Kmer twin(km.twin());

    char c;

    size_t j = 0;

    bool has_no_neighbor;
    bool self_loop = false;

    isIsolated = false;

    while (fwStepBBF(end, end, c, has_no_neighbor, l_ignored_km_tip)) {

        ++j;

        self_loop = (end == km);

        if (self_loop || (end == twin) || (end == last.twin())) break;

        fw_s.push_back(c);
        last = end;
    }

    if (!self_loop) {

        Kmer front(km);
        Kmer first(front);

        isIsolated = (j == 0) && has_no_neighbor;
        j = 0;

        while (bwStepBBF(front, front, c, has_no_neighbor, l_ignored_km_tip)) {

            ++j;

            if ((front == km) || (front == twin) || (front == first.twin())) break;

            bw_s.push_back(c);
            first = front;
        }

        isIsolated = isIsolated && (j == 0) && has_no_neighbor;

        reverse(bw_s.begin(), bw_s.end());
    }

    char tmp[MAX_KMER_SIZE];

    km.toString(tmp);

    s.reserve(k_ + fw_s.size() + bw_s.size());
    s.append(bw_s);
    s.append(tmp);
    s.append(fw_s);

    return bw_s.size();
}

template<typename U, typename G>
bool CompactedDBG<U, G>::bwStepBBF(const Kmer km, Kmer& front, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, const bool check_fp_cand) const {

    char km_tmp[MAX_KMER_SIZE];

    int found_fp_bw = 0;

    bool neigh_bw[4] = {false, false, false, false};
    uint64_t hashes_bw[4];

    RepHash rep_h(k_), rep_h_cpy;

    front.toString(km_tmp);
    rep_h.init(km_tmp);

    rep_h_cpy = rep_h;
    rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[0]);
    hashes_bw[0] = rep_h_cpy.hash();

    rep_h_cpy = rep_h;
    rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[1]);
    hashes_bw[1] = rep_h_cpy.hash();

    rep_h_cpy = rep_h;
    rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[2]);
    hashes_bw[2] = rep_h_cpy.hash();

    rep_h_cpy = rep_h;
    rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[3]);
    hashes_bw[3] = rep_h_cpy.hash();

    std::memmove(km_tmp + 1, km_tmp, (k_ - 1) * sizeof(char));

    const uint64_t it_min_h_bw = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

    size_t nb_neigh = bf.contains(hashes_bw, it_min_h_bw, neigh_bw, 2 * (static_cast<size_t>(check_fp_cand) + 1));
    size_t j = static_cast<size_t>(neigh_bw[1]) + (2ULL & (static_cast<size_t>(!neigh_bw[2]) - 1)) + (3ULL & (static_cast<size_t>(!neigh_bw[3]) - 1));

    if (check_fp_cand && (nb_neigh >= 2)){

        size_t j_tmp = 0;

        for (size_t i = 0; i < 4; ++i) {

            if (neigh_bw[i]){

                char dummy;
                bool has_no_neighbor_tmp = false;

                Kmer km_fp(front.backwardBase(alpha[i]));

                bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                neigh_bw[i] = has_no_neighbor_tmp && fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);
                found_fp_bw += neigh_bw[i];
                j_tmp += (i - j_tmp) & (static_cast<size_t>(neigh_bw[i]) - 1);
            }
        }

        if (found_fp_bw != 0){

            if ((nb_neigh - found_fp_bw) != 0){

                j = j_tmp;
                nb_neigh -= found_fp_bw;
            }
            else found_fp_bw = 0;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    if (check_fp_cand){

        int found_fp_fw = 0;

        bool neigh_fw[4] = {false, false, false, false};
        uint64_t hashes_fw[4];

        const Kmer bw(front.backwardBase(alpha[j]));

        bw.toString(km_tmp);
        rep_h.init(km_tmp);

        rep_h_cpy = rep_h;
        rep_h_cpy.updateFW(km_tmp[0], alpha[0]);
        hashes_fw[0] = rep_h_cpy.hash();

        rep_h_cpy = rep_h;
        rep_h_cpy.updateFW(km_tmp[0], alpha[1]);
        hashes_fw[1] = rep_h_cpy.hash();

        rep_h_cpy = rep_h;
        rep_h_cpy.updateFW(km_tmp[0], alpha[2]);
        hashes_fw[2] = rep_h_cpy.hash();

        rep_h_cpy = rep_h;
        rep_h_cpy.updateFW(km_tmp[0], alpha[3]);
        hashes_fw[3] = rep_h_cpy.hash();

        std::memmove(km_tmp, km_tmp + 1, (k_ - 1) * sizeof(char));

        const uint64_t it_min_h_fw = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

        nb_neigh = bf.contains(hashes_fw, it_min_h_fw, neigh_fw, 4);

        if (nb_neigh >= 2){

            for (size_t i = 0; i < 4; ++i) {

                if (neigh_fw[i]){

                    char dummy;
                    bool has_no_neighbor_tmp = false;

                    Kmer km_fp(bw.forwardBase(alpha[i]));

                    fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    neigh_fw[i] = has_no_neighbor_tmp && bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (neigh_fw[i] && (km_fp == km)){

                        found_fp_fw = 0;
                        break;
                    }
                    else found_fp_fw += (neigh_fw[i] && (km_fp != km));
                }
            }

            if (found_fp_fw){

                if (nb_neigh - found_fp_fw) nb_neigh -= found_fp_fw;
                else found_fp_fw = 0;
            }
        }

        if (nb_neigh != 1) return false;

        if (bw != km) {
            // exactly one k-mer in bw, character used is c

            for (size_t i = 0; (i < 4) && (found_fp_fw != 0); ++i) {

                if (neigh_fw[i]){

                    l_ignored_km_tip.push_back(bw.forwardBase(alpha[i]).rep());

                    --found_fp_fw;
                }
            }

             for (size_t i = 0; (i < 4) && (found_fp_bw != 0); ++i) {

                if (neigh_bw[i]){

                    l_ignored_km_tip.push_back(front.backwardBase(alpha[i]).rep());

                    --found_fp_bw;
                }
            }

            front = bw;
            c = alpha[j];

            return true;
        }

        return false;
    }

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::fwStepBBF(const Kmer km, Kmer& end, char& c, bool& has_no_neighbor, vector<Kmer>& l_ignored_km_tip, const bool check_fp_cand) const {

    char km_tmp[MAX_KMER_SIZE];

    int found_fp_fw = 0;

    bool neigh_fw[4] = {false, false, false, false};
    uint64_t hashes_fw[4];

    RepHash rep_h(k_), rep_h_cpy;

    end.toString(km_tmp);
    rep_h.init(km_tmp);

    rep_h_cpy = rep_h;
    rep_h_cpy.updateFW(km_tmp[0], alpha[0]);
    hashes_fw[0] = rep_h_cpy.hash();

    rep_h_cpy = rep_h;
    rep_h_cpy.updateFW(km_tmp[0], alpha[1]);
    hashes_fw[1] = rep_h_cpy.hash();

    rep_h_cpy = rep_h;
    rep_h_cpy.updateFW(km_tmp[0], alpha[2]);
    hashes_fw[2] = rep_h_cpy.hash();

    rep_h_cpy = rep_h;
    rep_h_cpy.updateFW(km_tmp[0], alpha[3]);
    hashes_fw[3] = rep_h_cpy.hash();

    std::memmove(km_tmp, km_tmp + 1, (k_ - 1) * sizeof(char));

    const uint64_t it_min_h_fw = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

    size_t nb_neigh = bf.contains(hashes_fw, it_min_h_fw, neigh_fw, 2 * (static_cast<size_t>(check_fp_cand) + 1));
    size_t j = static_cast<size_t>(neigh_fw[1]) + (2ULL & (static_cast<size_t>(!neigh_fw[2]) - 1)) + (3ULL & (static_cast<size_t>(!neigh_fw[3]) - 1));

    if (check_fp_cand && (nb_neigh >= 2)){

        size_t j_tmp = 0;

        for (size_t i = 0; i < 4; ++i) {

            if (neigh_fw[i]){

                char dummy;
                bool has_no_neighbor_tmp = false;

                Kmer km_fp(end.forwardBase(alpha[i]));

                fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                neigh_fw[i] = has_no_neighbor_tmp && bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                found_fp_fw += neigh_fw[i];
                j_tmp += (i - j_tmp) & (static_cast<size_t>(neigh_fw[i]) - 1);
            }
        }

        if (found_fp_fw){

            if (nb_neigh - found_fp_fw){

                j = j_tmp;
                nb_neigh -= found_fp_fw;
            }
            else found_fp_fw = 0;
        }
    }

    if (nb_neigh != 1) {

        has_no_neighbor = (nb_neigh == 0);
        return false;
    }
    else has_no_neighbor = false;

    if (check_fp_cand){

        int found_fp_bw = 0;

        bool neigh_bw[4] = {false, false, false, false};
        uint64_t hashes_bw[4];

        const Kmer fw(end.forwardBase(alpha[j]));

        fw.toString(km_tmp);
        rep_h.init(km_tmp);

        rep_h_cpy = rep_h;
        rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[0]);
        hashes_bw[0] = rep_h_cpy.hash();

        rep_h_cpy = rep_h;
        rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[1]);
        hashes_bw[1] = rep_h_cpy.hash();

        rep_h_cpy = rep_h;
        rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[2]);
        hashes_bw[2] = rep_h_cpy.hash();

        rep_h_cpy = rep_h;
        rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[3]);
        hashes_bw[3] = rep_h_cpy.hash();

        std::memmove(km_tmp + 1, km_tmp, (k_ - 1) * sizeof(char));

        const uint64_t it_min_h_bw = minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash();

        nb_neigh = bf.contains(hashes_bw, it_min_h_bw, neigh_bw, 4);

        if (nb_neigh >= 2){

            for (size_t i = 0; i < 4; ++i) {

                if (neigh_bw[i]){

                    char dummy;
                    bool has_no_neighbor_tmp = false;

                    Kmer km_fp(fw.backwardBase(alpha[i]));

                    bwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    neigh_bw[i] = has_no_neighbor_tmp && fwStepBBF(km_fp, km_fp, dummy, has_no_neighbor_tmp, l_ignored_km_tip, false);

                    if (neigh_bw[i] && (km_fp == km)){

                        found_fp_bw = 0;
                        break;
                    }
                    else found_fp_bw += (neigh_bw[i] && (km_fp != km));
                }
            }

            if (found_fp_bw){

                if (nb_neigh - found_fp_bw) nb_neigh -= found_fp_bw;
                else found_fp_bw = 0;
            }
        }

        if (nb_neigh != 1) return false;

        if (fw != km) {

            for (size_t i = 0; (i < 4) && (found_fp_bw != 0); ++i) {

                if (neigh_bw[i]){

                    l_ignored_km_tip.push_back(fw.backwardBase(alpha[i]).rep());

                    --found_fp_bw;
                }
            }

            for (size_t i = 0; (i < 4) && (found_fp_fw != 0); ++i) {

                if (neigh_fw[i]){

                    l_ignored_km_tip.push_back(end.forwardBase(alpha[i]).rep());

                    --found_fp_fw;
                }
            }

            end = fw;
            c = alpha[j];

            return true;
        }

        return false;
    }

    return true;
}

// use:  cc = cm.findUnitig(km,s,pos)
// pre:  s[pos,pos+k-1] is the kmer km
// post: cc contains either the reference to the unitig position
//       or empty if none found
template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const Kmer& km, const char* s, const size_t pos) {

    // need to check if we find it right away, need to treat this common case
    UnitigMap<U, G> um = find(km);

    if (!um.isEmpty && !um.isShort && !um.isAbundant){

        const size_t mask = static_cast<size_t>(um.strand) - 1;

        um.len = v_unitigs[um.pos_unitig]->getSeq().jump(s, pos, um.dist + ((k_ - 1) & mask), !um.strand) - k_ + 1;
        um.dist = um.dist - ((um.len - 1) & mask);
    }

    return um;
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const char* s, const size_t pos, const size_t len) {

    if ((len < k_) || (pos > len - k_)) return UnitigMap<U, G>();

    for (size_t i = pos; i != pos + k_; ++i){

        if (!isDNA(s[i])) return UnitigMap<U, G>();
    }

    // need to check if we find it right away, need to treat this common case
    UnitigMap<U, G> um = find(Kmer(s + pos));

    if (!um.isEmpty && !um.isShort && !um.isAbundant){

        const size_t mask = static_cast<size_t>(um.strand) - 1;

        um.len = v_unitigs[um.pos_unitig]->getSeq().jump(s, pos, um.dist + ((k_ - 1) & mask), !um.strand) - k_ + 1;
        um.dist = um.dist - ((um.len - 1) & mask);
    }

    return um;
}

template<typename U, typename G>
const_UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const char* s, const size_t pos, const size_t len) const {

    if ((len < k_) || (pos > len - k_)) return const_UnitigMap<U, G>();

    for (size_t i = pos; i != pos + k_; ++i){

        if (!isDNA(s[i])) return const_UnitigMap<U, G>();
    }

    // need to check if we find it right away, need to treat this common case
    const_UnitigMap<U, G> um = find(Kmer(s + pos));

    if (!um.isEmpty && !um.isShort && !um.isAbundant){

        const size_t mask = static_cast<size_t>(um.strand) - 1;

        um.len = v_unitigs[um.pos_unitig]->getSeq().jump(s, pos, um.dist + ((k_ - 1) & mask), !um.strand) - k_ + 1;
        um.dist = um.dist - ((um.len - 1) & mask);
    }

    return um;
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const char* s, const size_t pos, const size_t len, const minHashIterator<RepHash>& it_min) {

    if ((len < k_) || (pos > len - k_)) return UnitigMap<U, G>();

    for (size_t i = pos; i != pos + k_; ++i){

        if (!isDNA(s[i])) return UnitigMap<U, G>();
    }

    // need to check if we find it right away, need to treat this common case
    UnitigMap<U, G> um = find(s, pos, it_min);

    if (!um.isEmpty && !um.isShort && !um.isAbundant){

        const size_t mask = static_cast<size_t>(um.strand) - 1;

        um.len = v_unitigs[um.pos_unitig]->getSeq().jump(s, pos, um.dist + ((k_ - 1) & mask), !um.strand) - k_ + 1;
        um.dist = um.dist - ((um.len - 1) & mask);
    }

    return um;
}

template<typename U, typename G>
const_UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const char* s, const size_t pos, const size_t len, const minHashIterator<RepHash>& it_min) const {

    if ((len < k_) || (pos > len - k_)) return const_UnitigMap<U, G>();

    for (size_t i = pos; i != pos + k_; ++i){

        if (!isDNA(s[i])) return const_UnitigMap<U, G>();
    }

    // need to check if we find it right away, need to treat this common case
    const_UnitigMap<U, G> um = find(s, pos, it_min);

    if (!um.isEmpty && !um.isShort && !um.isAbundant){

        const size_t mask = static_cast<size_t>(um.strand) - 1;

        um.len = v_unitigs[um.pos_unitig]->getSeq().jump(s, pos, um.dist + ((k_ - 1) & mask), !um.strand) - k_ + 1;
        um.dist = um.dist - ((um.len - 1) & mask);
    }

    return um;
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::findUnitig(const Kmer& km, const char* s, const size_t pos, const preAllocMinHashIterator<RepHash>& it_min_h) {

    // need to check if we find it right away, need to treat this common case
    UnitigMap<U, G> um = find(km, it_min_h);

    if (!um.isEmpty && !um.isShort && !um.isAbundant){

        const size_t mask = static_cast<size_t>(um.strand) - 1;

        um.len = v_unitigs[um.pos_unitig]->getSeq().jump(s, pos, um.dist + ((k_ - 1) & mask), !um.strand) - k_ + 1;
        um.dist = um.dist - ((um.len - 1) & mask);
    }

    return um;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::addUnitig(const string& str_unitig, const size_t id_unitig){

    int pos;

    const size_t len = str_unitig.size();
    size_t pos_id_unitig = id_unitig << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    bool isShort = false, isAbundant = false, isForbidden = false;

    const char* c_str = str_unitig.c_str();

    char km_tmp[MAX_KMER_SIZE];

    Kmer km_rep;

    if (len == k_){ // Unitig to add is short, maybe abundant as well

        isShort = true;

        pos_id_unitig |= MASK_CONTIG_TYPE;

        km_rep = Kmer(c_str).rep();
        km_rep.toString(km_tmp);
        c_str = km_tmp;
    }

    minHashIterator<RepHash> it_min(c_str, len, k_, g_, RepHash(), true), it_min_end;
    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

        //If current minimizer was not seen before
        if ((last_pos_min < it_min.getPosition()) || isForbidden){

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(c_str + min_h_res.pos).rep(); //Get the minimizer to insert
                std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                packed_tiny_vector* pck_tinyv = &(p.first.getVector());
                size_t flag = p.first.getVectorSize();
                size_t v_sz = pck_tinyv->size(flag);

                pos = min_h_res.pos;
                pos_id_unitig = (pos_id_unitig & mask) | ((size_t) pos);

                if (!isShort){

                    mhr = min_h_res;

                    while ((v_sz >= max_abundance_lim) || ((v_sz > 0) && (((*pck_tinyv)(v_sz-1, flag) & mask) == mask))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if (((*pck_tinyv)(v_sz-1, flag) & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if (((*pck_tinyv)(v_sz-1, flag) & MASK_CONTIG_ID) == MASK_CONTIG_ID) (*pck_tinyv)(v_sz-1, flag) |= MASK_CONTIG_TYPE;
                                else p.first.getVectorSize() = (flag = pck_tinyv->push_back(mask, flag));
                            }

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(c_str + mhr.pos).rep();

                            p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                            pck_tinyv = &(p.first.getVector());
                            flag = p.first.getVectorSize();
                            v_sz = pck_tinyv->size(flag);
                        }
                        else break;
                    }
                }

                packed_tiny_vector& v = p.first.getVector();
                uint8_t& flag_v = p.first.getVectorSize();

                if (v_sz == 0) flag_v = v.push_back(pos_id_unitig, flag_v); //Newly created vector, just push unitig ID
                else if (isShort && (v_sz >= min_abundance_lim)){ //The minimizer (is/might be) too abundant

                    isShort = false;
                    isAbundant = true;

                    it_min = it_min_end;

                    break;
                }
                else if ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID){ //If minimizer is abundant or crowded

                    if ((v_sz == 1) || (v(v_sz-2, flag_v) != pos_id_unitig)) flag_v = v.insert(pos_id_unitig, v_sz-1, flag_v);
                }
                else if (v(v_sz-1, flag_v) != pos_id_unitig) flag_v = v.push_back(pos_id_unitig, flag_v);

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }

    if (isAbundant){

        if (id_unitig == km_unitigs.size()) km_unitigs.push_back(km_rep);
        else km_unitigs.set(id_unitig, km_rep);

        deleteUnitig_<is_void<U>::value>(true, false, id_unitig, false);

        if (id_unitig == km_unitigs.size() - 1) km_unitigs.resize(km_unitigs.size() - 1);

        it_min = minHashIterator<RepHash>(c_str, len, k_, g_, RepHash(), true);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

            if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(c_str + min_h_res.pos).rep(); //Get the minimizer to insert
                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();
                    const size_t v_sz = v.size(flag_v);

                    if ((v_sz > 0) && ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v(v_sz-1, flag_v)++;
                    else flag_v = v.push_back(MASK_CONTIG_ID + 1, flag_v);

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }

        h_kmers_ccov.insert(km_rep, CompressedCoverage_t<U>(1));
    }
    else if (isShort){

        if (id_unitig == km_unitigs.size()) km_unitigs.push_back(km_rep);
        else km_unitigs.set(id_unitig, km_rep);
    }
    else if (id_unitig == v_unitigs.size()) v_unitigs.push_back(new Unitig<U>(c_str)); //Push unitig to list of unitigs
    else v_unitigs[id_unitig] = new Unitig<U>(c_str);

    return isAbundant;
}

template<typename U, typename G>
void CompactedDBG<U, G>::moveToAbundant() {

    MinimizerIndex::iterator it = hmap_min_unitigs.begin();
    MinimizerIndex::iterator it_end = hmap_min_unitigs.end();

    vector<pair<size_t, string>> to_move;

    char km_tmp[MAX_KMER_SIZE];

    Kmer km_rep;

    to_move.reserve(16);

    while (it != it_end){

        const packed_tiny_vector* pck_tinyv = &(it.getVector());
        const size_t flag = it.getVectorSize();
        const size_t v_sz = pck_tinyv->size(flag);

        if (v_sz >= min_abundance_lim){

            for (size_t i = 0; i < v_sz; ++i){

                if (static_cast<bool>((*pck_tinyv)(i, flag) | MASK_CONTIG_TYPE)){

                    size_t unitig_id = (*pck_tinyv)(i, flag) >> 32;
                    string unitig_str = km_unitigs.getKmer(unitig_id).toString();

                    to_move.push_back({std::move(unitig_id), std::move(unitig_str)});
                }
            }

            if (min_abundance_lim > (v_sz - to_move.size())){

                for (const auto& p : to_move) deleteUnitig_(true, false, p.first, p.second);

                for (const auto& p : to_move){

                    const char* c_str = p.second.c_str();
                    const size_t len = p.second.length();

                    km_rep = Kmer(c_str).rep();
                    c_str = km_tmp;

                    km_rep.toString(km_tmp);

                    minHashIterator<RepHash> it_min = minHashIterator<RepHash>(c_str, len, k_, g_, RepHash(), true), it_min_end;

                    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

                        if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                            while (it_it_min != it_it_min_end){

                                const minHashResult& min_h_res = *it_it_min;
                                Minimizer minz_rep = Minimizer(c_str + min_h_res.pos).rep(); //Get the minimizer to insert

                                std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                                packed_tiny_vector& v = p.first.getVector();
                                uint8_t& flag_v = p.first.getVectorSize();
                                const size_t v_sz = v.size(flag_v);

                                if ((v_sz > 0) && ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v(v_sz-1, flag_v)++;
                                else flag_v = v.push_back(MASK_CONTIG_ID + 1, flag_v);

                                last_pos_min = min_h_res.pos;
                                ++it_it_min;
                            }
                        }
                    }

                    h_kmers_ccov.insert(km_rep, CompressedCoverage_t<U>(1));
                }
            }

            to_move.clear();
        }

        ++it;
    }
}

template<typename U, typename G>
bool CompactedDBG<U, G>::addUnitig(const string& str_unitig, const size_t id_unitig, SpinLock& lck_unitig, SpinLock& lck_kmer){

    int pos;

    const size_t len = str_unitig.size();
    size_t pos_id_unitig = id_unitig << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    bool isShort = false, isAbundant = false, isForbidden = false;

    const char* c_str = str_unitig.c_str();

    char km_tmp[MAX_KMER_SIZE];

    Kmer km_rep;

    if (len == k_){ // Unitig to add is short, maybe abundant as well

        isShort = true;

        pos_id_unitig |= MASK_CONTIG_TYPE;

        km_rep = Kmer(c_str).rep();
        km_rep.toString(km_tmp);
        c_str = km_tmp;
    }

    minHashIterator<RepHash> it_min(c_str, len, k_, g_, RepHash(), true), it_min_end;
    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

        //If current minimizer was not seen before
        if ((last_pos_min < it_min.getPosition()) || isForbidden){

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(c_str + min_h_res.pos).rep(); //Get the minimizer to insert

                pos = min_h_res.pos;
                pos_id_unitig = (pos_id_unitig & mask) | ((size_t) pos);

                std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert_p(minz_rep, packed_tiny_vector(), 0);
                packed_tiny_vector* pck_tinyv = &(p.first.getVector());
                size_t flag = p.first.getVectorSize();
                size_t v_sz = pck_tinyv->size(flag);

                if (!isShort){

                    mhr = min_h_res;

                    while ((v_sz >= max_abundance_lim) || ((v_sz > 0) && (((*pck_tinyv)(v_sz-1, flag) & mask) == mask))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if (((*pck_tinyv)(v_sz-1, flag) & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if (((*pck_tinyv)(v_sz-1, flag) & MASK_CONTIG_ID) == MASK_CONTIG_ID) (*pck_tinyv)(v_sz-1, flag) |= MASK_CONTIG_TYPE;
                                else p.first.getVectorSize() = (flag = pck_tinyv->push_back(mask, flag));
                            }

                            hmap_min_unitigs.release_p(p.first);

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(c_str + mhr.pos).rep();

                            p = hmap_min_unitigs.insert_p(minz_rep, packed_tiny_vector(), 0);
                            pck_tinyv = &(p.first.getVector());
                            flag = p.first.getVectorSize();
                            v_sz = pck_tinyv->size(flag);
                        }
                        else break;
                    }
                }

                packed_tiny_vector& v = p.first.getVector();
                uint8_t& flag_v = p.first.getVectorSize();

                if (v_sz == 0) flag_v = v.push_back(pos_id_unitig, flag_v); //Newly created vector, just push unitig ID
                else if ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID){ //If minimizer is abundant or crowded

                    if ((v_sz == 1) || (v(v_sz-2, flag_v) != pos_id_unitig)) flag_v = v.insert(pos_id_unitig, v_sz-1, flag_v);
                }
                else if (v(v_sz-1, flag_v) != pos_id_unitig) flag_v = v.push_back(pos_id_unitig, flag_v);

                hmap_min_unitigs.release_p(p.first);

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }

    if (isShort){

        lck_kmer.acquire();

        if (id_unitig >= km_unitigs.size()) km_unitigs.resize(id_unitig + 1);

        km_unitigs.set(id_unitig, km_rep);

        lck_kmer.release();
    }
    else {

        lck_unitig.acquire();

        if (id_unitig >= v_unitigs.size()) v_unitigs.resize(id_unitig + 1);

        v_unitigs[id_unitig] = new Unitig<U>(c_str); //Push unitig to list of unitigs

        lck_unitig.release();
    }

    return isAbundant;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::addUnitig(const string& str_unitig, const size_t id_unitig, const size_t id_unitig_r, const size_t is_short_r){

    int pos;

    const size_t pos_id_unitig_r = (id_unitig_r << 32) | (is_short_r ? MASK_CONTIG_TYPE : 0);

    const size_t len = str_unitig.size();
    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    size_t pos_id_unitig = id_unitig << 32;
    size_t i = 0;

    bool isShort = false, isAbundant = false, isForbidden = false;

    const char* c_str = str_unitig.c_str();

    char km_tmp[MAX_KMER_SIZE];

    Kmer km_rep;

    if (len == k_){ // Unitig to add is short, maybe abundant as well

        isShort = true;

        pos_id_unitig |= MASK_CONTIG_TYPE;

        km_rep = Kmer(c_str).rep();
        km_rep.toString(km_tmp);
        c_str = km_tmp;
    }

    minHashIterator<RepHash> it_min(c_str, len, k_, g_, RepHash(), true), it_min_end;
    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

        //If current minimizer was not seen before
        if ((last_pos_min < it_min.getPosition()) || isForbidden){

            minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(c_str + min_h_res.pos).rep(); //Get the minimizer to insert
                std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                packed_tiny_vector* pck_tinyv = &(p.first.getVector());
                size_t flag = p.first.getVectorSize();
                size_t v_sz = pck_tinyv->size(flag);

                pos = min_h_res.pos;
                pos_id_unitig = (pos_id_unitig & mask) | ((size_t) pos);

                for (i = 0; i < v_sz; ++i){

                    if (((*pck_tinyv)(i, flag) & mask) == pos_id_unitig_r){

                        (*pck_tinyv)(i, flag) = pos_id_unitig;
                        break;
                    }
                }

                if (!isShort){

                    mhr = min_h_res;

                    while ((i == v_sz) && ((v_sz >= max_abundance_lim) || ((v_sz > 0) && (((*pck_tinyv)(v_sz-1, flag) & mask) == mask)))){

                        mhr_tmp = it_min.getNewMin(mhr);
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            if (((*pck_tinyv)(v_sz-1, flag) & mask) != mask){ // Minimizer was never signaled before as overcrowded

                                // If minimizer bin already contains abundant k-mer, just set flag for unitig overcrowding
                                if (((*pck_tinyv)(v_sz-1, flag) & MASK_CONTIG_ID) == MASK_CONTIG_ID) (*pck_tinyv)(v_sz-1, flag) |= MASK_CONTIG_TYPE;
                                else p.first.getVectorSize() = (flag = pck_tinyv->push_back(mask, flag));
                            }

                            mhr = mhr_tmp;
                            minz_rep = Minimizer(c_str + mhr.pos).rep();

                            p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                            pck_tinyv = &(p.first.getVector());
                            flag = p.first.getVectorSize();
                            v_sz = pck_tinyv->size(flag);

                            for (i = 0; i < v_sz; ++i){

                                if (((*pck_tinyv)(i, flag) & mask) == pos_id_unitig_r){

                                    (*pck_tinyv)(i, flag) = pos_id_unitig;
                                    break;
                                }
                            }
                        }
                        else {

                            i = v_sz;
                            break;
                        }
                    }
                }

                if (i == v_sz){

                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();

                    if (v_sz == 0) flag_v = v.push_back(pos_id_unitig, flag_v); //Newly created vector, just push unitig ID
                    else if (isShort && (v_sz >= min_abundance_lim)){ //The minimizer (is/might be) too abundant

                        isShort = false;
                        isAbundant = true;

                        it_min = it_min_end;

                        break;
                    }
                    else if ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID){ //If minimizer is abundant or crowded

                        if ((v_sz == 1) || (v(v_sz-2, flag_v) != pos_id_unitig)) flag_v = v.insert(pos_id_unitig, v_sz-1, flag_v);
                    }
                    else if (v(v_sz-1, flag_v) != pos_id_unitig) flag_v = v.push_back(pos_id_unitig, flag_v);
                }

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }

    if (isAbundant){

        if (id_unitig == km_unitigs.size()) km_unitigs.push_back(km_rep);
        else km_unitigs.set(id_unitig, km_rep);

        deleteUnitig_<is_void<U>::value>(true, false, id_unitig, false);

        if (id_unitig == km_unitigs.size() - 1) km_unitigs.resize(km_unitigs.size() - 1);

        it_min = minHashIterator<RepHash>(c_str, len, k_, g_, RepHash(), true);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){

            if (last_pos_min < it_min.getPosition()){ //If current minimizer was not seen before

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(c_str + min_h_res.pos).rep(); //Get the minimizer to insert
                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);
                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();
                    const size_t v_sz = v.size(flag_v);

                    if ((v_sz > 0) && ((v(v_sz-1, flag_v) & MASK_CONTIG_ID) == MASK_CONTIG_ID)) v(v_sz-1, flag_v)++;
                    else flag_v = v.push_back(MASK_CONTIG_ID + 1, flag_v);

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }

        h_kmers_ccov.insert(km_rep, CompressedCoverage_t<U>(1));
    }
    else if (isShort){

        if (id_unitig == km_unitigs.size()) km_unitigs.push_back(km_rep);
        else km_unitigs.set(id_unitig, km_rep);
    }
    else if (id_unitig == v_unitigs.size()) v_unitigs.push_back(new Unitig<U>(c_str)); //Push unitig to list of unitigs
    else v_unitigs[id_unitig] = new Unitig<U>(c_str);

    return isAbundant;
}

template<typename U, typename G>
void CompactedDBG<U, G>::swapUnitigs(const bool isShort, const size_t id_a, const size_t id_b){

    size_t shift_id_unitig_a = id_a << 32;
    size_t shift_id_unitig_b = id_b << 32;

    const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

    string str;

    size_t h_sz = 0;

    // Swap the unitig pointers in v_unitigs
    if (isShort){

        km_unitigs.swap(id_a, id_b);

        shift_id_unitig_a |= MASK_CONTIG_TYPE;
        shift_id_unitig_b |= MASK_CONTIG_TYPE;

        str = km_unitigs.getKmer(id_a).toString();
        h_sz = 4;
    }
    else {

        std::swap(v_unitigs[id_a], v_unitigs[id_b]);

        str = v_unitigs[id_a]->getSeq().toString();
        h_sz = (v_unitigs[id_a]->length() / 4) + (v_unitigs[id_b]->length() / 4);
    }

    MinimizerHashTable<uint8_t> minimizers(h_sz);

    bool isForbidden = false;

    // Swap the unitig IDs in the minimizer hash table
    const char* s = str.c_str();

    size_t len = str.size();

    minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;
    minHashResult mhr, mhr_tmp;

    for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig

        if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer is found in unitig

            minHashResultIterator<RepHash> it_it_min(*it_min), it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(s + min_h_res.pos).rep();
                MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                if (it_h != hmap_min_unitigs.end()){

                    if (minimizers.find(minz_rep) == minimizers.end()){

                        minimizers.insert(minz_rep, 0);

                        packed_tiny_vector& v = it_h.getVector();
                        const uint8_t flag_v = it_h.getVectorSize();
                        const size_t v_sz = v.size(flag_v);

                        for (size_t i = 0; i < v_sz; ++i){

                             // Swap the unitig ids but do not change positions;
                            if ((v(i, flag_v) & mask) == shift_id_unitig_b) v(i, flag_v) = shift_id_unitig_a | (v(i, flag_v) & MASK_CONTIG_POS);
                            else if ((v(i, flag_v) & mask) == shift_id_unitig_a) v(i, flag_v) = shift_id_unitig_b | (v(i, flag_v) & MASK_CONTIG_POS);
                        }
                    }

                    packed_tiny_vector* v = &(it_h.getVector());
                    uint8_t flag_v = it_h.getVectorSize();
                    size_t v_sz = v->size(flag_v);

                    mhr = min_h_res;

                    while (((*v)(v_sz-1, flag_v) & mask) == mask){

                        mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            minz_rep = Minimizer(s + mhr_tmp.pos).rep();
                            it_h = hmap_min_unitigs.find(minz_rep);

                            if (it_h == hmap_min_unitigs.end()) break;

                            if (minimizers.find(minz_rep) == minimizers.end()){

                                minimizers.insert(minz_rep, 0);

                                packed_tiny_vector& v_tmp = it_h.getVector();
                                const uint8_t flag_v_tmp = it_h.getVectorSize();
                                const size_t v_sz_tmp = v_tmp.size(flag_v_tmp);

                                for (size_t i = 0; i < v_sz_tmp; ++i){

                                     // Swap the unitig ids but do not change positions;
                                    if ((v_tmp(i, flag_v_tmp) & mask) == shift_id_unitig_b){

                                        v_tmp(i, flag_v_tmp) = shift_id_unitig_a | (v_tmp(i, flag_v_tmp) & MASK_CONTIG_POS);
                                    }
                                    else if ((v_tmp(i, flag_v_tmp) & mask) == shift_id_unitig_a){

                                        v_tmp(i, flag_v_tmp) = shift_id_unitig_b | (v_tmp(i, flag_v_tmp) & MASK_CONTIG_POS);
                                    }
                                }
                            }

                            mhr = mhr_tmp;
                            v = &(it_h.getVector());
                            flag_v = it_h.getVectorSize();
                            v_sz = v->size(flag_v);
                        }
                        else break;
                    }
                }

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }

    str = isShort ? km_unitigs.getKmer(id_b).toString() : v_unitigs[id_b]->getSeq().toString();
    s = str.c_str();
    len = str.size();

    isForbidden = false;

    minHashIterator<RepHash> it_min2(s, len, k_, g_, RepHash(), true);

    for (int64_t last_pos_min = -1; it_min2 != it_min_end; ++it_min2){ // Iterate over minimizers of unitig

        if ((last_pos_min < it_min2.getPosition()) || isForbidden){ // If a new minimizer is found in unitig

            minHashResultIterator<RepHash> it_it_min(*it_min2), it_it_min_end;
            isForbidden = false;

            while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig

                const minHashResult& min_h_res = *it_it_min;
                Minimizer minz_rep = Minimizer(s + min_h_res.pos).rep();
                MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                if (it_h != hmap_min_unitigs.end()){

                    if (minimizers.find(minz_rep) == minimizers.end()){

                        minimizers.insert(minz_rep, 0);

                        packed_tiny_vector& v = it_h.getVector();
                        const uint8_t flag_v = it_h.getVectorSize();
                        const size_t v_sz = v.size(flag_v);

                        for (size_t i = 0; i < v_sz; ++i){

                            if ((v(i, flag_v) & mask) == shift_id_unitig_a) v(i, flag_v) = shift_id_unitig_b | (v(i, flag_v) & MASK_CONTIG_POS);
                        }
                    }

                    packed_tiny_vector* v = &(it_h.getVector());
                    uint8_t flag_v = it_h.getVectorSize();
                    size_t v_sz = v->size(flag_v);

                    mhr = min_h_res;

                    while (((*v)(v_sz-1, flag_v) & mask) == mask){

                        mhr_tmp = it_min2.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                        isForbidden = true;

                        if (mhr_tmp.hash != mhr.hash){

                            minz_rep = Minimizer(s + mhr_tmp.pos).rep();
                            it_h = hmap_min_unitigs.find(minz_rep);

                            if (it_h == hmap_min_unitigs.end()) break;

                            if (minimizers.find(minz_rep) == minimizers.end()){

                                minimizers.insert(minz_rep, 0);

                                packed_tiny_vector& v_tmp = it_h.getVector();
                                const uint8_t flag_v_tmp = it_h.getVectorSize();
                                const size_t v_sz_tmp = v_tmp.size(flag_v_tmp);

                                for (size_t i = 0; i < v_sz_tmp; ++i){

                                    if ((v_tmp(i, flag_v_tmp) & mask) == shift_id_unitig_a) {

                                        v_tmp(i, flag_v_tmp) = shift_id_unitig_b | (v_tmp(i, flag_v_tmp) & MASK_CONTIG_POS);
                                    }
                                }
                            }

                            mhr = mhr_tmp;
                            v = &(it_h.getVector());
                            flag_v = it_h.getVectorSize();
                            v_sz = v->size(flag_v);
                        }
                        else break;
                    }
                }

                last_pos_min = min_h_res.pos;
                ++it_it_min;
            }
        }
    }
}

// Input sequence must not contain duplicated k-mers
// If it does, use CompactedDBG<U, G>::add().
template<typename U, typename G>
bool CompactedDBG<U, G>::mergeUnitig(const string& seq, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::mergeUnitig(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::mergeUnitig(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = km_unitigs.size();

    size_t nb_curr_pred = 0;
    size_t nb_curr_succ = 0;
    size_t nb_prev_pred = 0xFFFFFFFFFFFFFFFF;
    size_t nb_prev_succ = 0xFFFFFFFFFFFFFFFF;

    size_t added = 0;
    size_t split_before = 0, split_after = 0;

    bool prev_found = true;

    string curr_unitig;

    vector<Kmer> v_joins;

    KmerHashTable<vector<size_t>> kht;

    const char* str_seq = seq.c_str();

    const size_t str_seq_len = seq.length();

    auto add_graph_function = [&](const string& unitig){

        const char* str_unitig = unitig.c_str();
        const size_t len_unitig = unitig.length();

        if (len_unitig == k_){

            if (!addUnitig(str_unitig, v_kmers_sz)){

                km_unitigs.setFull(v_kmers_sz);

                ++v_kmers_sz;
                ++added;
            }
            else h_kmers_ccov.find(Kmer(str_unitig).rep())->ccov.setFull();
        }
        else {

            addUnitig(str_unitig, nxt_pos_insert_v_unitigs);

            v_unitigs[nxt_pos_insert_v_unitigs]->getCov().setFull();

            ++nxt_pos_insert_v_unitigs;
            ++added;

            v_joins.push_back(Kmer(str_unitig + len_unitig - k_));
        }

        v_joins.push_back(Kmer(str_unitig));
    };

    auto add_split_function = [&](){

        KmerHashTable<vector<size_t>>::iterator it(kht.begin());
        const KmerHashTable<vector<size_t>>::iterator it_end(kht.end());

        vector<pair<int,int>> sp;

        while (it != it_end){

            ++split_before;
            ++split_after;

            const_UnitigMap<U, G> um(find(it.getKey(), true));

            vector<size_t>& split_v = *it;

            sort(split_v.begin(), split_v.end());

            size_t prev_split_pos = 0;

            for (const auto pos : split_v){

                if (pos != prev_split_pos){

                    sp.push_back({prev_split_pos, pos});

                    v_joins.push_back(um.getUnitigKmer(prev_split_pos));
                    if (prev_split_pos != pos - 1) v_joins.push_back(um.getUnitigKmer(pos - 1));

                    prev_split_pos = pos;

                    ++split_after;
                }
            }

            sp.push_back({prev_split_pos, um.size - k_ + 1});

            v_joins.push_back(um.getUnitigKmer(prev_split_pos));
            if (prev_split_pos != um.size - k_) v_joins.push_back(um.getUnitigKmer(um.size - k_));

            extractUnitig_<is_void<U>::value>(um.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

            sp.clear();
            split_v.clear();

            ++it;
        }

        kht.clear();
    };

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        UnitigMap<U, G> cm = findUnitig(p.first, str_seq, p.second);

        if (cm.isEmpty){

            vector<UnitigMap<U, G>> um_bw(findPredecessors(p.first));
            vector<UnitigMap<U, G>> um_fw(findSuccessors(p.first));

            nb_curr_pred = 0;
            nb_curr_succ = 0;

            for (auto& um : um_bw){

                if (!um.isEmpty){

                    ++nb_curr_pred;

                    if (!um.isAbundant && !um.isShort){

                        um.dist += um.strand;
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) kht.insert(um.getUnitigHead(), vector<size_t>()).first->push_back(um.dist);
                    }
                }
            }

            for (auto& um : um_fw){

                if (!um.isEmpty){

                    ++nb_curr_succ;

                    if (!um.isAbundant && !um.isShort){

                        um.dist += !um.strand;
                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) kht.insert(um.getUnitigHead(), vector<size_t>()).first->push_back(um.dist);
                    }
                }
            }

            if (prev_found){ // Previous k-mer was found in the graph => current k-mer (not found in the graph) starts a new unitig

                prev_found = false;
                curr_unitig = p.first.toString();
            }
            else if ((nb_prev_succ == 0) && (nb_curr_pred == 0)) curr_unitig.push_back(str_seq[p.second + k_ - 1]);
            else {

                add_graph_function(curr_unitig);

                curr_unitig = p.first.toString();
            }

            nb_prev_succ = nb_curr_succ;
            nb_prev_pred = nb_curr_pred;

            ++it_km;
        }
        else {

            if ((p.second == 0) && !cm.isAbundant && !cm.isShort){

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)) kht.insert(cm.getUnitigHead(), vector<size_t>()).first->push_back(cm.dist);
            }
            else if ((p.second + cm.len == str_seq_len - k_ + 1) && !cm.isAbundant && !cm.isShort){

                cm.dist += cm.len;

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)) kht.insert(cm.getUnitigHead(), vector<size_t>()).first->push_back(cm.dist);
            }

            it_km += cm.len;

            if (!prev_found){

                add_graph_function(curr_unitig);

                prev_found = true;
            }
        }
    }

    add_split_function();

    if (!prev_found) add_graph_function(curr_unitig);

    if (nxt_pos_insert_v_unitigs < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert_v_unitigs);
    if (v_kmers_sz < km_unitigs.size()) km_unitigs.resize(v_kmers_sz);

    const size_t joined = joinUnitigs_<is_void<U>::value>(&v_joins);

    if (verbose){

        cout << "CompactedDBG::mergeUnitig(): Added " << added << " new unitigs to the graph." << endl;
        cout << "CompactedDBG::mergeUnitig(): Split " << split_before << " unitigs into " << split_after << " new unitigs." << endl;
        cout << "CompactedDBG::mergeUnitig(): Joined " << joined << " unitigs from the graph." << endl;
    }

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::annotateSplitUnitig(const string& seq, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::annotateSplitUnitig(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::annotateSplitUnitig(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_kmers_sz = km_unitigs.size();

    size_t nb_curr_pred = 0;
    size_t nb_curr_succ = 0;
    size_t nb_prev_pred = 0xFFFFFFFFFFFFFFFF;
    size_t nb_prev_succ = 0xFFFFFFFFFFFFFFFF;

    bool prev_found = true;

    string curr_unitig;

    const char* str_seq = seq.c_str();

    const size_t str_seq_len = seq.length();

    auto add_graph_function = [&](const string& unitig){

        const char* str_unitig = unitig.c_str();
        const size_t len_unitig = unitig.length();

        if (len_unitig == k_){

            if (!addUnitig(str_unitig, v_kmers_sz)) km_unitigs.setFull(v_kmers_sz++);
            else h_kmers_ccov.find(Kmer(str_unitig).rep())->ccov.setFull();
        }
        else {

            addUnitig(str_unitig, nxt_pos_insert_v_unitigs);

            v_unitigs[nxt_pos_insert_v_unitigs++]->getCov().setFull();
        }
    };

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        UnitigMap<U, G> cm = findUnitig(p.first, str_seq, p.second);

        if (cm.isEmpty){

            vector<UnitigMap<U, G>> um_bw(findPredecessors(p.first));
            vector<UnitigMap<U, G>> um_fw(findSuccessors(p.first));

            nb_curr_pred = 0;
            nb_curr_succ = 0;

            for (auto& um : um_bw){

                if (!um.isEmpty){

                    ++nb_curr_pred;

                    if (!um.isAbundant && !um.isShort){

                        um.dist += um.strand;

                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) unmapRead(um);
                    }
                }
            }

            for (auto& um : um_fw){

                if (!um.isEmpty){

                    ++nb_curr_succ;

                    if (!um.isAbundant && !um.isShort){

                        um.dist += !um.strand;

                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) unmapRead(um);
                    }
                }
            }

            if (prev_found){ // Previous k-mer was found in the graph => current k-mer (not found in the graph) starts a new unitig

                prev_found = false;
                curr_unitig = p.first.toString();
            }
            else if (((nb_prev_succ == 0) && (nb_curr_pred == 0))) curr_unitig.push_back(str_seq[p.second + k_ - 1]);
            else {

                add_graph_function(curr_unitig);

                curr_unitig = p.first.toString();
            }

            nb_prev_succ = nb_curr_succ;
            nb_prev_pred = nb_curr_pred;

            ++it_km;
        }
        else {

            if ((p.second == 0) && !cm.isAbundant && !cm.isShort){

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;
                    unmapRead(cm);
                    cm.len = len;
                }
            }
            else if ((p.second + cm.len == str_seq_len - k_ + 1) && !cm.isAbundant && !cm.isShort){

                cm.dist += cm.len;

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;
                    unmapRead(cm);
                    cm.len = len;
                }
            }

            it_km += cm.len;

            if (!prev_found){

                add_graph_function(curr_unitig);

                prev_found = true;
            }
        }
    }

    if (!prev_found) add_graph_function(curr_unitig);

    return true;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::annotateSplitUnitig(const string& seq, LockGraph& lck_g, const bool verbose){

    if (invalid){

        cerr << "CompactedDBG::annotateSplitUnitig(): Graph is invalid and no sequence can be added to it" << endl;
        return false;
    }

    if (seq.length() < k_){

        cerr << "CompactedDBG::annotateSplitUnitig(): Input sequence length cannot be less than k = " << k_ << endl;
        return false;
    }

    size_t nb_curr_pred = 0;
    size_t nb_curr_succ = 0;
    size_t nb_prev_pred = 0xFFFFFFFFFFFFFFFF;
    size_t nb_prev_succ = 0xFFFFFFFFFFFFFFFF;

    bool prev_found = true;

    string curr_unitig;

    const char* str_seq = seq.c_str();

    const size_t str_seq_len = seq.length();

    auto add_graph_function = [&](const string& unitig){

        const char* str_unitig = unitig.c_str();
        const size_t len_unitig = unitig.length();

        lck_g.acquire_writer();

        if (len_unitig == k_){

            if (!addUnitig(str_unitig, km_unitigs.size())) km_unitigs.setFull(km_unitigs.size() - 1);
            else h_kmers_ccov.find(Kmer(str_unitig).rep())->ccov.setFull();
        }
        else {

            addUnitig(str_unitig, v_unitigs.size());

            v_unitigs[v_unitigs.size() - 1]->getCov().setFull();
        }

        lck_g.release_writer();
    };

    for (KmerIterator it_km(str_seq), it_km_end; it_km != it_km_end;) { //non-ACGT char. are discarded

        const std::pair<Kmer, int>& p = *it_km;

        lck_g.acquire_reader();

        UnitigMap<U, G> cm = findUnitig(p.first, str_seq, p.second);

        if (cm.isEmpty){

            vector<UnitigMap<U, G>> um_bw(findPredecessors(p.first));
            vector<UnitigMap<U, G>> um_fw(findSuccessors(p.first));

            nb_curr_pred = 0;
            nb_curr_succ = 0;

            for (auto& um : um_bw){

                if (!um.isEmpty){

                    ++nb_curr_pred;

                    if (!um.isAbundant && !um.isShort){

                        um.dist += um.strand;

                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) unmapRead(um, lck_g);
                    }
                }
            }

            for (auto& um : um_fw){

                if (!um.isEmpty){

                    ++nb_curr_succ;

                    if (!um.isAbundant && !um.isShort){

                        um.dist += !um.strand;

                        if ((um.dist != 0) && (um.dist != um.size - k_ + 1)) unmapRead(um, lck_g);
                    }
                }
            }

            lck_g.release_reader();

            if (prev_found){ // Previous k-mer was found in the graph => current k-mer (not found in the graph) starts a new unitig

                prev_found = false;
                curr_unitig = p.first.toString();
            }
            else if (((nb_prev_succ == 0) && (nb_curr_pred == 0))) curr_unitig.push_back(str_seq[p.second + k_ - 1]);
            else {

                add_graph_function(curr_unitig);

                curr_unitig = p.first.toString();
            }

            nb_prev_succ = nb_curr_succ;
            nb_prev_pred = nb_curr_pred;

            ++it_km;
        }
        else {

            if ((p.second == 0) && !cm.isAbundant && !cm.isShort){

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;
                    unmapRead(cm, lck_g);
                    cm.len = len;
                }
            }
            else if ((p.second + cm.len == str_seq_len - k_ + 1) && !cm.isAbundant && !cm.isShort){

                cm.dist += cm.len;

                if ((cm.dist != 0) && (cm.dist != cm.size - k_ + 1)){

                    const size_t len = cm.len;

                    cm.len = 1;
                    unmapRead(cm, lck_g);
                    cm.len = len;
                }
            }

            lck_g.release_reader();

            it_km += cm.len;

            if (!prev_found){

                add_graph_function(curr_unitig);

                prev_found = true;
            }
        }
    }

    if (!prev_found) add_graph_function(curr_unitig);

    return true;
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, void>::type CompactedDBG<U, G>::deleteUnitig_(const bool isShort, const bool isAbundant,
                                                                                const size_t id_unitig, const bool delete_data){

    if (isAbundant){

        char km_str[MAX_KMER_SIZE];

        typename h_kmers_ccov_t::iterator it_h = h_kmers_ccov.find(id_unitig);

        const Kmer km = it_h.getKey();

        km.toString(km_str);

        minHashIterator<RepHash> it_min(km_str, k_, k_, g_, RepHash(), true), it_min_end;

        if (delete_data) it_h->getData()->clear(UnitigMap<U, G>(id_unitig, 0, 1, k_, isShort, isAbundant, true, this));

        it_h->ccov.clear();
        h_kmers_ccov.erase(it_h);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(km_str + min_h_res.pos).rep(); // Get canonical minimizer
                    MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVector();
                        uint8_t& flag_v = it_h.getVectorSize();

                        const size_t last_pos_v = v.size(flag_v) - 1;

                        v(last_pos_v, flag_v)--;

                        if (((v(last_pos_v, flag_v) & RESERVED_ID) == 0) && ((v(last_pos_v, flag_v) & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_unitigs.erase(minz_rep);
                            else flag_v = v.remove(v.size(flag_v) - 1, flag_v);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
    else {

        bool isForbidden = false;

        size_t pos_id_unitig = id_unitig << 32;
        const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

        string str;

        if (isShort){

            str = km_unitigs.getKmer(id_unitig).toString();
            pos_id_unitig |= MASK_CONTIG_TYPE;
        }
        else str = v_unitigs[id_unitig]->getSeq().toString();

        const char* s = str.c_str();

        const size_t len = str.size();

        minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;
        minHashResult mhr, mhr_tmp;

        // The unitig is deleted but its space in the unitig vector is not because:
        // 1 - It would change indices in the minimizer hash table
        if (isShort){

            if (delete_data) km_unitigs.getData(id_unitig)->clear(UnitigMap<U, G>(id_unitig, 0, 1, len, isShort, isAbundant, true, this));

            Kmer empty_km;

            empty_km.set_deleted();

            km_unitigs.set(id_unitig, empty_km, 0);
        }
        else {

            if (delete_data) v_unitigs[id_unitig]->getData()->clear(UnitigMap<U, G>(id_unitig, 0, len - k_ + 1, len, isShort, isAbundant, true, this));

            delete v_unitigs[id_unitig];
            v_unitigs[id_unitig] = nullptr;
        }

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
                isForbidden = false;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(s + min_h_res.pos).rep(); // Get canonical minimizer
                    MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    mhr = min_h_res;

                    while (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVector();
                        uint8_t& flag_v = it_h.getVectorSize();
                        size_t i = 0, v_sz = v.size(flag_v);

                        while (i < v_sz){

                            if ((v(i, flag_v) & mask) == pos_id_unitig){

                                flag_v = v.remove(i, flag_v);
                                --v_sz;
                            }
                            else ++i;
                        }

                        it_h = hmap_min_unitigs.end();

                        if (v.size(flag_v) == 0) hmap_min_unitigs.erase(minz_rep);
                        else if (!isShort && ((v(v_sz-1, flag_v) & mask) == mask)){ //Minimizer bin is overcrowded


                            mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                            isForbidden = true;

                            if (mhr_tmp.hash != mhr.hash){

                                mhr = mhr_tmp;
                                minz_rep = Minimizer(s + mhr.pos).rep();
                                it_h = hmap_min_unitigs.find(minz_rep);
                            }
                            else break;
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, void>::type CompactedDBG<U, G>::deleteUnitig_( const bool isShort, const bool isAbundant,
                                                                                const size_t id_unitig, const bool delete_data){

    if (isAbundant){

        char km_str[MAX_KMER_SIZE];

        typename h_kmers_ccov_t::iterator it_h = h_kmers_ccov.find(id_unitig);

        const Kmer km = it_h.getKey();

        km.toString(km_str);

        minHashIterator<RepHash> it_min(km_str, k_, k_, g_, RepHash(), true), it_min_end;

        it_h->ccov.clear();
        h_kmers_ccov.erase(it_h);

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(km_str + min_h_res.pos).rep(); // Get canonical minimizer
                    MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVector();
                        uint8_t& flag_v = it_h.getVectorSize();

                        const size_t last_pos_v = v.size(flag_v) - 1;

                        v(last_pos_v, flag_v)--;

                        if (((v(last_pos_v, flag_v) & RESERVED_ID) == 0) && ((v(last_pos_v, flag_v) & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_unitigs.erase(minz_rep);
                            else flag_v = v.remove(v.size(flag_v) - 1, flag_v);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
    else {

        bool isForbidden = false;

        size_t pos_id_unitig = id_unitig << 32;
        const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

        string str;

        if (isShort){

            str = km_unitigs.getKmer(id_unitig).toString();
            pos_id_unitig |= MASK_CONTIG_TYPE;
        }
        else str = v_unitigs[id_unitig]->getSeq().toString();

        const char* s = str.c_str();

        const size_t len = str.size();

        minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;
        minHashResult mhr, mhr_tmp;

        // The unitig is deleted but its space in the unitig vector is not because:
        // 1 - It would change indices in the minimizer hash table
        if (isShort){

            Kmer empty_km;

            empty_km.set_deleted();
            km_unitigs.set(id_unitig, empty_km, 0);
        }
        else {

            delete v_unitigs[id_unitig];
            v_unitigs[id_unitig] = nullptr;
        }

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
                isForbidden = false;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(s + min_h_res.pos).rep(); // Get canonical minimizer
                    MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    mhr = min_h_res;

                    while (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVector();
                        uint8_t& flag_v = it_h.getVectorSize();
                        size_t i = 0, v_sz = v.size(flag_v);

                        while (i < v_sz){

                            if ((v(i, flag_v) & mask) == pos_id_unitig){

                                flag_v = v.remove(i, flag_v);
                                --v_sz;
                            }
                            else ++i;
                        }

                        it_h = hmap_min_unitigs.end();

                        if (v.size(flag_v) == 0) hmap_min_unitigs.erase(minz_rep);
                        else if (!isShort && ((v(v_sz-1, flag_v) & mask) == mask)){ //Minimizer bin is overcrowded


                            mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                            isForbidden = true;

                            if (mhr_tmp.hash != mhr.hash){

                                mhr = mhr_tmp;
                                minz_rep = Minimizer(s + mhr.pos).rep();
                                it_h = hmap_min_unitigs.find(minz_rep);
                            }
                            else break;
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::deleteUnitig_(const bool isShort, const bool isAbundant, const size_t id_unitig, const string& str){

    const char* s = str.c_str();
    const size_t len = str.size();

    if (isAbundant){

        minHashIterator<RepHash> it_min(s, k_, k_, g_, RepHash(), true), it_min_end;

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if (last_pos_min < it_min.getPosition()){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(s + min_h_res.pos).rep(); // Get canonical minimizer
                    MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    if (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVector();
                        uint8_t& flag_v = it_h.getVectorSize();

                        const size_t last_pos_v = v.size(flag_v) - 1;

                        v(last_pos_v, flag_v)--;

                        if (((v(last_pos_v, flag_v) & RESERVED_ID) == 0) && ((v(last_pos_v, flag_v) & MASK_CONTIG_TYPE) != MASK_CONTIG_TYPE)){

                            if (last_pos_v == 0) hmap_min_unitigs.erase(it_h);
                            else flag_v = v.remove(v.size(flag_v) - 1, flag_v);
                        }
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
    else {

        bool isForbidden = false;

        size_t pos_id_unitig = id_unitig << 32;
        const size_t mask = MASK_CONTIG_ID | MASK_CONTIG_TYPE;

        if (isShort) pos_id_unitig |= MASK_CONTIG_TYPE;

        minHashIterator<RepHash> it_min(s, len, k_, g_, RepHash(), true), it_min_end;
        minHashResult mhr, mhr_tmp;

        for (int64_t last_pos_min = -1; it_min != it_min_end; ++it_min){ // Iterate over minimizers of unitig to delete

            if ((last_pos_min < it_min.getPosition()) || isForbidden){ // If a new minimizer hash is found in unitig to delete

                minHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;
                isForbidden = false;

                while (it_it_min != it_it_min_end){ // Iterate over minimizers of current k-mer in unitig to delete

                    const minHashResult& min_h_res = *it_it_min;
                    Minimizer minz_rep = Minimizer(s + min_h_res.pos).rep(); // Get canonical minimizer
                    MinimizerIndex::iterator it_h = hmap_min_unitigs.find(minz_rep);

                    mhr = min_h_res;

                    while (it_h != hmap_min_unitigs.end()){ // If the minimizer is found

                        packed_tiny_vector& v = it_h.getVector();
                        uint8_t& flag_v = it_h.getVectorSize();
                        size_t i = 0, v_sz = v.size(flag_v);

                        while (i < v_sz){

                            if ((v(i, flag_v) & mask) == pos_id_unitig){

                                flag_v = v.remove(i, flag_v);
                                --v_sz;
                            }
                            else ++i;
                        }

                        if (v.size(flag_v) == 0) it_h = hmap_min_unitigs.end();
                        else if (!isShort && ((v(v_sz-1, flag_v) & mask) == mask)){ //Minimizer bin is overcrowded

                            mhr_tmp = it_min.getNewMin(mhr); //Recompute a new (different) minimizer for current k-mer
                            isForbidden = true;

                            if (mhr_tmp.hash != mhr.hash){

                                mhr = mhr_tmp;
                                minz_rep = Minimizer(s + mhr.pos).rep();
                                it_h = hmap_min_unitigs.find(minz_rep);
                            }
                            else {

                                it_h = hmap_min_unitigs.end();
                                break;
                            }
                        }
                        else it_h = hmap_min_unitigs.end();
                    }

                    last_pos_min = min_h_res.pos;
                    ++it_it_min;
                }
            }
        }
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, bool>::type CompactedDBG<U, G>::extractUnitig_(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs,
                                                                               size_t& v_unitigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp){

    bool deleted = true;

    if (!sp.empty()){

        const Unitig<U>* unitig = v_unitigs[pos_v_unitigs];

        //const pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        //const size_t totalcoverage = unitig->coveragesum - lowpair.second;
        //const size_t ccov_size = unitig->ccov.size();

        const string str = unitig->getSeq().toString();

        size_t i = 0;

        UnitigMap<U, G> um(pos_v_unitigs, 0, 0, unitig->length(), false, false, true, this);

        vector<Unitig<U>> v_data(sp.size());

        deleted = false;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit, ++i) { //Iterate over created split unitigs

            um.dist = sit->first;
            um.len = sit->second - um.dist;

            if (um.len == 1){

                const string split_str = str.substr(um.dist, um.len + k_ - 1);

                um.strand = (split_str <= reverse_complement(split_str));
            }
            else um.strand = true;

            v_data[i] = std::move(um.splitData(sit+1 == sp.end()));
        }

        i = 0;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit, ++i) { //Iterate over created split unitigs

            const size_t len = sit->second - sit->first;
            //const uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage
            const uint64_t cov_tmp = len * CompressedCoverage::getFullCoverage();

            string split_str = str.substr(sit->first, len + k_ - 1); // Split unitig sequence

            if (len == 1){

                const string split_str_rev = reverse_complement(split_str);

                if (split_str > split_str_rev) split_str = split_str_rev;

                if (addUnitig(split_str, v_kmers_sz)){

                    CompressedCoverage_t<U>& cc_t = *h_kmers_ccov.find(Kmer(split_str.c_str()).rep());

                    cc_t.ccov.setFull();

                    *(cc_t.getData()) = std::move(*(v_data[i].getData()));
                }
                else {

                    km_unitigs.setFull(v_kmers_sz); // We don't care about the coverage per k-mer anymore

                    *(km_unitigs.getData(v_kmers_sz)) = std::move(*(v_data[i].getData()));

                    ++v_kmers_sz;
                }
            }
            else {

                addUnitig(split_str, nxt_pos_insert_v_unitigs);

                //v_unitigs[nxt_pos_insert_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[nxt_pos_insert_v_unitigs]->getCov() = CompressedCoverage(v_unitigs[nxt_pos_insert_v_unitigs]->getSeq().size() - k_ + 1, true);
                //v_unitigs[nxt_pos_insert_v_unitigs]->coveragesum = cov_tmp;

                *(v_unitigs[nxt_pos_insert_v_unitigs]->getData()) = std::move(*(v_data[i].getData()));

                ++nxt_pos_insert_v_unitigs;
            }
        }
    }

    --nxt_pos_insert_v_unitigs; //Position of the last unitig in the vector which is not NULL

    if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

        swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

        // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
        // created in this method, do not try to split it again
        if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) ++pos_v_unitigs;
        else --v_unitigs_sz;
    }
    else --v_unitigs_sz;

    deleteUnitig_<false>(false, false, nxt_pos_insert_v_unitigs);

    return deleted;
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, bool>::type CompactedDBG<U, G>::extractUnitig_(size_t& pos_v_unitigs, size_t& nxt_pos_insert_v_unitigs,
                                                                              size_t& v_unitigs_sz, size_t& v_kmers_sz, const vector<pair<int,int>>& sp){

    bool deleted = true;

    if (!sp.empty()){

        const Unitig<U>* unitig = v_unitigs[pos_v_unitigs];

        //const pair<size_t, size_t> lowpair = unitig->ccov.lowCoverageInfo();

        //const size_t totalcoverage = unitig->coveragesum - lowpair.second;
        //const size_t ccov_size = unitig->ccov.size();

        const string str = unitig->getSeq().toString();

        deleted = false;

        for (vector<pair<int,int>>::const_iterator sit = sp.begin(); sit != sp.end(); ++sit) { //Iterate over created split unitigs

            const size_t len = sit->second - sit->first;
            //const uint64_t cov_tmp = (totalcoverage * len) / (ccov_size - lowpair.first); // Split unitig coverage
            const uint64_t cov_tmp = len * CompressedCoverage::getFullCoverage();
            const string split_str = str.substr(sit->first, len + k_ - 1); // Split unitig sequence

            if (len == 1){

                if (addUnitig(split_str, v_kmers_sz, pos_v_unitigs, false)) h_kmers_ccov.find(Kmer(split_str.c_str()).rep())->ccov.setFull();
                else km_unitigs.setFull(v_kmers_sz++); // We don't care about the coverage per k-mer anymore}
            }
            else {

                addUnitig(split_str, nxt_pos_insert_v_unitigs, pos_v_unitigs, false);

                //v_unitigs[nxt_pos_insert_v_unitigs]->initializeCoverage(true); //We don't care about the coverage per k-mer anymore
                v_unitigs[nxt_pos_insert_v_unitigs]->getCov() = CompressedCoverage(v_unitigs[nxt_pos_insert_v_unitigs]->getSeq().size() - k_ + 1, true);
                //v_unitigs[nxt_pos_insert_v_unitigs]->coveragesum = cov_tmp;

                ++nxt_pos_insert_v_unitigs;
            }
        }
    }

    --nxt_pos_insert_v_unitigs; //Position of the last unitig in the vector which is not NULL

    if (pos_v_unitigs != nxt_pos_insert_v_unitigs){ // Do not proceed to swap if swap positions are the same

        swapUnitigs(false, pos_v_unitigs, nxt_pos_insert_v_unitigs); // Swap unitigs

        // If the swapped unitig, previously in position nxt_pos_insert, was a split unitig
        // created in this method, do not try to split it again
        if (nxt_pos_insert_v_unitigs >= v_unitigs_sz) ++pos_v_unitigs;
        else --v_unitigs_sz;
    }
    else --v_unitigs_sz;

    deleteUnitig_<true>(false, false, nxt_pos_insert_v_unitigs);

    return deleted;
}

template<typename U, typename G>
UnitigMap<U, G> CompactedDBG<U, G>::find(const Kmer& km, const preAllocMinHashIterator<RepHash>& it_min_h) {

    const Kmer km_twin = km.twin();
    const Kmer& km_rep = km < km_twin ? km : km_twin;

    bool isShort;

    size_t unitig_id, unitig_id_pos_tmp, len;

    int64_t pos_match;

    const int diff = k_ - g_;

    preAllocMinHashIterator<RepHash> it_min(it_min_h, k_);
    preAllocMinHashResultIterator<RepHash> it_it_min = *it_min, it_it_min_end;

    minHashResult mhr, mhr_tmp;

    while (it_it_min != it_it_min_end){

        const minHashResult& min_h_res = *it_it_min;
        Minimizer minz = Minimizer(it_min.s + min_h_res.pos).rep();
        MinimizerIndex::const_iterator it = hmap_min_unitigs.find(minz);

        mhr = min_h_res;

        while (it != hmap_min_unitigs.end()){ // If the minimizer is found

            const packed_tiny_vector& v = it.getVector();
            const uint8_t flag_v = it.getVectorSize();
            const int v_sz = v.size(flag_v);

            it = hmap_min_unitigs.end();

            for (size_t i = 0; i < v_sz; ++i){

                unitig_id = v(i, flag_v) >> 32;

                if (unitig_id == RESERVED_ID){

                    if ((v(i, flag_v) & RESERVED_ID) != 0){ //This minimizer has abundant k-mers

                        typename h_kmers_ccov_t::const_iterator it_km = h_kmers_ccov.find(km_rep);

                        if (it_km != h_kmers_ccov.end()){

                            return UnitigMap<U, G>(it_km.getHash(), 0, 1, k_, false, true, km == km_rep, this);
                        }
                    }

                    if ((v(i, flag_v) & MASK_CONTIG_TYPE) == MASK_CONTIG_TYPE){ //This minimizer is unitig overcrowded

                        mhr_tmp = it_min.getNewMin(mhr);

                        if (mhr_tmp.hash != mhr.hash){

                            mhr = mhr_tmp;
                            minz = Minimizer(it_min.s + mhr.pos).rep();
                            it = hmap_min_unitigs.find(minz);
                        }
                    }
                }
                else {

                    isShort = (v(i, flag_v) & MASK_CONTIG_TYPE) != 0;
                    unitig_id_pos_tmp = v(i, flag_v) & MASK_CONTIG_POS;

                    if (isShort){

                        const Kmer km_unitig = km_unitigs.getKmer(unitig_id);

                        if (min_h_res.pos == unitig_id_pos_tmp){

                            if (km_unitig == km_rep){

                                return UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, true, this);
                            }
                        }
                        else if ((min_h_res.pos == diff - unitig_id_pos_tmp) && (km_unitig == km_rep)){

                            return UnitigMap<U, G>(unitig_id, 0, 1, k_, true, false, false, this);
                        }
                    }
                    else {

                        pos_match = unitig_id_pos_tmp - min_h_res.pos;
                        len = v_unitigs[unitig_id]->length() - k_;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km)){

                            return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, true, this);
                        }

                        pos_match = unitig_id_pos_tmp - diff + min_h_res.pos;

                        if ((pos_match >= 0) && (pos_match <= len) && v_unitigs[unitig_id]->getSeq().compareKmer(pos_match, k_, km_twin)){

                            return UnitigMap<U, G>(unitig_id, pos_match, 1, len + k_, false, false, false, this);
                        }
                    }
                }
            }
        }

        ++it_it_min;
    }

    return UnitigMap<U, G>();
}

// pre: Some k-mers in the unitigs might have a coverage which is less than CompressedCoverage::getFullCoverage()
//      Those k-mers must be deleted from the unitigs.
// post: All unitigs have a per k-mer coverage of CompressedCoverage::getFullCoverage(). The graph is not
//       necessarily compacted after calling this function.
template<typename U, typename G>
pair<size_t, size_t> CompactedDBG<U, G>::extractAllUnitigs() {

    size_t i;
    size_t split = 0, deleted = 0;
    size_t v_kmers_sz = km_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t nxt_pos_insert = v_unitigs.size();

	size_t delete1 = 0, delete2 = 0, delete3 = 0;

    for (typename h_kmers_ccov_t::iterator it(h_kmers_ccov.begin()); it != h_kmers_ccov.end(); ++it) {

        if (!it->ccov.isFull()){

            deleteUnitig_<is_void<U>::value>(false, true, it.getHash());
            ++deleted;
        }
    }

    for (i = 0; i < v_kmers_sz;) {

        if (!km_unitigs.isFull(i)) {

            --v_kmers_sz;

            if (i != v_kmers_sz) swapUnitigs(true, i, v_kmers_sz);

            deleteUnitig_<is_void<U>::value>(true, false, v_kmers_sz);

            ++deleted;
        }
        else ++i;
    }

    for (i = 0; i < v_unitigs_sz;) { // Iterate over unitigs created so far

        if (!v_unitigs[i]->getCov().isFull()) { //Coverage not full, unitig must be splitted

            vector<pair<int,int>> sp = v_unitigs[i]->getCov().splittingVector();

            if (extractUnitig_<is_void<U>::value>(i, nxt_pos_insert, v_unitigs_sz, v_kmers_sz, sp)) ++deleted;
            else {

                ++split;
                sp.clear();
            }
        }
        else ++i;
    }

    if (nxt_pos_insert < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert);
    if (v_kmers_sz < km_unitigs.size()) km_unitigs.resize(v_kmers_sz);

    return {split, deleted};
}

// pre: Some k-mers in the unitigs might have a coverage which is less than CompressedCoverage::getFullCoverage().
//      Unitigs must be split at the positions matching those low-coverage k-mers (not deleted).
// post: All unitigs have a per k-mer coverage of CompressedCoverage::getFullCoverage(). The graph is not
//       necessarily compacted after calling this function.
template<typename U, typename G>
pair<size_t, size_t> CompactedDBG<U, G>::splitAllUnitigs() {

    pair<size_t, size_t> p = {0, 0};

    size_t v_kmers_sz = km_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t nxt_pos_insert = v_unitigs.size();

    const size_t cov_full = CompressedCoverage::getFullCoverage();

    for (size_t i = 0; i < v_unitigs_sz;) { // Iterate over unitigs created so far

        const CompressedCoverage& ccov = v_unitigs[i]->getCov();

        if (!ccov.isFull()) { //Coverage not full, unitig must be splitted

            size_t prev_split_pos = 0;

            vector<pair<int,int>> sp;

            for (size_t pos = 0; pos < ccov.size(); ++pos){

                if ((ccov.covAt(pos) != cov_full) && (pos != prev_split_pos)){

                    sp.push_back({prev_split_pos, pos});
                    ++(p.second);

                    prev_split_pos = pos;
                }
            }

            sp.push_back({prev_split_pos, ccov.size()});

            ++(p.second);
            ++(p.first);

            extractUnitig_<is_void<U>::value>(i, nxt_pos_insert, v_unitigs_sz, v_kmers_sz, sp);
        }
        else ++i;
    }

    if (nxt_pos_insert < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert);
    if (v_kmers_sz < km_unitigs.size()) km_unitigs.resize(v_kmers_sz);

    return p;
}

template<typename U, typename G>
pair<size_t, size_t> CompactedDBG<U, G>::getSplitInfoAllUnitigs() const {

    pair<size_t, size_t> p = {0, 0};

    const size_t cov_full = CompressedCoverage::getFullCoverage();

    for (size_t i = 0; i < v_unitigs.size(); ++i) { // Iterate over unitigs created so far

        const CompressedCoverage& ccov = v_unitigs[i]->getCov();

        if (!ccov.isFull()) { //Coverage not full, unitig must be splitted

            size_t prev_split_pos = 0;

            for (size_t pos = 0; pos < ccov.size(); ++pos){

                if ((ccov.covAt(pos) != cov_full) && (pos != prev_split_pos)){

                    ++(p.second);

                    prev_split_pos = pos;
                }
            }

            ++(p.first);
            ++(p.second);
        }
    }

    return p;
}

template<typename U, typename G>
void CompactedDBG<U, G>::createJoinHT(vector<Kmer>* v_joins, KmerHashTable<Kmer>& joins, const size_t nb_threads) const {

    const size_t v_unitigs_size = v_unitigs.size();
    const size_t v_kmers_size = km_unitigs.size();
    const size_t chunk = 1024;

    if (v_joins == nullptr){

        if (nb_threads == 1){

            for (typename h_kmers_ccov_t::const_iterator it_ccov(h_kmers_ccov.begin()); it_ccov != h_kmers_ccov.end(); ++it_ccov) {

                const Kmer tail(it_ccov.getKey());
                const Kmer head_twin(tail.twin());

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(it_ccov.getHash(), 0, 1, k_, false, true, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
            }

            for (size_t i = 0; i != v_kmers_size; ++i) {

                const Kmer tail = km_unitigs.getKmer(i);
                const Kmer head_twin = tail.twin();

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(i, 0, 1, k_, true, false, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
            }

            for (size_t i = 0; i != v_unitigs_size; ++i) {

                const CompressedSequence& seq = v_unitigs[i]->getSeq();

                const Kmer head_twin(seq.getKmer(0).twin());
                const Kmer tail(seq.getKmer(seq.size() - k_));

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(i, 0, 1, seq.size(), false, false, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail);
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin);
            }
        }
        else {

            //Hybrid_SpinLockRW_MCS<> lck(nb_threads);
            SpinLockRW lck;

            auto worker_v_abundant = [chunk, &joins, &lck, this](typename h_kmers_ccov_t::const_iterator* l_it_ccov){

                typename h_kmers_ccov_t::const_iterator& it_ccov = *l_it_ccov;

                // We need to deal with the tail of long unitigs
                for (size_t i = 0; (it_ccov != h_kmers_ccov.end()) && (i < chunk); ++i, ++it_ccov) {

                    const Kmer tail(it_ccov.getKey());
                    const Kmer head_twin(tail.twin());

                    Kmer fw, bw;

                    const const_UnitigMap<U, G> cm(it_ccov.getHash(), 0, 1, k_, false, true, true, this);

                    lck.acquire_reader();

                    const bool joins_tail = (joins.find(tail) == joins.end());
                    const bool joins_head = (joins.find(head_twin) == joins.end());

                    lck.release_reader();

                    if (joins_tail && checkJoin(tail, cm, fw)){

                        lck.acquire_writer();
                        joins.insert(fw.twin(), tail);
                        lck.release_writer();
                    }

                    if (joins_head && checkJoin(head_twin, cm, bw)){

                        lck.acquire_writer();
                        joins.insert(bw.twin(), head_twin);
                        lck.release_writer();
                    }
                }
            };

            auto worker_v_kmers = [&joins, &lck, this](const size_t idx_a, const size_t idx_b){

                for (size_t i = idx_a; i != idx_b; ++i) {

                    Kmer fw, bw;

                    const Kmer tail = km_unitigs.getKmer(i);
                    const Kmer head_twin = tail.twin();

                    const const_UnitigMap<U, G> cm(i, 0, 1, k_, true, false, true, this);

                    lck.acquire_reader();

                    const bool joins_tail = (joins.find(tail) == joins.end());
                    const bool joins_head = (joins.find(head_twin) == joins.end());

                    lck.release_reader();

                    if (joins_tail && checkJoin(tail, cm, fw)){

                        lck.acquire_writer();
                        joins.insert(fw.twin(), tail);
                        lck.release_writer();
                    }

                    if (joins_head && checkJoin(head_twin, cm, bw)){

                        lck.acquire_writer();
                        joins.insert(bw.twin(), head_twin);
                        lck.release_writer();
                    }
                }
            };

            auto worker_v_unitigs = [&joins, &lck, this](   typename vector<Unitig<U>*>::const_iterator a,
                                                            typename vector<Unitig<U>*>::const_iterator b){

                for (size_t i = a - v_unitigs.begin(), end = b - v_unitigs.begin(); i != end; ++i) {

                    const CompressedSequence& seq = v_unitigs[i]->getSeq();

                    const const_UnitigMap<U, G> cm(i, 0, 1, seq.size(), false, false, true, this);

                    const Kmer head_twin(seq.getKmer(0).twin());
                    const Kmer tail(seq.getKmer(seq.size() - k_));

                    Kmer fw, bw;

                    lck.acquire_reader();

                    const bool joins_tail = (joins.find(tail) == joins.end());
                    const bool joins_head = (joins.find(head_twin) == joins.end());

                    lck.release_reader();

                    if (joins_tail && checkJoin(tail, cm, fw)){

                        lck.acquire_writer();
                        joins.insert(fw.twin(), tail);
                        lck.release_writer();
                    }

                    if (joins_head && checkJoin(head_twin, cm, bw)){

                        lck.acquire_writer();
                        joins.insert(bw.twin(), head_twin);
                        lck.release_writer();
                    }
                }
            };

            {
                typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(), it_end = h_kmers_ccov.end();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            typename h_kmers_ccov_t::const_iterator l_it;

                            bool stop;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it);

                                    l_it = it;

                                    for (size_t i = 0; (it != it_end) && (i < chunk); ++i, ++it){}

                                    stop = (l_it == it_end) && (it == it_end);
                                }

                                if (!stop) worker_v_abundant(&l_it);
                                else return;
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            {

                size_t it_km = 0, it_km_end = km_unitigs.size();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it_km;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            size_t l_it_km, l_it_km_end;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it_km);

                                    if (it_km == it_km_end) return;

                                    l_it_km = it_km;
                                    it_km = min(it_km + chunk, it_km_end);
                                    l_it_km_end = it_km;
                                }

                                worker_v_kmers(l_it_km, l_it_km_end);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            {
                auto it_unitig = v_unitigs.begin();
                auto it_unitig_end = v_unitigs.end();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it_unitig;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            auto l_it_unitig = v_unitigs.begin();
                            auto l_it_unitig_end = v_unitigs.end();

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it_unitig);

                                    if (it_unitig == it_unitig_end) return;

                                    l_it_unitig = it_unitig;
                                    l_it_unitig_end = it_unitig;

                                    if (distance(l_it_unitig, it_unitig_end) >= chunk){

                                        advance(l_it_unitig_end, chunk);

                                        it_unitig = l_it_unitig_end;
                                    }
                                    else {

                                        it_unitig = it_unitig_end;
                                        l_it_unitig_end = it_unitig_end;
                                    }
                                }

                                worker_v_unitigs(l_it_unitig, l_it_unitig_end);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }
        }
    }
    else {

        Kmer fw;

        for (auto km : *v_joins){

            const const_UnitigMap<U, G> cm = find(km, true);

            if (!cm.isEmpty){

                if (!cm.isShort && !cm.isAbundant){

                    if ((cm.dist == 0 && cm.strand) || (cm.dist != 0 && !cm.strand)) km = km.twin();
                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km);
                }
                else {

                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km);

                    km = km.twin();

                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km);
                }
            }
        }
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::createJoinHT(vector<Kmer>* v_joins, KmerHashTable<char>& joins, const size_t nb_threads) const {

    const size_t v_unitigs_size = v_unitigs.size();
    const size_t v_kmers_size = km_unitigs.size();

    const size_t chunk = 1024;

    if (v_joins == nullptr){

        if (nb_threads == 1){

            for (typename h_kmers_ccov_t::const_iterator it_ccov(h_kmers_ccov.begin()); it_ccov != h_kmers_ccov.end(); ++it_ccov) {

                const Kmer tail(it_ccov.getKey());
                const Kmer head_twin(tail.twin());

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(it_ccov.getHash(), 0, 1, k_, false, true, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail.getChar(0));
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin.getChar(0));
            }

            for (size_t i = 0; i != v_kmers_size; ++i) {

                const Kmer tail = km_unitigs.getKmer(i);
                const Kmer head_twin = tail.twin();

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(i, 0, 1, k_, true, false, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail.getChar(0));
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin.getChar(0));
            }

            for (size_t i = 0; i != v_unitigs_size; ++i) {

                const CompressedSequence& seq = v_unitigs[i]->getSeq();

                const Kmer head_twin(seq.getKmer(0).twin());
                const Kmer tail(seq.getKmer(seq.size() - k_));

                Kmer fw, bw;

                const const_UnitigMap<U, G> cm(i, 0, 1, seq.size(), false, false, true, this);

                if ((joins.find(tail) == joins.end()) && checkJoin(tail, cm, fw)) joins.insert(fw.twin(), tail.getChar(0));
                if ((joins.find(head_twin) == joins.end()) && checkJoin(head_twin, cm, bw)) joins.insert(bw.twin(), head_twin.getChar(0));
            }
        }
        else {

            //Hybrid_SpinLockRW_MCS<> lck(nb_threads);
            SpinLockRW lck;

            auto worker_v_abundant = [chunk, &joins, &lck, this](typename h_kmers_ccov_t::const_iterator* l_it_ccov){

                typename h_kmers_ccov_t::const_iterator& it_ccov = *l_it_ccov;

                // We need to deal with the tail of long unitigs
                for (size_t i = 0; (it_ccov != h_kmers_ccov.end()) && (i < chunk); ++i, ++it_ccov) {

                    const Kmer tail(it_ccov.getKey());
                    const Kmer head_twin(tail.twin());

                    Kmer fw, bw;

                    const const_UnitigMap<U, G> cm(it_ccov.getHash(), 0, 1, k_, false, true, true, this);

                    lck.acquire_reader();

                    const bool joins_tail = (joins.find(tail) == joins.end());
                    const bool joins_head = (joins.find(head_twin) == joins.end());

                    lck.release_reader();

                    if (joins_tail && checkJoin(tail, cm, fw)){

                        lck.acquire_writer();
                        joins.insert(fw.twin(), tail.getChar(0));
                        lck.release_writer();
                    }

                    if (joins_head && checkJoin(head_twin, cm, bw)){

                        lck.acquire_writer();
                        joins.insert(bw.twin(), head_twin.getChar(0));
                        lck.release_writer();
                    }
                }
            };

            auto worker_v_kmers = [&joins, &lck, this](const size_t idx_a, const size_t idx_b){

                for (size_t i = idx_a; i != idx_b; ++i) {

                    Kmer fw, bw;

                    const Kmer tail = km_unitigs.getKmer(i);
                    const Kmer head_twin = tail.twin();

                    const const_UnitigMap<U, G> cm(i, 0, 1, k_, true, false, true, this);

                    lck.acquire_reader();

                    const bool joins_tail = (joins.find(tail) == joins.end());
                    const bool joins_head = (joins.find(head_twin) == joins.end());

                    lck.release_reader();

                    if (joins_tail && checkJoin(tail, cm, fw)){

                        lck.acquire_writer();
                        joins.insert(fw.twin(), tail.getChar(0));
                        lck.release_writer();
                    }

                    if (joins_head && checkJoin(head_twin, cm, bw)){

                        lck.acquire_writer();
                        joins.insert(bw.twin(), head_twin.getChar(0));
                        lck.release_writer();
                    }
                }
            };

            auto worker_v_unitigs = [&joins, &lck, this](   typename vector<Unitig<U>*>::const_iterator a,
                                                            typename vector<Unitig<U>*>::const_iterator b){

                for (size_t i = a - v_unitigs.begin(), end = b - v_unitigs.begin(); i != end; ++i) {

                    const CompressedSequence& seq = v_unitigs[i]->getSeq();

                    const const_UnitigMap<U, G> cm(i, 0, 1, seq.size(), false, false, true, this);

                    const Kmer head_twin(seq.getKmer(0).twin());
                    const Kmer tail(seq.getKmer(seq.size() - k_));

                    Kmer fw, bw;

                    lck.acquire_reader();

                    const bool joins_tail = (joins.find(tail) == joins.end());
                    const bool joins_head = (joins.find(head_twin) == joins.end());

                    lck.release_reader();

                    if (joins_tail && checkJoin(tail, cm, fw)){

                        lck.acquire_writer();
                        joins.insert(fw.twin(), tail.getChar(0));
                        lck.release_writer();
                    }

                    if (joins_head && checkJoin(head_twin, cm, bw)){

                        lck.acquire_writer();
                        joins.insert(bw.twin(), head_twin.getChar(0));
                        lck.release_writer();
                    }
                }
            };

            {
                typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(), it_end = h_kmers_ccov.end();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            typename h_kmers_ccov_t::const_iterator l_it;

                            bool stop;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it);

                                    l_it = it;

                                    for (size_t i = 0; (it != it_end) && (i < chunk); ++i, ++it){}

                                    stop = (l_it == it_end) && (it == it_end);
                                }

                                if (!stop) worker_v_abundant(&l_it);
                                else return;
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            {

                size_t it_km = 0, it_km_end = km_unitigs.size();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it_km;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            size_t l_it_km, l_it_km_end;

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it_km);

                                    if (it_km == it_km_end) return;

                                    l_it_km = it_km;
                                    it_km = min(it_km + chunk, it_km_end);
                                    l_it_km_end = it_km;
                                }

                                worker_v_kmers(l_it_km, l_it_km_end);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }

            {
                auto it_unitig = v_unitigs.begin();
                auto it_unitig_end = v_unitigs.end();

                vector<thread> workers; // need to keep track of threads so we can join them

                mutex mutex_it_unitig;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            auto l_it_unitig = v_unitigs.begin();
                            auto l_it_unitig_end = v_unitigs.end();

                            while (true) {

                                {
                                    unique_lock<mutex> lock(mutex_it_unitig);

                                    if (it_unitig == it_unitig_end) return;

                                    l_it_unitig = it_unitig;
                                    l_it_unitig_end = it_unitig;

                                    if (distance(l_it_unitig, it_unitig_end) >= chunk){

                                        advance(l_it_unitig_end, chunk);

                                        it_unitig = l_it_unitig_end;
                                    }
                                    else {

                                        it_unitig = it_unitig_end;
                                        l_it_unitig_end = it_unitig_end;
                                    }
                                }

                                worker_v_unitigs(l_it_unitig, l_it_unitig_end);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }
        }
    }
    else {

        Kmer fw;

        for (auto km : *v_joins){

            const const_UnitigMap<U, G> cm = find(km, true);

            if (!cm.isEmpty){

                if (!cm.isShort && !cm.isAbundant){

                    if ((cm.dist == 0 && cm.strand) || (cm.dist != 0 && !cm.strand)) km = km.twin();
                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km.getChar(0));
                }
                else {

                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km.getChar(0));

                    km = km.twin();

                    if (checkJoin(km, cm, fw)) joins.insert(fw.twin(), km.getChar(0));
                }
            }
        }
    }
}

// use:  joined = mapper.joinUnitigs()
// pre:  no short unitigs exist in sUnitigs.
// post: all unitigs that could be connected have been connected
//       joined is the number of joined unitigs
template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, size_t>::type CompactedDBG<U, G>::joinUnitigs_(vector<Kmer>* v_joins, const size_t nb_threads) {

    size_t joined = 0;
    size_t cov_full = CompressedCoverage::getFullCoverage();
    size_t v_unitigs_size = v_unitigs.size();
    size_t v_kmers_size = km_unitigs.size();

    // a and b are candidates for joining
    /*KmerHashTable<Kmer> joins;

    createJoinHT(v_joins, joins, nb_threads);

    if (v_joins != nullptr) v_joins->clear();

    for (KmerHashTable<Kmer>::iterator it(joins.begin()); it != joins.end(); ++it) {

        const Kmer head(*it);
        const Kmer tail(it.getKey().twin());*/

    //-------------------------------------------------------------------------------
    KmerHashTable<char> joins;

    createJoinHT(v_joins, joins, nb_threads);

    if (v_joins != nullptr) v_joins->clear();

    for (KmerHashTable<char>::iterator it(joins.begin()); it != joins.end(); ++it) {

        const Kmer tail(it.getKey().twin());
        const Kmer head(tail.backwardBase(*it));
    //-------------------------------------------------------------------------------

        UnitigMap<U, G> cmHead(find(head, true));
        UnitigMap<U, G> cmTail(find(tail, true));

        if (!cmHead.isEmpty && !cmTail.isEmpty) {

            Kmer cmHead_head, cmTail_head;

            if (cmHead.isShort) cmHead_head = km_unitigs.getKmer(cmHead.pos_unitig);
            else if (cmHead.isAbundant) cmHead_head = h_kmers_ccov.find(cmHead.pos_unitig).getKey();
            else cmHead_head = v_unitigs[cmHead.pos_unitig]->getSeq().getKmer(0);

            if (cmTail.isShort) cmTail_head = km_unitigs.getKmer(cmTail.pos_unitig);
            else if (cmTail.isAbundant) cmTail_head = h_kmers_ccov.find(cmTail.pos_unitig).getKey();
            else cmTail_head = v_unitigs[cmTail.pos_unitig]->getSeq().getKmer(0);

            if (cmHead_head != cmTail_head) { // can't join a sequence with itself, either hairPin, loop or mobius loop

                // both kmers are still end-kmers
                bool headDir;
                bool len_k_head = cmHead.isShort || cmHead.isAbundant;

                if (len_k_head && (head == cmHead_head)) headDir = true;
                else if (!len_k_head && (head == v_unitigs[cmHead.pos_unitig]->getSeq().getKmer(v_unitigs[cmHead.pos_unitig]->numKmers()-1))) headDir = true;
                else if (head.twin() == cmHead_head) headDir = false;
                else continue;

                bool tailDir;
                bool len_k_tail = cmTail.isShort || cmTail.isAbundant;

                if (tail == cmTail_head) tailDir = true;
                else if (len_k_tail){
                    if (tail.twin() == cmTail_head) tailDir = false;
                    else continue;
                }
                else if (tail.twin() == v_unitigs[cmTail.pos_unitig]->getSeq().getKmer(v_unitigs[cmTail.pos_unitig]->numKmers()-1)) tailDir = false;
                else continue;

                //Compute join sequence
                string joinSeq;

                joinSeq.reserve((len_k_head ? k_ : cmHead.size) + (len_k_tail ? k_ : cmTail.size) - k_ + 1);

                if (headDir) joinSeq = len_k_head ? cmHead_head.toString() : v_unitigs[cmHead.pos_unitig]->getSeq().toString();
                else joinSeq = len_k_head ? cmHead_head.twin().toString() : v_unitigs[cmHead.pos_unitig]->getSeq().rev().toString();

                if (tailDir) joinSeq.append(len_k_tail ? cmTail_head.toString() : v_unitigs[cmTail.pos_unitig]->getSeq().toString(), k_ - 1, string::npos);
                else joinSeq.append(len_k_tail ? cmTail_head.twin().toString() : v_unitigs[cmTail.pos_unitig]->getSeq().rev().toString(), k_ - 1, string::npos);

                //Compute new coverage
                /*uint64_t covsum;

                if (len_k_head){

                    CompressedCoverage& ccov = cmHead.isShort ? v_kmers[cmHead.pos_unitig].second.ccov : h_kmers_ccov.find(cmHead.pos_unitig)->ccov;
                    covsum = ccov.covAt(0);
                }
                else covsum = v_unitigs[cmHead.pos_unitig]->coveragesum;

                if (len_k_tail){

                    CompressedCoverage& ccov = cmTail.isShort ? v_kmers[cmTail.pos_unitig].second.ccov : h_kmers_ccov.find(cmTail.pos_unitig)->ccov;
                    covsum += ccov.covAt(0);
                }
                else covsum += v_unitigs[cmTail.pos_unitig]->coveragesum;*/

                Unitig<U> data_tmp; //Store temporarily the new merged data
                Unitig<U>* unitig; //New unitig

                cmHead.strand = headDir;
                cmHead.dist = 0;
                cmHead.len = cmHead.size - k_ + 1;

                cmTail.strand = tailDir;
                cmTail.dist = 0;
                cmTail.len = cmTail.size - k_ + 1;

                data_tmp.getData()->concat(cmHead, cmTail);

                cmHead.getData()->clear(cmHead);
                cmTail.getData()->clear(cmTail);

                if (cmHead.isShort || cmHead.isAbundant){

                    if (cmHead.isShort){ //If head is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmHead.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmHead.pos_unitig, v_kmers_size);

                            // If the last unitig of the vector used for the swap was the tail
                            if (cmTail.isShort && (v_kmers_size == cmTail.pos_unitig)) cmTail.pos_unitig = cmHead.pos_unitig;
                        }

                        deleteUnitig_<false>(true, false, v_kmers_size);
                    }
                    else if (cmHead.isAbundant) deleteUnitig_<false>(false, true, cmHead.pos_unitig);
                }

                if (cmTail.isShort || cmTail.isAbundant){

                    if (cmTail.isShort){ //If tail is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmTail.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmTail.pos_unitig, v_kmers_size);

                            if (cmHead.isShort && (v_kmers_size == cmHead.pos_unitig)) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig_<false>(true, false, v_kmers_size);
                    }
                    else if (cmTail.isAbundant) deleteUnitig_<false>(false, true, cmTail.pos_unitig);
                }

                if (len_k_head && len_k_tail){

                    addUnitig(joinSeq, v_unitigs_size);
                    unitig = v_unitigs[v_unitigs_size];
                    ++v_unitigs_size;
                }
                else if (len_k_head){

                    deleteUnitig_<false>(false, false, cmTail.pos_unitig);
                    addUnitig(joinSeq, cmTail.pos_unitig);
                    unitig = v_unitigs[cmTail.pos_unitig];
                }
                else {

                    if (!len_k_tail){

                        --v_unitigs_size;

                        if (cmTail.pos_unitig != v_unitigs_size){

                            swapUnitigs(false, cmTail.pos_unitig, v_unitigs_size);

                            if (v_unitigs_size == cmHead.pos_unitig) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig_<false>(false, false, v_unitigs_size);
                    }

                    deleteUnitig_<false>(false, false, cmHead.pos_unitig);
                    addUnitig(joinSeq, cmHead.pos_unitig);
                    unitig = v_unitigs[cmHead.pos_unitig];
                }

                //unitig->coveragesum = covsum;
                //if (covsum >= cov_full * unitig->numKmers()) unitig->ccov.setFull();
                unitig->getCov().setFull();

                *(unitig->getData()) = std::move(*(data_tmp.getData()));

                ++joined;
            }
        }
    }

    if (v_unitigs_size < v_unitigs.size()) v_unitigs.resize(v_unitigs_size);
    if (v_kmers_size < km_unitigs.size()) km_unitigs.resize(v_kmers_size);

    return joined;
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, size_t>::type CompactedDBG<U, G>::joinUnitigs_(vector<Kmer>* v_joins, const size_t nb_threads) {

    size_t joined = 0;
    size_t cov_full = CompressedCoverage::getFullCoverage();
    size_t v_unitigs_size = v_unitigs.size();
    size_t v_kmers_size = km_unitigs.size();

    // a and b are candidates for joining
    /*KmerHashTable<Kmer> joins;

    createJoinHT(v_joins, joins, nb_threads);

    if (v_joins != nullptr) v_joins->clear();

    for (KmerHashTable<Kmer>::iterator it = joins.begin(); it != joins.end(); ++it) {

        const Kmer head = *it;
        const Kmer tail = it.getKey().twin();*/

    //-------------------------------------------------------------------------------
    KmerHashTable<char> joins;

    createJoinHT(v_joins, joins, nb_threads);

    if (v_joins != nullptr) v_joins->clear();

    for (KmerHashTable<char>::iterator it(joins.begin()); it != joins.end(); ++it) {

        const Kmer tail(it.getKey().twin());
        const Kmer head(tail.backwardBase(*it));
    //-------------------------------------------------------------------------------

        UnitigMap<U, G> cmHead = find(head, true);
        UnitigMap<U, G> cmTail = find(tail, true);

        if (!cmHead.isEmpty && !cmTail.isEmpty) {

            Kmer cmHead_head, cmTail_head;

            if (cmHead.isShort) cmHead_head = km_unitigs.getKmer(cmHead.pos_unitig);
            else if (cmHead.isAbundant) cmHead_head = h_kmers_ccov.find(cmHead.pos_unitig).getKey();
            else cmHead_head = v_unitigs[cmHead.pos_unitig]->getSeq().getKmer(0);

            if (cmTail.isShort) cmTail_head = km_unitigs.getKmer(cmTail.pos_unitig);
            else if (cmTail.isAbundant) cmTail_head = h_kmers_ccov.find(cmTail.pos_unitig).getKey();
            else cmTail_head = v_unitigs[cmTail.pos_unitig]->getSeq().getKmer(0);

            if (cmHead_head != cmTail_head) { // can't join a sequence with itself, either hairPin, loop or mobius loop

                // both kmers are still end-kmers
                bool headDir;
                bool len_k_head = cmHead.isShort || cmHead.isAbundant;

                if (len_k_head && (head == cmHead_head)) headDir = true;
                else if (!len_k_head && (head == v_unitigs[cmHead.pos_unitig]->getSeq().getKmer(v_unitigs[cmHead.pos_unitig]->numKmers()-1))) headDir = true;
                else if (head.twin() == cmHead_head) headDir = false;
                else continue;

                bool tailDir;
                bool len_k_tail = cmTail.isShort || cmTail.isAbundant;

                if (tail == cmTail_head) tailDir = true;
                else if (len_k_tail){
                    if (tail.twin() == cmTail_head) tailDir = false;
                    else continue;
                }
                else if (tail.twin() == v_unitigs[cmTail.pos_unitig]->getSeq().getKmer(v_unitigs[cmTail.pos_unitig]->numKmers()-1)) tailDir = false;
                else continue;

                //Compute join sequence
                string joinSeq;

                joinSeq.reserve((len_k_head ? k_ : cmHead.size) + (len_k_tail ? k_ : cmTail.size) - k_ + 1);

                if (headDir) joinSeq = len_k_head ? cmHead_head.toString() : v_unitigs[cmHead.pos_unitig]->getSeq().toString();
                else joinSeq = len_k_head ? cmHead_head.twin().toString() : v_unitigs[cmHead.pos_unitig]->getSeq().rev().toString();

                if (tailDir) joinSeq.append(len_k_tail ? cmTail_head.toString() : v_unitigs[cmTail.pos_unitig]->getSeq().toString(), k_ - 1, string::npos);
                else joinSeq.append(len_k_tail ? cmTail_head.twin().toString() : v_unitigs[cmTail.pos_unitig]->getSeq().rev().toString(), k_ - 1, string::npos);

                //Compute new coverage
                /*uint64_t covsum;

                if (len_k_head){

                    CompressedCoverage& ccov = cmHead.isShort ? v_kmers[cmHead.pos_unitig].second.ccov : h_kmers_ccov.find(cmHead.pos_unitig)->ccov;
                    covsum = ccov.covAt(0);
                }
                else covsum = v_unitigs[cmHead.pos_unitig]->coveragesum;

                if (len_k_tail){

                    CompressedCoverage& ccov = cmTail.isShort ? v_kmers[cmTail.pos_unitig].second.ccov : h_kmers_ccov.find(cmTail.pos_unitig)->ccov;
                    covsum += ccov.covAt(0);
                }
                else covsum += v_unitigs[cmTail.pos_unitig]->coveragesum;*/

                Unitig<U>* unitig; //New unitig

                cmTail.strand = tailDir;
                cmHead.strand = headDir;

                if (cmHead.isShort || cmHead.isAbundant){

                    if (cmHead.isShort){ //If head is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmHead.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmHead.pos_unitig, v_kmers_size);

                            // If the last unitig of the vector used for the swap was the tail
                            if (cmTail.isShort && (v_kmers_size == cmTail.pos_unitig)) cmTail.pos_unitig = cmHead.pos_unitig;
                        }

                        deleteUnitig_<true>(true, false, v_kmers_size);
                    }
                    else if (cmHead.isAbundant) deleteUnitig_<true>(false, true, cmHead.pos_unitig);
                }

                if (cmTail.isShort || cmTail.isAbundant){

                    if (cmTail.isShort){ //If tail is a short unitig, swap and delete it

                        --v_kmers_size;

                        if (cmTail.pos_unitig != v_kmers_size){

                            swapUnitigs(true, cmTail.pos_unitig, v_kmers_size);

                            if (cmHead.isShort && (v_kmers_size == cmHead.pos_unitig)) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig_<true>(true, false, v_kmers_size);
                    }
                    else if (cmTail.isAbundant) deleteUnitig_<true>(false, true, cmTail.pos_unitig);
                }

                if (len_k_head && len_k_tail){

                    addUnitig(joinSeq, v_unitigs_size);
                    unitig = v_unitigs[v_unitigs_size];
                    ++v_unitigs_size;
                }
                else if (len_k_head){

                    deleteUnitig_<true>(false, false, cmTail.pos_unitig);
                    addUnitig(joinSeq, cmTail.pos_unitig);
                    unitig = v_unitigs[cmTail.pos_unitig];
                }
                else {

                    if (!len_k_tail){

                        --v_unitigs_size;

                        if (cmTail.pos_unitig != v_unitigs_size){

                            swapUnitigs(false, cmTail.pos_unitig, v_unitigs_size);

                            if (v_unitigs_size == cmHead.pos_unitig) cmHead.pos_unitig = cmTail.pos_unitig;
                        }

                        deleteUnitig_<true>(false, false, v_unitigs_size);
                    }

                    deleteUnitig_<true>(false, false, cmHead.pos_unitig);
                    addUnitig(joinSeq, cmHead.pos_unitig);
                    unitig = v_unitigs[cmHead.pos_unitig];
                }

                //unitig->coveragesum = covsum;
                //if (covsum >= cov_full * unitig->numKmers()) unitig->ccov.setFull();
                unitig->getCov().setFull();

                ++joined;
            }
        }
    }

    if (v_unitigs_size < v_unitigs.size()) v_unitigs.resize(v_unitigs_size);
    if (v_kmers_size < km_unitigs.size()) km_unitigs.resize(v_kmers_size);

    return joined;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::checkJoin(const Kmer& a, const const_UnitigMap<U, G>& cm_a, Kmer& b) const {

    size_t i, j, count_succ;

    vector<const_UnitigMap<U, G>> v_um(findSuccessors(a, 2, true));

    for (i = 0, count_succ = 0; i != 4; ++i){

        if (!v_um[i].isEmpty){

            ++count_succ;
            j = i;
        }
    }

    if (count_succ == 1) {

        Kmer cand_head, ac_head;

        const Kmer fw_cand(a.forwardBase(alpha[j]));

        const const_UnitigMap<U, G> cm_cand(v_um[j]);

        if (cm_cand.isShort) cand_head = km_unitigs.getKmer(cm_cand.pos_unitig);
        else if (cm_cand.isAbundant) cand_head = h_kmers_ccov.find(cm_cand.pos_unitig).getKey();
        else cand_head = v_unitigs[cm_cand.pos_unitig]->getSeq().getKmer(0);

        if (cm_a.isShort) ac_head = km_unitigs.getKmer(cm_a.pos_unitig);
        else if (cm_a.isAbundant) ac_head = h_kmers_ccov.find(cm_a.pos_unitig).getKey();
        else ac_head = v_unitigs[cm_a.pos_unitig]->getSeq().getKmer(0);

        if (cand_head != ac_head) {

            v_um = findSuccessors(fw_cand.twin(), 2, true);

            for (i = 0, count_succ = 0; i != 4; ++i) count_succ += !v_um[i].isEmpty;

            if (count_succ == 1) {

                b = fw_cand;
                return true;
            }
        }
    }

    return false;
}

template<typename U, typename G>
void CompactedDBG<U, G>::check_fp_tips(KmerHashTable<bool>& ignored_km_tips){

    uint64_t nb_real_short_tips = 0;

    size_t nxt_pos_insert_v_unitigs = v_unitigs.size();
    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = km_unitigs.size();

    vector<pair<int,int>> sp;

    for (KmerHashTable<bool>::iterator it(ignored_km_tips.begin()); it != ignored_km_tips.end(); ++it) {

        Kmer km(it.getKey());

        UnitigMap<U, G> cm(find(km, true)); // Check if the (short) tip actually exists

        if (!cm.isEmpty){ // IF the tip exists

            ++nb_real_short_tips;

            bool not_found = true;

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                UnitigMap<U, G> cm_bw(find(km.backwardBase(alpha[i])));

                if (!cm_bw.isEmpty && !cm_bw.isAbundant && !cm_bw.isShort){

                    //if (cm_bw.strand) ++cm_bw.dist;
                    cm_bw.dist += cm_bw.strand;

                    if ((cm_bw.dist != 0) && (cm_bw.dist != cm_bw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_bw.dist));
                        sp.push_back(make_pair(cm_bw.dist, cm_bw.size - k_ + 1));

                        extractUnitig_<is_void<U>::value>(cm_bw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                        sp.clear();
                    }

                    not_found = false;
                }
            }

            for (size_t i = 0; (i < 4) && not_found; ++i) {

                UnitigMap<U, G> cm_fw(find(km.forwardBase(alpha[i])));

                if (!cm_fw.isEmpty && !cm_fw.isAbundant && !cm_fw.isShort){

                    //if (!cm_fw.strand) ++cm_fw.dist;
                    cm_fw.dist += !cm_fw.strand;

                    if ((cm_fw.dist != 0) && (cm_fw.dist != cm_fw.size - k_ + 1)){

                        sp.push_back(make_pair(0, cm_fw.dist));
                        sp.push_back(make_pair(cm_fw.dist, cm_fw.size - k_ + 1));

                        extractUnitig_<is_void<U>::value>(cm_fw.pos_unitig, nxt_pos_insert_v_unitigs, v_unitigs_sz, v_kmers_sz, sp);

                        sp.clear();
                    }

                    not_found = false;
                }
            }
        }
    }

    if (nxt_pos_insert_v_unitigs < v_unitigs.size()) v_unitigs.resize(nxt_pos_insert_v_unitigs);
    if (v_kmers_sz < km_unitigs.size()) km_unitigs.resize(v_kmers_sz);
}

template<typename U, typename G>
size_t CompactedDBG<U, G>::removeUnitigs(bool rmIsolated, bool clipTips, vector<Kmer>& v){

    if (!rmIsolated && !clipTips) return 0;

    const bool rm_and_clip = rmIsolated && clipTips;

    size_t v_unitigs_sz = v_unitigs.size();
    size_t v_kmers_sz = km_unitigs.size();
    size_t removed = 0;
    size_t i;

    int64_t j;

    const int lim = (clipTips ? 1 : 0);

    int nb_pred, nb_succ;

    Kmer km;

    Unitig<U>* unitig = nullptr;

    for (j = 0; j < v_unitigs_sz; ++j) {

        unitig = v_unitigs[j];

        if (unitig->numKmers() < k_){

            const Kmer head(unitig->getSeq().getKmer(0));

            nb_pred = 0;

            for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

                if (!find(head.backwardBase(alpha[i]), true).isEmpty){

                    ++nb_pred;
                    if (clipTips) km = head.backwardBase(alpha[i]);
                }
            }

            if (nb_pred <= lim){

                const Kmer tail(unitig->getSeq().getKmer(unitig->length() - k_));

                nb_succ = 0;

                for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                    if (!find(tail.forwardBase(alpha[i]), true).isEmpty){

                        ++nb_succ;
                        if (clipTips) km = tail.forwardBase(alpha[i]);
                    }
                }

                if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                    ++removed;
                    --v_unitigs_sz;

                    if (j != v_unitigs_sz){

                        swapUnitigs(false, j, v_unitigs_sz),
                        --j;
                    }

                    if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
                }
            }
        }
    }

    for (j = 0; j < v_kmers_sz; ++j) {

        const Kmer km_unitig = km_unitigs.getKmer(j);

        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            if (!find(km_unitig.backwardBase(alpha[i]), true).isEmpty){

                ++nb_pred;

                if (clipTips) km = km_unitig.backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                if (!find(km_unitig.forwardBase(alpha[i]), true).isEmpty){

                    ++nb_succ;

                    if (clipTips) km = km_unitig.forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))) { //Unitig is isolated

                ++removed;
                --v_kmers_sz;

                if (j != v_kmers_sz){

                    swapUnitigs(true, j, v_kmers_sz),
                    --j;
                }

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (typename h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it) {

        nb_pred = 0;

        for (i = 0; (i != 4) && (nb_pred <= lim); ++i) {

            if (!find(it.getKey().backwardBase(alpha[i]), true).isEmpty){

                ++nb_pred;
                if (clipTips) km = it.getKey().backwardBase(alpha[i]);
            }
        }

        if (nb_pred <= lim){

            nb_succ = 0;

            for (i = 0; (i != 4) && (nb_succ <= lim); ++i) {

                if (!find(it.getKey().forwardBase(alpha[i]), true).isEmpty){

                    ++nb_succ;
                    if (clipTips) km = it.getKey().forwardBase(alpha[i]);
                }
            }

            if ((rm_and_clip && ((nb_pred + nb_succ) <= lim)) || (!rm_and_clip && ((nb_pred + nb_succ) == lim))){

                ++removed;

                *it = CompressedCoverage_t<U>();

                if (clipTips && ((nb_pred + nb_succ) == lim)) v.push_back(km);
            }
        }
    }

    for (j = v_unitigs_sz; j < v_unitigs.size(); ++j) deleteUnitig_<is_void<U>::value>(false, false, j);
    v_unitigs.resize(v_unitigs_sz);

    for (j = v_kmers_sz; j < km_unitigs.size(); ++j) deleteUnitig_<is_void<U>::value>(true, false, j);
    km_unitigs.resize(v_kmers_sz);

    for (typename h_kmers_ccov_t::iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it){

        if (it->ccov.size() == 0) deleteUnitig_<is_void<U>::value>(false, true, it.getHash());
    }

    return removed;
}

template<typename U, typename G>
void CompactedDBG<U, G>::writeFASTA(const string& graphfilename) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = km_unitigs.size();

    size_t i = 0;

    ofstream graphfile;
    ostream graph(0);

    graphfile.open(graphfilename.c_str());
    graph.rdbuf(graphfile.rdbuf());
    //graph.sync_with_stdio(false);

    assert(!graphfile.fail());

    for (size_t j = 0; j < v_unitigs_sz; ++j, ++i) graph << ">" << i << "\n" << v_unitigs[j]->getSeq().toString() << "\n";
    for (size_t j = 0; j < v_kmers_sz; ++j, ++i) graph << ">" << i << "\n" << km_unitigs.getKmer(j).toString() << "\n";

    for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); it != h_kmers_ccov.end(); ++it, ++i) {

        graph << ">" << i << "\n" << it.getKey().toString() << "\n";
    }

    graphfile.close();
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    size_t labelA = 1;

    for (const auto& unitig : *this){

        const string seq(unitig.referenceUnitigToString());

        graph.write_sequence(std::to_string(labelA), seq.size(), seq, unitig.getData()->serialize(unitig));

        if (unitig.isAbundant) idmap.insert(Kmer(unitig.referenceUnitigToString().c_str()), labelA);

        ++labelA;
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    size_t labelA = 1;

    for (const auto& unitig : *this){

        const string seq(unitig.referenceUnitigToString());

        graph.write_sequence(std::to_string(labelA), seq.size(), seq, "");

        if (unitig.isAbundant) idmap.insert(Kmer(unitig.referenceUnitigToString().c_str()), labelA);

        ++labelA;
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::writeGFA(const string& graphfilename, const size_t nb_threads) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = km_unitigs.size();

    size_t labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    const string header_tag("BV:Z:" + string(BFG_VERSION) + "\t" + "KL:Z:" + to_string(k_) + "\t" + "ML:Z:" + to_string(g_));

    KmerHashTable<size_t> idmap(h_kmers_ccov.size());

    GFA_Parser graph(graphfilename);

    graph.open_write(1, header_tag);

    writeGFA_sequence_<is_void<U>::value>(graph, idmap);

    if (nb_threads == 1){

        for (labelA = 1; labelA <= v_unitigs_sz; labelA++) {

            const Unitig<U>* unitig = v_unitigs[labelA - 1];
            const Kmer head = unitig->getSeq().getKmer(0);
            const Kmer tail = unitig->getSeq().getKmer(unitig->length() - k_);

            const vector<const_UnitigMap<U, G>> pred = findPredecessors(head, true);
            const vector<const_UnitigMap<U, G>> succ = findSuccessors(tail, 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }

        for (labelA = v_unitigs_sz + 1; labelA <= v_kmers_sz + v_unitigs_sz; labelA++) {

            const Kmer km_unitig = km_unitigs.getKmer(labelA - v_unitigs_sz - 1);

            const vector<const_UnitigMap<U, G>> pred = findPredecessors(km_unitig, true);
            const vector<const_UnitigMap<U, G>> succ = findSuccessors(km_unitig, 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }

        for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

            labelA = *it;

            const vector<const_UnitigMap<U, G>> pred = findPredecessors(it.getKey(), true);
            const vector<const_UnitigMap<U, G>> succ = findSuccessors(it.getKey(), 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const string slabelA = std::to_string(labelA);
                    const string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }
    }
    else {

        const size_t chunk_size = 1024;

        auto worker_v_unitigs = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Unitig<U>* unitig = v_unitigs[labelA - 1];

                const Kmer head = unitig->getSeq().getKmer(0);
                const Kmer tail = unitig->getSeq().getKmer(unitig->length() - k_);

                const vector<const_UnitigMap<U, G>> pred = this->findPredecessors(head, true);
                const vector<const_UnitigMap<U, G>> succ = this->findSuccessors(tail, 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_kmers = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Kmer km_unitig = km_unitigs.getKmer(labelA - v_unitigs_sz - 1);

                const vector<const_UnitigMap<U, G>> pred = this->findPredecessors(km_unitig, true);
                const vector<const_UnitigMap<U, G>> succ = this->findSuccessors(km_unitig, 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, false), make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(labelA, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_abundant = [v_unitigs_sz, chunk_size, &idmap, this](  KmerHashTable<size_t>::iterator* l_it,
                                                                            vector<pair<pair<size_t, bool>, pair<size_t, bool>>>* v_out){

            KmerHashTable<size_t>::iterator& it = *l_it;

            // We need to deal with the tail of long unitigs
            for (size_t i = 0; (it != idmap.end()) && (i < chunk_size); ++i, ++it) {

                const vector<const_UnitigMap<U, G>> pred = this->findPredecessors(it.getKey(), true);
                const vector<const_UnitigMap<U, G>> succ = this->findSuccessors(it.getKey(), 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(*it, false), make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(make_pair(make_pair(*it, true), make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        {
            atomic<size_t> label(1);

            vector<vector<pair<pair<size_t, bool>, pair<size_t, bool>>>> v_out(nb_threads);
            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            const size_t old_labelA = label.fetch_add(chunk_size);

                            if (old_labelA <= v_unitigs_sz){

                                if (old_labelA + chunk_size <= v_unitigs_sz) worker_v_unitigs(old_labelA, old_labelA + chunk_size, &v_out[t]);
                                else worker_v_unitigs(old_labelA, v_unitigs_sz + 1, &v_out[t]);

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const string slabelA = std::to_string(p.first.first);
                                        const string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        {
            const size_t v_kmers_unitigs_sz = v_kmers_sz + v_unitigs_sz;

            atomic<size_t> label(v_unitigs_sz + 1);

            vector<vector<pair<pair<size_t, bool>, pair<size_t, bool>>>> v_out(nb_threads);
            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            const size_t old_labelA = label.fetch_add(chunk_size);

                            if (old_labelA <= v_kmers_unitigs_sz){

                                if (old_labelA + chunk_size <= v_kmers_unitigs_sz) worker_v_kmers(old_labelA, old_labelA + chunk_size, &v_out[t]);
                                else worker_v_kmers(old_labelA, v_kmers_unitigs_sz + 1, &v_out[t]);

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const string slabelA = std::to_string(p.first.first);
                                        const string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        {
            KmerHashTable<size_t>::iterator it = idmap.begin(), it_end = idmap.end();

            vector<vector<pair<pair<size_t, bool>, pair<size_t, bool>>>> v_out(nb_threads);
            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file, mutex_it;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        KmerHashTable<size_t>::iterator l_it;

                        bool stop;

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_it);

                                l_it = it;

                                for (size_t i = 0; (it != it_end) && (i < chunk_size); ++i, ++it){}

                                stop = (l_it == it_end) && (it == it_end);
                            }

                            if (!stop){

                                worker_v_abundant(&l_it, &v_out[t]);

                                {
                                    unique_lock<mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const string slabelA = std::to_string(p.first.first);
                                        const string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }
    }

    graph.close();
}

template<typename U, typename G>
void CompactedDBG<U, G>::readGFA(const string& graphfilename, const size_t nb_threads) {

    size_t graph_file_id = 0;

    bool new_file_opened = false;

    GFA_Parser graph(graphfilename);

    graph.open_read();

    GFA_Parser::GFA_line r = graph.read(graph_file_id, new_file_opened, true);

    if (nb_threads == 1){

        while ((r.first != nullptr) || (r.second != nullptr)){

            if (r.first != nullptr) addUnitig(r.first->seq, (r.first->seq.length() == k_) ? km_unitigs.size() : v_unitigs.size());

            r = graph.read(graph_file_id, new_file_opened, true);
        }
    }
    else {

        const size_t block_sz = 1024;

        std::atomic<size_t> v_kmers_sz;
        std::atomic<size_t> v_unitigs_sz;

        bool is_first = true;
        bool stop = false;

        SpinLock lck_unitig, lck_kmer;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        v_kmers_sz = 0;
        v_unitigs_sz = 0;

        hmap_min_unitigs.init_threads();

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    vector<string> seq;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) return;

                            seq.clear();

                            for (size_t i = 0; (i < block_sz) && !stop; ++i){

                                if (!is_first) r = graph.read(graph_file_id, new_file_opened, true);
                                if (r.first != nullptr) seq.push_back(r.first->seq);

                                stop = ((r.first == nullptr) && (r.second == nullptr));
                                is_first = false;
                            }
                        }

                        for (const auto& s : seq) addUnitig(s, (s.length() == k_) ? v_kmers_sz++ : v_unitigs_sz++, lck_unitig, lck_kmer);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        hmap_min_unitigs.release_threads();
        moveToAbundant();
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::readFASTA(const string& graphfilename, const size_t nb_threads) {

    size_t graph_file_id = 0;

    FastqFile ff(vector<string>(1, graphfilename));

    string seq;

    if (nb_threads == 1){

        while (ff.read_next(seq, graph_file_id) != -1) addUnitig(seq, (seq.length() == k_) ? km_unitigs.size() : v_unitigs.size());
    }
    else {

        const size_t block_sz = 1024;

        bool stop = false;

        std::atomic<size_t> v_kmers_sz;
        std::atomic<size_t> v_unitigs_sz;

        SpinLock lck_unitig, lck_kmer;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        v_kmers_sz = 0;
        v_unitigs_sz = 0;

        hmap_min_unitigs.init_threads();

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    vector<string> v_seq;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) return;

                            v_seq.clear();

                            for (size_t i = 0; (i < block_sz) && !stop; ++i){

                                stop = (ff.read_next(seq, graph_file_id) == -1);

                                if (!stop && !seq.empty()) v_seq.push_back(seq);
                            }
                        }

                        for (const auto& s : v_seq) addUnitig(s, (s.length() == k_) ? v_kmers_sz++ : v_unitigs_sz++, lck_unitig, lck_kmer);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        hmap_min_unitigs.release_threads();
        moveToAbundant();
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::mapRead(const const_UnitigMap<U, G>& um) {

    if (um.isEmpty) return; // nothing maps, move on

    if (um.isShort) km_unitigs.cover(um.pos_unitig);
    else if (um.isAbundant) h_kmers_ccov.find(um.pos_unitig)->ccov.cover(um.dist, um.dist + um.len - 1);
    else v_unitigs[um.pos_unitig]->getCov().cover(um.dist, um.dist + um.len - 1);
}

template<typename U, typename G>
void CompactedDBG<U, G>::mapRead(const const_UnitigMap<U, G>& um, LockGraph& lck_g) {

    if (um.isEmpty) return; // nothing maps, move on

    if (um.isShort) km_unitigs.cover_thread_safe(um.pos_unitig);
    else {

        size_t lock_unitig_id = um.pos_unitig;

        lock_unitig_id += v_unitigs.size() & (static_cast<size_t>(!um.isShort) - 1);
        lock_unitig_id += (v_unitigs.size() + km_unitigs.size()) & (static_cast<size_t>(!um.isAbundant) - 1);

        lck_g.lock_unitig(lock_unitig_id);

        if (um.isAbundant) h_kmers_ccov.find(um.pos_unitig)->ccov.cover(um.dist, um.dist + um.len - 1);
        else v_unitigs[um.pos_unitig]->getCov().cover(um.dist, um.dist + um.len - 1);

        lck_g.unlock_unitig(lock_unitig_id);
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::unmapRead(const const_UnitigMap<U, G>& um) {

    if (um.isEmpty) return; // nothing maps, move on

    if (um.isShort) km_unitigs.uncover(um.pos_unitig);
    else if (um.isAbundant) h_kmers_ccov.find(um.pos_unitig)->ccov.uncover(um.dist, um.dist + um.len - 1);
    else v_unitigs[um.pos_unitig]->getCov().uncover(um.dist, um.dist + um.len - 1);
}

template<typename U, typename G>
void CompactedDBG<U, G>::unmapRead(const const_UnitigMap<U, G>& um, LockGraph& lck_g) {

    if (um.isEmpty) return; // nothing maps, move on

    if (um.isShort) km_unitigs.uncover_thread_safe(um.pos_unitig);
    else {

        size_t lock_unitig_id = um.pos_unitig;

        lock_unitig_id += v_unitigs.size() & (static_cast<size_t>(!um.isShort) - 1);
        lock_unitig_id += (v_unitigs.size() + km_unitigs.size()) & (static_cast<size_t>(!um.isAbundant) - 1);

        lck_g.lock_unitig(lock_unitig_id);

        if (um.isAbundant) h_kmers_ccov.find(um.pos_unitig)->ccov.uncover(um.dist, um.dist + um.len - 1);
        else v_unitigs[um.pos_unitig]->getCov().uncover(um.dist, um.dist + um.len - 1);

        lck_g.unlock_unitig(lock_unitig_id);
    }
}

template<typename U, typename G>
vector<Kmer> CompactedDBG<U, G>::extractMercyKmers(BlockedBloomFilter& bf_uniq_km, const size_t nb_threads, const bool verbose) {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = km_unitigs.size();

    size_t i, j;

    char km_tmp[MAX_KMER_SIZE];

    KmerHashTable<uint8_t> tips;

    vector<Kmer> v_out;

    for (typename h_kmers_ccov_t::iterator it_ccov = h_kmers_ccov.begin(); it_ccov != h_kmers_ccov.end(); ++it_ccov) {

        const Kmer km = it_ccov.getKey().rep();

        vector<UnitigMap<U, G>> v_um = findPredecessors(km, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        v_um = findSuccessors(km, true);

        for (j = 0; (j != 4) && v_um[j].isEmpty; ++j){}

        if ((i == 4) && (j == 4)) tips.insert(km, 3);
        else if (j == 4) tips.insert(km, 2);
        else if (i == 4) tips.insert(km, 1);
    }

    for (size_t it_v_km = 0; it_v_km != v_kmers_sz; ++it_v_km) {

        const Kmer km = km_unitigs.getKmer(it_v_km).rep();

        vector<UnitigMap<U, G>> v_um = findPredecessors(km, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        v_um = findSuccessors(km, true);

        for (j = 0; (j != 4) && v_um[j].isEmpty; ++j){}

        if ((i == 4) && (j == 4)) tips.insert(km, 3);
        else if (j == 4) tips.insert(km, 2);
        else if (i == 4) tips.insert(km, 1);
    }

    for (size_t it_v_unitig = 0; it_v_unitig != v_unitigs_sz; ++it_v_unitig) {

        const CompressedSequence& seq = v_unitigs[it_v_unitig]->getSeq();

        const Kmer head = seq.getKmer(0);
        const Kmer tail = seq.getKmer(seq.size() - k_);

        vector<UnitigMap<U, G>> v_um = findPredecessors(head, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        if (i == 4){

            const Kmer head_rep = head.rep();

            if (head == head_rep) tips.insert(head, 1);
            else tips.insert(head_rep, 2);
        }

        v_um = findSuccessors(tail, true);

        for (i = 0; (i != 4) && v_um[i].isEmpty; ++i){}

        if (i == 4){

            const Kmer tail_rep = tail.rep();

            if (tail == tail_rep) tips.insert(tail, 2);
            else tips.insert(tail_rep, 1);
        }
    }

    for (typename KmerHashTable<uint8_t>::iterator it_a = tips.begin(); it_a != tips.end(); ++it_a) {

        if ((*it_a == 1) || (*it_a == 3)){ // Corresponding k-mer has no predecessor in the graph

            bool pres_neigh_bw[4] = {false, false, false, false};
            uint64_t hashes_bw[4];

            const Kmer km_a = it_a.getKey();

            km_a.toString(km_tmp); // Get the k-mer
            RepHash rep_h(k_), rep_h_cpy;

            rep_h.init(km_tmp);

            for (i = 0; i != 4; ++i) {

                rep_h_cpy = rep_h;
                rep_h_cpy.updateBW(km_tmp[k_ - 1], alpha[i]);

                hashes_bw[i] = rep_h_cpy.hash(); // Prepare the hash of its predecessor
            }

            std::memmove(km_tmp + 1, km_tmp, (k_ - 1) * sizeof(char));

            // Query the MBBF for all possible predecessors
            bf_uniq_km.contains(hashes_bw, minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash(), pres_neigh_bw, 4);

            for (i = 0; i != 4; ++i) {

                if (pres_neigh_bw[i]){ // A potential mercy k-mer has been found

                    const Kmer km_mercy = km_a.backwardBase(alpha[i]);

                    for (j = 0; j != 4; ++j) {

                        const Kmer km_b = km_mercy.backwardBase(alpha[j]); // Possible predecessor of this mercy k-mer
                        const Kmer km_b_rep = km_b.rep();

                        typename KmerHashTable<uint8_t>::iterator it_b = tips.find(km_b_rep); // Query the tips with predecessor

                        if (it_b != tips.end()){ // Possible predecessor of this mercy k-mer exists

                            if ((km_b == km_b_rep) && ((*it_b == 2) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 2 ? 0 : 1);
                                *it_a = (*it_a == 1 ? 0 : 2);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy); //Push the mercy k-mer found
                            }
                            else if ((km_b != km_b_rep) && ((*it_b == 1) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 1 ? 0 : 2);
                                *it_a = (*it_a == 1 ? 0 : 2);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy);  //Push the mercy k-mer found
                            }
                        }
                    }
                }
            }
        }

        if ((*it_a == 2) || (*it_a == 3)){ // Corresponding k-mer has no successor in the graph

            bool pres_neigh_fw[4] = {false, false, false, false};
            uint64_t hashes_fw[4];

            const Kmer km_a = it_a.getKey();

            km_a.toString(km_tmp); // Get the k-mer
            RepHash rep_h(k_), rep_h_cpy; // Prepare its hash

            rep_h.init(km_tmp);

            for (i = 0; i != 4; ++i) {

                rep_h_cpy = rep_h;
                rep_h_cpy.updateFW(km_tmp[0], alpha[i]);

                hashes_fw[i] = rep_h_cpy.hash(); // Prepare the hash of its successor
            }

            std::memmove(km_tmp, km_tmp + 1, (k_ - 1) * sizeof(char));

            // Query the MBBF for all possible predecessors
            bf_uniq_km.contains(hashes_fw, minHashKmer<RepHash>(km_tmp, k_, g_, RepHash(), true).getHash(), pres_neigh_fw, 4);

            for (i = 0; i != 4; ++i) {

                if (pres_neigh_fw[i]){ // A potential mercy k-mer has been found

                    const Kmer km_mercy = km_a.forwardBase(alpha[i]);

                    for (j = 0; j != 4; ++j) {

                        const Kmer km_b = km_mercy.forwardBase(alpha[j]); // Possible successor of this mercy k-mer
                        const Kmer km_b_rep = km_b.rep();

                        typename KmerHashTable<uint8_t>::iterator it_b = tips.find(km_b_rep); // Query the tips with successor

                        if (it_b != tips.end()){ // Possible successor of this mercy k-mer exists

                            if ((km_b == km_b_rep) && ((*it_b == 1) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 1 ? 0 : 2);
                                *it_a = (*it_a == 2 ? 0 : 1);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy); //Push the mercy k-mer found
                            }
                            else if ((km_b != km_b_rep) && ((*it_b == 2) || (*it_b == 3))){

                                // Those tips can't be use anymore in these directions
                                *it_b = (*it_b == 2 ? 0 : 1);
                                *it_a = (*it_a == 2 ? 0 : 1);
                                i = 3; j = 3; // Break the loops

                                v_out.push_back(km_mercy);  //Push the mercy k-mer found
                            }
                        }
                    }
                }
            }
        }
    }

    if (verbose) cout << "CompactedDBG::extractMercyKmers(): " << v_out.size() << " k-mers extracted" << endl;

    return v_out;
}

template<typename U, typename G>
size_t CompactedDBG<U, G>::joinTips(string filename_MBBF_uniq_kmers, const size_t nb_threads, const bool verbose) {

    if (invalid){

        cerr << "CompactedDBG::joinTips(): Graph is invalid and tips cannot be joined" << endl;
        return 0;
    }

    BlockedBloomFilter mbbf;

    FILE* f_mbbf;

    if ((f_mbbf = fopen(filename_MBBF_uniq_kmers.c_str(), "rb")) == NULL){

        cerr << "CompactedDBG::joinTips(): Minimizer Blocked Bloom filter file of unique k-mers cannot be opened" << endl;
        return 0;
    }

    mbbf.ReadBloomFilter(f_mbbf);
    fclose(f_mbbf);

    vector<Kmer> v_mercy_km = extractMercyKmers(mbbf, nb_threads, verbose);

    for (const auto& km_mercy : v_mercy_km) addUnitig(km_mercy.rep().toString(), km_unitigs.size());

    size_t nb_join = joinUnitigs_<is_void<U>::value>(&v_mercy_km, nb_threads);

    if (verbose) cout << "CompactedDBG<U, G>::joinTips(): " << nb_join << " unitigs have been joined using mercy k-mers" << endl;

    return nb_join;
}

template<typename U, typename G>
void CompactedDBG<U, G>::setKmerGmerLength(const int kmer_length, const int minimizer_length){

    invalid = false;

    if (kmer_length <= 2){

        cerr << "CompactedDBG::CompactedDBG(): Length k of k-mers cannot be less than 3" << endl;
        invalid = true;
    }

    if (kmer_length >= MAX_KMER_SIZE){

        cerr << "CompactedDBG::CompactedDBG(): Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        invalid = true;
    }

    if (minimizer_length == 0){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot be equal to 0" << endl;
        invalid = true;
    }

    if (minimizer_length >= MAX_GMER_SIZE){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        invalid = true;
    }

    if ((minimizer_length > 0) && (minimizer_length > kmer_length - 2)){

        cerr << "CompactedDBG::CompactedDBG(): Length g of minimizers cannot exceed k - 2" << endl;
        invalid = true;
    }

    if (!invalid){

        k_ = kmer_length;

        if (minimizer_length >= 0) g_ = minimizer_length;
        else if (kmer_length >= 15) g_ = k_ - DEFAULT_G_DEC1;
        else if (kmer_length >= 7) g_ = k_ - DEFAULT_G_DEC2;
        else g_ = k_ - 2;

        Kmer::set_k(k_);
        Minimizer::set_g(g_);
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::setFullCoverage(const size_t cov) const {

    CompressedCoverage::setFullCoverage(cov);
    KmerCovIndex<U>::setFullCoverage(cov);
}

template<typename U, typename G>
void CompactedDBG<U, G>::print() const {

    cout << "CompactedDBG::print(): v_unitigs.size() = " << v_unitigs.size() << endl;
    cout << "CompactedDBG::print(): v_kmers.size() = " << km_unitigs.size() << endl;
    cout << "CompactedDBG::print(): h_kmers_ccov.size() = " << h_kmers_ccov.size() << endl;
    cout << "CompactedDBG::print(): hmap_min_unitigs.size() = " << hmap_min_unitigs.size() << endl;
}

#endif
