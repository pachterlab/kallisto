#import "MetalUtils.h"
#import "MetalProcessReads.h"
#include "common.h"
#include "KmerIndex.h"
#include "PlaintextWriter.h"
#include "GPUIndexFormat.h"
#include "BenchmarkStats.h"

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <ctime>

// k must be 32 so each Kmer fits in a single uint64_t
static_assert(MAX_KMER_SIZE == 32, "MAX_KMER_SIZE must be 32 for Metal code");

using namespace std;

std::string argv_to_string(int argc, char *argv[]) {
    std::string res;
    for (int i = 0; i < argc; ++i) {
        res += argv[i];
        if (i + 1 < argc) res += " ";
    }
    return res;
}

std::string get_local_time() {
    time_t rawtime;
    struct tm* timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    std::string ret(asctime(timeinfo));
    return ret.substr(0, ret.size() - 1);
}

void usageEMMetal(bool valid_input = true) {
    if (valid_input) {
        cout << "metalkallisto " << KALLISTO_VERSION << endl
             << "Pseudoalignment pipeline on Apple Metal (M-series)" << endl << endl;
    }
    cout << "Usage: metalkallisto quant [arguments] FASTQ-files" << endl << endl
         << "Required arguments:" << endl
         << "-i, --index=STRING            Kallisto index file" << endl
         << "-o, --output-dir=STRING       Output directory" << endl << endl
         << "Optional arguments:" << endl
         << "-t, --threads=INT             Number of threads (default: 1)" << endl
         << "    --verbose                 Print progress every 1M reads" << endl
         << "    --single                  Single-end reads" << endl
         << "    --fr-stranded             Strand-specific reads (FR)" << endl
         << "    --rf-stranded             Strand-specific reads (RF)" << endl
         << "-l, --fragment-length=FLOAT   Estimated average fragment length (single-end)" << endl
         << "-s, --sd=FLOAT                Estimated standard deviation of fragment length" << endl
         << "-n, --iterations=INT          Number of EM iterations (default: 500)" << endl
         << "    --metallib=STRING         Path to kallisto.metallib (overrides compiled-in default)" << endl;
}

void ParseOptionsEMMetal(int argc, char **argv, ProgramOptions& opt, std::string& metallib_path) {
    int verbose_flag  = 0;
    int single_flag   = 0;
    int strand_FR_flag = 0;
    int strand_RF_flag = 0;

    const char *opt_string = "t:i:l:s:o:n:b:";
    static struct option long_options[] = {
        {"verbose",          no_argument,       &verbose_flag,   1},
        {"single",           no_argument,       &single_flag,    1},
        {"fr-stranded",      no_argument,       &strand_FR_flag, 1},
        {"rf-stranded",      no_argument,       &strand_RF_flag, 1},
        {"threads",          required_argument, 0, 't'},
        {"index",            required_argument, 0, 'i'},
        {"fragment-length",  required_argument, 0, 'l'},
        {"sd",               required_argument, 0, 's'},
        {"output-dir",       required_argument, 0, 'o'},
        {"iterations",       required_argument, 0, 'n'},
        {"bootstrap-samples",required_argument, 0, 'b'},
        {"metallib",         required_argument, 0, 'm'},
        {0, 0, 0, 0}
    };

    int c, option_index = 0;
    while (true) {
        c = getopt_long(argc, argv, opt_string, long_options, &option_index);
        if (c == -1) break;
        switch (c) {
        case 0: break;
        case 't': stringstream(optarg) >> opt.threads;   break;
        case 'i': opt.index  = optarg; break;
        case 'l': stringstream(optarg) >> opt.fld;       break;
        case 's': stringstream(optarg) >> opt.sd;        break;
        case 'o': opt.output = optarg; break;
        case 'n': stringstream(optarg) >> opt.iterations; break;
        case 'b': stringstream(optarg) >> opt.bootstrap; break;
        case 'm': metallib_path = optarg; break;
        default: break;
        }
    }

    for (int i = optind; i < argc; ++i)
        opt.files.push_back(argv[i]);

    if (verbose_flag)    opt.verbose = true;
    if (single_flag)     opt.single_end = true;
    if (strand_FR_flag) { opt.strand_specific = true; opt.strand = ProgramOptions::StrandType::FR; }
    if (strand_RF_flag) { opt.strand_specific = true; opt.strand = ProgramOptions::StrandType::RF; }
}

int main(int argc, char** argv) {
    std::cout.sync_with_stdio(false);
    setvbuf(stdout, NULL, _IOFBF, 1048576);

    if (argc < 2) {
        usageEMMetal();
        return 1;
    }

    std::string cmd(argv[1]);

    if (cmd == "quant") {
        if (argc == 2) { usageEMMetal(); return 0; }

        auto start_time = get_local_time();
        ProgramOptions opt;
        std::string metallib_path;

        ParseOptionsEMMetal(argc - 1, argv + 1, opt, metallib_path);

        if (opt.index.empty() || opt.output.empty()) {
            std::cerr << "Error: -i (index) and -o (output-dir) are required.\n";
            usageEMMetal(false);
            return 1;
        }

        // Initialize Metal context
        try {
            MetalContext::get().init(metallib_path);
        } catch (const std::exception& e) {
            std::cerr << "Error initializing Metal: " << e.what() << std::endl;
            return 1;
        }

        int64_t num_processed = 0, num_pseudoaligned = 0, num_unique = 0;
        size_t  num_trans = 0;
        size_t  index_version = 0;
        int     index_k = 31;
        MetalRunStats run_stats;

        if (is_gpu_index(opt.index)) {
            GPUIndex gpu_idx;
            auto t0 = std::chrono::high_resolution_clock::now();
            if (!gpu_idx.load(opt.index)) {
                std::cerr << "Error: failed to load GPU index " << opt.index << std::endl;
                return 1;
            }
            Kmer::set_k(gpu_idx.k);
            auto t1 = std::chrono::high_resolution_clock::now();
            g_benchmark_stats.setup_index_load_ms +=
                std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
            num_trans     = gpu_idx.num_transcripts;
            index_version = GPU_INDEX_FORMAT_VERSION;
            index_k       = (int)gpu_idx.k;
            metal_run(opt, gpu_idx, &run_stats);
        } else {
            auto t0 = std::chrono::high_resolution_clock::now();
            KmerIndex index(opt);
            index.load(opt, true, true, /*gpuMode=*/true);
            auto t1 = std::chrono::high_resolution_clock::now();
            g_benchmark_stats.setup_index_load_ms +=
                std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
            num_trans     = index.num_trans;
            index_version = index.INDEX_VERSION;
            index_k       = (int)index.k;
            metal_run(opt, index, &run_stats);
        }

        num_processed    = run_stats.num_processed;
        num_pseudoaligned = run_stats.num_pseudoaligned;
        num_unique        = run_stats.num_unique;

        std::string call = argv_to_string(argc, argv);
        plaintext_aux(
            opt.output + "/run_info.json",
            std::to_string(num_trans),
            std::to_string(0),
            std::to_string(num_processed),
            std::to_string(num_pseudoaligned),
            std::to_string(num_unique),
            KALLISTO_VERSION,
            std::to_string(index_version),
            std::to_string(index_k),
            start_time,
            call,
            "");

        std::cerr << "\n[quant] processed " << num_processed << " reads, "
                  << num_pseudoaligned << " reads pseudoaligned\n";
    } else {
        std::cerr << "Unknown command: " << cmd << std::endl;
        usageEMMetal();
        return 1;
    }

    return 0;
}
