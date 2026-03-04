#include "common.h"
#include "KmerIndex.h"
#include "EMAlgorithm.h"
#include "PlaintextWriter.h"
#include "GPUProcessReads.cuh"
#include "GPUPipeline.cuh"
#include "GPUIndexFormat.h"
// #include "Kmer.hpp"

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <ctime>

// KSEQ_INIT(gzFile, gzread)

// Ensure MAX_KMER_SIZE is 32 so each Kmer fits in a single uint64_t
static_assert(MAX_KMER_SIZE == 32, "MAX_KMER_SIZE must be 32 for GPU code");

using namespace std;


std::string argv_to_string(int argc, char *argv[]) {
  std::string res;
  for (int i = 0; i < argc; ++i) {
    res += argv[i];
    if (i + 1 < argc) {
      res += " ";
    }
  }

  return res;
}

std::string get_local_time() {
  time_t rawtime;
  struct tm * timeinfo;

  time( &rawtime );
  timeinfo = localtime( &rawtime );
  std::string ret(asctime(timeinfo));

  // chomp off the newline
  return ret.substr(0, ret.size() - 1);
}

void usageEMGPU(bool valid_input = true) {
  if (valid_input) {

  cout << "gpukallisto " << KALLISTO_VERSION << endl
       << "Computes equivalence classes for reads and quantifies abundances" << endl << endl;
  }
  //      "----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|"
  cout << "Usage: gpukallisto quant [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              quantification" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "    --verbose                 Print out progress information every 1M proccessed reads" << endl
       << "    --build-gpu-index         When using kallisto index, also build .gpuidx for faster future runs" << endl;
}

void usageConvertIndex() {
  cout << "Usage: gpukallisto convert-index -i <kallisto-index> -o <gpu-index>" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Input kallisto index (.idx)" << endl
       << "-o, --output=STRING            Output GPU index (.gpuidx)" << endl << endl
       << "Converts a kallisto index to the lightweight GPU index format for faster loading." << endl;
}


bool CheckOptionsEMGPU(ProgramOptions& opt, bool emonly = false) {

  return true;
}


void ParseOptionsEMGPU(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int plaintext_flag = 0;
  int write_index_flag = 0;
  int single_flag = 0;
  int single_overhang_flag = 0;
  int strand_FR_flag = 0;
  int strand_RF_flag = 0;
  int bias_flag = 0;
  int pbam_flag = 0;
  int gbam_flag = 0;
  int fusion_flag = 0;
  int build_gpu_index_flag = 0;
  int em_only_flag = 0;

  const char *opt_string = "t:i:l:s:o:n:m:d:b:g:c:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"build-gpu-index", no_argument, &build_gpu_index_flag, 1},
    {"em-only", no_argument, &em_only_flag, 1},
    {"plaintext", no_argument, &plaintext_flag, 1},
    {"write-index", no_argument, &write_index_flag, 1},
    {"single", no_argument, &single_flag, 1},
    {"single-overhang", no_argument, &single_overhang_flag, 1},
    {"fr-stranded", no_argument, &strand_FR_flag, 1},
    {"rf-stranded", no_argument, &strand_RF_flag, 1},
    {"bias", no_argument, &bias_flag, 1},
    {"pseudobam", no_argument, &pbam_flag, 1},
    {"genomebam", no_argument, &gbam_flag, 1},
    {"fusion", no_argument, &fusion_flag, 1},
    {"seed", required_argument, 0, 'd'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"index", required_argument, 0, 'i'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {"bootstrap-samples", required_argument, 0, 'b'},
    {"gtf", required_argument, 0, 'g'},
    {"chromosomes", required_argument, 0, 'c'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'l': {
      stringstream(optarg) >> opt.fld;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.sd;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    case 'n': {
      stringstream(optarg) >> opt.iterations;
      break;
    }
    case 'b': {
      stringstream(optarg) >> opt.bootstrap;
      break;
    }
    case 'g': {
      stringstream(optarg) >> opt.gtfFile;
      break;
    }
    case 'c': {
      stringstream(optarg) >> opt.chromFile;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.seed;
      break;
    }

    default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }

  if (verbose_flag) {
    std::cout << "Verbose mode enabled" << std::endl;
    opt.verbose = true;
  }

  if (plaintext_flag) {
    opt.plaintext = true;
  }

  if (write_index_flag) {
    opt.write_index = true;
  }

  if (single_flag) {
    opt.single_end = true;
  }

  if (single_overhang_flag) {
    opt.single_overhang = true;
  }

  if (strand_FR_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::FR;
  }
  
  if (strand_RF_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::RF;
  }

  if (bias_flag) {
    opt.bias = true;
  }

  if (pbam_flag) {
    opt.pseudobam = true;
  }

  if (gbam_flag) {
    opt.pseudobam = true;
    opt.genomebam = true;    
  }

  if (fusion_flag) {
    opt.fusion = true;
  }

  if (build_gpu_index_flag) {
    opt.build_gpu_index = true;
  }

  if (em_only_flag) {
    opt.em_only = true;
  }
}


int main(int argc, char** argv) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);


  if (argc < 2) {
    usageEMGPU();
    exit(1);
  } else {
    auto start_time(get_local_time());
    ProgramOptions opt;
    string cmd(argv[1]);

    if (cmd == "convert-index") {
      if (argc < 6) {
        usageConvertIndex();
        return 1;
      }
      string kallisto_idx, gpuidx_out;
      const char* opt_string = "i:o:";
      static struct option long_options[] = {
        {"index", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {0, 0, 0, 0}
      };
      int c;
      int option_index = 0;
      while ((c = getopt_long(argc - 1, argv + 1, opt_string, long_options, &option_index)) != -1) {
        if (c == 'i') kallisto_idx = optarg;
        else if (c == 'o') gpuidx_out = optarg;
      }
      if (kallisto_idx.empty() || gpuidx_out.empty()) {
        cerr << "Error: -i and -o are required for convert-index" << endl;
        usageConvertIndex();
        return 1;
      }
      if (!convert_kallisto_to_gpu_index(kallisto_idx, gpuidx_out, opt)) {
        return 1;
      }
      return 0;
    } else if (cmd == "quant") {
      if (argc==2) {
        usageEMGPU();
        return 0;
      }
      ParseOptionsEMGPU(argc-1,argv+1,opt);
      if (!CheckOptionsEMGPU(opt)) {
        cerr << endl;
        usageEMGPU(false);
        exit(1);
      } else {
        int64_t num_processed = 0;
        int64_t num_pseudoaligned = 0;
        int64_t num_unique = 0;
        size_t num_trans = 0;
        size_t index_version = 0;
        int index_k = 31;
        GPURunStats run_stats;

        if (is_gpu_index(opt.index)) {
          GPUIndex gpu_idx;
          auto setup_start = std::chrono::high_resolution_clock::now();
          if (!gpu_idx.load(opt.index)) {
            cerr << "Error: failed to load GPU index " << opt.index << endl;
            return 1;
          }
          Kmer::set_k(gpu_idx.k);  // Required for read k-mer extraction (DeviceKmerLoader::run uses Kmer::k)
          auto setup_end = std::chrono::high_resolution_clock::now();
          auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
          g_benchmark_stats.setup_index_load_ms += setup_duration.count() / 1000.0;
          num_trans = gpu_idx.num_transcripts;
          index_version = GPU_INDEX_FORMAT_VERSION;
          index_k = static_cast<int>(gpu_idx.k);
          gpu_run(opt, gpu_idx, &run_stats);
          num_processed = run_stats.num_processed;
          num_pseudoaligned = run_stats.num_pseudoaligned;
          num_unique = run_stats.num_unique;
        } else {
          auto setup_start = std::chrono::high_resolution_clock::now();
          KmerIndex index(opt);
          index.load(opt, true, true, /*gpuMode=*/true);
          auto setup_end = std::chrono::high_resolution_clock::now();
          auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);
          g_benchmark_stats.setup_index_load_ms += setup_duration.count() / 1000.0;
          num_trans = index.num_trans;
          index_version = index.INDEX_VERSION;
          index_k = static_cast<int>(index.k);
          gpu_run(opt, index, &run_stats);
          num_processed = run_stats.num_processed;
          num_pseudoaligned = run_stats.num_pseudoaligned;
          num_unique = run_stats.num_unique;

          if (opt.build_gpu_index) {
            std::string gpuidx_path = opt.index + ".gpuidx";
            std::cerr << "[index] Building GPU index: " << gpuidx_path << std::endl;
            GPUIndex gpu_idx;
            gpu_idx.k = index.k;
            gpu_idx.num_transcripts = index.num_trans;
            gpu_idx.num_ecs = index.ecmapinv.size();
            gpu_idx.num_contigs = index.dbg.size();
            gpu_idx.target_names_ = index.target_names_;
            gpu_idx.target_lens_ = index.target_lens_;
            gpu_idx.ecmapinv = index.ecmapinv;
            std::vector<SparseVector<uint32_t>> vals;
            for (const auto& contig : index.dbg) {
              auto n = contig.getData();
              gpu_idx.contig_sequences.push_back(contig.referenceUnitigToString());
              n->ec.get_vals(vals);
              std::vector<ECBlock> blocks;
              int j = 0;
              size_t contigpos = 0;
              while (contigpos < contig.len) {
                auto mc = n->ec.get_block_at(contigpos);
                const auto& val = vals[j];
                const Roaring& trs = val.getIndices();
                auto ec_it = index.ecmapinv.find(trs);
                if (ec_it != index.ecmapinv.end()) {
                  blocks.push_back({static_cast<uint32_t>(mc.first),
                                   static_cast<uint32_t>(mc.second),
                                   ec_it->second});
                }
                contigpos = mc.second;
                ++j;
              }
              gpu_idx.contig_ec_blocks.push_back(std::move(blocks));
            }
            gpu_idx.write(gpuidx_path);
          }
        }

        std::string call = argv_to_string(argc, argv);
        plaintext_aux(
            opt.output + "/run_info.json",
            std::string(std::to_string(num_trans)),
            std::string(std::to_string(0)),
            std::string(std::to_string(num_processed)),
            std::string(std::to_string(num_pseudoaligned)),
            std::string(std::to_string(num_unique)),
            KALLISTO_VERSION,
            std::string(std::to_string(index_version)),
            std::string(std::to_string(index_k)),
            start_time,
            call,
            "");
      }
    } 
  }
}
