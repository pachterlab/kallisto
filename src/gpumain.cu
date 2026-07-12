#include "common.h"
#include "KmerIndex.h"
#include "EMAlgorithm.h"
#include "PlaintextWriter.h"
#include "GPUProcessReads.cuh"
#include "GPUPipeline.cuh"
#include "GPUIndexFormat.h"
#include "GPUBusProcess.cuh"
#include "BusOptionsParser.h"
#include <sys/stat.h>
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

void usageBusGPU() {
  cout << "gpukallisto " << KALLISTO_VERSION << endl
       << "Generates BUS files for single-cell sequencing using GPU acceleration" << endl << endl
       << "Usage: gpukallisto bus [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index (.idx or .gpuidx)" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl
       << "-x, --technology=STRING       Single-cell technology used (e.g. 10xv3, dropseq) or"
       << " bc:umi:cdna" << endl << endl
       << "Optional arguments:" << endl
       << "-l, --list                    List the supported technologies and exit" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "    --unstranded              Disable strand-specific filtering (currently required for"
       << " v1)" << endl
       << "    --verbose                 Print verbose progress" << endl << endl
       << "Notes (v1 scope):" << endl
       << "  - Only 2-file presets/custom -x with one BC piece on file 0, one UMI piece on file 0,"
       << " and one cDNA piece on file 1 are supported." << endl
       << "  - --unstranded is required (presets that default to fr/rf-stranded must be passed"
       << " --unstranded)." << endl;
}

static bool file_exists(const std::string& fn) {
  struct stat st;
  return stat(fn.c_str(), &st) == 0;
}

bool CheckOptionsBusGPU(ProgramOptions& opt) {
  bool ret = true;

  if (opt.index.empty()) {
    cerr << "Error: kallisto index file missing" << endl;
    ret = false;
  } else if (!file_exists(opt.index)) {
    cerr << "Error: kallisto index file not found " << opt.index << endl;
    ret = false;
  }

  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory with -o" << endl;
    ret = false;
  } else {
    struct stat st;
    auto s = stat(opt.output.c_str(), &st);
    if (s == 0) {
      if (!S_ISDIR(st.st_mode)) {
        cerr << "Error: " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      }
    } else {
#ifdef _WIN32
      if (mkdir(opt.output.c_str()) != 0) {
#else
      if (mkdir(opt.output.c_str(), 0777) != 0) {
#endif
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  }

  if (opt.technology.empty()) {
    cerr << "Error: -x/--technology is required" << endl;
    ret = false;
    return ret;
  }

  // Reject unsupported v1 features early.
  if (opt.batch_mode) {
    cerr << "Error: -B/--batch is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.bam) {
    cerr << "Error: --bam is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.long_read) {
    cerr << "Error: --long is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.aa) {
    cerr << "Error: --aa is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.genomebam || opt.pseudobam) {
    cerr << "Error: --genomebam/--bam is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.input_interleaved_nfiles != 0) {
    cerr << "Error: --inleaved is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.record_batch_bus_barcode) {
    cerr << "Error: --batch-barcodes is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (!opt.tagsequence.empty()) {
    cerr << "Error: -T/--tag is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.num) {
    cerr << "Error: --num is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.max_num_reads != 0) {
    cerr << "Error: -N/--numReads is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.do_union) {
    cerr << "Error: --union is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  if (opt.no_jump) {
    cerr << "Error: --no-jump is not supported in gpukallisto bus (v1)" << endl;
    ret = false;
  }
  // Single end bus runs use file 0 for BC/UMI, file 1 for cDNA, so we treat
  // it as paired-end (two-file) input internally.
  opt.single_end = false;

  // Files: must be exactly 2 input fastq files (one batch).
  if (opt.files.empty()) {
    cerr << "Error: missing read files (expected exactly 2)" << endl;
    ret = false;
  } else if (opt.files.size() != 2) {
    cerr << "Error: gpukallisto bus (v1) requires exactly 2 input files (got "
         << opt.files.size() << ")" << endl;
    ret = false;
  } else {
    for (const auto& fn : opt.files) {
      if (!file_exists(fn)) {
        cerr << "Error: file not found " << fn << endl;
        ret = false;
      }
    }
  }
  if (!ret) return ret;

  // Apply preset technology / parse custom -x.
  auto& busopt = opt.busOptions;
  busopt = BUSOptions{};
  busopt.nfiles = 1;
  busopt.keep_fastq_comments = false;
  busopt.paired = false;
  busopt.long_read = false;
  busopt.unmapped = false;
  busopt.error_rate = 0.0;
  busopt.threshold = 0.8;
  busopt.aa = false;

  std::vector<std::string> errs;
  ProgramOptions::StrandType preset_strand = ProgramOptions::StrandType::None;
  bool ok = bus_parser::ApplyBusPresetTechnology(opt, preset_strand, errs);
  if (!ok) {
    for (const auto& e : errs) cerr << e << endl;
    return false;
  }

  // v1 strand: --unstranded must be set (or be the technology's default).
  bool unstranded = (opt.strand == ProgramOptions::StrandType::None);
  if (!unstranded) {
    cerr << "Error: gpukallisto bus (v1) only supports --unstranded; pass"
         << " --unstranded to override the technology default" << endl;
    ret = false;
  }
  if (preset_strand != ProgramOptions::StrandType::None && !unstranded) {
    cerr << "Error: technology " << opt.technology
         << " defaults to a stranded mode; --unstranded is required for v1" << endl;
    ret = false;
  }
  // Force unstranded behavior for downstream bookkeeping.
  opt.strand_specific = false;
  opt.strand = ProgramOptions::StrandType::None;

  // v1 shape: 2 files, 1 BC piece on file 0, 1 UMI piece on file 0, 1 cDNA piece on file 1.
  if (busopt.nfiles != 2) {
    cerr << "Error: gpukallisto bus (v1) requires a 2-file technology (got nfiles="
         << busopt.nfiles << ")" << endl;
    ret = false;
  }
  if (busopt.bc.size() != 1 || busopt.bc[0].fileno != 0 || busopt.bc[0].start < 0 ||
      busopt.bc[0].stop <= busopt.bc[0].start) {
    cerr << "Error: gpukallisto bus (v1) requires exactly one BC piece on file 0" << endl;
    ret = false;
  }
  if (busopt.umi.size() != 1 || busopt.umi[0].fileno != 0 || busopt.umi[0].start < 0 ||
      busopt.umi[0].stop <= busopt.umi[0].start) {
    cerr << "Error: gpukallisto bus (v1) requires exactly one UMI piece on file 0" << endl;
    ret = false;
  }
  if (busopt.seq.size() != 1 || busopt.seq[0].fileno != 1) {
    cerr << "Error: gpukallisto bus (v1) requires exactly one cDNA piece on file 1" << endl;
    ret = false;
  }
  if (busopt.paired) {
    cerr << "Error: gpukallisto bus (v1) does not support paired cDNA reads" << endl;
    ret = false;
  }

  if (opt.threads <= 0) opt.threads = 1;

  opt.bus_mode = true;

  return ret;
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

    if (cmd == "bus") {
      if (argc == 2) {
        usageBusGPU();
        return 0;
      }
      bus_parser::ParseOptionsBus(argc - 1, argv + 1, opt);
      if (!CheckOptionsBusGPU(opt)) {
        cerr << endl;
        usageBusGPU();
        return 1;
      }
      GPURunStats run_stats;
      if (is_gpu_index(opt.index)) {
        GPUIndex gpu_idx;
        auto setup_start = std::chrono::high_resolution_clock::now();
        if (!gpu_idx.load(opt.index)) {
          cerr << "Error: failed to load GPU index " << opt.index << endl;
          return 1;
        }
        Kmer::set_k(gpu_idx.k);
        auto setup_end = std::chrono::high_resolution_clock::now();
        g_benchmark_stats.setup_index_load_ms +=
            std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start).count() / 1000.0;
        gpu_bus_run(opt, gpu_idx, start_time, argc, argv, &run_stats);
      } else {
        auto setup_start = std::chrono::high_resolution_clock::now();
        KmerIndex index(opt);
        index.load(opt, true, true, /*gpuMode=*/true);
        auto setup_end = std::chrono::high_resolution_clock::now();
        g_benchmark_stats.setup_index_load_ms +=
            std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start).count() / 1000.0;
        gpu_bus_run(opt, index, start_time, argc, argv, &run_stats);
      }
      return 0;
    } else if (cmd == "convert-index") {
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
