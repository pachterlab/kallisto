#include <getopt.h>
#include <sys/stat.h>
#include <time.h>

#include <CompactedDBG.hpp>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "Bootstrap.h"
#include "BusOptionsParser.h"
#include "EMAlgorithm.h"
#include "GeneModel.h"
#include "H5Writer.h"
#include "Inspect.h"
#include "KmerIndex.h"
#include "MinCollector.h"
#include "PlaintextWriter.h"
#include "ProcessReads.h"
#include "common.h"
#include "weights.h"

// #define ERROR_STR "\033[1mError:\033[0m"
#define ERROR_STR "Error:"


int my_mkdir(const char *path, mode_t mode) {
#ifdef _WIN64
  return mkdir(path);
#else
  return mkdir(path, mode);
#endif
}

bool checkFileExists(std::string fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}

void ParseOptionsIndex(int argc, char **argv, ProgramOptions &opt) {
  int verbose_flag = 0;
  int make_unique_flag = 0;
  int aa_flag = 0;
  int distinguish_flag = 0;
  int skip_index_flag = 0;
  const char *opt_string = "i:k:m:e:t:d:T:";
  static struct option long_options[] = {// long args
                                         {"verbose", no_argument, &verbose_flag, 1},
                                         {"make-unique", no_argument, &make_unique_flag, 1},
                                         {"aa", no_argument, &aa_flag, 1},
                                         {"skip-index", no_argument, &skip_index_flag, 1},
                                         {"distinguish", no_argument, &distinguish_flag, 1},
                                         // short args
                                         {"index", required_argument, 0, 'i'},
                                         {"kmer-size", required_argument, 0, 'k'},
                                         {"min-size", required_argument, 0, 'm'},
                                         {"ec-max-size", required_argument, 0, 'e'},
                                         {"threads", required_argument, 0, 't'},
                                         {"tmp", required_argument, 0, 'T'},
                                         {"d-list", required_argument, 0, 'd'},
                                         {"d-list-overhang", required_argument, 0,
                                          'D'},  // Do we have to have a one-letter flag as well?
                                         {0, 0, 0, 0}};
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 0:
        break;
      case 'i': {
        opt.index = optarg;
        break;
      }
      case 'k': {
        std::stringstream(optarg) >> opt.k;
        break;
      }
      case 'm': {
        std::stringstream(optarg) >> opt.g;
        break;
      }
      case 'e': {
        std::string max_size_str;
        std::stringstream(optarg) >> max_size_str;
        std::transform(max_size_str.begin(), max_size_str.end(), max_size_str.begin(), ::toupper);
        if (max_size_str == "AUTO") {
          opt.max_ec_size = 0;
        } else {
          std::stringstream(max_size_str) >> opt.max_ec_size;
        }
        break;
      }
      case 't': {
        std::stringstream(optarg) >> opt.threads;
        break;
      }
      case 'T': {
        std::stringstream(optarg) >> opt.tmp_dir;
        break;
      }
      case 'd': {
        std::string d_list;
        std::stringstream(optarg) >> d_list;
        std::stringstream ss(d_list);
        std::string filename;
        while (std::getline(ss, filename, ',')) {
          opt.d_list.push_back(filename);
        }
        break;
      }
      case 'D': {
        std::stringstream(optarg) >> opt.d_list_overhang;
        break;
      }
      default:
        break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
  if (make_unique_flag) {
    opt.make_unique = true;
  }
  if (aa_flag) {
    opt.aa = true;
    if (!opt.d_list.empty() && opt.d_list_overhang < 3) {
      // std::cerr << "[index] WARNING --d-list-overhang set to < 3; should be >= 3 with --aa" << std::endl;
      std::cerr << "[index] --d-list-overhang was set to 3 (with --aa, the d-list overhang must be >= 3)"
           << std::endl;
      opt.d_list_overhang = 3;
    }
  }
  if (distinguish_flag) {
    opt.distinguish = true;
  }

  for (int i = optind; i < argc; i++) {
    opt.transfasta.push_back(argv[i]);
  }
}

void ParseOptionsInspect(int argc, char **argv, ProgramOptions &opt) {
  const char *opt_string = "G:g:b:t:f:";

  int para_flag = 0;
  static struct option long_options[] = {// long args
                                         {"gfa", required_argument, 0, 'G'},
                                         {"gtf", required_argument, 0, 'g'},
                                         {"bed", required_argument, 0, 'b'},
                                         {"threads", required_argument, 0, 't'},
                                         {"fasta", required_argument, 0, 'f'},
                                         {"paranoid", no_argument, &para_flag, 1},
                                         {0, 0, 0, 0}};

  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 0:
        break;
      case 'G': {
        opt.gfa = optarg;
        break;
      }
      case 'b': {
        std::stringstream(optarg) >> opt.bedFile;
        break;
      }
      case 'g': {
        std::stringstream(optarg) >> opt.gtfFile;
        break;
      }
      case 't': {
        std::stringstream(optarg) >> opt.threads;
        if (opt.threads <= 0) {
          opt.threads = 1;
        }
        break;
      }
      case 'f': {
        opt.transcriptsFile = optarg;
        break;
      }
      default:
        break;
    }
  }
  opt.index = argv[optind];

  if (para_flag) {
    opt.inspect_thorough = true;
  }
}

void ParseOptionsEM(int argc, char **argv, ProgramOptions &opt) {
  int verbose_flag = 0;
  int plaintext_flag = 0;
  int write_index_flag = 0;
  int single_flag = 0;
  int long_read_flag = 0;
  int single_overhang_flag = 0;
  int strand_FR_flag = 0;
  int strand_RF_flag = 0;
  int bias_flag = 0;
  int pbam_flag = 0;
  int gbam_flag = 0;
  int fusion_flag = 0;
  int do_union_flag = 0;
  int no_jump_flag = 0;

  const char *opt_string = "t:i:l:P:s:o:n:m:d:b:g:c:p:";
  static struct option long_options[] = {// long args
                                         {"verbose", no_argument, &verbose_flag, 1},
                                         {"plaintext", no_argument, &plaintext_flag, 1},
                                         {"write-index", no_argument, &write_index_flag, 1},
                                         {"single", no_argument, &single_flag, 1},
                                         {"long", no_argument, &long_read_flag, 1},
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
                                         {"platform", required_argument, 0, 'P'},
                                         {"sd", required_argument, 0, 's'},
                                         {"output-dir", required_argument, 0, 'o'},
                                         {"iterations", required_argument, 0, 'n'},
                                         {"min-range", required_argument, 0, 'm'},
                                         {"bootstrap-samples", required_argument, 0, 'b'},
                                         {"gtf", required_argument, 0, 'g'},
                                         {"chromosomes", required_argument, 0, 'c'},
                                         {"priors", required_argument, 0, 'p'},
                                         {"union", no_argument, &do_union_flag, 1},
                                         {"no-jump", no_argument, &no_jump_flag, 1},
                                         {0, 0, 0, 0}};
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 0:
        break;
      case 't': {
        std::stringstream(optarg) >> opt.threads;
        break;
      }
      case 'i': {
        opt.index = optarg;
        break;
      }
      case 'l': {
        std::stringstream(optarg) >> opt.fld;
        break;
      }
      case 'P': {
        std::stringstream(optarg) >> opt.platform;
        break;
      }
      case 's': {
        std::stringstream(optarg) >> opt.sd;
        break;
      }
      case 'o': {
        opt.output = optarg;
        break;
      }
      case 'n': {
        std::stringstream(optarg) >> opt.iterations;
        break;
      }
      case 'm': {
        std::stringstream(optarg) >> opt.min_range;
      }
      case 'b': {
        std::stringstream(optarg) >> opt.bootstrap;
        break;
      }
      case 'g': {
        std::stringstream(optarg) >> opt.gtfFile;
        break;
      }
      case 'c': {
        std::stringstream(optarg) >> opt.chromFile;
        break;
      }
      case 'd': {
        std::stringstream(optarg) >> opt.seed;
        break;
      }
      case 'p': {
        opt.priors = optarg;
        break;
      }

      default:
        break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }

  if (verbose_flag) {
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

  if (long_read_flag) {
    opt.long_read = true;
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

  if (do_union_flag) {
    opt.do_union = true;
  }

  if (no_jump_flag) {
    opt.no_jump = true;
  }
}

void ParseOptionsTCCQuant(int argc, char **argv, ProgramOptions &opt) {
  const char *opt_string = "o:i:T:e:f:P:l:s:t:g:G:b:d:p:";
  int matrix_to_files = 0;
  int matrix_to_directories = 0;
  int plaintext_flag = 0;
  int long_read_flag = 0;
  static struct option long_options[] = {
      {"plaintext", no_argument, &plaintext_flag, 1},
      {"matrix-to-files", no_argument, &matrix_to_files, 1},
      {"matrix-to-directories", no_argument, &matrix_to_directories, 1},
      {"index", required_argument, 0, 'i'},
      {"txnames", required_argument, 0, 'T'},
      {"threads", required_argument, 0, 't'},
      {"fragment-file", required_argument, 0, 'f'},
      {"long", no_argument, &long_read_flag, 1},
      {"platform", required_argument, 0, 'P'},
      {"fragment-length", required_argument, 0, 'l'},
      {"sd", required_argument, 0, 's'},
      {"output-dir", required_argument, 0, 'o'},
      {"ec-file", required_argument, 0, 'e'},
      {"genemap", required_argument, 0, 'g'},
      {"gtf", required_argument, 0, 'G'},
      {"bootstrap-samples", required_argument, 0, 'b'},
      {"seed", required_argument, 0, 'd'},
      {"priors", required_argument, 0, 'p'},
      {0, 0, 0, 0}};
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 0:
        break;
      case 't': {
        std::stringstream(optarg) >> opt.threads;
        break;
      }
      case 'f': {
        std::stringstream(optarg) >> opt.fldFile;
        break;
      }
      case 'P': {
        std::stringstream(optarg) >> opt.platform;
        std::transform(opt.platform.begin(), opt.platform.end(), opt.platform.begin(), ::toupper);
        break;
      }
      case 'l': {
        std::stringstream(optarg) >> opt.fld;
        break;
      }
      case 's': {
        std::stringstream(optarg) >> opt.sd;
        break;
      }
      case 'o': {
        opt.output = optarg;
        break;
      }
      case 'i': {
        opt.index = optarg;
        break;
      }
      case 'e': {
        opt.ecFile = optarg;
        break;
      }
      case 'g': {
        opt.genemap = optarg;
        break;
      }
      case 'G': {
        opt.gtfFile = optarg;
        break;
      }
      case 'b': {
        std::stringstream(optarg) >> opt.bootstrap;
        break;
      }
      case 'd': {
        std::stringstream(optarg) >> opt.seed;
        break;
      }
      case 'T': {
        opt.transcriptsFile = optarg;
        break;
      }
      case 'p': {
        opt.priors = optarg;
        break;
      }
      default:
        break;
    }
  }
  if (long_read_flag) {
    opt.long_read = true;
  }
  if (matrix_to_files) {
    opt.matrix_to_files = true;
  }
  if (matrix_to_directories) {
    opt.matrix_to_files = true;  // This must be true if matrix_to_directories is true
    opt.matrix_to_directories = true;
  }
  if (plaintext_flag) {
    opt.plaintext = true;
  }
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
  if (opt.files.size() > 0) {
    opt.tccFile = opt.files[0];
  }
}

void ListSingleCellTechnologies() { bus_parser::ListSingleCellTechnologies(); }

void ParseOptionsBus(int argc, char **argv, ProgramOptions &opt) {
  bus_parser::ParseOptionsBus(argc, argv, opt);
}

bool ParseTechnology(const std::string &techstr, BUSOptions &busopt,
                     std::vector<std::string> &errorList) {
  return bus_parser::ParseTechnology(techstr, busopt, errorList);
}

void ParseOptionsH5Dump(int argc, char **argv, ProgramOptions &opt) {
  int peek_flag = 0;
  const char *opt_string = "o:";
  static struct option long_options[] = {// long args
                                         {"peek", no_argument, &peek_flag, 1},
                                         // short args
                                         {"output-dir", required_argument, 0, 'o'},
                                         {0, 0, 0, 0}};

  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
      case 0:
        break;
      case 'o': {
        opt.output = optarg;
        break;
      }
      default:
        break;
    }
  }

  // there should only be one thing here...
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }

  if (peek_flag) {
    opt.peek = true;
  }
}

bool CheckOptionsBus(ProgramOptions &opt) {
  // Initialize BUS options (to ensure no variables are uninitialized for good measure); BUS options
  // shouldn't be intialized before this
  opt.busOptions.nfiles = 1;
  opt.busOptions.keep_fastq_comments = false;
  opt.busOptions.paired = false;
  opt.busOptions.long_read = false;
  opt.busOptions.unmapped = false;
  opt.busOptions.error_rate = 0.0;
  opt.busOptions.threshold = 0.8;
  opt.busOptions.aa = false;

  bool ret = true;

  std::cerr << std::endl;

// check bam options
#ifdef NO_HTSLIB
  if (opt.bam || opt.genomebam || opt.pseudobam) {
    std::cerr << ERROR_STR << " in order to use BAM, must compile with BAM option enabled" << std::endl;
    ret = false;
  }
#endif

  // check index
  if (opt.index.empty()) {
    std::cerr << ERROR_STR << " kallisto index file missing" << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " kallisto index file not found " << opt.index << std::endl;
      ret = false;
    }
  }

  if (opt.max_num_reads < 0) {
    std::cerr << ERROR_STR << " --numReads must be a positive number" << std::endl;
    ret = false;
  }

  if (opt.long_read && !(0 < opt.threshold && opt.threshold < 1)) {
    std::cerr << "Threshold not in (0,1). Setting default threshold for unmapped kmers to 0.8"
              << std::endl;
    opt.threshold = 0.8;
  }

  if (opt.do_union && (opt.long_read || opt.aa)) {
    std::cerr << "--union is not compatible with this mode" << std::endl;
    ret = false;
  }
  if (opt.no_jump && (opt.long_read || opt.aa)) {
    std::cerr << "--no-jump is not compatible with this mode" << std::endl;
    ret = false;
  }

  if (opt.long_read) {  // opt.error_rate <= 0) {
    // hiding for release, not used for this version
    // std::cerr << "No sequencing error-rate: invalid error-rate; must be greater than zero" <<
    // std::endl; ret = false;
  }

  // check files
  bool read_stdin = opt.files.size() == 1 && !opt.batch_mode && !opt.bam &&
                    opt.files[0] == "-";  // - means read from stdin
  if (opt.files.size() == 0 && !opt.batch_mode) {
    std::cerr << ERROR_STR << " Missing read files" << std::endl;
    ret = false;
  } else if (!read_stdin) {
    struct stat stFileInfo;
    for (auto &fn : opt.files) {
      auto intStat = stat(fn.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << ERROR_STR << " file not found " << fn << std::endl;
        ret = false;
      }
    }
  }

  if (opt.input_interleaved_nfiles != 0) {
    if (opt.files.size() > 1) {
      std::cerr << ERROR_STR << " interleaved input cannot consist of more than one input" << std::endl;
      ret = false;
    }
    if (opt.bam) {
      std::cerr << ERROR_STR << " interleaved input is not compatible with the bam option" << std::endl;
      ret = false;
    }
    if (opt.batch_mode) {
      std::cerr << ERROR_STR << " interleaved input cannot be specified with a batch file" << std::endl;
      ret = false;
    }
  }

  if (opt.output.empty()) {
    std::cerr << "Error: need to specify output directory " << opt.output << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
        ret = false;
      }
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  }

  if (opt.threads <= 0) {
    std::cerr << "Error: invalid number of threads " << opt.threads << std::endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      std::cerr << "Warning: you asked for " << opt.threads << ", but only " << n
           << " cores on the machine" << std::endl;
      opt.threads = n;
    }
  }

  ProgramOptions::StrandType strand = ProgramOptions::StrandType::None;

  bool do_bulk_mode = (opt.technology == "BULK");
  if (do_bulk_mode) {
    opt.technology = "";
  }
  if (opt.technology.empty()) {  // kallisto pseudo
    // check for read files
    if (!opt.batch_mode) {
      if (ret && !do_bulk_mode) {
        std::cerr << "Error: the technology must be specified via -x, use \"bulk\" for regular RNA-seq "
                "reads"
             << std::endl;
        ret = false;
      }
      opt.batch_mode = true;
      if (ret && !opt.single_end && opt.files.size() % 2 != 0 &&
          opt.input_interleaved_nfiles == 0) {
        std::cerr << "Error: paired-end mode requires an even number of input files" << std::endl;
        ret = false;
      } else if (ret) {
        int i = 0;
        int sample_id = 0;
        while (i < opt.files.size()) {
          opt.batch_ids.push_back("batch" + std::to_string(sample_id));
          std::string f1, f2;
          struct stat stFileInfo;
          f1 = opt.files[i];
          if (opt.single_end) {
            opt.batch_files.push_back({f1});
            auto intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0 && !read_stdin) {
              std::cerr << ERROR_STR << " file not found " << f1 << std::endl;
              ret = false;
            }
          } else if (opt.input_interleaved_nfiles != 0) {
            f2 = "interleaved";
            opt.batch_files.push_back({f1, f2});
            opt.input_interleaved_nfiles = 2;
            auto intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0 && !read_stdin) {
              std::cerr << ERROR_STR << " file not found " << f1 << std::endl;
              ret = false;
            }
          } else {
            i++;
            f2 = opt.files[i];
            opt.batch_files.push_back({f1, f2});
            auto intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              std::cerr << ERROR_STR << " file not found " << f1 << std::endl;
              ret = false;
            }
            intstat = stat(f2.c_str(), &stFileInfo);
            if (intstat != 0) {
              std::cerr << ERROR_STR << " file not found " << f2 << std::endl;
              ret = false;
            }
          }
          i++;
          sample_id++;
        }
      }
    } else if (ret) {
      std::cerr << "[bus] will try running read files supplied in batch file" << std::endl;
      if (!opt.single_end) {
        std::cerr
            << "[bus] --paired ignored; single/paired-end is inferred from number of files supplied"
            << std::endl;
      }
      if (opt.files.size() != 0) {
        std::cerr << ERROR_STR << " cannot specify batch mode and supply read files" << std::endl;
        ret = false;
      } else {
        // check for batch files
        struct stat stFileInfo;
        auto intstat = stat(opt.batch_file_name.c_str(), &stFileInfo);
        if (intstat != 0) {
          std::cerr << ERROR_STR << " file not found " << opt.batch_file_name << std::endl;
          ret = false;
        }
        // open the file, parse and fill the batch_files values
        std::ifstream bfile(opt.batch_file_name);
        std::string line;
        std::string id, f1, f2;
        bool read_first_batch_file_line = false;
        while (std::getline(bfile, line)) {
          if (line.size() == 0) {
            continue;
          }
          std::stringstream ss(line);
          ss >> id;
          if (id[0] == '#') {
            continue;
          }
          opt.batch_ids.push_back(id);
          ss >> f1 >> f2;
          if (!read_first_batch_file_line) {
            if (f2.empty()) {
              opt.single_end = true;
            } else {
              opt.single_end = false;
            }
            read_first_batch_file_line = true;
          }
          if (opt.single_end) {
            opt.batch_files.push_back({f1});
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              std::cerr << ERROR_STR << " file not found " << f1 << std::endl;
              ret = false;
            }
            if (!f2.empty()) {
              std::cerr << ERROR_STR << " batch file malformatted" << std::endl;
              ret = false;
              break;
            }
          } else {
            opt.batch_files.push_back({f1, f2});
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              std::cerr << ERROR_STR << " file not found " << f1 << std::endl;
              ret = false;
            }
            if (f2.empty()) {
              std::cerr << ERROR_STR << " batch file malformatted" << std::endl;
              ret = false;
              break;
            }
            intstat = stat(f2.c_str(), &stFileInfo);
            if (intstat != 0) {
              std::cerr << ERROR_STR << " file not found " << f2 << std::endl;
              ret = false;
            }
          }
          f1.clear();
          f2.clear();
        }
      }
    }
    if (opt.pseudobam) {
      std::cerr << "Error: Pseudobam not supported yet in this mode" << std::endl;
      ret = false;
    }
    if (opt.bam) {
      std::cerr << "Error: --bam not supported in this mode" << std::endl;
      ret = false;
    }
    if (!opt.tagsequence.empty()) {
      std::cerr << "Error: --tag not supported in this mode" << std::endl;
      ret = false;
    }
    if (opt.strand_specific && opt.strand == ProgramOptions::StrandType::None) {
      opt.strand_specific = false;
    }
    if (ret) {  // Set up BUS options
      auto &busopt = opt.busOptions;
      busopt.seq.push_back(BUSOptionSubstr(0, 0, 0));
      busopt.umi.resize(0);
      busopt.bc.resize(0);
      busopt.aa = opt.aa;
      busopt.keep_fastq_comments = false;
      if (opt.single_end) {
        busopt.nfiles = 1;
        busopt.paired = false;
      } else {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
        busopt.paired = true;
      }
      if (opt.long_read) {
        busopt.long_read = true;
        busopt.paired = false;
        busopt.error_rate = opt.error_rate;
      }
      busopt.umi.push_back(BUSOptionSubstr(-1, -1, -1));
    }
    return ret;
  } else {  // User supplied -x (technology option)
    int nfiles_per_batch = 0;
    if (opt.batch_mode) {  // using -x with batch
      std::cerr << "[bus] will try running read files supplied in batch file" << std::endl;
      if (opt.files.size() != 0) {
        std::cerr << ERROR_STR << " cannot specify batch mode and supply read files" << std::endl;
        ret = false;
      }
      struct stat stFileInfo;
      auto intstat = stat(opt.batch_file_name.c_str(), &stFileInfo);
      if (intstat != 0) {
        std::cerr << ERROR_STR << " file not found " << opt.batch_file_name << std::endl;
        ret = false;
      }
      // open the file, parse and fill the batch_files values
      std::ifstream bfile(opt.batch_file_name);
      std::string line;
      std::string id;
      std::vector<std::string> f;
      bool read_first_batch_file_line = false;
      while (std::getline(bfile, line)) {
        if (line.size() == 0) {
          continue;
        }
        std::stringstream ss(line);
        ss >> id;
        if (id[0] == '#') {
          continue;
        }
        opt.batch_ids.push_back(id);
        std::string s;
        while (ss >> s) {
          f.push_back(s);
          opt.files.push_back(s);
          intstat = stat(s.c_str(), &stFileInfo);
          if (intstat != 0) {
            std::cerr << ERROR_STR << " file not found " << s << std::endl;
            ret = false;
          }
        }
        opt.batch_files.push_back(f);
        if (!read_first_batch_file_line) {
          read_first_batch_file_line = true;
          nfiles_per_batch = f.size();
        }
        if (nfiles_per_batch != f.size()) {
          std::cerr << ERROR_STR << " batch file malformatted" << std::endl;
          ret = false;
          break;
        }
        f.clear();
      }
    }
    if (!ret)
      return false;
    auto &busopt = opt.busOptions;
    busopt.aa = opt.aa;

    busopt.long_read = opt.long_read;
    busopt.threshold = opt.threshold;
    busopt.paired = false;
    busopt.keep_fastq_comments = false;
    if (opt.bam) {
      busopt.nfiles = 1;  // Note: only 10xV2 has been tested
    }
    {
      std::vector<std::string> errorList;
      bool ok = bus_parser::ApplyBusPresetTechnology(opt, strand, errorList);
      if (!ok) {
        for (const auto &e : errorList) std::cerr << e << std::endl;
        ret = false;
      }
    }
    if (opt.bam) {
      busopt.nfiles = 1;  // Make sure it's still one file for BAM
    }
    if (nfiles_per_batch != 0 && nfiles_per_batch != busopt.nfiles) {
      std::cerr << "Wrong number of files per batch for technology: " << opt.technology << std::endl;
    }
    if (opt.batch_mode) {
      opt.files.clear();  // Don't need these
    }
  }

  if (opt.tagsequence.empty() && (opt.technology == "SMARTSEQ3")) {
    opt.tagsequence = "ATTGCGCAATG";
    std::cerr << "[bus] Using " << opt.tagsequence << " as UMI tag sequence" << std::endl;
  }

  if (opt.strand_specific && opt.strand == ProgramOptions::StrandType::None) {
    opt.strand_specific = false;  // User specified --unstranded
  } else if (!opt.strand_specific && ret) {
    opt.strand = strand;
    if (opt.strand == ProgramOptions::StrandType::FR) {
      std::cerr << "[bus] Note: Strand option was not specified; setting it to --fr-stranded for "
              "specified technology"
           << std::endl;
      opt.strand_specific = true;
    } else if (opt.strand == ProgramOptions::StrandType::RF) {
      std::cerr << "[bus] Note: Strand option was not specified; setting it to --rf-stranded for "
              "specified technology"
           << std::endl;
      opt.strand_specific = true;
    } else {
      std::cerr << "[bus] Note: Strand option was not specified; setting it to --unstranded for "
              "specified technology"
           << std::endl;
    }
  }

  if (!opt.tagsequence.empty()) {
    opt.busOptions.umi[0].start += opt.tagsequence.length();
    if (opt.busOptions.umi[0].start >= opt.busOptions.umi[0].stop) {
      std::cerr << "Error: Tag sequence must be shorter than UMI sequence" << std::endl;
      ret = false;
    }
  }

  if (!opt.single_end && !opt.busOptions.paired && !opt.long_read) {
    std::cerr << "Error: Paired reads are not compatible with the specified technology" << std::endl;
    ret = false;
  }

  if (opt.strand_specific) {
    if (opt.busOptions.seq.size() != 1 &&
        !(opt.busOptions.seq.size() == 2 && opt.busOptions.paired)) {
      std::cerr << "Error: Strand-specific read processing is only supported for technologies with a "
              "single CDNA read file or paired-end reads"
           << std::endl;
      ret = false;
    }
  }

  if (opt.genomebam) {
    if (opt.busOptions.seq.size() != 1 &&
        !(opt.busOptions.seq.size() == 2 && opt.busOptions.paired)) {
      std::cerr << "Error: BAM output is currently only supported for technologies with a single CDNA "
              "read file or paired-end reads"
           << std::endl;
      ret = false;
    }
    if (!opt.gtfFile.empty()) {
      if (!checkFileExists(opt.gtfFile)) {
        std::cerr << "Error: GTF file " << opt.gtfFile << " does not exist" << std::endl;
        ret = false;
      }
    } else {
      std::cerr << "Error: need GTF file for genome alignment" << std::endl;
      ret = false;
    }
    if (!opt.chromFile.empty()) {
      if (!checkFileExists(opt.chromFile)) {
        std::cerr << "Error: Chromosome file not found: " << opt.chromFile << std::endl;
        ret = false;
      }
    }
  }

  if (ret && !opt.bam && !(opt.input_interleaved_nfiles != 0) &&
      opt.files.size() % opt.busOptions.nfiles != 0) {
    std::cerr << "Error: Number of files (" << opt.files.size()
         << ") does not match number of input files required by "
         << "technology " << opt.technology << " (" << opt.busOptions.nfiles << ")" << std::endl;
    ret = false;
  }
  if (opt.input_interleaved_nfiles != 0) {
    opt.input_interleaved_nfiles = opt.busOptions.nfiles;
    for (int i = 1; i < opt.busOptions.nfiles; i++) {
      opt.files.push_back("interleaved");
    }
  }

  if (opt.bam && opt.num) {
    std::cerr << "Warning: --bam option was used, so --num option will be ignored" << std::endl;
  }

  return ret;
}

bool CheckOptionsIndex(ProgramOptions &opt) {
  bool ret = true;

  if (opt.threads <= 0) {
    std::cerr << "Error: invalid number of threads " << opt.threads << std::endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      std::cerr << "Warning: you asked for " << opt.threads << ", but only " << n
           << " cores on the machine" << std::endl;
      opt.threads = n;
    }
  }

  if (opt.k <= 1 || opt.k >= MAX_KMER_SIZE) {
    std::cerr << "Error: invalid k-mer length " << opt.k << ", minimum is 3 and maximum is "
         << (MAX_KMER_SIZE - 1) << std::endl;
    ret = false;
  }

  if (opt.k % 2 == 0) {
    std::cerr << "Error: k needs to be an odd number" << std::endl;
    ret = false;
  }

  if (opt.transfasta.empty()) {
    std::cerr << "Error: no FASTA files specified" << std::endl;
    ret = false;
  } else {
    for (auto &fasta : opt.transfasta) {
      // we want to generate the index, check k, index and transfasta
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << "Error: FASTA file not found " << fasta << std::endl;
        ret = false;
      }
    }
  }

  if (!opt.d_list.empty()) {
    for (auto &fasta : opt.d_list) {
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << "Error: D-list FASTA file not found \"" << fasta << "\"" << std::endl;
        ret = false;
      }
    }
  }

  if (opt.index.empty()) {
    std::cerr << "Error: need to specify kallisto index name" << std::endl;
    ret = false;
  }

  if (opt.g != 0) {
    if (opt.g <= 2 || opt.g > opt.k - 2) {
      std::cerr << "Error: invalid minimizer size " << opt.g << ", minimum is 3 and maximum is k - 2"
           << std::endl;
      ret = false;
    }
  }
  if (opt.max_ec_size < -1) {
    std::cerr << "Error: invalid max ec size " << opt.max_ec_size << std::endl;
    ret = false;
  }

  return ret;
}

bool CheckOptionsEM(ProgramOptions &opt, bool emonly = false) {
  bool ret = true;

#ifdef NO_HTSLIB
  if (opt.bam || opt.genomebam || opt.pseudobam) {
    std::cerr << ERROR_STR << " in order to use BAM, must compile with BAM option enabled" << std::endl;
    ret = false;
  }
#endif

  std::cerr << std::endl;
  // check for index
  if (!emonly) {
    if (opt.index.empty()) {
      std::cerr << ERROR_STR << " kallisto index file missing" << std::endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      auto intStat = stat(opt.index.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << ERROR_STR << " kallisto index file not found " << opt.index << std::endl;
        ret = false;
      }
    }
  }

  // check for read files
  if (!emonly) {
    if (opt.files.size() == 0) {
      std::cerr << ERROR_STR << " Missing read files" << std::endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      for (auto &fn : opt.files) {
        auto intStat = stat(fn.c_str(), &stFileInfo);
        if (intStat != 0) {
          std::cerr << ERROR_STR << " file not found " << fn << std::endl;
          ret = false;
        }
      }
    }

    /*
    if (opt.strand_specific && !opt.single_end && !opt.long_read) {
      std::cerr << "Error: strand-specific mode requires single-end mode" << std::endl;
      ret = false;
    }*/

    if (!opt.single_end && !opt.long_read) {
      if (opt.files.size() % 2 != 0) {
        std::cerr << "Error: paired-end mode requires an even number of input files" << std::endl
             << "       (use --single for processing single-end reads)" << std::endl;
        ret = false;
      }
    }
  }

  if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
    std::cerr << "Error: cannot supply mean/sd without supplying both -l and -s" << std::endl;
    ret = false;
  }

  if (opt.do_union && (opt.long_read || opt.aa)) {
    std::cerr << "--union is not compatible with this mode" << std::endl;
    ret = false;
  }
  if (opt.no_jump && (opt.long_read || opt.aa)) {
    std::cerr << "--no-jump is not compatible with this mode" << std::endl;
    ret = false;
  }

  if ((opt.single_end && !opt.long_read) && (opt.fld == 0.0 || opt.sd == 0.0)) {
    std::cerr << "Error: fragment length mean and sd must be supplied for single-end reads using -l and "
            "-s"
         << std::endl;
    ret = false;
  } else if (opt.fld == 0.0 && ret) {
    // In the future, if we have single-end data we should require this
    // argument
    std::cerr << "[quant] fragment length distribution will be estimated from the data" << std::endl;
  } else if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
    std::cerr << "[quant] fragment length distribution is truncated gaussian with mean = " << opt.fld
         << ", sd = " << opt.sd << std::endl;
  }

  if (!opt.single_end && (opt.fld > 0.0 && opt.sd > 0.0)) {
    std::cerr << "[~warn] you specified using a gaussian but have paired end data" << std::endl;
    std::cerr << "[~warn] we suggest omitting these parameters and let us estimate the distribution "
            "from data"
         << std::endl;
  }

  if (opt.fld < 0.0) {
    std::cerr << "Error: invalid value for mean fragment length " << opt.fld << std::endl;
    ret = false;
  }

  if (opt.sd < 0.0) {
    std::cerr << "Error: invalid value for fragment length standard deviation " << opt.sd << std::endl;
    ret = false;
  }

  if (opt.iterations <= 0) {
    std::cerr << "Error: invalid number of iterations " << opt.iterations << std::endl;
    ret = false;
  }

  if (opt.min_range <= 0) {
    std::cerr << "Error: invalid value for minimum range " << opt.min_range << std::endl;
    ret = false;
  }

  if (opt.genomebam) {
    if (!opt.gtfFile.empty()) {
      struct stat stFileInfo;
      auto intStat = stat(opt.gtfFile.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << "Error: GTF file " << opt.gtfFile << " does not exist" << std::endl;
        ret = false;
      }
    } else {
      std::cerr << "Error: need GTF file for genome alignment" << std::endl;
      ret = false;
    }

    if (!opt.chromFile.empty()) {
      struct stat stFileInfo;
      auto intStat = stat(opt.chromFile.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << "Error: Chromosome file not found: " << opt.chromFile << std::endl;
        ret = false;
      }
    }
  }

  if (opt.output.empty()) {
    std::cerr << "Error: need to specify output directory " << opt.output << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        std::cerr << "Error: file " << opt.output << " exists and is not a directory" << std::endl;
        ret = false;
      } else if (emonly) {
        // check for directory/counts.txt
        struct stat stCountInfo;
        auto intcountstat = stat((opt.output + "/counts.txt").c_str(), &stCountInfo);
        if (intcountstat != 0) {
          std::cerr << "Error: could not find file " << opt.output << "/counts.txt" << std::endl;
          ret = false;
        }

        // check for directory/index.saved
        struct stat stIndexInfo;
        auto intindexstat = stat((opt.output + "/index.saved").c_str(), &stIndexInfo);
        if (intindexstat != 0) {
          std::cerr << "Error: could not find index " << opt.output << "/index.saved" << std::endl;
          ret = false;
        }
        opt.index = (opt.output + "/index.saved");
      }
    } else {
      if (emonly) {
        std::cerr << "Error: output directory needs to exist, run quant first" << std::endl;
        ret = false;
      } else {
        // create directory
        if (my_mkdir(opt.output.c_str(), 0777) == -1) {
          std::cerr << "Error: could not create directory " << opt.output << std::endl;
          ret = false;
        }
      }
    }
  }

  if (opt.threads <= 0) {
    std::cerr << "Error: invalid number of threads " << opt.threads << std::endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      std::cerr << "Warning: you asked for " << opt.threads << ", but only " << n
           << " cores on the machine" << std::endl;
    }
    if (opt.threads > 1 && opt.pseudobam) {
      // std::cerr << "Error: pseudobam is not compatible with running on many threads."<< std::endl;
      // ret = false;
    }
  }

  if (opt.bootstrap < 0) {
    std::cerr << "Error: number of bootstrap samples must be a non-negative integer." << std::endl;
    ret = false;
  }

#ifndef USE_HDF5
  if (opt.bootstrap > 0 && !opt.plaintext) {
    std::cerr << "Warning: kallisto was not compiled with HDF5 support so no bootstrapping" << std::endl
         << "will be performed. Run quant with --plaintext option or recompile with" << std::endl
         << "HDF5 support to obtain bootstrap estimates." << std::endl;
    opt.bootstrap = 0;
  }
#endif
  return ret;
}

bool CheckOptionsTCCQuant(ProgramOptions &opt) {
  bool ret = true;

  std::cerr << std::endl;
  // check for index
  if (opt.index.empty() && opt.transcriptsFile.empty()) {
    std::cerr << ERROR_STR << " either a kallisto index file or a transcripts file need to be supplied"
         << std::endl;
    ret = false;
  } else if (!opt.index.empty() && !opt.transcriptsFile.empty()) {
    std::cerr << ERROR_STR << " cannot supply both a kallisto index file and a transcripts file" << std::endl;
    ret = false;
  } else if (!opt.index.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " kallisto index file not found " << opt.index << std::endl;
      ret = false;
    }
  } else if (!opt.transcriptsFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.transcriptsFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " transcripts file not found " << opt.transcriptsFile << std::endl;
      ret = false;
    }
  }

  if (opt.tccFile.empty()) {
    std::cerr << ERROR_STR << " transcript-compatibility counts file missing" << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.tccFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " transcript-compatibility counts file not found " << opt.tccFile
           << std::endl;
      ret = false;
    }
  }

  if (opt.ecFile.empty() && !opt.transcriptsFile.empty()) {
    std::cerr << ERROR_STR << " equivalence class file must be supplied if transcripts file is supplied "
         << opt.tccFile << std::endl;
    ret = false;
  }

  if (!opt.ecFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.ecFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " equivalence class file not found " << opt.ecFile << std::endl;
      ret = false;
    }
  }

  if (!opt.fldFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.fldFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " fragment length distribution file not found " << opt.fldFile << std::endl;
      ret = false;
    }
  }

  if (!opt.genemap.empty() && !opt.gtfFile.empty()) {
    std::cerr << ERROR_STR << " Cannot supply both --genemap and --gtf" << std::endl;
    ret = false;
  }

  if (!opt.gtfFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.gtfFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " GTF file not found " << opt.gtfFile << std::endl;
      ret = false;
    }
  }

  if (!opt.genemap.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.genemap.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << ERROR_STR << " file for mapping transcripts to genes not found " << opt.genemap
           << std::endl;
      ret = false;
    }
  }

  if ((opt.fld != 0.0 || opt.sd != 0.0) && !opt.fldFile.empty()) {
    std::cerr << ERROR_STR
         << " cannot supply mean or sd while also supplying a fragment length distribution file"
         << std::endl;
    ret = false;
  }

  if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
    std::cerr << ERROR_STR << " cannot supply mean/sd without supplying both -l and -s" << std::endl;
    ret = false;
  }

  if (opt.index.empty() && (!opt.fldFile.empty() || opt.fld != 0.0 || opt.sd != 0.0)) {
    std::cerr << ERROR_STR
         << " cannot supply fragment length information unless a kallisto index is provided"
         << std::endl;
    ret = false;
  }

  if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
    std::cerr << "[tcc] fragment length distribution is truncated gaussian with mean = " << opt.fld
         << ", sd = " << opt.sd << std::endl;
  }

  if (opt.fld < 0.0) {
    std::cerr << ERROR_STR << " invalid value for mean fragment length " << opt.fld << std::endl;
    ret = false;
  }

  if (opt.sd < 0.0) {
    std::cerr << ERROR_STR << " invalid value for fragment length standard deviation " << opt.sd << std::endl;
    ret = false;
  }

  if (opt.index.empty() && (opt.matrix_to_files || opt.matrix_to_directories)) {
    std::cerr << ERROR_STR << " cannot get abundance tsv files unless a kallisto index is provided"
         << std::endl;
    ret = false;
  }

  if (opt.output.empty()) {
    std::cerr << ERROR_STR << " need to specify output directory " << opt.output << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        std::cerr << ERROR_STR << " file " << opt.output << " exists and is not a directory" << std::endl;
        ret = false;
      }
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        std::cerr << ERROR_STR << " could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  }

  if (opt.threads <= 0) {
    std::cerr << ERROR_STR << " invalid number of threads " << opt.threads << std::endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      std::cerr << "Warning: you asked for " << opt.threads << ", but only " << n
           << " cores on the machine" << std::endl;
    }
  }

  if (opt.bootstrap < 0) {
    std::cerr << "Error: number of bootstrap samples must be a non-negative integer." << std::endl;
    ret = false;
  }

  return ret;
}

bool CheckOptionsInspect(ProgramOptions &opt) {
  bool ret = true;

  if (opt.threads <= 0) {
    std::cerr << "Error: invalid number of threads " << opt.threads << std::endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      std::cerr << "Warning: you asked for " << opt.threads << ", but only " << n
           << " cores on the machine" << std::endl;
      opt.threads = n;
    }
  }

  // check for index
  if (opt.index.empty()) {
    std::cerr << "Error: kallisto index file missing" << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << "Error: kallisto index file not found " << opt.index << std::endl;
      ret = false;
    }
  }

  if (!opt.bedFile.empty() || !opt.gtfFile.empty()) {
    opt.pseudobam = true;
    opt.genomebam = true;
  }

  if (opt.genomebam) {
    struct stat stFileInfo;
    auto intStat = stat(opt.gtfFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      std::cerr << "Error: GTF file not found " << opt.gtfFile << std::endl;
      ret = false;
    }

    if (!opt.chromFile.empty()) {
      struct stat stFileInfo;
      auto intStat = stat(opt.chromFile.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << "Error: Chromosome file not found " << opt.chromFile << std::endl;
        ret = false;
      }
    }

    if (opt.bedFile.empty()) {
      opt.bedFile = opt.index + ".bed";
    }
  }

  return ret;
}

bool CheckOptionsH5Dump(ProgramOptions &opt) {
  bool ret = true;
  if (!opt.peek) {
    if (opt.output.size() == 0) {
      std::cerr << "Error: You must specify an output directory." << std::endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      auto intStat = stat(opt.output.c_str(), &stFileInfo);
      if (intStat == 0) {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode)) {
          std::cerr << "Error: tried to open " << opt.output << " but another file already exists there"
               << std::endl;
          ret = false;
        }
      } else if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        std::cerr << "Error: could not create directory " << opt.output << std::endl;
        ret = false;
      }
    }
  } else {
    if (opt.output.size() > 0) {
      std::cerr << "Error: Cannot specify output directory and '--peek'. Please specify only one."
           << std::endl;
      ret = false;
    }
  }

  if (opt.files.size() == 0) {
    std::cerr << "Error: Missing H5 files" << std::endl;
    ret = false;
  } else if (opt.files.size() > 1) {
    std::cerr << "Error: Please specify only one H5 file" << std::endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    for (auto &fn : opt.files) {
      auto intStat = stat(fn.c_str(), &stFileInfo);
      if (intStat != 0) {
        std::cerr << "Error: H5 file not found " << fn << std::endl;
        ret = false;
      }
    }
  }

  return ret;
}

void PrintCite() {
  std::cout << "When using this program in your research, please cite" << std::endl
       << std::endl
       << "  Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L." << std::endl
       << "  Near-optimal probabilistic RNA-seq quantification, " << std::endl
       << "  Nature Biotechnology 34, 525-527(2016), doi:10.1038/nbt.3519" << std::endl
       << std::endl;
}

void PrintVersion() {
  std::cout << "kallisto, version " << KALLISTO_VERSION;
#if KALLISTO_CPP_VERSION < 201703L
  std::cout << ".";  // Just an internal use hidden "marker" so we know what C++ version kallisto
                     // was compiled in
#endif
  std::cout << std::endl;
}

void usage() {
  std::cout << "kallisto " << KALLISTO_VERSION << std::endl
       << std::endl
       << "Usage: kallisto <CMD> [arguments] .." << std::endl
       << std::endl
       << "Where <CMD> can be one of:" << std::endl
       << std::endl
       << "    index         Builds a kallisto index " << std::endl
       << "    quant         Runs the quantification algorithm " << std::endl
       << "    quant-tcc     Runs quantification on transcript-compatibility counts" << std::endl
       << "    bus           Generate BUS files for single-cell data " << std::endl
       << "    h5dump        Converts HDF5-formatted results to plaintext" << std::endl
       << "    inspect       Inspects and gives information about an index" << std::endl
       << "    version       Prints version information" << std::endl
       << "    cite          Prints citation information" << std::endl
       << std::endl
       << "Running kallisto <CMD> without arguments prints usage information for <CMD>" << std::endl
       << std::endl;
}

void usageBus() {
  std::cout << "kallisto " << KALLISTO_VERSION << std::endl
       << "Generates BUS files for single-cell sequencing" << std::endl
       << std::endl
       << "Usage: kallisto bus [arguments] FASTQ-files" << std::endl
       << std::endl
       << "Required arguments:" << std::endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << std::endl
       << "                              pseudoalignment" << std::endl
       << "-o, --output-dir=STRING       Directory to write output to" << std::endl
       << std::endl
       << "Optional arguments:" << std::endl
       << "-x, --technology=STRING       Single-cell technology used " << std::endl
       << "-l, --list                    List all single-cell technologies supported" << std::endl
       << "-B, --batch=FILE              Process files listed in FILE" << std::endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << std::endl
       << "-b, --bam                     Input file is a BAM file" << std::endl
       << "-n, --num                     Output number of read in flag column (incompatible with "
          "--bam)"
       << std::endl
       << "-N, --numReads                Maximum number of reads to process from supplied input"
       << std::endl
       << "-T, --tag=STRING              5′ tag sequence to identify UMI reads for certain "
          "technologies"
       << std::endl
       << "    --fr-stranded             Strand specific reads for UMI-tagged reads, first read "
          "forward"
       << std::endl
       << "    --rf-stranded             Strand specific reads for UMI-tagged reads, first read "
          "reverse"
       << std::endl
       << "    --unstranded              Treat all read as non-strand-specific" << std::endl
       << "    --paired                  Treat reads as paired" << std::endl
       << "    --long                    Treat reads as long"
       << std::endl
       //<< "    --error_rate              Estimated error rate of long reads (required for --long)"
       //<< std::endl
       << "    --threshold               Threshold for rate of unmapped kmers per read"
       << std::endl
       //<< "    --unmapped              Computed ratio of unmapped kmers for first 1M
       //reads and output unmapped reads" << std::endl
       << "    --aa                      Align to index generated from a FASTA-file containing "
          "amino acid sequences"
       << std::endl
       << "    --inleaved                Specifies that input is an interleaved FASTQ file" << std::endl
       << "    --batch-barcodes          Records both batch and extracted barcode in BUS file"
       << std::endl
       #ifndef NO_HTSLIB
       << "    --genomebam               Project pseudoalignments to genome sorted BAM file" << std::endl
       << "-g, --gtf                     GTF file for transcriptome information" << std::endl
       << "                              (required for --genomebam)" << std::endl
       << "-c, --chromosomes             Tab separated file with chromosome names and lengths" << std::endl
       << "                              (optional for --genomebam, but recommended)" << std::endl
       #endif // NO_HTSLIB
       << "    --verbose                 Print out progress information every 1M proccessed reads"
       << std::endl;
}

void usageIndex() {
  std::cout << "kallisto " << KALLISTO_VERSION << std::endl
       << "Builds a kallisto index" << std::endl
       << std::endl
       << "Usage: kallisto index [arguments] FASTA-files" << std::endl
       << std::endl
       << "Required argument:" << std::endl
       << "-i, --index=STRING          Filename for the kallisto index to be constructed " << std::endl
       << std::endl
       << "Optional argument:" << std::endl
       << "-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: "
       << (MAX_KMER_SIZE - 1) << ")" << std::endl
       << "-t, --threads=INT           Number of threads to use (default: 1)" << std::endl
       << "-d, --d-list=STRING         Path to a FASTA-file containing sequences to mask from "
          "quantification"
       << std::endl
       << "    --make-unique           Replace repeated target names with unique names" << std::endl
       << "    --aa                    Generate index from a FASTA-file containing amino acid "
          "sequences"
       << std::endl
       << "    --distinguish           Generate index where sequences are distinguished by the "
          "sequence name"
       << std::endl
       << "-T, --tmp=STRING            Temporary directory (default: tmp)" << std::endl
       << "-m, --min-size=INT          Length of minimizers (default: automatically chosen)" << std::endl
       << "-e, --ec-max-size=INT       Maximum number of targets in an equivalence class (default: "
          "no maximum)"
       << std::endl
       << std::endl;
}

void usageh5dump() {
  std::cout << "kallisto " << KALLISTO_VERSION << std::endl
       << "Converts HDF5-formatted results to plaintext" << std::endl
       << std::endl
       << "Usage:  kallisto h5dump [arguments] abundance.h5" << std::endl
       << std::endl
       << "Required argument:" << std::endl
       << "-o, --output-dir=STRING       Directory to write output to" << std::endl
       << std::endl;
}

void usageInspect() {
  std::cout << "kallisto " << KALLISTO_VERSION << std::endl
       << std::endl
       << "Usage: kallisto inspect INDEX-file" << std::endl
       << std::endl
       << "Optional arguments:" << std::endl
       << "-t                      Number of threads" << std::endl
       << std::endl;
}

void usageEM(bool valid_input = true) {
  if (valid_input) {
    std::cout << "kallisto " << KALLISTO_VERSION << std::endl
         << "Computes equivalence classes for reads and quantifies abundances" << std::endl
         << std::endl;
  }
  //      "----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|"
  std::cout << "Usage: kallisto quant [arguments] FASTQ-files" << std::endl
       << std::endl
       << "Required arguments:" << std::endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << std::endl
       << "                              quantification" << std::endl
       << "-o, --output-dir=STRING       Directory to write output to" << std::endl
       << std::endl
       << "Optional arguments:" << std::endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << std::endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << std::endl
       << "    --plaintext               Output plaintext instead of HDF5" << std::endl
       << "    --single                  Quantify single-end reads" << std::endl
       << "    --single-overhang         Include reads where unobserved rest of fragment is" << std::endl
       << "                              predicted to lie outside a transcript" << std::endl
       << "    --fr-stranded             Strand specific reads, first read forward" << std::endl
       << "    --rf-stranded             Strand specific reads, first read reverse" << std::endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << std::endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << std::endl
       << "                              (default: -l, -s values are estimated from paired" << std::endl
       << "                               end data, but are required when using --single)" << std::endl
       << "-p, --priors                  Priors for the EM algorithm, either as raw counts or as"
       << std::endl
       << "                              probabilities. Pseudocounts are added to raw reads to"
       << std::endl
       << "                              prevent zero valued priors. Supplied in the same order"
       << std::endl
       << "                              as the transcripts in the transcriptome" << std::endl
       #ifndef NO_HTSLIB
       << "    --pseudobam               Save pseudoalignments to transcriptome to BAM file" << std::endl
       << "    --genomebam               Project pseudoalignments to genome sorted BAM file" << std::endl
       << "-g, --gtf                     GTF file for transcriptome information" << std::endl
       << "                              (required for --genomebam)" << std::endl
       << "-c, --chromosomes             Tab separated file with chromosome names and lengths" << std::endl
       << "                              (optional for --genomebam, but recommended)" << std::endl
       #endif // NO_HTSLIB
       << "-t, --threads=INT             Number of threads to use (default: 1)" << std::endl
       << "    --verbose                 Print out progress information every 1M proccessed reads"
       << std::endl;
}

void usageEMOnly() {
  std::cout << "kallisto " << KALLISTO_VERSION << std::endl
       << "Computes equivalence classes for reads and quantifies abundance" << std::endl
       << std::endl
       << "Usage: kallisto quant-only [arguments]" << std::endl
       << std::endl
       << "Required argument:" << std::endl
       << "-o, --output-dir=STRING       Directory to write output to" << std::endl
       << std::endl
       << "Optional arguments:" << std::endl
       << "-l, --fragment-length=DOUBLE  Estimated fragment length (default: value is estimated "
          "from the input data)"
       << std::endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << std::endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << std::endl
       << "    --plaintext               Output plaintext instead of HDF5" << std::endl
       << std::endl;
}

void usageTCCQuant(bool valid_input = true) {
  if (valid_input) {
    std::cout << "kallisto " << KALLISTO_VERSION << std::endl
         << "Quantifies abundance from pre-computed transcript-compatibility counts" << std::endl
         << std::endl;
  }

  std::cout << "Usage: kallisto quant-tcc [arguments] transcript-compatibility-counts-file" << std::endl
       << std::endl
       << "Required arguments:" << std::endl
       << "-o, --output-dir=STRING       Directory to write output to" << std::endl
       << std::endl
       << "Optional arguments:" << std::endl
       << "-i, --index=STRING            Filename for the kallisto index to be used" << std::endl
       << "                              (required if file with names of transcripts not supplied)"
       << std::endl
       << "-T, --txnames=STRING          File with names of transcripts" << std::endl
       << "                              (required if index file not supplied)" << std::endl
       << "-e, --ec-file=FILE            File containing equivalence classes" << std::endl
       << "                              (default: equivalence classes are taken from the index)"
       << std::endl
       << "-f, --fragment-file=FILE      File containing fragment length distribution" << std::endl
       << "                              (default: effective length normalization is not performed)"
       << std::endl
       << "--long                        Use version of EM for long reads " << std::endl
       << "-P, --platform.               [PacBio or ONT] used for sequencing " << std::endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << std::endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << std::endl
       << "                              (note: -l, -s values only should be supplied when" << std::endl
       << "                               effective length normalization needs to be performed"
       << std::endl
       << "                               but --fragment-file is not specified)" << std::endl
       << "-p, --priors                  Priors for the EM algorithm, either as raw counts or as"
       << std::endl
       << "                              probabilities. Pseudocounts are added to raw reads to"
       << std::endl
       << "                              prevent zero valued priors. Supplied in the same order"
       << std::endl
       << "                              as the transcripts in the transcriptome" << std::endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << std::endl
       << "-g, --genemap                 File for mapping transcripts to genes" << std::endl
       << "                              (required for obtaining gene-level abundances)" << std::endl
       << "-G, --gtf=FILE                GTF file for transcriptome information" << std::endl
       << "                              (can be used instead of --genemap for obtaining "
          "gene-level abundances)"
       << std::endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << std::endl
       << "    --matrix-to-files         Reorganize matrix output into abundance tsv files" << std::endl
       << "    --matrix-to-directories   Reorganize matrix output into abundance tsv files across "
          "multiple directories"
       << std::endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << std::endl
       << "    --plaintext               Output plaintext only, not HDF5" << std::endl;
}

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
  struct tm *timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  std::string ret(asctime(timeinfo));

  // chomp off the newline
  return ret.substr(0, ret.size() - 1);
}

int main(int argc, char *argv[]) {
  std::cout.sync_with_stdio(false);
  setvbuf(stdout, NULL, _IOFBF, 1048576);

  if (argc < 2) {
    usage();
    exit(1);
  } else {
    auto start_time(get_local_time());
    ProgramOptions opt;
   std::string cmd(argv[1]);
    if (cmd == "version") {
      PrintVersion();
    } else if (cmd == "cite") {
      PrintCite();
    } else if (cmd == "index") {
      std::cerr << std::endl;
      if (argc == 2) {
        usageIndex();
        return 0;
      }
      ParseOptionsIndex(argc - 1, argv + 1, opt);
      if (!CheckOptionsIndex(opt)) {
        usageIndex();
        exit(1);
      } else {
        // create an index
        opt.tmp_dir = (opt.tmp_dir.empty() ? "tmp" : opt.tmp_dir);
        Kmer::set_k(opt.k);
        KmerIndex index(opt);

        std::ofstream out;
        out.open(opt.index, std::ios::out | std::ios::binary);
        if (opt.distinguish) {
          index.BuildDistinguishingGraph(opt, out);
        } else {
          index.BuildTranscripts(opt, out);
        }
        index.write(out, opt.threads);
        rmdir(opt.tmp_dir.c_str());  // Remove temp directory if non-empty
      }
      std::cerr << std::endl;
    } else if (cmd == "inspect") {
      if (argc == 2) {
        usageInspect();
        return 0;
      }
      ParseOptionsInspect(argc - 1, argv + 1, opt);
      if (!CheckOptionsInspect(opt)) {
        usageInspect();
        exit(1);
      } else {
        KmerIndex index(opt);
        index.load(opt);
        InspectIndex(index, opt);
      }
    } else if (cmd == "bus") {
      if (argc == 2) {
        usageBus();
        return 0;
      }

      ParseOptionsBus(argc - 1, argv + 1, opt);
      int64_t num_processed = 0;
      if (!CheckOptionsBus(opt)) {
        usageBus();
        exit(1);
      }
      int64_t num_pseudoaligned = 0;
      int64_t num_unique = 0;
      uint32_t bclen = 0;
      uint32_t umilen = 0;
      std::string prefix = opt.output + "/matrix";
      std::string ecfilename = prefix + ".ec";
      std::string cellnamesfilename = prefix + ".cells";
      std::string genelistname = prefix + ".genelist.txt";
      std::string busbarcodelistname = prefix + ".sample.barcodes";
      bool batch_mode = opt.batch_mode;

      if (!batch_mode) {
        opt.bus_mode = true;  // bus_mode = !batch_mode (however, if either is true, means we're
                              // writing BUS file)
        opt.single_end = false;
      }

      KmerIndex index(opt);
      index.load(opt);

      if (opt.long_read) {
        double error_rate_threshold_tmp = ((1.0 / opt.error_rate - 2 * index.k) * opt.error_rate);
        // std::cerr << "Suggested threshold for novel reads to " << error_rate_threshold_tmp <<
        // std::endl;
        if (1 > opt.threshold && opt.threshold > 0) {
          // std::cerr << "Using supplied threshold " << opt.threshold << std::endl;
        } else if (0 < error_rate_threshold_tmp && error_rate_threshold_tmp < 1) {
          opt.threshold = error_rate_threshold_tmp;
          std::cerr << "Using computed threshold " << opt.threshold << std::endl;
        } else {
          opt.threshold = .8;
          std::cerr << "Supplied and computed threshold are invalid, using default value of "
                    << opt.threshold << std::endl;
        }
        /***
        int suggested_k = int((1.0/opt.error_rate)/2.0) - 1;
        if (suggested_k % 2 == 0) {
           std::cerr << "Suggested kmer length for error rate is: " << suggested_k-1 << std::endl;
        } else {
           std::cerr << "Suggested kmer length for error rate is: " << suggested_k << std::endl;
        }
        ***/
      }

      bool guessChromosomes = true;
      Transcriptome model;  // empty
      if (opt.genomebam) {
        if (!opt.chromFile.empty()) {
          guessChromosomes = false;
          model.loadChromosomes(opt.chromFile);
        } else {
          guessChromosomes = true;
        }
      }
      if (!opt.gtfFile.empty()) {
        model.parseGTF(opt.gtfFile, index, opt, guessChromosomes);
      }

      MinCollector collection(index, opt);
      MasterProcessor MP(index, opt, collection, model);
      if (batch_mode) {
        num_processed = ProcessBatchReads(MP, opt);
        writeCellIds(cellnamesfilename, opt.batch_ids);
        // Write out fake barcodes that identify each cell
        if ((opt.batch_mode && opt.technology.empty()) || opt.record_batch_bus_barcode) {
          std::vector<std::string> fake_bcs;
          for (size_t j = 0; j < MP.batch_id_mapping.size(); j++) {
            fake_bcs.push_back(binaryToString(MP.batch_id_mapping[j], BUSFORMAT_FAKE_BARCODE_LEN));
          }
          writeCellIds(busbarcodelistname, fake_bcs);
        }
        // Write out index if necessary (basically, when no UMIs or when paired-end)
        if (!opt.single_end || opt.technology.empty() || opt.busOptions.paired ||
            opt.busOptions.umi[0].fileno == -1) {
          index.write((opt.output + "/index.saved"), false, opt.threads);
        }
        // Write out fragment length distributions if reads paired-end or long:
        if (!opt.single_end || opt.long_read) {
          std::remove((opt.output + "/flens.txt").c_str());
          std::ofstream flensout_f((opt.output + "/flens.txt"));
          for (size_t id = 0; id < opt.batch_ids.size(); id++) {
            if (opt.long_read) {
              // Should I be using batchFlens?
              std::vector<uint32_t> fld_lr = MP.batchFlens_lr[id];
              std::vector<uint32_t> fld_lr_c = MP.batchFlens_lr_c[id];
              for (size_t i = 0; i < fld_lr.size(); ++i) {
                if (i != 0) {
                  flensout_f << " ";
                }
                if (fld_lr_c[i] > 0.5) {
                  // Good results with comment below.
                  // flensout_f << std::fabs((double)fld_lr[i] / (double)fld_lr_c[i] -
                  // index.k);//index.target_lens_[i] - (double)fld_lr[i] / (double)fld_lr_c[i] -
                  // k); // take mean of recorded uniquely aligning read lengths
                  flensout_f << std::fabs(((double)fld_lr[i] / (double)fld_lr_c[i]) - index.k);
                } else {
                  flensout_f << std::fabs(index.target_lens_[i] -
                                          index.k);  // index.target_lens_[i]);
                }
              }
              flensout_f << "\n";
            } else {
              std::vector<uint32_t> fld = MP.batchFlens[id];
              for (size_t i = 0; i < fld.size(); ++i) {
                if (i != 0) {
                  flensout_f << " ";
                }
                flensout_f << fld[i];
              }
              flensout_f << "\n";
            }
          }
          flensout_f.close();

          if (opt.unmapped) {
            std::ofstream um_f((opt.output + "/unmapped_ratio.txt"));
            for (size_t id = 0; id < opt.batch_ids.size(); id++) {
              std::vector<double> unmapped_l = MP.tc.unmapped_list;
              for (size_t i = 0; i < unmapped_l.size(); ++i) {
                if (i != 0) {
                  um_f << ",";
                }
                um_f << unmapped_l[i];
              }
              um_f << ",";
            }
            um_f.close();
          }
        }
      } else {
        num_processed = ProcessBUSReads(MP, opt);
        for (int i = 0; i <= 32; i++) {
          if (MP.bus_bc_len[i] > MP.bus_bc_len[bclen]) {
            bclen = i;
          }
          if (MP.bus_umi_len[i] > MP.bus_umi_len[umilen]) {
            umilen = i;
          }
        }

        bool write = false;
        // hack, open the bus file and write over the values in there.
        if (opt.busOptions.getBCLength() == 0) {
          if (bclen > 0) {
            write = true;
          }
        } else {
          bclen = opt.busOptions.getBCLength();
        }
        if (opt.busOptions.getUMILength() == 0) {
          if (umilen > 0) {
            write = true;
          }
        } else {
          umilen = opt.busOptions.getUMILength();
        }

        if (write) {
          std::FILE *fp = std::fopen((opt.output + "/output.bus").c_str(), "r+b");
          if (fp != nullptr) {
            std::fseek(fp, 8, SEEK_SET);  // skip magicstd::string and version
            // write to uint32_t values
            std::fwrite(&bclen, sizeof(bclen), 1, fp);
            std::fwrite(&umilen, sizeof(umilen), 1, fp);
            std::fclose(fp);
            fp = nullptr;
          }
        }
        std::vector<uint32_t> fld;
        std::vector<uint32_t> fld_c;
        if (opt.busOptions.paired && !opt.long_read) {
          fld = collection.flens;  // copy
          collection.compute_mean_frag_lens_trunc();
          // Write out index:
          index.write((opt.output + "/index.saved"), false, opt.threads);
          // Write out fragment length distribution:
          std::remove((opt.output + "/flens.txt").c_str());
          std::ofstream flensout_f((opt.output + "/flens.txt"));
          for (size_t i = 0; i < fld.size(); ++i) {
            if (i != 0) {
              flensout_f << " ";
            }
            flensout_f << fld[i];
          }
          flensout_f << "\n";
          flensout_f.close();
        } else if (opt.long_read) {
          opt.busOptions.threshold = opt.threshold;
          fld = collection.flens_lr;      // copy
          fld_c = collection.flens_lr_c;  // copy
          // Write out index:
          index.write((opt.output + "/index.saved"), false, opt.threads);
          // Write out fragment length distribution:
          std::remove((opt.output + "/flens.txt").c_str());
          std::ofstream flensout_f((opt.output + "/flens.txt"));
          for (size_t i = 0; i < fld.size(); ++i) {
            if (i > 0.0001) {
              flensout_f << " ";
            }
            if (fld_c[i] > 0.5) {
              flensout_f << std::fabs(
                  (double)fld[i] / (double)fld_c[i] -
                  index.k);  // index.target_lens_[i] - (double)fld[i] / (double)fld_c[i] - k); //
                             // take mean of recorded uniquely aligning read lengths
            } else {
              flensout_f << std::fabs(index.target_lens_[i] -
                                      index.k);  // index.target_lens_[i] - 31);
            }
          }
          flensout_f << "\n";
          flensout_f.close();

          if (opt.unmapped) {
            std::ofstream um_f((opt.output + "/unmapped_ratio.txt"));
            std::vector<double> unmapped_l = collection.unmapped_list;
            for (size_t i = 0; i < unmapped_l.size(); ++i) {
              if (i != 0) {
                um_f << ",";
              }
              um_f << unmapped_l[i];
            }
            um_f.close();
          }
        } else if (opt.busOptions.umi[0].fileno == -1) {
          // Write out index:
          index.write((opt.output + "/index.saved"), false, opt.threads);
        }
      }
      writeECList(ecfilename, index);

      // write transcript names
      std::ofstream transout_f((opt.output + "/transcripts.txt"));
      for (size_t i = 0; i < index.onlist_sequences.cardinality(); i++) {
        transout_f << index.target_names_[i] << "\n";
      }
      transout_f.close();

      for (const auto &elem : index.ecmapinv) {
        if (elem.first.cardinality() == 1) {
          num_unique += collection.counts[elem.second];
        }
        num_pseudoaligned += collection.counts[elem.second];
      }

      // write json file
      std::string call = argv_to_string(argc, argv);
      plaintext_aux(opt.output + "/run_info.json",
                    std::string(std::to_string(index.onlist_sequences.cardinality())),
                    std::string(std::to_string(0)), std::string(std::to_string(num_processed)),
                    std::string(std::to_string(num_pseudoaligned)),
                    std::string(std::to_string(num_unique)), KALLISTO_VERSION,
                    std::string(std::to_string(index.INDEX_VERSION)),
                    std::string(std::to_string(index.k)), start_time, call,
                    opt.aa ? std::to_string(collection.cardinality_clashes) : "");

      if (opt.pseudobam) {
        std::vector<double> fl_means(index.target_lens_.size(), 0.0);
        EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
#ifndef NO_HTSLIB
        MP.processAln(em, false);
#endif
      }

      std::cerr << std::endl;
      if (!opt.gtfFile.empty()) {
        // write out gene info
        writeGeneList(genelistname, model);
      }
      if (opt.max_num_reads != 0 && num_processed < opt.max_num_reads) {
        std::cerr << "Note: Number of reads processed is less than --numReads: "
                  << opt.max_num_reads << ", returning 1" << std::endl;
        return 1;
      }
      if (num_pseudoaligned == 0) {
        exit(1);  // exit with error
      }
    } else if (cmd == "merge") {
      std::cerr << "Deprecated: `kallisto merge` is deprecated. See `kallisto bus`." << std::endl;
    } else if (cmd == "quant") {
      if (argc == 2) {
        usageEM();
        return 0;
      }
      ParseOptionsEM(argc - 1, argv + 1, opt);
      if (!CheckOptionsEM(opt)) {
        std::cerr << std::endl;
        usageEM(false);
        exit(1);
      } else {
        // run the em algorithm
        KmerIndex index(opt);
        index.load(opt);
        if (opt.fusion) {
          // need full transcript sequences
          index.loadTranscriptSequences();
        }

        bool guessChromosomes = false;
        Transcriptome model;
        if (opt.genomebam) {
          if (!opt.chromFile.empty()) {
            model.loadChromosomes(opt.chromFile);
          } else {
            guessChromosomes = true;
          }
          model.parseGTF(opt.gtfFile, index, opt, guessChromosomes);
        }

        int64_t num_processed = 0;
        int64_t num_pseudoaligned = 0;
        int64_t num_unique = 0;

        MinCollector collection(index, opt);
        MasterProcessor MP(index, opt, collection, model);
        num_processed = ProcessReads(MP, opt);

        // save modified index for future use
        if (opt.write_index) {
          index.write((opt.output + "/index.saved"), false, opt.threads);
        }

        // if mean FL not provided, estimate
        std::vector<uint32_t> fld;
        if (opt.fld == 0.0) {
          fld = collection.flens;  // copy
          collection.compute_mean_frag_lens_trunc();
        } else {
          auto mean_fl = (opt.fld > 0.0) ? opt.fld : collection.get_mean_frag_len();
          auto sd_fl = opt.sd;
          collection.init_mean_fl_trunc(mean_fl, sd_fl);
          // fld.resize(MAX_FRAG_LEN,0); // no obersvations
          fld = trunc_gaussian_counts(0, MAX_FRAG_LEN, mean_fl, sd_fl, 10000);
        }

        std::vector<int32_t> preBias(4096, 1);
        if (opt.bias) {
          preBias = collection.bias5;  // copy
        }

        auto fl_means = get_frag_len_means(index.target_lens_, collection.mean_fl_trunc);

        EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
        if (opt.priors != "") {
          std::vector<double> priors = EMAlgorithm::read_priors(opt.priors);
          em.set_priors(priors);
          priors.clear();
        }
        em.run(10000, 50, true, opt.bias);

        std::string call = argv_to_string(argc, argv);

#ifdef USE_HDF5
        H5Writer writer;
        if (!opt.plaintext) {
          std::vector<int> fld_(std::begin(fld), std::end(fld));
          writer.init(opt.output + "/abundance.h5", opt.bootstrap, num_processed, fld_, preBias,
                      em.post_bias_, 6, index.INDEX_VERSION, call, start_time);
          std::vector<int> lengths_2(std::begin(index.target_lens_), std::end(index.target_lens_));
          writer.write_main(em, index.target_names_, lengths_2);
        }
#endif  // USE_HDF5

        for (const auto &elem : index.ecmapinv) {
          if (elem.first.cardinality() == 1) {
            num_unique += collection.counts[elem.second];
          }
          num_pseudoaligned += collection.counts[elem.second];
        }

        if (num_pseudoaligned == 0) {
          std::cerr << "[~warn] Warning, zero reads pseudoaligned check your input files and index"
                    << std::endl;
        }

        plaintext_aux(
            opt.output + "/run_info.json",
            std::string(std::to_string(index.onlist_sequences.cardinality())),
            std::string(std::to_string(opt.bootstrap)), std::string(std::to_string(num_processed)),
            std::string(std::to_string(num_pseudoaligned)), std::string(std::to_string(num_unique)),
            KALLISTO_VERSION, std::string(std::to_string(index.INDEX_VERSION)),
            std::string(std::to_string(index.k)), start_time, call,
            opt.aa ? std::to_string(collection.cardinality_clashes) : "");

        plaintext_writer(opt.output + "/abundance.tsv", em.target_names_, em.alpha_, em.eff_lens_,
                         index.target_lens_);

        {
          std::ofstream ecf(opt.output + "/counts_cpu.txt");
          for (const auto& it : index.ecmapinv) {
            int ec_id = it.second;
            ecf << ec_id << "\t" << collection.counts[ec_id] << "\t";
            bool first = true;
            for (uint32_t tr : it.first) {
              if (!first) ecf << ",";
              ecf << tr;
              first = false;
            }
            ecf << "\n";
          }
          ecf.close();
        }

        if (opt.bootstrap > 0 && num_pseudoaligned == 0) {
          // this happens if nothing aligns, then we write an empty bootstrap file
          for (int b = 0; b < opt.bootstrap; b++) {
            if (!opt.plaintext) {
#ifdef USE_HDF5
              writer.write_bootstrap(em, b);  // em is empty
#endif
            } else {
              plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".tsv",
                               em.target_names_, em.alpha_, em.eff_lens_,
                               index.target_lens_);  // re-use empty input
            }
          }
        } else if (opt.bootstrap > 0 && num_pseudoaligned > 0) {
          auto B = opt.bootstrap;
          std::mt19937_64 rand;
          rand.seed(opt.seed);

          std::vector<size_t> seeds;
          for (auto s = 0; s < B; ++s) {
            seeds.push_back(rand());
          }

          if (opt.threads > 1) {
            auto n_threads = opt.threads;
            if (opt.threads > opt.bootstrap) {
              std::cerr << "[~warn] number of threads (" << opt.threads
                   << ") greater than number of bootstraps" << std::endl
                   << "[~warn] (cont'd) updating threads to number of bootstraps " << opt.bootstrap
                   << std::endl;
              n_threads = opt.bootstrap;
            }
#ifdef USE_HDF5
            BootstrapThreadPool pool(opt.threads, seeds, collection.counts, index, collection,
                                     em.eff_lens_, opt, &writer, fl_means);
#endif
          } else {
            for (auto b = 0; b < B; ++b) {
              Bootstrap bs(collection.counts, index, collection, em.eff_lens_, seeds[b], fl_means,
                           opt);
              std::cerr << "[bstrp] running EM for the bootstrap: " << b + 1 << "\r";
              auto res = bs.run_em();

              if (!opt.plaintext) {
#ifdef USE_HDF5
                writer.write_bootstrap(res, b);
#endif
              } else {
                plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".tsv",
                                 em.target_names_, res.alpha_, em.eff_lens_, index.target_lens_);
              }
            }
          }

          std::cerr << std::endl;
        }

        if (opt.pseudobam) {
#ifndef NO_HTSLIB
          MP.processAln(em, true);
#endif
        }

        std::cerr << std::endl;
        if (num_pseudoaligned == 0) {
          exit(1);  // exit with error
        }
      }
    } else if (cmd == "quant-only") {
      std::cerr << "[deprecated] quant-only option has been deprecated." << std::endl;
      exit(0);
    } else if (cmd == "quant-tcc") {
      if (argc == 2) {
        usageTCCQuant();
        return 0;
      }
      ParseOptionsTCCQuant(argc - 1, argv + 1, opt);
      if (!CheckOptionsTCCQuant(opt)) {
        usageTCCQuant();
        exit(1);
      } else {
        KmerIndex index(opt);
        index.load(opt, false, false);  // skip the k-mer map and the D-list
        MinCollector collection(index, opt);
        std::vector<std::vector<std::pair<uint32_t, uint32_t>>> batchCounts;  // Stores TCCs
        std::ifstream in((opt.tccFile));
        bool firstline = true;
        bool isMatrixFile = false;
        size_t nrow = 0, ncol = 0, nlines = 0;
        int i = 0, prevRow = 0, prevCol = 0;
        if (in.is_open()) {  // Read in TCC file
          std::string line;
          while (getline(in, line)) {
            if (firstline) {  // First line decides whether the file is in matrix format
              firstline = false;
              if (line.rfind("%%MatrixMarket", 0) == 0) {
                std::cerr
                    << "[tcc] Parsing transcript-compatibility counts (TCC) file as a matrix file"
                    << std::endl;
                isMatrixFile = true;
                while (getline(in, line) && line.rfind("%", 0) == 0) {
                }  // Skip header comments
                std::stringstream ss(line);
                ss >> nrow >> ncol >> nlines;
                std::cerr << "[tcc] Matrix dimensions: " << pretty_num(nrow) << " x "
                          << pretty_num(ncol) << std::endl;
                batchCounts.assign(nrow, {});
                if (opt.bootstrap > 0) {
                  if (opt.plaintext) {
                    std::cerr << "[tcc] Bootstrapping will be performed and outputted as plaintext"
                         << std::endl;
                  } else {
#ifdef USE_HDF5
                    std::cerr << "[tcc] Bootstrapping will be performed and outputted as HDF5" << std::endl;
#else
                    std::cerr << "[tcc] Bootstrapping will NOT be performed because kallisto was "
                            "compiled without HDF5 support. ";
                    std::cerr << "Bootstrapping can be outputted as plaintext if the --plaintext option "
                            "is supplied."
                         << std::endl;
#endif
                  }
                }
                continue;
              } else {
                std::cerr
                    << "[tcc] Transcript-compatibility counts (TCC) file is not in matrix format; "
                    << "it will not be parsed as a matrix file" << std::endl;
                isMatrixFile = false;
                batchCounts.assign(1, {});
              }
            }
            std::stringstream ss(line);
            int row, col, val;
            if (isMatrixFile) {
              if (i >= nlines) {
                std::cerr << "[tcc] Warning: TCC matrix file contains additional lines which will "
                             "not be read; "
                          << "only " << pretty_num(nlines)
                          << " entries, as specified on the first line, will be read." << std::endl;
                break;
              }
              ss >> row >> col >> val;
              if (row > nrow || col > ncol) {
                std::cerr << "Error: TCC matrix file is malformed; "
                          << "row numbers or column numbers exceed the dimensions of the matrix."
                          << std::endl;
                exit(1);
              }
            } else {
              ss >> col >> val;
              col += 1;  // Because we'll consider non-matrix TCC files to be zero-indexed
              row = 1;
              nrow = 1;
              if (ncol < col) {
                ncol = col;
              }
            }
            if (row <= 0 || col <= 0) {
              std::cerr << "Error: Invalid indices in TCC file." << std::endl;
              exit(1);
            }
            if (row < prevRow || (row == prevRow && col <= prevCol)) {
              std::cerr << "Error: TCC file is not sorted." << std::endl;
              exit(1);
            }
            prevRow = row;
            prevCol = col;
            auto &bc = batchCounts[row - 1];
            bc.push_back({col - 1, val});
            i++;
          }
        } else {
          std::cerr << "Error: could not open file " << opt.tccFile << std::endl;
          exit(1);
        }
        if (isMatrixFile && i < nlines) {
          std::cerr << "Error: Found only " << pretty_num(i) << " entries in TCC matrix file, "
                    << "expected " << pretty_num(nlines) << std::endl;
          exit(1);
        }

        std::vector<std::vector<std::pair<int32_t, double>>> Abundance_mat, Abundance_mat_gene,
            TPM_mat, TPM_mat_gene, EffLen_mat;
        Abundance_mat.resize(nrow, {});
        Abundance_mat_gene.resize(nrow, {});
        TPM_mat.resize(nrow, {});
        TPM_mat_gene.resize(nrow, {});
        std::vector<std::pair<double, double>> FLD_mat;
        FLD_mat.resize(nrow, {});
        EffLen_mat.resize(nrow, {});
        std::vector<std::vector<uint32_t>> FLDs;
        Transcriptome model;  // empty model

        if (index.onlist_sequences.cardinality() > 0) {
          std::ofstream transout_f((opt.output + "/transcripts.txt"));
          for (size_t i = 0; i < index.onlist_sequences.cardinality(); i++) {
            transout_f << index.target_names_[i] << "\n";
          }
          transout_f.close();
        }
        std::string prefix = opt.output + "/matrix";
        std::string abtsvfilename = opt.output + "/abundance.tsv";
        std::string abtsvprefix = opt.output + "/abundance";
        std::string abtsvprefixh5 = opt.output + "/abundance";
        std::string bootstrapprefix = opt.output + "/bs_abundance";
        std::string gene_abtsvfilename = opt.output + "/abundance.gene.tsv";
        std::string genelistname = opt.output + "/genes.txt";
        std::string abfilename = prefix + ".abundance.mtx";
        std::string abtpmfilename = prefix + ".abundance.tpm.mtx";
        std::string efflenmtxfilename = prefix + ".efflens.mtx";
        std::string gene_abfilename = prefix + ".abundance.gene.mtx";
        std::string gene_abtpmfilename = prefix + ".abundance.gene.tpm.mtx";
        std::string fldfilename = prefix + ".fld.tsv";

        size_t num_trans = index.target_names_.size();

        const bool calcEffLen = !opt.fldFile.empty() || opt.fld != 0.0;
        if (calcEffLen && !opt.fldFile.empty() &&
            (!opt.long_read ||
             opt.platform == "PACBIO")) {  // Parse supplied fragment length distribution file
          std::ifstream infld((opt.fldFile));
          if (infld.is_open()) {
            std::string line;
            while (getline(infld, line)) {
              if (line.empty() || line.rfind("#", 0) == 0) {
                continue;  // Ignore empty lines or lines that begin with a #
              }
              std::vector<uint32_t> tmp_vec;
              std::stringstream ss(line);
              while (ss.good()) {
                std::string tmp_val;
                getline(ss, tmp_val, ' ');
                int tmp_val_num = std::atoi(tmp_val.c_str());
                if (tmp_val_num < 0) {
                  std::cerr << "Error: Fragment length distribution file contains invalid value: "
                            << tmp_val_num << std::endl;
                  exit(1);
                }
                tmp_vec.push_back(tmp_val_num);
              }
              if (tmp_vec.size() != MAX_FRAG_LEN && tmp_vec.size() != index.target_lens_.size()) {
                std::cerr << "Error: Fragment length distribution file contains a line with "
                          << tmp_vec.size() << " values; expected: " << MAX_FRAG_LEN << std::endl;
                exit(1);
              }
              FLDs.push_back(tmp_vec);
            }
            if (FLDs.size() != 1 &&
                (FLDs.size() != nrow && FLDs.size() != index.target_lens_.size())) {
              std::cerr << "Error: Fragment length distribution file contains " << FLDs.size()
                        << " valid lines; expected: " << nrow << std::endl;
              exit(1);
            }
          } else {
            std::cerr << "Error: could not open file " << opt.fldFile << std::endl;
            exit(1);
          }
        }

        if (!opt.genemap.empty()) {  // Parse supplied gene mapping file
          model.parseGeneMap(opt.genemap, index, opt);
        } else if (!opt.gtfFile.empty()) {
          model.parseGTF(opt.gtfFile, index, opt, true);
        }
        const bool gene_level_counting = !opt.genemap.empty() || !opt.gtfFile.empty();

        std::cerr << "[quant] Running EM algorithm..." << std::endl;

        std::vector<double> priors;
        if (opt.priors != "") {
          priors = EMAlgorithm::read_priors(opt.priors);
        }
        auto EM_lambda = [&](int id) {
          std::cerr << "[quant] Processing sample/cell " << std::to_string(id) << std::endl;
          MinCollector collection(index, opt);
          const auto &bc = batchCounts[id];
          for (const auto &p : bc) {
            collection.counts[p.first] = p.second;
          }

          std::vector<double> fl_means;
          if (calcEffLen && (!opt.long_read || !(opt.platform == "ONT"))) {
            if (opt.fld != 0.0) {
              collection.init_mean_fl_trunc(opt.fld, opt.sd);
            } else {
              std::vector<uint32_t> fld;
              if (FLDs.size() ==
                  1) {  // Only one distribution supplied; will use this for everything
                fld = FLDs[0];
              } else {
                fld = FLDs[id];
              }
              collection.flens = fld;
              collection.compute_mean_frag_lens_trunc(false);
            }
            fl_means = get_frag_len_means(index.target_lens_, collection.mean_fl_trunc);
            double mean_fl = collection.get_mean_frag_len(true);
            double sd_fl = collection.get_sd_frag_len();
            FLD_mat[id] = {mean_fl, sd_fl};
          } else {  // Set mean fragment lengths to be target lengths so that effective lengths are
                    // all one
            fl_means = std::vector<double>(index.target_lens_.begin(), index.target_lens_.end());
          }

          EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
          if (opt.priors != "") {
            em.set_priors(priors);
          }
          em.run(10000, 50, false, false);

          if (isMatrixFile) {  // Update abundances matrix
            auto &ab_m = Abundance_mat[id];
            auto &tpm_m = TPM_mat[id];
            auto &elen_m = EffLen_mat[id];
            std::vector<double> gc;
            std::vector<double> gc_tpm;
            if (gene_level_counting) {
              gc.assign(model.genes.size(), 0.0);
              gc_tpm.assign(model.genes.size(), 0.0);
            }
            auto tpm = counts_to_tpm(em.alpha_, em.eff_lens_);
            for (int i = 0; i < em.alpha_.size(); i++) {
              if (em.alpha_[i] > 0.0) {
                ab_m.push_back({i, em.alpha_[i]});
                tpm_m.push_back({i, tpm[i]});
                if (calcEffLen) {
                  elen_m.push_back({i, em.eff_lens_[i]});
                }
                if (gene_level_counting) {
                  int g_id = model.transcripts[i].gene_id;
                  if (g_id != -1) {
                    gc[g_id] += em.alpha_[i];
                    gc_tpm[g_id] += tpm[i];
                  }
                }
              }
            }
            if (gene_level_counting) {
              auto &ab_m_g = Abundance_mat_gene[id];
              auto &tpm_m_g = TPM_mat_gene[id];
              for (int i = 0; i < gc.size(); i++) {
                if (gc[i] > 0.0) {
                  ab_m_g.push_back({i, gc[i]});
                  tpm_m_g.push_back({i, gc_tpm[i]});
                }
              }
            }
            // Write plaintext abundance for current row in matrix
            if (opt.matrix_to_files) {
              std::string output_dir = "";

              if (opt.matrix_to_directories) {
                output_dir = opt.output + "/abundance_" + std::to_string(id + 1);
                struct stat stFileInfo;
                auto intStat = stat(output_dir.c_str(), &stFileInfo);
                if (intStat == 0) {
                  // file/dir exits
                  if (!S_ISDIR(stFileInfo.st_mode)) {
                    std::cerr << "Error: file " << output_dir << " exists and is not a directory"
                         << std::endl;
                    exit(1);
                  }
                } else {
                  // create directory
                  if (my_mkdir(output_dir.c_str(), 0777) == -1) {
                    std::cerr << "Error: could not create directory " << output_dir << std::endl;
                    exit(1);
                  }
                }
                output_dir = output_dir + "/";
                abtsvprefix = "abundance";
                abtsvprefixh5 = "abundance";
                bootstrapprefix = "bs_abundance";
              }
              std::string id_suffix = (output_dir.empty() ? ("_" + std::to_string(id + 1)) : "");

              plaintext_writer(output_dir + abtsvprefix + id_suffix + ".tsv", em.target_names_,
                               em.alpha_, em.eff_lens_, index.target_lens_);
              if (gene_level_counting) {
                plaintext_writer_gene(output_dir + abtsvprefix + ".gene" + id_suffix + ".tsv",
                                      em.target_names_, em.alpha_, em.eff_lens_, model);
              }
#ifdef USE_HDF5
              H5Writer writer;
              std::vector<uint32_t> fld_i;
              if (!opt.plaintext) {
                fld_i = collection.flens;
                if (opt.fld != 0.0) {
                  fld_i = trunc_gaussian_counts(0, MAX_FRAG_LEN, opt.fld, opt.sd, 10000);
                }
                std::vector<int> fld_(std::begin(fld_i), std::end(fld_i));
                std::vector<int32_t> preBias(4096, 1);
                writer.init(output_dir + abtsvprefixh5 + id_suffix + ".h5", opt.bootstrap, 0, fld_,
                            preBias, em.post_bias_, 6, index.INDEX_VERSION, "", "");
                std::vector<int> lengths_2(std::begin(index.target_lens_),
                                           std::end(index.target_lens_));
                writer.write_main(em, index.target_names_, lengths_2);
              }
#endif  // USE_HDF5
        // Write out bootstraps if requested:
              if (opt.bootstrap > 0) {
                if (ab_m.size() == 0) {  // No abundances (e.g. nothing aligned), write empty
                  for (int b = 0; b < opt.bootstrap; b++) {
                    if (!opt.plaintext) {
#ifdef USE_HDF5
                      writer.write_bootstrap(em, b);  // em is empty
#endif                                                // USE_HDF5
                    } else {
                      plaintext_writer(output_dir + bootstrapprefix + id_suffix + "_" +
                                           std::to_string(b) + ".tsv",
                                       em.target_names_, em.alpha_, em.eff_lens_,
                                       index.target_lens_);
                      if (gene_level_counting) {
                        plaintext_writer_gene(output_dir + bootstrapprefix + ".gene" + id_suffix +
                                                  "_" + std::to_string(b) + ".tsv",
                                              em.target_names_, em.alpha_, em.eff_lens_, model);
                      }
                    }
                  }
                } else {
                  auto B = opt.bootstrap;
                  std::mt19937_64 rand;
                  rand.seed(opt.seed);
                  std::vector<size_t> seeds;
                  for (auto s = 0; s < B; ++s) {
                    seeds.push_back(rand());
                  }
                  for (auto b = 0; b < B; ++b) {
                    Bootstrap bs(collection.counts, index, collection, em.eff_lens_, seeds[b],
                                 fl_means, opt);
                    auto res = bs.run_em();
                    if (!opt.plaintext) {
#ifdef USE_HDF5
                      writer.write_bootstrap(res, b);
#endif  // USE_HDF5
                    } else {
                      plaintext_writer(output_dir + bootstrapprefix + id_suffix + "_" +
                                           std::to_string(b) + ".tsv",
                                       em.target_names_, res.alpha_, em.eff_lens_,
                                       index.target_lens_);
                      if (gene_level_counting) {
                        plaintext_writer_gene(output_dir + bootstrapprefix + ".gene" + id_suffix +
                                                  "_" + std::to_string(b) + ".tsv",
                                              em.target_names_, res.alpha_, em.eff_lens_, model);
                      }
                    }
                  }
                }
              }     // End bootstrapping
            }       // End matrix-to-files
          } else {  // Write plaintext abundances (for non-matrix files)
            plaintext_writer(abtsvfilename, em.target_names_, em.alpha_, em.eff_lens_,
                             index.target_lens_);
            if (gene_level_counting) {
              plaintext_writer_gene(gene_abtsvfilename, em.target_names_, em.alpha_, em.eff_lens_,
                                    model);
            }
          }

          if (opt.bootstrap > 0 && !isMatrixFile) {
            auto B = opt.bootstrap;
            std::mt19937_64 rand;
            rand.seed(opt.seed);
            std::vector<size_t> seeds;
            for (auto s = 0; s < B; ++s) {
              seeds.push_back(rand());
            }
            std::cerr << std::endl;
            for (auto b = 0; b < B; ++b) {
              Bootstrap bs(collection.counts, index, collection, em.eff_lens_, seeds[b], fl_means,
                           opt);
              std::cerr << "[bstrp] running EM for the bootstrap: " << b + 1 << "\r";
              auto res = bs.run_em();
              plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".tsv",
                               em.target_names_, res.alpha_, em.eff_lens_, index.target_lens_);
            }
            std::cerr << std::endl;
          }
        };  // end of EM_lambda

        std::vector<std::thread> workers;
        int num_ids = nrow;
        int id = 0;
        while (id < num_ids) {
          workers.clear();
          int nt = std::min(opt.threads, (num_ids - id));
          int first_id = id;
          for (int i = 0; i < nt; i++, id++) {
            workers.emplace_back(std::thread(EM_lambda, id));
          }
          for (int i = 0; i < nt; i++) {
            workers[i].join();
          }
        }
        std::cerr << " done" << std::endl;

        std::cerr << std::endl;

        if (isMatrixFile) {
          writeSparseBatchMatrix(abfilename, Abundance_mat, num_trans);
          writeSparseBatchMatrix(abtpmfilename, TPM_mat, num_trans);
          if (calcEffLen) {
            writeSparseBatchMatrix(efflenmtxfilename, EffLen_mat, num_trans);
          }
          if (gene_level_counting) {
            writeSparseBatchMatrix(gene_abfilename, Abundance_mat_gene, model.genes.size());
            writeSparseBatchMatrix(gene_abtpmfilename, TPM_mat_gene, model.genes.size());
            writeGeneList(genelistname, model, true);
          }
        }
        if (calcEffLen) {
          writeFLD(fldfilename, FLD_mat);
          std::ofstream translens_f((opt.output + "/transcript_lengths.txt"));
          for (size_t i = 0; i < num_trans; i++) {
            translens_f << index.target_names_[i] << " " << index.target_lens_[i] << "\n";
          }
          translens_f.close();
        }
      }
    } else if (cmd == "pseudo") {
      std::cerr << "Deprecated: `kallisto pseudo` is deprecated. See `kallisto bus`." << std::endl;
    } else if (cmd == "h5dump") {
      if (argc == 2) {
        usageh5dump();
        exit(1);
      }

      ParseOptionsH5Dump(argc - 1, argv + 1, opt);
      if (!CheckOptionsH5Dump(opt)) {
        usageh5dump();
        exit(1);
      }
#ifdef USE_HDF5
      H5Converter h5conv(opt.files[0], opt.output);
      if (!opt.peek) {
        h5conv.write_aux();
        h5conv.convert();
      }
#endif
    } else {
      std::cerr << "Error: invalid command " << cmd << std::endl;
      usage();
      exit(1);
    }
  }

  fflush(stdout);

  return 0;
}
