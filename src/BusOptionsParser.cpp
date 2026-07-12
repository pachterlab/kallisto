#include "BusOptionsParser.h"

#include <getopt.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace bus_parser {

void ListSingleCellTechnologies() {
  std::cout << "List of supported single-cell technologies" << std::endl
            << std::endl
            << "short name       description" << std::endl
            << "----------       -----------" << std::endl
            << "10xv1            10x version 1 chemistry" << std::endl
            << "10xv2            10x version 2 chemistry" << std::endl
            << "10xv3            10x version 3 chemistry" << std::endl
            << "10xv4            10x version 4 chemistry" << std::endl
            << "Bulk             Bulk RNA-seq" << std::endl
            << "ParseV3          Parse Evercode V3" << std::endl
            << "SmartSeq2        Smart-seq2 (multiplexed)" << std::endl
            << "BDWTA            BD Rhapsody WTA" << std::endl
            << "CELSeq           CEL-Seq" << std::endl
            << "CELSeq2          CEL-Seq version 2" << std::endl
            << "DropSeq          DropSeq" << std::endl
            << "inDropsv1        inDrops version 1 chemistry" << std::endl
            << "inDropsv2        inDrops version 2 chemistry" << std::endl
            << "inDropsv3        inDrops version 3 chemistry" << std::endl
            << "MATQSEQ          MATQ-SEQ" << std::endl
            << "PETRISEQ         PETRI-SEQ" << std::endl
            << "SCRBSeq          SCRB-Seq" << std::endl
            << "SmartSeq3        Smart-seq3" << std::endl
            << "SPLiT-seq        SPLiT-seq" << std::endl
            << "STORM-seq        STORM-seq" << std::endl
            << "SureCell         SureCell for ddSEQ" << std::endl
            << "VASA-seq         VASA-seq" << std::endl
            << "Visium           10x Visium Spatial Transcriptomics" << std::endl
            << std::endl;
}

void ParseOptionsBus(int argc, char **argv, ProgramOptions &opt) {
  int verbose_flag = 0;
  int gbam_flag = 0;
  int paired_end_flag = 0;
  int long_read_flag = 0;
  int unmapped_flag = 0;
  int aa_flag = 0;
  int strand_FR_flag = 0;
  int strand_RF_flag = 0;
  int unstranded_flag = 0;
  int interleaved_flag = 0;
  int batch_barcodes_flag = 0;
  int dfk_onlist_flag = 0;
  int do_union_flag = 0;
  int no_jump_flag = 0;

  const char *opt_string = "i:o:x:t:lbng:c:T:P:r:e:B:N:";
  static struct option long_options[] = {{"verbose", no_argument, &verbose_flag, 1},
                                         {"dfk-onlist", no_argument, &dfk_onlist_flag, 1},
                                         {"index", required_argument, 0, 'i'},
                                         {"output-dir", required_argument, 0, 'o'},
                                         {"technology", required_argument, 0, 'x'},
                                         {"list", no_argument, 0, 'l'},
                                         {"batch", required_argument, 0, 'B'},
                                         {"threads", required_argument, 0, 't'},
                                         {"bam", no_argument, 0, 'b'},
                                         {"num", no_argument, 0, 'n'},
                                         {"genomebam", no_argument, &gbam_flag, 1},
                                         {"gtf", required_argument, 0, 'g'},
                                         {"chromosomes", required_argument, 0, 'c'},
                                         {"tag", required_argument, 0, 'T'},
                                         {"fr-stranded", no_argument, &strand_FR_flag, 1},
                                         {"rf-stranded", no_argument, &strand_RF_flag, 1},
                                         {"unstranded", no_argument, &unstranded_flag, 1},
                                         {"paired", no_argument, &paired_end_flag, 1},
                                         {"long", no_argument, &long_read_flag, 1},
                                         {"platform", required_argument, 0, 'P'},
                                         {"threshold", required_argument, 0, 'r'},
                                         {"unmapped", no_argument, &unmapped_flag, 1},
                                         {"aa", no_argument, &aa_flag, 1},
                                         {"inleaved", no_argument, &interleaved_flag, 1},
                                         {"numReads", required_argument, 0, 'N'},
                                         {"batch-barcodes", no_argument, &batch_barcodes_flag, 1},
                                         {"union", no_argument, &do_union_flag, 1},
                                         {"no-jump", no_argument, &no_jump_flag, 1},
                                         {0, 0, 0, 0}};

  int list_flag = 0;
  int c;
  int option_index = 0;

  optind = 1;
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
      case 'l': {
        list_flag = 1;
        break;
      }
      case 'o': {
        opt.output = optarg;
        break;
      }
      case 'x': {
        opt.technology = optarg;
        std::transform(opt.technology.begin(), opt.technology.end(), opt.technology.begin(),
                       ::toupper);
        break;
      }
      case 't': {
        std::stringstream(optarg) >> opt.threads;
        break;
      }
      case 'b': {
        opt.bam = true;
        break;
      }
      case 'B': {
        opt.batch_mode = true;
        opt.batch_file_name = optarg;
        break;
      }
      case 'n': {
        opt.num = true;
        break;
      }
      case 'P': {
        std::stringstream(optarg) >> opt.platform;
        std::transform(opt.platform.begin(), opt.platform.end(), opt.platform.begin(), ::toupper);
        break;
      }
      case 'e': {
        std::stringstream(optarg) >> opt.error_rate;
        break;
      }
      case 'r': {
        std::stringstream(optarg) >> opt.threshold;
        break;
      }
      case 'N': {
        std::stringstream(optarg) >> opt.max_num_reads;
        if (opt.max_num_reads == 0) {
          opt.max_num_reads = -1;
        }
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
      case 'T': {
        std::stringstream(optarg) >> opt.tagsequence;
        break;
      }
      default:
        break;
    }
  }

  if (list_flag) {
    ListSingleCellTechnologies();
    exit(1);
  }

  if (opt.technology.find('%') !=
      std::string::npos) {  // Process technology strings of format -x bc:umi:cdna%strand%parity
    std::string first = opt.technology.substr(opt.technology.find("%") + 1);
    if (first.length() >= 7 && first.substr(0, 7) == "FORWARD") {
      opt.strand_specific = true;
      opt.strand = ProgramOptions::StrandType::FR;
    } else if (first.length() >= 7 && first.substr(0, 7) == "REVERSE") {
      opt.strand_specific = true;
      opt.strand = ProgramOptions::StrandType::RF;
    }
    if (first.find('%') != std::string::npos) {
      std::string second = first.substr(first.find("%") + 1);
      if (second.length() >= 6 && second.substr(0, 6) == "PAIRED") {
        opt.single_end = false;
        paired_end_flag = true;
      } else {
        opt.single_end = true;
      }
    }
    opt.technology = opt.technology.substr(0, opt.technology.find("%"));
  }

  if (verbose_flag) {
    opt.verbose = true;
  }

  if (gbam_flag) {
    opt.pseudobam = true;
    opt.genomebam = true;
  }

  if (strand_FR_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::FR;
  }

  if (strand_RF_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::RF;
  }

  if (unstranded_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::None;
  }

  if (paired_end_flag) {
    opt.single_end = false;
  } else {
    opt.single_end = true;
  }

  if (long_read_flag) {
    opt.long_read = true;
  }

  if (unmapped_flag) {
    opt.unmapped = true;
  }

  if (interleaved_flag) {
    opt.input_interleaved_nfiles = 1;
  }

  if (batch_barcodes_flag) {
    opt.record_batch_bus_barcode = true;
  }

  if (dfk_onlist_flag) {
    opt.dfk_onlist = true;
  }

  if (do_union_flag) {
    opt.do_union = true;
  }

  if (no_jump_flag) {
    opt.no_jump = true;
  }

  opt.single_overhang = true;

  if (aa_flag) {
    opt.aa = true;
    opt.dfk_onlist = true;
    opt.single_end = true;
    if (paired_end_flag) {
      std::cerr << "[bus] --paired ignored; --aa only supports single-end reads" << std::endl;
    }
  }

  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
}

bool ParseTechnology(const std::string &techstr, BUSOptions &busopt,
                     std::vector<std::string> &errorList) {
  auto i1 = techstr.find(':');
  if (i1 == std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), none found: \"" +
                        techstr + "\"");
    return false;
  }
  auto i2 = techstr.find(':', i1 + 1);
  if (i2 == std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), only one found: \"" +
                        techstr + "\"");
    return false;
  }
  auto ip = techstr.find(':', i2 + 1);
  if (ip != std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), three found: \"" +
                        techstr + "\"");
    return false;
  }
  auto bcstr = techstr.substr(0, i1);
  auto umistr = techstr.substr(i1 + 1, i2 - i1 - 1);
  auto seqstr = techstr.substr(i2 + 1);

  int maxnf = 0;

  auto convert_commas_to_vector = [&](const std::string &s,
                                      std::vector<BUSOptionSubstr> &v) -> bool {
    std::vector<int> vv;
    v.clear();
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, ',')) {
      try {
        int i = stoi(t);
        vv.push_back(i);
      } catch (std::invalid_argument &e) {
        errorList.push_back("Error: converting to int: \"" + t + "\"");
        return false;
      }
    }

    int nv = vv.size();
    if (nv % 3 == 0) {
      for (int i = 0; i + 2 < nv; i += 3) {
        int f = vv[i];
        int a = vv[i + 1];
        int b = vv[i + 2];
        if (f < -1) {
          errorList.push_back("Error: invalid file number (" + std::to_string(f) + ")  " + s);
          return false;
        }
        if (a < 0 && f != -1) {
          errorList.push_back("Error: invalid start (" + std::to_string(a) + ")  " + s);
          return false;
        }
        if (b != 0 && b <= a && f != -1) {
          errorList.push_back("Error: invalid stop (" + std::to_string(b) + ") has to be after start (" +
                              std::to_string(a) + ")  " + s);
          return false;
        }
        v.push_back(BUSOptionSubstr(f, a, b));
        if (f > maxnf) {
          maxnf = f;
        }
      }
    } else {
      errorList.push_back("Error: number of values has to be multiple of 3 " + s);
      return false;
    }

    busopt.nfiles = maxnf + 1;
    return true;
  };

  std::vector<BUSOptionSubstr> v;
  if (!convert_commas_to_vector(bcstr, v)) {
    return false;
  }
  if (v.empty()) {
    errorList.push_back("Error: empty barcode list " + bcstr);
    return false;
  }
  busopt.bc = std::move(v);

  if (umistr == "RX" || busopt.keep_fastq_comments) {
    busopt.keep_fastq_comments = true;
    v.push_back(BUSOptionSubstr(-1, -1, -1));
  } else if (!convert_commas_to_vector(umistr, v)) {
    return false;
  }
  if (v.empty()) {
    errorList.push_back("Error: empty UMI list " + umistr);
    return false;
  }
  busopt.umi = std::move(v);

  if (!convert_commas_to_vector(seqstr, v)) {
    return false;
  }
  if (v.empty()) {
    errorList.push_back("Error: empty sequence list " + bcstr);
    return false;
  }

  busopt.seq = std::move(v);

  return true;
}

bool ApplyBusPresetTechnology(ProgramOptions &opt,
                              ProgramOptions::StrandType &strand_out,
                              std::vector<std::string> &errorList) {
  auto &busopt = opt.busOptions;
  strand_out = ProgramOptions::StrandType::None;

  if (opt.technology == "10XV2") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 16, 26));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 16));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "10XV3") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 16, 28));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 16));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "PARSEV3") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(0, 0, 0));
    busopt.bc.push_back(BUSOptionSubstr(1, 10, 18));
    busopt.bc.push_back(BUSOptionSubstr(1, 30, 38));
    busopt.bc.push_back(BUSOptionSubstr(1, 50, 58));
    busopt.umi.push_back(BUSOptionSubstr(1, 0, 10));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "MATQSEQ") {
    busopt.nfiles = 1;
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 8));
    busopt.seq.push_back(BUSOptionSubstr(0, 8, 0));
    busopt.umi.push_back(BUSOptionSubstr(-1, -1, -1));
    strand_out = ProgramOptions::StrandType::None;
    opt.strand_specific = false;
  } else if (opt.technology == "PETRISEQ") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 17));
    busopt.bc.push_back(BUSOptionSubstr(0, 7, 14));
    busopt.bc.push_back(BUSOptionSubstr(0, 29, 36));
    busopt.bc.push_back(BUSOptionSubstr(0, 50, 58));
    busopt.umi.push_back(BUSOptionSubstr(0, 0, 7));
    strand_out = ProgramOptions::StrandType::None;
    opt.strand_specific = false;
  } else if (opt.technology == "VISIUM") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 16, 28));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 16));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "10XV1") {
    busopt.nfiles = 3;
    busopt.seq.push_back(BUSOptionSubstr(2, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(1, 0, 10));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 14));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "SURECELL") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 51, 59));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 6));
    busopt.bc.push_back(BUSOptionSubstr(0, 21, 27));
    busopt.bc.push_back(BUSOptionSubstr(0, 42, 48));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "DROPSEQ") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 12, 20));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 12));
    opt.strand_specific = true;
  } else if (opt.technology == "INDROPSV1") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 42, 48));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 11));
    busopt.bc.push_back(BUSOptionSubstr(0, 30, 38));
  } else if (opt.technology == "INDROPSV2") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(0, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(1, 42, 48));
    busopt.bc.push_back(BUSOptionSubstr(1, 0, 11));
    busopt.bc.push_back(BUSOptionSubstr(1, 30, 38));
  } else if (opt.technology == "INDROPSV3") {
    busopt.nfiles = 3;
    busopt.seq.push_back(BUSOptionSubstr(2, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(1, 8, 14));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 8));
    busopt.bc.push_back(BUSOptionSubstr(1, 0, 8));
  } else if (opt.technology == "CELSEQ") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 8, 12));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 8));
    strand_out = ProgramOptions::StrandType::FR;
  } else if (opt.technology == "CELSEQ2") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 0, 6));
    busopt.bc.push_back(BUSOptionSubstr(0, 6, 12));
    strand_out = ProgramOptions::StrandType::FR;
  } else if (opt.technology == "SPLIT-SEQ") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(0, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(1, 0, 10));
    busopt.bc.push_back(BUSOptionSubstr(1, 10, 18));
    busopt.bc.push_back(BUSOptionSubstr(1, 48, 56));
    busopt.bc.push_back(BUSOptionSubstr(1, 78, 86));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "STORM-SEQ") {
    busopt.nfiles = 2;
    busopt.bc.push_back(BUSOptionSubstr(-1, -1, -1));
    busopt.umi.push_back(BUSOptionSubstr(1, 0, 8));
    busopt.seq.push_back(BUSOptionSubstr(0, 0, 0));
    busopt.seq.push_back(BUSOptionSubstr(1, 14, 0));
    busopt.paired = true;
    strand_out = ProgramOptions::StrandType::RF;
    opt.strand_specific = true;
  } else if (opt.technology == "SCRBSEQ") {
    busopt.nfiles = 2;
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(0, 6, 16));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 6));
  } else if (opt.technology == "SMARTSEQ3") {
    busopt.nfiles = 4;
    busopt.seq.push_back(BUSOptionSubstr(2, 22, 0));
    busopt.seq.push_back(BUSOptionSubstr(3, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(2, 0, 19));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 0));
    busopt.bc.push_back(BUSOptionSubstr(1, 0, 0));
    busopt.paired = true;
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "SMARTSEQ2") {
    busopt.nfiles = 3;
    busopt.seq.push_back(BUSOptionSubstr(2, 0, 0));
    busopt.umi.push_back(BUSOptionSubstr(-1, -1, -1));
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 0));
    busopt.bc.push_back(BUSOptionSubstr(1, 0, 0));
    if (!opt.single_end) {
      busopt.nfiles++;
      busopt.seq.push_back(BUSOptionSubstr(3, 0, 0));
      if (!opt.long_read) {
        busopt.paired = true;
      }
    }
  } else if (opt.technology == "BDWTA") {
    busopt.nfiles = 2;
    busopt.bc.push_back(BUSOptionSubstr(0, 0, 9));
    busopt.bc.push_back(BUSOptionSubstr(0, 9 + 12, 9 + 12 + 9));
    busopt.bc.push_back(BUSOptionSubstr(0, 9 + 12 + 9 + 13, 9 + 12 + 9 + 13 + 9));
    busopt.umi.push_back(
        BUSOptionSubstr(0, 9 + 12 + 9 + 13 + 9, 9 + 12 + 9 + 13 + 9 + 8));
    busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else if (opt.technology == "VASA-SEQ") {
    busopt.nfiles = 1;
    busopt.bc.push_back(BUSOptionSubstr(0, 6, 14));
    busopt.umi.push_back(BUSOptionSubstr(0, 0, 6));
    busopt.seq.push_back(BUSOptionSubstr(0, 14, 0));
    strand_out = ProgramOptions::StrandType::FR;
    opt.strand_specific = true;
  } else {
    bool valid = ParseTechnology(opt.technology, busopt, errorList);

    if (busopt.seq.size() == 2 && !opt.single_end && !opt.long_read) {
      busopt.paired = true;
    }

    if (!valid) {
      errorList.push_back("Unable to create technology: " + opt.technology);
      return false;
    }
  }

  return true;
}

}  // namespace bus_parser
