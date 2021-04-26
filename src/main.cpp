#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
#include <thread>
#include <time.h>
#include <algorithm>

#include <cstdio>

#include <zlib.h>

#include "common.h"
#include "ProcessReads.h"
#include "KmerIndex.h"
#include "Kmer.hpp"
#include "MinCollector.h"
#include "EMAlgorithm.h"
#include "weights.h"
#include "Inspect.h"
#include "Bootstrap.h"
#include "H5Writer.h"
#include "PlaintextWriter.h"
#include "GeneModel.h"
#include "Merge.h"


//#define ERROR_STR "\033[1mError:\033[0m"
#define ERROR_STR "Error:"

using namespace std;


int my_mkdir(const char *path, mode_t mode) {
  #ifdef _WIN64
  return mkdir(path);
  #else
  return mkdir(path,mode);
  #endif
}

bool checkFileExists(std::string fn) {
  struct stat stFileInfo;
  auto intStat = stat(fn.c_str(), &stFileInfo);
  return intStat == 0;
}


void ParseOptionsIndex(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int make_unique_flag = 0;
  const char *opt_string = "i:k:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"make-unique", no_argument, &make_unique_flag, 1},
    // short args
    {"index", required_argument, 0, 'i'},
    {"kmer-size", required_argument, 0, 'k'},
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
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.k;
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
  if (make_unique_flag) {
    opt.make_unique = true;
  }

  for (int i = optind; i < argc; i++) {
    opt.transfasta.push_back(argv[i]);
  }
}

void ParseOptionsInspect(int argc, char **argv, ProgramOptions& opt) {


  const char *opt_string = "G:g:b:";

  int para_flag = 0;
  static struct option long_options[] = {
    // long args
    {"gfa", required_argument, 0, 'G'},
    {"gtf", required_argument, 0, 'g'},
    {"bed", required_argument, 0, 'b'},
    {"paranoid", no_argument, &para_flag, 1},
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
    case 'G': {
      opt.gfa = optarg;
      break;
    }
    case 'b': {
      stringstream(optarg) >> opt.bedFile;
      break;
    }
    case 'g': {
      stringstream(optarg) >> opt.gtfFile;
      break;
    }
    default: break;
    }
  }
  opt.index = argv[optind];

  if (para_flag) {
    opt.inspect_thorough = true;
  }
}



void ParseOptionsEM(int argc, char **argv, ProgramOptions& opt) {
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

  const char *opt_string = "t:i:l:s:o:n:m:d:b:g:c:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
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
    {"min-range", required_argument, 0, 'm'},
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
    case 'm': {
      stringstream(optarg) >> opt.min_range;
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
}

void ParseOptionsEMOnly(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int plaintext_flag = 0;

  const char *opt_string = "t:s:l:s:o:n:m:d:b:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"plaintext", no_argument, &plaintext_flag, 1},
    {"seed", required_argument, 0, 'd'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {"min-range", required_argument, 0, 'm'},
    {"bootstrap-samples", required_argument, 0, 'b'},
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
    case 'm': {
      stringstream(optarg) >> opt.min_range;
    }
    case 'b': {
      stringstream(optarg) >> opt.bootstrap;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.seed;
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }

  if (plaintext_flag) {
    opt.plaintext = true;
  }
}

void ParseOptionsTCCQuant(int argc, char **argv, ProgramOptions& opt) {
  const char *opt_string = "o:i:T:e:f:l:s:t:g:G:b:d:";
  static struct option long_options[] = {
    {"index", required_argument, 0, 'i'},
    {"txnames", required_argument, 0, 'T'},
    {"threads", required_argument, 0, 't'},
    {"fragment-file", required_argument, 0, 'f'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {"ec-file", required_argument, 0, 'e'},
    {"genemap", required_argument, 0, 'g'},
    {"gtf", required_argument, 0, 'G'},
    {"bootstrap-samples", required_argument, 0, 'b'},
    {"seed", required_argument, 0, 'd'},
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
    case 'f': {
      stringstream(optarg) >> opt.fldFile;
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
      stringstream(optarg) >> opt.bootstrap;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.seed;
      break;
    }
    case 'T': {
      opt.transcriptsFile = optarg;
      break;
    }
    default: break;
    }
  }
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
  if (opt.files.size() > 0) {
    opt.tccFile = opt.files[0];
  }
}

void ParseOptionsPseudo(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int single_flag = 0;
  int strand_flag = 0;
  int pbam_flag = 0;
  int gbam_flag = 0;
  int umi_flag = 0;
  int quant_flag = 0;
  int bus_flag = 0;

  const char *opt_string = "t:i:l:s:o:b:u:g:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"single", no_argument, &single_flag, 1},
    //{"strand-specific", no_argument, &strand_flag, 1},
    {"pseudobam", no_argument, &pbam_flag, 1},
    {"quant", no_argument, &quant_flag, 1},
    {"bus", no_argument, &bus_flag, 1},
    {"umi", no_argument, &umi_flag, 'u'},
    {"batch", required_argument, 0, 'b'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"gtf", required_argument, 0, 'g'},
    {"index", required_argument, 0, 'i'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
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
    case 'g': {
      stringstream(optarg) >> opt.gtfFile;
      break;
    }
    case 'b': {
      opt.batch_mode = true;
      opt.batch_file_name = optarg;
      break;
    }
    default: break;
    }
  }
  
  if (umi_flag) {
    opt.umi = true;
    opt.single_end = true; // UMI implies single-end reads
  }
  
  if (bus_flag) {
    opt.batch_bus = true;
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
 
  if (verbose_flag) {
    opt.verbose = true;
  }

  if (single_flag) {
    opt.single_end = true;
    opt.single_overhang = true;
  }
  
  if (quant_flag) {
    opt.pseudo_quant = true;
  }

  if (strand_flag) {
    opt.strand_specific = true;
  }

  if (pbam_flag) {
    opt.pseudobam = true;
  }
}


void ParseOptionsMerge(int argc, char **argv, ProgramOptions& opt) {

  const char *opt_string = "i:o:";
  static struct option long_options[] = {
    {"index", required_argument, 0, 'i'},
    {"output-dir", required_argument, 0, 'o'},
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
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    default: break;
    }
  }
  
  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
 
}

void ListSingleCellTechnologies() {
  //todo, figure this out
  cout << "List of supported single-cell technologies" << endl << endl 
  << "short name       description" << endl
  << "----------       -----------" << endl
  << "10xv1            10x version 1 chemistry" << endl
  << "10xv2            10x version 2 chemistry" << endl
  << "10xv3            10x version 3 chemistry" << endl
  << "BDWTA            BD Rhapsody WTA" << endl
  << "CELSeq           CEL-Seq" << endl
  << "CELSeq2          CEL-Seq version 2" << endl
  << "DropSeq          DropSeq" << endl
  << "inDropsv1        inDrops version 1 chemistry" << endl
  << "inDropsv2        inDrops version 2 chemistry" << endl
  << "inDropsv3        inDrops version 3 chemistry" << endl
  << "SCRBSeq          SCRB-Seq" << endl
  << "SmartSeq3        Smart-seq3" << endl
  << "SureCell         SureCell for ddSEQ" << endl
  << endl;
 }

void ParseOptionsBus(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int gbam_flag = 0;
  int paired_end_flag = 0;
  int strand_FR_flag = 0;
  int strand_RF_flag = 0;

  const char *opt_string = "i:o:x:t:lbng:c:T:";
  static struct option long_options[] = {
    {"verbose", no_argument, &verbose_flag, 1},
    {"index", required_argument, 0, 'i'},
    {"output-dir", required_argument, 0, 'o'},
    {"technology", required_argument, 0, 'x'},
    {"list", no_argument, 0, 'l'},
    {"threads", required_argument, 0, 't'},
    {"bam", no_argument, 0, 'b'},
    {"num", no_argument, 0, 'n'},
    {"genomebam", no_argument, &gbam_flag, 1},
    {"gtf", required_argument, 0, 'g'},
    {"chromosomes", required_argument, 0, 'c'},
    {"tag", required_argument, 0, 'T'},
    {"fr-stranded", no_argument, &strand_FR_flag, 1},
    {"rf-stranded", no_argument, &strand_RF_flag, 1},
    {"paired", no_argument, &paired_end_flag, 1},
    {0,0,0,0}
  };

  int list_flag = 0;
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
      std::transform(opt.technology.begin(), opt.technology.end(),opt.technology.begin(), ::toupper);
      break;
    }
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'b': {
      opt.bam = true;
      break;
    }
    case 'n': {
      opt.num = true;
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
    case 'T': {
      stringstream(optarg) >> opt.tagsequence;
      break;
    }
    default: break;
    }
  }

  if (list_flag) {
    ListSingleCellTechnologies();
    exit(1);
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

  if (paired_end_flag) {
    opt.single_end = false;
  } else {
    opt.single_end = true;
  }
  
  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
 
}

bool ParseTechnology(const std::string &techstr, BUSOptions& busopt, std::vector<std::string> &errorList) {
  auto i1 = techstr.find(':');
  if (i1 == std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), none found: \"" + techstr + "\"");    
    return false;
  }
  auto i2 = techstr.find(':', i1+1);
  if (i2 == std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), only one found: \"" + techstr + "\"");    
    return false;
  }
  auto ip = techstr.find(':', i2+1);
  if (ip != std::string::npos) {
    errorList.push_back("Error: technology string must contain two colons (:), three found: \"" + techstr + "\"");    
    return false;
  }
  auto bcstr = techstr.substr(0, i1);
  auto umistr = techstr.substr(i1+1,i2-i1-1);
  auto seqstr = techstr.substr(i2+1);
  


  int maxnf = 0;

  auto convert_commas_to_vector = [&](const std::string &s, std::vector<BUSOptionSubstr> &v) -> bool {
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
      for (int i = 0; i+2 < nv; i+=3) {
        int f = vv[i];
        int a = vv[i+1];
        int b = vv[i+2];
        if (f < 0) {
          errorList.push_back("Error: invalid file number (" + to_string(f) + ")  " + s);
        }
        if (a <  0) {
          errorList.push_back("Error: invalid start (" + to_string(a) + ")  " + s);
        }
        if (b != 0 && b <= a) {
          errorList.push_back("Error: invalid stop (" + to_string(b) + ") has to be after start (" + to_string(a) + ")  " + s);
        }
        v.push_back(BUSOptionSubstr(f,a,b));
        if (f > maxnf) {
          maxnf = f;
        }
      }
    } else {
      errorList.push_back("Error: number of values has to be multiple of 3 " + s);
      return false;
    }

    busopt.nfiles = maxnf+1;
    return true;
  };

  

  std::vector<BUSOptionSubstr> v;
  if (!convert_commas_to_vector(bcstr,v)) {
    return false;
  }
  if (v.empty()) {
    errorList.push_back("Error: empty barcode list " + bcstr);
    return false;
  }
  busopt.bc = std::move(v);

  if (!convert_commas_to_vector(umistr, v)) {
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



void ParseOptionsH5Dump(int argc, char **argv, ProgramOptions& opt) {
  int peek_flag = 0;
  const char *opt_string = "o:";
  static struct option long_options[] = {
    // long args
    {"peek", no_argument, &peek_flag, 1},
    // short args
    {"output-dir", required_argument, 0, 'o'},
    {0,0,0,0}
  };

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
    default: break;
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

bool CheckOptionsBus(ProgramOptions& opt) {
  bool ret = true;

  cerr << endl;

  // check index
  if (opt.index.empty()) {
    cerr << ERROR_STR << " kallisto index file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
      ret = false;
    }
  }

  // check files
  if (opt.files.size() == 0) {
    cerr << ERROR_STR << " Missing read files" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    for (auto& fn : opt.files) {
      auto intStat = stat(fn.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << ERROR_STR << " file not found " << fn << endl;
        ret = false;
      }
    }
  }

  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      } 
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }    
  }

  if (opt.technology.empty()) {
    cerr << "Error: need to specify technology to use" << endl;
    ret = false;
  } else {
    auto& busopt = opt.busOptions;
   
    if (opt.bam) { // Note: only 10xV2 has been tested
      busopt.nfiles = 1;
      busopt.paired = false;
      if (opt.technology == "10XV2") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0)); // second file, entire string
        busopt.umi.push_back(BUSOptionSubstr(0,16,26)); // first file [16:26]
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
      } else if (opt.technology == "10XV3") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,16,28));
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
//      } else if (opt.technology == "10XV1") {

      } else if (opt.technology == "SURECELL") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,18,26));
        busopt.bc.push_back(BUSOptionSubstr(0,0,18));
      } else if (opt.technology == "DROPSEQ") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,12,20));
        busopt.bc.push_back(BUSOptionSubstr(0,0,12));
      } else if (opt.technology == "INDROPSV1") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,42,48));
        busopt.bc.push_back(BUSOptionSubstr(0,0,11));
        busopt.bc.push_back(BUSOptionSubstr(0,30,38));  
      } else if (opt.technology == "INDROPSV2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(0,0,0));
        busopt.umi.push_back(BUSOptionSubstr(1,42,48));
        busopt.bc.push_back(BUSOptionSubstr(1,0,11));
        busopt.bc.push_back(BUSOptionSubstr(1,30,38));  
      } else if (opt.technology == "INDROPSV3") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi.push_back(BUSOptionSubstr(1,8,14));
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
        busopt.bc.push_back(BUSOptionSubstr(1,0,8));
      } else if (opt.technology == "CELSEQ") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,8,12));
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
      } else if (opt.technology == "CELSEQ2") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,0,6));
        busopt.bc.push_back(BUSOptionSubstr(0,6,12));
      } else if (opt.technology == "SCRBSEQ") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,6,16));
        busopt.bc.push_back(BUSOptionSubstr(0,0,6));
      } else if (opt.technology == "INDROPSV3") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,0,6));
        busopt.bc.push_back(BUSOptionSubstr(0,6,16));
      } else if (opt.technology == "SMARTSEQ3") {
        busopt.nfiles = 4;
        busopt.seq.push_back(BUSOptionSubstr(2,22,0));
        busopt.seq.push_back(BUSOptionSubstr(3,0,0));
        busopt.umi.push_back(BUSOptionSubstr(2,0,19));
        busopt.bc.push_back(BUSOptionSubstr(0,0,0));
        busopt.bc.push_back(BUSOptionSubstr(1,0,0));
        busopt.paired = true;
      } else if (opt.technology == "BDWTA") {
        busopt.nfiles = 2;
        busopt.bc.push_back(BUSOptionSubstr(0, 0, 9)); // bc1 CLS1
        // busopt.bc.push_back(BUSOptionSubstr(0,9,9+12)); // linker
        busopt.bc.push_back(BUSOptionSubstr(0, 9 + 12, 9 + 12 + 9)); // bc2 CLS2
        // busopt.bc.push_back(BUSOptionSubstr(0,9+12+9,9+12+9+13)); // linker
        busopt.bc.push_back(BUSOptionSubstr(0, 9 + 12 + 9 + 13, 9 + 12 + 9 + 13 + 9));          // bc3 CLS3
        busopt.umi.push_back(BUSOptionSubstr(0, 9 + 12 + 9 + 13 + 9, 9 + 12 + 9 + 13 + 9 + 8)); // umi
        busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
      } else {
        vector<int> files;
        vector<BUSOptionSubstr> values;
        vector<BUSOptionSubstr> bcValues;
        vector<std::string> errorList;
        //bool invalid = ParseTechnology(opt.technology, values, files, errorList, bcValues);
        bool valid = ParseTechnology(opt.technology, busopt, errorList);
        
        if (busopt.seq.size() == 2 && !opt.single_end) {
          busopt.paired = true;
        }
        
        if(!valid) {
          /*
          busopt.nfiles = files.size(); 
          for(int i = 0; i < bcValues.size(); i++) {
            busopt.bc.push_back(bcValues[i]);
          }
          busopt.umi.push_back(values[0]);
          busopt.seq.push_back(values[1]);
          */
        } else {
          for(int j = 0; j < errorList.size(); j++) {
            cerr << errorList[j] << endl;
          }
          cerr << "Unable to create technology: " << opt.technology << endl;
          ret = false;
        }
      }
    } else {
      busopt.paired = false;
      if (opt.technology == "10XV2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0)); // second file, entire string
        busopt.umi.push_back(BUSOptionSubstr(0,16,26)); // first file [16:26]
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
      } else if (opt.technology == "10XV3") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,16,28));
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
      } else if (opt.technology == "10XV1") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi.push_back(BUSOptionSubstr(1,0,10));
        busopt.bc.push_back(BUSOptionSubstr(0,0,14));
      } else if (opt.technology == "SURECELL") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,51,59));
        busopt.bc.push_back(BUSOptionSubstr(0,0,6));
        busopt.bc.push_back(BUSOptionSubstr(0,21,27));
        busopt.bc.push_back(BUSOptionSubstr(0,42,48));
      } else if (opt.technology == "DROPSEQ") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,12,20));
        busopt.bc.push_back(BUSOptionSubstr(0,0,12));
      } else if (opt.technology == "INDROPSV1") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,42,48));
        busopt.bc.push_back(BUSOptionSubstr(0,0,11));
        busopt.bc.push_back(BUSOptionSubstr(0,30,38));  
      } else if (opt.technology == "INDROPSV2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(0,0,0));
        busopt.umi.push_back(BUSOptionSubstr(1,42,48));
        busopt.bc.push_back(BUSOptionSubstr(1,0,11));
        busopt.bc.push_back(BUSOptionSubstr(1,30,38));  
      } else if (opt.technology == "INDROPSV3") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi.push_back(BUSOptionSubstr(1,8,14));
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
        busopt.bc.push_back(BUSOptionSubstr(1,0,8));
      } else if (opt.technology == "CELSEQ") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,8,12));
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
      } else if (opt.technology == "CELSEQ2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,0,6));
        busopt.bc.push_back(BUSOptionSubstr(0,6,12));
      } else if (opt.technology == "SCRBSEQ") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi.push_back(BUSOptionSubstr(0,6,16));
        busopt.bc.push_back(BUSOptionSubstr(0,0,6));
      } else if (opt.technology == "INDROPSV3") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi.push_back(BUSOptionSubstr(1,8,14));
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
        busopt.bc.push_back(BUSOptionSubstr(1,0,8));
      } else if (opt.technology == "SMARTSEQ3") {
        busopt.nfiles = 4;
        busopt.seq.push_back(BUSOptionSubstr(2,22,0));
        busopt.seq.push_back(BUSOptionSubstr(3,0,0));
        busopt.umi.push_back(BUSOptionSubstr(2,0,19));
        busopt.bc.push_back(BUSOptionSubstr(0,0,0));
        busopt.bc.push_back(BUSOptionSubstr(1,0,0));
        busopt.paired = true;
      } else if (opt.technology == "BDWTA") {
        busopt.nfiles = 2;
        busopt.bc.push_back(BUSOptionSubstr(0, 0, 9)); // bc1 CLS1
        // busopt.bc.push_back(BUSOptionSubstr(0,9,9+12)); // linker
        busopt.bc.push_back(BUSOptionSubstr(0, 9 + 12, 9 + 12 + 9)); // bc2 CLS2
        // busopt.bc.push_back(BUSOptionSubstr(0,9+12+9,9+12+9+13)); // linker
        busopt.bc.push_back(BUSOptionSubstr(0, 9 + 12 + 9 + 13, 9 + 12 + 9 + 13 + 9));          // bc3 CLS3
        busopt.umi.push_back(BUSOptionSubstr(0, 9 + 12 + 9 + 13 + 9, 9 + 12 + 9 + 13 + 9 + 8)); // umi
        busopt.seq.push_back(BUSOptionSubstr(1, 0, 0));
      } else {
        vector<int> files;
        vector<BUSOptionSubstr> values;
        vector<BUSOptionSubstr> bcValues;
        vector<std::string> errorList;        
        //bool invalid = ParseTechnology(opt.technology, values, files, errorList, bcValues);
        bool valid = ParseTechnology(opt.technology, busopt, errorList);
        
        if (busopt.seq.size() == 2 && !opt.single_end) {
          busopt.paired = true;
        }
        
        if(valid) {
          /*
          busopt.nfiles = files.size(); 
          for(int i = 0; i < bcValues.size(); i++) {
            busopt.bc.push_back(bcValues[i]);
          }
          busopt.umi.push_back(values[0]);
          busopt.seq.push_back(values[1]);
          */
        } else {
          for(int j = 0; j < errorList.size(); j++) {
            cerr << errorList[j] << endl;
          }
          cerr << "Unable to create technology: " << opt.technology << endl;
          ret = false;
        }
      }
    }
  }
  
  if (opt.tagsequence.empty() && (opt.technology == "SMARTSEQ3")) {
    opt.tagsequence = "ATTGCGCAATG";
    cerr << "[bus] Using " <<  opt.tagsequence << " as UMI tag sequence" << endl;
  }
  
  if (!opt.tagsequence.empty()) {
    opt.busOptions.umi[0].start += opt.tagsequence.length();
    if (opt.busOptions.umi[0].start >= opt.busOptions.umi[0].stop) {
      cerr << "Error: Tag sequence must be shorter than UMI sequence" << endl;
      ret = false;
    }
  }
  
  if (!opt.single_end && !opt.busOptions.paired) {
      cerr << "Error: Paired reads are not compatible with the specified technology" << endl;
      ret = false;
  }

  if (opt.strand_specific) {
    if (opt.busOptions.seq.size() != 1 && !(opt.busOptions.seq.size() == 2 && opt.busOptions.paired)) {
      cerr << "Error: Strand-specific read processing is only supported for technologies with a single CDNA read file or paired-end reads" << endl;
      ret = false;
    }
  }

  if (opt.genomebam) {
    if (opt.busOptions.seq.size() != 1 && !(opt.busOptions.seq.size() == 2 && opt.busOptions.paired)) {
      cerr << "Error: BAM output is currently only supported for technologies with a single CDNA read file or paired-end reads" << endl;
      ret = false;
    }
    if (!opt.gtfFile.empty()) {
      if (!checkFileExists(opt.gtfFile)) {
        cerr << "Error: GTF file " << opt.gtfFile << " does not exist" << endl;
        ret = false;
      }
    } else {
      cerr << "Error: need GTF file for genome alignment" << endl;
      ret = false;
    }
    if (!opt.chromFile.empty()) {
      if (!checkFileExists(opt.chromFile)) {
        cerr << "Error: Chromosome file not found: " << opt.chromFile << endl;
        ret = false;
      }
    }
  }

  

  if (ret && !opt.bam && opt.files.size() %  opt.busOptions.nfiles != 0) {
    cerr << "Error: Number of files (" << opt.files.size() << ") does not match number of input files required by "
    << "technology " << opt.technology << " (" << opt.busOptions.nfiles << ")" << endl;
    ret = false;
  }

  if (opt.bam && opt.num) {
    cerr << "Warning: --bam option was used, so --num option will be ignored" << endl;
  }

  return ret;
}

bool CheckOptionsIndex(ProgramOptions& opt) {

  bool ret = true;

  if (opt.k <= 1 || opt.k >= Kmer::MAX_K) {
    cerr << "Error: invalid k-mer length " << opt.k << ", minimum is 3 and  maximum is " << (Kmer::MAX_K -1) << endl;
    ret = false;
  }

  if (opt.k % 2 == 0) {
    cerr << "Error: k needs to be an odd number" << endl;
    ret = false;
  }

  if (opt.transfasta.empty()) {
    cerr << "Error: no FASTA files specified" << endl;
    ret = false;
  } else {

    for (auto& fasta : opt.transfasta) {
      // we want to generate the index, check k, index and transfasta
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: FASTA file not found " << fasta << endl;
        ret = false;
      }
    }
  }

  if (opt.index.empty()) {
    cerr << "Error: need to specify kallisto index name" << endl;
    ret = false;
  }

  return ret;
}

bool CheckOptionsEM(ProgramOptions& opt, bool emonly = false) {

  bool ret = true;

  cerr << endl;
  // check for index
  if (!emonly) {
    if (opt.index.empty()) {
      cerr << ERROR_STR << " kallisto index file missing" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      auto intStat = stat(opt.index.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
        ret = false;
      }
    }
  }

  // check for read files
  if (!emonly) {
    if (opt.files.size() == 0) {
      cerr << ERROR_STR << " Missing read files" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      for (auto& fn : opt.files) {
        auto intStat = stat(fn.c_str(), &stFileInfo);
        if (intStat != 0) {
          cerr << ERROR_STR << " file not found " << fn << endl;
          ret = false;
        }
      }
    }

    /*
    if (opt.strand_specific && !opt.single_end) {
      cerr << "Error: strand-specific mode requires single-end mode" << endl;
      ret = false;
    }*/

    if (!opt.single_end) {
      if (opt.files.size() % 2 != 0) {
        cerr << "Error: paired-end mode requires an even number of input files" << endl
            << "       (use --single for processing single-end reads)" << endl;
        ret = false;
      }
    }    
  }

  if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
    cerr << "Error: cannot supply mean/sd without supplying both -l and -s" << endl;
    ret = false;
  }

  if (opt.single_end && (opt.fld == 0.0 || opt.sd == 0.0)) {
    cerr << "Error: fragment length mean and sd must be supplied for single-end reads using -l and -s" << endl;
    ret = false;
  } else if (opt.fld == 0.0 && ret) {
    // In the future, if we have single-end data we should require this
    // argument
    cerr << "[quant] fragment length distribution will be estimated from the data" << endl;
  } else if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
    cerr << "[quant] fragment length distribution is truncated gaussian with mean = " <<
      opt.fld << ", sd = " << opt.sd << endl;
  }

  if (!opt.single_end && (opt.fld > 0.0 && opt.sd > 0.0)) {
    cerr << "[~warn] you specified using a gaussian but have paired end data" << endl;
    cerr << "[~warn] we suggest omitting these parameters and let us estimate the distribution from data" << endl;
  }

  if (opt.fld < 0.0) {
    cerr << "Error: invalid value for mean fragment length " << opt.fld << endl;
    ret = false;
  }

  if (opt.sd < 0.0) {
    cerr << "Error: invalid value for fragment length standard deviation " << opt.sd << endl;
    ret = false;
  }

  if (opt.iterations <= 0) {
    cerr << "Error: invalid number of iterations " << opt.iterations << endl;
    ret = false;
  }

  if (opt.min_range <= 0) {
    cerr << "Error: invalid value for minimum range " << opt.min_range << endl;
    ret = false;
  }
  
  if (opt.genomebam) {
    if (!opt.gtfFile.empty()) {
      struct stat stFileInfo;
      auto intStat = stat(opt.gtfFile.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: GTF file " << opt.gtfFile << " does not exist" << endl;
        ret = false;
      }
    } else {
      cerr << "Error: need GTF file for genome alignment" << endl;
      ret = false;
    }

    if (!opt.chromFile.empty()) {
      struct stat stFileInfo;
      auto intStat = stat(opt.chromFile.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: Chromosome file not found: " << opt.chromFile << endl;
        ret = false;
      }
    }
    
  }


  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      } else if (emonly) {
        // check for directory/counts.txt
        struct stat stCountInfo;
        auto intcountstat = stat((opt.output + "/counts.txt" ).c_str(), &stCountInfo);
        if (intcountstat != 0) {
          cerr << "Error: could not find file " << opt.output << "/counts.txt" << endl;
          ret = false;
        }

        // check for directory/index.saved
        struct stat stIndexInfo;
        auto intindexstat = stat((opt.output + "/index.saved").c_str(), &stIndexInfo);
        if (intindexstat != 0) {
          cerr << "Error: could not find index " << opt.output << "/index.saved" << endl;
          ret = false;
        }
        opt.index = (opt.output + "/index.saved");
      }
    } else {
      if (emonly) {
        cerr << "Error: output directory needs to exist, run quant first" << endl;
        ret = false;
      } else {
        // create directory
        if (my_mkdir(opt.output.c_str(), 0777) == -1) {
          cerr << "Error: could not create directory " << opt.output << endl;
          ret = false;
        }
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }
    if (opt.threads > 1 && opt.pseudobam) {
      //cerr << "Error: pseudobam is not compatible with running on many threads."<< endl;
      //ret = false;
    }
  }

  if (opt.bootstrap < 0) {
    cerr << "Error: number of bootstrap samples must be a non-negative integer." << endl;
    ret = false;
  }

  #ifndef USE_HDF5  
  if (opt.bootstrap > 0 && !opt.plaintext) {
    cerr << "Warning: kallisto was not compiled with HDF5 support so no bootstrapping" << endl
         << "will be performed. Run quant with --plaintext option or recompile with" << endl
         << "HDF5 support to obtain bootstrap estimates." << endl;
    opt.bootstrap = 0;
  }
  #endif
  return ret;
}

bool CheckOptionsTCCQuant(ProgramOptions& opt) {
  
  bool ret = true;
  
  cerr << endl;
  // check for index
  if (opt.index.empty() && opt.transcriptsFile.empty()) {
    cerr << ERROR_STR << " either a kallisto index file or a transcripts file need to be supplied" << endl;
    ret = false;
  } else if (!opt.index.empty() && !opt.transcriptsFile.empty()) {
    cerr << ERROR_STR << " cannot supply both a kallisto index file and a transcripts file" << endl;
    ret = false;
  } else if (!opt.index.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
      ret = false;
    }
  } else if (!opt.transcriptsFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.transcriptsFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " transcripts file not found " << opt.transcriptsFile << endl;
      ret = false;
    }
  }

  if (opt.tccFile.empty()) {
    cerr << ERROR_STR << " transcript-compatibility counts file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.tccFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " transcript-compatibility counts file not found " << opt.tccFile << endl;
      ret = false;
    }
  }
  
  if (opt.ecFile.empty() && !opt.transcriptsFile.empty()) {
    cerr << ERROR_STR << " equivalence class file must be supplied if transcripts file is supplied " << opt.tccFile << endl;
    ret = false;
  }

  if (!opt.ecFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.ecFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " equivalence class file not found " << opt.ecFile << endl;
      ret = false;
    }
  }
  
  if (!opt.fldFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.fldFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " fragment length distribution file not found " << opt.fldFile << endl;
      ret = false;
    }
  }
  
  if (!opt.genemap.empty() && !opt.gtfFile.empty()) {
    cerr << ERROR_STR << " Cannot supply both --genemap and --gtf" << endl;
    ret = false;
  }
  
  if (!opt.gtfFile.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.gtfFile.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " GTF file not found " << opt.gtfFile << endl;
      ret = false;
    }
  }
  
  if (!opt.genemap.empty()) {
    struct stat stFileInfo;
    auto intStat = stat(opt.genemap.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " file for mapping transcripts to genes not found " << opt.genemap << endl;
      ret = false;
    }
  }
  
  if ((opt.fld != 0.0 || opt.sd != 0.0) && !opt.fldFile.empty()) {
    cerr << ERROR_STR <<" cannot supply mean or sd while also supplying a fragment length distribution file" << endl;
    ret = false;
  }
  
  if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
    cerr << ERROR_STR << " cannot supply mean/sd without supplying both -l and -s" << endl;
    ret = false;
  }
  
  if (opt.index.empty() && (!opt.fldFile.empty() || opt.fld != 0.0 || opt.sd != 0.0)) {
    cerr << ERROR_STR << " cannot supply fragment length information unless a kallisto index is provided" << endl;
    ret = false;
  }
  
  if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
    cerr << "[tcc] fragment length distribution is truncated gaussian with mean = " <<
      opt.fld << ", sd = " << opt.sd << endl;
  }
  
  if (opt.fld < 0.0) {
    cerr << ERROR_STR << " invalid value for mean fragment length " << opt.fld << endl;
    ret = false;
  }
  
  if (opt.sd < 0.0) {
    cerr << ERROR_STR << " invalid value for fragment length standard deviation " << opt.sd << endl;
    ret = false;
  }

  if (opt.output.empty()) {
    cerr << ERROR_STR << " need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << ERROR_STR << " file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      }
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << ERROR_STR << " could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << ERROR_STR << " invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }
  }
  
  if (opt.bootstrap < 0) {
    cerr << "Error: number of bootstrap samples must be a non-negative integer." << endl;
    ret = false;
  }

   return ret;
}

bool CheckOptionsMerge(ProgramOptions& opt) {

  bool ret = true;

  cerr << endl;
  // check for index
  if (opt.index.empty()) {
    cerr << ERROR_STR << " kallisto index file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
      ret = false;
    }
  }

    
  if (opt.files.size() == 0) {
    cerr << ERROR_STR << " Missing input directory to merge" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;      
    for (auto& fn : opt.files) {        
      auto intStat = stat(fn.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << ERROR_STR << " input directory not found " << fn << endl;
        ret = false;
      } else {
        if (!S_ISDIR(stFileInfo.st_mode)) {
          cerr << "Error: file " << fn << " exists but is not a directory" << endl;
          ret = false;
        }

        
        if (!checkFileExists(fn + "/matrix.ec")) {
          cerr << "Error: file " << fn << "/matrix.ec was not found, check that it was run in batch mode" << endl;
          ret = false;
        }
        if (!checkFileExists(fn + "/matrix.cells")) {
          cerr << "Error: file " << fn << "/matrix.cells was not found, check that it was run in batch mode" << endl;
          ret = false;
        }
        if (!checkFileExists(fn + "/matrix.tcc.mtx")) {
          cerr << "Error: file " << fn << "/matrix.tcc.mtx was not found, check that it was run in batch mode" << endl;
          ret = false;
        }
      }
    }
  }


  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      }

      auto it = std::find(opt.files.begin(), opt.files.end(), opt.output);
      if (it != opt.files.end()) {
        cerr << "Error: output directory cannot be part of input directory " << opt.output << endl;
        ret = false;
      }
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  }

  
  return ret;
}

bool CheckOptionsPseudo(ProgramOptions& opt) {

  bool ret = true;

  cerr << endl;
  // check for index
  if (opt.index.empty()) {
    cerr << ERROR_STR << " kallisto index file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
      ret = false;
    }
  }

  if (opt.pseudo_quant) {
    if (!opt.batch_mode) {
      cerr << ERROR_STR << " --quant can only be run with in batch mode with --batch" << endl;
      ret = false;
    }
    if (opt.batch_bus) {
      cerr << ERROR_STR << " --quant cannot be run with --bus" << endl;
      ret = false;
    }
  }
  
  if (opt.batch_bus && opt.umi) {
    cerr << ERROR_STR << " UMI cannot be run with --bus" << endl;
    ret = false;   
  }

  // check for read files
  if (!opt.batch_mode) {
    if (opt.umi) {
      cerr << ERROR_STR << " UMI must be run in batch mode, use --batch option" << endl;
      ret = false;      
    }
    if (opt.batch_bus) {
      cerr << ERROR_STR << " --bus must be run in batch mode, use --batch option" << endl;
      ret = false;      
    }
    
    if (opt.files.size() == 0) {
      cerr << ERROR_STR << " Missing read files" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;      
      for (auto& fn : opt.files) {        
        auto intStat = stat(fn.c_str(), &stFileInfo);
        if (intStat != 0) {
          cerr << ERROR_STR << " file not found " << fn << endl;
          ret = false;
        }
      }
    }
  } else {
    if (opt.files.size() != 0) {
      cerr << ERROR_STR << " cannot specify batch mode and supply read files" << endl;
      ret = false;
    } else {
      // check for batch files
      if (opt.batch_mode) {
        struct stat stFileInfo;
        auto intstat = stat(opt.batch_file_name.c_str(), &stFileInfo);
        if (intstat != 0) {
          cerr << ERROR_STR << " file not found " << opt.batch_file_name << endl;
          ret = false;
        }
        // open the file, parse and fill the batch_files values
        std::ifstream bfile(opt.batch_file_name);
        std::string line;
        std::string id,f1,f2;
        while (std::getline(bfile,line)) {
          if (line.size() == 0) {
            continue;
          }
          std::stringstream ss(line);
          ss >> id;
          if (id[0] == '#') {
            continue;
          }
          opt.batch_ids.push_back(id);
          if (opt.single_end && !opt.umi) {
            ss >> f1;
            opt.batch_files.push_back({f1});
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f1 << endl;
              ret = false;
            }
          } else {
            ss >> f1 >> f2;
            if (!opt.umi) {
              opt.batch_files.push_back({f1,f2});
            } else {
              opt.umi_files.push_back(f1);
              opt.batch_files.push_back({f2});
            }
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f1 << endl;
              ret = false;
            }
            intstat = stat(f2.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f2 << endl;
              ret = false;
            }
          }
        }
      }
    }
  }


  /*
  if (opt.strand_specific && !opt.single_end) {
    cerr << "Error: strand-specific mode requires single-end mode" << endl;
    ret = false;
  }*/

  if (!opt.single_end) {
    if (opt.files.size() % 2 != 0) {
      cerr << "Error: paired-end mode requires an even number of input files" << endl
           << "       (use --single for processing single-end reads)" << endl;
      ret = false;
    }
  }
  
  if (opt.umi) {
    opt.single_end = true;
    if (opt.fld != 0.0 || opt.sd != 0.0) {
      cerr << "[~warn] you supplied fragment length information for UMI data which will be ignored" << endl;
    }
  } else {
    if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
      cerr << "Error: cannot supply mean/sd without supplying both -l and -s" << endl;
      ret = false;
    }

    if (opt.single_end && (opt.fld == 0.0 || opt.sd == 0.0)) {
      cerr << "Error: fragment length mean and sd must be supplied for single-end reads using -l and -s" << endl;
      ret = false;
    } else if (opt.fld == 0.0 && ret) {
      // In the future, if we have single-end data we should require this
      // argument
      cerr << "[quant] fragment length distribution will be estimated from the data" << endl;
    } else if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
      cerr << "[quant] fragment length distribution is truncated gaussian with mean = " <<
        opt.fld << ", sd = " << opt.sd << endl;
    }

    if (!opt.single_end && (opt.fld > 0.0 && opt.sd > 0.0)) {
      cerr << "[~warn] you specified using a gaussian but have paired end data" << endl;
      cerr << "[~warn] we suggest omitting these parameters and let us estimate the distribution from data" << endl;
    }
  }

  if (opt.fld < 0.0) {
    cerr << "Error: invalid value for mean fragment length " << opt.fld << endl;
    ret = false;
  }

  if (opt.sd < 0.0) {
    cerr << "Error: invalid value for fragment length standard deviation " << opt.sd << endl;
    ret = false;
  }

  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      }
    } else {
      // create directory
      if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "[~warn]  you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }
    if (opt.threads > 1 && opt.pseudobam) {
      cerr << "Error: pseudobam is not compatible with running on many threads."<< endl;
      ret = false;
    }
  }

  if (opt.pseudobam) {
    cerr << "Error: Pseudobam not supported yet in pseudo mode, use quant mode to obtain BAM files" << endl;
    ret = false;
  }
  
  if (opt.batch_bus && ret) { // Set up BUS options
    auto& busopt = opt.busOptions;
    busopt.seq.push_back(BUSOptionSubstr(0,0,0));
    if (opt.single_end) {
      busopt.nfiles = 1;
    } else {
      busopt.nfiles = 2;
      busopt.seq.push_back(BUSOptionSubstr(1,0,0));
      busopt.paired = true;
    }
  }

  return ret;
}


bool CheckOptionsInspect(ProgramOptions& opt) {

  bool ret = true;
  // check for index
  if (opt.index.empty()) {
    cerr << "Error: kallisto index file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << "Error: kallisto index file not found " << opt.index << endl;
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
      cerr << "Error: GTF file not found " << opt.gtfFile << endl;
      ret = false;
    }

    if (!opt.chromFile.empty()) {
      struct stat stFileInfo;
      auto intStat = stat(opt.chromFile.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: Chromosome file not found " << opt.chromFile << endl;
        ret = false;
      }
    }

    if (opt.bedFile.empty()) {
      opt.bedFile = opt.index + ".bed";
    }
  }

  return ret;
}

bool CheckOptionsH5Dump(ProgramOptions& opt) {
  bool ret = true;
  if (!opt.peek) {
    if ( opt.output.size() == 0) {
      cerr << "Error: You must specify an output directory." << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      auto intStat = stat(opt.output.c_str(), &stFileInfo);
      if (intStat == 0) {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode)) {
          cerr << "Error: tried to open " << opt.output << " but another file already exists there" << endl;
          ret = false;
        }
      } else if (my_mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  } else {
    if (opt.output.size() > 0) {
      cerr << "Error: Cannot specify output directory and '--peek'. Please specify only one." << endl;
      ret = false;
    }
  }

  if (opt.files.size() == 0) {
    cerr << "Error: Missing H5 files" << endl;
    ret = false;
  } else if (opt.files.size() > 1) {
    cerr << "Error: Please specify only one H5 file" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    for (auto& fn : opt.files) {
      auto intStat = stat(fn.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: H5 file not found " << fn << endl;
        ret = false;
      }
    }
  }

  return ret;
}

void PrintCite() {
  cout << "When using this program in your research, please cite" << endl << endl
       << "  Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L." << endl
       << "  Near-optimal probabilistic RNA-seq quantification, "<< endl
       << "  Nature Biotechnology 34, 525-527(2016), doi:10.1038/nbt.3519" << endl
       << endl;
}

void PrintVersion() {
  cout << "kallisto, version " << 	KALLISTO_VERSION << endl;
}

void usage() {
  cout << "kallisto " << KALLISTO_VERSION << endl << endl
       << "Usage: kallisto <CMD> [arguments] .." << endl << endl
       << "Where <CMD> can be one of:" << endl << endl
       << "    index         Builds a kallisto index "<< endl
       << "    quant         Runs the quantification algorithm " << endl
       << "    quant-tcc     Runs quantification on transcript-compatibility counts" << endl
       << "    bus           Generate BUS files for single-cell data " << endl
       << "    pseudo        Runs the pseudoalignment step " << endl
       << "    merge         Merges several batch runs " << endl
       << "    h5dump        Converts HDF5-formatted results to plaintext" << endl
       << "    inspect       Inspects and gives information about an index" << endl 
       << "    version       Prints version information" << endl
       << "    cite          Prints citation information" << endl << endl
       << "Running kallisto <CMD> without arguments prints usage information for <CMD>"<< endl << endl;
}

void usageBus() {
  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Generates BUS files for single-cell sequencing" << endl << endl
       << "Usage: kallisto bus [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              pseudoalignment" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl 
       << "-x, --technology=STRING       Single-cell technology used " << endl << endl
       << "Optional arguments:" << endl
       << "-l, --list                    List all single-cell technologies supported" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "-b, --bam                     Input file is a BAM file" << endl
       << "-n, --num                     Output number of read in flag column (incompatible with --bam)" << endl
       << "-T, --tag=STRING              5′ tag sequence to identify UMI reads for certain technologies" << endl
       << "    --fr-stranded             Strand specific reads for UMI-tagged reads, first read forward" << endl
       << "    --rf-stranded             Strand specific reads for UMI-tagged reads, first read reverse" << endl
       << "    --paired                  Treat reads as paired" << endl
       << "    --verbose                 Print out progress information every 1M proccessed reads" << endl;
}

void usageIndex() {
  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Builds a kallisto index" << endl << endl
       << "Usage: kallisto index [arguments] FASTA-files" << endl << endl
       << "Required argument:" << endl
       << "-i, --index=STRING          Filename for the kallisto index to be constructed " << endl << endl
       << "Optional argument:" << endl
       << "-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: " << (Kmer::MAX_K-1) << ")" << endl
       << "    --make-unique           Replace repeated target names with unique names" << endl
       << endl;

}

void usageh5dump() {
  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Converts HDF5-formatted results to plaintext" << endl << endl
       << "Usage:  kallisto h5dump [arguments] abundance.h5" << endl << endl
       << "Required argument:" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl;
}

void usageInspect() {
  cout << "kallisto " << KALLISTO_VERSION << endl << endl
       << "Usage: kallisto inspect INDEX-file" << endl << endl
       << "Optional arguments:" << endl
       << "-G, --gfa=STRING        Filename for GFA output of T-DBG" << endl
       << "-g, --gtf=STRING        Filename for GTF file" << endl
       << "-b, --bed=STRING        Filename for BED output (default: index + \".bed\")" << endl << endl;
 
}

void usageEM(bool valid_input = true) {
  if (valid_input) {

  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Computes equivalence classes for reads and quantifies abundances" << endl << endl;
  }
  //      "----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|"
  cout << "Usage: kallisto quant [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              quantification" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "    --bias                    Perform sequence based bias correction" << endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << endl
       << "    --plaintext               Output plaintext instead of HDF5" << endl
       << "    --fusion                  Search for fusions for Pizzly" << endl
       << "    --single                  Quantify single-end reads" << endl
       << "    --single-overhang         Include reads where unobserved rest of fragment is" << endl
       << "                              predicted to lie outside a transcript" << endl
       << "    --fr-stranded             Strand specific reads, first read forward" << endl
       << "    --rf-stranded             Strand specific reads, first read reverse" << endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << endl
       << "                              (default: -l, -s values are estimated from paired" << endl
       << "                               end data, but are required when using --single)" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "    --pseudobam               Save pseudoalignments to transcriptome to BAM file" << endl
       << "    --genomebam               Project pseudoalignments to genome sorted BAM file" << endl
       << "-g, --gtf                     GTF file for transcriptome information" << endl
       << "                              (required for --genomebam)" << endl
       << "-c, --chromosomes             Tab separated file with chromosome names and lengths" << endl
       << "                              (optional for --genomebam, but recommended)" << endl
       << "    --verbose                 Print out progress information every 1M proccessed reads" << endl;

}

void usagePseudo(bool valid_input = true) {
  if (valid_input) {
    cout << "kallisto " << KALLISTO_VERSION << endl
         << "Computes equivalence classes for reads and quantifies abundances" << endl << endl;
  }

  cout << "Usage: kallisto pseudo [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              pseudoalignment" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "-u  --umi                     First file in pair is a UMI file" << endl
       << "-b  --batch=FILE              Process files listed in FILE" << endl
       << "    --quant                   Quantify using EM algorithm (only in batch mode)" << endl
       << "    --bus                     Output a BUS file" << endl
       << "    --single                  Quantify single-end reads" << endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << endl
       << "                              (default: -l, -s values are estimated from paired" << endl
       << "                               end data, but are required when using --single)" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl;
//       << "    --pseudobam               Output pseudoalignments in SAM format to stdout" << endl;

}

void usageMerge(bool valid_input = true) {
  if (valid_input) {
    cout << "kallisto " << KALLISTO_VERSION << endl
         << "Computes equivalence classes for reads and quantifies abundances" << endl << endl;
  }

  cout << "Usage: kallisto merge [arguments] ouput-directories" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              pseudoalignment" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl;
}

void usageEMOnly() {
  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Computes equivalence classes for reads and quantifies abundance" << endl << endl
       << "Usage: kallisto quant-only [arguments]" << endl << endl
       << "Required argument:" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "-l, --fragment-length=DOUBLE  Estimated fragment length (default: value is estimated from the input data)" << endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << endl
       << "    --plaintext               Output plaintext instead of HDF5" << endl << endl;
}

void usageTCCQuant(bool valid_input = true) {
  if (valid_input) {
    cout << "kallisto " << KALLISTO_VERSION << endl
         << "Quantifies abundance from pre-computed transcript-compatibility counts" << endl << endl;
  }
  
  cout << "Usage: kallisto quant-tcc [arguments] transcript-compatibility-counts-file" << endl << endl
       << "Required arguments:" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used" << endl
       << "                              (required if file with names of transcripts not supplied)" << endl
       << "-T, --txnames=STRING          File with names of transcripts" << endl
       << "                              (required if index file not supplied)" << endl
       << "-e, --ec-file=FILE            File containing equivalence classes" << endl
       << "                              (default: equivalence classes are taken from the index)" << endl
       << "-f, --fragment-file=FILE      File containing fragment length distribution" << endl
       << "                              (default: effective length normalization is not performed)" << endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << endl
       << "                              (note: -l, -s values only should be supplied when" << endl
       << "                               effective length normalization needs to be performed" << endl
       << "                               but --fragment-file is not specified)" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "-g, --genemap                 File for mapping transcripts to genes" << endl
       << "                              (required for obtaining gene-level abundances)" << endl
       << "-G, --gtf=FILE                GTF file for transcriptome information" << endl
       << "                              (can be used instead of --genemap for obtaining gene-level abundances)" << endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << endl;
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
  struct tm * timeinfo;

  time( &rawtime );
  timeinfo = localtime( &rawtime );
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
    string cmd(argv[1]);
    if (cmd == "version") {
      PrintVersion();
    } else if (cmd == "cite") {
      PrintCite();
    } else if (cmd == "index") {
      cerr << endl;
      if (argc==2) {
        usageIndex();
        return 0;
      }
      ParseOptionsIndex(argc-1,argv+1,opt);
      if (!CheckOptionsIndex(opt)) {
        usageIndex();
        exit(1);
      } else {
        // create an index
        Kmer::set_k(opt.k);
        KmerIndex index(opt);
        index.BuildTranscripts(opt);
        index.write(opt.index);
      }
      cerr << endl;
    } else if (cmd == "inspect") {
      if (argc==2) {
        usageInspect();
        return 0;
      }
      ParseOptionsInspect(argc-1, argv+1, opt);
      if (!CheckOptionsInspect(opt)) {
        usageInspect();
        exit(1);
      } else {
        KmerIndex index(opt);
        index.load(opt);
        InspectIndex(index,opt);
      }
    } else if (cmd == "bus") {
      if (argc ==2) {
        usageBus();        
        return 0;
      }
      ParseOptionsBus(argc-1, argv+1,opt);
      if (!CheckOptionsBus(opt)) {
        usageBus();
        exit(1);
      } else {
        int num_trans, index_version;
        int64_t num_processed = 0;
        int64_t num_pseudoaligned = 0;
        int64_t num_unique = 0;

        opt.bus_mode = true;
        opt.single_end = false;
        KmerIndex index(opt);
        index.load(opt);


        bool guessChromosomes = false;
        Transcriptome model; // empty
        if (opt.genomebam) {
          if (!opt.chromFile.empty()) {
            model.loadChromosomes(opt.chromFile);
          } else {
            guessChromosomes = true;
          }
          model.parseGTF(opt.gtfFile, index, opt, guessChromosomes);
        }


        MinCollector collection(index, opt); 
        MasterProcessor MP(index, opt, collection, model);
        num_processed = ProcessBUSReads(MP, opt);

        uint32_t bclen = 0;
        uint32_t umilen = 0;

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
          std::FILE* fp = std::fopen((opt.output + "/output.bus").c_str(), "r+b");
          if (fp != nullptr) {
            std::fseek(fp,8,SEEK_SET); // skip magic string and version
            // write to uint32_t values
            std::fwrite(&bclen, sizeof(bclen),1,fp);
            std::fwrite(&umilen, sizeof(umilen),1,fp);
            std::fclose(fp);
            fp = nullptr;
          }
        }


        writeECList(opt.output + "/matrix.ec", index);

        // write transcript names
        std::ofstream transout_f((opt.output + "/transcripts.txt"));
        for (const auto &t : index.target_names_) {
          transout_f << t << "\n";
        }
        transout_f.close();
        
        std::vector<int> fld;
        if (opt.busOptions.paired && !opt.tagsequence.empty()) {
          fld = collection.flens; // copy
          collection.compute_mean_frag_lens_trunc();
          // Write out index:
          index.write((opt.output + "/index.saved"), false);
          // Write out fragment length distribution:
          std::ofstream flensout_f((opt.output + "/flens.txt"));
          for ( size_t i = 0 ; i < fld.size(); ++i ) {
            if (i != 0) {
              flensout_f << " ";
            }
            flensout_f << fld[i];
          }
          flensout_f << "\n";
          flensout_f.close();
        }

        // gather stats
        num_unique = 0;
        for (int i = 0; i < index.num_trans; i++) {
          num_unique += collection.counts[i];          
        }
        num_pseudoaligned = 0;
        for (int i = 0; i < collection.counts.size(); i++) {
          num_pseudoaligned += collection.counts[i];
        }
        
        // write json file
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
            start_time,
            call);

        if (opt.pseudobam) {
          std::vector<double> fl_means(index.target_lens_.size(),0.0);
          EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
          MP.processAln(em, false);
        }


        cerr << endl;
        if (num_pseudoaligned == 0) {
          exit(1); // exit with error
        }
      }
    } else if (cmd == "merge") {
      if (argc == 2) {
        usageMerge();
        return 0;
      }
      ParseOptionsMerge(argc -1, argv + 1, opt);
      if (!CheckOptionsMerge(opt)) {
        usageMerge();
        exit(1);        
      } else {
        int num_trans, index_version;
        int64_t num_processed, num_pseudoaligned, num_unique;
  
        bool b = MergeBatchDirectories(opt, num_trans, num_processed, num_pseudoaligned, num_unique, index_version);        
        if (!b) {
          exit(1);
        }
        // write json file

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
            start_time,
            call);

      }
    } else if (cmd == "quant") {
      if (argc==2) {
        usageEM();
        return 0;
      }
      ParseOptionsEM(argc-1,argv+1,opt);
      if (!CheckOptionsEM(opt)) {
        cerr << endl;
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
          //model.loadTranscriptome(index, in, opt);
        }


        
        int64_t num_processed = 0;
        int64_t num_pseudoaligned = 0;
        int64_t num_unique = 0;

      
        MinCollector collection(index, opt);        
        MasterProcessor MP(index, opt, collection, model);
        num_processed = ProcessReads(MP, opt);

        // save modified index for future use
        if (opt.write_index) {
          index.write((opt.output + "/index.saved"), false);
        }

        // if mean FL not provided, estimate
        std::vector<int> fld;
        if (opt.fld == 0.0) {
          fld = collection.flens; // copy
          collection.compute_mean_frag_lens_trunc();
        } else {
          auto mean_fl = (opt.fld > 0.0) ? opt.fld : collection.get_mean_frag_len();
          auto sd_fl = opt.sd;
          collection.init_mean_fl_trunc( mean_fl, sd_fl );
          //fld.resize(MAX_FRAG_LEN,0); // no obersvations
          fld = trunc_gaussian_counts(0, MAX_FRAG_LEN, mean_fl, sd_fl, 10000);

          // for (size_t i = 0; i < collection.mean_fl_trunc.size(); ++i) {
          //   cout << "--- " << i << '\t' << collection.mean_fl_trunc[i] << endl;
          // }
        }

        std::vector<int> preBias(4096,1);
        if (opt.bias) {
          preBias = collection.bias5; // copy
        }

        auto fl_means = get_frag_len_means(index.target_lens_, collection.mean_fl_trunc);

        /*for (int i = 0; i < collection.bias3.size(); i++) {
          std::cout << i << "\t" << collection.bias3[i] << "\t" << collection.bias5[i] << "\n";
          }*/

        EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
        em.run(10000, 50, true, opt.bias);

        std::string call = argv_to_string(argc, argv);


        #ifdef USE_HDF5
        H5Writer writer;
        if (!opt.plaintext) {
          writer.init(opt.output + "/abundance.h5", opt.bootstrap, num_processed, fld, preBias, em.post_bias_, 6,
              index.INDEX_VERSION, call, start_time);
          writer.write_main(em, index.target_names_, index.target_lens_);
        }
        #endif // USE_HDF5

        for (int i = 0; i < index.num_trans; i++) {
          num_unique += collection.counts[i];          
        }
        for (int i = 0; i < collection.counts.size(); i++) {
          num_pseudoaligned += collection.counts[i];
        }

        if (num_pseudoaligned == 0) {
          std::cerr << "[~warn] Warning, zero reads pseudoaligned check your input files and index" << std::endl;
        }

        plaintext_aux(
            opt.output + "/run_info.json",
            std::string(std::to_string(index.num_trans)),
            std::string(std::to_string(opt.bootstrap)),
            std::string(std::to_string(num_processed)),
            std::string(std::to_string(num_pseudoaligned)),
            std::string(std::to_string(num_unique)),
            KALLISTO_VERSION,
            std::string(std::to_string(index.INDEX_VERSION)),
            start_time,
            call);

        plaintext_writer(opt.output + "/abundance.tsv", em.target_names_,
            em.alpha_, em.eff_lens_, index.target_lens_);

        if (opt.bootstrap > 0 && num_pseudoaligned == 0) {
          // this happens if nothing aligns, then we write an empty bootstrap file
          for (int b = 0; b < opt.bootstrap; b++) {
            if (!opt.plaintext) {
              #ifdef USE_HDF5
              writer.write_bootstrap(em, b); // em is empty
              #endif
            } else {
              plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".tsv",
                    em.target_names_, em.alpha_, em.eff_lens_, index.target_lens_); // re-use empty input
            }
          }
        } else if (opt.bootstrap > 0 && num_pseudoaligned > 0) {
          auto B = opt.bootstrap;
          std::mt19937_64 rand;
          rand.seed( opt.seed );

          std::vector<size_t> seeds;
          for (auto s = 0; s < B; ++s) {
            seeds.push_back( rand() );
          }

          if (opt.threads > 1) {
            auto n_threads = opt.threads;
            if (opt.threads > opt.bootstrap) {
              cerr
                << "[~warn] number of threads (" << opt.threads <<
                ") greater than number of bootstraps" << endl
                << "[~warn] (cont'd) updating threads to number of bootstraps "
                << opt.bootstrap << endl;
              n_threads = opt.bootstrap;
            }
            #ifdef USE_HDF5            
            BootstrapThreadPool pool(opt.threads, seeds, collection.counts, index,
                collection, em.eff_lens_, opt, &writer, fl_means);
            #endif
          } else {
            for (auto b = 0; b < B; ++b) {
              Bootstrap bs(collection.counts, index, collection, em.eff_lens_, seeds[b], fl_means, opt);
              cerr << "[bstrp] running EM for the bootstrap: " << b + 1 << "\r";
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

          cerr << endl;
        }

        if (opt.pseudobam) {
        
          MP.processAln(em, true);
        }
        

        cerr << endl;
        if (num_pseudoaligned == 0) {
          exit(1); // exit with error
        }
      }
    } else if (cmd == "quant-only") {
      if (argc==2) {
        usageEMOnly();
        return 0;
      }
      ParseOptionsEMOnly(argc-1,argv+1,opt);
      if (!CheckOptionsEM(opt, true)) {
        usageEMOnly();
        exit(1);
      } else {
        // run the em algorithm
        KmerIndex index(opt);
        index.load(opt, false); // skip the k-mer map
        MinCollector collection(index, opt);
        collection.loadCounts(opt);

        std::vector<int> fld;
        // if mean FL not provided, estimate
        if (opt.fld == 0.0) {
          collection.compute_mean_frag_lens_trunc();
          fld = collection.flens;
        } else {
          auto mean_fl = (opt.fld > 0.0) ? opt.fld : collection.get_mean_frag_len();
          auto sd_fl = opt.sd;
          collection.init_mean_fl_trunc( mean_fl, sd_fl );
          //cout << collection.mean_fl_trunc.size() << endl;
          //fld.resize(MAX_FRAG_LEN,0);
          fld = trunc_gaussian_counts(0, MAX_FRAG_LEN, mean_fl, sd_fl, 10000);
          // for (size_t i = 0; i < collection.mean_fl_trunc.size(); ++i) {
          //   cout << "--- " << i << '\t' << collection.mean_fl_trunc[i] << endl;
          // }
        }

        std::vector<int> preBias(4096,1); // default
        if (opt.bias) {
          // fetch the observed bias
          preBias = collection.bias5; // copy
        }

        auto fl_means = get_frag_len_means(index.target_lens_, collection.mean_fl_trunc);

        EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
        em.run(10000, 50, true, opt.bias);

        std::string call = argv_to_string(argc, argv);
        //H5Writer writer;

        if (!opt.plaintext) {
          // setting num_processed to 0 because quant-only is for debugging/special ops
          /* writer.init(opt.output + "/abundance.h5", opt.bootstrap, 0, fld, preBias, em.post_bias_, 6,
              index.INDEX_VERSION, call, start_time);
          writer.write_main(em, index.target_names_, index.target_lens_);*/
        } else {
          plaintext_aux(
              opt.output + "/run_info.json",
              std::string(std::to_string(index.num_trans)),
              std::string(std::to_string(opt.bootstrap)),
              std::string(std::to_string(0)),
              std::string(std::to_string(0)),
              std::string(std::to_string(0)),
              KALLISTO_VERSION,
              std::string(std::to_string(index.INDEX_VERSION)),
              start_time,
              call);

          plaintext_writer(opt.output + "/abundance.tsv", em.target_names_,
              em.alpha_, em.eff_lens_, index.target_lens_);
        }

        int64_t num_pseudoaligned =0;

        for (int i = 0; i < collection.counts.size(); i++) {
          num_pseudoaligned += collection.counts[i];
        }

        if (num_pseudoaligned == 0) {
          std::cerr << "[~warn] Warning, zero reads pseudoaligned check your input files and index" << std::endl;
        }

        if (opt.bootstrap > 0 && num_pseudoaligned == 0) {
          // this happens if nothing aligns, then we write an empty bootstrap file
          for (int b = 0; b < opt.bootstrap; b++) {
            if (!opt.plaintext) {
              //writer.write_bootstrap(em, b); // em is empty
            } else {
              plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".tsv",
                    em.target_names_, em.alpha_, em.eff_lens_, index.target_lens_); // re-use empty input
            }
          }
        } else if (opt.bootstrap > 0 && num_pseudoaligned > 0) {
          std::cerr << "Bootstrapping!" << std::endl;
          auto B = opt.bootstrap;
          std::mt19937_64 rand;
          rand.seed( opt.seed );

          std::vector<size_t> seeds;
          for (auto s = 0; s < B; ++s) {
            seeds.push_back( rand() );
          }

          if (opt.threads > 1) {
            auto n_threads = opt.threads;
            if (opt.threads > opt.bootstrap) {
              cerr << "[btstrp] Warning: number of threads (" << opt.threads <<
                ") greater than number of bootstraps." << endl
                << "[btstrp] (cont'd): Updating threads to number of bootstraps "
                << opt.bootstrap << endl;
              n_threads = opt.bootstrap;
            }

            /*BootstrapThreadPool pool(n_threads, seeds, collection.counts, index,
                collection, em.eff_lens_, opt, writer, fl_means);*/
          } else {
            for (auto b = 0; b < B; ++b) {
              Bootstrap bs(collection.counts, index, collection, em.eff_lens_, seeds[b], fl_means, opt);
              cerr << "[bstrp] running EM for the bootstrap: " << b + 1 << "\r";
              auto res = bs.run_em();

              if (!opt.plaintext) {
                //writer.write_bootstrap(res, b);
              } else {
                plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".tsv",
                    em.target_names_, res.alpha_, em.eff_lens_, index.target_lens_);
              }
            }
          }
        }
        cerr << endl;

        if (num_pseudoaligned == 0) {
          exit(1);
        }
      }
    } else if (cmd == "quant-tcc") {
      if (argc==2) {
        usageTCCQuant();
        return 0;
      }
      ParseOptionsTCCQuant(argc-1,argv+1,opt);
      if (!CheckOptionsTCCQuant(opt)) {
        usageTCCQuant();
        exit(1);
      } else {
        KmerIndex index(opt);
        index.load(opt, false); // skip the k-mer map
        size_t num_ecs = index.ecmap.size();
        MinCollector collection(index, opt);
        std::vector<std::vector<std::pair<int32_t, int32_t>>> batchCounts; // Stores TCCs
        std::ifstream in((opt.tccFile));
        bool firstline = true;
        bool isMatrixFile = false;
        size_t nrow = 0, ncol = 0, nlines = 0;
        int i = 0, prevRow = 0, prevCol = 0;
        if (in.is_open()) { // Read in TCC file
          std::string line;
          while (getline(in, line)) {
            if (firstline) { // First line decides whether the file is in matrix format
              firstline = false;
              if (line.rfind("%%MatrixMarket", 0) == 0) {
                std::cerr << "[tcc] Parsing transcript-compatibility counts (TCC) file as a matrix file" << std::endl;
                isMatrixFile = true;
                while (getline(in, line) && line.rfind("%", 0) == 0) { } // Skip header comments
                std::stringstream ss(line);
                ss >> nrow >> ncol >> nlines;
                std::cerr << "[tcc] Matrix dimensions: " << pretty_num(nrow) << " x " << pretty_num(ncol) << std::endl;
                batchCounts.assign(nrow, {});
                if (opt.bootstrap > 0) {
                  cerr << "[tcc] Warning: Bootstrap can only be used for non-matrix files; will not be performed" << endl;
                }
                continue;
              } else {
                std::cerr << "[tcc] Transcript-compatibility counts (TCC) file is not in matrix format; "
                          << "it will not be parsed as a matrix file" << std::endl;
                isMatrixFile = false;
                batchCounts.assign(1, {});
              }
            }
            std::stringstream ss(line);
            int row, col, val;
            if (isMatrixFile) {
              if (i >= nlines) {
                std::cerr << "[tcc] Warning: TCC matrix file contains additional lines which will not be read; "
                          << "only " << pretty_num(nlines) << " entries, as specified on the first line, will be read." << std::endl;
                break;
              }
              ss >> row >> col >> val;
              if (row > nrow || col > ncol) {
                std::cerr << "Error: TCC matrix file is malformed; " 
                          << "row numbers or column numbers exceed the dimensions of the matrix." << std::endl;
                exit(1);
              }
            } else {
              ss >> col >> val;
              col += 1; // Because we'll consider non-matrix TCC files to be zero-indexed
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
            auto &bc = batchCounts[row-1];
            bc.push_back({col-1, val});
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
        if ((isMatrixFile && ncol != num_ecs) || (!isMatrixFile && ncol > num_ecs)) {
          std::cerr << "Error: number of equivalence classes does not match index. Found "
                    << pretty_num(ncol) << ", expected " << pretty_num(num_ecs) << std::endl;
          exit(1);
        }
        
        std::vector<std::vector<std::pair<int32_t, double>>> Abundance_mat, Abundance_mat_gene, TPM_mat, TPM_mat_gene;
        Abundance_mat.resize(nrow, {});
        Abundance_mat_gene.resize(nrow, {});
        TPM_mat.resize(nrow, {});
        TPM_mat_gene.resize(nrow, {});
        std::vector<std::pair<double, double>> FLD_mat;
        FLD_mat.resize(nrow, {});
        std::vector<std::vector<int>> FLDs;
        Transcriptome model; // empty model

        std::ofstream transout_f((opt.output + "/transcripts.txt"));
        for (const auto &t : index.target_names_) {
          transout_f << t << "\n";
        }
        transout_f.close();
        std::string prefix = opt.output + "/matrix";
        std::string abtsvfilename = opt.output + "/abundance.tsv";
        std::string gene_abtsvfilename = opt.output + "/abundance.gene.tsv";
        std::string genelistname = opt.output + "/genelist.txt";
        std::string abfilename = prefix + ".abundance.mtx";
        std::string abtpmfilename = prefix + ".abundance.tpm.mtx";
        std::string gene_abfilename = prefix + ".abundance.gene.mtx";
        std::string gene_abtpmfilename = prefix + ".abundance.gene.tpm.mtx";
        std::string fldfilename = prefix + ".fld.tsv";

        const bool calcEffLen = !opt.fldFile.empty() || opt.fld != 0.0;
        if (calcEffLen && !opt.fldFile.empty()) { // Parse supplied fragment length distribution file
          std::ifstream infld((opt.fldFile));
          if (infld.is_open()) {
            std::string line;
            while (getline(infld, line)) {
              if (line.empty() || line.rfind("#", 0) == 0) {
                continue; // Ignore empty lines or lines that begin with a #
              }
              std::vector<int> tmp_vec;
              std::stringstream ss(line);
              while(ss.good()) {
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
              if (tmp_vec.size() != MAX_FRAG_LEN) {
                std::cerr << "Error: Fragment length distribution file contains a line with " 
                          << tmp_vec.size() << " values; expected: " << MAX_FRAG_LEN << std::endl;
                exit(1);
              }
              FLDs.push_back(tmp_vec);
            }
            if (FLDs.size() != 1 && FLDs.size() != nrow) {
              std::cerr << "Error: Fragment length distribution file contains " 
                        << FLDs.size() << " valid lines; expected: " << nrow << std::endl;
              exit(1);
            }
          } else {
            std::cerr << "Error: could not open file " << opt.fldFile << std::endl;
            exit(1);
          }
        }

        if (!opt.genemap.empty()) { // Parse supplied gene mapping file
          model.parseGeneMap(opt.genemap, index, opt);
        } else if (!opt.gtfFile.empty()) {
          model.parseGTF(opt.gtfFile, index, opt, true);
        }
        const bool gene_level_counting = !opt.genemap.empty() || !opt.gtfFile.empty();

        std::cerr << "[quant] Running EM algorithm..."; std::cerr.flush();
        auto EM_lambda = [&](int id) {
          MinCollector collection(index, opt);
          collection.counts.assign(index.ecmap.size(), 0);
          const auto& bc = batchCounts[id];
          for (const auto &p : bc) {
            collection.counts[p.first] = p.second;
          }
          
          std::vector<double> fl_means;
          if (calcEffLen) {
            if (opt.fld != 0.0) {
              collection.init_mean_fl_trunc(opt.fld, opt.sd);
            } else {
              std::vector<int> fld;
              if (FLDs.size() == 1) { // Only one distribution supplied; will use this for everything
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
          } else { // Set mean fragment lengths to be target lengths so that effective lengths are all one
            fl_means = std::vector<double>(index.target_lens_.begin(), index.target_lens_.end());
          }

          EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
          em.run(10000, 50, false, false);
          
          if (isMatrixFile) { // Update abundances matrix
            auto &ab_m = Abundance_mat[id];
            auto &tpm_m = TPM_mat[id];
            std::vector<double> gc;
            std::vector<double> gc_tpm;
            if (gene_level_counting) {
              gc.assign(model.genes.size(), 0.0);
              gc_tpm.assign(model.genes.size(), 0.0);
            }
            auto tpm = counts_to_tpm(em.alpha_, em.eff_lens_);
            for (int i = 0; i < em.alpha_.size(); i++) {
              if (em.alpha_[i] > 0.0) {
                ab_m.push_back({i,em.alpha_[i]});
                tpm_m.push_back({i,tpm[i]});
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
          } else { // Write plaintext abundances
            plaintext_writer(abtsvfilename, em.target_names_,
                             em.alpha_, em.eff_lens_, index.target_lens_);
            if (gene_level_counting) {
              plaintext_writer_gene(gene_abtsvfilename, em.target_names_,
                                     em.alpha_, em.eff_lens_, model);
            }
          }
          
          if (opt.bootstrap > 0 && !isMatrixFile) {
            auto B = opt.bootstrap;
            std::mt19937_64 rand;
            rand.seed( opt.seed );
            std::vector<size_t> seeds;
            for (auto s = 0; s < B; ++s) {
              seeds.push_back( rand() );
            }
            cerr << endl;
            for (auto b = 0; b < B; ++b) {
              Bootstrap bs(collection.counts, index, collection, em.eff_lens_, seeds[b], fl_means, opt);
              cerr << "[bstrp] running EM for the bootstrap: " << b + 1 << "\r";
              auto res = bs.run_em();
              plaintext_writer(opt.output + "/bs_abundance_" + std::to_string(b) + ".tsv",
                               em.target_names_, res.alpha_, em.eff_lens_, index.target_lens_);
            }
            cerr << endl;
          }
        }; // end of EM_lambda
        
        std::vector<std::thread> workers;
        int num_ids = nrow;
        int id = 0;
        while (id < num_ids) {
          workers.clear();
          int nt = std::min(opt.threads, (num_ids - id));
          int first_id = id;
          for (int i = 0; i < nt; i++,id++) {
            workers.emplace_back(std::thread(EM_lambda, id));
          }
          for (int i = 0; i < nt; i++) {
            workers[i].join();
          }
        }
        std::cerr << " done" << std::endl;
        
        cerr << endl;
        

        if (isMatrixFile) {
          writeSparseBatchMatrix(abfilename, Abundance_mat, index.num_trans);
          writeSparseBatchMatrix(abtpmfilename, TPM_mat, index.num_trans);
          if (gene_level_counting) {
            writeSparseBatchMatrix(gene_abfilename, Abundance_mat_gene, model.genes.size());
            writeSparseBatchMatrix(gene_abtpmfilename, TPM_mat_gene, model.genes.size());
            writeGeneList(genelistname, model);
          }
        }
        if (calcEffLen) {
          writeFLD(fldfilename, FLD_mat);
        }
      }
    } else if (cmd == "pseudo") {
      if (argc==2) {
        usagePseudo();
        return 0;
      }
      ParseOptionsPseudo(argc-1,argv+1,opt);
      if (!CheckOptionsPseudo(opt)) {
        cerr << endl;
        usagePseudo(false);
        exit(1);
      } else {
        // pseudoalign the reads
        KmerIndex index(opt);
        index.load(opt);

        MinCollector collection(index, opt);
        int64_t num_processed = 0;
        int64_t num_pseudoaligned = 0;
        int64_t num_unique = 0;

        Transcriptome model; // empty model
        if (!opt.gtfFile.empty()) {          
          model.parseGTF(opt.gtfFile, index, opt, true);
        }
        MasterProcessor MP(index, opt, collection, model);
        if (!opt.batch_mode) {
          num_processed = ProcessReads(MP, opt);
          collection.write((opt.output + "/pseudoalignments"));
        } else {          
          num_processed = ProcessBatchReads(MP,opt);
        }

        std::string call = argv_to_string(argc, argv);


        for (int id = 0; id < MP.batchCounts.size(); id++) {
          const auto &cc = MP.batchCounts[id];
          for (const auto &p : cc) {
            if (p.first < index.num_trans) {
              num_unique += p.second;
            }
            num_pseudoaligned += p.second;
          }
        }
        
        std::ofstream transout_f((opt.output + "/transcripts.txt"));
        for (const auto &t : index.target_names_) {
          transout_f << t << "\n";
        }
        transout_f.close();

        plaintext_aux(
            opt.output + "/run_info.json",
            std::string(std::to_string(index.num_trans)),
            std::string(std::to_string(0)), // no bootstraps in pseudo
            std::string(std::to_string(num_processed)),
            std::string(std::to_string(num_pseudoaligned)),
            std::string(std::to_string(num_unique)),
            KALLISTO_VERSION,
            std::string(std::to_string(index.INDEX_VERSION)),
            start_time,
            call);

        
        
        std::vector<std::vector<std::pair<int32_t, double>>> Abundance_mat;
        std::vector<std::pair<double, double>> FLD_mat;
          
        if (opt.pseudo_quant) {
          int n_batch_files = opt.batch_files.size();
          Abundance_mat.resize(n_batch_files, {});
          FLD_mat.resize(n_batch_files, {});

          std::cerr << "[quant] Running EM algorithm for each cell .."; std::cerr.flush();
          
          auto EM_lambda = [&](int id) {          
            MinCollector collection(index, opt);
            collection.flens = MP.batchFlens[id];
            collection.counts.assign(index.ecmap.size(), 0);
            const auto& bc = MP.batchCounts[id];
            for (const auto &p : bc) {
              collection.counts[p.first] = p.second;
            }
            // if mean FL not provided, estimate
            std::vector<int> fld;
            if (opt.fld == 0.0) {
              fld = collection.flens; // copy
              collection.compute_mean_frag_lens_trunc(false);
            } else {
              auto mean_fl = (opt.fld > 0.0) ? opt.fld : collection.get_mean_frag_len(true);
              if (mean_fl == std::numeric_limits<double>::max()) {
                std::cerr << "Couldn't estimate fragment length for batch file number " << id << std::endl;
                return;
              }
              auto sd_fl = opt.sd;
              collection.init_mean_fl_trunc( mean_fl, sd_fl );
              fld = trunc_gaussian_counts(0, MAX_FRAG_LEN, mean_fl, sd_fl, 10000);
            }
            std::vector<int> preBias(4096,1);
            if (opt.bias) {
              //preBias = collection.bias5; // copy
            }

            auto fl_means = get_frag_len_means(index.target_lens_, collection.mean_fl_trunc);

            EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
            em.run(10000, 50, false, opt.bias);

            auto &ab_m = Abundance_mat[id];
            for (int i = 0; i < em.alpha_.size(); i++) {
              if (em.alpha_[i] > 0.0) {
                ab_m.push_back({i,em.alpha_[i]});
              }
            }

            double mean_fl = collection.get_mean_frag_len(true);
            double sd_fl = collection.get_sd_frag_len();
            FLD_mat[id] = {mean_fl, sd_fl};
          }; // end of EM_lambda

          std::vector<std::thread> workers;
          int num_ids = opt.batch_ids.size();
          int id =0;
          while (id < num_ids) {
            workers.clear();
            int nt = std::min(opt.threads, (num_ids - id));
            int first_id = id;
            for (int i = 0; i < nt; i++,id++) {
              workers.emplace_back(std::thread(EM_lambda, id));
              //workers.emplace_back(std::thread(ReadProcessor(index, opt, tc, *this, id,i)));
            }
            
            for (int i = 0; i < nt; i++) {
              workers[i].join();
            }
          }

          std::cerr << " done" << std::endl;


        }
        cerr << endl;

        std::string prefix = opt.output + "/matrix";
        std::string ecfilename = prefix + ".ec";
        std::string tccfilename = prefix + ".tcc.mtx";
        std::string abfilename = prefix + ".abundance.mtx";
        std::string cellnamesfilename = prefix + ".cells";
        std::string fldfilename = prefix + ".fld.tsv";
        std::string genelistname = prefix + ".genelist.txt";
        std::string genecountname = prefix + ".genes.mtx";
        std::string busbarcodelistname = prefix + ".barcodes";
        std::string busoutputname = opt.output + "/output.bus";

        writeECList(ecfilename, index);
        writeCellIds(cellnamesfilename, opt.batch_ids);
        if (opt.batch_bus) {
          // Write BUS file
          writeBUSMatrix(busoutputname, MP.batchCounts, index.ecmap.size());
          if (!MP.batchCounts.empty()) {
            // Write out fake barcodes that identify each cell
            std::vector<std::string> fake_bcs;
            fake_bcs.reserve(MP.batchCounts.size());
            for (size_t j = 0; j < MP.batchCounts.size(); j++) {
              fake_bcs.push_back(binaryToString(j, BUSFORMAT_FAKE_BARCODE_LEN));
            }
            writeCellIds(busbarcodelistname, fake_bcs);
          }
          // Write out index:
          index.write((opt.output + "/index.saved"), false);
          // Write out fragment length distributions if reads paired-end:
          if (!opt.single_end) {
            std::ofstream flensout_f((opt.output + "/flens.txt"));
            for (size_t id = 0; id < opt.batch_ids.size(); id++) {
              std::vector<int> fld = MP.batchFlens[id];
              for ( size_t i = 0 ; i < fld.size(); ++i ) {
                if (i != 0) {
                  flensout_f << " ";
                }
                flensout_f << fld[i];
              }
              flensout_f << "\n";
            }
            flensout_f.close();
          }
        } else {
          writeSparseBatchMatrix(tccfilename, MP.batchCounts, index.ecmap.size());
        }
        if (opt.pseudo_quant) {
          writeSparseBatchMatrix(abfilename, Abundance_mat, index.num_trans);
          writeFLD(fldfilename, FLD_mat);
        }
        if (!opt.gtfFile.empty()) {
          // write out gene info
          std::vector<std::vector<std::pair<int32_t, double>>> geneCounts;
          geneCounts.assign(MP.batchCounts.size(), {});
          
          std::unordered_set<int> gene_ids;
          gene_ids.reserve(100);
          int n_batch_files = opt.batch_files.size();
          std::vector<double> gc;

          for (int id = 0; id < n_batch_files; id++) {
            auto& sgc = geneCounts[id];
            gc.assign(model.genes.size(), 0.0);  
            const auto& bc = MP.batchCounts[id];
            for (auto &p : bc) {
              int ec = p.first;
              if (ec < 0) {
                continue; 
              }
              if (ec < index.num_trans) {
                int g_id = model.transcripts[ec].gene_id;
                if (g_id != -1) {
                  gc[g_id] += p.second;
                }
              } else {
                gene_ids.clear();
                for (auto t : index.ecmap[ec]) {
                  int g_id = model.transcripts[t].gene_id;
                  if (g_id != -1) {
                    gene_ids.insert(g_id);
                  }
                }
                if (!gene_ids.empty()) {
                  double n_genes = gene_ids.size();
                  for (auto &g_id : gene_ids) {
                    gc[g_id] += p.second / n_genes;
                  }
                }
              }
            }

            for (int j = 0; j < gc.size(); j++) {
              if (gc[j] > 0.0) {
                sgc.push_back({j, gc[j]});                
              }
            }
          }


         

          writeGeneList(genelistname, model);
          writeSparseBatchMatrix(genecountname, geneCounts, model.genes.size());
        }


        if (opt.pseudobam) {       
          std::vector<double> fl_means(index.target_lens_.size(),0.0);
          EMAlgorithm em(collection.counts, index, collection, fl_means, opt);
          MP.processAln(em, false);
        }
      }

      
    } else if (cmd == "h5dump") {

      if (argc == 2) {
        usageh5dump();
        exit(1);
      }

      ParseOptionsH5Dump(argc-1,argv+1,opt);
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
    }  else {
      cerr << "Error: invalid command " << cmd << endl;
      usage();
      exit(1);
    }

  }

  fflush(stdout);

  return 0;
}
