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

void ParseOptionsPseudo(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int single_flag = 0;
  int strand_flag = 0;
  int pbam_flag = 0;
  int gbam_flag = 0;
  int umi_flag = 0;
  int quant_flag = 0;

  const char *opt_string = "t:i:l:s:o:b:u:g:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"single", no_argument, &single_flag, 1},
    //{"strand-specific", no_argument, &strand_flag, 1},
    {"pseudobam", no_argument, &pbam_flag, 1},
    {"quant", no_argument, &quant_flag, 1},
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
  << "CELSeq           CEL-Seq" << endl
  << "CELSeq2          CEL-Seq version 2" << endl
  << "DropSeq          DropSeq" << endl
  << "inDropsv1        inDrops version 1 chemistry" << endl
  << "inDropsv2        inDrops version 2 chemistry" << endl
  << "inDropsv3        inDrops version 3 chemistry" << endl
  << "SCRBSeq          SCRB-Seq" << endl
  << "SureCell         SureCell for ddSEQ" << endl
  << endl;
 }

void ParseOptionsBus(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int gbam_flag = 0;

  const char *opt_string = "i:o:x:t:lbng:c:";
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
  if (v.size() != 1) {
    errorList.push_back("Error: only a single UMI list allow " + umistr);
    return false;
  }
  busopt.umi = std::move(v.front());


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
      if (opt.technology == "10XV2") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0)); // second file, entire string
        busopt.umi = BUSOptionSubstr(0,16,26); // first file [16:26]
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
      } else if (opt.technology == "10XV3") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,16,28);
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
//      } else if (opt.technology == "10XV1") {

      } else if (opt.technology == "SURECELL") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,18,26);
        busopt.bc.push_back(BUSOptionSubstr(0,0,18));
      } else if (opt.technology == "DROPSEQ") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,12,20);
        busopt.bc.push_back(BUSOptionSubstr(0,0,12));
      } else if (opt.technology == "INDROPSV1") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,42,48);
        busopt.bc.push_back(BUSOptionSubstr(0,0,11));
        busopt.bc.push_back(BUSOptionSubstr(0,30,38));  
      } else if (opt.technology == "INDROPSV2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(0,0,0));
        busopt.umi = BUSOptionSubstr(1,42,48);
        busopt.bc.push_back(BUSOptionSubstr(1,0,11));
        busopt.bc.push_back(BUSOptionSubstr(1,30,38));  
      } else if (opt.technology == "INDROPSV3") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi = BUSOptionSubstr(1,8,14);
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
        busopt.bc.push_back(BUSOptionSubstr(1,0,8));
      } else if (opt.technology == "CELSEQ") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,8,12);
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
      } else if (opt.technology == "CELSEQ2") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,0,6);
        busopt.bc.push_back(BUSOptionSubstr(0,6,12));
      } else if (opt.technology == "SCRBSEQ") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,6,16);
        busopt.bc.push_back(BUSOptionSubstr(0,0,6));
      } else if (opt.technology == "INDROPSV3") {
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,0,6);
        busopt.bc.push_back(BUSOptionSubstr(0,6,16));
      } else {
        vector<int> files;
        vector<BUSOptionSubstr> values;
        vector<BUSOptionSubstr> bcValues;
        vector<std::string> errorList;
        //bool invalid = ParseTechnology(opt.technology, values, files, errorList, bcValues);
        bool valid = ParseTechnology(opt.technology, busopt, errorList);
        
        if(!valid) {
          /*
          busopt.nfiles = files.size(); 
          for(int i = 0; i < bcValues.size(); i++) {
            busopt.bc.push_back(bcValues[i]);
          }
          busopt.umi = values[0];
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
      if (opt.technology == "10XV2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0)); // second file, entire string
        busopt.umi = BUSOptionSubstr(0,16,26); // first file [16:26]
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
      } else if (opt.technology == "10XV3") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,16,28);
        busopt.bc.push_back(BUSOptionSubstr(0,0,16));
      } else if (opt.technology == "10XV1") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi = BUSOptionSubstr(1,0,10);
        busopt.bc.push_back(BUSOptionSubstr(0,0,14));
      } else if (opt.technology == "SURECELL") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,51,59);
        busopt.bc.push_back(BUSOptionSubstr(0,0,6));
        busopt.bc.push_back(BUSOptionSubstr(0,21,27));
        busopt.bc.push_back(BUSOptionSubstr(0,42,48));
      } else if (opt.technology == "DROPSEQ") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,12,20);
        busopt.bc.push_back(BUSOptionSubstr(0,0,12));
      } else if (opt.technology == "INDROPSV1") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,42,48);
        busopt.bc.push_back(BUSOptionSubstr(0,0,11));
        busopt.bc.push_back(BUSOptionSubstr(0,30,38));  
      } else if (opt.technology == "INDROPSV2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(0,0,0));
        busopt.umi = BUSOptionSubstr(1,42,48);
        busopt.bc.push_back(BUSOptionSubstr(1,0,11));
        busopt.bc.push_back(BUSOptionSubstr(1,30,38));  
      } else if (opt.technology == "INDROPSV3") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi = BUSOptionSubstr(1,8,14);
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
        busopt.bc.push_back(BUSOptionSubstr(1,0,8));
      } else if (opt.technology == "CELSEQ") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,8,12);
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
      } else if (opt.technology == "CELSEQ2") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,0,6);
        busopt.bc.push_back(BUSOptionSubstr(0,6,12));
      } else if (opt.technology == "SCRBSEQ") {
        busopt.nfiles = 2;
        busopt.seq.push_back(BUSOptionSubstr(1,0,0));
        busopt.umi = BUSOptionSubstr(0,6,16);
        busopt.bc.push_back(BUSOptionSubstr(0,0,6));
      } else if (opt.technology == "INDROPSV3") {
        busopt.nfiles = 3;
        busopt.seq.push_back(BUSOptionSubstr(2,0,0));
        busopt.umi = BUSOptionSubstr(1,8,14);
        busopt.bc.push_back(BUSOptionSubstr(0,0,8));
        busopt.bc.push_back(BUSOptionSubstr(1,0,8));
      } else {
        vector<int> files;
        vector<BUSOptionSubstr> values;
        vector<BUSOptionSubstr> bcValues;
        vector<std::string> errorList;        
        //bool invalid = ParseTechnology(opt.technology, values, files, errorList, bcValues);
        bool valid = ParseTechnology(opt.technology, busopt, errorList);
        
        
        if(valid) {
          /*
          busopt.nfiles = files.size(); 
          for(int i = 0; i < bcValues.size(); i++) {
            busopt.bc.push_back(bcValues[i]);
          }
          busopt.umi = values[0];
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

  if (opt.genomebam) {
    if (opt.busOptions.seq.size() != 1) {
      cerr << "Error: BAM output is currently only supported for technologies with a single CDNA read file" << endl;
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
  }

  // check for read files
  if (!opt.batch_mode) {
    if (opt.umi) {
      cerr << ERROR_STR << " UMI must be run in batch mode, use --batch option" << endl;
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

        writeECList(ecfilename, index);
        writeCellIds(cellnamesfilename, opt.batch_ids);
        writeSparseBatchMatrix(tccfilename, MP.batchCounts, index.ecmap.size());
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
