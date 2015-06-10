#ifndef KALLISTO_PROCESSREADS_H
#define KALLISTO_PROCESSREADS_H

#include <zlib.h>
#include "kseq.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

#include <iostream>
#include <fstream>
#include "MinCollector.h"


#include "common.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif




//void printVector(const std::vector<int>& v, std::ostream& o);

//bool isSubset(const std::vector<int>& x, const std::vector<int>& y);
void ProcessReads(KmerIndex& index, const ProgramOptions& opt, MinCollector& tc);


#endif // KALLISTO_PROCESSREADS_H
