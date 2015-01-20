#ifndef KALLISTO_KMERINDEX_H
#define KALLISTO_KMERINDEX_H

#include <vector>

#include "common.h"


struct KmerIndex
{
	KmerIndex(const ProgramOptions& opt) {}
	~KmerIndex() {}

	void match(const char *s, int l, std::vector<int> & v) const {}

	
	
	
	
};

#endif // KALLISTO_KMERINDEX_H
