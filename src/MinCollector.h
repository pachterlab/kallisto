#ifndef KALLISTO_MINCOLLECTOR_H
#define KALLISTO_MINCOLLECTOR_H

#include "common.h"
#include <iostream>
#include <vector>

struct MinCollector {
	
	MinCollector(const ProgramOptions& opt) {}

	void collect(const std::vector<int>& v) {}
	void write(std::ostream& o) {
		o << "TEST";
	}
	
};

#endif // KALLISTO_MINCOLLECTOR_H
