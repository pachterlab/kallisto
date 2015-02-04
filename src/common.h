#ifndef KALLISTO_COMMON_H
#define KALLISTO_COMMON_H


#define KALLISTO_VERSION "0.1"


#include <string>
#include <vector>




struct ProgramOptions {
	bool verbose;
	int threads;
	std::string index;
	int k;
	int iterations;
	std::string output;
	size_t seed;
	std::string transfasta;
	std::vector<std::string> files;

ProgramOptions() : verbose(false), seed(0), threads(1), k(21), iterations(500)  {}
};

#endif // KALLISTO_COMMON_H
