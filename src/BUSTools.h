#ifndef KALLISTO_BUSTOOLS_H
#define KALLISTO_BUSTOOLS_H

#include <fstream>
#include <vector>
#include "BUSData.h"

void writeBUSData(std::ofstream &out, const std::vector<BUSData> &bv);

#endif // KALLISTO_BUSTOOLS_H