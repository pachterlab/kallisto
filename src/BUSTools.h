#ifndef KALLISTO_BUSTOOLS_H
#define KALLISTO_BUSTOOLS_H

#include <fstream>
#include <vector>
#include "BUSData.h"

void writeBUSData(std::ofstream &out, const std::vector<BUSData> &bv);
void writeBUSHeader(std::ofstream &out, int bclen, int umilen);

#endif // KALLISTO_BUSTOOLS_H