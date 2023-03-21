#ifndef KALLISTO_BUSTOOLS_H
#define KALLISTO_BUSTOOLS_H

#include <fstream>
#include <vector>
#include "MinCollector.h"
#include "BUSData.h"

const uint32_t BUSFORMAT_FAKE_BARCODE_LEN = 16;

size_t writeBUSData(std::ofstream &out, const std::vector<BUSData> &bv, MinCollector* tc = nullptr);
void writeBUSHeader(std::ofstream &out, int bclen, int umilen);
void writeBUSMatrix(const std::string &filename,
                    std::vector<std::vector<std::pair<int,int>>> &data, int cols, std::vector<int> mapping);

#endif // KALLISTO_BUSTOOLS_H
