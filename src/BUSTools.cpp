#include <vector>
#include <fstream>
#include "BUSData.h"

void writeBUSHeader(std::ofstream &out, int bclen, int umilen) {
  out.write((char*)(&bclen), sizeof(bclen));
  out.write((char*)(&umilen), sizeof(umilen));
}

void writeBUSData(std::ofstream &out, const std::vector<BUSData> &bv) {
  for (const auto &b : bv) {
    if (b.ec != -1) { // maybe keep non-mapping reads ?!?
      out.write((char*)(&b), sizeof(b));
    }
  }
}