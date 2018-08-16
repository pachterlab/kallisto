#include <vector>
#include <fstream>
#include "BUSData.h"

void writeBUSData(std::ofstream &out, const std::vector<BUSData> &bv) {
  for (const auto &b : bv) {
    if (b.ec != -1) { // maybe keep non-mapping reads ?!?
      out.write((char*)(&b), sizeof(b));
    }
  }
}