#include <vector>
#include <fstream>
#include "BUSData.h"

void writeBUSHeader(std::ofstream &out, int bclen, int umilen) {
  out.write("BUS\0", 4);
  out.write((char*)(&BUSFORMAT_VERSION), sizeof(BUSFORMAT_VERSION));
  out.write((char*)(&bclen), sizeof(bclen));
  out.write((char*)(&umilen), sizeof(umilen));
  std::string header_text = "BUS file produced by kallisto";
  uint32_t len = header_text.size();
  out.write((char*)(&len),sizeof(len));  
  out.write(header_text.c_str(), len);
}

void writeBUSData(std::ofstream &out, const std::vector<BUSData> &bv) {
  for (const auto &b : bv) {
    if (b.ec != -1) { // maybe keep non-mapping reads ?!?
      out.write((char*)(&b), sizeof(b));
    }
  }
}