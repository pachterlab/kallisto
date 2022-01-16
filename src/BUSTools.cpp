#include <vector>
#include <fstream>
#include "BUSTools.h"

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

void writeBUSMatrix(const std::string &filename,
                    std::vector<std::vector<std::pair<int,int>>> &data, int cols) {
  std::ofstream of;
  of.open(filename.c_str(), std::ios::out | std::ios::binary);
  writeBUSHeader(of, BUSFORMAT_FAKE_BARCODE_LEN, 1);
  if (!data.empty()) { // reduntant  
    for (size_t j = 0; j < data.size(); j++) {
      const auto &v = data[j];
      for (size_t i = 0; i < v.size(); i++) {
        if (v[i].second != 0 && v[i].first != -1) {
          BUSData b;
          b.barcode = j;
          b.flags = 0;
          b.UMI = -1;
          b.ec = v[i].first;
          b.count = v[i].second;
          of.write((char*)(&b), sizeof(b));
        }
      }
    }
  }
  of.close();
}