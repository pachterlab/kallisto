#ifndef KALLISTO_BUSDATA_H
#define KALLISTO_BUSDATA_H

#include <string>
#include <vector>
#include <stdint.h>

struct BUSTranscript {
  std::string name;
  uint32_t transcriptLength;
  BUSTranscript() : transcriptLength(0) {}  
};


struct BUSHeader {
  std::string text;
  std::vector<BUSTranscript> transcripts;
  std::vector<std::vector<int32_t>> ecs;
};

struct BUSData {
  uint32_t barcode;
  uint32_t UMI;
  int32_t ec;
  uint32_t count;
  uint32_t flags;

  BUSData() : barcode(0), UMI(0), ec(-1), count(0), flags(0) {}
};

uint32_t stringToBinary(const std::string &s, uint32_t &flag);
uint32_t stringToBinary(const char* s, const size_t len, uint32_t &flag);

#endif // KALLISTO_BUSDATA_H