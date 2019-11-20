#ifndef KALLISTO_BUSDATA_H
#define KALLISTO_BUSDATA_H

#include <string>
#include <vector>
#include <stdint.h>

const uint32_t BUSFORMAT_VERSION = 1;

struct BUSTranscript {
  std::string name;
  uint32_t transcriptLength;
  BUSTranscript() : transcriptLength(0) {}  
};


struct BUSHeader {
  std::string text;
  std::vector<BUSTranscript> transcripts;
  std::vector<std::vector<int32_t>> ecs;
  uint32_t version;
  uint32_t bclen;
  uint32_t umilen;
  BUSHeader() : version(0), bclen(0), umilen(0) {}

};

struct BUSData {
  uint64_t barcode;
  uint64_t UMI;
  int32_t ec;
  uint32_t count;
  uint32_t flags;  
  uint32_t pad;
  BUSData() : barcode(0), UMI(0), ec(-1), count(0), flags(0), pad(0) {}
};

uint64_t stringToBinary(const std::string &s, uint32_t &flag);
uint64_t stringToBinary(const char* s, const size_t len, uint32_t &flag);
std::string binaryToString(uint64_t x, size_t len);
#endif // KALLISTO_BUSDATA_H