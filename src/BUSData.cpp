#include "BUSData.h"


uint32_t stringToBinary(const std::string &s, uint32_t &flag) {
  return stringToBinary(s.c_str(), s.size(), flag);
}

uint32_t stringToBinary(const char* s, const size_t len, uint32_t &flag) {
  uint32_t r = 0;
  flag = 0;
  int numN = 0;
  size_t posN = 0;
  size_t k = len;
  if (k > 16) {
    k = 16;
  }
  for (size_t i = 0; i < k; ++i) {
    uint32_t x = ((*s) & 4) >> 1;
    if (((*s) & 3) == 2) {
      if (numN == 0) {
        posN = i;
      }
      ++numN;
    }
    r |= ((x + ((x ^ (*s & 2)) >>1)) << (2*i));
    s++;
  }
  if (numN>0) {
    if (numN > 3) {
      numN = 3;      
    }    
    flag = (numN & 3) | (posN & 15) << 2;
  }
  return r;
}