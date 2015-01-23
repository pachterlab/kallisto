#ifndef BFG_KMER_HPP
#define BFG_KMER_HPP

#ifndef MAX_KMER_SIZE
 #define MAX_KMER_SIZE 32
#endif

#include <stdio.h>
#include <stdint.h>
#include <cassert>
#include <cstring>
#include <string>

#include "hash.hpp"


//for debug
void int2bin(uint32_t a, char *buffer, int buf_size);


/* Short description: 
 *  - Store kmer strings by using 2 bits per base instead of 8 
 *  - Easily return reverse complements of kmers, e.g. TTGG -> CCAA
 *  - Easily compare kmers
 *  - Provide hash of kmers
 *  - Get last and next kmer, e.g. ACGT -> CGTT or ACGT -> AACGT
 *  */
class Kmer {
 public:

  Kmer();
  Kmer(const Kmer& o);
  explicit Kmer(const char *s); 
  

  Kmer& operator=(const Kmer& o);
  
  void set_deleted();

  bool operator<(const Kmer& o) const;

  bool operator==(const Kmer& o) const;

  bool operator!=(const Kmer& o) const {
    return !(*this == o);
  }

  void set_kmer(const char *s);

  uint64_t hash() const;

  

  Kmer twin() const;
  Kmer rep() const;

  Kmer getLink(const size_t index) const;

  Kmer forwardBase(const char b) const;

  Kmer backwardBase(const char b) const;
  
  void printBinary() const;
  
  void toString(char * s) const;
  std::string toString() const;

  // static functions
  static void set_k(unsigned int _k);


  static const unsigned int MAX_K = MAX_KMER_SIZE;
  static unsigned int k;

 private:
  static unsigned int k_bytes;
  //  static unsigned int k_longs;
  static unsigned int k_modmask;

  // data fields
  union {
    uint8_t bytes[MAX_K/4];
    //uint32_t longs[MAX_K/16];
  };


  // private functions
  void shiftForward(int shift);
  
  void shiftBackward(int shift);

};


struct KmerHash {
  size_t operator()(const Kmer &km) const {
    return km.hash();
  }
};


#endif // BFG_KMER_HPP
