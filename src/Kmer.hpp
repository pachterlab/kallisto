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

  void set_empty();
  void set_deleted();


  bool operator<(const Kmer& o) const;

  // use:  b = (km1 == km2);
  // pre:
  // post: b is true <==> the DNA strings in km1 and km2 are equal
  inline bool operator==(const Kmer& o) const {
    for (size_t i = 0; i < MAX_K/32; i++) {
      if (longs[i] != o.longs[i]) {
        return false;
      }
    }
    return true;
    //  return memcmp(bytes,o.bytes,MAX_K/4)==0;
  }

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

  std::string getBinary() const;

  void toString(char *s) const;
  std::string toString() const;

  // static functions
  static void set_k(unsigned int _k);


  static const unsigned int MAX_K = MAX_KMER_SIZE;
  static unsigned int k;

 private:
  static unsigned int k_bytes;
  static unsigned int k_longs;
  static unsigned int k_modmask; // int?

  // data fields
  union {
    uint8_t bytes[MAX_K/4];
    uint64_t longs[MAX_K/32];
  };
  // By default MAX_K == 64 so the union uses 16 bytes
  // However sizeof(Kmer) == 24
  // Are the 8 extra bytes alignment?

  // private functions
//void shiftForward(int shift);

//void shiftBackward(int shift);

};


struct KmerHash {
  size_t operator()(const Kmer& km) const {
    return km.hash();
  }
};


#endif // BFG_KMER_HPP
