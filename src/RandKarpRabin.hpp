#ifndef RANDKARPRABIN_H
#define RANDKARPRABIN_H

#include <stdint.h>
#include <cassert>

struct RandKarpRabin {
	RepHash(int _k) : k(_k) //, charmask (31)
  {
    k = k & 63;
    seed(0);
  }

  RepHash() : k(0) //, charmask(31)
  {
    seed(0);
  }

  void seed(int s) {}

	void init(const char* _s) {
    const unsigned char* s = (const unsigned char*) _s;
    for (size_t i = 0; i < k; i++) {
			// finish
    }
  }

	inline void update(const unsigned char out, const unsigned char in) {}


	void setK(int _k) {}

	inline uint64_t charmask(unsigned char x) {
    return (x&6)>>1;
  }

  inline uint64_t twinmask(unsigned char x) {
    return ((x^4)&6)>>1;
  }

	size_t k;
	uint64_t a = 14199033468433001729ULL;
	uint64_t b = 12279310770576267508ULL;	
}

#endif // RANDKARPRABIN_H
