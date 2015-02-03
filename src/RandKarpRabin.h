#ifndef RANDKARPRABIN_H
#define RANDKARPRABIN_H

#include <stdint.h>
#include <cassert>

struct RandKarpRabin {
	RandKarpRabin(int _k) : k(_k) //, charmask (31)
  {
		setK(k&63);
    seed(0);
  }

  RandKarpRabin() : k(0) //, charmask(31)
  {
    seed(0);
  }

	uint64_t hash() const {
		return (h ^ ht);
	}

  void seed(int s) {}

	void init(const char* _s) {
    const unsigned char* s = (const unsigned char*) _s;
		uint64_t sh = 1ULL;
    for (size_t i = 0; i < k; i++) {
			h <<= 2;
			h |= charmask(s[i]);
			ht |= twinmask(s[i])*sh;
			sh <<= 2;
    }
  }

	inline void update(const unsigned char out, const unsigned char in) {
		h = (h << 2) & bitmask;
		h |= charmask(in);
		ht >>= 2;
		ht |= twinmask(out) << (2*(k-1));
	}


	void setK(int _k) {
		h=0;
		ht=0;
		k = _k;
		bitmask = (1ULL << (2*k))-1;
	}

	inline uint64_t charmask(unsigned char x) {
    return (x&6)>>1;
  }

  inline uint64_t twinmask(unsigned char x) {
    return ((x^4)&6)>>1;
  }

	size_t k;
	uint64_t a = 14199033468433001729ULL;
	uint64_t b = 12279310770576267508ULL;
	uint64_t bitmask;
	uint64_t h;
	uint64_t ht;
};

#endif // RANDKARPRABIN_H
