#ifndef SIMPLEKMERHASH_H
#define SIMPLEKMERHASH_H

#include <stdint.h>
#include <cassert>

extern  const uint64_t twin_table[256];

struct SimpleKmerHash {

	SimpleKmerHash(int _k) {
		assert(_k == Kmer::k && "k has to match");
	}

	uint64_t hash() const {
		return km.hash() ^ tw.hash();
	}

	void seed(int s) {}

	void init(const char* s) {
		km.set_kmer(s);
		tw = km.twin();
	}

	inline void update(const unsigned char out, const unsigned char in) {
		km.forwardBase(in);
		tw.backwardBase((unsigned char)(twin_table[out]));
	}

	void setK(int _k) {
		assert(_k == Kmer::k && "k has to match");
	}
	
	Kmer km,tw;
};

#endif // SIMPLEKMERHASH_H
