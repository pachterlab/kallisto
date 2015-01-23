#include <stdint.h>
#include <cstring>
#include "hash.hpp"

uint64_t inline _rotl64(uint64_t value, int8_t amount) {
  return ((value) << (amount)) | ((value) >> (64 - (amount)));
}

uint32_t SuperFastHash (const char * data, int len) {
uint32_t hash = len, tmp;
int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}




//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

inline uint64_t getblock ( const uint64_t * p, int i )
{
  return p[i];
}

//----------
// Block mix - combine the key bits with the hash bits and scramble everything

inline void bmix64 ( uint64_t & h1, uint64_t & h2, uint64_t & k1, uint64_t & k2, uint64_t & c1, uint64_t & c2 )
{
  k1 *= c1; 
  k1  = _rotl64(k1,23); 
  k1 *= c2;
  h1 ^= k1;
  h1 += h2;

  h2 = _rotl64(h2,41);

  k2 *= c2; 
  k2  = _rotl64(k2,23);
  k2 *= c1;
  h2 ^= k2;
  h2 += h1;

  h1 = h1*3+0x52dce729;
  h2 = h2*3+0x38495ab5;

  c1 = c1*5+0x7b7d159c;
  c2 = c2*5+0x6bce6396;
}

//----------
// Finalization mix - avalanches all bits to within 0.05% bias

inline uint64_t fmix64 ( uint64_t k )
{
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccd;
  k ^= k >> 33;
  k *= 0xc4ceb9fe1a85ec53;
  k ^= k >> 33;

  return k;
}

void MurmurHash3_x64_128 ( const void * key, const int len, const uint32_t seed, void * out )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 16;

  uint64_t h1 = 0x9368e53c2f6af274 ^ seed;
  uint64_t h2 = 0x586dcd208f7cd3fd ^ seed;

  uint64_t c1 = 0x87c37b91114253d5;
  uint64_t c2 = 0x4cf5ad432745937f;

  //----------
  // body

  const uint64_t * blocks = (const uint64_t *)(data);

  for(int i = 0; i < nblocks; i++)
    {
      uint64_t k1 = getblock(blocks,i*2+0);
      uint64_t k2 = getblock(blocks,i*2+1);

      bmix64(h1,h2,k1,k2,c1,c2);
    }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*16);

  uint64_t k1 = 0;
  uint64_t k2 = 0;

  switch(len & 15)
    {
    case 15: k2 ^= uint64_t(tail[14]) << 48;
    case 14: k2 ^= uint64_t(tail[13]) << 40;
    case 13: k2 ^= uint64_t(tail[12]) << 32;
    case 12: k2 ^= uint64_t(tail[11]) << 24;
    case 11: k2 ^= uint64_t(tail[10]) << 16;
    case 10: k2 ^= uint64_t(tail[ 9]) << 8;
    case  9: k2 ^= uint64_t(tail[ 8]) << 0;

    case  8: k1 ^= uint64_t(tail[ 7]) << 56;
    case  7: k1 ^= uint64_t(tail[ 6]) << 48;
    case  6: k1 ^= uint64_t(tail[ 5]) << 40;
    case  5: k1 ^= uint64_t(tail[ 4]) << 32;
    case  4: k1 ^= uint64_t(tail[ 3]) << 24;
    case  3: k1 ^= uint64_t(tail[ 2]) << 16;
    case  2: k1 ^= uint64_t(tail[ 1]) << 8;
    case  1: k1 ^= uint64_t(tail[ 0]) << 0;
      bmix64(h1,h2,k1,k2,c1,c2);
    };

  //----------
  // finalization

  h2 ^= len;

  h1 += h2;
  h2 += h1;

  h1 = fmix64(h1);
  h2 = fmix64(h2);

  h1 += h2;
  h2 += h1;

  ((uint64_t*)out)[0] = h1;
  ((uint64_t*)out)[1] = h2;
}

//-----------------------------------------------------------------------------
// If we need a smaller hash value, it's faster to just use a portion of the 
// 128-bit hash

void MurmurHash3_x64_32 ( const void * key, int len, uint32_t seed, void * out )
{
  uint32_t temp[4];

  MurmurHash3_x64_128(key,len,seed,temp);

  *(uint32_t*)out = temp[0];
}

//----------

void MurmurHash3_x64_64 ( const void * key, int len, uint32_t seed, void * out )
{
  uint64_t temp[2];

  MurmurHash3_x64_128(key,len,seed,temp);

  *(uint64_t*)out = temp[0];
} 

//-----------------------------------------------------------------------------
