#ifndef SPARSEVECTOR_HPP
#define SPARSEVECTOR_HPP

#include <vector>
#include <iostream>
#include "roaring.hh"

template <class T>
class SparseVector { // This class is stored in a BlockArray and associates transcript IDs with positions along unitig
public:
  SparseVector(bool init = false);
  SparseVector(SparseVector<T> &&arg);
  SparseVector(const SparseVector<T> &arg);
  ~SparseVector();
  
  void insert(size_t i, const T &elem); // Inserts transcript ID i along with its position info: elem (note: a single transcript can be associated with multiple positions along the unitig)
  void insert(const std::pair<size_t, T> &elem);
  Roaring remove(size_t i);
  void clear();
  
  const Roaring& getIndices() const; // Get all transcript IDs stored in this object	
  // Populates elems with (index, element) tuples
  void getElements(std::vector<std::pair<uint32_t, Roaring> > &elems) const;
  Roaring get(size_t i, bool getOne=false) const;
  bool contains(size_t i) const; // Whether the transcript with id i exists in this object
  bool isEmpty() const;
  
  SparseVector<T>& operator=(SparseVector<T> &&other);
  SparseVector<T>& operator=(const SparseVector<T> &other);
  // checks if r is equal between two objects (only r is considered, nothing else):
  bool operator==(const SparseVector<T>& other) const;
  char operator[] (size_t i); // Returns strandedness: (pos & 0x7FFFFFFF) == pos, for transcript id: i; could also return 2 if ambiguous (e.g. for tx i, a + part and a - part exist in this block)
  const char operator[] (size_t i) const;
  
  // Serialization/Deserialization
  void serialize(std::ostream& out) const;
  void deserialize(std::istream& in, bool small=true);
  
  void runOptimize();
  size_t cardinality() const;
  
private:
  Roaring r; // Set of transcripts in this data structure {tx A, tx B, tx C, ...}
  uint8_t flag; // How the storage should work see below:
  // flag=0 means uninitialized / low-memory / no member in the union activated
  // flag=1 means store strand+positional info in posinfo struct (high memory)
  // flag=2 means store strand info in a char array and don't store positional info
  // flag=3 is like flag=2 except store strand info is stored as individual bits in a 64-bit integer (lowest memory but only works if cardinality of set r <= 64)
  // flag=4 means strand+positional info grows dynamically (we can add stuff into it like a vector; this is only for insertion/serialization purposes, not for querying purposes)
  
  // TODO: See if we can just make the following one array (instead of two) in order to optimize alignment
  struct posinfo { // 16 byte data structure containing two arrays (see below):
    uint32_t* v; // For each element (num elements = r.cardinality()), if second most significant bit (MSB) is 0 [since first MSB is the strand], the transcript has a single position+strand which is contained in the remaining 30 bits (note: the third MSB is excluded)
    uint32_t* a; // If an element in v's second MSB is 1, the remaining bits in v (after the third MSB) contains the offset w.r.t. the pointer a. a+o is an array of the the multiple positions+strands with the first element being the the number of elements in the a+o array (aka the positions/strands span from (a+o)[1] to (a+o)[1]+(a+o)[0], inclusive, where o=offset*sizeof(uint32_t)).
  }; // What about the third MSB in v's element? It's 1 if the transcript is strand-ambiguous (i.e. both + and - exist in that transcript's array a), 0 otherwise
  
  union {
    posinfo arr; // activated when flag=1
    char* tinyarr; // activated when flag=2
    uint64_t tinybits; // activated when flag=3
    std::vector<Roaring>* v; // activated when flag=4 (pointer to vector to avoid alignment issues since vector is 24 bytes)
  };
};

#include "SparseVector.tcc"

#endif
