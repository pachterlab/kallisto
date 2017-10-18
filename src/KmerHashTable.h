#ifndef KALLISTO_KMERHASHTABLE_H
#define KALLISTO_KMERHASHTABLE_H

#include <utility>
#include <string>
#include <iterator>

/*#include <iostream> // debug
	using namespace std;*/


template<typename T, typename Hash = KmerHash>
struct KmerHashTable {
  using value_type = std::pair<Kmer, T>;
  using key_type = Kmer;
  using mapped_type = T;

  Hash hasher;
  value_type *table;
  size_t size_, pop;
  value_type empty;
  value_type deleted;


// ---- iterator ----

  template<bool is_const_iterator = true>
  class iterator_ : public std::iterator<std::forward_iterator_tag, value_type> {
   public:

    typedef typename std::conditional<is_const_iterator, const KmerHashTable *, KmerHashTable *>::type DataStructurePointerType;
    typedef typename std::conditional<is_const_iterator, const value_type&, value_type&>::type ValueReferenceType;
    typedef typename std::conditional<is_const_iterator, const value_type *, value_type *>::type ValuePointerType;


    DataStructurePointerType ht;
    size_t h;

    iterator_() : ht(nullptr), h(0) {}
    iterator_(DataStructurePointerType ht_) : ht(ht_), h(ht_->size_) {}
    iterator_(DataStructurePointerType ht_, size_t h_) :  ht(ht_), h(h_) {}
    iterator_(const iterator_<false>& o) : ht(o.ht), h(o.h) {}
    iterator_& operator=(const iterator_& o) {ht=o.ht; h=o.h; return *this;}

    ValueReferenceType operator*() const {return ht->table[h];}
    ValuePointerType operator->() const {return &(ht->table[h]);}

    void find_first() {
      h = 0;
      if (ht->table != nullptr && ht->size_>0) {
        Kmer& km = ht->table[h].first;
        if (km == ht->empty.first || km == ht->deleted.first) {
          operator++();
        }
      }
    }

    iterator_ operator++(int) {
      const iterator_ old(*this);
      ++(*this);
      return old;
    }

    iterator_& operator++() {
      if (h == ht->size_) {
        return *this;
      }
      ++h;
      for (; h < ht->size_; ++h) {
        Kmer& km = ht->table[h].first;
        if (km != ht->empty.first && km != ht->deleted.first) {
          break;
        }
      }
      return *this;
    }
    bool operator==(const iterator_ &o) const {return (ht->table == o.ht->table) && (h == o.h);}
    bool operator!=(const iterator_ &o) const {return !(this->operator==(o));}
    friend class iterator_<true>;
  };

  typedef iterator_<true> const_iterator;
  typedef iterator_<false> iterator;


  // --- hash table


  KmerHashTable(const Hash& h = Hash() ) : hasher(h), table(nullptr), size_(0), pop(0) {
    empty.first.set_empty();
    deleted.first.set_deleted();
    init_table(1024);
  }

  KmerHashTable(size_t sz, const Hash& h = Hash() ) : hasher(h), table(nullptr), size_(0), pop(0) {
    empty.first.set_empty();
    deleted.first.set_deleted();
    init_table((size_t) (1.2*sz));
  }

  ~KmerHashTable() {
    clear_table();
  }

  void clear_table() {
    if (table != nullptr) {
      delete[] table;
      table = nullptr;
    }
    size_ = 0;
    pop  = 0;
  }

  size_t size() const {
    return pop;
  }

  void clear() {
    std::fill(table, table+size_, empty);
    pop = 0;
  }

  void init_table(size_t sz) {
    clear_table();
    size_ = rndup(sz);
    //cerr << "init table of size " << size_ << endl;
    table = new value_type[size_];
    std::fill(table, table+size_, empty);
  }

  iterator find(const Kmer& key) {
    size_t h = hasher(key) & (size_-1);

    for (;; h =  (h+1!=size_ ? h+1 : 0)) {
      if (table[h].first == empty.first) {
        // empty slot, not in table
        return iterator(this);
      } else if (table[h].first == key) {
        // same key, found
        return iterator(this, h);
      } // if it is deleted, we still have to continue
    }
  }

  const_iterator find(const Kmer& key) const {

    size_t h = hasher(key) & (size_-1);

    for (;; h =  (h+1!=size_ ? h+1 : 0)) {
      if (table[h].first == empty.first) {
        // empty slot, not in table
        return const_iterator(this);
      } else if (table[h].first == key) {
        // same key, found
        return const_iterator(this, h);
      }
    }
  }


  iterator erase(const_iterator pos) {
    if (pos == this->end()) {
      return this->end();
    }
    size_t h = pos.h;
    table[h] = deleted;
    --pop;
    return ++iterator(this, h); // return pointer to next element
  }

  size_t erase(const Kmer& km) {
    const_iterator pos = find(km);
    size_t oldpop = pop;
    if (pos != this->end()) {
      erase(pos);
    }
    return oldpop-pop;
  }

  std::pair<iterator,bool> insert(const value_type& val) {
    //cerr << "inserting " << val.first.toString() << " = " << val.second << endl;
    if ((pop + (pop>>4))> size_) { // if more than 80% full
      //cerr << "-- triggered resize--" << endl;
      reserve(2*size_, false);
    }

    size_t h = hasher(val.first) & (size_-1);
    //cerr << " hash value = " << h << endl;
    for (;; h = (h+1!=size_ ? h+1 : 0)) {
      //cerr << "  lookup at " << h << endl;
      if (table[h].first == empty.first || table[h].first == deleted.first) {
        //cerr << "   found empty slot" << endl;
        // empty slot, insert here
        table[h] = val;
        ++pop; // new table
        return {iterator(this, h), true};
      } else if (table[h].first == val.first) {
        // same key, update value
        //cerr << "   found key already here " << table[h].first.toString() << " = " << table[h].second <<  endl;
        return {iterator(this, h), false};
      }
    }

  }

  void reserve(size_t sz, bool expand = true) {
    if (expand) {
      sz = sz + (sz>>4); // what we actually need for insert
    }
    if (sz <= size_) {
      return;
    }

    value_type *old_table = table;
    size_t old_size_ = size_;


    size_ = rndup(sz);
    pop = 0;

    table = new value_type[size_];
    std::fill(table, table+size_, empty);
    for (size_t i = 0; i < old_size_; i++) {
      if (old_table[i].first != empty.first && old_table[i].first != deleted.first) {
        insert(old_table[i]);
      }
    }
    delete[] old_table;
    old_table = nullptr;

  }

  size_t rndup(size_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
  }

  iterator begin() {
    iterator it(this);
    it.find_first();
    return it;
  }

  const_iterator begin() const {
    const_iterator it(this);
    it.find_first();
    return it;
  }

  iterator end() {
    return iterator(this);
  }

  const_iterator end() const {
    return const_iterator(this);
  }





};

#endif // KALLISTO_KMERHASHTABLE_H
