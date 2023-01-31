template <class T>
SparseVector<T>::SparseVector(bool init) {
  r = Roaring();
  flag = 0;
  if (init) { // Initialize it for insertion/serialization purposes
    flag = 4;
    new (&v) std::vector<Roaring>*;
    v = new std::vector<Roaring>();
  }
}

template <class T>
SparseVector<T>::~SparseVector() {
  switch (flag) {
  case 1:
    if (arr.a != nullptr) {
      delete[] arr.a;
      arr.a = nullptr;
    }
    if (arr.v != nullptr) {
      delete[] arr.v;
      arr.v = nullptr;
    }
    break;
  case 2:
    if (tinyarr != nullptr) {
      delete[] tinyarr;
      tinyarr = nullptr;
    }
    break;
  case 4:
    if (v != nullptr) {
      v->clear();
      //(*v).~vector();
      delete v;
      v = nullptr;
    }
    break;
  }
}

template <class T>
SparseVector<T>::SparseVector(SparseVector<T> &&arg) {
  flag = arg.flag;
  arg.flag = 0; 
  switch (flag) {
  case 1:
    new (&arr) posinfo;
    arr.v = arg.arr.v;
    arr.a = arg.arr.a;
    arg.arr.v = nullptr;
    arg.arr.a = nullptr;
    break;
  case 2:
    new (&tinyarr) auto(arg.tinyarr);
    tinyarr = arg.tinyarr;
    arg.tinyarr = nullptr;
    break;
  case 3:
    new (&tinybits) auto(arg.tinybits);
    tinybits = arg.tinybits;
    break;
  case 4:
    new (&v) std::vector<Roaring>*;
    v = arg.v;
    arg.v = nullptr;
    break;
  }
  r = std::move(arg.r);
};

template <class T>
SparseVector<T>::SparseVector(const SparseVector<T> &arg) {
  flag = arg.flag;
  if (flag == 1) {
    new (&arr) posinfo;
    arr.v = new uint32_t[arg.cardinality()];
    size_t num_elems = 0;
    for (size_t i = 0; i < arg.cardinality(); i++) {
      arr.v[i] = arg.arr.v[i];
      if ((arg.arr.v[i] | 0x40000000) == arg.arr.v[i]) { // Second MSB is 1
        uint32_t offset = arg.arr.v[i] & ~(0x60000000);
        num_elems += 1; // To account for first element of the array storing the size
        num_elems += arg.arr.a[offset]; // The size of the remaining elements of the array
      }
    }
    if (num_elems != 0) {
      arr.a = new uint32_t[num_elems];
      for (size_t i = 0; i < num_elems; i++) {
        arr.a[i] = arg.arr.a[i];
      } } else {
        arr.a = nullptr;
      }
  } else if (flag == 2) {
    new (&tinyarr) auto(arg.tinyarr);
    tinyarr = new char[arg.cardinality()];
    for (size_t i = 0; i < arg.cardinality(); i++) {
      tinyarr[i] = arg.tinyarr[i];
    }
  } else if (flag == 3) {
    new (&tinybits) auto(arg.tinybits);
    tinybits = arg.tinybits;
  } else if (flag == 4) {
    new (&v) std::vector<Roaring>*;
    v = new std::vector<Roaring>();
    if (arg.v != nullptr) {
      v->reserve(arg.v->size());
      for (size_t i = 0; i < arg.v->size(); i++) {
        v->push_back((*(arg.v))[i]);
      }
    }
  }
  r = arg.r;
};

template <class T>
SparseVector<T>& SparseVector<T>::operator=(SparseVector<T> &&other) {
  if (this == &other) {
    return *this;
  }
  flag = other.flag;
  other.flag = 0;
  switch (flag) {
  case 1:
    new (&arr) posinfo;
    arr.v = other.arr.v;
    arr.a = other.arr.a;
    other.arr.v = nullptr;
    other.arr.a = nullptr;
    break;
  case 2:
    new (&tinyarr) auto(other.tinyarr);
    tinyarr = other.tinyarr;
    other.tinyarr = nullptr;
    break;
  case 3:
    new (&tinybits) auto(other.tinybits);
    tinybits = other.tinybits;
    break;
  case 4:
    new (&v) std::vector<Roaring>*;
    v = other.v;
    other.v = nullptr;
    break;
  }
  r = std::move(other.r);
  return *this;
}

template <class T>
SparseVector<T>& SparseVector<T>::operator=(const SparseVector<T> &other) {
  flag = other.flag;
  switch (flag) {
  case 1:
    new (&arr) posinfo;
    arr.v = new uint32_t[other.cardinality()];
    {
      size_t num_elems = 0;
      for (size_t i = 0; i < other.cardinality(); i++) {
        arr.v[i] = other.arr.v[i];
        if ((other.arr.v[i] | 0x40000000) == other.arr.v[i]) { // Second MSB is 1
          uint32_t offset = other.arr.v[i] & ~(0x60000000);
          num_elems += 1; // To account for first element of the array storing the size
          num_elems += other.arr.a[offset]; // The size of the remaining elements of the array
        }
      }
      if (num_elems != 0) {
        arr.a = new uint32_t[num_elems];
        for (size_t i = 0; i < num_elems; i++) {
          arr.a[i] = other.arr.a[i];
        }
      } else {
        arr.a = nullptr;
      }
    }
    
    break;
  case 2:
    new (&tinyarr) auto(other.tinyarr);
    tinyarr = other.tinyarr;
    tinyarr = new char[other.cardinality()];
    for (size_t i = 0; i < other.cardinality(); i++) {
      tinyarr[i] = other.tinyarr[i];
    }
    break;
  case 3:
    new (&tinybits) auto(other.tinybits);
    tinybits = other.tinybits;
    break;
  case 4:
    new (&v) std::vector<Roaring>*;
    v = new std::vector<Roaring>();
    if (other.v != nullptr) {
      v->reserve(other.v->size());
      for (size_t i = 0; i < other.v->size(); i++) {
        v->push_back((*(other.v))[i]);
      }
    }
    break;
  }
  r = other.r;
  return *this;
}

// Warning:
// Worst case performance is O(N).
// It is recommended to insert items in order of ascending indices.
template <class T>
void SparseVector<T>::insert(size_t i, const T &elem) {
  if (flag == 0) {
    flag = 4;
    new (&v) std::vector<Roaring>*;
    if (v == nullptr) v = new std::vector<Roaring>();
  } else if (flag != 4) {
    throw std::runtime_error("Invalid call to insert() in SparseVector.");
  }
  size_t idx;
  if (r.contains(i)) {
    idx = r.rank(i) - 1;
    (*v)[idx].add(elem);
  } else {
    idx = r.rank(i);
    r.add(i);
    Roaring new_elem;
    new_elem.add(elem);
    if (v->size() == idx) {
      v->push_back(new_elem);
    } else {
      v->emplace(v->begin() + idx, new_elem);
    }
  }
}

template <class T>
void SparseVector<T>::insert(const std::pair<size_t, T> &elem) {
  insert(elem.first, elem.second);
}

// Warning:
// Worst case performance is O(N).
// It is recommended to remove items in order of descending indices.
template <class T>
Roaring SparseVector<T>::remove(size_t i) {
  if (flag != 4) {
    throw std::runtime_error("Invalid call to remove() in SparseVector.");
  }
  Roaring t;
  if (r.contains(i)) {
    size_t idx = r.rank(i) - 1;
    t = (*v)[idx];
    r.remove(i);
    v->erase(v->begin()+idx);
  }
  return t;
}

template <class T>
void SparseVector<T>::clear() {
  if (flag == 4) {
    if (v != nullptr) v->clear();
  }
  if (flag == 2) {
    if (tinyarr != nullptr) {
      delete[] tinyarr;
      tinyarr = nullptr;
    }
  }
  if (flag == 1) {
    if (arr.a != nullptr) {
      delete[] arr.a;
      arr.a = nullptr;
    }
    if (arr.v != nullptr) {
      delete[] arr.v;
      arr.v = nullptr;
    }
  }
  r = Roaring();
  flag = 0;
}

template <class T>
const Roaring& SparseVector<T>::getIndices() const {
  return r;
}

template <class T>
void SparseVector<T>::getElements(std::vector<std::pair<uint32_t, Roaring> > &elems) const {
  if (flag != 4) {
    throw std::runtime_error("Invalid call to getElements() in SparseVector.");
  }
  if (r.isEmpty()) return;
  elems.reserve(r.cardinality());
  uint32_t i = 0;
  for (const auto &idx : r) {
    elems.push_back(std::pair<uint32_t, Roaring>(idx, (*v)[i]));
    ++i;
  }
}

template <class T>
Roaring SparseVector<T>::get(size_t i, bool getOne) const {
  if (!(flag == 1)) {
    throw std::runtime_error("Invalid call to get() in SparseVector.");
  }
  if (flag == 1 && r.contains(i)) {
    Roaring x;
    i = r.rank(i)-1;
    if ((arr.v[i] | 0x40000000) == arr.v[i]) { // Second MSB is 1
        uint32_t offset = arr.v[i] & ~(0x60000000); // Mask out second and third MSB
        uint32_t arr_size = arr.a[offset];
        for (size_t n = 0; n < arr_size; n++) {
          x.add(arr.a[offset+n+1]);
          if (getOne) {
            break;
          }
        }
    } else { // Second MSB is 0
        x.add(arr.v[i] & ~(0x60000000)); // Mask out second and third MSB
    }
    return x;
  }
  throw std::invalid_argument("Index not present in SparseVector.");
}

template <class T>
bool SparseVector<T>::contains(size_t i) const {
  return r.contains(i);
}

template <class T>
bool SparseVector<T>::isEmpty() const {
  return r.isEmpty();
}

template <class T>
char SparseVector<T>::operator[] (size_t i) {
  if (r.contains(i)) {
    if (flag == 2) {
      return (tinyarr[r.rank(i) - 1]);
    } else if (flag == 3) {
      return (bool)((tinybits & ((uint64_t)1 << (uint64_t)(r.rank(i)-1))) != 0);
    } else if (flag == 4) {
      return ((*v)[r.rank(i) - 1].minimum() & 0x7FFFFFFF) == (*v)[r.rank(i) - 1].minimum();
    } else if (flag == 1) {
      i = r.rank(i)-1;
      if ((arr.v[i] | 0x20000000) == arr.v[i]) { // Third MSB is set to 1 (aka ambiguous strand)
        return 2; 
      }
      if ((arr.v[i] | 0x40000000) != arr.v[i]) { // Second MSB not set to 1
        return (arr.v[i] & 0x7FFFFFFF) == arr.v[i];
      } else { // Second MSB set to 1
        uint32_t offset = arr.v[i] & ~(0x60000000); // Mask out second and third MSB
        return (arr.a[offset+1] & 0x7FFFFFFF) == arr.a[offset+1];
      }
    } else {
      throw std::runtime_error("Invalid call to operator[] in SparseVector.");
    }
  }
  throw std::invalid_argument("Index not present in SparseVector.");
}

template <class T>
const char SparseVector<T>::operator[] (size_t i) const {
  if (r.contains(i)) {
    if (flag == 2) {
      return (tinyarr[r.rank(i) - 1]);
    } else if (flag == 3) {
      return (bool)((tinybits & ((uint64_t)1 << (uint64_t)(r.rank(i)-1))) != 0);
    } else if (flag == 4) {
      return ((*v)[r.rank(i) - 1].minimum() & 0x7FFFFFFF) == (*v)[r.rank(i) - 1].minimum();
    } else if (flag == 1) {
      i = r.rank(i)-1;
      if ((arr.v[i] | 0x20000000) == arr.v[i]) { // Third MSB is set to 1 (aka ambiguous strand)
        return 2;
      }
      if ((arr.v[i] | 0x40000000) != arr.v[i]) { // Second MSB not set to 1
        return (arr.v[i] & 0x7FFFFFFF) == arr.v[i];
      } else { // Second MSB set to 1
        uint32_t offset = arr.v[i] & ~(0x60000000); // Mask out second and third MSB
        return (arr.a[offset+1] & 0x7FFFFFFF) == arr.a[offset+1];
      }
    } else {
      throw std::runtime_error("Invalid call to operator[] in SparseVector.");
    }
  }
  throw std::invalid_argument("Index not present in SparseVector.");
}


template <class T>
bool SparseVector<T>::operator==(const SparseVector<T>& other) const {
  return r == other.r; 
}

template <class T>
void SparseVector<T>::serialize(std::ostream& out) const {
  if (flag != 4) {
    throw std::runtime_error("Invalid call to serialize() in SparseVector.");
  }
  // Write Roaring
  size_t tmp_size = r.getSizeInBytes(false);
  char *buffer = new char[tmp_size];
  r.write(buffer, false);
  out.write((char *)&tmp_size, sizeof(tmp_size));
  out.write(buffer, tmp_size);
  delete[] buffer;
  // Write vector
  tmp_size = 0;
  if (v != nullptr) tmp_size = v->size();
  out.write((char *)&tmp_size, sizeof(tmp_size));
  if (v == nullptr) return;
  for (auto p : *v) {
    p.runOptimize();
    tmp_size = p.getSizeInBytes(false);
    buffer = new char[tmp_size];
    p.write(buffer, false);
    out.write((char *)&tmp_size, sizeof(tmp_size));
    out.write(buffer, tmp_size);
    delete[] buffer;
  }
}

template <class T>
void SparseVector<T>::deserialize(std::istream& in, bool small) {
  if (flag == 4) {
    throw std::runtime_error("Invalid call to deserialize() in SparseVector.");
  }
  size_t tmp_size;
  in.read((char *)&tmp_size, sizeof(tmp_size)); // Roaring size
  char* buffer = new char[tmp_size];
  in.read(buffer, tmp_size);
  r = Roaring::read(buffer, false);
  delete[] buffer;
  size_t v_size;
  in.read((char *)&v_size, sizeof(v_size)); // Number of elements (aka number of transcripts in set)
  char* tinyarr_ = nullptr;
  uint64_t tinybits_ = 0;
  if (small) {
    tinyarr_ = new char[v_size];
    if (v_size > 64) { // large enough so might as well use tinyarr instead of tinybits
      flag = 2;
      new (&tinyarr) char*;
      tinyarr = tinyarr_; // copy pointer
    }
  } else { // Store everything: strand+position
    if (flag == 0) {
      flag = 1;
      new (&arr) posinfo;
      arr.a = nullptr;
      arr.v = nullptr;
    }
    arr.v = new uint32_t[v_size];
  }
  assert(r.cardinality() == v_size);
  size_t offset = 0; // offset (only for flag=1)
  std::vector<uint32_t> arr_a_vec; // Temporary vector to store contents that will be transferred to arr.a (only for flag=1)
  for (size_t i = 0; i < v_size; ++i) {
    in.read((char *)&tmp_size, sizeof(tmp_size)); // Roaring size
    buffer = new char[tmp_size];
    in.read(buffer, tmp_size);
    if (!small) {
      Roaring x = Roaring::read(buffer, false);
      if (x.cardinality() == 1) { // Single item (pos/strand) in this transcript's set
        arr.v[i] = (uint32_t)x.minimum();
      } else { // Multiple items in this transcript's set
        arr.v[i] = (uint32_t)(offset); // Set offset
        arr.v[i] |= 0x40000000; // Set second MSB to 1 (denoting multiple items in transcript's set)
        if (((x.minimum() & 0x7FFFFFFF) == x.minimum()) != ((x.maximum() & 0x7FFFFFFF) == x.maximum())) { // ambiguous strand
          arr.v[i] |= 0x20000000; // Set third MSB to 1 (denoting strand ambiguity)
        }
        arr_a_vec.push_back(x.cardinality());
        offset++; // Since first element is the size (number of elements)
        for (uint32_t p : x) {
          offset++;
          arr_a_vec.push_back(p);
        }
      }
    } else {
      Roaring x = Roaring::read(buffer, false);
      if (flag != 2) {
        if ((x.minimum() & 0x7FFFFFFF) == x.minimum()) {
          tinybits_ |= (uint64_t)((uint64_t)1 << i);
        }
        tinyarr_[i] = (char)((x.minimum() & 0x7FFFFFFF) == x.minimum());
        if (((x.minimum() & 0x7FFFFFFF) == x.minimum()) != ((x.maximum() & 0x7FFFFFFF) == x.maximum())) { // ambiguous strand; don't store in compressed format
          flag = 2;
          new (&tinyarr) char*;
          tinyarr = tinyarr_; // copy pointer
        }
      }
      if (flag == 2) {
        tinyarr_[i] = (char)(((x.minimum() & 0x7FFFFFFF) == x.minimum()));
        if (((x.minimum() & 0x7FFFFFFF) == x.minimum()) != ((x.maximum() & 0x7FFFFFFF) == x.maximum())) tinyarr_[i] = 2; // ambiguous
      }
    }
    delete[] buffer;
  }
  // For small case with comppressed 64-bit representation, set flag and do transfer
  if (small && flag != 2) {
		      flag = 3;
		      new (&tinybits) uint64_t;
		      tinybits = tinybits_;
		      delete[] tinyarr_; // free memory
		      tinyarr_ = nullptr;
  }
  // Transfer from arr_a_vec to arr.a (for flag=1)
  if (!small && !arr_a_vec.empty()) {
    arr.a = new uint32_t[arr_a_vec.size()];
    size_t i = 0;
    for (auto x : arr_a_vec) {
      arr.a[i++] = x;
    }
  }
}

template <class T>
void SparseVector<T>::runOptimize() {
  r.runOptimize(); // Note: Roarings in v not optimized here
}

template <class T>
size_t SparseVector<T>::cardinality() const {
  return r.cardinality();
}
