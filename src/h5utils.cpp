#include "h5utils.h"

// allocate a contiguous block of memory, dependent on the largest string
char* vec_to_ptr(const std::vector<std::string>& v) {
  size_t max_len = 0;
  for (auto& x : v) {
    if (x.size() > max_len) {
      max_len = x.size();
    }
  }
  max_len += 1;

  // allocate a contiguous block of memory
  char *pool = new char[max_len * v.size()];
  char *ptr = pool;

  for (size_t i = 0; i < v.size(); ++i, ptr += max_len) {
    strcpy(ptr, v[i].c_str());
  }

  return pool;
}

const double* vec_to_ptr(const std::vector<double>& v) {
  return v.data();
}

const int* vec_to_ptr(const std::vector<int>& v) {
  return v.data();
}

hid_t get_datatype_id(const std::vector<std::string>& v) {
  size_t max_len = 0;
  for (auto& x : v) {
    if (x.size() > max_len) {
      max_len = x.size();
    }
  }
  max_len += 1;

  hid_t datatype_id = H5Tcopy (H5T_C_S1);
  /* herr_t status = H5Tset_size (datatype_id, H5T_VARIABLE); */
  /* assert( status >= 0 ); */
  herr_t status = H5Tset_size (datatype_id, max_len);
  assert( status >= 0 );

  return datatype_id;
}

hid_t get_datatype_id(const std::vector<double>& v) {
  v.size(); // shutup, compiler
  return H5T_NATIVE_DOUBLE;
}

hid_t get_datatype_id(const std::vector<int>& v) {
  v.size(); // shutup, compiler
  return H5T_NATIVE_INT;
}
