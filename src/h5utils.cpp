#include "h5utils.h"

const char** vec_to_ptr(const std::vector<std::string>& v) {
  const char** ret;
  ret = new const char*[v.size()];
  for (size_t i = 0; i < v.size(); ++i) {
    ret[i] = v[i].c_str();
  }
  return ret;
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
  herr_t status = H5Tset_size (datatype_id, H5T_VARIABLE);
  assert( status >= 0 );
  /* herr_t status = H5Tset_size (datatype_id, max_len); */
  H5Tset_strpad(datatype_id, H5T_STR_NULLTERM);

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
