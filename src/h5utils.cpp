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
  hid_t datatype_id = H5Tcopy (H5T_C_S1);
  H5Tset_size (datatype_id, H5T_VARIABLE);
  v.size(); // shutup, compiler
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
