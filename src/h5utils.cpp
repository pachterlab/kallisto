#ifdef USE_HDF5
#include "h5utils.h"

// allocate a contiguous block of memory, dependent on the largest string
char* vec_to_ptr(const std::vector<std::string>& v) {
  size_t max_len = 0;
  for (auto& x : v) {
    if (x.size() > max_len) {
      max_len = x.size();
    }
  }
  // account for 0 terminator
  max_len += 1;
  // allocate a contiguous block of memory
  char *pool = new char[max_len * v.size()];
  memset(pool,0,max_len * v.size());
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

  // account for 0 terminator
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

void read_vector(
    hid_t dataset_id,
    hid_t datatype_id,
    hid_t dataspace_id,
    std::vector<std::string>& out) {

  hsize_t size;
  hsize_t dims[1];

  H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

  size = H5Tget_size(datatype_id);

  char *pool = new char[ size * dims[0] ];
  char *ptr = pool;
  char *buf = new char[size];

  H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, pool);
  out.reserve(dims[0]);

  for (size_t i = 0; i < dims[0]; ++i, ptr += size) {
    memcpy(buf, ptr, size);
    out.push_back( buf );
  }

  delete [] pool;
  delete [] buf;
}

void read_vector(
    hid_t dataset_id,
    hid_t datatype_id,
    hid_t dataspace_id,
    std::vector<int>& out) {

  hsize_t size;
  hsize_t dims[1];

  H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

  size = H5Tget_size(datatype_id);

  int *pool = new int[ dims[0] ];

  H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, pool);
  out.reserve(dims[0]);

  for (size_t i = 0; i < dims[0]; ++i) {
    out.push_back( pool[i] );
  }

  delete [] pool;
}

void read_vector(
    hid_t dataset_id,
    hid_t datatype_id,
    hid_t dataspace_id,
    std::vector<double>& out) {

  hsize_t size;
  hsize_t dims[1];

  H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

  size = H5Tget_size(datatype_id);

  double *pool = new double[ dims[0] ];

  H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, pool);
  out.reserve(dims[0]);

  for (size_t i = 0; i < dims[0]; ++i) {
    out.push_back( pool[i] );
  }

  delete [] pool;
}
#endif // USE_HDF5