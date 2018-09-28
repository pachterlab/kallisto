#ifndef KALLISTO_H5_UTILS
#define KALLISTO_H5_UTILS

#include <assert.h>

#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include <iostream>

#include "hdf5.h"

#ifdef _WIN64
typedef unsigned int uint;
#endif

// begin: writing utils
// XXX: remember to cleanup result!
char* vec_to_ptr(const std::vector<std::string>& v);

const double* vec_to_ptr(const std::vector<double>& v);

const int* vec_to_ptr(const std::vector<int>& v);

hid_t get_datatype_id(const std::vector<std::string>& v);

hid_t get_datatype_id(const std::vector<double>& v);

hid_t get_datatype_id(const std::vector<int>& v);

// str_vec: a vector of string to be written out
// group_id: a group_id which has already been opened
// dataset_name: the to write out to
// release_type: if 'true', release the datatype and ptr
// compression_level: the level of compression (6 seems reasonable)
//
// return: the status of H5Dwrite (last H5 operation)
template <typename T>
herr_t vector_to_h5(
    const std::vector<T>& str_vec,
    hid_t group_id,
    const std::string& dataset_name,
    bool release_type,
    uint compression_level = 6
    ) {
  herr_t status;

  hsize_t dims[1] = {str_vec.size()};

  // create the propery which allows for compression
  hid_t prop_id = H5Pcreate(H5P_DATASET_CREATE);
  // chunk size is same size as vector
  status = H5Pset_chunk(prop_id, 1, dims);
  assert( status >= 0 );

  status = H5Pset_deflate(prop_id, compression_level);
  assert( status >= 0 );

  // create the data type
  hid_t datatype_id = get_datatype_id(str_vec);

  // create the dataspace
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

  // create the dataset
  hid_t dataset_id = H5Dcreate(group_id, dataset_name.c_str(), datatype_id,
      dataspace_id, H5P_DEFAULT, prop_id, H5P_DEFAULT);

  // get the ptrs from the string and write out
  auto ptr = vec_to_ptr(str_vec);
  status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
      ptr);
  assert( status >= 0 );

  status = H5Pclose(prop_id);
  assert( status >= 0 );
  status = H5Dclose(dataset_id);
  assert( status >= 0 );
  status = H5Sclose(dataspace_id);
  assert( status >= 0 );
  if (release_type) {
    status = H5Tclose(datatype_id);
    assert( status >= 0 );
    delete [] ptr;
  }

  return status;
}

// end: writing utils

// begin: reading utils

// pre: out is an empty vector
// post: out contains the data read in from HDF5
void read_vector(hid_t dataset_id, hid_t datatype_id, hid_t dataspace_id,
    std::vector<std::string>& out);

void read_vector(hid_t dataset_id, hid_t datatype_id, hid_t dataspace_id,
    std::vector<int>& out);

void read_vector(hid_t dataset_id, hid_t datatype_id, hid_t dataspace_id,
    std::vector<double>& out);

template <typename T>
void read_dataset(hid_t group_id,
    const std::string& dset_name,
    std::vector<T>& out) {

  hid_t dataset_id;
  hid_t datatype_id;
  hid_t dataspace_id;
  hid_t prop_id;

  H5Z_filter_t filter_type;

  unsigned int flags;
  unsigned int filter_info;

  size_t nelem = 0;

  herr_t status;

  dataset_id = H5Dopen(group_id, dset_name.c_str(), H5P_DEFAULT);
  prop_id = H5Dget_create_plist(dataset_id);
  filter_type = H5Pget_filter(prop_id, 0, &flags, &nelem, NULL, 0, NULL,
                &filter_info);

  datatype_id = H5Dget_type(dataset_id);

  dataspace_id = H5Dget_space(dataset_id);
  read_vector(dataset_id, datatype_id, dataspace_id, out);

  status = H5Pclose(prop_id);
  assert(status >= 0);
  status = H5Dclose(dataset_id);
  assert(status >= 0);
  status = H5Sclose(dataspace_id);
  assert(status >= 0);
}

// end: reading utils

#endif // KALLISTO_H5_UTILS
