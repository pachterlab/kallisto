#include "H5Writer.h"

H5Writer::H5Writer(const std::string& fname, int num_bootstrap, uint compression) :
  num_bootstrap_(num_bootstrap)
{
  compression_ = compression;
  file_id_ = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  root_ = H5Gopen(file_id_, "/", H5P_DEFAULT);
  aux_ = H5Gcreate(file_id_, "/aux", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (num_bootstrap_ > 0) {
    bs_ = H5Gcreate(file_id_, "/bootstrap", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  std::vector<int> n_bs;
  n_bs.push_back(num_bootstrap);
  vector_to_h5(n_bs, aux_, "num_bootstrap", false, compression_);
}

H5Writer::~H5Writer() {
  if (num_bootstrap_ > 0) {
    H5Gclose(bs_);
  }
  H5Gclose(aux_);
  H5Gclose(root_);
  H5Fclose(file_id_);
}

void H5Writer::write_main(const EMAlgorithm& em,
    const std::vector<std::string>& targ_ids,
    const std::vector<int>& lengths) {
  vector_to_h5(em.alpha_, root_, "est_counts", false, compression_);

  vector_to_h5(targ_ids, aux_, "ids", true, compression_);
  vector_to_h5(em.eff_lens_, aux_, "eff_lengths", false, compression_);
  vector_to_h5(lengths, aux_, "lengths", false, compression_);
}

void H5Writer::write_bootstrap(const EMAlgorithm& em, int bs_id) {
  std::string bs_id_str("bs" + std::to_string( bs_id ));
  vector_to_h5(em.alpha_, bs_, bs_id_str.c_str(), false, compression_);
}

/**********************************************************************/

H5Reader::H5Reader(const std::string& h5_fname, const std::string& out_dir) {
  file_id_ = H5Fopen(h5_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  root_ = H5Gopen(file_id_, "/", H5P_DEFAULT);
  aux_ = H5Gopen(file_id_, "/aux", H5P_DEFAULT);

  std::vector<int> n_bs_vec;
  read_dataset(aux_, "num_bootstrap", n_bs_vec);
  n_bs_ = n_bs_vec[0];

  // TODO: create output directory
}

H5Reader::~H5Reader() {
  H5Gclose(aux_);
  H5Gclose(root_);
  H5Fclose(file_id_);
}
