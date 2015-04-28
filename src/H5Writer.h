#ifndef KALLISTO_H5WRITER_H
#define KALLISTO_H5WRITER_H

#include "EMAlgorithm.h"

#include "h5utils.h"

class H5Writer {
  public:
    H5Writer(const std::string& fname, int num_bootstrap, uint compression);
    ~H5Writer();

    void write_main(const EMAlgorithm& em,
        const std::vector<std::string>& targ_ids,
        const std::vector<int>& lengths);

    void write_bootstrap(const EMAlgorithm& em, int bs_id);


  private:
      int num_bootstrap_;
      uint compression_;

      hid_t file_id_;
      hid_t root_;
      hid_t aux_;
      hid_t bs_;
};

class H5Reader {
  public:
    H5Reader(const std::string& h5_fname, const std::string& out_dir);
    ~H5Reader();

    void rw_bootstraps();

  private:
    void rw_main();

    /* EMAlgorithm main_em_; */
    std::vector<int>& lengths_;

    hid_t file_id_;
    hid_t root_;
    hid_t aux_;
    hid_t bs_;

    int n_bs_;
}
#endif // KALLISTO_H5WRITER_H
