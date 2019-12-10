#ifdef USE_HDF5
#ifndef KALLISTO_H5WRITER_H
#define KALLISTO_H5WRITER_H

#include "EMAlgorithm.h"

#include "h5utils.h"
#include "PlaintextWriter.h"
#include "Bootstrap.h"

class H5Writer : public BootstrapWriter {
  public:
    H5Writer() : primed_(false) {}
    ~H5Writer();

    virtual void init(const std::string& fname, int num_bootstrap, int num_processed,
      const std::vector<int>& fld, const std::vector<int>& preBias, const std::vector<double>& postBias, uint compression, size_t index_version,
      const std::string& shell_call, const std::string& start_time);

    virtual void write_main(const EMAlgorithm& em,
        const std::vector<std::string>& targ_ids,
        const std::vector<int>& lengths);

    virtual void write_bootstrap(const EMAlgorithm& em, int bs_id);

  private:
    bool primed_;

    int num_bootstrap_;
    uint compression_;

    hid_t file_id_;
    hid_t root_;
    hid_t aux_;
    hid_t bs_;
};

class H5Converter {
  public:
    // assumes 'out_dir' is already
    H5Converter(const std::string& h5_fname, const std::string& out_dir);
    ~H5Converter();

    void write_aux();
    void convert();

  private:
    void rw_from_counts(hid_t group_id, const std::string& count_name,
        const std::string& out_fname);

    std::string out_dir_;

    // run info
    std::string kal_version_;
    size_t idx_version_;
    std::string start_time_;
    std::string call_;

    // auxilary
    std::vector<std::string> targ_ids_;
    std::vector<int> lengths_;
    std::vector<double> eff_lengths_;

    // buffers used for every single read/write
    std::vector<int> kal_id_;
    std::vector<double> alpha_buf_;
    std::vector<double> tpm_buf_;

    hid_t file_id_;
    hid_t root_;
    hid_t aux_;
    hid_t bs_;

    int n_bs_;
    int n_proc_;
    int n_paln_;
    size_t n_targs_;
};

#endif // KALLISTO_H5WRITER_H
#endif // USE_HDF5