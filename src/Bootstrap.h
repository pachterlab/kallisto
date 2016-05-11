#ifndef KALLISTO_BOOTSTRAP_H
#define KALLISTO_BOOTSTRAP_H

#include <mutex>
#include <thread>

#include "KmerIndex.h"
#include "MinCollector.h"
#include "weights.h"
#include "EMAlgorithm.h"
#include "Multinomial.hpp"
#include "H5Writer.h"

class Bootstrap {
    // needs:
    // - "true" counts
    // - ecmap
    // - target_names
    // - eff_lens
public:
  Bootstrap(const std::vector<int>& true_counts,
            const KmerIndex& index,
            const MinCollector& tc,
            const std::vector<double>& eff_lens,
            size_t seed,
            const std::vector<double>& mean_fls,
            const ProgramOptions& opt) :
    index_(index),
    tc_(tc),
    eff_lens_(eff_lens),
    seed_(seed),
    mult_(true_counts, seed_),
    mean_fls_(mean_fls),
    opt(opt)
    {}

  // EM Algorithm generates a sample from the Multinomial, then returns
  // an "EMAlgorithm" that has already run the EM as well as compute the
        // rho values
  EMAlgorithm run_em();

private:
  const KmerIndex& index_;
  const MinCollector& tc_;
  const std::vector<double>& eff_lens_;
  size_t seed_;
  Multinomial mult_;
  const std::vector<double>& mean_fls_;
  const ProgramOptions& opt;
};

class BootstrapThreadPool {
  friend class BootstrapWorker;

  public:
    BootstrapThreadPool(
        size_t n_threads,
        std::vector<size_t> seeds,
        const std::vector<int>& true_counts,
        const KmerIndex& index,
        const MinCollector& tc,
        const std::vector<double>& eff_lens,
        const ProgramOptions& p_opts,
        H5Writer& h5writer,
        const std::vector<double>& mean_fls
        );

    size_t num_threads() {return n_threads_;}

    ~BootstrapThreadPool();
  private:
    std::vector<size_t> seeds_;
    size_t n_threads_;

    std::vector<std::thread> threads_;
    std::mutex seeds_mutex_;
    std::mutex write_lock_;

    size_t n_complete_;

    // things to run bootstrap
    const std::vector<int> true_counts_;
    const KmerIndex& index_;
    const MinCollector& tc_;
    const std::vector<double>& eff_lens_;
    const ProgramOptions& opt_;
    H5Writer& writer_;
    const std::vector<double>& mean_fls_;
};

class BootstrapWorker {
  public:
    BootstrapWorker(BootstrapThreadPool& pool, size_t thread_id) :
      pool_(pool),
      thread_id_(thread_id)
    {}

  void operator()();

  private:
    BootstrapThreadPool& pool_;
    size_t thread_id_;
};

#endif
