#ifndef KALLISTO_BOOTSTRAP_H
#define KALLISTO_BOOTSTRAP_H

#include <mutex>
#include <thread>

#include "KmerIndex.h"
#include "MinCollector.h"
#include "weights.h"
#include "EMAlgorithm.h"
#include "Multinomial.hpp"

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
            double mean,
            size_t seed) :
    index_(index),
    tc_(tc),
    eff_lens_(eff_lens),
    mean_fl(mean),
    seed_(seed),
    mult_(true_counts, seed_)
    {}

  // EM Algorithm generates a sample from the Multinomial, then returns
  // an "EMAlgorithm" that has already run the EM as well as compute the
        // rho values
  EMAlgorithm run_em(const EMAlgorithm& em_start);

private:
  const KmerIndex& index_;
  const MinCollector& tc_;
  const std::vector<double>& eff_lens_;
  double mean_fl;
  size_t seed_;
  Multinomial mult_;
};

class BootstrapThreadPool {
  friend class BootstrapWorker;

  public:
    BootstrapThreadPool(
        size_t n_threads,
        std::vector<size_t> seeds
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
