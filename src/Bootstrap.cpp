#include "Bootstrap.h"
#include "PlaintextWriter.h"

EMAlgorithm Bootstrap::run_em() {
    auto counts = mult_.sample();
    EMAlgorithm em(counts, index_, tc_, mean_fls_, opt);

    //em.set_start(em_start);
    em.run(10000, 50, false, false);
    /* em.compute_rho(); */

    return em;
}

BootstrapThreadPool::BootstrapThreadPool(
    size_t n_threads,
    std::vector<size_t> seeds,
    const std::vector<int>& true_counts,
    const KmerIndex& index,
    const MinCollector& tc,
    const std::vector<double>& eff_lens,
    const ProgramOptions& p_opts,
    BootstrapWriter  *bswriter,
    const std::vector<double>& mean_fls
    ) :
  n_threads_(n_threads),
  seeds_(seeds),
  n_complete_(0),
  true_counts_(true_counts),
  index_(index),
  tc_(tc),
  eff_lens_(eff_lens),
  opt_(p_opts),
  writer_(bswriter),
  mean_fls_(mean_fls)
{
  for (size_t i = 0; i < n_threads_; ++i) {
    threads_.push_back( std::thread(BootstrapWorker(*this, i)) );
  }
}

BootstrapThreadPool::~BootstrapThreadPool() {
  for (size_t i = 0; i < n_threads_; ++i) {
    threads_[i].join();
  }
}

void BootstrapWorker::operator() (){
  while (true) {
    size_t cur_seed;
    size_t cur_id;

    // acquire a seed
    {
      std::unique_lock<std::mutex> lock(pool_.seeds_mutex_);

      if (pool_.seeds_.empty()) {
        // no more bootstraps to perform, this thread is done
        return;
      }

      cur_id = pool_.seeds_.size() - 1;
      cur_seed = pool_.seeds_.back();
      pool_.seeds_.pop_back();
      // std::cout << "cur seed from thread (" << thread_id_ << "): " <<
      //   cur_seed <<  " id: " << cur_id << std::endl;
    } // release lock

    Bootstrap bs(pool_.true_counts_,
        pool_.index_,
        pool_.tc_,
        pool_.eff_lens_,
        cur_seed,
        pool_.mean_fls_, pool_.opt_);

    auto res = bs.run_em();

    if (!pool_.opt_.plaintext) {
      std::unique_lock<std::mutex> lock(pool_.write_lock_);
      ++pool_.n_complete_;
      std::cerr << "[bstrp] number of EM bootstraps complete: " << pool_.n_complete_ << "\r";
      pool_.writer_->write_bootstrap(res, cur_id);
      // release write lock
    } else {
      // can write out plaintext in parallel
      plaintext_writer(pool_.opt_.output + "/bs_abundance_" +
          std::to_string(cur_id) + ".tsv",
          pool_.index_.target_names_, res.alpha_,
          pool_.eff_lens_, pool_.index_.target_lens_);
    }
  }
}
