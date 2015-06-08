#include "Bootstrap.h"
// #include "weights.h"
// #include "EMAlgorithm.h"

EMAlgorithm Bootstrap::run_em(const EMAlgorithm& em_start) {
    auto counts = mult_.sample();
    EMAlgorithm em(counts, index_, tc_, mean_fl);

    //em.set_start(em_start);
    em.run(10000, 50, false, false);
    /* em.compute_rho(); */

    return em;
}

BootstrapThreadPool::BootstrapThreadPool(
    size_t n_threads,
    std::vector<size_t> seeds) :
  n_threads_(n_threads),
  seeds_(seeds),
  n_complete_(0)
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
      std::cout << "cur seed from thread (" << thread_id_ << "): " <<
        cur_seed <<  " id: " << cur_id << std::endl;
    } // release lock


    // write out data
    {
      std::unique_lock<std::mutex> lock(pool_.write_lock_);
      ++pool_.n_complete_;
      //std::cerr << "[bstrp] bootstraps complete: " << pool_.n_complete_ << "\r";

    } // release write lock
  }

}
