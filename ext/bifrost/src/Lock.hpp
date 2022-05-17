#ifndef BIFROST_SPINLOCK_HPP
#define BIFROST_SPINLOCK_HPP

#include <vector>

#include "rw_spin_lock.h"

class LockGraph : public SpinLockRW {

    public:

        LockGraph(const size_t nb_locks_min) : nb_locks(rndup(nb_locks_min)), mask_nb_locks(nb_locks-1) {

            spinlocks_unitigs = std::vector<SpinLock>(nb_locks);
        }

        BFG_INLINE void clear() {

            spinlocks_unitigs.clear();
        }

        BFG_INLINE void lock_unitig(const size_t unitig_id){

             spinlocks_unitigs[unitig_id & mask_nb_locks].acquire();
        }

        BFG_INLINE void unlock_unitig(const size_t unitig_id){

            spinlocks_unitigs[unitig_id & mask_nb_locks].release();
        }

    private :

        const size_t nb_locks;
        const size_t mask_nb_locks;

        std::vector<SpinLock> spinlocks_unitigs;
};

/*class LockGraph : public SpinLockRW_MCS {

    public:

        LockGraph(const size_t nb_threads) :    SpinLockRW_MCS(nb_threads), mask_nb_locks(rndup(nb_threads * nb_locks_per_thread) - 1),
                                                spinlocks_unitigs(mask_nb_locks + 1) {}

        BFG_INLINE void clear() {

            SpinLockRW_MCS::clear();
            spinlocks_unitigs.clear();
        }

        BFG_INLINE void lock_unitig(const size_t unitig_id){

             spinlocks_unitigs[unitig_id & mask_nb_locks].acquire();
        }

        BFG_INLINE void unlock_unitig(const size_t unitig_id){

            spinlocks_unitigs[unitig_id & mask_nb_locks].release();
        }

    private :

        const size_t nb_locks_per_thread = 1024;
        const size_t mask_nb_locks;

        std::vector<SpinLock> spinlocks_unitigs;
};*/

/*class LockGraph : public Hybrid_SpinLockRW_MCS<> {

    public:

        LockGraph(const size_t nb_threads) :    Hybrid_SpinLockRW_MCS(nb_threads),
                                                mask_nb_locks(rndup(nb_threads * nb_locks_per_thread) - 1),
                                                spinlocks_unitigs(mask_nb_locks + 1) {}

        BFG_INLINE void clear() {

            spinlocks_unitigs.clear();
        }

        BFG_INLINE void lock_unitig(const size_t unitig_id){

             spinlocks_unitigs[unitig_id & mask_nb_locks].acquire();
        }

        BFG_INLINE void unlock_unitig(const size_t unitig_id){

            spinlocks_unitigs[unitig_id & mask_nb_locks].release();
        }

    private :

        const size_t nb_locks_per_thread = 1024;
        const size_t mask_nb_locks;

        std::vector<SpinLock> spinlocks_unitigs;
};*/

#endif
