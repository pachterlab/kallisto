#ifndef BIFROST_RW_SPINLOCK_HPP
#define BIFROST_RW_SPINLOCK_HPP

#include <thread>
#include <atomic>

//#if defined(__SSE2__)
//#include <emmintrin.h>
//#endif

#include "Common.hpp"

#define HAS_WRITER 0x80000000UL
#define HAS_WRITER_WAITING 0x40000000UL
#define MASK_READER 0xBFFFFFFFUL

#define RETRY_THRESHOLD 100
#define RETRY_THRESHOLD_MAX 1023

/*class SpinLock {

    public:

        SpinLock() : _bits(0), padding{0} {}

        BFG_INLINE void acquire(){

            int retry = 0;

            while (true){

                uint32_t prev_bits = _bits;

                if ((prev_bits == 0) && _bits.compare_exchange_weak(prev_bits, 0xFFFFFFFFUL)) return;

                if (++retry > RETRY_THRESHOLD){

                    retry = 0;
                    this_thread::yield();
                }
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }
        }

        BFG_INLINE void release() {

            _bits = 0;
        }

    private:

        std::atomic<uint32_t> _bits;
        const int padding[15];
};*/

class SpinLock {

    public:

        SpinLock() : padding{0} {

            lock.clear();
        }

        BFG_INLINE void acquire() {

            while (lock.test_and_set(std::memory_order_acquire));
        }

        BFG_INLINE void release() {

            lock.clear(std::memory_order_release);
        }

    private:

        std::atomic_flag lock;
        const char padding[63];
};

class SpinLockRW {

    public:

        SpinLockRW() : _bits(0), padding{0} {}

        BFG_INLINE void acquire_reader() {

            int retry = 0;

            while (true) {

                uint32_t prev_bits = _bits;

                if (prev_bits < HAS_WRITER_WAITING) {

                    uint32_t new_bits = prev_bits + 1;

                    if (_bits.compare_exchange_weak(prev_bits, new_bits)) return;
                }

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }
        }

        BFG_INLINE void release_reader(){

            _bits.fetch_sub(1);
        }

        BFG_INLINE void acquire_writer(){

            int retry = 0;

            while (true){

                uint32_t prev_bits = _bits;

                if (((prev_bits & MASK_READER) == 0) && _bits.compare_exchange_weak(prev_bits, HAS_WRITER)) return;
                if ((prev_bits & HAS_WRITER_WAITING) == 0) _bits.fetch_or(HAS_WRITER_WAITING);

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }
        }

        BFG_INLINE void release_writer() {

            _bits = 0;
        }

        BFG_INLINE void release_writer_acquire_reader() {

            _bits = 1;
        }

        BFG_INLINE void release_all() {

            _bits = 0;
        }

    private:

        std::atomic<uint32_t> _bits;
        const int padding[15];
};

class SpinLockRW_MCS {

    public:

        SpinLockRW_MCS(const size_t nb_readers) :   writer(nullptr), lock_pool(nullptr), it_lock_pool(0),
                                                    load_lock_pool(0), mask_it(rndup(nb_readers + 1) - 1),
                                                    padding1{0}, padding2{0}, padding3{0}, padding4{0} {

            if (nb_readers <= std::thread::hardware_concurrency()){

                lock_pool = new Lock[mask_it + 1];
                lock_pool[0].is_locked = false;
            }
        }

        ~SpinLockRW_MCS() {

            clear();
        }

        BFG_INLINE void clear() {

            if (lock_pool != nullptr){

                delete[] lock_pool;
                lock_pool = nullptr;
            }

            writer = nullptr;

            it_lock_pool = 0;
            load_lock_pool = 0;
        }

        BFG_INLINE void acquire_reader() {

            uint_fast32_t retry = 0;

            const size_t prev_reader_id = it_lock_pool.fetch_add(1) & mask_it;
            const size_t new_reader_id = (prev_reader_id + 1) & mask_it;

            while (lock_pool[prev_reader_id].is_locked){

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }

            ++load_lock_pool;

            lock_pool[prev_reader_id].is_locked = true;
            lock_pool[new_reader_id].is_locked = false;
        }

        BFG_INLINE void release_reader() {

            --load_lock_pool;
        }

        BFG_INLINE void acquire_writer() {

            uint_fast32_t retry = 0;

            const size_t prev_reader_id = it_lock_pool.fetch_add(1) & mask_it;
            const size_t new_reader_id = (prev_reader_id + 1) & mask_it;

            while (lock_pool[prev_reader_id].is_locked){

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }

            while (load_lock_pool){

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }

            lock_pool[prev_reader_id].is_locked = true;

            writer = lock_pool + new_reader_id;
        }

        BFG_INLINE void release_writer() {

            writer->is_locked = false;
        }

        BFG_INLINE void release_writer_acquire_reader() {

            ++load_lock_pool;

            writer->is_locked = false;
        }

    private:

        struct Lock {

            std::atomic<bool> is_locked;
            const int padding[15];

            Lock() : is_locked(true), padding{0} {}
        };

        Lock* writer;
        const int padding1[14];
        Lock* lock_pool;
        const int padding2[14];
        const size_t mask_it;
        const int padding3[14];
        std::atomic<size_t> it_lock_pool;
        const int padding4[14];
        std::atomic<size_t> load_lock_pool;
};

template<size_t N = 4, size_t M = 256>
class Hybrid_SpinLockRW_MCS {

    public:

        Hybrid_SpinLockRW_MCS(const size_t nb_threads_) :   load(0), id(0), id_w(0), nb_threads(std::min(nb_threads_, N)), flag_full(0x1ULL << (nb_threads + 1)),
                                                            mask_id(M - 1), mask_full((0x1ULL << (nb_threads + 2)) - 1), mask_bits(mask_full - 0x1ULL - flag_full),
                                                            padding1{0}, padding2{0}, padding3{0}, padding4{0}, padding5{0}, padding6{0}, padding7{0} {

            static_assert(N <= 62, "Hybrid_SpinLockRW_MCS(): Number of reader/writer per block cannot be more than 62");

            if (nb_threads <= std::thread::hardware_concurrency()) lcks[mask_id].bits = mask_full;
            else {

                cerr << "Hybrid_SpinLockRW_MCS(): Number of threads required is greater than number of threads possible on this machine (" <<
                std::thread::hardware_concurrency() << ")" << endl;
            }
        }

        BFG_INLINE void acquire_reader() {

            uint_fast32_t retry = 0;

            const size_t l_id = id.fetch_add(1);
            const size_t new_ = (l_id / nb_threads) & mask_id;
            const size_t prev_ = (new_ - 1) & mask_id;
            const size_t mod_ = l_id % nb_threads;
            const size_t pos_ = 0x2ULL << mod_;

            while (lcks[prev_].bits != mask_full){

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }

            ++load;

            // If element is full of readers only, just unlock it and reset previous element in list
            if (((lcks[new_].bits.fetch_or(pos_) | pos_) & mask_bits) == mask_bits){

                lcks[prev_].bits = 0;
                lcks[new_].bits = mask_full;
            }
            // Else if element is full but as writers waiting, set flag_full to unlock writers
            else if (mod_ == (nb_threads - 1)) lcks[new_].bits |= flag_full;
        }

        BFG_INLINE void release_reader() {

            --load;
        }

        BFG_INLINE void acquire_writer() {

            uint_fast32_t retry = 0;

            const size_t l_id = id.fetch_add(1);
            const size_t div_ = l_id / nb_threads;
            const size_t mod_ = l_id % nb_threads;
            const size_t new_ = div_ & mask_id;
            const size_t prev_ = (new_ - 1) & mask_id;
            const size_t pos_ = 0x2ULL << mod_;
            const size_t mask_ = mask_full - ((pos_ << 1) - 1);
            const size_t end_ = div_ * nb_threads + nb_threads;

            size_t start_ = (l_id + 1);

            // Wait while previous element is not ready -> LOCAL SPIN
            // Skipping this step and just spinning on *load* would make all writers access the same ressource at the same moment
            while (lcks[prev_].bits != mask_full){

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }

            // Element is full but as writers waiting so unlocks them
            if (mod_ == (nb_threads - 1)) lcks[new_].bits |= flag_full;

            // Wait that the current element is filled in with the max readers and writers possible -> LOCAL SPIN
            while (~lcks[new_].bits & mask_){

                if (++retry > RETRY_THRESHOLD){

                    this_thread::yield();

                    if ((retry & RETRY_THRESHOLD_MAX) == 0) {

                        while ((start_ < end_) && !id.compare_exchange_strong(start_, start_ + 1));

                        if (start_ < end_){

                            const size_t mod_r = start_ % nb_threads;

                            ++start_;
                            lcks[new_].bits |= (0x2ULL << mod_r);

                            if (mod_r == (nb_threads - 1)) lcks[new_].bits |= flag_full;
                        }
                    }
                }
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }

            while (load){ // Wait that all readers are done reading -> GLOBAL SPIN

                if (++retry > RETRY_THRESHOLD) this_thread::yield();
                //#if defined(__SSE2__)
                //else _mm_pause();
                //#endif
            }

            id_w = l_id; // Set the current ticket number that is processed
        }

        BFG_INLINE void release_writer() {

            const size_t curr_ = (id_w / nb_threads) & mask_id;
            const size_t pos_ = 0x2ULL << (id_w % nb_threads);
            const size_t prev_bits = lcks[curr_].bits.fetch_or(pos_);

            if ((prev_bits | pos_) == (mask_full - 0x1ULL)){

                const size_t prev_ = (curr_ - 1) & mask_id;

                lcks[prev_].bits = 0;
                lcks[curr_].bits = mask_full;
            }
        }

    private:

        struct Lock {

            std::atomic<size_t> bits;
            const int padding[14];

            Lock() : bits(0), padding{0} {}
        };

        Lock lcks[M];
        atomic<size_t> load;
        const int padding1[14];
        atomic<size_t> id;
        const int padding2[14];
        size_t id_w;
        const int padding3[14];
        const size_t nb_threads;
        const int padding4[14];
        const size_t flag_full;
        const int padding5[14];
        const size_t mask_id;
        const int padding6[14];
        const size_t mask_full;
        const int padding7[14];
        const size_t mask_bits;
};

#endif
