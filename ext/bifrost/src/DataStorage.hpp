#ifndef BIFROST_DATA_STORAGE_HPP
#define BIFROST_DATA_STORAGE_HPP

#include <atomic>

#include "ColorSet.hpp"
#include "CompactedDBG.hpp"

#define BFG_COLOREDCDBG_FORMAT_VERSION 2

template<typename Unitig_data_t> class ColoredCDBG;
template<typename Unitig_data_t> class DataAccessor;
template<typename Unitig_data_t> class DataStorage;

template<typename U> using UnitigColorMap = UnitigMap<DataAccessor<U>, DataStorage<U>>;
template<typename U> using const_UnitigColorMap = const_UnitigMap<DataAccessor<U>, DataStorage<U>>;

namespace std
{
    template<>
    struct hash<pair<Kmer, size_t>> {

        size_t operator()(pair<Kmer, size_t> const& p) const {

            return (2 * p.second + 1) * p.first.hash();
        }
    };
}

template<typename Unitig_data_t = void>
class DataStorage {

    template<typename U> friend class ColoredCDBG;
    template<typename U> friend class DataAccessor;

    public:

        typedef Unitig_data_t U;

        DataStorage();
        DataStorage(const size_t nb_seeds_, const size_t sz_cs_, const vector<string>& color_names_);
        DataStorage(const DataStorage& o);
        DataStorage(DataStorage&& o);

        ~DataStorage();

        DataStorage<U>& operator=(const DataStorage& o);
        DataStorage<U>& operator=(DataStorage&& o);

        void clear();

        const U* getData(const const_UnitigColorMap<U>& um) const;
        U* getData(const UnitigColorMap<U>& um);

        const UnitigColors* getUnitigColors(const const_UnitigColorMap<U>& um) const;
        UnitigColors* getUnitigColors(const UnitigColorMap<U>& um);

        UnitigColors getSubUnitigColors(const const_UnitigColorMap<U>& um) const;

        vector<string> getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const;

        bool write(const string& prefix_output_filename, const bool verbose = false) const;
        bool read(const string& filename_colors, const size_t nb_threads = 1, const bool verbose = false);

        bool addUnitigColors(const UnitigColorMap<U>& um_dest, const const_UnitigColorMap<U>& um_src);
        UnitigColors joinUnitigColors(const const_UnitigColorMap<U>& um_dest, const const_UnitigColorMap<U>& um_src) const;

        inline size_t getNbColors() const { return color_names.size(); }

        size_t getUnitigColorsSize(const size_t nb_threads = 1) const;

        uint64_t getHash(const UnitigColorMap<U>& um) const;

        pair<DataAccessor<U>, pair<UnitigColors*, U*>> insert(const UnitigColorMap<U>& um, const bool force_overflow = false);
        pair<DataAccessor<U>, pair<UnitigColors*, U*>> insert(const Kmer head_unitig, const size_t unitig_sz, const bool force_overflow = false);

        void remove(const UnitigColorMap<U>& um_dest);

        void resize(const double growth = 0.1);

    private:

        void releaseMemory();

        pair<DataAccessor<U>, UnitigColors*> insert_(const Kmer head_unitig, const size_t unitig_sz, const bool force_overflow = false);

        size_t nb_seeds;

        size_t nb_cs;
        size_t sz_cs;
        size_t pos_empty_cs;

        size_t sz_shared_cs;

        uint64_t seeds[256];

        UnitigColors* color_sets;
        UnitigColors::SharedUnitigColors* shared_color_sets;

        atomic<uint64_t>* unitig_cs_link;

        U* data;

        unordered_map<pair<Kmer, size_t>, size_t> overflow;

        mutable mutex mutex_overflow;

        vector<string> color_names;
};

#endif
