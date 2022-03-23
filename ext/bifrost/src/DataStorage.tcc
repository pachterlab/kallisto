#ifndef BIFROST_DATA_STORAGE_TCC
#define BIFROST_DATA_STORAGE_TCC

template<typename U>
DataStorage<U>::DataStorage() : color_sets(nullptr), shared_color_sets(nullptr), unitig_cs_link(nullptr), data(nullptr),
                                nb_seeds(0), nb_cs(0), sz_cs(0), sz_shared_cs(0), pos_empty_cs(0) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    for (size_t i = 0; i != 256; ++i) seeds[i] = distribution(generator); //Initialize the hash function seeds
}

template<typename U>
DataStorage<U>::DataStorage(const size_t nb_seeds_, const size_t sz_cs_, const vector<string>& color_names_) :
                            color_sets(nullptr), shared_color_sets(nullptr), unitig_cs_link(nullptr), data(nullptr),
                            nb_seeds(nb_seeds_), nb_cs(sz_cs_), sz_cs(sz_cs_), sz_shared_cs(0), pos_empty_cs(0),
                            color_names(color_names_) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    const size_t sz_unitig_cs_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    for (size_t i = 0; i != 256; ++i) seeds[i] = distribution(generator); //Initialize the hash function seeds

    color_sets = new UnitigColors[sz_cs];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];
    data = new U[sz_cs];

    for (size_t i = 0; i != sz_unitig_cs_link; ++i) unitig_cs_link[i] = 0;
}

template<>
inline DataStorage<void>::DataStorage(const size_t nb_seeds_, const size_t sz_cs_, const vector<string>& color_names_) :
                                color_sets(nullptr), shared_color_sets(nullptr), unitig_cs_link(nullptr), data(nullptr),
                                nb_seeds(nb_seeds_), nb_cs(sz_cs_), sz_cs(sz_cs_), sz_shared_cs(0), pos_empty_cs(0),
                                color_names(color_names_) {

    std::random_device rd; //Seed
    std::default_random_engine generator(rd()); //Random number generator
    std::uniform_int_distribution<long long unsigned> distribution(0, 0xFFFFFFFFFFFFFFFF); //Distribution on which to apply the generator

    const size_t sz_unitig_cs_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    for (int i = 0; i < 256; ++i) seeds[i] = distribution(generator); //Initialize the hash function seeds

    color_sets = new UnitigColors[sz_cs];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];

    for (size_t i = 0; i != sz_unitig_cs_link; ++i) unitig_cs_link[i] = 0;
}

template<typename U>
DataStorage<U>::DataStorage(const DataStorage& o) : color_sets(nullptr), shared_color_sets(nullptr), unitig_cs_link(nullptr),
                                                    data(nullptr), overflow(o.overflow), nb_seeds(o.nb_seeds), nb_cs(o.nb_cs),
                                                    sz_cs(o.sz_cs), sz_shared_cs(o.sz_shared_cs), pos_empty_cs(o.pos_empty_cs),
                                                    color_names(o.color_names) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.shared_color_sets != nullptr) && (o.sz_shared_cs != 0)){

        shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];

        copy(o.shared_color_sets, o.shared_color_sets + sz_shared_cs, shared_color_sets);
    }

    if ((o.color_sets != nullptr) && (o.sz_cs != 0)){

        color_sets = new UnitigColors[sz_cs];

        for (size_t i = 0; i != sz_cs; ++i) new (&color_sets[i]) UnitigColors(o.color_sets[i], o.shared_color_sets, shared_color_sets);

        const size_t sz_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.sz_link[i].load();
    }

    if ((o.data != nullptr) && (o.sz_cs != 0)){

        data = new U[sz_cs];

        copy(o.data, o.data + sz_cs, data);
    }
}

template<>
inline DataStorage<void>::DataStorage(const DataStorage& o) :   color_sets(nullptr), shared_color_sets(nullptr), unitig_cs_link(nullptr),
                                                                data(nullptr), overflow(o.overflow), nb_seeds(o.nb_seeds), nb_cs(o.nb_cs),
                                                                sz_cs(o.sz_cs), sz_shared_cs(o.sz_shared_cs), pos_empty_cs(o.pos_empty_cs),
                                                                color_names(o.color_names) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.shared_color_sets != nullptr) && (o.sz_shared_cs != 0)){

        shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];

        copy(o.shared_color_sets, o.shared_color_sets + sz_shared_cs, shared_color_sets);
    }

    if ((o.color_sets != nullptr) && (o.sz_cs != 0)){

        color_sets = new UnitigColors[sz_cs];

        for (size_t i = 0; i != sz_cs; ++i) new (&color_sets[i]) UnitigColors(o.color_sets[i], o.shared_color_sets, shared_color_sets);

        const size_t sz_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.unitig_cs_link[i].load();
    }
}

template<typename U>
DataStorage<U>::DataStorage(DataStorage&& o) :  color_sets(o.color_sets), shared_color_sets(o.shared_color_sets), data(o.data),
                                                unitig_cs_link(o.unitig_cs_link), overflow(move(o.overflow)), nb_seeds(o.nb_seeds),
                                                nb_cs(o.nb_cs), sz_cs(o.sz_cs), sz_shared_cs(o.sz_shared_cs), pos_empty_cs(o.pos_empty_cs),
                                                color_names(move(o.color_names)) {

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    o.color_sets = nullptr;
    o.shared_color_sets = nullptr;
    o.unitig_cs_link = nullptr;
    o.data = nullptr;

    o.clear();
}

template<typename U>
void DataStorage<U>::clear() {

    nb_seeds = 0;
    nb_cs = 0;
    sz_cs = 0;
    sz_shared_cs = 0;
    pos_empty_cs = 0;

    releaseMemory();
}

template<typename U>
void DataStorage<U>::releaseMemory() {

    if (color_sets != nullptr){

        delete[] color_sets;
        color_sets = nullptr;
    }

    if (shared_color_sets != nullptr){

        delete[] shared_color_sets;
        shared_color_sets = nullptr;
    }

    if (unitig_cs_link != nullptr){

        delete[] unitig_cs_link;
        unitig_cs_link = nullptr;
    }

    if (data != nullptr){

        delete[] data;
        data = nullptr;
    }

    color_names.clear();
    overflow.clear();
}

template<>
inline void DataStorage<void>::releaseMemory() {

    if (color_sets != nullptr){

        delete[] color_sets;
        color_sets = nullptr;
    }

    if (shared_color_sets != nullptr){

        delete[] shared_color_sets;
        shared_color_sets = nullptr;
    }

    if (unitig_cs_link != nullptr){

        delete[] unitig_cs_link;
        unitig_cs_link = nullptr;
    }

    data = nullptr;

    color_names.clear();
    overflow.clear();
}

template<typename U>
DataStorage<U>::~DataStorage() {

    releaseMemory();
}

template<typename U>
DataStorage<U>& DataStorage<U>::operator=(const DataStorage& o) {

    releaseMemory();

    nb_seeds = o.nb_seeds;
    nb_cs = o.nb_cs;
    sz_cs = o.sz_cs;
    sz_shared_cs = o.sz_shared_cs;
    pos_empty_cs = o.pos_empty_cs;

    color_names = o.color_names;
    overflow = o.overflow;

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.shared_color_sets != nullptr) && (o.sz_shared_cs != 0)){

        shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];

        copy(o.shared_color_sets, o.shared_color_sets + sz_shared_cs, shared_color_sets);
    }

    if ((o.color_sets != nullptr) && (o.sz_cs != 0)){

        color_sets = new UnitigColors[sz_cs];

        for (size_t i = 0; i != sz_cs; ++i) new (&color_sets[i]) UnitigColors(o.color_sets[i], o.shared_color_sets, shared_color_sets);

        const size_t sz_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];
        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.unitig_cs_link[i].load();
    }

    if ((o.data != nullptr) && (o.sz_cs != 0)){

        data = new U[sz_cs];

        copy(o.data, o.data + sz_cs, data);
    }

    return *this;
}

template<>
inline DataStorage<void>& DataStorage<void>::operator=(const DataStorage& o) {

    releaseMemory();

    nb_seeds = o.nb_seeds;
    nb_cs = o.nb_cs;
    sz_cs = o.sz_cs;
    sz_shared_cs = o.sz_shared_cs;
    pos_empty_cs = o.pos_empty_cs;

    color_names = o.color_names;
    overflow = o.overflow;

    memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

    if ((o.shared_color_sets != nullptr) && (o.sz_shared_cs != 0)){

        shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];

        copy(o.shared_color_sets, o.shared_color_sets + sz_shared_cs, shared_color_sets);
    }

    if ((o.color_sets != nullptr) && (o.sz_cs != 0)){

        color_sets = new UnitigColors[sz_cs];

        for (size_t i = 0; i != sz_cs; ++i) new (&color_sets[i]) UnitigColors(o.color_sets[i], o.shared_color_sets, shared_color_sets);

        const size_t sz_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

        unitig_cs_link = new atomic<uint64_t>[sz_link];

        for (size_t i = 0; i != sz_link; ++i) unitig_cs_link[i] = o.unitig_cs_link[i].load();
    }

    return *this;
}

template<typename U>
DataStorage<U>& DataStorage<U>::operator=(DataStorage&& o) {

    if (this != &o) {

        releaseMemory();

        nb_seeds = o.nb_seeds;
        nb_cs = o.nb_cs;
        sz_cs = o.sz_cs;
        sz_shared_cs = o.sz_shared_cs;
        pos_empty_cs = o.pos_empty_cs;

        color_names = move(o.color_names);

        overflow = move(o.overflow);

        color_sets = o.color_sets;
        shared_color_sets = o.shared_color_sets;
        unitig_cs_link = o.unitig_cs_link;
        data = o.data;

        memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

        o.color_sets = nullptr;
        o.shared_color_sets = nullptr;
        o.unitig_cs_link = nullptr;
        o.data = nullptr;

        o.clear();
    }

    return *this;
}

template<>
inline DataStorage<void>& DataStorage<void>::operator=(DataStorage&& o) {

    if (this != &o) {

        releaseMemory();

        nb_seeds = o.nb_seeds;
        nb_cs = o.nb_cs;
        sz_cs = o.sz_cs;
        sz_shared_cs = o.sz_shared_cs;
        pos_empty_cs = o.pos_empty_cs;

        color_names = move(o.color_names);

        overflow = move(o.overflow);

        color_sets = o.color_sets;
        shared_color_sets = o.shared_color_sets;
        unitig_cs_link = o.unitig_cs_link;

        memcpy(seeds, o.seeds, 256 * sizeof(uint64_t));

        o.color_sets = nullptr;
        o.shared_color_sets = nullptr;
        o.unitig_cs_link = nullptr;

        o.clear();
    }

    return *this;
}

template<typename U>
const UnitigColors* DataStorage<U>::getUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            unique_lock<mutex> lock(mutex_overflow);

            unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it = overflow.find({head, um.size});
            if (it != overflow.end()) return &color_sets[it->second];
        }
        else return &(color_sets[head.hash(seeds[da_id - 1]) % nb_cs]);
    }

    return nullptr;
}

template<typename U>
UnitigColors* DataStorage<U>::getUnitigColors(const UnitigColorMap<U>& um) {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            unique_lock<mutex> lock(mutex_overflow);

            unordered_map<pair<Kmer, size_t>, size_t>::iterator it = overflow.find({head, um.size});
            if (it != overflow.end()) return &color_sets[it->second];
        }
        else return &(color_sets[head.hash(seeds[da_id - 1]) % nb_cs]);
    }

    return nullptr;
}

template<typename U>
const U* DataStorage<U>::getData(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (data != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            unique_lock<mutex> lock(mutex_overflow);

            unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it = overflow.find({head, um.size});
            if (it != overflow.end()) return &data[it->second];
        }
        else return &(data[head.hash(seeds[da_id - 1]) % nb_cs]);
    }

    return nullptr;
}

template<> inline const void* DataStorage<void>::getData(const const_UnitigColorMap<void>& um) const {

    return nullptr;
}

template<typename U>
U* DataStorage<U>::getData(const UnitigColorMap<U>& um) {

    if (!um.isEmpty && (data != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            unique_lock<mutex> lock(mutex_overflow);

            unordered_map<pair<Kmer, size_t>, size_t>::iterator it = overflow.find({head, um.size});
            if (it != overflow.end()) return &data[it->second];
        }
        else return &(data[head.hash(seeds[da_id - 1]) % nb_cs]);
    }

    return nullptr;
}

template<> inline void* DataStorage<void>::getData(const UnitigColorMap<void>& um) {

    return nullptr;
}

template<typename U>
UnitigColors DataStorage<U>::getSubUnitigColors(const const_UnitigColorMap<U>& um) const {

    UnitigColors new_cs;

    if (!um.isEmpty && (color_sets != nullptr)){

        const UnitigColors* cs = getUnitigColors(um);

        if (cs != nullptr){

            const size_t nb_colors = um.getGraph()->getData()->color_names.size();

            UnitigColorMap<U> um_tmp(0, 1, um.len + um.getGraph()->getK() - 1, um.strand);

            if ((um.len * nb_colors * 16) < cs->size(um)){

                const size_t end = um.dist + um.len;
                const size_t um_km_sz = um.size - um.getGraph()->getK() + 1;

                for (size_t colorID = 0; colorID < nb_colors; ++colorID){

                    for (size_t km_dist = um.dist; km_dist < end; ++km_dist){

                        if (cs->contains(colorID * um_km_sz + km_dist)){

                            um_tmp.dist = um.strand ? (km_dist - um.dist) : (um.dist + um.len - km_dist - 1);
                            new_cs.add(um_tmp, colorID);
                        }
                    }
                }
            }
            else {

                UnitigColors::const_iterator it(cs->begin(um));
                const UnitigColors::const_iterator it_end(cs->end());

                while (it != it_end){

                    const size_t km_dist = it.getKmerPosition();

                    um_tmp.dist = um.strand ? (km_dist - um.dist) : (um.dist + um.len - km_dist - 1);
                    new_cs.add(um_tmp, it.getColorID());

                    ++it;
                }
            }
        }
    }

    return new_cs;
}

template<typename U>
vector<string> DataStorage<U>::getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const {

    vector<string> v_out;

    if (!um.isEmpty && (color_sets != nullptr)){

        const UnitigColors* cs = getUnitigColors(um);

        if (cs != nullptr){

            UnitigColors::const_iterator it(cs->begin(um));
            const UnitigColors::const_iterator it_end(cs->end());

            for (; it != it_end; it.nextColor()) v_out.push_back(color_names[it.getColorID()]);
        }
    }

    return v_out;
}

template<typename U>
bool DataStorage<U>::write(const string& prefix_output_filename, const bool verbose) const {

    if (verbose) cout << endl << "DataStorage::write(): Writing colors to disk" << endl;

    const string out = prefix_output_filename + ".bfg_colors";

    FILE* fp = fopen(out.c_str(), "wb");

    if (fp == NULL) {

        cerr << "DataStorage::write(): Could not open file " << out << " for writing color sets" << endl;
        return false;
    }
    else {

        fclose(fp);

        if (std::remove(out.c_str()) != 0) cerr << "DataStorage::write(): Could not remove temporary file " << out << endl;
    }

    ofstream colorsfile_out;
    ostream colors_out(nullptr);

    colorsfile_out.open(out.c_str(), ios_base::out | ios_base::binary);
    colors_out.rdbuf(colorsfile_out.rdbuf());
    //colors_out.sync_with_stdio(false);

    const size_t format_version = BFG_COLOREDCDBG_FORMAT_VERSION;
    const size_t overflow_sz = overflow.size();
    const size_t nb_colors = color_names.size();

    const size_t block_sz = 1024;

    const char nl = '\n';

    streampos pos_f_cs;

    vector<streampos> v_pos_f_cs;

    //Write the file format version number
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&format_version), sizeof(size_t));
    //Write number of different seeds for hash function
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_seeds), sizeof(size_t));
   //Write number of colors in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_colors), sizeof(size_t));
    //Write number of color sets in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_cs), sizeof(size_t));
    //Write number of elements allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_shared_cs), sizeof(size_t));
    //Write number of (kmer, color set) overflowing
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&overflow_sz), sizeof(size_t));
    //Write the hash function seeds of the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(seeds), nb_seeds * sizeof(uint64_t));
    // Write the block size of the color sets
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&block_sz), sizeof(size_t));

    pos_f_cs = colors_out.tellp();

    const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
    const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);

    for (size_t i = 0; (i < nb_pos_shared_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_pos_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of non-shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_colors) && colors_out.good(); ++i){
        //Write the color names of the graph
        colors_out.write(color_names[i].c_str(), color_names[i].size() * sizeof(char));
        colors_out.write(&nl, sizeof(char));
    }

    for (uint64_t i = 0, j = ((sz_cs >> 6) + ((sz_cs & 0x3F) != 0)), e; (i != j) && colors_out.good(); ++i){

        e = unitig_cs_link[i].load();
        colors_out.write(reinterpret_cast<const char*>(&e), sizeof(uint64_t));
    }

    for (size_t i = 0; (i < sz_shared_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        if (shared_color_sets[i].first.write(colors_out)){

            colors_out.write(reinterpret_cast<const char*>(&(shared_color_sets[i].second)), sizeof(size_t));
        }
    }

    for (size_t i = 0; (i < sz_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        color_sets[i].write(colors_out, false); //Write the color sets
    }

    unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it(overflow.begin());
    const unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it_end(overflow.end());

    for (; (it != it_end) && colors_out.good(); ++it){

        it->first.first.write(colors_out); // Write the k-mer

        colors_out.write(reinterpret_cast<const char*>(&(it->first.second)), sizeof(size_t)); // Write the unitig length
        colors_out.write(reinterpret_cast<const char*>(&(it->second)), sizeof(size_t)); // Write the position
    }

    if (colors_out.good()){

        colors_out.seekp(pos_f_cs); // Re-position cursor to array of position of shared color sets at the beginning

        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&v_pos_f_cs[0]), v_pos_f_cs.size() * sizeof(streampos));
    }

    const bool ret = colors_out.good();

    colorsfile_out.close();

    return ret;
}

template<>
inline bool DataStorage<void>::write(const string& prefix_output_filename, const bool verbose) const {

    if (verbose) cout << endl << "DataStorage::write(): Writing colors to disk" << endl;

    const string out = prefix_output_filename + ".bfg_colors";

    FILE* fp = fopen(out.c_str(), "wb");

    if (fp == NULL) {

        cerr << "DataStorage::write(): Could not open file " << out << " for writing color sets" << endl;
        return false;
    }
    else {

        fclose(fp);

        if (std::remove(out.c_str()) != 0) cerr << "DataStorage::write(): Could not remove temporary file " << out << endl;
    }

    ofstream colorsfile_out;
    ostream colors_out(nullptr);

    colorsfile_out.open(out.c_str(), ios_base::out | ios_base::binary);
    colors_out.rdbuf(colorsfile_out.rdbuf());
    //colors_out.sync_with_stdio(false);

    const size_t format_version = BFG_COLOREDCDBG_FORMAT_VERSION;
    const size_t overflow_sz = overflow.size();
    const size_t nb_colors = color_names.size();

    const size_t block_sz = 1024;

    const char nl = '\n';

    streampos pos_f_cs;

    vector<streampos> v_pos_f_cs;

    //Write the file format version number
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&format_version), sizeof(size_t));
    //Write number of different seeds for hash function
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_seeds), sizeof(size_t));
   //Write number of colors in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_colors), sizeof(size_t));
    //Write number of color sets in the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&nb_cs), sizeof(size_t));
    //Write number of elements allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&sz_shared_cs), sizeof(size_t));
    //Write number of (kmer, color set) overflowing
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&overflow_sz), sizeof(size_t));
    //Write the hash function seeds of the graph
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(seeds), nb_seeds * sizeof(uint64_t));
    // Write the block size of the color sets
    if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&block_sz), sizeof(size_t));

    pos_f_cs = colors_out.tellp();

    const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
    const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);

    for (size_t i = 0; (i < nb_pos_shared_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_pos_cs) && colors_out.good(); ++i){
        // Reserve space to write positions in file of non-shared colorsets blocks
        colors_out.write(reinterpret_cast<const char*>(&pos_f_cs), sizeof(streampos));
    }

    for (size_t i = 0; (i < nb_colors) && colors_out.good(); ++i){
        //Write the color names of the graph
        colors_out.write(color_names[i].c_str(), color_names[i].size() * sizeof(char));
        colors_out.write(&nl, sizeof(char));
    }

    for (uint64_t i = 0, j = ((sz_cs >> 6) + ((sz_cs & 0x3F) != 0)), e; (i != j) && colors_out.good(); ++i){

        e = unitig_cs_link[i].load();
        colors_out.write(reinterpret_cast<const char*>(&e), sizeof(uint64_t));
    }

    for (size_t i = 0; (i < sz_shared_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        if (shared_color_sets[i].first.write(colors_out)){

            colors_out.write(reinterpret_cast<const char*>(&(shared_color_sets[i].second)), sizeof(size_t));
        }
    }

    for (size_t i = 0; (i < sz_cs) && colors_out.good(); ++i){

        if (i % block_sz == 0) v_pos_f_cs.push_back(colors_out.tellp());

        color_sets[i].write(colors_out, false); //Write the color sets
    }

    unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it(overflow.begin());
    const unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it_end(overflow.end());

    for (; (it != it_end) && colors_out.good(); ++it){

        it->first.first.write(colors_out); // Write the k-mer

        colors_out.write(reinterpret_cast<const char*>(&(it->first.second)), sizeof(size_t)); // Write the unitig length
        colors_out.write(reinterpret_cast<const char*>(&(it->second)), sizeof(size_t)); // Write the position
    }

    if (colors_out.good()){

        colors_out.seekp(pos_f_cs); // Re-position cursor to array of position of shared color sets at the beginning

        if (colors_out.good()) colors_out.write(reinterpret_cast<const char*>(&v_pos_f_cs[0]), v_pos_f_cs.size() * sizeof(streampos));
    }

    const bool ret = colors_out.good();

    colorsfile_out.close();

    return ret;
}

template<typename U>
bool DataStorage<U>::read(const string& filename_colors, const size_t nb_threads, const bool verbose) {

    if (verbose) cout << endl << "DataStorage::read(): Reading color sets from disk" << endl;

    FILE* fp = fopen(filename_colors.c_str(), "rb");

    if (fp == NULL) {

        cerr << "DataStorage::read(): Could not open file " << filename_colors << " for reading color sets" << endl;
        return false;
    }
    else fclose(fp);

    auto readSharedColorSets = [](UnitigColors::SharedUnitigColors* shared_color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i){

            if (shared_color_sets[i].first.read(colors_in)){

                colors_in.read(reinterpret_cast<char*>(&(shared_color_sets[i].second)), sizeof(size_t));
            }
        }
    };

    auto readColorSets = [](UnitigColors* color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i) color_sets[i].read(colors_in);
    };

    Kmer km;

    size_t format_version, overflow_sz, nb_colors;

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    clear();

    //Read the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Read number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Read number of colors in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_colors), sizeof(size_t));
    //Read number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_cs), sizeof(size_t));
    //Read number of elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_shared_cs), sizeof(size_t));
    //Read number of (kmer, color set) overflowing
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&overflow_sz), sizeof(size_t));

    if (nb_seeds >= 256){

        cerr << "DataStorage::read(): Does not support more than 255 hash seeds" << endl;
        return false;
    }

    const size_t sz_unitig_cs_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    overflow = unordered_map<pair<Kmer, size_t>, size_t>(overflow_sz);

    color_sets = new UnitigColors[sz_cs];
    shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];
    data = new U[sz_cs];

    //Read the hash function seeds of the graph
    colors_in.read(reinterpret_cast<char*>(seeds), nb_seeds * sizeof(uint64_t));

    if (format_version == 1){

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
        readColorSets(color_sets, colors_in, sz_cs);
    }
    else {

        size_t block_sz = 0;

        streampos* pos_f_cs = nullptr;

        //Read the hash function seeds of the graph
        if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&block_sz), sizeof(size_t));

        const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
        const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);
        const size_t pos_f_cs_sz = nb_pos_shared_cs + nb_pos_cs;

        if (pos_f_cs_sz != 0){

            pos_f_cs = new streampos[pos_f_cs_sz];

            if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(pos_f_cs), pos_f_cs_sz * sizeof(streampos));
        }

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        if ((nb_threads == 1) || (pos_f_cs_sz == 0)) {

            readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
            readColorSets(color_sets, colors_in, sz_cs);
        }
        else {

            streampos colors_in_pos = colors_in.tellg();

            colorsfile_in.close();

            mutex m_colors_in_pos;

            vector<thread> workers; // need to keep track of threads so we can join them

            std::atomic<size_t> i;

            i = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            const size_t l_i = i++;

                            if (l_i >= nb_pos_shared_cs){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);
                            readSharedColorSets(shared_color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_shared_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            workers.clear();

            i = nb_pos_shared_cs;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            size_t l_i = i++;

                            if (l_i >= pos_f_cs_sz){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);

                            l_i -= nb_pos_shared_cs;

                            readColorSets(color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
            colors_in.rdbuf(colorsfile_in.rdbuf());
            colors_in.seekg(colors_in_pos);
        }

        if (pos_f_cs != nullptr) delete[] pos_f_cs;
    }

    for (size_t i = 0, sz, pos; (i < overflow_sz) && colors_in.good(); ++i){

        km.read(colors_in);

        colors_in.read(reinterpret_cast<char*>(&sz), sizeof(size_t));
        colors_in.read(reinterpret_cast<char*>(&pos), sizeof(size_t));

        overflow.insert({{km, sz}, pos});
    }

    const bool ret = colors_in.good();

    colorsfile_in.close();

    return ret;
}

template<>
inline bool DataStorage<void>::read(const string& filename_colors, const size_t nb_threads, const bool verbose) {

    if (verbose) cout << endl << "DataStorage::read(): Reading color sets from disk" << endl;

    FILE* fp = fopen(filename_colors.c_str(), "rb");

    if (fp == NULL) {

        cerr << "DataStorage::read(): Could not open file " << filename_colors << " for reading color sets" << endl;
        return false;
    }
    else fclose(fp);

    auto readSharedColorSets = [](UnitigColors::SharedUnitigColors* shared_color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i){

            if (shared_color_sets[i].first.read(colors_in)){

                colors_in.read(reinterpret_cast<char*>(&(shared_color_sets[i].second)), sizeof(size_t));
            }
        }
    };

    auto readColorSets = [](UnitigColors* color_sets, istream& colors_in, const size_t sz){

        for (size_t i = 0; (i < sz) && colors_in.good(); ++i) color_sets[i].read(colors_in);
    };

    Kmer km;

    size_t format_version, overflow_sz, nb_colors;

    ifstream colorsfile_in;
    istream colors_in(nullptr);

    colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
    colors_in.rdbuf(colorsfile_in.rdbuf());

    clear();

    //Read the file format version number
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&format_version), sizeof(size_t));
    //Read number of different seeds for hash function
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_seeds), sizeof(size_t));
    //Read number of colors in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_colors), sizeof(size_t));
    //Read number of color sets in the graph
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&nb_cs), sizeof(size_t));
    //Read number of elements allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_cs), sizeof(size_t));
    //Write number of SharedUnitigColors allocated
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&sz_shared_cs), sizeof(size_t));
    //Read number of (kmer, color set) overflowing
    if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&overflow_sz), sizeof(size_t));

    if (nb_seeds >= 256){

        cerr << "DataStorage::read(): Does not support more than 255 hash seeds" << endl;
        return false;
    }

    const size_t sz_unitig_cs_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    overflow = unordered_map<pair<Kmer, size_t>, size_t>(overflow_sz);

    color_sets = new UnitigColors[sz_cs];
    shared_color_sets = new UnitigColors::SharedUnitigColors[sz_shared_cs];
    unitig_cs_link = new atomic<uint64_t>[sz_unitig_cs_link];

    //Read the hash function seeds of the graph
    colors_in.read(reinterpret_cast<char*>(seeds), nb_seeds * sizeof(uint64_t));

    if (format_version == 1){

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
        readColorSets(color_sets, colors_in, sz_cs);
    }
    else {

        size_t block_sz = 0;

        streampos* pos_f_cs = nullptr;

        //Read the hash function seeds of the graph
        if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(&block_sz), sizeof(size_t));

        const size_t nb_pos_shared_cs = (sz_shared_cs / block_sz) + static_cast<size_t>((sz_shared_cs % block_sz) != 0);
        const size_t nb_pos_cs = (sz_cs / block_sz) + static_cast<size_t>((sz_cs % block_sz) != 0);
        const size_t pos_f_cs_sz = nb_pos_shared_cs + nb_pos_cs;

        if (pos_f_cs_sz != 0){

            pos_f_cs = new streampos[pos_f_cs_sz];

            if (colors_in.good()) colors_in.read(reinterpret_cast<char*>(pos_f_cs), pos_f_cs_sz * sizeof(streampos));
        }

        for (size_t i = 0; (i < nb_colors) && colors_in.good(); ++i){
            //Read the hash function seeds of the graph
            color_names.push_back(string());
            getline(colors_in, color_names[i]);
        }

        for (uint64_t i = 0, e; (i != sz_unitig_cs_link) && colors_in.good(); ++i){

            colors_in.read(reinterpret_cast<char*>(&e), sizeof(uint64_t));
            unitig_cs_link[i] = e;
        }

        if ((nb_threads == 1) || (pos_f_cs_sz == 0)) {

            readSharedColorSets(shared_color_sets, colors_in, sz_shared_cs);
            readColorSets(color_sets, colors_in, sz_cs);
        }
        else {

            streampos colors_in_pos = colors_in.tellg();

            colorsfile_in.close();

            mutex m_colors_in_pos;

            vector<thread> workers; // need to keep track of threads so we can join them

            std::atomic<size_t> i;

            i = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            const size_t l_i = i++;

                            if (l_i >= nb_pos_shared_cs){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);
                            readSharedColorSets(shared_color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_shared_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            workers.clear();

            i = nb_pos_shared_cs;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        ifstream colorsfile_in_t;
                        istream colors_in_t(nullptr);

                        colorsfile_in_t.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
                        colors_in_t.rdbuf(colorsfile_in_t.rdbuf());

                        while (true) {

                            size_t l_i = i++;

                            if (l_i >= pos_f_cs_sz){

                                const streampos colors_in_t_pos = colors_in_t.tellg();

                                {
                                    unique_lock<mutex> lock(m_colors_in_pos);

                                    colors_in_pos = max(colors_in_pos, colors_in_t_pos);
                                }

                                colorsfile_in_t.close();

                                break;
                            }

                            colors_in_t.seekg(pos_f_cs[l_i]);

                            l_i -= nb_pos_shared_cs;

                            readColorSets(color_sets + (l_i * block_sz), colors_in_t, min(block_sz, sz_cs - (l_i * block_sz)));
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();

            colorsfile_in.open(filename_colors.c_str(), ios_base::in | ios_base::binary);
            colors_in.rdbuf(colorsfile_in.rdbuf());
            colors_in.seekg(colors_in_pos);
        }

        if (pos_f_cs != nullptr) delete[] pos_f_cs;
    }

    for (size_t i = 0, sz, pos; (i < overflow_sz) && colors_in.good(); ++i){

        km.read(colors_in);

        colors_in.read(reinterpret_cast<char*>(&sz), sizeof(size_t));
        colors_in.read(reinterpret_cast<char*>(&pos), sizeof(size_t));

        overflow.insert({{km, sz}, pos});
    }

    const bool ret = colors_in.good();

    colorsfile_in.close();

    return ret;
}

template<typename U>
bool DataStorage<U>::addUnitigColors(const UnitigColorMap<U>& um_dest, const const_UnitigColorMap<U>& um_src) {

    if (!um_dest.isEmpty && !um_src.isEmpty && (um_src.len == um_dest.len) && (color_sets != nullptr)){

        UnitigColors* cs_dest = getUnitigColors(um_dest);

        const UnitigColors* cs_src = um_src.getGraph()->getData()->getUnitigColors(um_src);

        if ((cs_dest != nullptr) && (cs_src != nullptr)){

            const size_t nb_colors_src = um_src.getGraph()->getData()->color_names.size();
            const size_t nb_colors_dest = um_dest.getGraph()->getData()->color_names.size() - nb_colors_src;

            UnitigColorMap<U> um_tmp(0, 1, um_dest.size, um_dest.strand);

            if ((um_src.len * nb_colors_src * 16) < cs_src->size(um_src)){

                const size_t pos_end_src = um_src.dist + um_src.len;
                const size_t sz_km_src = um_src.size - Kmer::k + 1;

                for (size_t colorID = 0; colorID < nb_colors_src; ++colorID){

                    for (size_t km_dist = um_src.dist; km_dist < pos_end_src; ++km_dist){

                        if (cs_src->contains(colorID * sz_km_src + km_dist)){

                            if (um_dest.strand != um_src.strand) um_tmp.dist = (um_src.len - 1 - (km_dist - um_src.dist)) + um_dest.dist;
                            else um_tmp.dist = (km_dist - um_src.dist) + um_dest.dist;

                            cs_dest->add(um_tmp, colorID + nb_colors_dest);
                        }
                    }
                }
            }
            else {

                UnitigColors::const_iterator it(cs_src->begin(um_src));
                const UnitigColors::const_iterator it_end(cs_src->end());

                while (it != it_end){

                    const size_t km_dist = it.getKmerPosition();

                    if (um_dest.strand != um_src.strand) um_tmp.dist = (um_src.len - 1 - (km_dist - um_src.dist)) + um_dest.dist;
                    else um_tmp.dist = (km_dist - um_src.dist) + um_dest.dist;

                    cs_dest->add(um_tmp, it.getColorID() + nb_colors_dest);

                    ++it;
                }
            }
        }
    }

    return false;
}

template<typename U>
UnitigColors DataStorage<U>::joinUnitigColors(const const_UnitigColorMap<U>& um_dest, const const_UnitigColorMap<U>& um_src) const {

    UnitigColors new_cs;

    if (!um_dest.isEmpty && !um_src.isEmpty && (color_sets != nullptr)){

        const UnitigColors* color_set_dest = um_dest.getGraph()->getData()->getUnitigColors(um_dest);
        const UnitigColors* color_set_src = um_src.getGraph()->getData()->getUnitigColors(um_src);

        if ((color_set_dest != nullptr) || (color_set_src != nullptr)){

            size_t prev_color_id = 0xffffffffffffffff;
            size_t prev_km_dist = 0xffffffffffffffff;

            const size_t k = um_dest.getGraph()->getK();
            const size_t um_dest_km_sz = um_dest.size - k + 1;
            const size_t um_src_km_sz = um_src.size - k + 1;

            if (color_set_dest != nullptr){

                UnitigColors csd_rev;

                UnitigColorMap<U> new_um_dest(0, 0, um_dest.size + um_src.size - k + 1, um_dest.strand);

                if (!um_dest.strand){

                    csd_rev = color_set_dest->reverse(um_dest);
                    color_set_dest = &csd_rev;
                }

                UnitigColors::const_iterator it(color_set_dest->begin(0, um_dest_km_sz, um_dest_km_sz));
                const UnitigColors::const_iterator it_end(color_set_dest->end());

                if (it != it_end){

                    prev_km_dist = it.getKmerPosition();
                    prev_color_id = it.getColorID();

                    new_um_dest.dist = prev_km_dist;
                    new_um_dest.len = 1;

                    ++it;
                }

                // Insert colors layer by layer
                for (; it != it_end; ++it){

                    const size_t km_dist = it.getKmerPosition();
                    const size_t color_id = it.getColorID();

                    if ((color_id != prev_color_id) || (km_dist != prev_km_dist + 1)){

                        new_cs.add(new_um_dest, prev_color_id);

                        new_um_dest.dist = km_dist;
                        new_um_dest.len = 1;
                    }
                    else ++(new_um_dest.len);

                    prev_color_id = color_id;
                    prev_km_dist = km_dist;
                }

                if (new_um_dest.dist + new_um_dest.len != 0) new_cs.add(new_um_dest, prev_color_id);
            }

            if (color_set_src != nullptr){

                UnitigColors css_rev;

                UnitigColorMap<U> new_um_src(0, 0, um_dest.size + um_src.size - k + 1, um_src.strand);

                if (!um_src.strand){

                    css_rev = color_set_src->reverse(um_src);
                    color_set_src = &css_rev;
                }

                UnitigColors::const_iterator it(color_set_src->begin(0, um_src_km_sz, um_src_km_sz));
                const UnitigColors::const_iterator it_end(color_set_src->end());

                if (it != it_end){

                    prev_km_dist = it.getKmerPosition();
                    prev_color_id = it.getColorID();

                    new_um_src.dist = prev_km_dist + um_dest.size - k + 1;
                    new_um_src.len = 1;

                    ++it;
                }

                // Insert colors layer by layer
                for (; it != it_end; ++it){

                    const size_t km_dist = it.getKmerPosition();
                    const size_t color_id = it.getColorID();

                    if ((color_id != prev_color_id) || (km_dist != prev_km_dist + 1)){

                        new_cs.add(new_um_src, prev_color_id);

                        new_um_src.dist = km_dist + um_dest.size - k + 1;
                        new_um_src.len = 1;
                    }
                    else ++(new_um_src.len);

                    prev_color_id = color_id;
                    prev_km_dist = km_dist;
                }

                if (new_um_src.dist + new_um_src.len != 0) new_cs.add(new_um_src, prev_color_id);
            }
        }
    }

    return new_cs;
}

template<typename U>
void DataStorage<U>::remove(const UnitigColorMap<U>& um) {

    if (!um.isEmpty && (color_sets != nullptr) && (data != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            unique_lock<mutex> lock(mutex_overflow);

            unordered_map<pair<Kmer, size_t>, size_t>::iterator it = overflow.find({head, um.size});

            if (it != overflow.end()){

                const uint64_t bit = 1ULL << (it->second & 0x3F);
                const uint64_t mask = 0xFFFFFFFFFFFFFFFFULL - bit;

                if ((unitig_cs_link[it->second >> 6].fetch_and(mask) & bit) == bit){

                    color_sets[it->second].clear();
                    data[it->second].clear(um);

                    overflow.erase(it);
                }
            }
        }
        else {

            const uint64_t h_v = head.hash(seeds[da_id - 1]) % nb_cs;
            const uint64_t bit = 1ULL << (h_v & 0x3F);
            const uint64_t mask = 0xFFFFFFFFFFFFFFFFULL - bit;

            if ((unitig_cs_link[h_v >> 6].fetch_and(mask) & bit) == bit){

                color_sets[h_v].clear();
                data[h_v].clear(um);
            }
        }
    }
}

template<>
inline void DataStorage<void>::remove(const UnitigColorMap<void>& um) {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            unique_lock<mutex> lock(mutex_overflow);

            unordered_map<pair<Kmer, size_t>, size_t>::iterator it = overflow.find({head, um.size});

            if (it != overflow.end()){

                const uint64_t bit = 1ULL << (it->second & 0x3F);
                const uint64_t mask = 0xFFFFFFFFFFFFFFFFULL - bit;

                if ((unitig_cs_link[it->second >> 6].fetch_and(mask) & bit) == bit){

                    color_sets[it->second].clear();

                    overflow.erase(it);
                }
            }
        }
        else {

            const uint64_t h_v = head.hash(seeds[da_id - 1]) % nb_cs;
            const uint64_t bit = 1ULL << (h_v & 0x3F);
            const uint64_t mask = 0xFFFFFFFFFFFFFFFFULL - bit;

            if ((unitig_cs_link[h_v >> 6].fetch_and(mask) & bit) == bit) color_sets[h_v].clear();
        }
    }
}

template<typename U>
size_t DataStorage<U>::getUnitigColorsSize(const size_t nb_threads) const {

    if (color_sets != nullptr){

        atomic<size_t> sz_in_bytes(0);

        auto worker_function = [&sz_in_bytes, this](const size_t idx_start, const size_t idx_end){

            size_t cpt = 0;

            for (size_t i = idx_start; i < idx_end; ++i) cpt += color_sets[i].getSizeInBytes();

            sz_in_bytes += cpt;
        };

        vector<thread> workers;

        const size_t load_per_thread = nb_cs / nb_threads;

        size_t t = 0;

        for (; t < (nb_threads - 1) * load_per_thread; t += load_per_thread) workers.emplace_back(worker_function, t, t + load_per_thread);

        workers.emplace_back(worker_function, t, nb_cs);

        for (auto& t : workers) t.join();

        return sz_in_bytes;
    }

    return 0;
}

template<typename U>
uint64_t DataStorage<U>::getHash(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (color_sets != nullptr)){

        const Kmer head = um.getUnitigHead();
        const uint8_t da_id = um.getData()->get();

        if (da_id == 0){

            unique_lock<mutex> lock(mutex_overflow);

            unordered_map<pair<Kmer, size_t>, size_t>::const_iterator it = overflow.find({head, um.size});
            if (it != overflow.end()) return it->second;
        }
        else return head.hash(seeds[da_id - 1]) % nb_cs;
    }

    return 0;
}

template<typename U>
pair<DataAccessor<U>, UnitigColors*> DataStorage<U>::insert_(const Kmer head, const size_t unitig_sz, const bool force_overflow) {

    if (color_sets == nullptr) return {DataAccessor<U>(0), nullptr};

    uint64_t i, pos, id_link_mod;

    for (i = (force_overflow ? nb_seeds : 0); i < nb_seeds; ++i){

        pos = head.hash(seeds[i]) % nb_cs; // Hash to which we can possibly put our colorset for current kmer
        id_link_mod = 1ULL << (pos & 0x3F);

        if ((unitig_cs_link[pos >> 6].fetch_or(id_link_mod) & id_link_mod) == 0) break;
    }

    if (i == nb_seeds){ // IF we couldn't find a hash matching an unoccupied color set for current k-mer

        unique_lock<mutex> lock(mutex_overflow);

        size_t j = 0;

        id_link_mod = 1ULL << (pos_empty_cs & 0x3F);

        while ((j != sz_cs) && ((unitig_cs_link[pos_empty_cs >> 6].fetch_or(id_link_mod) & id_link_mod) != 0)){

            ++j;
            ++pos_empty_cs;

            if (pos_empty_cs == sz_cs) pos_empty_cs = 0;

            id_link_mod = 1ULL << (pos_empty_cs & 0x3F);
        }

        if (j == sz_cs){ // No empty slot could be found,

            resize();

            j = 0;
            id_link_mod = 1ULL << (pos_empty_cs & 0x3F);

            while ((j != sz_cs) && ((unitig_cs_link[pos_empty_cs >> 6].fetch_or(id_link_mod) & id_link_mod) != 0)){

                ++j;
                ++pos_empty_cs;

                if (pos_empty_cs == sz_cs) pos_empty_cs = 0;

                id_link_mod = 1ULL << (pos_empty_cs & 0x3F);
            }
        }

        pos = pos_empty_cs;

        overflow.insert({{head, unitig_sz}, pos}); // Insertion

    }

    return {DataAccessor<U>(static_cast<uint8_t>(i == nb_seeds ? 0 : i + 1)), &color_sets[pos]};
}

template<typename U>
pair<DataAccessor<U>, pair<UnitigColors*, U*>> DataStorage<U>::insert(const Kmer head, const size_t unitig_sz, const bool force_overflow) {

    if ((color_sets == nullptr) && (data == nullptr)) return {DataAccessor<U>(0), {nullptr, nullptr}};

    const pair<DataAccessor<U>, UnitigColors*> p(insert_(head, unitig_sz, force_overflow));

    return {p.first, {p.second, &data[p.second - color_sets]}};
}

template<>
inline pair<DataAccessor<void>, pair<UnitigColors*, void*>> DataStorage<void>::insert(const Kmer head, const size_t unitig_sz, const bool force_overflow) {

    if ((color_sets == nullptr) && (data == nullptr)) return {DataAccessor<void>(0), {nullptr, nullptr}};

    const pair<DataAccessor<void>, UnitigColors*> p(insert_(head, unitig_sz, force_overflow));

    return {p.first, {p.second, nullptr}};
}

template<typename U>
pair<DataAccessor<U>, pair<UnitigColors*, U*>> DataStorage<U>::insert(const UnitigColorMap<U>& um, const bool force_overflow) {

    if ((color_sets == nullptr) && (data == nullptr)) return {DataAccessor<U>(0), {nullptr, nullptr}};

    const Kmer head(um.getMappedHead());
    const pair<DataAccessor<U>, UnitigColors*> p(insert_(head, um.len + um.getGraph()->getK() - 1, force_overflow));

    return {p.first, {p.second, &data[p.second - color_sets]}};
}

template<>
inline pair<DataAccessor<void>, pair<UnitigColors*, void*>> DataStorage<void>::insert(const UnitigColorMap<void>& um, const bool force_overflow) {

    if ((color_sets == nullptr) && (data == nullptr)) return {DataAccessor<void>(0), {nullptr, nullptr}};

    const Kmer head(um.getMappedHead());
    const pair<DataAccessor<void>, UnitigColors*> p(insert_(head, um.len + um.getGraph()->getK() - 1, force_overflow));

    return {p.first, {p.second, nullptr}};
}

template<typename U>
void DataStorage<U>::resize(const double growth) {

    UnitigColors* old_color_sets = color_sets;
    atomic<uint64_t>* old_unitig_cs_link = unitig_cs_link;
    U* old_data = data;

    const size_t old_sz_cs = sz_cs;
    const size_t old_sz_link = (old_sz_cs >> 6) + ((old_sz_cs & 0x3F) != 0);

    sz_cs += sz_cs * growth;

    const size_t sz_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    // Reallocate UnitigColors
    color_sets = new UnitigColors[sz_cs];

    move(old_color_sets, old_color_sets + old_sz_cs, color_sets);
    delete[] old_color_sets;

    // Reallocate UnitigColors occupancy
    unitig_cs_link = new atomic<uint64_t>[sz_link];

    for (size_t i = 0; i != old_sz_link; ++i) unitig_cs_link[i] = old_unitig_cs_link[i].load();
    for (size_t i = old_sz_link; i != sz_link; ++i) unitig_cs_link[i] = 0;

    delete[] old_unitig_cs_link;

    // Reallocate data
    data = new U[sz_cs];

    move(old_data, old_data + old_sz_cs, data);
    delete[] old_data;
}

template<>
inline void DataStorage<void>::resize(const double growth) {

    UnitigColors* old_color_sets = color_sets;
    atomic<uint64_t>* old_unitig_cs_link = unitig_cs_link;

    const size_t old_sz_cs = sz_cs;
    const size_t old_sz_link = (old_sz_cs >> 6) + ((old_sz_cs & 0x3F) != 0);

    sz_cs += sz_cs * growth;

    const size_t sz_link = (sz_cs >> 6) + ((sz_cs & 0x3F) != 0);

    color_sets = new UnitigColors[sz_cs];

    move(old_color_sets, old_color_sets + old_sz_cs, color_sets);
    delete[] old_color_sets;

    unitig_cs_link = new atomic<uint64_t>[sz_link];

    for (size_t i = 0; i != old_sz_link; ++i) unitig_cs_link[i] = old_unitig_cs_link[i].load();
    for (size_t i = old_sz_link; i != sz_link; ++i) unitig_cs_link[i] = 0;

    delete[] old_unitig_cs_link;
}

#endif
