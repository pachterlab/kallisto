#ifndef BIFROST_HASHID_TCC
#define BIFROST_HASHID_TCC

template<typename U>
DataAccessor<U>::DataAccessor(const uint8_t id) : da_id(id) {}

template<typename U>
void DataAccessor<U>::clear(const UnitigColorMap<U>& um) {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr){

            ds->remove(um);

            da_id = 0;
        }
    }
}

template<typename U>
const U* DataAccessor<U>::getData(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        const DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr) return ds->getData(um);
    }

    return nullptr;
}

template<>
inline const void* DataAccessor<void>::getData(const const_UnitigColorMap<void>& um) const {

    return nullptr;
}

template<typename U>
U* DataAccessor<U>::getData(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr) return ds->getData(um);
    }

    return nullptr;
}

template<>
inline void* DataAccessor<void>::getData(const UnitigColorMap<void>& um) const {

    return nullptr;
}

template<typename U>
const UnitigColors* DataAccessor<U>::getUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        const DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr) return ds->getUnitigColors(um);
    }

    return nullptr;
}

template<typename U>
UnitigColors* DataAccessor<U>::getUnitigColors(const UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr) return ds->getUnitigColors(um);
    }

    return nullptr;
}

template<typename U>
UnitigColors DataAccessor<U>::getSubUnitigColors(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        const DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr) return ds->getSubUnitigColors(um);
    }

    return UnitigColors();
}

template<typename U>
vector<string> DataAccessor<U>::getSubUnitigColorNames(const const_UnitigColorMap<U>& um) const {

    if (!um.isEmpty && (um.getGraph() != nullptr)){

        const DataStorage<U>* ds = um.getGraph()->getData();

        if (ds != nullptr) return ds->getSubUnitigColorNames(um);
    }

    return vector<string>();
}

template<typename U>
void DataAccessor<U>::concat(const UnitigColorMap<U>& um_dest, const UnitigColorMap<U>& um_src){

    DataStorage<U>* ds = um_dest.getGraph()->getData();

    DataAccessor<U>* da_dest = um_dest.getData();
    UnitigColors* uc_dest = da_dest->getUnitigColors(um_dest);
    U* data_dest = da_dest->getData(um_dest);

    DataAccessor<U>* da_src = um_src.getData();
    UnitigColors* uc_src = da_src->getUnitigColors(um_src);
    U* data_src = da_src->getData(um_src);

    if ((uc_dest != nullptr) || (uc_src != nullptr)){ // If a colorset exists for um_dest

        const Kmer head = um_dest.getUnitigHead();
        const Kmer new_head = um_dest.strand ? head : um_dest.getUnitigTail().twin();

        UnitigColors* uc_ptr = uc_dest;
        U* data_ptr = data_dest;

        UnitigColors uc_ = ds->joinUnitigColors(um_dest, um_src); // Join the color sets
        U data_;

        if ((data_dest != nullptr) || (data_src != nullptr)) data_.concat(um_dest, um_src);

        if ((uc_ptr == nullptr) || (head != new_head) || (da_dest->get() == 0)){
             // If (da_dest->get() == 0, UnitigColors of um_dest cannot be recycled 'cause new unitig length will be different from the one of um_dest
            const pair<DataAccessor<U>, pair<UnitigColors*, U*>> p  = ds->insert(new_head, um_dest.size + um_src.size - um_dest.getGraph()->getK() + 1);

            uc_ptr = p.second.first;
            data_ptr = p.second.second;
            *this = p.first;
        }
        else {

            *this = *da_dest; // UnitigColors of um_dest is recycled
            *da_dest = DataAccessor<U>(0);
        }

        *uc_ptr = move(uc_);
        *data_ptr = move(data_);
    }
}

template<>
inline void DataAccessor<void>::concat(const UnitigColorMap<void>& um_dest, const UnitigColorMap<void>& um_src){

    DataStorage<void>* ds = um_dest.getGraph()->getData();

    DataAccessor<void>* da_dest = um_dest.getData();
    UnitigColors* uc_dest = da_dest->getUnitigColors(um_dest);

    DataAccessor<void>* da_src = um_src.getData();
    UnitigColors* uc_src = da_src->getUnitigColors(um_src);

    if ((uc_dest != nullptr) || (uc_src != nullptr)){ // If a UnitigColors exists for at least one of the two unitigs

        const Kmer head = um_dest.getUnitigHead(); // Get head k-mer of reference unitig of um_dest
        const Kmer new_head = um_dest.strand ? head : um_dest.getUnitigTail().twin(); // Compute new head k-mer after concatenation

        UnitigColors uc_ = ds->joinUnitigColors(um_dest, um_src); // Join the UnitigColors of the two reference unitigs

        UnitigColors* uc_ptr = uc_dest;

        // If reference unitig of um_dest has no UnitigColors associated
        // OR if the new head k-mer (after concat.) is not the same has the current head k-mer (before concat.)
        // OR if UnitigColors associated with reference unitig of um_dest is in the overflow system of DataStorage (da_dest->get() == 0)
        if ((uc_ptr == nullptr) || (head != new_head) || (da_dest->get() == 0)){
            //
            // If (da_dest->get() == 0, UnitigColors of um_dest cannot be recycled 'cause new unitig length will be different from the one of um_dest
            const pair<DataAccessor<void>, pair<UnitigColors*, void*>> p  = ds->insert(new_head, um_dest.size + um_src.size - um_dest.getGraph()->getK() + 1);

            uc_ptr = p.second.first;
            *this = p.first;
        }
        else {

            *this = *da_dest; // UnitigColors associated with reference unitig of um_dest is recycled
            // Set the DataAccessor of um_dest to 0 (means inserted into overflow system of DataStorage)
            // to ensure the (recycled) UnitigColors is not deleted after concatenation
            *da_dest = DataAccessor<void>(0);
        }

        *uc_ptr = move(uc_);
    }
}

template<typename U>
void DataAccessor<U>::merge(const UnitigColorMap<U>& um_dest, const const_UnitigColorMap<U>& um_src){

    DataStorage<U>* ds = um_dest.getGraph()->getData(); // Get DataStorage where the UnitigColors and data are stored

    DataAccessor<U>* da_dest = um_dest.getData(); // Get DataAccessor stored with reference unitig of um_dest
    UnitigColors* uc_dest = da_dest->getUnitigColors(um_dest); // Get UnitigColors associated with reference unitig of um_dest
    U* data_dest = da_dest->getData(um_dest); // Get data associated with reference unitig of um_dest

    const DataAccessor<U>* da_src = um_src.getData(); // Get DataAccessor stored with reference unitig of um_src
    const UnitigColors* uc_src = da_src->getUnitigColors(um_src); // Get UnitigColors associated with reference unitig of um_src
    const U* data_src = da_src->getData(um_src); // Get data associated with reference unitig of um_src

    // If reference unitig of um_dest has no UnitigColors but reference unitig of um_src has one
    if ((uc_dest == nullptr) && (uc_src != nullptr)){

        UnitigColorMap<U> um(um_dest);

        um.dist = 0;
        um.len = um.size - um_dest.getGraph()->getK() + 1;
        um.strand = true;

        // Insert new UnitigColors and data for reference unitig of um_dest
        const pair<DataAccessor<U>, pair<UnitigColors*, U*>> p  = ds->insert(um);

        uc_dest = p.second.first; // Set the new UnitigColors
        data_dest = p.second.second; // Set the new data
        *this = p.first; // Set the new DataAccessor to locate new UnitigColors and data
    }

    // Merge colors for k-mers with positions matching mapping given in um_src
    if ((uc_dest != nullptr) && (uc_src != nullptr)) ds->addUnitigColors(um_dest, um_src);
    // Merge data for k-mers with positions matching mapping given in um_src
    if ((data_dest != nullptr) && (data_src != nullptr)) data_dest->merge(um_dest, um_src);
}

template<>
inline void DataAccessor<void>::merge(const UnitigColorMap<void>& um_dest, const const_UnitigColorMap<void>& um_src){

    DataStorage<void>* ds = um_dest.getGraph()->getData(); // Get DataStorage where the UnitigColors are stored

    UnitigColors* uc_dest = um_dest.getData()->getUnitigColors(um_dest); // Get UnitigColors associated with reference unitig of um_dest
    const UnitigColors* uc_src = um_src.getData()->getUnitigColors(um_src); // Get UnitigColors associated with reference unitig of um_src

    // If reference unitig of um_dest has no UnitigColors but reference unitig of um_src has one
    if ((uc_dest == nullptr) && (uc_src != nullptr)){

        UnitigColorMap<void> um(um_dest);

        um.dist = 0;
        um.len = um.size - um.getGraph()->getK() + 1;
        um.strand = true;

        // Insert new UnitigColors for reference unitig of um_dest
        const pair<DataAccessor<void>, pair<UnitigColors*, void*>> p  = ds->insert(um);

        uc_dest = p.second.first; // Set the new UnitigColors
        *this = p.first; // Set the new DataAccessor to locate new UnitigColors
    }

    // Merge colors for k-mers with positions matching mapping given in um_src
    if ((uc_dest != nullptr) && (uc_src != nullptr)) ds->addUnitigColors(um_dest, um_src);
}

template<typename U>
void DataAccessor<U>::extract(const UnitigColorMap<U>& um_src, const bool last_extraction) {

    DataStorage<U>* ds = um_src.getGraph()->getData(); // Get DataStorage where the UnitigColors and data are stored

    const pair<DataAccessor<U>, pair<UnitigColors*, U*>> p  = ds->insert(um_src); // Insert new UnitigColors + data associated with that mapping

    if (ds->getUnitigColors(um_src) != nullptr){ // If reference unitig of um_src has a UnitigColors associated

        UnitigColors new_cs = ds->getSubUnitigColors(um_src); // Extract colors for k-mer positions matching mapping of um_src

        if (!new_cs.isEmpty()) *(p.second.first) = move(new_cs); // Move extracted colors to new (inserted) UnitigColors
    }

    U* new_data = p.second.second; // Get new data slot

    if (new_data != nullptr) new_data->extract(p.second.first, um_src, last_extraction); // Extract new data

    *this = p.first; // Set the new DataAccessor to locate new UnitigColors and new data
}

template<>
inline void DataAccessor<void>::extract(const UnitigColorMap<void>& um_src, const bool last_extraction) {

    DataStorage<void>* ds = um_src.getGraph()->getData(); // Get DataStorage where the UnitigColors are stored

    if (ds->getUnitigColors(um_src) != nullptr){ // If reference unitig of um_src has a UnitigColors associated

        UnitigColors new_cs = ds->getSubUnitigColors(um_src); // Extract colors for k-mers with positions matching mapping given in um_src

        if (!new_cs.isEmpty()){ // If some colors were extracted

            const pair<DataAccessor<void>, pair<UnitigColors*, void*>> p  = ds->insert(um_src); // Insert new UnitigColors associated with that mapping

            *(p.second.first) = move(new_cs); // Move extracted colors to new (inserted) UnitigColors
            *this = p.first; // Set the new DataAccessor to locate new UnitigColors
        }
    }
}

template<typename U>
string DataAccessor<U>::serialize(const const_UnitigColorMap<U>& um_src) const {

    string da_str("DA:Z:" + std::to_string(da_id));

    const U* data_src = um_src.getData()->getData(um_src);

    if (data_src != nullptr) {

        const string data_str(data_src->serialize(um_src));

        if (!data_str.empty()) da_str += string('\t' + data_str);
    }

    return da_str;
}

template<>
inline string DataAccessor<void>::serialize(const const_UnitigColorMap<void>& um_src) const {

    return string("DA:Z:" + std::to_string(da_id));
}

#endif
