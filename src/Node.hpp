#ifndef NODE_HPP
#define NODE_HPP

#include <vector>

#include <CompactedDBG.hpp>

// Unitig to transcript
struct u2t {
    // Transcript ID
    int tr_id;
    // Position of unitig within transcript
    int pos;
    // Strandedness
    bool sense;
    u2t() {};
    u2t(int tr_id, int pos, bool sense) : tr_id(tr_id), pos(pos), sense(sense) {};
};

class Node: public CDBG_Data_t<Node> {
    public:
        // Unitig ID
        int id;
        // Length of unitig
        int len;
        // Equivalence class of unitig
        std::vector<int> ec;
        std::vector<u2t> transcripts;

    void initialize_ec(int len) {
        ec = std::vector<int>(len, -1);
    }

    void concat(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_dest = um_dest.getData();
        Node* data_src = um_src.getData();

        ec.clear();
        ec = data_dest->ec;
        ec.reserve(um_dest.size + um_src.size);
        for (const auto& e : data_src->ec) {
            ec.push_back(e);
        }
    }

    void merge(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_src = um_src.getData();

        ec.reserve(um_dest.size + um_src.size);
        for (const auto& e : data_src->ec) {
            ec.push_back(e);
        }
    }

    void clear(const UnitigMap<Node>& um_dest) {
        ec.clear();
        transcripts.clear();
    }

    void extract(const UnitigMap<Node>& um_src, bool last_extraction) {
        ec.clear();
        ec.reserve(um_src.size);

        Node* data = um_src.getData();
        for (size_t i = um_src.dist; i < um_src.size; ++i) {
            ec.push_back(data->ec[i]);
        }
    }
};

#endif
