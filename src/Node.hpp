#ifndef NODE_HPP
#define NODE_HPP

#include <vector>
#include <unordered_map>
#include <set>

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
        // Equivalence class of unitig
        std::vector<int> ec;
        std::unordered_map<int, std::vector<u2t> > transcripts;

    void initialize_ec(int len) {
        ec = std::vector<int>(len, -1);
    }

    void concat(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_dest = um_dest.getData();
        Node* data_src = um_src.getData();

        ec = data_dest->ec;
        ec.reserve(data_dest->ec.size() + data_src->ec.size());
        for (const auto& e : data_src->ec) {
            ec.push_back(e);
        }

        transcripts = data_dest->transcripts;
        for (const auto& kv : data_src->transcripts) {
            transcripts[kv.first] = kv.second;
        }
    }

    void merge(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_src = um_src.getData();

        ec.reserve(ec.size() + data_src->ec.size());
        for (const auto& e : data_src->ec) {
            ec.push_back(e);
        }

        for (const auto& kv : data_src->transcripts) {
            transcripts[kv.first] = kv.second;
        }
    }

    void clear(const UnitigMap<Node>& um_dest) {
        ec.clear();
        transcripts.clear();
    }

    void extract(const UnitigMap<Node>& um_src, bool last_extraction) {
        Node* data = um_src.getData();
        data->ec.reserve(um_src.size);

        std::set<int> ecs;
        for (size_t i = um_src.dist; i < um.dist + um_src.len; ++i) {
            data->ec.push_back(data->ec[i]);
            ecs.insert(data->ec[i]);
        }

        for (const auto& i : ecs) {
            data->transcripts[i] = transcripts[i];
        }
    }
};

#endif
