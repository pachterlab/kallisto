#ifndef NODE_HPP
#define NODE_HPP

#include <vector>
#include <unordered_map>
#include <set>

#include <CompactedDBG.hpp>

// Unitig to transcript
struct u2t {
    // Transcript ID
    uint32_t tr_id;
    // Position of unitig within transcript
    uint32_t pos;
    // Strandedness
    bool sense;
    u2t() {};
    u2t(int tr_id, int pos, bool sense) : tr_id(tr_id), pos(pos), sense(sense) {};
};

class Node: public CDBG_Data_t<Node> {
    public:
        int id;
        int ec;
        std::vector<u2t> transcripts;

    Node() : id(-1), ec(-1) {}

    void concat(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_dest = um_dest.getData();
        Node* data_src = um_src.getData();

        if (data_dest->ec != data_src->ec) {
            std::cout << "warning: concatenating two unitigs with corresponding";
            std::cout << " to different ECs" << std::endl;
        }
        ec = data_dest->ec;

        transcripts = data_dest->transcripts;
        transcripts.reserve(data_dest->transcripts.size() + data_src->transcripts.size());
        for (auto& trans : data_src->transcripts) {
            transcripts.push_back(trans);
            // Pos denotes where the unitig begins within the transcript. Since
            // we're prepending um_dest to um_src, the new unitig begins
            // um_dest.len + k - 1 base pairs earlier in the transcript
            transcripts[transcripts.size()-1].pos -= (um_dest.len + um_dest.getGraph()->getK() - 1);
        }
    }

    void merge(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_src = um_src.getData();

        if (ec != data_src->ec) {
            std::cout << "warning: merging two unitigs with corresponding";
            std::cout << " to different ECs" << std::endl;
        }


        transcripts.reserve(transcripts.size() + data_src->transcripts.size());
        for (auto& trans : data_src->transcripts) {
            transcripts.push_back(trans);
            // Pos denotes where the unitig begins within the transcript. Since
            // we're prepending um_dest to um_src, the new unitig begins
            // um_dest.len + k - 1 base pairs earlier in the transcript
            transcripts[transcripts.size()-1].pos -= (um_dest.len + um_dest.getGraph()->getK() - 1);
        }
    }

    void clear(const UnitigMap<Node>& um_dest) {
        transcripts.clear();
    }

    void extract(const UnitigMap<Node>& um_src, bool last_extraction) {
        if (ec == -1) {
          return;
        }

        Node* data = um_src.getData();
        data->ec = ec;
        data->transcripts = transcripts;

        for (auto tr : data->transcripts) {
            // We're cutting um_src.dist from the beginning of the unitig so
            // the new unitig begins um_src.dist further into back.
            tr.pos += um_src.dist;
        }
    }
};

#endif
