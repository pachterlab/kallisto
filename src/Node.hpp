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
        int ec;
        std::vector<u2t> transcripts;

        /*
        inline int get_pos() const {
            return (_pos & 0x0FFFFFFF);
        }

        inline int is_fw() const {
            return (_pos & 0xF0000000) == 0;
        }

        inline int get_dist(bool fwd) const {
            if (is_fw() == fwd) {
                return (len - 1 - get_pos());
            }
            return get_pos();
        }
        */

};


#endif
