#ifndef NODE_HPP
#define NODE_HPP

#include <vector>
#include <unordered_map>
#include <set>

#include <CompactedDBG.hpp>
#include <BlockArray.hpp>

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
        // Mosaic Equivalence Class:
        // Each kmer in the unitig can have a different equivalence class
        BlockArray<uint32_t> ec;
        // EC : [transcripts]
        // Each equivalence class that is part of this unitig's mosaic
        // equivalence class can have a different list of transcripits
        // associated with it
        std::unordered_map<uint32_t, std::vector<u2t> > transcripts;

    Node() : id(-1) {}

    // Returns [j, k), j<=i, k>=i, where j-1 is the last kmer to have a
    // different EC than kmer i, and k is the first kmer to have a different
    // EC than i.
    std::pair<size_t, size_t> get_mc_contig(size_t i) const {
        return ec.get_block_at(i);
    }

    void concat(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_dest = um_dest.getData();
        Node* data_src = um_src.getData();

        ec = data_dest->ec;
        ec.reserve(ec.size() + data_src->transcripts.size());
        for (const auto& e : data_src->ec) {
            ec.insert(e.lb + um_dest.len, e.ub + um_dest.len, e.val);
        }

        // Pos denotes where the unitig begins within the transcript. Since
        // we're prepending um_dest to um_src, the new unitig begins
        // um_dest.len + k - 1 base pairs earlier in the transcript
        int offset = (um_dest.len + um_dest.getGraph()->getK() - 1);
        transcripts = data_dest->transcripts;
        for (const auto& kv : data_src->transcripts) {
            transcripts[kv.first] = kv.second;
            for (auto tr : transcripts[kv.first]) {
                tr.pos -= offset;
            }
        }
    }

    void merge(const UnitigMap<Node>& um_dest, const UnitigMap<Node>& um_src) {
        Node* data_src = um_src.getData();

        ec.reserve(ec.size() + data_src->transcripts.size());
        for (const auto& e : data_src->ec) {
            ec.insert(e.lb + um_dest.len, e.ub + um_dest.len, e.val);
        }

        int offset = (um_dest.len + um_dest.getGraph()->getK() - 1);
        for (const auto& kv : data_src->transcripts) {
            transcripts[kv.first] = kv.second;
            for (auto tr : transcripts[kv.first]) {
                // Pos denotes where the unitig begins within the transcript. Since
                // we're prepending um_dest to um_src, the new unitig begins
                // um_dest.len + k - 1 base pairs earlier in the transcript
                tr.pos -= offset;
            }
        }
    }

    void clear(const UnitigMap<Node>& um_dest) {
        ec.clear();
        transcripts.clear();
    }

    void extract(const UnitigMap<Node>& um_src, bool last_extraction) {
        if (ec.size() == 0) {
          return;
        }
        Node* data = um_src.getData();

        ec = data->ec.get_slice(um_src.dist, um_src.dist + um_src.len);

        std::vector<uint32_t> ecs;
        ec.get_vals(ecs);

        for (uint32_t i : ecs) {
            data->transcripts[i] = transcripts[i];
            for (auto tr : data->transcripts[i]) {
                // We're cutting um_src.dist from the beginning of the unitig  so
                // the new unitig begins um_src.dist further into back.
                tr.pos += um_src.dist;
            }
        }
    }

    void serialize(std::ofstream& out) const {

        size_t tmp_size;

        // 1 Write id
        out.write((char *)&id, sizeof(id));

        // 2 Write mosaic equivalence class
        ec.serialize(out);

        // 3 Write transcripts for each equivalence class
        tmp_size = transcripts.size();
        out.write((char *)&tmp_size, sizeof(tmp_size));
        for (const auto& kv : transcripts) {
            // 3.1 Write EC
            out.write((char *)&kv.first, sizeof(kv.first));
            // 3.2 Write transcripts
            tmp_size = kv.second.size();
            out.write((char *)&tmp_size, sizeof(tmp_size));
            for (const auto& tr : kv.second) {
                out.write((char *)&tr.tr_id, sizeof(tr.tr_id));
                out.write((char *)&tr.pos, sizeof(tr.pos));
                out.write((char *)&tr.sense, sizeof(tr.sense));
            }
        }
    }

    void deserialize(std::ifstream& in) {

        size_t tmp_size;
        uint32_t tmp_uint;

        // 1 Read id
        in.read((char *)&id, sizeof(id));

        // 2 Read mosaic equivalence class
        ec.deserialize(in);

        // 3 Read transcripts for each equivalence class
        size_t n_transcripts;
        uint32_t tr_id, pos;
        bool sense;
        in.read((char *)&tmp_size, sizeof(tmp_size));
        transcripts.reserve(tmp_size);
        for (size_t i = 0; i < tmp_size; ++i) {
            // 3.1 Read EC
            in.read((char *)&tmp_uint, sizeof(tmp_uint));
            transcripts[tmp_uint] = std::vector<u2t>();
            // 3.2 Read transcripts
            in.read((char *)&n_transcripts, sizeof(n_transcripts));
            transcripts[tmp_uint].reserve(n_transcripts);
            for (size_t j = 0; j < n_transcripts; ++j) {
                in.read((char *)&tr_id, sizeof(tr_id));
                in.read((char *)&pos, sizeof(pos));
                in.read((char *)&sense, sizeof(sense));
                transcripts[tmp_uint].emplace_back(tr_id, pos, sense);
            }
        }
    }
};

#endif
