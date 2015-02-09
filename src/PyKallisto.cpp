#include "PyKallisto.h"

// counts gets overwritten with the counts per ecid
KmerIndex read_kal_py_index(const std::string& fname,
        const std::string& fasta,
        std::string& py_counts,
        std::vector<int>& counts) {

    KmerIndex kidx;

    // first read all the transcript names and lengths
    std::cout << "reading fasta" << std::endl;
    gzFile fp = gzopen(fasta.c_str(),"r");
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((l = kseq_read(seq)) > 0) {
      kidx.target_names_.push_back(seq->name.s);
      ++kidx.num_trans;
      kidx.trans_lens_.push_back(seq->seq.l);
    }
    kseq_destroy(seq);
    gzclose(fp);

    // next, read the python index

    std::cout << "reading python index" << std::endl;
    std::ifstream in;
    in.open(fname, std::ios::in);
    int max_ecid = -1;

    std::string line;
    while (std::getline(in, line)) {
        int ecid;
        std::istringstream iss(line);

        iss >> ecid;

        std::vector<int> ids;

        int tid;
        int cur_count;

        std::map<int, int> id_map;

        while (iss >> tid) {
            iss >> cur_count;

            if (tid % 2 == 0) {
                tid = tid >> 1;
            } else {
                tid = (tid - 1) >> 1;
            }

            auto search = id_map.find(tid);
            if (search != id_map.end()) {
                id_map[tid] += cur_count;
            } else {
                id_map.insert({tid, cur_count});
            }
        }

        for (auto kv : id_map) {
            ids.push_back(kv.first);
        }

        if (ids.size() == 0) {
            // empty equivalence class
            continue;
        }

        if (ecid > max_ecid) {
            max_ecid = ecid;
        }

        kidx.ecmap.insert({ecid, ids});

        // don't actually use inverse in py-em run
        // kidx.ecmapinv.insert({ids, ecid});
    }
    in.close();

    counts.clear();
    for (auto i = 0; i < max_ecid + 1; ++i) {
        counts.push_back(0);
    }

    std::ifstream py_counts_in;
    py_counts_in.open(py_counts, std::ios::in);

    // ignore first line
    std::getline(py_counts_in, line);

    while( std::getline(py_counts_in, line )) {
        int ecid;
        int cur_count;
        std::istringstream iss(line);

        iss >> ecid >> cur_count;
        counts[ecid] = cur_count;
    }
    py_counts_in.close();

    return kidx;
}
