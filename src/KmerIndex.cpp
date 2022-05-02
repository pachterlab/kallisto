#include <algorithm>
#include <random>
#include <sstream>
#include <ctype.h>
#include <zlib.h>
#include <unordered_set>
#include "kseq.h"
#include "KmerIndex.h"

#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread)
#endif

// helper functions
// pre: u is sorted
bool isUnique(const std::vector<int>& u) {
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] == u[j]) {
      return false;
    }
  }
  return true;
}

std::vector<int> unique(const std::vector<int>& u) {
  std::vector<int> v;
  v.reserve(u.size());
  v.push_back(u[0]);
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] != u[j]) {
      v.push_back(u[j]);
    }
  }
  return v;
}

const char Dna(int i) {
  static const char *dna = "ACGT";
  return dna[i & 0x03];
}

int hamming(const char *a, const char *b) {
  int h = 0;
  while (*a != 0 && *b != 0) {
    if (*a != *b) {
      h++;
    }
    a++;
    b++;
  }
  return h;
}

std::string revcomp(const std::string s) {
  std::string r(s);
  std::transform(s.rbegin(), s.rend(), r.begin(), [](char c) {
      switch(c) {
      case 'A': return 'T';
      case 'C': return 'G';
      case 'G': return 'C';
      case 'T': return 'A';
      default: return 'N';
      }
      return 'N';
    });
  return r;
}

void KmerIndex::BuildTranscripts(const ProgramOptions& opt) {
  // read input
  std::unordered_set<std::string> unique_names;
  int k = opt.k;
  for (auto& fasta : opt.transfasta) {
    std::cerr << "[build] loading fasta file " << fasta
              << std::endl;
  }
  std::cerr << "[build] k-mer length: " << k << std::endl;

  std::vector<std::string> seqs;

  // read fasta file
  gzFile fp = 0;
  kseq_t *seq;
  int l = 0;
  std::mt19937 gen(42);
  int countNonNucl = 0;
  int countUNuc = 0;
  int polyAcount = 0;

  for (auto& fasta : opt.transfasta) {
    fp = gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);
    while (true) {
      l = kseq_read(seq);
      if (l <= 0) {
        break;
      }
      seqs.emplace_back(seq->seq.s);
      std::string& str = *seqs.rbegin();
      auto n = str.size();
      for (auto i = 0; i < n; i++) {
        char c = str[i];
        c = ::toupper(c);
        if (c=='U') {
          str[i] = 'T';
          countUNuc++;
        } else if (c !='A' && c != 'C' && c != 'G' && c != 'T') {
          str[i] = Dna(gen()); // replace with pseudorandom string
          countNonNucl++;
        }
      }
      std::transform(str.begin(), str.end(),str.begin(), ::toupper);

      if (str.size() >= 10 && str.substr(str.size()-10,10) == "AAAAAAAAAA") {
        // clip off polyA tail
        //std::cerr << "[index] clipping off polyA tail" << std::endl;
        polyAcount++;
        int j;
        for (j = str.size()-1; j >= 0 && str[j] == 'A'; j--) {}
        str = str.substr(0,j+1);
      }

      target_lens_.push_back(seq->seq.l);
      std::string name(seq->name.s);
      size_t p = name.find(' ');
      if (p != std::string::npos) {
        name = name.substr(0,p);
      }

      if (unique_names.find(name) != unique_names.end()) {
        if (!opt.make_unique) {
          std::cerr << "Error: repeated name in FASTA file " << fasta << "\n" << name << "\n\n" << "Run with --make-unique to replace repeated names with unique names" << std::endl;
          exit(1);
        } else {
          for (int i = 1; ; i++) { // potential bug if you have more than 2^32 repeated names
            std::string new_name = name + "_" + std::to_string(i);
            if (unique_names.find(new_name) == unique_names.end()) {
              name = new_name;
              break;
            }
          }
        }
      }
      unique_names.insert(name);
      target_names_.push_back(name);

    }
    gzclose(fp);
    fp=0;
  }

  if (polyAcount > 0) {
    std::cerr << "[build] warning: clipped off poly-A tail (longer than 10)" << std::endl << "        from " << polyAcount << " target sequences" << std::endl;
  }

  if (countNonNucl > 0) {
    std::cerr << "[build] warning: replaced " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl << "        with pseudorandom nucleotides" << std::endl;
  }
  if (countUNuc > 0) {
    std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
  }

  num_trans = seqs.size();

  // for each target, create its own equivalence class
  for (int i = 0; i < seqs.size(); i++ ) {
    std::vector<int> single(1,i);
    //ecmap.insert({i,single});
    ecmap.push_back(single);
    ecmapinv.insert({single,i});
  }

  BuildDeBruijnGraph(opt, seqs);
  BuildEquivalenceClasses(opt, seqs);
  //BuildEdges(opt);

}

void KmerIndex::BuildDeBruijnGraph(const ProgramOptions& opt, const std::vector<std::string>& seqs) {

  std::string tmp_file = ".tmp_dbg.fasta";
  std::ofstream of(tmp_file);
  for (size_t i = 0; i < seqs.size(); ++i) {
    of << ">" << std::to_string(i) << std::endl;
    of << seqs[i] << std::endl;
  }
  of.close();

  CDBG_Build_opt c_opt;
  c_opt.k = k;
  c_opt.nb_threads = opt.threads;
  c_opt.build = true;
  c_opt.clipTips = false;
  c_opt.deleteIsolated = false;
  c_opt.useMercyKmers = true;
  c_opt.verbose = true;
  c_opt.filename_ref_in.push_back(tmp_file);

  dbg = CompactedDBG<Node>(k);
  dbg.build(c_opt);

  std::remove(tmp_file.c_str());
}

void KmerIndex::BuildEquivalenceClasses(const ProgramOptions& opt, const std::vector<std::string>& seqs) {

  std::cerr << "[build] creating equivalence classes ... "; std::cerr.flush();

  std::vector<std::vector<TRInfo>> trinfos(dbg.size());
  UnitigMap<Node> um;
  uint32_t running_id = 0;
  //std::cout << "Mapping target " << std::endl;
  for (size_t i = 0; i < seqs.size(); ++i) {
    const auto& seq = seqs[i];
    if (seq.size() < k) continue;

    int seqlen = seq.size() - k + 1; // number of k-mers
    const char *s = seq.c_str();
    //std::cout << "sequence number " << i << std::endl;
    size_t proc = 0;
    while (proc < seqlen) {
      um = dbg.findUnitig(seq.c_str(), proc, seq.size());
      if (um.isEmpty) {
        ++proc;
        continue;
      }
      if (um.getData()->id == -1) um.getData()->id = running_id++;

      TRInfo tr;

      tr.trid = i;
      tr.pos = proc;
      tr.sense = um.strand;
      tr.start = um.dist;
      tr.stop  = um.dist + um.len;
      /*
      if (tr.sense) {
        tr.start = um.dist;
        // tr = um[start, stop)
        tr.stop  = um.dist + um.len;
      } else {
        tr.start = um.dist + um.len;
        // tr = um[start, stop)
        tr.stop  = um.dist + 1;
      }
      */

      trinfos[um.getData()->id].push_back(tr);

      proc += um.len;
    }
  }

  PopulateMosaicECs(trinfos);

  std::cerr << " done" << std::endl;
  std::cerr << "[build] target de Bruijn graph has " << dbg.size() << " contigs and contains "  << dbg.nbKmers() << " k-mers " << std::endl;
  std::cerr << "[build] target de Bruijn graph contains " << ecmapinv.size() << " equivalence classes from " << seqs.size() << " sequences." << std::endl;
}

void KmerIndex::PopulateMosaicECs(std::vector<std::vector<TRInfo> >& trinfos) {

  for (const auto& um : dbg) {

    Node* n = um.getData();
    // Mosaic EC vector should always be the length of the number of kmers
    // in the corresponding unitig
    n->ec.resize(um.size - k + 1);

    // Find the overlaps
    std::vector<int> brpoints;
    for (const auto& x : trinfos[n->id]) {
      brpoints.push_back(x.start);
      brpoints.push_back(x.stop);
    }
    sort(brpoints.begin(), brpoints.end());
    assert(brpoints[0] == 0);
    assert(brpoints[brpoints.size()-1]==um.size());

    // Find unique break points
    if (!isUnique(brpoints)) {
      std::vector<int> u = unique(brpoints);
      swap(u,brpoints);
    }

    // Create a mosaic EC for the unitig, where each break point interval
    // corresponds to one set of transcripts and therefore an EC
    for (size_t i = 1; i < brpoints.size(); ++i) {
      std::vector<int> u;
      std::vector<u2t> transcripts;
      for (const auto& tr : trinfos[n->id]) {
        // If a transcript encompasses the full breakpoint interval
        if (tr.start <= brpoints[i-1] && tr.stop >= brpoints[i]) {
          u.push_back(tr.trid);
          transcripts.emplace_back(tr.trid, tr.pos, tr.sense);
        }
      }

      sort(u.begin(), u.end());
      if (!isUnique(u)){
        std::vector<int> v = unique(u);
        swap(u,v);
      }

      assert(!u.empty());

      auto search = ecmapinv.find(u);
      int ec = -1;
      if (search != ecmapinv.end()) {
        // See if we've created an EC for this set of transcripts yet
        ec = search->second;
      } else {
        // Otherwise, we create a new EC
        ec = ecmapinv.size();
        ecmapinv.insert({u,ec});
        ecmap.push_back(u);
      }

      // Assign mosaic EC to the corresponding part of unitig
      for (size_t j = brpoints[i-1]; j < brpoints[i]; ++j) {
        n->ec[j] = ec;
      }
      // Assign the transcripts to the EC
      n->transcripts[ec] = transcripts;
    }
  }
}

void KmerIndex::write(const std::string& index_out, bool writeKmerTable) {
  std::ofstream out;
  out.open(index_out, std::ios::out | std::ios::binary);

  if (!out.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index output file could not be opened!";
    exit(1);
  }

  // 1. write version
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. write k
  out.write((char *)&k, sizeof(k));

  // 3. write number of targets
  out.write((char *)&num_trans, sizeof(num_trans));

  // 4. write out target lengths
  for (int tlen : target_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  // 5. write fake k-mer size
  size_t kmap_size = 0;
  out.write((char *)&kmap_size, sizeof(kmap_size));

  // 6. write none of the kmer->ec values

  // 7. write number of equivalence classes
  size_t tmp_size;
  tmp_size = ecmap.size();
  out.write((char *)&tmp_size, sizeof(tmp_size));

  // 8. write out each equiv class
  for (int ec = 0; ec < ecmap.size(); ec++) {
    out.write((char *)&ec, sizeof(ec));
    auto& v = ecmap[ec];
    // 8.1 write out the size of equiv class
    tmp_size = v.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));
    // 8.2 write each member
    for (auto& val: v) {
      out.write((char *)&val, sizeof(val));
    }
  }

  // 9. Write out target ids
  // XXX: num_trans should equal to target_names_.size(), so don't need
  // to write out again.
  assert(num_trans == target_names_.size());
  for (auto& tid : target_names_) {
    // 9.1 write out how many bytes
    // XXX: Note: this doesn't actually encore the max targ id size.
    // might cause problems in the future
    // tmp_size = tid.size();
    tmp_size = strlen(tid.c_str());
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 9.2 write out the actual string
    out.write(tid.c_str(), tmp_size);
  }

  // 10. write out contigs
  if (writeKmerTable) {
    // Split the unitigs in the graph into vanilla kallisto contigs that have
    // non-mosaic ECs
    tmp_size = 0;
    for (const auto& um : dbg) {
      // The key in the transcripts unordered_map is the EC, and the number of
      // ECs in a unitig denote how many vanilla kallisto contigs they get split
      // into.
      tmp_size += um.getData()->transcripts.size();
    }

    out.write((char*)&tmp_size, sizeof(tmp_size));
    // Store the ECs in the same order as we write out the contigs
    std::vector<uint32_t> ecs;
    for (const auto& um : dbg) {
      // I beg forgiveness for the double pass.
      Node* n = um.getData();
      uint32_t curr_ec = n->ec[0];
      std::string ref = um.referenceUnitigToString();
      size_t start = 0;
      for (size_t i = 0; i < n->ec.size(); ++i) {
        if (curr_ec == n->ec[i]) continue;
        ecs.push_back(curr_ec);
        std::string contig = ref.substr(start, i-start+k);
        out.write((char*)&n->id, sizeof(n->id));
        tmp_size = strlen(contig.c_str());
        out.write((char*)&tmp_size, sizeof(tmp_size));
        out.write(contig.c_str(), tmp_size);

        // 10.1 write out transcript info
        tmp_size = n->transcripts[curr_ec].size();
        out.write((char*)&tmp_size, sizeof(tmp_size));
        for (const auto& info : n->transcripts[curr_ec]) {
          out.write((char*)&info.tr_id, sizeof(info.tr_id));
          out.write((char*)&info.pos, sizeof(info.pos));
          out.write((char*)&info.sense, sizeof(info.sense));
        }

        start = i;
        curr_ec = n->ec[i];
      }
      ecs.push_back(curr_ec);
      std::string contig = ref.substr(start, ref.length()-start);
      out.write((char*)&n->id, sizeof(n->id));
      tmp_size = strlen(contig.c_str());
      out.write((char*)&tmp_size, sizeof(tmp_size));
      out.write(contig.c_str(), tmp_size);

      // 10.1 write out transcript info for last mosaic EC in unitig
      tmp_size = n->transcripts[curr_ec].size();
      out.write((char*)&tmp_size, sizeof(tmp_size));
      for (const auto& info : n->transcripts[curr_ec]) {
        out.write((char*)&info.tr_id, sizeof(info.tr_id));
        out.write((char*)&info.pos, sizeof(info.pos));
        out.write((char*)&info.sense, sizeof(info.sense));
      }
    }

    // 11. write out ecs info
    for (uint32_t ec : ecs) {
      out.write((char*)&ec, sizeof(ec));
    }

  } else {
    // write empty dBG
    tmp_size = 0;
    out.write((char*)&tmp_size, sizeof(tmp_size));
  }

  out.flush();
  out.close();
}

bool KmerIndex::fwStep(Kmer km, Kmer& end) const {
  /*
  int j = -1;
  int fw_count = 0;
  for (int i = 0; i < 4; i++) {
    Kmer fw_rep = end.forwardBase(Dna(i)).rep();
    auto search = kmap.find(fw_rep);
    if (search != kmap.end()) {
      j = i;
      ++fw_count;
      if (fw_count > 1) {
        return false;
      }
    }
  }

  if (fw_count != 1) {
    return false;
  }

  Kmer fw = end.forwardBase(Dna(j));

  int bw_count = 0;
  for (int i = 0; i < 4; i++) {
    Kmer bw_rep = fw.backwardBase(Dna(i)).rep();
    if (kmap.find(bw_rep) != kmap.end()) {
      ++bw_count;
      if (bw_count > 1) {
        return false;
      }
    }
  }

  if (bw_count != 1) {
    return false;
  } else {
    if (fw != km) {
      end = fw;
      return true;
    } else {
      return false;
    }
  }
  */
  // TODO: Remove
  return false;
}

void KmerIndex::load(ProgramOptions& opt, bool loadKmerTable) {
  if (opt.index.empty() && !loadKmerTable) {
    // Make an index from transcript and EC files
    loadTranscriptsFromFile(opt);
    loadECsFromFile(opt);
    return;
  }

  std::string& index_in = opt.index;
  std::ifstream in;

  in.open(index_in, std::ios::in | std::ios::binary);

  if (!in.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index input file could not be opened!";
    exit(1);
  }

  // 1. read version
  size_t header_version = 0;
  in.read((char *)&header_version, sizeof(header_version));

  if (header_version != INDEX_VERSION) {
    std::cerr << "Error: incompatible indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
              << "Rerun with index to regenerate";
    exit(1);
  }

  // 2. read k
  in.read((char *)&k, sizeof(k));
  dbg = CompactedDBG<Node>(k);
  if (Kmer::k == 0) {
    //std::cerr << "[index] no k has been set, setting k = " << k << std::endl;
    Kmer::set_k(k);
    opt.k = k;
  } else if (Kmer::k == k) {
    //std::cerr << "[index] Kmer::k has been set and matches" << k << std::endl;
    opt.k = k;
  } else {
    std::cerr << "Error: Kmer::k was already set to = " << Kmer::k << std::endl
              << "       conflicts with value of k  = " << k << std::endl;
    exit(1);
  }

  // 3. read in number of targets
  in.read((char *)&num_trans, sizeof(num_trans));

  // 4. read in length of targets
  target_lens_.clear();
  target_lens_.reserve(num_trans);

  for (int i = 0; i < num_trans; i++) {
    int tlen;
    in.read((char *)&tlen, sizeof(tlen));
    target_lens_.push_back(tlen);
  }

  // 5. read number of k-mers
  size_t kmap_size;
  in.read((char *)&kmap_size, sizeof(kmap_size));

  std::cerr << "[index] k-mer length: " << k << std::endl;
  std::cerr << "[index] number of targets: " << pretty_num(num_trans)
    << std::endl;
  std::cerr << "[index] number of k-mers: " << pretty_num(kmap_size)
    << std::endl;

  // 6. read kmer->ec values
  Kmer tmp_kmer;
  KmerEntry tmp_val;
  for (size_t i = 0; i < kmap_size; ++i) {
    // Spool over kmap. Don't need it.
    in.read((char *)&tmp_kmer, sizeof(tmp_kmer));
    in.read((char *)&tmp_val, sizeof(tmp_val));
  }

  // 7. read number of equivalence classes
  size_t ecmap_size;
  in.read((char *)&ecmap_size, sizeof(ecmap_size));

  std::cerr << "[index] number of equivalence classes: "
    << pretty_num(ecmap_size) << std::endl;
  ecmap.resize(ecmap_size);
  int tmp_id;
  int tmp_ecval;
  size_t vec_size;

  // 8. read each equiv class
  for (size_t ec = 0; ec < ecmap_size; ++ec) {
    in.read((char *)&tmp_id, sizeof(tmp_id));

    // 8.1 read size of equiv class
    in.read((char *)&vec_size, sizeof(vec_size));

    // 8.2 read each member
    std::vector<int> tmp_vec;
    tmp_vec.reserve(vec_size);
    for (size_t j = 0; j < vec_size; ++j ) {
      in.read((char *)&tmp_ecval, sizeof(tmp_ecval));
      tmp_vec.push_back(tmp_ecval);
    }
    //ecmap.insert({tmp_id, tmp_vec});
    ecmap[tmp_id] = tmp_vec;
    ecmapinv.insert({tmp_vec, tmp_id});
  }
  std::cout << "read ECs" << std::endl;

  // 9. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

  std::cout << "read target ids" << std::endl;
  size_t tmp_size;
  size_t bufsz = 1024;
  char *buffer = new char[bufsz];
  for (auto i = 0; i < num_trans; ++i) {
    // 9.1 read in the size
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size +1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }

    // clear the buffer
    memset(buffer,0,bufsz);
    // 9.2 read in the character string
    in.read(buffer, tmp_size);

    target_names_.push_back(std::string(buffer));
  }
  std::cout << "read target names" << std::endl;

  // 10. read contigs
  size_t contig_size;
  in.read((char *)&contig_size, sizeof(contig_size));

  int c_len, c_id;
  std::string c_seq;

  UnitigMap<Node> um;
  std::vector<std::vector<u2t> > transcripts(contig_size);
  Kmer kmer;
  Node* data;
  // Keep track of the order the contigs are listed in the index in order to
  // match them up with the ECs later on
  // Keep track of the original strandedness of the contig in order to fix
  // transcript.sense if the strandedness changes. (Because the strandedness
  // of Bifrost unitigs is not guaranteed to be the same as the strandedness
  // of kallisto contigs).
  // TODO:
  // See if we can get away with not keeping all the contig sequences in this
  // vector by restructuring the index.
  std::vector<std::string> canonical_contigs;
  canonical_contigs.reserve(contig_size);
  for (size_t i = 0; i < contig_size; ++i) {
    in.read((char *)&c_id, sizeof(c_id));
    // TODO:
    // We never use c_len. Remove it.
    //in.read((char *)&c_len, sizeof(c_len));
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size + 1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }

    memset(buffer,0,bufsz);
    in.read(buffer, tmp_size);
    c_seq = std::string(buffer); // copy
    canonical_contigs.push_back(c_seq);

    //dbg.add(c_seq, false);

    // 10.1 read transcript info
    in.read((char*)&tmp_size, sizeof(tmp_size));
    int tr_id, pos;
    bool sense;
    for (size_t j = 0; j < tmp_size; ++j) {
      in.read((char*)&tr_id, sizeof(tr_id));
      in.read((char*)&pos, sizeof(pos));
      in.read((char*)&sense, sizeof(sense));
      transcripts[i].emplace_back(tr_id, pos, sense);
    }
  }

  // TODO:
  // We read all the contigs in the index into a vector and write them out into
  // a temporary fasta file, from which we build a Bifrost CompactedDBG.
  // In the future, we should implement a method in Bifrost that can take in
  // a vector of strings and create a CompactedDBG from it, so we don't have
  // to write to disk.
  std::string tmp_file = ".tmp_dbg.fasta";
  std::ofstream of(tmp_file);
  for (size_t i = 0; i < canonical_contigs.size(); ++i) {
    of << ">" << std::to_string(i) << std::endl;
    of << canonical_contigs[i] << std::endl;
  }
  of.close();

  CDBG_Build_opt c_opt;
  c_opt.k = k;
  c_opt.nb_threads = opt.threads;
  c_opt.build = true;
  c_opt.clipTips = false;
  c_opt.deleteIsolated = false;
  c_opt.useMercyKmers = true;
  c_opt.filename_ref_in.push_back(tmp_file);

  dbg.build(c_opt);

  std::remove(tmp_file.c_str());

  // 11. read ecs info
  // TODO:
  // Change to uint32_t
  int tmp_ec;
  size_t running_id = 0;
  for (size_t i = 0; i < canonical_contigs.size(); ++i) {
    in.read((char*)&tmp_ec, sizeof(tmp_ec));
    size_t proc = 0;
    while (proc < canonical_contigs[i].length() - k + 1) {
      um = dbg.findUnitig(canonical_contigs[i].c_str(), proc, canonical_contigs[i].length());
      data = um.getData();
      if (data->ec.size() < um.size - k + 1) {
          data->ec.resize(um.size - k + 1);
      }
      for (size_t j = um.dist; j < um.dist + um.len; ++j) {
        data->ec[j] = tmp_ec;
      }
      proc += um.len;
      data->transcripts[tmp_ec] = transcripts[i];
      data->id = running_id++;

      // With Bifrost we have no control over the strandedness of the contigs.
      // We iterate over all the contigs we read from the index and check whether
      // the corresponding Bifrost unitig has the same strandedness. If not, we
      // flip the transcript.sense for all transcripts belonging to that unitig.
      if (!um.isEmpty && !um.strand) {
        // If the canonical contig sequence is not a substring of the found
        // unitig sequence
        for (auto& tr: data->transcripts[tmp_ec]) {
          tr.sense = !tr.sense;
        }
      }
    }
  }

  int polychromatic = 0;
  for (auto& um : dbg) {
    Node* n = um.getData();
    uint32_t ec = n->ec[0];
    for (size_t i = 0; i < n->ec.size(); ++i) {
      if (n->ec[i] != ec) {
        n->monochrome = false;
        break;
        polychromatic++;
      }
    }
  }
  std::cout << "There are " << polychromatic << " polychromatic unitigs in the graph out of " << dbg.size() << " unitigs (" << float(polychromatic) / dbg.size() << ")." << std::endl;

  // delete the buffer
  delete[] buffer;
  buffer=nullptr;

  in.close();

  if (!opt.ecFile.empty()) {
    loadECsFromFile(opt);
  }
}

void KmerIndex::loadECsFromFile(const ProgramOptions& opt) {
  ecmap.clear();
  ecmapinv.clear();
  int i = 0;
  std::ifstream in((opt.ecFile));
  if (in.is_open()) {
    std::string line;
    while (getline(in, line)) {
      std::stringstream ss(line);
      int ec;
      std::string transcripts;
      ss >> ec >> transcripts;
      if (i != ec) {
        std::cerr << "Error: equivalence class file has a misplaced equivalence class."
                  << " Found " << ec << ", expected " << i << std::endl;
        exit(1);
      }
      std::vector<int> tmp_vec;
      std::stringstream ss2(transcripts);
      while(ss2.good()) {
        std::string tmp_ecval;
        getline(ss2, tmp_ecval, ',');
        int tmp_ecval_num = std::atoi(tmp_ecval.c_str());
        if (tmp_ecval_num < 0 || tmp_ecval_num >= num_trans) {
          std::cerr << "Error: equivalence class file has invalid value: "
                    << tmp_ecval << " in " << transcripts << std::endl;
          exit(1);
        }
        tmp_vec.push_back(tmp_ecval_num);
      }
      ecmap.push_back(tmp_vec); // copy
      ecmapinv.insert({std::move(tmp_vec), i}); // move
      i++;
    }
  } else {
    std::cerr << "Error: could not open file " << opt.ecFile << std::endl;
    exit(1);
  }
  std::cerr << "[index] number of equivalence classes loaded from file: "
            << pretty_num(ecmap.size()) << std::endl;
}

void KmerIndex::loadTranscriptsFromFile(const ProgramOptions& opt) {
  target_names_.clear();
  int i = 0;
  std::ifstream in((opt.transcriptsFile));
  if (in.is_open()) {
    std::string txp;
    while (in >> txp) {
      target_names_.push_back(txp);
      i++;
    }
  } else {
    std::cerr << "Error: could not open file " << opt.transcriptsFile << std::endl;
    exit(1);
  }
  num_trans = i;
  target_lens_.assign(num_trans, 0);
  std::cerr << "[index] number of targets loaded from file: "
            << pretty_num(num_trans) << std::endl;
}

int KmerIndex::mapPair(const char *s1, int l1, const char *s2, int l2, int ec) const {
  bool d1 = true;
  bool d2 = true;
  int p1 = -1;
  int p2 = -1;
  int c1 = -1;
  int c2 = -1;

  KmerIterator kit1(s1), kit_end;
  const_UnitigMap<Node> um1, um2;

  bool found1 = false;
  for (; kit1 != kit_end; ++kit1) {
    um1 = dbg.find(kit1->first);
    if (!um1.isEmpty) {
      found1 = true;
      if (um1.strand)
        p1 = um1.dist - kit1->second;
      else
        p1 = um1.dist + k + kit1->second;
      break;
    }
  }

  if (!found1) {
    return -1;
  }

  KmerIterator kit2(s2);
  bool found2 = false;

  for (; kit2 != kit_end; ++kit2) {
    um2 = dbg.find(kit2->first);
    if (!um2.isEmpty) {
      found2 = true;
      if (um2.strand)
        p2 = um2.dist - kit2->second;
      else
        p2 = um2.dist + k + kit2->second;
      break;
    }
  }

  if (!found2) {
    return -1;
  }

  // We want the reads to map within the same EC block on the same unitig
  if (!um1.isSameReferenceUnitig(um2) ||
      um1.getData()->ec[um1.dist] != um2.getData()->ec[um2.dist]) {
    return -1;
  }

  // Paired reads need to map to opposite strands
  if (!(um1.strand ^ um2.strand)) {
    //std::cerr << "Reads map to same strand " << s1 << "\t" << s2 << std::endl;
    return -1;
  }

  if (p1>p2) {
    return p1-p2;
  } else {
    return p2-p1;
  }
}

// use:  match(s,l,v)
// pre:  v is initialized
// post: v contains all equiv classes for the k-mers in s
void KmerIndex::match(const char *s, int l, std::vector<std::pair<const_UnitigMap<Node>, int>>& v) const{
  const Node* n;

  // TODO:
  // Rework KmerIndex::match() such that it uses the following type of logic
  // rather than the jumping logic below

  /*
  size_t proc = 0;
  while (proc < l - k + 1) {
    const_UnitigMap<Node> um = dbg.findUnitig(s, proc, l);
    if (um.isEmpty) {
      proc++;
      continue;
    }

    n = um.getData();
    uint32_t curr_ec = n->ec[um.dist];
    v.emplace_back(um, proc);
    // Add one entry to v for each EC that is part of the mosaic EC of the contig.
    for (size_t i = 0; i < um.len; ++i) {
      if (n->ec[um.dist + i] != curr_ec) {
        curr_ec = n->ec[um.dist + i];
        v.emplace_back(dbg.find(um.getUnitigKmer(um.dist + i)), proc + i);
      }
    }
    proc += um.len;
  }
  */

  KmerIterator kit(s), kit_end;
  bool backOff = false;
  int nextPos = 0; // nextPosition to check
  for (int i = 0;  kit != kit_end; ++i,++kit) {
    // need to check it
    const_UnitigMap<Node> um = dbg.find(kit->first);
    n = um.getData();

    int pos = kit->second;

    if (!um.isEmpty) {

      v.push_back({um, kit->second});

      // Find start and end of O.G. kallisto contig w.r.t. the bifrost-kallisto
      // unitig
      size_t contig_start = 0, contig_length = um.size - k + 1;
      if (!n->monochrome) {
        auto p = n->get_mc_contig(um.dist);
        contig_start += p.first;
        contig_length = p.second - contig_start;
      }
      // Looks like kallisto thinks that canonical kmer means forward strand?
      //bool forward = (um.strand == (kit->first == kit->first.rep()));
      bool forward = um.strand;
      int dist = (forward) ? (contig_length - 1 - (um.dist - contig_start)) : um.dist - contig_start;

      // see if we can skip ahead
      if (dist >= 2) {
        // where should we jump to?
        int nextPos = pos+dist; // default jump

        if (pos + dist >= l-k) {
          // if we can jump beyond the read, check the end
          nextPos = l-k;
        }

        // check next position
        KmerIterator kit2(kit);
        kit2 += nextPos-pos;
        if (kit2 != kit_end) {
          const_UnitigMap<Node> um2 = dbg.find(kit2->first);
          bool found2 = false;
          int  found2pos = pos+dist;
          if (um2.isEmpty) {
            found2=true;
            found2pos = pos;
          } else if (um.isSameReferenceUnitig(um2) &&
                     n->ec[um.dist] == um2.getData()->ec[um2.dist]) {
            // um and um2 are on the same unitig and also share the same EC
            found2=true;
            found2pos = pos+dist;
          }
          if (found2) {
            // great, a match (or nothing) see if we can move the k-mer forward
            if (found2pos >= l-k) {
              v.push_back({um, l-k}); // push back a fake position
              break; //
            } else {
              v.push_back({um, found2pos});
              kit = kit2; // move iterator to this new position
            }
          } else {
            // this is weird, let's try the middle k-mer
            bool foundMiddle = false;
            if (dist > 4) {
              int middlePos = (pos + nextPos)/2;
              int middleContig = -1;
              int found3pos = pos+dist;
              KmerIterator kit3(kit);
              kit3 += middlePos-pos;
              if (kit3 != kit_end) {
                const_UnitigMap<Node> um3 = dbg.find(kit3->first);
                if (!um3.isEmpty) {
                  if (um.isSameReferenceUnitig(um3) &&
                      n->ec[um.dist] == um3.getData()->ec[um3.dist]) {
                    foundMiddle = true;
                    found3pos = middlePos;
                  } else if (um2.isSameReferenceUnitig(um3) &&
                             um2.getData()->ec[um2.dist] == um3.getData()->ec[um3.dist]) {
                    foundMiddle = true;
                    found3pos = pos+dist;
                  }
                }

                if (foundMiddle) {
                  v.push_back({um3, found3pos});
                  if (nextPos >= l-k) {
                    break;
                  } else {
                    kit = kit2;
                  }
                }
              }
            }

            if (!foundMiddle) {
              ++kit;
              backOff = true;
              goto donejumping; // sue me Dijkstra!
            }
          }
        } else {
          // the sequence is messed up at this point, let's just take the match
          //v.push_back({dbGraph.ecs[val.contig], l-k});
          break;
        }
      }
    }

donejumping:

    if (backOff) {
      // backup plan, let's play it safe and search incrementally for the rest, until nextStop
      for (int j = 0; kit != kit_end; ++kit,++j) {
        if (j==skip) {
          j=0;
        }
        if (j==0) {
          // need to check it
          const_UnitigMap<Node> um4 = dbg.find(kit->first);
          if (!um4.isEmpty) {
            // if k-mer found
            v.push_back({um4, kit->second}); // add equivalence class, and position
          }
        }

        if (kit->second >= nextPos) {
          backOff = false;
          break; // break out of backoff for loop
        }
      }
    }
  }
}

std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, int p) const{
  const_UnitigMap<Node> um = dbg.find(km);
  if (!um.isEmpty) {
    return findPosition(tr, km, um, p);
  } else {
    return {-1,true};
  }
}

//use:  (pos,sense) = index.findPosition(tr,km,val,p)
//pre:  index.kmap[km] == val,
//      km is the p-th k-mer of a read
//      val.contig maps to tr
//post: km is found in position pos (1-based) on the sense/!sense strand of tr
std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, const_UnitigMap<Node>& um, int p) const {
  bool csense = um.strand;

  int trpos = -1;
  bool trsense = true;
  if (um.getData()->id < 0) {
    return {-1, true};
  }
  const Node* n = um.getData();
  const auto trs = n->transcripts.find(n->ec[um.dist]);
  for (const auto& x : trs->second) {
    if (x.tr_id == tr) {
      trpos = x.pos;
      trsense = x.sense;
      break;
    }
  }

  if (trpos == -1) {
    return {-1,true};
  }

  // TODO:
  // Check whether um.dist does the same thing as KmerEntry.getPos();
  // If something doesn't work, it's most likely this!
  if (trsense) {
    if (csense) {
      return {trpos + um.dist - p + 1, csense}; // 1-based, case I
    } else {
      return {trpos + um.dist + k + p, csense}; // 1-based, case III
    }
  } else {
    if (csense) {
      // UnitigMapBase.size is the length in base pairs
      return {trpos + (um.size - um.dist - 1) + k + p, !csense};  // 1-based, case IV
    } else {
      return {trpos + (um.size - um.dist) - p, !csense}; // 1-based, case II
    }
  }
}

// use:  res = intersect(ec,v)
// pre:  ec is in ecmap, v is a vector of valid targets
//       v is sorted in increasing order
// post: res contains the intersection  of ecmap[ec] and v sorted increasing
//       res is empty if ec is not in ecma
std::vector<int> KmerIndex::intersect(int ec, const std::vector<int>& v) const {
  std::vector<int> res;
  //auto search = ecmap.find(ec);
  if (ec < ecmap.size()) {
    //if (search != ecmap.end()) {
    //auto& u = search->second;
    auto& u = ecmap[ec];
    res.reserve(v.size());

    auto a = u.begin();
    auto b = v.begin();

    while (a != u.end() && b != v.end()) {
      if (*a < *b) {
        ++a;
      } else if (*b < *a) {
        ++b;
      } else {
        // match
        res.push_back(*a);
        ++a;
        ++b;
      }
    }
  }
  return res;
}

void KmerIndex::loadTranscriptSequences() const {
  if (target_seqs_loaded) {
    return;
  }

  std::vector<std::vector<std::pair<std::string, u2t> > > trans_contigs(num_trans);
  for (auto& um : dbg) {
    std::cout << um.referenceUnitigToString() << std::endl;
    const Node* data = um.getData();
    std::string um_seq = um.mappedSequenceToString();
    if (!um.strand) {
      um_seq = revcomp(um_seq);
    }
    const auto trs = data->transcripts.find(data->ec[um.dist]);
    for (const auto& tr : trs->second) {
      // TODO:
      // Verify that this method of choosing to reverse complement is legit
      if (um.strand ^ tr.sense) {
        trans_contigs[tr.tr_id].push_back({revcomp(um_seq), tr});
      } else {
        trans_contigs[tr.tr_id].push_back({um_seq, tr});
      }
    }
  }
  std::cout << "====================" << std::endl;
  Kmer km("AGATGAT");
  auto um = dbg.find(km);
  std::cout << um.getMappedHead().toString() << std::endl;
  std::cout << um.getUnitigHead().toString() << std::endl;
  std::cout << um.strand << std::endl;
  std::cout << um.dist << std::endl;
  std::cout << "---" << std::endl;

  km = Kmer("ATCATCT");
  um = dbg.find(km);
  std::cout << um.getMappedHead().toString() << std::endl;
  std::cout << um.getUnitigHead().toString() << std::endl;
  std::cout << um.strand << std::endl;
  std::cout << um.dist << std::endl;
  std::cout << (km == um.getMappedHead()) << std::endl;

  auto &target_seqs = const_cast<std::vector<std::string>&>(target_seqs_);

  for (size_t i = 0; i < trans_contigs.size(); ++i) {
    auto& v = trans_contigs[i];
    std::sort(v.begin(),
              v.end(),
              [](std::pair<std::string, u2t> a, std::pair<std::string, u2t> b) {
                return a.second.pos < b.second.pos;
              });

    std::string seq;
    seq.reserve(target_lens_[i]);
    for (const auto& pct : v) {
      int start = (pct.second.pos == 0) ? 0 : k-1;
      seq.append(pct.first.substr(start));
    }
    target_seqs.push_back(seq);
  }

  bool &t = const_cast<bool&>(target_seqs_loaded);
  t = true;//target_seqs_loaded = true;
  return;
}

void KmerIndex::clear() {
  /*
  kmap.clear_table();
  ecmap.resize(0);
  dbGraph.ecs.resize(0);
  dbGraph.contigs.resize(0);
  {
    std::unordered_map<std::vector<int>, int, SortedVectorHasher> empty;
    std::swap(ecmapinv, empty);
  }
  
  target_lens_.resize(0);
  target_names_.resize(0);
  target_seqs_.resize(0);
  */
}

void KmerIndex::writePseudoBamHeader(std::ostream &o) const {
  /*
  // write out header
  o << "@HD\tVN:1.0\n";
  for (int i = 0; i < num_trans; i++) {
    o << "@SQ\tSN:" << target_names_[i] << "\tLN:" << target_lens_[i] << "\n";
  }
  o << "@PG\tID:kallisto\tPN:kallisto\tVN:"<< KALLISTO_VERSION << "\n";
  o.flush();
  */
}
