#include <algorithm>
#include <random>
#include <sstream>
#include <ctype.h>
#include <zlib.h>
#include <unordered_set>
#include <functional>
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

std::string generate_tmp_file(std::string seed) {
  std::string base = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  std::string tmp_file = ".kallisto.";
  srand((unsigned int)std::hash<std::string>{}(seed));
  int pos;
  while(tmp_file.length() < 32) {
    pos = ((rand() % (base.size() - 1)));
    tmp_file += base.substr(pos, 1);
  }
  return tmp_file;
}

void KmerIndex::BuildTranscripts(const ProgramOptions& opt, std::ofstream& out) {
  // read input
  std::unordered_set<std::string> unique_names;
  int k = opt.k;
  for (auto& fasta : opt.transfasta) {
    std::cerr << "[build] loading fasta file " << fasta
              << std::endl;
  }
  std::cerr << "[build] k-mer length: " << k << std::endl;

  // Generate random file name
  std::string tmp_file = generate_tmp_file(opt.index);
  std::ofstream of(tmp_file);
  num_trans = 0;

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
      std::string str = seq->seq.s;
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
      of << ">" << num_trans++ << "\n" << str << std::endl;

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

  of.close();

  if (polyAcount > 0) {
    std::cerr << "[build] warning: clipped off poly-A tail (longer than 10)" << std::endl << "        from " << polyAcount << " target sequences" << std::endl;
  }

  if (countNonNucl > 0) {
    std::cerr << "[build] warning: replaced " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl << "        with pseudorandom nucleotides" << std::endl;
  }
  if (countUNuc > 0) {
    std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
  }

  BuildDeBruijnGraph(opt, tmp_file, out);
  BuildEquivalenceClasses(opt, tmp_file);

  /*
  for (auto& um : dbg) {
    auto* node = um.getData();
    node->ec.shrink_to_size();
    node->pos.shrink_to_size();
  }
  */

  std::remove(tmp_file.c_str());
}

void KmerIndex::BuildDeBruijnGraph(const ProgramOptions& opt, const std::string& tmp_file, std::ofstream& out) {

  CDBG_Build_opt c_opt;
  c_opt.k = k;
  c_opt.nb_threads = opt.threads;
  c_opt.build = true;
  c_opt.clipTips = false;
  c_opt.deleteIsolated = false;
  c_opt.verbose = true;
  c_opt.filename_ref_in.push_back(tmp_file);

  if (opt.g > 0) { // If minimizer length supplied, override the default
    c_opt.g = opt.g;
  } else { // Define minimizer length defaults
    int g = k-8;
    if (k <= 15) {
      g = k-2;
    } else if (k <= 17) {
      g = k-4;
    } else if (k <= 19) {
      g = k-6;
    }
    c_opt.g = g;
  }
  dbg = CompactedDBG<Node>(k, c_opt.g);
  dbg.build(c_opt);

  Offlist(opt.offlist);

  // Write graph to a temporary binary file and read it back into a static graph
  // TODO:
  // Only do this if graph exceeds certain size
  //dbg.to_static();
  //std::string tmp_bin = generate_tmp_file(tmp_file);
  //dbg.writeBinary(tmp_bin, opt.threads);
  //dbg = CompactedDBG<Node>(k, c_opt.g);
  //dbg.readBinary(tmp_bin, true, opt.threads);
  //std::remove(tmp_bin.c_str());

  // 1. write version
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. serialize dBG
  auto pos1 = out.tellp();
  // Write dummy size of graph
  size_t tmp_size = 1337;
  out.write((char *)&tmp_size, sizeof(tmp_size));
  bool res = dbg.writeBinary(out, opt.threads);

  if (res == 0) {
    std::cerr << "Error: could not write de Bruijn Graph to disk." << std::endl;
    exit(1);
  }

  auto pos2 = out.tellp();
  out.seekp(pos1);

  // Write real size of graph
  tmp_size = pos2 - pos1 - sizeof(tmp_size);
  out.write((char *)&tmp_size, sizeof(tmp_size));
  out.seekp(pos2);

  std::vector<Minimizer> minz;
  dbg.clearAndGetMinimizers(minz);
  boophf_t* mphf = new boophf_t(minz.size(), minz, opt.threads, 1.0, true, true);
  out.write((char *)&tmp_size, sizeof(tmp_size));
  mphf->save(out);
  pos1 = out.tellp();
  out.seekp(pos2);
  tmp_size = pos1 - pos2 - sizeof(tmp_size);
  out.write((char *)&tmp_size, sizeof(tmp_size));

  out.seekp(pos1);

  minz.clear();
  dbg.clear();

  std::ifstream infile;
  infile.open(opt.index, std::ios::in | std::ios::binary);
  std::istream in(0);
  in.rdbuf(infile.rdbuf());
  in.ignore(sizeof(INDEX_VERSION));
  in.ignore(sizeof(tmp_size));

  dbg = CompactedDBG<Node>(k, c_opt.g);
  dbg.readBinary(in, mphf);
  //infile.close();

  uint32_t running_id = 0;
  for (auto& um : dbg) {
    um.getData()->id = running_id++;
  }
}

void KmerIndex::Offlist(const std::string& path) {
  std::ifstream infile(path.c_str());
  if (!infile.good()) return;

  size_t count = 0;
  std::string line;
  std::stringstream buf;
  const_UnitigMap<Node> um;
  while (std::getline(infile, line)) {
    if (line[0] == '>' && buf.str() != "") {
      // Search for everything in buffer in graph and delete the corresponding
      // nodes in case of a hit.
      size_t pos = 0;
      std::string seq = buf.str();
      while (pos < seq.size() - k + 1) {

        um = dbg.findUnitig(seq.c_str(), pos, seq.size());

        if (um.isEmpty) {
          ++pos;
          continue;
        }

        pos += um.len;
        count += um.len;

        dbg.remove(um, true);
      }

      // Clear buffer.
      buf.str("");
    } else {
      buf << line;
    }
  }
  infile.close();
  std::cerr << "[build] Removed " << count << " blacklisted kmers from graph." << std::endl;
}

void KmerIndex::BuildEquivalenceClasses(const ProgramOptions& opt, const std::string& tmp_file) {

  std::cerr << "[build] creating equivalence classes ... " << std::endl;

  std::vector<std::vector<TRInfo> > trinfos(dbg.size());
  UnitigMap<Node> um;
  size_t EC_THRESHOLD = 250;
  size_t EC_SOFT_THRESHOLD = 800;
  size_t EC_MAX_N_ABOVE_THRESHOLD = 6000; // Thresholding ECs to size EC_THRESHOLD will only occur if we encounter >EC_MAX_N_ABOVE_THRESHOLD number of nodes that have size >EC_SOFT_THRESHOLD
  if (opt.max_ec_size > 0) { // If max EC size is supplied, override the default thresholds
    EC_THRESHOLD = opt.max_ec_size;
    EC_SOFT_THRESHOLD = EC_THRESHOLD;
    EC_MAX_N_ABOVE_THRESHOLD = 0;
  }
  uint32_t sense = 0x80000000, missense = 0;

  std::ifstream infile(tmp_file);
  std::string line;
  size_t j = 0;
  size_t n_above_threshold = 0;
  //for (size_t i = 0; i < seqs.size(); ++i) {
  while (std::getline(infile, line)) {
    if (line[0] == '>') continue;
    const auto& seq = line;
    if (seq.size() < k) { j++; continue; }

    int seqlen = seq.size() - k + 1; // number of k-mers
    size_t proc = 0;
    while (proc < seqlen) {
      um = dbg.findUnitig(seq.c_str(), proc, seq.size());

      if (um.isEmpty) {
        ++proc;
        continue;
      }

      proc += um.len;
      const Node* n = um.getData();
      if (trinfos[n->id].size() > EC_THRESHOLD) {
        if (trinfos[n->id].size() == EC_SOFT_THRESHOLD+1) {
          n_above_threshold++;
        }
        if (n_above_threshold > EC_MAX_N_ABOVE_THRESHOLD) {
          trinfos[n->id].clear();
          std::vector<TRInfo>().swap(trinfos[n->id]); // potentially free up memory
          TRInfo tr_discard;
          tr_discard.trid = std::numeric_limits<uint32_t>::max();
          trinfos[n->id].reserve(1);
          trinfos[n->id].push_back(tr_discard);
          continue;
        }
      } else if (trinfos[n->id].size() == 1 && trinfos[n->id][0].trid == std::numeric_limits<uint32_t>::max()) {
        continue;
      }
      TRInfo tr;

      tr.trid = j;
      tr.pos = proc | (!um.strand ? sense : missense);
      tr.start = um.dist;
      tr.stop  = um.dist + um.len;

      trinfos[n->id].reserve(trinfos[n->id].size()+1);
      trinfos[n->id].push_back(tr);
    }
    j++;
  }
  infile.close();

  // Threshold large ECs
  if (n_above_threshold > EC_MAX_N_ABOVE_THRESHOLD) {
    size_t n_removed = 0;
    for (auto& trinfo : trinfos) {
      if (trinfo.size() > EC_THRESHOLD) {
        trinfo.clear();
        std::vector<TRInfo>().swap(trinfo); // potentially free up memory
        ++n_removed;
      } else if (trinfo.size() == 1 && trinfo[0].trid == std::numeric_limits<uint32_t>::max()) {
        trinfo.clear();
        ++n_removed;
      }
    }
    std::cerr << "[build] discarded " << n_removed << " ECs larger than threshold." << std::endl;
  }

  PopulateMosaicECs(trinfos);

  std::cerr << "[build] target de Bruijn graph has k-mer length " << dbg.getK() << " and minimizer length "  << dbg.getG() << std::endl;
  std::cerr << "[build] target de Bruijn graph has " << dbg.size() << " contigs and contains "  << dbg.nbKmers() << " k-mers " << std::endl;
  //std::cerr << "[build] target de Bruijn graph contains " << ecmapinv.size() << " equivalence classes from " << seqs.size() << " sequences." << std::endl;
}

void KmerIndex::PopulateMosaicECs(std::vector<std::vector<TRInfo> >& trinfos) {

  for (const auto& um : dbg) {

    Node* n = um.getData();

    // Process empty ECs
    if (trinfos[n->id].size() == 0) {
      Roaring r;
      n->ec.insert(0, um.len, std::move(r));
      continue;
    }

    // Find the overlaps
    std::vector<int> brpoints;
    brpoints.reserve(2 * trinfos[n->id].size());
    for (const auto& x : trinfos[n->id]) {
      brpoints.push_back(x.start);
      brpoints.push_back(x.stop);
    }

    sort(brpoints.begin(), brpoints.end());
    assert(brpoints[0] == 0);
    assert(brpoints[brpoints.size()-1]==um.size-k+1);

    // Find unique break points
    if (!isUnique(brpoints)) {
      std::vector<int> u = unique(brpoints);
      swap(u,brpoints);
    }

    std::sort(trinfos[n->id].begin(), trinfos[n->id].end(),
              [](const TRInfo& lhs, const TRInfo& rhs) -> bool {
                return (lhs.trid < rhs.trid);
              });

    size_t j = 0;
    std::vector<uint32_t> pos;
    // Create a mosaic EC for the unitig, where each break point interval
    // corresponds to one set of transcripts and therefore an EC
    for (size_t i = 1; i < brpoints.size(); ++i) {

      Roaring u;

      for (const auto& tr : trinfos[n->id]) {
        // If a transcript encompasses the full breakpoint interval
        if (tr.start <= brpoints[i-1] && tr.stop >= brpoints[i]) {
          u.add(tr.trid);
          pos.reserve(pos.size()+1);
          pos.push_back(tr.pos);
        }
      }

      assert(!u.isEmpty());
      u.runOptimize();

      // Assign mosaic EC to the corresponding part of unitig
      n->ec.insert(brpoints[i-1], brpoints[i], std::move(u));
    }
    // Assign position and sense for all transcripts belonging to unitig
    n->pos = std::move(pos);
    std::vector<TRInfo>().swap(trinfos[n->id]); // potentially free up memory
  }
}

void KmerIndex::write(std::ofstream& out, int threads) {

  size_t tmp_size;

  if (!out.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index output file could not be opened!";
    exit(1);
  }


  // 3. serialize nodes
  tmp_size = dbg.size();
  out.write((char *)&tmp_size, sizeof(tmp_size));
  for (const auto& um : dbg) {
    // 3.1 write head kmer to associate unitig with node
    std::string kmer = um.getUnitigHead().toString();
    out.write(kmer.c_str(), strlen(kmer.c_str()));
    // 3.2 serialize node
    um.getData()->serialize(out);
  }

  // 4. write number of targets
  out.write((char *)&num_trans, sizeof(num_trans));

  // 5. write out target lengths
  for (int tlen : target_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  // 6. Write out target ids
  // XXX: num_trans should equal to target_names_.size(), so don't need
  // to write out again.
  assert(num_trans == target_names_.size());
  for (auto& tid : target_names_) {
    // 6.1 write out how many bytes
    tmp_size = strlen(tid.c_str());
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 6.2 write out the actual string
    out.write(tid.c_str(), tmp_size);
  }
}

void KmerIndex::write(const std::string& index_out, bool writeKmerTable, int threads) {

  std::ofstream out;
  out.open(index_out, std::ios::out | std::ios::binary);

  if (!out.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index output file could not be opened!";
    exit(1);
  }

  size_t tmp_size;

  // append is true if we have already written the INDEX_VERSION and the dBG
  // when constructing the dBG
  // 1. write version
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. serialize dBG
  if (writeKmerTable) {
    auto pos1 = out.tellp();
    // Write dummy size of graph
    tmp_size = 1337;
    out.write((char *)&tmp_size, sizeof(tmp_size));
    bool res = dbg.writeBinary(out, threads);

    if (res == 0) {
      std::cerr << "Error: could not write de Bruijn Graph to disk." << std::endl;
      exit(1);
    }

    auto pos2 = out.tellp();
    out.seekp(pos1);

    // Write real size of graph
    tmp_size = pos2 - pos1 - sizeof(tmp_size);
    out.write((char *)&tmp_size, sizeof(tmp_size));
    out.seekp(pos2);
  } else {
    tmp_size = 0;
    out.write((char *)&tmp_size, sizeof(tmp_size));
  }

  // 3. serialize nodes
  if (writeKmerTable) {
    tmp_size = dbg.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));
    for (const auto& um : dbg) {
      // 3.1 write head kmer to associate unitig with node
      std::string kmer = um.getUnitigHead().toString();
      out.write(kmer.c_str(), strlen(kmer.c_str()));
      // 3.2 serialize node
      um.getData()->serialize(out);
    }
  } else {
    tmp_size = 0;
    out.write((char *)&tmp_size, sizeof(tmp_size));
  }

  // 4. write number of targets
  out.write((char *)&num_trans, sizeof(num_trans));

  // 5. write out target lengths
  for (int tlen : target_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  // 6. Write out target ids
  // XXX: num_trans should equal to target_names_.size(), so don't need
  // to write out again.
  assert(num_trans == target_names_.size());
  for (auto& tid : target_names_) {
    // 6.1 write out how many bytes
    tmp_size = strlen(tid.c_str());
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 6.2 write out the actual string
    out.write(tid.c_str(), tmp_size);
  }

  out.flush();
  out.close();
}

void KmerIndex::load(ProgramOptions& opt, bool loadKmerTable) {

  if (opt.index.empty() && !loadKmerTable) {
    // Make an index from transcript and EC files
    loadTranscriptsFromFile(opt);
    loadECsFromFile(opt);
    return;
  }

  std::string& index_in = opt.index;
  std::ifstream in;//, in_minz;

  in.open(index_in, std::ios::in | std::ios::binary);
  //in_minz.open(index_in, std::ios::in | std::ios::binary);

  if (!in.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index input file could not be opened!";
    exit(1);
  }

  // 1. read version
  size_t header_version = 0;
  in.read((char *)&header_version, sizeof(header_version));
  //in_minz.ignore(sizeof(header_version));

  if (header_version != INDEX_VERSION) {
    std::cerr << "Error: incompatible indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
              << "Rerun with index to regenerate";
    exit(1);
  }

  // 2. deserialize dBG
  size_t tmp_size;
  in.read((char *)&tmp_size, sizeof(tmp_size));
  if (tmp_size > 0) {

    auto pos1 = in.tellg();
    in.ignore(tmp_size);
    boophf_t* mphf = new boophf_t();
    in.read((char *)&tmp_size, sizeof(tmp_size));
    mphf->load(in);
    auto pos2 = in.tellg();
    in.seekg(pos1);

    auto success = dbg.readBinary(in, mphf);

    in.seekg(pos2);

    //dbg.to_static();
    k = dbg.getK();
  }

  // 3. deserialize nodes
  Kmer kmer;
  UnitigMap<Node> um;
  size_t kmer_size = k * sizeof(char);
  char* buffer = new char[kmer_size];
  in.read((char *)&tmp_size, sizeof(tmp_size));
  for (size_t i = 0; i < tmp_size; ++i) {
    // 3.1 read head kmer
    memset(buffer, 0, kmer_size);
    in.read(buffer, kmer_size);
    kmer = Kmer(buffer);
    um = dbg.find(kmer);

    if (um.isEmpty) {
      std::cerr << "Error: Corrupted index; unitig not found: " << std::string(buffer) << std::endl;
      exit(1);
    }

    // 3.2 deserialize node
    um.getData()->deserialize(in);
  }
  delete[] buffer;
  buffer = nullptr;

  // 4. read number of targets
  in.read((char *)&num_trans, sizeof(num_trans));

  // 5. read out target lengths
  target_lens_.clear();
  target_lens_.reserve(num_trans);
  int tmp_int;
  for (size_t i = 0; i < num_trans; ++i) {
    in.read((char *)&tmp_int, sizeof(tmp_int));
    target_lens_.push_back(tmp_int);
  }

  // 6. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

  size_t bufsz = 1024;
  buffer = new char[bufsz];
  for (auto i = 0; i < num_trans; ++i) {
    // 6.1 read in the size
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size +1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }

    // clear the buffer
    memset(buffer,0,bufsz);
    // 6.2 read in the character string
    in.read(buffer, tmp_size);

    target_names_.push_back(std::string(buffer));
  }

  // delete the buffer
  delete[] buffer;
  buffer=nullptr;

  in.close();
  exit(0);
  if (!opt.ecFile.empty()) {
    loadECsFromFile(opt);
  }
}

void KmerIndex::loadECsFromFile(const ProgramOptions& opt) {
  ecmapinv.clear();
  int32_t i = 0;
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

      Roaring r;
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
        r.add(tmp_ecval_num);
      }
      ecmapinv.insert({std::move(r), i}); // move
      i++;
    }
  } else {
    std::cerr << "Error: could not open file " << opt.ecFile << std::endl;
    exit(1);
  }
  std::cerr << "[index] number of equivalence classes loaded from file: "
            << pretty_num(ecmapinv.size()) << std::endl;
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

int KmerIndex::mapPair(const char *s1, int l1, const char *s2, int l2) const {
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
      !(um1.getData()->ec[um1.dist] == um2.getData()->ec[um2.dist])) {
    return -1;
  }

  // Paired reads need to map to opposite strands
  if (!(um1.strand ^ um2.strand)) {
    //std::cerr << "Reads map to same strand " << s1 << "\t" << s2 << std::endl;
    return -1;
  }

  if (um1.getData()->get_mc_contig(um1.dist).second != um2.getData()->get_mc_contig(um2.dist).second) {
    return -1; // If the mc contigs for um1 and um2 are actually not the same (despite having the same color)
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
      auto p = n->get_mc_contig(um.dist);
      contig_start += p.first;
      contig_length = p.second - contig_start;

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
  uint32_t bitmask = 0x7FFFFFFF, rawpos;
  bool trsense = true;
  if (um.getData()->id == -1) {
    return {-1, true};
  }
  const Node* n = um.getData();
  auto ecs = n->ec.get_leading_vals(um.dist);
  size_t offset = 0;
  const Roaring& ec = ecs[ecs.size() - 1];

  for (size_t i = 0; i < ecs.size() - 1; ++i) {
    offset += ecs[i].cardinality();
  }

  uint32_t* trs = new uint32_t[ec.cardinality()];
  ec.toUint32Array(trs);
  for (size_t i = 0; i < ec.cardinality(); ++i) {
    if (trs[i] == tr) {
      rawpos = n->pos[offset + i];
      trpos = rawpos & bitmask;
      trsense = (rawpos == trpos);
      break;
    }
  }
  delete[] trs;
  trs = nullptr;

  if (trpos == -1) {
    return {-1,true};
  }

  // TODO:
  // Check whether um.dist does the same thing as KmerEntry.getPos();
  // If something doesn't work, it's most likely this!
  if (trsense) {
    if (csense) {
      return {(trpos - p - (um.size - um.dist - k)), csense}; // 1-based, case I
    } else {
      return {(trpos + p + k - (um.size - k + 1 - um.dist)), csense}; // 1-based, case III
    }
  } else {
    if (csense) {
      // UnitigMapBase.size is the length in base pairs
      return {trpos + (-um.dist - 1) + k + p, !csense};  // 1-based, case IV
    } else {
      return {trpos + (-um.dist) - p, !csense}; // 1-based, case II
    }
  }
}

// use:  res = intersect(ec,v)
// pre:  ec is an equivalence class, v is a vector of valid targets
//       v is sorted in increasing order
// post: res contains the intersection  of ec and v
Roaring KmerIndex::intersect(const Roaring& ec, const Roaring& v) const {
  Roaring res;
  if (ec.cardinality() == 0) {
    // If the EC has no transcripts, it has been filtered due to its size
    res = v;
  } else if (v.cardinality() == 0) {
    // If transcript vector is empty, it represents an EC that has been filtered
    // due to its size
    res = ec;
  } else {
    // Do an actual intersect
    res = ec & v;
  }
  return res;
}

void KmerIndex::loadTranscriptSequences() const {
  if (target_seqs_loaded) {
    return;
  }

  std::vector<std::vector<std::pair<std::string, u2t> > > trans_contigs(num_trans);
  for (auto& um : dbg) {
    const Node* data = um.getData();
    std::string um_seq = um.mappedSequenceToString();
    // XXX:
    // We also reverse complement in the following for-loop. Is this correct?
    if (!um.strand) {
      um_seq = revcomp(um_seq);
    }

    const Node* n = um.getData();
    auto ecs = n->ec.get_leading_vals(um.dist);
    size_t offset = 0;
    const Roaring& ec = ecs[ecs.size() - 1];
    for (size_t i = 0; i < ecs.size() - 1; ++i) {
      offset += ecs[i].cardinality();
    }

    uint32_t* trs = new uint32_t[ec.cardinality()];
    ec.toUint32Array(trs);
    for (size_t i = 0; i < ec.cardinality(); ++i) {
      bool sense = (n->pos[offset + i] & 0x7FFFFFFF) != n->pos[offset + i];
      u2t tr(trs[i], n->pos[offset + i]);
      // TODO:
      // Verify that this method of choosing to reverse complement is legit
      if (um.strand ^ sense) {
        trans_contigs[trs[i]].emplace_back(revcomp(um_seq), tr);
      } else {
        trans_contigs[trs[i]].emplace_back(um_seq, tr);
      }
    }
  }

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
  dbg.clear();

  {
    EcMapInv empty;
    std::swap(ecmapinv, empty);
  }

  target_lens_.resize(0);
  target_names_.resize(0);
  target_seqs_.resize(0);
}

void KmerIndex::writePseudoBamHeader(std::ostream &o) const {
  // write out header
  o << "@HD\tVN:1.0\n";
  for (int i = 0; i < num_trans; i++) {
    o << "@SQ\tSN:" << target_names_[i] << "\tLN:" << target_lens_[i] << "\n";
  }
  o << "@PG\tID:kallisto\tPN:kallisto\tVN:"<< KALLISTO_VERSION << "\n";
  o.flush();
}
