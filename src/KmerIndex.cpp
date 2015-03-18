#include "KmerIndex.h"
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>


using namespace seqan;

bool KmerIndex::loadSuffixArray(const ProgramOptions& opt ) {
  // open the fasta sequences
  
  index = TIndex(seqs);
  // try to open the file

  if (!open(getFibre(index, FibreSA()), (opt.index + ".sa").c_str(), OPEN_RDONLY | OPEN_QUIET)) {
    return false;
  }
  //setHaystack(finder, index); // finder searches index
  
  return true;
}

void KmerIndex::clearSuffixArray() {
  clear(index);
  //clear(finder);
}

void KmerIndex::BuildTranscripts(const ProgramOptions& opt) {
  const std::string& fasta = opt.transfasta;
  std::cerr << "Loading fasta file " << fasta
            << std::endl;
  std::cerr << "k: " << k << std::endl;
  

  SeqFileIn seqFileIn(fasta.c_str());

  // read fasta file
  seqan::StringSet<seqan::CharString> ids;
  readRecords(ids, seqs, seqFileIn);

  int transid = 0;
  

  if (!loadSuffixArray(opt)) {
    std::cerr << "Constructing index ... "; std::cerr.flush();  

    index = TIndex(seqs);
    indexRequire(index, EsaSA());
    // write fasta to disk
    SeqFileOut fastaIndex;
    if (!open(fastaIndex, (opt.index+".fa").c_str())) {
      std::cerr << "Error: could not open file " << opt.index << ".fa for writing" << std::endl;
      exit(1);
    }
    try {
      writeRecords(fastaIndex, ids, seqs);
    } catch (IOError const & e) {
      std::cerr << "Error: writing to file " << opt.index << ". " << e.what() << std::endl; 
    }
    close(fastaIndex);
    
    // write index to disk
    save(getFibre(index, FibreSA()),(opt.index + ".sa").c_str());

    std::cerr << " done " << std::endl;
  } else {
    std::cerr << "Found index" << std::endl;
  }

  TFinder finder(index);
  find(finder, "A");
  clear(finder);

  
  assert(length(seqs) == length(ids));
  num_trans = length(seqs);

  
  
  // for each transcript, create it's own equivalence class
  for (int i = 0; i < length(seqs); i++ ) {
    std::string name(toCString(ids[i]));
    size_t p = name.find(' ');
    if (p != std::string::npos) {
      name = name.substr(0,p);
    }
    target_names_.push_back(name);

    CharString seq = value(seqs,i);
    trans_lens_.push_back(length(seq));

    std::vector<int> single(1,i);
    ecmap.insert({i,single});
    ecmapinv.insert({single,i});
  }

  // traverse the suffix tree
  Iterator<TIndex, TopDown<ParentLinks<>>>::Type it(index);
  do {
    
    if (repLength(it) >= k) {
      //std::cout << representative(it) << std::endl;
      // process the k-mers
      CharString seq = representative(it);
      Kmer km(toCString(seq));
      Kmer rep = km.rep();
      if (kmap.find(rep) == kmap.end()) {
        // if we have never seen this k-mer before
        
        std::vector<int> ecv;

        for (auto x : getOccurrences(it)) {
          ecv.push_back(x.i1); // record transcript
        }

        // search for second part
        Kmer twin = km.twin(); // other k-mer
        clear(finder);
        while (find(finder, twin.toString().c_str())) {
          ecv.push_back(getSeqNo(position(finder)));
        }
        
        // common case
        if (ecv.size() == 1) {
          int ec = ecv[0];
          kmap.insert({rep, ec});
        } else {
          sort(ecv.begin(), ecv.end());
          std::vector<int> u;
          u.push_back(ecv[0]);
          for (int i = 1; i < ecv.size(); i++) {
            if (ecv[i-1] != ecv[i]) {
              u.push_back(ecv[i]);
            }
          }
          
          auto search = ecmapinv.find(u);
          if (search != ecmapinv.end()) {
            kmap.insert({rep, search->second});
          } else {
            int eqs_id = ecmapinv.size();
            ecmapinv.insert({u, eqs_id });
            ecmap.insert({eqs_id, u});
            kmap.insert({rep, eqs_id});
          }
        }
      }
    }
    // next step
    if (!goDown(it) || repLength(it) > k) {
      // if we can't go down or the sequence is too long
      do {
        if (!goRight(it)) {
          while (goUp(it) && !goRight(it)) {
            // go up the tree until you can go to the right
          }
        }
      } while (repLength(it) > k);
    }
  } while (!isRoot(it));
  
  /*
  // for each transcript
  for (int i = 0; i < length(seqs); i++) {
    // if it it long enough
    CharString seq = value(seqs,i);
    if (length(seq) >= k) {
      KmerIterator kit(toCString(seq)), kit_end;

      // for each k-mer add to map
      for(; kit != kit_end; ++kit) {
        Kmer km = kit->first;
        Kmer rep = km.rep();
        if (kmap.find(km) != kmap.end()) {
          // we've seen this k-mer before
          continue;
        }
        
        CharString kms(km.toString().c_str()); // hacky, but works

        std::vector<int> ecv;
        
        // find all instances of the k-mer
        clear(finder);
        while (find(finder, kms)) {
          ecv.push_back(getSeqNo(position(finder)));
        }
        // find all instances of twin
        reverseComplement(kms);
        clear(finder);
        while (find(finder,kms)) {
          ecv.push_back(getSeqNo(position(finder)));
        }

        // most common case
        if (ecv.size() == 1) {
          int ec = ecv[0];
          kmap.insert({rep, ec});
        } else {
          sort(ecv.begin(), ecv.end());
          std::vector<int> u;
          u.push_back(ecv[0]);
          for (int i = 1; i < ecv.size(); i++) {
            if (ecv[i-1] != ecv[i]) {
              u.push_back(ecv[i]);
            }
          }
          auto search = ecmapinv.find(u);
          if (search != ecmapinv.end()) {
            kmap.insert({rep, search->second});
          } else {
            int eqs_id = ecmapinv.size();
            ecmapinv.insert({u, eqs_id });
            ecmap.insert({eqs_id, u});
            kmap.insert({rep, eqs_id});
          }
        }
      }
    }
  }
    */

  // remove polyA close k-mers
  CharString polyA;
  resize(polyA, k, 'A');
  std::cerr << "searching for neighbors of " << polyA << std::endl;

  {
    auto search = kmap.find(Kmer(toCString(polyA)).rep());
    if (search != kmap.end()){
      kmap.erase(search);
    }
  }
  
  for (int i = 0; i < k; i++) {
    for (int a = 1; a < 4; a++) {
      CharString x(polyA);
      assignValue(x, i, Dna(a));
      {
        auto search = kmap.find(Kmer(toCString(x)).rep());
        if (search != kmap.end()){
          kmap.erase(search);
        }
      }
      /*
      for (int j = i+1; j < k; j++) {
        CharString y(x);
        for (int b = 1; b < 4; b++) {
          assignValue(y, j, Dna(b));
          {
            auto search = kmap.find(Kmer(toCString(y)).rep());
            if (search != kmap.end()){
              kmap.erase(search);
            }
          }
        }
      }
      */
      
    }
  }

  std::cerr << "Found " << num_trans << " transcripts"
            << std::endl
            << "Size of k-mer map " << kmap.size()<< std::endl;


 
  

  int eqs_id = num_trans;


  std::cerr << "Created " << ecmap.size() << " equivalence classes from " << num_trans << " transcripts" << std::endl;

  std::cerr << "K-mer map has " << kmap.size() << " k-mers and " << std::endl;
}


void KmerIndex::write(const std::string& index_out, bool writeKmerTable) {
  std::ofstream out;
  out.open(index_out, std::ios::out | std::ios::binary);

  if (!out.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index output file could not be opened!";
    exit(1);
  }

  // 1. write index
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. write k
  out.write((char *)&k, sizeof(k));

  // 3. write number of transcripts
  out.write((char *)&num_trans, sizeof(num_trans));

  // 4. write out transcript lengths
  for (int tlen : trans_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  size_t kmap_size = kmap.size();

  if (writeKmerTable) {
    // 5. write number of k-mers in map
    out.write((char *)&kmap_size, sizeof(kmap_size));

    // 6. write kmer->ec values
    for (auto& kv : kmap) {
      out.write((char *)&kv.first, sizeof(kv.first));
      out.write((char *)&kv.second, sizeof(kv.second));
    }
  } else {
    // 5. write fake k-mer size
    kmap_size = 0;
    out.write((char *)&kmap_size, sizeof(kmap_size));

    // 6. write none of the kmer->ec values
  }
  // 7. write number of equivalence classes
  size_t tmp_size;
  tmp_size = ecmap.size();
  out.write((char *)&tmp_size, sizeof(tmp_size));

  // 8. write out each equiv class
  for (auto& kv : ecmap) {
    out.write((char *)&kv.first, sizeof(kv.first));

    // 8.1 write out the size of equiv class
    tmp_size = kv.second.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));
    // 8.2 write each member
    for (auto& val: kv.second) {
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
    tmp_size = tid.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 9.2 write out the actual string
    out.write(tid.c_str(), tid.size());
  }

  out.flush();
  out.close();
}


void KmerIndex::load(ProgramOptions& opt, bool loadKmerTable) {


  
  std::string& index_in = opt.index;
  std::ifstream in;

  // load sequences
  SeqFileIn seqFileIn;

  if (!open(seqFileIn, (index_in + ".fa").c_str())){
    std::cerr << "Error: index fasta could not be opened " << std::endl;
    exit(1);
  }
  seqan::StringSet<seqan::CharString> ids; // not used
  readRecords(ids, seqs, seqFileIn);
  
  
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
    std::cerr << "Error: Incompatiple indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
              << "Rerun with index to regenerate!";
    exit(1);
  }

  // 2. read k
  in.read((char *)&k, sizeof(k));
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

  // 3. read number of transcripts
  in.read((char *)&num_trans, sizeof(num_trans));

  // 4. read number of transcripts
  trans_lens_.clear();
  trans_lens_.reserve(num_trans);

  for (int i = 0; i < num_trans; i++) {
    int tlen;
    in.read((char *)&tlen, sizeof(tlen));
    trans_lens_.push_back(tlen);
  }

  // 5. read number of k-mers
  size_t kmap_size;
  in.read((char *)&kmap_size, sizeof(kmap_size));

  std::cerr << "[index] k: " << k << std::endl;
  std::cerr << "[index] num_trans read: " << num_trans << std::endl;
  std::cerr << "[index] kmap size: " << kmap_size << std::endl;

  kmap.clear();
  if (loadKmerTable) {
    kmap.reserve(kmap_size);
  }

  // 6. read kmer->ec values
  Kmer tmp_kmer;
  int tmp_val;
  for (size_t i = 0; i < kmap_size; ++i) {
    in.read((char *)&tmp_kmer, sizeof(tmp_kmer));
    in.read((char *)&tmp_val, sizeof(tmp_val));

    if (loadKmerTable) {
      kmap.insert({tmp_kmer, tmp_val});
    }
  }

  // 7. read number of equivalence classes
  size_t ecmap_size;
  in.read((char *)&ecmap_size, sizeof(ecmap_size));

  std::cerr << "[index] ecmap size: " << ecmap_size << std::endl;

  int tmp_id;
  size_t vec_size;
  // 8. read each equiv class
  for (size_t i = 0; i < ecmap_size; ++i) {
    in.read((char *)&tmp_id, sizeof(tmp_id));

    // 8.1 read size of equiv class
    in.read((char *)&vec_size, sizeof(vec_size));

    // 8.2 read each member
    std::vector<int> tmp_vec;
    tmp_vec.reserve(vec_size);
    for (size_t j = 0; j < vec_size; ++j ) {
      in.read((char *)&tmp_val, sizeof(tmp_val));
      tmp_vec.push_back(tmp_val);
    }
    ecmap.insert({tmp_id, tmp_vec});
    ecmapinv.insert({tmp_vec, tmp_id});
  }

  // 9. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

  size_t tmp_size;
  char buffer[1024]; // if your target_name is longer than this, screw you.
  for (auto i = 0; i < num_trans; ++i) {
    // 9.1 read in the size
    in.read((char *)&tmp_size, sizeof(tmp_size));

    // 9.2 read in the character string
    in.read(buffer, tmp_size);

    std::string tmp_targ_id( buffer );
    target_names_.push_back(std::string( buffer ));

    // clear the buffer for next string
    memset(buffer,0,strlen(buffer));
  }

  in.close();
}

int KmerIndex::mapPair(const char* s1, int l1, const char* s2, int l2, int ec, TFinder& finder) const {
  bool d1 = true;
  bool d2 = true;
  int p1 = -1;
  int p2 = -1;

  
  KmerIterator kit1(s1), kit_end;

  bool found1 = false;
  while (kit1 != kit_end) {
    if (kmap.find(kit1->first.rep()) != kmap.end()) {
      found1 = true;
      break;
    }
  }

  if (!found1) {
    return -1;
  }

  clear(finder);
  //CharString c1(s1+kit1->second);
  // try to locate it in the suffix array
  //Prefix<CharString>::Type
  //CharString pre = prefix(c1,k);
  //setNeedle(finder, kit1->first.toString().c_str());
  while(find(finder, kit1->first.toString().c_str())) {
    if (getSeqNo(position(finder)) == ec) {
      p1 = getSeqOffset(position(finder)) - kit1->second;
      break;
    }
  }
  if (p1 < 0) {
    //reverseComplement(pre);
    clear(finder);
    while(find(finder, kit1->first.twin().toString().c_str())) {
      if (getSeqNo(position(finder)) == ec) {
        p1 = getSeqOffset(position(finder)) + k + kit1->second;
        d1 = false;
        break;
      }
    }
  }

  if (p1 == -1) {
    return -1;
  }


  KmerIterator kit2(s2);
  bool found2 = false;
  while (kit2 != kit_end) {
    if (kmap.find(kit2->first.rep()) != kmap.end()) {
      found2 = true;
      break;
    }
  }

  if (!found2) {
    return -1;
  }


  //CharString c2(s2+kit2->second);
  // try to locate it in the suffix array
  //Prefix<CharString>::Type pre = prefix(c2,k);
  clear(finder);
  while(find(finder, kit2->first.toString().c_str())) {
    if (getSeqNo(position(finder)) == ec) {
      p2 = getSeqOffset(position(finder)) - kit2->second;
      break;
    }
  }
  if (p2 < 0) {
    //reverseComplement(pre);
    clear(finder);
    while(find(finder, kit2->first.twin().toString().c_str())) {
      if (getSeqNo(position(finder)) == ec) {
        p2 = getSeqOffset(position(finder)) + k + kit2->second;
        d2 = false;
        break;
      }
    }
  }

  if (p2 == -1) {
    return -1;
  }

  if ((d1 && d2) || (!d1 && !d2)) {
    //std::cerr << "Reads map to same strand " << s1 << "\t" << s2 << std::endl;
    return -1;
  }

  if (p1>p2) {
    return p1-p2;
  } else {
    return p2-p1;
  }
  
}
