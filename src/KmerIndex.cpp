#include "KmerIndex.h"


#include <zlib.h>
#include "kseq.h"

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
  int k = opt.k;
  std::cerr << "[build] Loading fasta file " << opt.transfasta
            << std::endl;
  std::cerr << "[build] k: " << k << std::endl;


  std::vector<std::string> seqs;

  // read fasta file  
  gzFile fp = 0;
  kseq_t *seq;
  int l = 0;
  fp = gzopen(opt.transfasta.c_str(), "r");
  seq = kseq_init(fp);
  while (true) {
    l = kseq_read(seq);
    if (l <= 0) {
      break;
    }
    seqs.push_back(seq->seq.s);
    trans_lens_.push_back(seq->seq.l);
    std::string name(seq->name.s);
    size_t p = name.find(' ');
    if (p != std::string::npos) {
      name = name.substr(0,p);
    }
    target_names_.push_back(name);

  }
  gzclose(fp);
  fp=0;
  
  num_trans = seqs.size();
  
  // for each transcript, create it's own equivalence class
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
  

  std::cerr << "[build] Counting k-mers ... "; std::cerr.flush();
  // gather all k-mers
  for (int i = 0; i < seqs.size(); i++) {
    const char *s = seqs[i].c_str();
    KmerIterator kit(s),kit_end;
    for (; kit != kit_end; ++kit) {
      kmap.insert({kit->first.rep(), KmerEntry()}); // don't care about repeats
    }
  }
  std::cerr << "done." << std::endl;
  
  std::cerr << "[build] Building de Bruijn Graph ... "; std::cerr.flush();
  // find out how much we can skip ahead for each k-mer.
  for (auto& kv : kmap) {
    if (kv.second.contig == -1) {
      // ok we haven't processed the k-mer yet
      std::vector<Kmer> flist, blist;

      // iterate in forward direction
      Kmer km = kv.first;
      Kmer end = km;
      Kmer last = end;
      Kmer twin = km.twin();
      bool selfLoop = false;
      flist.push_back(km);

      while (fwStep(end,end)) {
        if (end == km) {
          // selfloop
          selfLoop = true;
          break;
        } else if (end == twin) {
          selfLoop = (flist.size() > 1); // hairpins are not loops
          // mobius loop
          break;
        } else if (end == last.twin()) {
          // hairpin
          break;
        }
        flist.push_back(end);
        last = end;
      }

      Kmer front = twin;
      Kmer first = front;

      if (!selfLoop) {
        while (fwStep(front,front)) {
          if (front == twin) {
            // selfloop
            selfLoop = true;
            break;
          } else if (front == km) {
            // mobius loop
            selfLoop = true;
            break;
          } else if (front == first.twin()) {
            // hairpin
            break;
          }
          blist.push_back(front);
          first = front;
        }
      }

      std::vector<Kmer> klist;
      for (auto it = blist.rbegin(); it != blist.rend(); ++it) {
        klist.push_back(it->twin());
      }
      for (auto x : flist) {
        klist.push_back(x);
      }


      Contig contig;
      contig.id = dbGraph.contigs.size();
      contig.length = klist.size();
      contig.seq = klist[0].toString();
      contig.seq.reserve(contig.length + k-1);


      for (int i = 0; i < klist.size(); i++) {
        Kmer x = klist[i];
        Kmer xr = x.rep();
        bool forward = (x==xr);
        auto it = kmap.find(xr);
        it->second = KmerEntry(contig.id, i, forward);
        if (i > 0) {
          contig.seq.push_back(x.toString()[k-1]);
        }
      }
      
      dbGraph.contigs.push_back(contig);
      dbGraph.ecs.push_back(-1);
    }
  }
  std::cerr << " finished building de Bruijn Graph found for " << kmap.size() << " kmers and with " << dbGraph.contigs.size() << " contigs" << std::endl;
  std::cerr << "[build] Creating equivalence classes ... "; std::cerr.flush();
}

void KmerIndex::BuildEquivalenceClasses(const ProgramOptions& opt, const std::vector<std::string>& seqs) {

  std::vector<std::vector<TRInfo>> trinfos(dbGraph.contigs.size());
  //std::cout << "Mapping transcript " << std::endl;
  for (int i = 0; i < seqs.size(); i++) {
    int seqlen = seqs[i].size() - k + 1; // number of k-mers
    const char *s = seqs[i].c_str();
    //std::cout << "sequence number " << i << std::endl;
    KmerIterator kit(s), kit_end;
    for (; kit != kit_end; ++kit) {
      Kmer x = kit->first;
      Kmer xr = x.rep();
      auto search = kmap.find(xr);
      bool forward = (x==xr);
      KmerEntry val = search->second;
      std::vector<TRInfo>& trinfo = trinfos[val.contig];
      Contig& contig = dbGraph.contigs[val.contig];

      

      TRInfo tr;
      tr.trid = i;
      int jump = kit->second;
      if (forward == val.isFw()) {
        tr.start = val.getPos();
        if (contig.length - tr.start > seqlen - kit->second) {
          // transcript stops
          tr.stop = tr.start + seqlen - kit->second;
          jump = seqlen;
        } else {
          tr.stop = contig.length;
          jump = kit->second + (tr.stop - tr.start)-1;
        }
      } else {
        tr.stop = val.getPos()+1;
        int stpos = tr.stop - (seqlen - kit->second);
        if (stpos > 0) {
          tr.start = stpos;
          jump = seqlen;
        } else {
          tr.start = 0;
          jump = kit->second + (tr.stop - tr.start) - 1;
        }
      }

      // debugging -->
      //std::cout << "covering seq" << std::endl << seqs[i].substr(kit->second, jump-kit->second +k) << std::endl;
      //std::cout << "id = " << tr.trid << ", (" << tr.start << ", " << tr.stop << ")" << std::endl;
      //std::cout << "contig seq" << std::endl;
      if (forward == val.isFw()) {
        //std::cout << contig.seq << std::endl;
        //assert(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start) == seqs[i].substr(kit->second, jump-kit->second +k) );
      } else {
        //std::cout << revcomp(contig.seq) << std::endl;
        //assert(revcomp(contig.seq) == seqs[i].substr(kit->second, jump-kit->second +k));
        //assert(revcomp(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start)) == seqs[i].substr(kit->second, jump-kit->second +k));
      }
      if (jump == seqlen) {
        //std::cout << std::string(k-1+(tr.stop-tr.start)-1,'-') << "^" << std::endl;
      }

      // <-- debugging
      
      trinfo.push_back(tr);
      kit.jumpTo(jump);
    }
  }

  
  FixSplitContigs(opt, trinfos);
  int perftr = 0;
  for (int i = 0; i < trinfos.size(); i++) {
    bool all = true;

    int contigLen = dbGraph.contigs[i].length;
    //std::cout << "contig " << i << ", length = " << contigLen << ", seq = " << dbGraph.contigs[i].seq << std::endl << "tr = ";
    for (auto x : trinfos[i]) {
      if (x.start!=0 || x.stop !=contigLen) {
        all = false;
      }
      //std::cout << "[" << x.trid << ",(" << x.start << ", " << x.stop << ")], " ;
    }
    //std::cout << std::endl;
    if (all) {
      perftr++;
    } 
  }
  std::cerr << "For " << dbGraph.contigs.size() << ", " << (dbGraph.contigs.size() - perftr) << " of them need to be split" << std::endl;
  // need to create the equivalence classes

  assert(dbGraph.contigs.size() == trinfos.size());
  // for each contig
  for (int i = 0; i < trinfos.size(); i++) {
    std::vector<int> u;
    for (auto x : trinfos[i]) {
      u.push_back(x.trid);
    }
    sort(u.begin(), u.end());
    if (!isUnique(u)){
      std::vector<int> v = unique(u);
      swap(u,v);
    }

    assert(!u.empty());

    auto search = ecmapinv.find(u);
    if (search != ecmapinv.end()) {
      // insert contig -> ec info
      dbGraph.ecs[i] = search->second;
    } else {
      int eqs_id = ecmapinv.size();
      ecmapinv.insert({u,eqs_id});
      ecmap.push_back(u);
      dbGraph.ecs[i] = eqs_id;
    }
  }
}

void KmerIndex::FixSplitContigs(const ProgramOptions& opt, std::vector<std::vector<TRInfo>>& trinfos) {

  int perftr = 0;
  int orig_size = trinfos.size();

  for (int i = 0; i < orig_size; i++) {
    bool all = true;

    int contigLen = dbGraph.contigs[i].length;
    //std::cout << "contig " << i << ", length = " << contigLen << ", seq = " << dbGraph.contigs[i].seq << std::endl << "tr = ";
    for (auto x : trinfos[i]) {
      if (x.start!=0 || x.stop !=contigLen) {
        all = false;
      }
      //std::cout << "[" << x.trid << ",(" << x.start << ", " << x.stop << ")], " ;
      assert(x.start < x.stop);
    }
    //std::cout << std::endl;

    if (all) {
      perftr++;
    } else {
      // break up equivalence classes
      // sort by start/stop
      std::vector<int> brpoints;
      for (auto& x : trinfos[i]) {
        brpoints.push_back(x.start);
        brpoints.push_back(x.stop);
      }
      sort(brpoints.begin(), brpoints.end());
      assert(brpoints[0] == 0);
      assert(brpoints[brpoints.size()-1]==contigLen);

      // find unique points
      if (!isUnique(brpoints)) {
        std::vector<int> u = unique(brpoints);
        swap(u,brpoints);
      }

      assert(!brpoints.empty());
      
      // copy sequence
      std::string seq = dbGraph.contigs[i].seq;
      // copy old trinfo
      std::vector<TRInfo> oldtrinfo = trinfos[i];
      
      for (int j = 1; j < brpoints.size(); j++) {
        assert(brpoints[j-1] < brpoints[j]);
        Contig newc;
        newc.seq = seq.substr(brpoints[j-1], brpoints[j]-brpoints[j-1]+k-1);
        newc.length = brpoints[j]-brpoints[j-1];

        if (j>1) {
          newc.id = dbGraph.contigs.size();
          dbGraph.contigs.push_back(newc);
          dbGraph.ecs.push_back(-1);
        } else {
          newc.id = i;
          dbGraph.contigs[i] = newc;
        }

        // repair k-mer mapping
        KmerIterator kit(newc.seq.c_str()), kit_end;
        for (; kit != kit_end; ++kit) {
          Kmer x = kit->first;
          Kmer xr = x.rep();
          auto search = kmap.find(xr);
          assert(search != kmap.end());
          bool forward = (x==xr);
          search->second = KmerEntry(newc.id, kit->second, forward);
        }

        // repair tr-info
        std::vector<TRInfo> newtrinfo;
        for (auto x : oldtrinfo) {
          if (!(x.stop <= brpoints[j-1] || x.start >= brpoints[j])) {
            TRInfo trinfo;
            trinfo.trid = x.trid;
            trinfo.start = 0;
            trinfo.stop = newc.length;
            newtrinfo.push_back(trinfo);
          }
        }
        if (j>1) {
          trinfos.push_back(newtrinfo);
        } else {
          trinfos[i] = newtrinfo;
        }
      }
    }
  }


  std::cerr << "For " << dbGraph.contigs.size() << ", " << (dbGraph.contigs.size() - perftr) << " of them need to be split" << std::endl;

  
}




/*
void dummy_funk() {
  const std::string& fasta = opt.transfasta;
  std::cerr << "[build] Loading fasta file " << fasta
            << std::endl;
  std::cerr << "[build] k: " << k << std::endl;


  SeqFileIn seqFileIn(fasta.c_str());

  // read fasta file
  seqan::StringSet<seqan::CharString> ids;
  readRecords(ids, seqs, seqFileIn);

  int transid = 0;

  // probably remove this -->
  
  if (!loadSuffixArray(opt)) {
    std::cerr << "[build] Constructing suffix array ... "; std::cerr.flush();

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
    } catch (IOError const& e) {
      std::cerr << "Error: writing to file " << opt.index << ". " << e.what() << std::endl;
    }
    close(fastaIndex);

    // write index to disk
    save(getFibre(index, FibreSA()),(opt.index + ".sa").c_str());

    std::cerr << " finished " << std::endl;
  } else {
    std::cerr << "[build] Found suffix array" << std::endl;
  }

  TFinder finder(index);
  find(finder, "A");
  clear(finder);
  * /
  // <-- 
  

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
    //ecmap.insert({i,single});
    ecmap.push_back(single);
    ecmapinv.insert({single,i});
  }


  // remove this <-- 
  // traverse the suffix tree
  / *
  int nk = 0;
  Iterator<TIndex, TopDown<ParentLinks<>>>::Type it(index);
  do {

    if (repLength(it) >= k) {
      nk++;

      // if (nk % 1000 == 0) {
      //   std::cerr << "."; std::cerr.flush();
      // }

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
          kmap.insert({rep, KmerEntry(ec)});
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
            kmap.insert({rep, KmerEntry(search->second)});
          } else {
            int eqs_id = ecmapinv.size();
            ecmapinv.insert({u, eqs_id });
            ecmap.push_back(u);
            //ecmap.insert({eqs_id, u});
            kmap.insert({rep, KmerEntry(eqs_id)});
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
  * /
  // <-- 

  std::cerr << " ... finished creating equivalence classes" << std::endl;


  / *
  // remove polyA close k-mers
  CharString polyA;
  resize(polyA, k, 'A');
  std::cerr << "[build] Removing all k-mers within hamming distance 1 of " << polyA << std::endl;

  {
    auto search = kmap.find(Kmer(toCString(polyA)).rep());
    if (search != kmap.end()) {
      kmap.erase(search);
    }
  }

  for (int i = 0; i < k; i++) {
    for (int a = 1; a < 4; a++) {
      CharString x(polyA);
      assignValue(x, i, Dna(a));
      {
        auto search = kmap.find(Kmer(toCString(x)).rep());
        if (search != kmap.end()) {
          kmap.erase(search);
        }
      }
      
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
      //
    }
  }
  * /

  

  std::cerr << "[build] Found " << num_trans << " transcripts"
            << std::endl;

  int eqs_id = num_trans;

  std::cerr << "[build] Created " << ecmap.size() << " equivalence classes from " << num_trans << " transcripts" << std::endl;

  std::cerr << "[build] K-mer map has " << kmap.size() << " k-mers" << std::endl;
}

*/

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
  //  for (auto& kv : ecmap) {
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
  assert(dbGraph.contigs.size() == dbGraph.ecs.size());
  tmp_size = dbGraph.contigs.size();
  out.write((char*)&tmp_size, sizeof(tmp_size));
  for (auto& contig : dbGraph.contigs) {
    out.write((char*)&contig.id, sizeof(contig.id));
    out.write((char*)&contig.length, sizeof(contig.length));
    tmp_size = strlen(contig.seq.c_str());
    out.write((char*)&tmp_size, sizeof(tmp_size));
    out.write(contig.seq.c_str(), tmp_size);
  }

  // 11. write out ecs info
  for (auto ec : dbGraph.ecs) {
    out.write((char*)&ec, sizeof(ec));
  }
  

  out.flush();
  out.close();
}

bool KmerIndex::fwStep(Kmer km, Kmer& end) const {
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

}

void KmerIndex::load(ProgramOptions& opt, bool loadKmerTable) {

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

  // 4. read length of transcripts
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
  KmerEntry tmp_val;
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

  // 9. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

  size_t tmp_size;
  size_t bufsz = 1024;
  char *buffer = new char[bufsz];
  for (auto i = 0; i < num_trans; ++i) {
    // 9.1 read in the size
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size > bufsz+1) {
      delete[] buffer;
      bufsz *= 2;
      buffer = new char[bufsz];
    }
    
    // clear the buffer 
    memset(buffer,0,bufsz);
    // 9.2 read in the character string
    in.read(buffer, tmp_size);

    /* std::string tmp_targ_id( buffer ); */
    target_names_.push_back(std::string( buffer ));
  }


  // 10. read contigs
  size_t contig_size;
  in.read((char *)&contig_size, sizeof(contig_size));
  dbGraph.contigs.clear();
  dbGraph.contigs.reserve(contig_size);
  for (auto i = 0; i < contig_size; i++) {
    Contig c;
    in.read((char *)&c.id, sizeof(c.id));
    in.read((char *)&c.length, sizeof(c.length));
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size > bufsz + 1) {
      delete[] buffer;
      bufsz *= 2;
      buffer = new char[bufsz];
    }

    memset(buffer,0,bufsz);
    in.read(buffer, tmp_size);
    c.seq = std::string(buffer); // copy
    dbGraph.contigs.push_back(c);
  }

  // 11. read ecs info
  dbGraph.ecs.clear();
  dbGraph.ecs.reserve(contig_size);
  int tmp_ec;
  for (auto i = 0; i < contig_size; i++) {
    in.read((char *)&tmp_ec, sizeof(tmp_ec));
    dbGraph.ecs.push_back(tmp_ec);
  }

  // delete the buffer
  delete[] buffer;
  buffer=nullptr;
  
  in.close();
}

int KmerIndex::mapPair(const char *s1, int l1, const char *s2, int l2, int ec) const {
  bool d1 = true;
  bool d2 = true;
  int p1 = -1;
  int p2 = -1;
  int c1 = -1;
  int c2 = -1;


  KmerIterator kit1(s1), kit_end;

  bool found1 = false;
  for (; kit1 != kit_end; ++kit1) {
    Kmer x = kit1->first;
    Kmer xr = x.rep();
    auto search = kmap.find(xr);
    bool forward = (x==xr);

    if (search != kmap.end()) {
      found1 = true;
      KmerEntry val = search->second;
      c1 = val.contig;
      if (forward == val.isFw()) {
        p1 = val.getPos() - kit1->second;
        d1 = true;
      } else {
        p1 = val.getPos() + k + kit1->second;
        d1 = false;
      }
      break;
    }
  }

  if (!found1) {
    return -1;
  }

  

  KmerIterator kit2(s2);
  bool found2 = false;

  for (; kit2 != kit_end; ++kit2) {
    Kmer x = kit2->first;
    Kmer xr = x.rep();
    auto search = kmap.find(xr);
    bool forward = (x==xr);

    if (search != kmap.end()) {
      found2 = true;
      KmerEntry val = search->second;
      c2 = val.contig;
      if (forward== val.isFw()) {
        p2 = val.getPos() - kit2->second;
        d2 = true;
      } else {
        p2 = val.getPos() + k + kit2->second;
        d2 = false;
      }
      break;
    }
  }

  if (!found2) {
    return -1;
  }

  if (c1 != c2) {
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

// use:  match(s,l,v)
// pre:  v is initialized
// post: v contains all equiv classes for the k-mers in s
void KmerIndex::match(const char *s, int l, std::vector<std::pair<int, int>>& v) const {
  KmerIterator kit(s), kit_end;
  bool backOff = false;
  int nextPos = 0; // nextPosition to check
  for (int i = 0;  kit != kit_end; ++i,++kit) {
    // need to check it
    Kmer rep = kit->first.rep();
    auto search = kmap.find(rep);
    int pos = kit->second;

    if (search != kmap.end()) {

      KmerEntry val = search->second;
      
      v.push_back({dbGraph.ecs[val.contig], kit->second});

      // see if we can skip ahead
      // bring thisback later
      bool forward = (kit->first == rep);
      int dist = dbGraph.contigs[val.contig].getDist(val, forward);


      //const int lastbp = 10;
      if (dist >= 2) {
        // where should we jump to?
        int nextPos = pos+dist; // default jump

        if (pos + dist >= l-k) {
          // if we can jump beyond the read, check the middle
          nextPos = l-k;
          if (pos < (l-k)/2) {
            // if we are in the first half of the read
            nextPos = (l-k)/2 + 1;
          } else {
            nextPos = l-k;
          }
        }

        // check next position
        KmerIterator kit2(kit);
        kit2.jumpTo(nextPos);
        if (kit2 != kit_end) {
          Kmer rep2 = kit2->first.rep();
          auto search2 = kmap.find(rep2);
          if (search2 == kmap.end() || (val.contig == search2->second.contig)) {
            // great, a match (or nothing) see if we can move the k-mer forward
            if (pos + dist >= l-k) {
              v.push_back({dbGraph.ecs[val.contig], l-k}); // push back a fake position
              break; //
            } else {
              v.push_back({dbGraph.ecs[val.contig], pos+dist});
              kit = kit2; // move iterator to this new position
            }
          } else {
            // this is weird, let's try the middle k-mer
            bool foundMiddle = false;
            if (dist > 4) {
              int middlePos = (pos + nextPos)/2;
              int middleContig = -1;
              KmerIterator kit3(kit);
              kit3.jumpTo(middlePos);
              if (kit3 != kit_end) {
                Kmer rep3 = kit3->first.rep();
                auto search3 = kmap.find(rep3);
                if (search3 != kmap.end()) {
                  middleContig = search3->second.contig;
                  if (middleContig == val.contig) {
                    foundMiddle = true;
                  } else if (middleContig == search2->second.contig) {
                    foundMiddle = true;
                  }
                }
              }

              if (foundMiddle) {
                v.push_back({dbGraph.ecs[middleContig], nextPos});
                if (nextPos >= l-k) {
                  break;
                } else {
                  kit = kit2; 
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
          v.push_back({dbGraph.ecs[val.contig], l-k});
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
          Kmer rep = kit->first.rep();
          auto search = kmap.find(rep);
          if (search != kmap.end()) {
            // if k-mer found
            v.push_back({dbGraph.ecs[search->second.contig], kit->second}); // add equivalence class, and position
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

// use:  res = intersect(ec,v)
// pre:  ec is in ecmap, v is a vector of valid transcripts
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
