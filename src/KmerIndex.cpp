#include "KmerIndex.h"
#include <algorithm>
#include <random>
#include <ctype.h>
#include <zlib.h>
#include <unordered_set>
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
  
  // for each target, create it's own equivalence class
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
  

  std::cerr << "[build] counting k-mers ... "; std::cerr.flush();
  // gather all k-mers
  for (int i = 0; i < seqs.size(); i++) {
    const char *s = seqs[i].c_str();
    KmerIterator kit(s),kit_end;
    for (; kit != kit_end; ++kit) {
      kmap.insert({kit->first.rep(), KmerEntry()}); // don't care about repeats
      //    Kmer rep = kit->first.rep();
      // std::cout << rep.toString() << "\t" << rep.hash() << std::endl;
    }
  }
  std::cerr << "done." << std::endl;
  
  std::cerr << "[build] building target de Bruijn graph ... "; std::cerr.flush();
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
        assert(it->second.contig==-1);
        it->second = KmerEntry(contig.id, contig.length, i, forward);
        if (i > 0) {
          contig.seq.push_back(x.toString()[k-1]);
        }
      }
      
      dbGraph.contigs.push_back(contig);
      dbGraph.ecs.push_back(-1);
    }
  }
  std::cerr << " done " << std::endl;

}

void KmerIndex::BuildEquivalenceClasses(const ProgramOptions& opt, const std::vector<std::string>& seqs) {
  std::cerr << "[build] creating equivalence classes ... "; std::cerr.flush();

  std::vector<std::vector<TRInfo>> trinfos(dbGraph.contigs.size());
  //std::cout << "Mapping target " << std::endl;
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
        tr.sense = true;
        tr.start = val.getPos();
        if (contig.length - tr.start > seqlen - kit->second) {
          // tartget stops
          tr.stop = tr.start + seqlen - kit->second;
          jump = seqlen;
        } else {
          tr.stop = contig.length;
          jump = kit->second + (tr.stop - tr.start)-1;
        }
      } else {
        tr.sense = false;
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

      // // debugging -->
      //std::cout << "covering seq" << std::endl << seqs[i].substr(kit->second, jump-kit->second +k) << std::endl;
      //std::cout << "id = " << tr.trid << ", (" << tr.start << ", " << tr.stop << ")" << std::endl;
      //std::cout << "contig seq" << std::endl;
      // if (forward == val.isFw()) {
      //   //std::cout << contig.seq << std::endl;
      //   assert(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start) == seqs[i].substr(kit->second, jump-kit->second +k) );
      // } else {
      //   //std::cout << revcomp(contig.seq) << std::endl;
      //   assert(revcomp(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start)) == seqs[i].substr(kit->second, jump-kit->second +k));
      // }
      // if (jump == seqlen) {
      //   //std::cout << std::string(k-1+(tr.stop-tr.start)-1,'-') << "^" << std::endl;
      // }

      // // <-- debugging
      
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
  //std::cerr << "For " << dbGraph.contigs.size() << ", " << (dbGraph.contigs.size() - perftr) << " of them need to be split" << std::endl;
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
    int ec = -1;
    if (search != ecmapinv.end()) {
      // insert contig -> ec info
      ec = search->second;
    } else {
      ec = ecmapinv.size();
      ecmapinv.insert({u,ec});
      ecmap.push_back(u);
    }
    dbGraph.ecs[i] = ec;
    assert(ec != -1);
    
    // record the transc
    Contig& contig = dbGraph.contigs[i];
    contig.ec = ec;
    // correct ec of all k-mers in contig
    /*
    KmerIterator kit(contig.seq.c_str()), kit_end;
    for (; kit != kit_end; ++kit) {
      auto ksearch = kmap.find(kit->first.rep());
      //ksearch->second.ec = ec;
    }
    */
  }

  // map transcripts to contigs
  //std::cout << std::endl;
  for (int i = 0; i < seqs.size(); i++) {
    int seqlen = seqs[i].size() - k + 1; // number of k-mers
    // debugging
    std::string stmp;
    const char *s = seqs[i].c_str();
    ////std::cout << "sequence number " << i << std::endl;
    //std::cout << ">" << target_names_[i] << std::endl;
    //std::cout << seqs[i] << std::endl;
    KmerIterator kit(s), kit_end;
    for (; kit != kit_end; ++kit) {
      Kmer x = kit->first;
      //std::cout << "position = " << kit->second << ", mapping " << x.toString() << std::endl;
      Kmer xr = x.rep();
      auto search = kmap.find(xr);
      bool forward = (x==xr);
      KmerEntry val = search->second;
      Contig& contig = dbGraph.contigs[val.contig];

      ContigToTranscript info;
      info.trid = i;
      info.pos = kit->second;
      info.sense = (forward == val.isFw());
      int jump = kit->second + contig.length-1;
      //std::cout << "mapped to contig " << val.contig << ", len = " << contig.length <<  ", pos = " << val.getPos() << ", sense = " << info.sense << std::endl;
      contig.transcripts.push_back(info);
      // debugging
      if (info.sense) {

        if (info.pos == 0) {
          stmp.append(contig.seq);
          //std::cout << contig.seq << std::endl;
        } else {
          stmp.append(contig.seq.substr(k-1));
          //std::cout << contig.seq.substr(k-1) << std::endl;
        }
      } else {
        std::string r = revcomp(contig.seq);
        if (info.pos == 0) {
          stmp.append(r);
          //std::cout << r << std::endl;
        } else {
          stmp.append(r.substr(k-1));
          //std::cout << r.substr(k-1) << std::endl;
        }
      }
      //std::cout << stmp << std::endl;
      //std::cout << "covering seq" << std::endl << seqs[i].substr(kit->second, jump-kit->second +k) << std::endl;
      //std::cout << "id = " << tr.trid << ", (" << tr.start << ", " << tr.stop << ")" << std::endl;
      //std::cout << "contig seq" << std::endl;
      // if (forward == val.isFw()) {
      //   //std::cout << contig.seq << std::endl;
      //   assert(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start) == seqs[i].substr(kit->second, jump-kit->second +k) );
      // } else {
      //   //std::cout << revcomp(contig.seq) << std::endl;
      //   assert(revcomp(contig.seq.substr(tr.start, k-1 + tr.stop-tr.start)) == seqs[i].substr(kit->second, jump-kit->second +k));
      // }
      // if (jump == seqlen) {
      //   //std::cout << std::string(k-1+(tr.stop-tr.start)-1,'-') << "^" << std::endl;
      // }

      // // <-- debugging
      
      kit.jumpTo(jump);
    }
    if (seqlen > 0 && seqs[i] != stmp) {
      /*std::cout << ">" << target_names_[i] << std::endl
                << seqs[i] << std::endl
                << stmp << std::endl;*/
      assert(false);
      
    }
  }

  // double check the contigs
  for (auto &c : dbGraph.contigs) {
    for (auto info : c.transcripts) {
      std::string r;
      if (info.sense) {
        r = c.seq;
      } else {
        r = revcomp(c.seq);
      }
      assert(r == seqs[info.trid].substr(info.pos,r.size()));
    }
  }

  
  std::cerr << " done" << std::endl;
  std::cerr << "[build] target de Bruijn graph has " << dbGraph.contigs.size() << " contigs and contains "  << kmap.size() << " k-mers " << std::endl;


  

  
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
          search->second = KmerEntry(newc.id, newc.length,  kit->second, forward);
        }

        // repair tr-info
        std::vector<TRInfo> newtrinfo;
        for (auto x : oldtrinfo) {
          if (!(x.stop <= brpoints[j-1] || x.start >= brpoints[j])) {
            TRInfo trinfo;
            trinfo.sense = x.sense;
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


  //std::cerr << "For " << dbGraph.contigs.size() << ", " << (dbGraph.contigs.size() - perftr) << " of them need to be split" << std::endl;

  
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
  if (writeKmerTable) {
    assert(dbGraph.contigs.size() == dbGraph.ecs.size());
    tmp_size = dbGraph.contigs.size();
    out.write((char*)&tmp_size, sizeof(tmp_size));
    for (auto& contig : dbGraph.contigs) {
      out.write((char*)&contig.id, sizeof(contig.id));
      out.write((char*)&contig.length, sizeof(contig.length));
      tmp_size = strlen(contig.seq.c_str());
      out.write((char*)&tmp_size, sizeof(tmp_size));
      out.write(contig.seq.c_str(), tmp_size);

      // 10.1 write out transcript info
      tmp_size = contig.transcripts.size();
      out.write((char*)&tmp_size, sizeof(tmp_size));
      for (auto& info : contig.transcripts) {
        out.write((char*)&info.trid, sizeof(info.trid));
        out.write((char*)&info.pos, sizeof(info.pos));
        out.write((char*)&info.sense, sizeof(info.sense));
      }
    }
    
    // 11. write out ecs info
    for (auto ec : dbGraph.ecs) {
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
    std::cerr << "Error: incompatible indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
              << "Rerun with index to regenerate";
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

  kmap.clear();
  if (loadKmerTable) {
    kmap.reserve(kmap_size,true);
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

  // 9. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

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

    if (tmp_size + 1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }

    memset(buffer,0,bufsz);
    in.read(buffer, tmp_size);
    c.seq = std::string(buffer); // copy
    
    // 10.1 read transcript info
    in.read((char*)&tmp_size, sizeof(tmp_size));
    c.transcripts.clear();
    c.transcripts.reserve(tmp_size);

    for (auto j = 0; j < tmp_size; j++) {
      ContigToTranscript info;
      in.read((char*)&info.trid, sizeof(info.trid));
      in.read((char*)&info.pos, sizeof(info.pos));
      in.read((char*)&info.sense, sizeof(info.sense));
      c.transcripts.push_back(info);
    }

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
void KmerIndex::match(const char *s, int l, std::vector<std::pair<KmerEntry, int>>& v) const {
  KmerIterator kit(s), kit_end;
  bool backOff = false;
  int nextPos = 0; // nextPosition to check
  for (int i = 0;  kit != kit_end; ++i,++kit) {
    // need to check it
    auto search = kmap.find(kit->first.rep());
    int pos = kit->second;

    if (search != kmap.end()) {

      KmerEntry val = search->second;
      
      v.push_back({val, kit->second});

      // see if we can skip ahead
      // bring thisback later
      bool forward = (kit->first == search->first);
      int dist = val.getDist(forward);


      //const int lastbp = 10;
      if (dist >= 2) {
        // where should we jump to?
        int nextPos = pos+dist; // default jump

        if (pos + dist >= l-k) {
          // if we can jump beyond the read, check the end
          nextPos = l-k;
        }

        // check next position
        KmerIterator kit2(kit);
        kit2.jumpTo(nextPos);
        if (kit2 != kit_end) {
          Kmer rep2 = kit2->first.rep();
          auto search2 = kmap.find(rep2);
          bool found2 = false;
          int  found2pos = pos+dist;
          if (search2 == kmap.end()) {
            found2=true;
            found2pos = pos;
          } else if (val.contig == search2->second.contig) {
            found2=true;
            found2pos = pos+dist;
          }
          if (found2) {
            // great, a match (or nothing) see if we can move the k-mer forward
            if (found2pos >= l-k) {
              v.push_back({val, l-k}); // push back a fake position
              break; //
            } else {
              v.push_back({val, found2pos});
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
              kit3.jumpTo(middlePos);
              KmerEntry val3;
              if (kit3 != kit_end) {
                Kmer rep3 = kit3->first.rep();
                auto search3 = kmap.find(rep3);
                if (search3 != kmap.end()) {
                  middleContig = search3->second.contig;
                  if (middleContig == val.contig) {
                    foundMiddle = true;
                    found3pos = middlePos;
                  } else if (middleContig == search2->second.contig) {
                    foundMiddle = true;
                    found3pos = pos+dist;
                  }
                }


                if (foundMiddle) {
                  v.push_back({search3->second, found3pos});
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
          Kmer rep = kit->first.rep();
          auto search = kmap.find(rep);
          if (search != kmap.end()) {
            // if k-mer found
            v.push_back({search->second, kit->second}); // add equivalence class, and position
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

std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, int p) const {
  auto it = kmap.find(km.rep());
  if (it != kmap.end()) {
    KmerEntry val = it->second;
    return findPosition(tr, km, val, p);
  } else {
    return {-1,true};
  }
}

//use:  (pos,sense) = index.findPosition(tr,km,val,p)
//pre:  index.kmap[km] == val,
//      km is the p-th k-mer of a read
//      val.contig maps to tr
//post: km is found in position pos (1-based) on the sense/!sense strand of tr
std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, KmerEntry val, int p) const {
  bool fw = (km == km.rep());
  bool csense = (fw == val.isFw());

  int trpos = -1;
  bool trsense = true;
  if (val.contig < 0) {
    return {-1, true};
  }
  const Contig &c = dbGraph.contigs[val.contig];
  for (auto x : c.transcripts) {
    if (x.trid == tr) {
      trpos = x.pos;
      trsense = x.sense;
      break;
    }
  }

  if (trpos == -1) {
    return {-1,true};
  }


  if (trsense) {
    if (csense) {
      return {trpos + val.getPos() - p + 1, csense}; // 1-based, case I
    } else {
      return {trpos + val.getPos() + k + p, csense}; // 1-based, case III
    }
  } else {
    if (csense) {
      return {trpos + (c.length - val.getPos() -1) + k + p, !csense};  // 1-based, case IV
    } else {
      return {trpos + (c.length - val.getPos())  - p, !csense}; // 1-based, case II
    }
  }
}

/*
// use:  r = matchEnd(s,l,v,p)
// pre:  v is initialized, p>=0
// post: v contains all equiv classes for the k-mer in s, as
//       well as the best match for s[p:]
//       if r is false then no match was found
bool KmerIndex::matchEnd(const char *s, int l, std::vector<std::pair<int,int>> &v, int maxPos) const {
  // kmer-iterator checks for N's and out of bounds
  KmerIterator kit(s+maxPos), kit_end;
  if (kit != kit_end && kit->second == 0) {
    Kmer last = kit->first;
    auto search = kmap.find(last.rep());

    if (search == kmap.end()) {
      return false; // shouldn't happen
    }

    KmerEntry val = search->second;
    bool forward = (kit->first == search->first);
    int dist = val.getDist(forward);
    int pos = maxPos + dist + 1; // move 1 past the end of the contig
    
    const char* sLeft = s + pos + (k-1);
    int sizeleft = l -  (pos + k-1); // characters left in sleft
    if (sizeleft <= 0) {
      return false; // nothing left to match
    }

    // figure out end k-mer
    const Contig& root = dbGraph.contigs[val.contig];
    Kmer end; // last k-mer
    bool readFw = (forward == val.isFw());
    if (readFw) {
      end = Kmer(root.seq.c_str() + root.length-1);
    } else {
      end = Kmer(root.seq.c_str()).twin();
    }

    
    int bestContig = -1;
    int bestDist = sizeleft;
    int numBest = 0;
    for (int i = 0; i < 4; i++) {
      Kmer x = end.forwardBase(Dna(i));
      Kmer xr = x.rep();
      auto searchx = kmap.find(xr);
      if (searchx != kmap.end()) {
        KmerEntry valx = searchx->second;
        const Contig& branch = dbGraph.contigs[valx.contig];
        if (branch.length < sizeleft) {
          return false; // haven't implemented graph walks yet
        }
        std::string cs;
        bool contigFw = (x==xr) && valx.isFw();
        // todo: get rid of this string copy
        if (valx.getPos() == 0 && contigFw) {
          cs = branch.seq.substr(k-1);
        } else if (valx.getPos() == branch.length-1 && !contigFw) {
          cs = revcomp(branch.seq).substr(k-1);
        }
        int hdist = hamming(sLeft,cs.c_str());
        if (bestDist >= hdist) {
          numBest++;
          if (bestDist > hdist) {
            bestDist = hdist;
            numBest = 1;
          }
          bestContig = valx.contig;
        }
      }
    }

    if (numBest == 1 && bestDist < 10) {
      v.push_back({dbGraph.ecs[bestContig], l-k});
      return true;
    } else {
      return false;
    }
  } else {
    return false;
  }
}
*/


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


  
  std::vector<std::vector<std::pair<int, ContigToTranscript>>> trans_contigs(num_trans);
  for (auto &c : dbGraph.contigs) {
    for (auto &ct : c.transcripts) {
      trans_contigs[ct.trid].push_back({c.id, ct});
    }
  }

  auto &target_seqs = const_cast<std::vector<std::string>&>(target_seqs_);
  
  for (int i = 0; i < trans_contigs.size(); i++) {
    auto &v = trans_contigs[i];
    std::sort(v.begin(), v.end(), [](std::pair<int,ContigToTranscript> a, std::pair<int,ContigToTranscript> b) {
        return a.second.pos < b.second.pos;
      });

    std::string seq;
    seq.reserve(target_lens_[i]);

    for (auto &pct : v) {
      auto ct = pct.second;
      int start = (ct.pos==0) ? 0 : k-1;
      const auto& contig = dbGraph.contigs[pct.first];
      if (ct.sense) {
        seq.append(contig.seq.substr(start));
      } else {
        seq.append(revcomp(contig.seq).substr(start));
      }
    }
    target_seqs.push_back(seq);
  }

  bool &t = const_cast<bool&>(target_seqs_loaded);
  t = true;//target_seqs_loaded = true;
  return;
}

void KmerIndex::clear() {
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
