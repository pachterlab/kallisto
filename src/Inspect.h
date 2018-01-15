#ifndef KALLISTO_INSPECTINDEX_H
#define KALLISTO_INSPECTINDEX_H

#include <iostream>
#include <unordered_set>
#include "KmerIndex.h"
#include "GeneModel.h"

using namespace std;


struct ECStruct {
  int ec;
  int chr;
  int start;
  int stop;
  std::vector<std::pair<int,int>> start_lens;
  std::vector<int> tlist;
};

std::vector<ECStruct> merge_contigs(std::vector<ECStruct> ecv) {
  if (ecv.size() <= 1) {
    return (ecv);
  }
  std::vector<ECStruct> out;
  std::sort(ecv.begin(), ecv.end(), [&](const ECStruct& a, const ECStruct& b) { return a.start < b.start;});

  // check if overlapping

  int a = 0, b = 0;
  while (b <= ecv.size()) {
    assert(a <= b);
    if (a == ecv.size()) {
      break;
    }

    // ecv[a:b] can be merged, see if ecv[a:(b+1)] can be
    if (b < ecv.size()) {
      if (a == b) {
        b++;
        continue; // ok, trivial to merge empty set
      }

      if (ecv[b-1].stop <= ecv[b].start) {           
        if (b+1 == ecv.size() || ecv[b].stop <= ecv[b+1].start) {
          b++;
          continue;
        }        
      }
    }
    // ok, push back ecv[a:b]
    if (a+1 == b) {
      out.push_back(ecv[a]);
    } else {
      std::unordered_set<int> tset;
      ECStruct ecs;
      ecs.ec = ecv[a].ec;
      ecs.chr = ecv[a].ec;
      ecs.start = ecv[a].start;
      ecs.stop = ecv[b-1].stop;
      int pos,len;
      for (int i = a; i < b; i++) {
        // len = 0;        
        pos = ecv[i].start - ecs.start;
        auto& sp = ecs.start_lens;
        for (auto x : ecv[i].start_lens) {
          sp.push_back({pos+x.first, x.second});
          // len += x.second;
        }       
        tset.insert(ecv[i].tlist.begin(), ecv[i].tlist.end());
      }
      for (auto t : tset) {
        ecs.tlist.push_back(t);
      }
      std::sort(ecs.tlist.begin(), ecs.tlist.end());
      out.push_back(ecs);
    }
    // 
    a = b;
  }

  return out;
}

void printVector(const vector<int>& v) {
  cout << "[";
  int i = 0;
  for (auto x : v) {
    if (i>0) {
      cout << ", ";
    }
    cout << x;
    i++;
  }
  cout << "]";
}

void printHisto(const unordered_map<int,int>& m, const string& header) {
  cout << header << "\n";
  int mn = std::numeric_limits<int>::max();
  int mx = 0;

  for (auto kv : m) {
    mn = min(mn,kv.first);
    mx = max(mx,kv.first);
  }

  for (int i = mn; i <= mx; i++) {
    auto search = m.find(i);
    if (search == m.end()) {
      cout << i << "\t0\n";
    } else {
      cout << i << "\t" << search->second << "\n";
    }
  }
}

void InspectIndex(const KmerIndex& index, const ProgramOptions& opt) {

  std::string gfa = opt.gfa;
  std::string bed = opt.bedFile;

  static const char *dna = "ACGT";
  auto Dna = [](int i) {return dna[i & 0x03];};

  int k = index.k;
  cout << "#[inspect] Index version number = " << index.INDEX_VERSION << endl;
  cout << "#[inspect] k = " << index.k << endl;;
  cout << "#[inspect] number of targets = " << index.num_trans << endl;

  cout << "#[inspect] number of equivalence classes = " << index.ecmap.size() << endl;


  if (index.ecmap.size() != index.ecmapinv.size()) {
    cout << "Error: sizes do not match. ecmap.size = " << index.ecmap.size()
         << ", ecmapinv.size = " << index.ecmapinv.size() << endl;
    exit(1);
  }

  if (index.dbGraph.ecs.size() != index.dbGraph.contigs.size()) {
    cout << "Error: sizes do not match. ecs.size = " << index.dbGraph.ecs.size()
         << ", contigs.size = " << index.dbGraph.contigs.size() << endl;
    exit(1);
  }

  cout << "#[inspect] number of contigs = " << index.dbGraph.contigs.size() << endl;
  

  unordered_map<int,int> echisto;

  //for (auto& ecv : index.ecmap) {
  for (int ec = 0; ec < index.ecmap.size(); ec++) {
    const vector<int>& v = index.ecmap[ec];
    ++echisto[v.size()];

    if (opt.inspect_thorough) {
      if (v.empty()) {
        cout << "Error: ec = " << ec  << " is empty!" << endl;
        exit(1);
      }
      for (int i = 0; i < v.size(); i++) {
        if (v[i] < 0 || v[i] >= index.num_trans) {
          cout << "Error: ec = " << ec  << " has invalid target id " << v[i] << endl;
          exit(1);
        }

        if (i > 0 && v[i] == v[i-1]) {
          cout << "Error: ec = " << ec  << " has repeated target id " << v[i] << endl;
          exit(1);
        }

        if (i > 0 && v[i] < v[i-1]) {
          cout << "Error: ec = " << ec  << " is not sorted!" << endl;
          exit(1);
        }
      }

      auto search = index.ecmapinv.find(v);
      if (search == index.ecmapinv.end()) {
        cout << "Error: could not find inverse for " << ec << endl;
        exit(1);
      } else {
        if (search->second != ec) {
          cout << "Error: inverse incorrect for ecmap -> ecmapinv,  ecv.first = "
              << ec <<  ", ecmapinv[ecv.second] = " << search->second << endl;
          exit(1);
        }
      }
    }
  }

  if (opt.inspect_thorough) {
    for (auto& eiv : index.ecmapinv) {
      //auto search = index.ecmap.find(eiv.second);
      //if (search == index.ecmap.end()) {
      if (eiv.second < 0 || eiv.second >= index.ecmap.size()) {
        cout << "Error: could not find inverse for ";
        printVector(eiv.first);
        cout << ", ecid = " << eiv.second << endl;
        exit(1);
      } else {
        auto &v = index.ecmap[eiv.second];
        if (v != eiv.first) {
          cout << "Error: inverse incorrect for ecmapinv -> ecmap,  eiv.first = ";
          printVector(eiv.first);
          cout <<  ", ecmap[eiv.second] = ";
          printVector(v);
          cout << endl;
          exit(1);
        }
      }
    }
  }

  cout << "#[inspect] Number of k-mers in index = " << index.kmap.size() << endl;
  unordered_map<int,int> kmhisto;

  for (auto& kv : index.kmap) {
    int id = kv.second.contig;
    int pos = kv.second.getPos();
    int fw = kv.second.isFw();

    if (id < 0 || id >= index.dbGraph.contigs.size()) {
      cerr << "Kmer " << kv.first.toString() << " mapped to contig " << id << ", which is not in the de Bruijn Graph" << endl;
      exit(1);
    } else {
      ++kmhisto[index.ecmap[index.dbGraph.ecs[id]].size()];
    }

    if (opt.inspect_thorough) {
      const Contig& c = index.dbGraph.contigs[id];
      const char* s = c.seq.c_str();
      Kmer x = Kmer(s+pos);
      Kmer xr = x.rep();

      bool bad = (fw != (x==xr)) || (xr != kv.first);
      if (bad) {
        cerr << "Kmer " << kv.first.toString() << " mapped to contig " << id << ", pos = " << pos << ", on " << (fw ? "forward" : "reverse") << " strand" << endl;
        cerr << "seq = " << c.seq << endl;
        cerr << "x  = " << x.toString() << endl;
        cerr << "xr = " << xr.toString() << endl;
        exit(1);
      }
    }
  }

  if (opt.inspect_thorough) {
    for (int i = 0; i < index.dbGraph.contigs.size(); i++) {
      const Contig& c = index.dbGraph.contigs[i];

      if (c.seq.size() != c.length + k-1) {
        cerr << "Length and string dont match " << endl << "seq = " << c.seq << " (length = " << c.seq.size() << "), c.length = " << c.length << endl;
        exit(1);
      }


      const char *s = c.seq.c_str();
      KmerIterator kit(s), kit_end;
      for (; kit != kit_end; ++kit) {
        Kmer x = kit->first;
        Kmer xr = x.rep();
        auto search = index.kmap.find(xr);
        if (search == index.kmap.end()) {
          cerr << "could not find kmer " << x.toString() << " in map " << endl << "seq = " << c.seq << ", pos = " << kit->second << endl;
          exit(1);
        }

        KmerEntry val = search->second;
        if (val.contig != i /*|| val.ec != index.dbGraph.ecs[i]*/ || val.contig_length != c.length || val.getPos() != kit->second || val.isFw() != (x==xr)) {
          cerr << "mismatch " << x.toString() << " in map " << endl << "id = " << i << ", ec = " << index.dbGraph.ecs[i] << ", length = " << c.length << ", seq = " << c.seq << ", pos = " << kit->second << endl;
          cerr << "val = " << val.contig << /* ", ec = " << val.ec << */ ", length = " << val.contig_length << ", pos = (" << val.getPos() << ", " << (val.isFw() ? "forward" :  "reverse") << ")" << endl;
          exit(1);
        }
      }
      
    }
  }

  if (gfa.empty() && bed.empty()) {
    printHisto(echisto, "#EC.size\tNum.targets");
    cout << endl << endl;

    printHisto(kmhisto, "#EC.size\tNum.kmers");
  }


  if (!gfa.empty()) {
    std::ofstream out;
    out.open(gfa);
    out << "H\tVN:Z:1.0\n";
    int i = 0;
    for (auto& c : index.dbGraph.contigs) {
      out << "S\t" << i << "\t" << c.seq << "\tXT:S:";
      for (int j = 0; j < c.transcripts.size(); j++) {
        auto &ct = c.transcripts[j];
        if (j > 0) {
          out << ",";
        }
        out << index.target_names_[ct.trid];
      }
      out << "\n";
      i++;
    }

    const auto& kmap = index.kmap;
    i = 0;
    for (auto& c : index.dbGraph.contigs) {
      auto& seq = c.seq;

      Kmer last(seq.c_str() + seq.size()-k);
      for (int j = 0; j < 4; j++) {
        Kmer after = last.forwardBase(Dna(j));
        auto search = kmap.find(after.rep());
        if (search != kmap.end()) {
          KmerEntry val = search->second;
          // check if + or -
          bool strand = val.isFw() == (after == after.rep());
          out << "L\t" << i << "\t+\t" << val.contig
              << "\t" << (strand ? '+' : '-') << "\t" << (k-1) << "M\n";
        }
      }

      // enumerate bw links
      Kmer first(seq.c_str());
      for (int j = 0; j < 4; j++) {
        Kmer before = first.backwardBase(Dna(j));
        auto search = kmap.find(before.rep());
        if (search != kmap.end()) {
          KmerEntry val = search->second;
          // check if + or -
          bool strand = val.isFw() == (before == before.rep());
          out << "L\t" << i << "\t-\t"
              << val.contig << "\t" << (strand ? '-' : '+')
              << "\t" << (k-1) << "M\n";
        }
      }

      i++;
    }

    out.flush();
    
    out.close();
  }

  if (!bed.empty()) {
    // export bed track with TCC information    
    bool guessChromosomes = false;
    Transcriptome model;
    if (opt.genomebam) {          
      if (!opt.chromFile.empty()) {
        model.loadChromosomes(opt.chromFile);
      } else {
        guessChromosomes = true;
      }          
      model.parseGTF(opt.gtfFile, index, opt, guessChromosomes);
      //model.loadTranscriptome(index, in, opt);
    }

    std::ofstream out;
    out.open(bed);

    out << "track name=\"Kallisto \" gffTags=\"on\"\n";
    std::unordered_map<TranscriptAlignment, std::vector<int>> cmap;
    cmap.reserve(100);
    std::vector<std::unordered_map<int, std::vector<ECStruct>>> ec_chrom(index.ecmap.size());

    for (const auto& c : index.dbGraph.contigs) {
      cmap.clear();
      // structure for TRaln
      TranscriptAlignment tra;
      
      int len = c.length;
      int cid = c.id;
      for (const auto& ct : c.transcripts) {
        // ct.trid, ct.pos, ct.sense
        model.translateTrPosition(ct.trid, ct.pos, len, ct.sense, tra);
        cmap[tra].push_back(ct.trid);
      }

      for (const auto& cp : cmap) {
        const auto& tra = cp.first;
        const auto& tlist = cp.second;      
        if (tra.chr != -1) {          
          ECStruct ecs;
          ecs.chr = tra.chr;
          ecs.ec = index.dbGraph.ecs[c.id];
          ecs.start = tra.chrpos;
          int pos = 0;
          for (uint32_t cig : tra.cigar) {
            int len = cig >> BAM_CIGAR_SHIFT;
            int type = cig & 0xF;
            if (type == BAM_CMATCH) {            
              ecs.start_lens.push_back({pos, len});
              pos += len;
            } else if (type == BAM_CREF_SKIP) {
              pos += len;
            } else {
              assert(false);
            }
          }
          ecs.stop = tra.chrpos + pos;
          ecs.tlist = tlist;    
          
          ec_chrom[ecs.ec][tra.chr].push_back(ecs);
        }
      }
    }
    
    int num_ecs = ec_chrom.size();
    for (int ec = 0; ec < num_ecs; ec++) {
      auto& chrmap = ec_chrom[ec];
      if (chrmap.empty()) {
        continue;
      }
      bool single_chrom = (chrmap.size() == 1);
      for (auto& x : chrmap) {
        int chr = x.first;
        auto& ecv = x.second;

        auto blocklist = merge_contigs(ecv);

        for (auto &ecs : blocklist) {
          int i;
          out << model.chr[chr].name << "\t" // chrom
              << ecs.start << "\t" << ecs.stop << "\t"; // start stop
          out << "Name=" << ec << ";Transcripts=";
          i = 0;
          for (const auto& t : index.ecmap[ec]) {
            bool ischr = (std::find(ecs.tlist.begin(), ecs.tlist.end(), t) != ecs.tlist.end());
            if (i++ > 0) {
              out << "%0A";
            }            
            out << index.target_names_[t];
            if (!ischr) {
              int tchr = model.transcripts[t].chr;
              if (tchr != -1) {
                if ( tchr != chr) {
                  out << "%20(" << model.chr[tchr].name << ")";
                } else {
                  out << "%20(*)";
                }
              } else {
                out << "%20(\?\?\?)";
              }
            }

          }
          out  << ";\t" <<  (int)(1000.0 / index.ecmap[ec].size()) << "\t" // score
               << ".\t" // strand
               << ecs.start << "\t" << ecs.stop // thick part
               << "\t0\t"; // color name
          const auto& sp = ecs.start_lens;
          out << sp.size() << "\t";
          i = 0;
          for (const auto& x : sp) {
            if (i++ > 0) {
              out << ',';
            }
            out << x.second;
          }
          out << "\t";
          i = 0;
          for (const auto& x : sp) {
            if (i++ > 0) {
              out << ',';
            }
            out << x.first;
          }
          out << "\n";
        }
      }
    }
  }
}

#endif // KALLISTO_INSPECTINDEX_H
