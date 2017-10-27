#include "GeneModel.h"
#include <assert.h>
#include <iostream>
#include <htslib/sam.h> // needed for CIGAR ops

char strandToChar(bool s) {
  return (s) ? '+' : '-';  
}



bool charToStrand(char c) {
  switch (c) {
    case '+':
      return true; break;
    case '-':
      return false; break;
    default:
      return true; break;
  }
}


bool Transcriptome::translateTrPosition(const int tr, const int pos, const int rlen, bool strand, TranscriptAlignment &aln) const {
  const TranscriptModel& model = transcripts[tr];
  aln.chr = model.chr;
  aln.cigar.clear();
  aln.chrpos = -1;

  aln.strand = (strand == model.strand);
  int trpos;
  int rpos = 0; // how many bp of read have been matched

  int n_exons = model.exons.size();
  
  if (model.strand) {
    trpos = pos;
    if (trpos < 0) { 
      // read starts a bit before the transcript, softclip to fit
      int softclip = -trpos;
      aln.cigar.push_back((softclip << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP);
      rpos += softclip;
      aln.chrpos = model.start;
    }
    for (int i = 0; i < n_exons; i++) {
      auto& exon = model.exons[i];
      int len = exon.stop - exon.start;      
      if (trpos < len) {
        if (rpos == 0) {
          aln.chrpos = exon.start + trpos; // left-most position
        }
        if (trpos + rlen <= len) {          
          aln.cigar.push_back(((rlen-rpos)<< BAM_CIGAR_SHIFT) | BAM_CMATCH); // rlen-trpos of M
          rpos = rlen;
          break; // end of the read
        } else {
          int mlen = 0;
          if (trpos < 0) { // match begins before this exon, extends over it
            mlen = len;            
          } else {     
            mlen = len - trpos;     
          }
          aln.cigar.push_back((mlen << BAM_CIGAR_SHIFT) | BAM_CMATCH);
          if (i +1 < n_exons) {
            aln.cigar.push_back(((model.exons[i+1].start-exon.stop) << BAM_CIGAR_SHIFT) | BAM_CREF_SKIP); // insert a bunch of N's
          }
          rpos += mlen;
        }
      }
      trpos -= len;
    }
    if (rpos < rlen) {
      aln.cigar.push_back(((rlen-rpos) << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP);
    }
  } else {
    trpos = model.length - pos - rlen; // counting from back from right most transcript position
    if (trpos < 0) {
      int softclip = -trpos;      
      aln.cigar.push_back((softclip << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP); // soft clip
      rpos += softclip;
      aln.chrpos = model.start;
    }
    for (int i = n_exons-1; i >= 0; i--) {
      auto& exon = model.exons[i];
      int len = exon.stop - exon.start;
      if (trpos < len) {
        if (rpos == 0) {
          aln.chrpos = exon.start + trpos;
        }
        if (trpos + rlen <= len) {
          aln.cigar.push_back(((rlen-rpos) << BAM_CIGAR_SHIFT) | BAM_CMATCH);
          rpos = rlen;
          break;
        } else {
          int mlen = 0;
          if (trpos < 0) { // match begins before this exon, extends over it
            mlen = len;            
          } else {     
            mlen = len - trpos;     
          }
          aln.cigar.push_back((mlen << BAM_CIGAR_SHIFT) | BAM_CMATCH);
          if (i > 0) {
            aln.cigar.push_back(((model.exons[i-1].start - exon.stop) << BAM_CIGAR_SHIFT) | BAM_CREF_SKIP);
          }
          rpos += mlen;
        }
      }
      trpos -= len;
    }
    if (rpos < rlen) {
      aln.cigar.push_back(((rlen-rpos) << BAM_CIGAR_SHIFT) | BAM_CSOFT_CLIP);
    }
  }


  return true;
}


/*
bool Transcriptome::translateTrPosition(const std::string &tr, const int _trpos, std::string &chr, int& chrpos, std::string &gene_id) {
  auto gid_it = trxToGeneId.find(tr);
  if (gid_it == trxToGeneId.end()) {
    return false;
  }
  gene_id = gid_it->second;
  auto g_it = genes.find(gene_id);
  if (g_it == genes.end()) {
    return false;
  }

  auto t_it = g_it->second.transcripts.find(tr);
  if (t_it == g_it->second.transcripts.end()) {
    return false;
  }

  auto &trxmodel = t_it->second;
  chr = trxmodel.chr;    
  

  bool fw = trxmodel.strand;
  int trpos = _trpos;
  chrpos = trxmodel.start;
  for (auto& exon : trxmodel.exons) {
    int len = (exon.strand ) ? (exon.stop - exon.start) : (exon.start - exon.stop);
    assert(len > 0);
    if (trpos < len) {
      // maps to this exon
      chrpos = exon.start;
      if (exon.strand) {
        chrpos += trpos;
      } else {
        chrpos -= trpos;
      }
      break;
    } else {
      trpos -= len;
    }
  }
  if (trpos > 0) {
    // goes beyond last exon, map to end
    chrpos = trxmodel.stop + trpos - 1;
  }
  return true;
}
*/

/*void writeTranscriptome(Transcriptome &transcriptome, std::ostream &out) {
  for(const auto &gene : transcriptome.genes) {
    const auto &gene_id = gene.first;
    const auto &model = gene.second;
    
    out << "GENE" << "\t" << model.id  << "\t"
        << model.name << "\t"
        << model.chr << "\t"
        << typeToString(model.type) << "\t"
        << strandToChar(model.strand) << "\t"
        << model.start << "\t"
        << model.stop << "\t";
    
    bool firsttr = true;
    for (auto &trlist : model.transcripts) {
      if (!firsttr) {
        out << ";";
      } else {
        firsttr = false;
      }
      out << trlist.first;
    }
    out << "\n";
  }
  for(const auto &gene : transcriptome.genes) {
    const auto &gene_id = gene.first;
    const auto &model = gene.second;
    for (auto &trlist : model.transcripts) {
      const auto &tr_id = trlist.first;
      const auto &tr = trlist.second;
      
      out << "TRANSCRIPT" << "\t" 
          << tr.id << "\t"
          << gene_id << "\t"
          << tr.chr << "\t"
          << typeToString(tr.type) << "\t"
          << strandToChar(tr.strand) << "\t"
          << tr.start << "\t"
          << tr.stop << "\t";
      
      bool firstex = true;
      for(const auto &exon : tr.exons) {
        if (!firstex) {
          out << ";";
        } else {
          firstex = false;
        }
        out << exon.start << "," << exon.stop;
      }
      out << "\n";
    }
  }
}
*/


void Transcriptome::loadTranscriptome(const KmerIndex& index, std::istream &in, const ProgramOptions& options) {
  std::string line;
  std::string segment;
  std::string type;
  std::string gtype, ttype;
  std::string gene_id, tr_id, chr_id;

  
  for (int i = 0; i < index.num_trans; i++) {
    TranscriptModel tr;
    tr.id = 1;
    tr.chr = -1;
    tr.gene_id = -1;
    transcripts.push_back(std::move(tr));
    trNameToId.insert({index.target_names_[i], i});
  }


  auto chrLookup = [&,this](const std::string chr) -> int {
    auto cit = chrNameToId.find(chr_id);
    if (cit != chrNameToId.end()) {      
      int c = cit->second;
      assert(c >= 0);
      //assert(c < chr.size());
      return c;
    } else {
      return -1;
    }  
  };


  //std::unordered_map<std::string,int>
  while (in.good()) {
    in >> type;
    if (type == "GENE") {
      GeneModel model;
      char strand = '?';
      in >> model.name >> model.commonName >> chr_id >> gtype >> strand >> model.start >> model.stop;
      model.chr = chrLookup(chr_id);
      if (model.chr == -1) {
        std::cerr << "Error: chromosome " << chr_id << " not definede before reference" << std::endl;
        assert(false);
      }      
      model.id = transcripts.size();      
      model.strand = charToStrand(strand);
      
      geneNameToId.insert({model.name, model.id});
      in >> line; // rest of transcripts, ignore for now
      genes.push_back(model);
      
    } else if (type == "TRANSCRIPT") {
      TranscriptModel model;
      char strand = '?';
      in >> tr_id >> gene_id >> chr_id >> ttype >> strand >> model.start >> model.stop;
      auto it = trNameToId.find(tr_id);
      if (it != trNameToId.end()) {
        model.id = it->second;        
      } else {
        std::cerr << "Warning: transcript " << tr_id << " is defined in GTF but not in FASTA" << std::endl;
        std::getline(in,line);
        continue;
        
      }


      model.chr = chrLookup(chr_id);
      if (model.chr == -1) {
        std::cerr << "Error: chromosome " << chr_id << " not definede before reference" << std::endl;        
        assert(false);
      }      

      model.strand = charToStrand(strand);    
      model.name = tr_id;
      
      in >> line;
      std::stringstream strline(line);
      while (std::getline(strline, segment, ';')) {
        ExonModel ex;
        ex.chr = model.chr;
        ex.strand = model.strand;
        int p = segment.find(',');
        ex.start = std::stoi(segment.substr(0,p));
        ex.stop = std::stoi(segment.substr(p+1));
        model.exons.push_back(std::move(ex));
      }
      model.length = 0;
      for (auto & exon : model.exons) {
        model.length += exon.stop - exon.start;
      }
      std::getline(in, line); 
      int id = model.id;
      if (transcripts[id].chr == -1) {
        transcripts[id] = std::move(model);
      }

    } else if (type == "CHROMOSOME") { 
      int id = chr.size();
      Chromosome c;      
      in >> c.name >> c.len;
      chr.push_back(c);
      chrNameToId.insert({c.name,id});
      std::getline(in,line);
    } else {
      std::getline(in, line); // read until end of line
    }
  }
}


/*void parseFasta(Transcriptome &transcriptome, const std::string &fasta_fn) {
  seqan::SeqFileIn seqFileIn(seqan::toCString(fasta_fn));
  seqan::CharString id;
  seqan::CharString seq;

  while(!atEnd(seqFileIn)) {
    seqan::readRecord(id,seq,seqFileIn);
    std::string name = std::string(seqan::toCString(id));
    size_t sp = name.find(' ');
    size_t pipe = name.find('|');
    seqan::toUpper(seq);
    transcriptome.seqs.insert({name.substr(0,std::min(sp,pipe)), std::move(seq)});
  }  
}

void parseGTF(Transcriptome &transcriptome, const std::string &gtf_fn, const ProgramOptions& options) {
  seqan::GffFileIn gtf(gtf_fn.c_str());
  seqan::GffRecord record;

  int n = 0;
  while(!seqan::atEnd(gtf)) {
    n++;
    seqan::readRecord(record,gtf);
    if (record.type == "gene") {
      GeneModel model;
      std::string gene_version;
      for (int i = 0; i < length(record.tagNames); i++) {
        if (record.tagNames[i] == "gene_id") {
          model.id = seqan::toCString(record.tagValues[i]);
        }
        if (record.tagNames[i] == "gene_name") {
          model.name = seqan::toCString(record.tagValues[i]);
        }
        if (record.tagNames[i] == "gene_biotype" || record.tagNames[i] == "gene_type") {
          std::string val = std::string(seqan::toCString(record.tagValues[i]));
          //auto &val = record.tagValues[i];
          if (val == "protein_coding") {
            model.type = BioType::PROTEIN;
          } else if (val.find("pseudogene") != std::string::npos) {
            model.type = BioType::PSEUDO;
          } else {
            model.type = BioType::OTHER;
          }
        }
        if (record.tagNames[i] == "gene_version") {
          gene_version = std::string(seqan::toCString(record.tagValues[i]));
        }
      }
      if (!gene_version.empty() && model.id.find('.') == std::string::npos) {        
        model.id += "." + gene_version;
      }
      if (model.name.empty()) {
        model.name = model.id;
      }
      if (options.ignoreProtein) {
        model.type = BioType::PROTEIN;
      }
      model.chr = seqan::toCString(record.ref);
      model.start = record.beginPos;
      model.stop = record.endPos;
      if (record.strand == '+') {
        model.strand = Strandedness::FORWARD;
      } else if (record.strand == '-') {
        model.strand = Strandedness::REVERSE;
      } else {
        model.strand = Strandedness::UNKNOWN;
      }
      transcriptome.genes.insert({model.id,std::move(model)});
    } else if (record.type == "transcript") {
      TranscriptModel model;
      bool bioTypeSet = false;
      model.type = BioType::OTHER; 
      std::string gene_id, gene_version, txp_version;
      for (int i = 0; i < length(record.tagNames); i++) {
        if (record.tagNames[i] == "gene_id") {
          gene_id = seqan::toCString(record.tagValues[i]);
        } else if (record.tagNames[i] == "transcript_id") {
          model.id = seqan::toCString(record.tagValues[i]);
        }
        if (record.tagNames[i] == "transcript_biotype" || record.tagNames[i] == "transcript_type") {
          std::string val = std::string(seqan::toCString(record.tagValues[i]));
          //auto &val = record.tagValues[i];
          if (val == "protein_coding") {
            model.type = BioType::PROTEIN;
          } else if (val.find("pseudogene") != std::string::npos) {
            model.type = BioType::PSEUDO;
          } 
          bioTypeSet = true;
        }   
        if (record.tagNames[i] == "gene_version") {
          gene_version = std::string(seqan::toCString(record.tagValues[i]));
        } 
        if (record.tagNames[i] == "transcript_version") {
          txp_version = std::string(seqan::toCString(record.tagValues[i]));
        }    
      }
      if (!gene_version.empty() && gene_id.find('.') == std::string::npos) {        
        gene_id += "." + gene_version;
      }
      if (!txp_version.empty() && model.id.find('.') == std::string::npos) {        
        model.id += "." + txp_version;
      }
      if (!bioTypeSet) {
        std::string source = seqan::toCString(record.source);
        // we need this for Ensembl versions 76 and below, 
        // transcripts don't have transcript_[bio]type set but store this info in the source name ?!?
        if (source == "protein_coding") {
          model.type = BioType::PROTEIN;
        } else if (source.find("pseudogene") != std::string::npos) {
          model.type = BioType::PSEUDO;
        }
      }
      if (options.ignoreProtein) {
        model.type = BioType::PROTEIN;
      }
      model.chr = seqan::toCString(record.ref);
      model.start = record.beginPos;
      model.stop = record.endPos;
      if (record.strand == '+') {
        model.strand = Strandedness::FORWARD;
      } else if (record.strand == '-') {
        model.strand = Strandedness::REVERSE;
      } else {
        model.strand = Strandedness::UNKNOWN;
      }
      assert(!gene_id.empty());
      auto it = transcriptome.genes.find(gene_id);
      assert(transcriptome.trxToGeneId.find(model.id) == transcriptome.trxToGeneId.end());
      transcriptome.trxToGeneId.insert({model.id,gene_id});
      assert(it != transcriptome.genes.end());
      it->second.transcripts.insert({model.id,model});      
    } else if (record.type == "exon") {
      ExonModel model;
      std::string trx_id;
      std::string gene_id;
      std::string txp_version, gene_version;
      for (int i = 0; i < length(record.tagNames); i++) {
         if (record.tagNames[i] == "gene_id") {
          gene_id = seqan::toCString(record.tagValues[i]);
        } else if (record.tagNames[i] == "transcript_id") {
          trx_id = seqan::toCString(record.tagValues[i]);
        }
        if (record.tagNames[i] == "gene_version") {
          gene_version = std::string(seqan::toCString(record.tagValues[i]));
        } 
        if (record.tagNames[i] == "transcript_version") {
          txp_version = std::string(seqan::toCString(record.tagValues[i]));
        }  
      }
      if (!gene_version.empty() && gene_id.find('.') == std::string::npos) {        
        gene_id += "." + gene_version;
      }
      if (!txp_version.empty() && trx_id.find('.') == std::string::npos) {        
        trx_id += "." + txp_version;
      }
      model.chr = seqan::toCString(record.ref);
      model.start = record.beginPos;
      model.stop = record.endPos;
      if (record.strand == '+') {
        model.strand = Strandedness::FORWARD;
      } else if (record.strand == '-') {
        model.strand = Strandedness::REVERSE;
      } else {
        model.strand = Strandedness::UNKNOWN;
      }

      auto g_it = transcriptome.genes.find(gene_id);
      assert(g_it != transcriptome.genes.end());
      auto t_it = g_it->second.transcripts.find(trx_id);
      assert(t_it != g_it->second.transcripts.end());
      t_it->second.exons.push_back(std::move(model));
    }
  }
  std::cerr << "GTF file contains " << transcriptome.genes.size() << " genes and " << transcriptome.trxToGeneId.size() << " transcripts" << std::endl;

  
}

*/