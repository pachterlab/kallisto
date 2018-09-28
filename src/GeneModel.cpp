#include "GeneModel.h"
#include <assert.h>
#include <iostream>
#include <htslib/sam.h> // needed for CIGAR ops
#include <zlib.h>

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

int chrLookup(const Transcriptome& model, const std::string chr) {
  auto cit = model.chrNameToId.find(chr);
  if (cit != model.chrNameToId.end()) {      
    int c = cit->second;
    assert(c >= 0);
    //assert(c < chr.size());
    return c;
  } else {
    return -1;
  }  
};

bool Transcriptome::translateTrPosition(const int tr, const int pos, const int rlen, bool strand, TranscriptAlignment &aln) const {
  const TranscriptModel& model = transcripts[tr];
  if (model.chr == -1) {
    return false;
  }
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

void Transcriptome::loadChromosomes(const std::string &chrom_fn) {
  std::ifstream in(chrom_fn);
  std::string type;
  while (in.good()) {
    int id = chr.size();
    Chromosome c;      
    c.len = -1;
    in >> c.name >> c.len;
    if (!c.name.empty() && c.len >= 0) {
      chr.push_back(c);
      chrNameToId.insert({c.name,id});
    }
    std::getline(in,type); // clear end of line   
  }
}

/*
void Transcriptome::loadTranscriptome(const KmerIndex& index, std::istream &in, const ProgramOptions& options) {
  std::string line;
  std::string segment;
  std::string type;
  std::string gtype, ttype;
  std::string gene_id, tr_id, chr_id;

  
  for (int i = 0; i < index.num_trans; i++) {
    TranscriptModel tr;
    tr.id = i;
    tr.chr = -1;
    tr.gene_id = -1;
    tr.length = -1;
    tr.strand = true;
    transcripts.push_back(std::move(tr));
    trNameToId.insert({index.target_names_[i], i});
  }


  

  int tr_extras = 0;
  //std::unordered_map<std::string,int>
  while (in.good()) {
    in >> type;
    if (type == "GENE") {
      GeneModel model;
      char strand = '?';
      in >> model.name >> model.commonName >> chr_id >> gtype >> strand >> model.start >> model.stop;
      model.chr = chrLookup(*this, chr_id);
      if (model.chr == -1) {
        std::cerr << "Error: chromosome " << chr_id << " not defined before reference" << std::endl;
        assert(false);
      }      
      model.id =  genes.size();
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
        tr_extras++;
        //std::cerr << "Warning: transcript " << tr_id << " is defined in GTF but not in FASTA" << std::endl;
        std::getline(in,line);
        continue;        
      }

      model.chr = chrLookup(*this, chr_id);
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

  if (tr_extras > 0) {
    std::cerr << "Warning: " << tr_extras << " transcript defined in GTF but not found in FASTA" << std::endl;
  }

  int tr_bad = 0;
  for (int i = 0; i < index.num_trans; i++) {
    const auto& tr = transcripts[i];
    if (tr.chr == -1) {
      tr_bad++;
    }
  }
  if (tr_bad > 0) {
    std::cerr << "Warning: there were " << tr_bad << " transcripts out of " 
              << index.num_trans << " that are missing in the GTF file " << std::endl;
  }

}
*/


int Transcriptome::addGTFLine(const std::string &line, const KmerIndex& index, bool guessChromosomes) {

  if(line.empty() || line[0] == '#') {
    return 0;
  }
  int p = 0, t=0; 
  // read chr
  std::string schr;
  t = line.find('\t',p);
  schr = line.substr(p,t-p);
  t = line.find('\t',t+1); // skip annotation source
  p = t+1;
  t = line.find('\t',p);
  std::string typestr = line.substr(p,t-p);
  enum Type {GENE, TRANSCRIPT, EXON, OTHER};
  Type type;
  if (typestr == "gene") {
    type = Type::GENE;
    // need id, name, chr, coding, strand, start, stop, transcripts
  } else if (typestr == "transcript") {
    type = Type::TRANSCRIPT;
    //id, gene, chr, coding, strand, start, stop, exons: start,stop; (increasing for +, decreasing for -)
  } else if (typestr == "exon") {
    type = Type::EXON;
  } else {
    type = Type::OTHER;
    return 0;
  }
  p = t+1;
  t = line.find('\t',p);
  int start = std::stoi(line.substr(p,t-p))-1;
  p = t+1;
  t = line.find('\t',p);
  int stop = std::stoi(line.substr(p,t-p));
  t = line.find('\t',t+1);
  p = t+1; // skip score
  t = line.find('\t',p);
  char strand = line[p];
  p = t+1;
  t = line.find('\t',p); // phase, ignore
  p = t+1;

  GeneModel gmodel;
  TranscriptModel tmodel;
  ExonModel emodel;

  int ichr = chrLookup(*this, schr);
  if (ichr == -1) {
    if (guessChromosomes) {
      // just add a new chr on the fly
      int i = chr.size();
      Chromosome c;
      c.name = schr;
      c.len = 536870911; // maximum that bai can index :(
      chr.push_back(c);
      chrNameToId.insert({schr, i});
    } else {
      return 1; // couldn't find chrom
    }
  }

  if (type == Type::GENE) {
    gmodel.id = genes.size(); // ?? ever used?
    gmodel.chr = ichr;
    gmodel.start = start;
    gmodel.stop = stop;
    gmodel.strand = (strand == '+') ? true : false;
    gmodel.id = -1; // figure out later
  } else if (type == Type::TRANSCRIPT) {
    tmodel.chr = ichr;
    tmodel.start = start;
    tmodel.stop = stop;
    tmodel.strand = (strand == '+') ? true : false;    
  } else if (type == Type::EXON) {
    emodel.chr = ichr;
    emodel.start = start;
    emodel.stop = stop;
    emodel.strand = (strand == '+') ? true : false;
  }


  // line[p:] contains the additional fields as 'key "value";'
  // line[p:t] is key, line[t+2:]
  int s = 0;
  std::string key,value;
  int keycount = 0;
  std::string gversion, tversion, gene_name, transcript_name;
  while (p != std::string::npos) {    
    if ((t = line.find('"',p))== std::string::npos) {
      break;
    }
    if ((s = line.find('"',t+1)) == std::string::npos) {
      break;
    }  
    assert(line[t-1] == ' ');
    key.assign(line.substr(p,t-p-1));
    assert(s > t);
    value.assign(line.substr(t+1,s-t-1));

    // common values
    if (key == "gene_id") {
      keycount++;
      gene_name = std::move(value);
    } else if (key == "gene_version") {
      keycount++;
      gversion = std::move(value);
    }


    if (type == Type::GENE) {
      if (key == "gene_name") {
        keycount++;
        gmodel.commonName = std::move(value);
      }

      if (keycount == 3) {
        break;
      }
    } else {
      if (key == "transcript_id") { 
        keycount++;
        transcript_name = std::move(value);        
      } else if (key == "transcript_version") {
        keycount++;
        tversion = std::move(value);
      } 

      if (type == Type::TRANSCRIPT) {
        if (keycount == 4) {
          break;
        }
      } else if (type == Type::EXON) {
        if (keycount == 4) {
          break;
        }
      }
    }

    assert(s +1 < line.size());
    assert(line[s+1] == ';');
    if ((p = line.find(' ',s)) != std::string::npos) {
      p++;
      if (p >= line.size()) {
        break;
      }
    }
    
  }

  if (type == Type::GENE) {
    assert(!gene_name.empty());
    if (!gversion.empty() && gmodel.name.find('.') == std::string::npos) {
      gmodel.name += "." + gversion;
    }
    // add to the transcriptome model
    geneNameToId.insert({gmodel.name, gmodel.id});
    genes.push_back(std::move(gmodel));
  } else if (type == Type::TRANSCRIPT) {
    tmodel.name = transcript_name;
    assert(!tmodel.name.empty());
    
    // canonical name includes version number
    if (!tversion.empty() && tmodel.name.find('.') == std::string::npos) {
      tmodel.name += "." + tversion;
    }
    auto it2 = trNameToId.find(tmodel.name); 
    if (it2 == trNameToId.end()) {
      tmodel.name = transcript_name; // try without version number
      it2 = trNameToId.find(tmodel.name);
    }
    
    
    if (it2 != trNameToId.end()) {
      tmodel.id = it2->second;
      tmodel.length = index.target_lens_[tmodel.id];
    } else {
      // transcript not found, do nothing
      return 2;
    }

    std::string tmp_genename = gene_name; // copy
    if (!gversion.empty() && gmodel.name.find('.') == std::string::npos) {
      gene_name += "." + gversion;
    }
    auto it = geneNameToId.find(gene_name);
    if (it == geneNameToId.end()) {
      gene_name = tmp_genename;
      it = geneNameToId.find(gene_name);
    }

    if (it != geneNameToId.end()) {
      tmodel.gene_id = it->second;
    }

    int id = tmodel.id;
    if (transcripts[id].chr == -1) {
      transcripts[id] = std::move(tmodel);
    }
  } else if (type == Type::EXON) {
    std::string tmp_transcript_name = transcript_name;
    if (!tversion.empty() && transcript_name.find('.') == std::string::npos) {
      transcript_name += "." + tversion;
    }

    auto it = trNameToId.find(transcript_name);
    if (it == trNameToId.end()) {
      transcript_name = tmp_transcript_name;
      it = trNameToId.find(transcript_name);
    }

    if (it != trNameToId.end()) {
      auto& tm = transcripts[it->second];
      if (tm.chr != -1) {
        tm.exons.push_back(std::move(emodel));
      }
    }
  }
  


  return 0;
}

void Transcriptome::parseGTF(const std::string &gtf_fn, const KmerIndex& index, const ProgramOptions& options, bool guessChromosomes) {

  for (int i = 0; i < index.num_trans; i++) {
    TranscriptModel tr;
    tr.id = i;
    tr.chr = -1;
    tr.gene_id = -1;
    tr.strand = true;
    transcripts.push_back(std::move(tr));
    trNameToId.insert({index.target_names_[i], i});
  }


  int num_chrom_missing = 0, num_trans_missing = 0;
  int buf_size = 1<<22; // line is < 4M
  char* buf = new char[buf_size+1];
  buf[buf_size] = 0;
  int buf_end = 0;
  int bufread = 0;
  int pos = 0;
  char *pos_end = buf;
  gzFile file = gzopen(gtf_fn.c_str(), "r");
  std::string line;
  if (!file) {
    std::cerr << "Error: could not open file " << gtf_fn << std::endl;    
    return;
  }
  
  bool done_reading = false;
  while (!done_reading) {
    // buf[0:pos] contains a previously unprocessed line
    int to_read = (buf+buf_size) - pos_end;
    bufread = gzread(file, pos_end, to_read);
    buf_end = (pos_end - buf) + bufread;
    if (bufread < to_read) {
      if (gzeof(file)) {
        done_reading = true;
      }
    }
    while (true) {
      pos_end = std::strchr(buf+pos, '\n');
      if (pos_end == nullptr) {
        pos_end = buf + buf_end+1;
        if (!done_reading) {
          break;
        }
      }

      // line goes from buf[pos:pos_end]
      line.assign(buf+pos, pos_end - (buf+pos));
      int r = addGTFLine(line, index, guessChromosomes);
      if (r == 1) {
        num_chrom_missing++;
      } else if (r==2) {
        num_trans_missing++;
      }

      pos = pos_end - buf + 1;
      if (pos >= buf_end) {
        break;
      }
      assert(0 <= pos);
      assert(pos < buf_size);
    }
    
    // buf[pos:bufend] contains now new line,
    if (done_reading) { 
      break;
    }
    int leftover = buf_end - pos;
    if (leftover > 0) {
      std::memmove(buf, buf+pos,leftover);
    }
    size_t toset = (buf_size - leftover);
    std::memset(buf+leftover, 0, toset);
    pos_end = buf + leftover;
    pos = 0;

  }
  

  
  delete[] buf;
  gzclose(file);

  if (num_chrom_missing > 0) {
    std::cerr << "Warning: could not find chromosomes for " << num_chrom_missing << " transcripts" << std::endl;
  } 
  if (num_trans_missing > 0) {
    std::cerr << "Warning: " << num_trans_missing << " transcripts were defined in GTF file, but not in the index" << std::endl;
  }
}
