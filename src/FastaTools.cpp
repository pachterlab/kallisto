//
//  FastaTools.cpp
//  TopHat
//
//  Created by Harold Pimentel on 10/27/11.
//

#include "FastaTools.h"

FastaReader::FastaReader()
{
    isPrimed_ = false;
}

FastaReader::FastaReader(std::string fname)
{
    isPrimed_ = false;
    init(fname);
}

FastaReader::~FastaReader()
{
    ifstream_.close();
}

void FastaReader::init(std::string fname)
{
    if (isPrimed_) {
        std::cerr << "Warning: object has already FastaReader has already been "
            << "initialized with file: " << fname_ << std::endl;
        return;
    }
    std::ios::sync_with_stdio(false); //to speed up slow iostream reading
    fname_ = fname;
    ifstream_.open(fname_.c_str(), std::ios::in);
    if (!ifstream_.good()) {
        std::cerr << "ERROR: Could not open file " << fname_ <<
            " in FastaReader" << std::endl;
        exit(1);
    }

    // Check the first character to see if it is valid
    char c = ifstream_.peek();
    if (c != '>')
    {
        std::cerr << "ERROR: Invalid format for FASTA file. Begins with a '" <<
            c << "'instead of a '>'" << std::endl;
        exit(1);
    }

    isPrimed_ = true;
}


bool FastaReader::good() const
{
    return ifstream_.good() && !ifstream_.eof();
}

//       Up to caller to allocate memory.
//       Only deallocates memory when there are no more records left
bool FastaReader::next(FastaRecord& rec) {
    if (!isPrimed_)
    {
        std::cerr << "ERROR: Stream has not been primed (FastaReader)"
            << std::endl;
        exit(1);
    }
    // Get the entire first line and description
    //ifstream_.getline(line_buf_, LINE_BUF_SIZE);
    if (ifstream_.eof() || !std::getline(ifstream_, line_buf_)) {
      rec.clear();
      return false;
    }

    if (line_buf_.empty() || !good()) {
        rec.clear();
        return false;
    }
    if (line_buf_.length()>0 && line_buf_[0]!='>') {
      std::cerr << "ERROR: no FASTA record start found (FastaReader)" << std::endl;
      exit(1);
    }
    size_t sp_pos = line_buf_.find(' ');
    if (sp_pos != std::string::npos) {
       rec.id_=line_buf_.substr(1, sp_pos-1);
       rec.desc_=line_buf_.substr(sp_pos+1);
    }
    else {
       rec.id_=line_buf_.substr(1);
       rec.desc_.clear();
    }
    rec.seq_.clear();

    // Read until you see another ">"
    while (ifstream_.peek() != '>') {
        //ifstream_ >> cur_line >> std::ws;
       if (std::getline(ifstream_, line_buf_))
           rec.seq_ += line_buf_;
       else {
           break; // if ifstream_.good() && !ifstream_.eof() &&
       }
    }
    return true;
}


FastaWriter::FastaWriter()
{
    isPrimed_ = false;
}

FastaWriter::FastaWriter(std::string fname)
{
    isPrimed_ = false;
    init(fname);
}

FastaWriter::~FastaWriter()
{
    ofstream_.close();
}

void FastaWriter::init(std::string fname)
{
    if (isPrimed_)
    {
        std::cerr << "Warning: Cannot allocate FastaWriter to file '" <<
            fname << "'. It has already been allocated to file '"
            << fname_ << "'" << std::endl;
        return;
    }

    ofstream_.open(fname.c_str(), std::ios::out);
    if (!ofstream_.good())
    {
        std::cerr << "ERROR: Could not open " << fname << " for writing in "
            << "FastaWriter" << std::endl;
        exit(1);
    }

    fname_ = fname;

    isPrimed_ = true;
}

void FastaWriter::write(FastaRecord& rec, size_t column_size)
{
    if (rec.seq_.length() == 0)
        return; //don't write empty records
    ofstream_ << ">" << rec.id_; //<< std::endl;
    if (rec.desc_.length()) {
      ofstream_ << " " << rec.desc_;
      }
    ofstream_ << std::endl;

    // iterate throught the string and print out the string
    size_t start = 0;
    while (start < rec.seq_.length())
    {
        ofstream_ << rec.seq_.substr(start, column_size) << std::endl;
        start += column_size;
    }
}
