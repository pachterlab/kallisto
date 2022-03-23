#include "FASTX_Parser.hpp"

FastqFile::FastqFile() : kseq(NULL), file_no(0) { fnit = fnames.end(); }

FastqFile::FastqFile(const vector<string> files) : kseq(NULL), fnames(files), file_no(0) {

    fnit = fnames.begin();
    fp = gzopen(fnit->c_str(), "r");
    kseq = kseq_init(fp);

    //std::ios::sync_with_stdio(false);
}

FastqFile::FastqFile(FastqFile&& o) : fp(o.fp), kseq(o.kseq), fnames(o.fnames), file_no(o.file_no) {

    fnit = fnames.begin();

    while (*fnit != *(o.fnit)) ++fnit;

    o.kseq = NULL;
}

FastqFile& FastqFile::operator=(FastqFile&& o){

    if (this != &o) {

        close();

        fp = o.fp;
        kseq = o.kseq;
        fnames = o.fnames;
        file_no = o.file_no;

        fnit = fnames.begin();

        while (*fnit != *(o.fnit)) ++fnit;

        o.kseq = NULL;
    }

    return *this;
}

FastqFile::~FastqFile() {

    close();
}

void FastqFile::reopen() {

    close();

    fnit = fnames.begin();
    fp = gzopen(fnit->c_str(), "r");
    kseq = kseq_init(fp);
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int FastqFile::read_next(char* read, size_t* read_len, string &seq, size_t* seq_len, unsigned int* file_id, char* qual) {

    int r;

    if ((r = kseq_read(kseq)) >= 0) {

        memcpy(read, kseq->name.s, kseq->name.l + 1); // 0-terminated string
        *read_len = kseq->name.l;
        seq.assign(kseq->seq.s);
        *seq_len = kseq->seq.l;

        if (qual != NULL) memcpy(qual, kseq->qual.s, kseq->qual.l + 1); // 0-terminated string
        if (file_id != NULL) *file_id = file_no / 2;
    }
    else if (r == -1) {

        open_next();

        if (fnit != fnames.end()) return read_next(read, read_len, seq, seq_len, file_id, qual);
        return -1;
    }

    return r;
}

// returns >=0 (length of seq), -1 end of last file, -2 truncated quality string
int FastqFile::read_next(string &seq, size_t& id) {

    const int r = kseq_read(kseq);

    if (r >= 0) seq.assign(kseq->seq.s);
    else if (r == -1) {

        open_next();

        if (fnit != fnames.end()){

            id = file_no;

            return read_next(seq, id);
        }
    }

    return r;
}

// next_file_opened == true indicates the next file to read in the list has been opened:
// first read entry of the file has not been read yet, method returns 0.
// Else, method returns >= 0 (length of seq), -1 end of last file, -2 truncated quality string.
int FastqFile::read_next(string &seq, size_t& id, bool& next_file_opened) {

    const int r = kseq_read(kseq);

    next_file_opened = false;

    if (r >= 0) seq.assign(kseq->seq.s);
    else if (r == -1) {

        open_next();

        if (fnit != fnames.end()){

            id = file_no;
            next_file_opened = true;

            return 0;
        }
    }

    return r;
}

int FastqFile::read_next(stringstream& ss, size_t& id, bool& next_file_opened) {

    const int r = kseq_read(kseq);

    next_file_opened = false;

    if (r >= 0) ss << kseq->seq.s;
    else if (r == -1) {

        open_next();

        if (fnit != fnames.end()){

            id = file_no;
            next_file_opened = true;

            return 0;
        }
    }

    return r;
}

int FastqFile::read_next() {

    const int r = kseq_read(kseq);

    if (r == -1) {

        open_next();

        if (fnit != fnames.end()) return read_next();
    }

    return r;
}

vector<string>::const_iterator FastqFile::open_next() {

    if (fnit != fnames.end()) {
        // close current file
        kseq_destroy(kseq);
        gzclose(fp);

        kseq = NULL;

        // get next file
        ++fnit;
        ++file_no;

        if (fnit != fnames.end()) {

            fp = gzopen(fnit->c_str(), "r");
            kseq = kseq_init(fp);
        }
    }

    return fnit;
}

void FastqFile::close() {

    if (kseq != NULL) {

        kseq_destroy(kseq);
        gzclose(fp);

        fnit = fnames.end();
        kseq = NULL;
    }
}
