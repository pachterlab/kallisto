#ifndef BIFROST_GFA_PARSER_HPP
#define BIFROST_GFA_PARSER_HPP

#include <string>
#include <cstring>
#include <vector>
#include <sys/stat.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

using namespace std;

class GFA_Parser {

    struct Sequence {

        string id;
        string seq;
        size_t len;

        vector<string> tags;

        Sequence() : seq("*"), len(0) {};
        Sequence(const string& id_, const string& seq_, const size_t len_) :  id(id_), seq(seq_), len(len_) {};

        inline void clear() {

            id.clear();
            tags.clear();

            seq = "*";
            len = 0;
        }
    };

    struct Edge {

        string edge_id;

        string vertexA_id;
        size_t pos_start_overlapA;
        size_t pos_end_overlapA;
        bool strand_overlapA;

        string vertexB_id;
        size_t pos_start_overlapB;
        size_t pos_end_overlapB;
        bool strand_overlapB;

        Edge() :    edge_id("*"),
                    vertexA_id(), pos_start_overlapA(0), pos_end_overlapA(0), strand_overlapA(true),
                    vertexB_id(), pos_start_overlapB(0), pos_end_overlapB(0), strand_overlapB(true) {};

        Edge(const string vertexA_id_, const size_t pos_start_overlapA_, const size_t pos_end_overlapA_, const bool strand_overlapA_,
             const string vertexB_id_, const size_t pos_start_overlapB_, const size_t pos_end_overlapB_, const bool strand_overlapB_,
             const string edge_id_ = "*") : edge_id(edge_id_),
             vertexA_id(vertexA_id_), pos_start_overlapA(pos_start_overlapA_), pos_end_overlapA(pos_end_overlapA_), strand_overlapA(strand_overlapA_),
             vertexB_id(vertexB_id_), pos_start_overlapB(pos_start_overlapB_), pos_end_overlapB(pos_end_overlapB_), strand_overlapB(strand_overlapB_) {};

        inline void clear() {

            vertexA_id.clear();
            vertexB_id.clear();

            edge_id = "*";

            pos_start_overlapA = 0;
            pos_end_overlapA = 0;
            pos_start_overlapB = 0;
            pos_end_overlapB = 0;

            strand_overlapA = true;
            strand_overlapB = true;
        }
    };

    public:

        typedef pair<const Sequence*, const Edge*> GFA_line;

        GFA_Parser();
        GFA_Parser(const string& filename);
        GFA_Parser(const vector<string>& filenames);

        ~GFA_Parser();

        GFA_Parser(GFA_Parser&& o);
        GFA_Parser& operator=(GFA_Parser&& o);

        bool open_write(const size_t version_GFA = 1, const string tags_line_header = "");
        bool open_read();

        void close();

        bool write_sequence(const string& id, const size_t len, const string seq = "*", const string tags_line = "");
        bool write_edge(const string vertexA_id, const size_t pos_start_overlapA, const size_t pos_end_overlapA, const bool strand_overlapA,
                        const string vertexB_id, const size_t pos_start_overlapB, const size_t pos_end_overlapB, const bool strand_overlapB,
                        const string edge_id = "*");

        GFA_line read(size_t& file_id);
        GFA_line read(size_t& file_id, bool& new_file_opened, const bool skip_edges = false);

    private:

        bool open(const size_t idx_filename);

        vector<string> graph_filenames;

        ifstream* graphfile_in;
        istream graph_in;

        ofstream* graphfile_out;
        ostream graph_out;

        size_t v_gfa;
        size_t file_no;

        char buffer_stream[8192];

        bool file_open_write;
        bool file_open_read;

        Sequence s;
        Edge e;
};

#endif
