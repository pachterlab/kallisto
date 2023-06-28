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

#include "zstr.hpp"


class GFA_Parser {

    struct Sequence {

        std::string id;
        std::string seq;
        size_t len;

        std::vector<std::string> tags;

        Sequence() : seq("*"), len(0) {};
        Sequence(const std::string& id_, const std::string& seq_, const size_t len_) :  id(id_), seq(seq_), len(len_) {};

        inline void clear() {

            id.clear();
            tags.clear();

            seq = "*";
            len = 0;
        }
    };

    struct Edge {

        std::string edge_id;

        std::string vertexA_id;
        size_t pos_start_overlapA;
        size_t pos_end_overlapA;
        bool strand_overlapA;

        std::string vertexB_id;
        size_t pos_start_overlapB;
        size_t pos_end_overlapB;
        bool strand_overlapB;

        Edge() :    edge_id("*"),
                    vertexA_id(), pos_start_overlapA(0), pos_end_overlapA(0), strand_overlapA(true),
                    vertexB_id(), pos_start_overlapB(0), pos_end_overlapB(0), strand_overlapB(true) {};

        Edge(const std::string vertexA_id_, const size_t pos_start_overlapA_, const size_t pos_end_overlapA_, const bool strand_overlapA_,
             const std::string vertexB_id_, const size_t pos_start_overlapB_, const size_t pos_end_overlapB_, const bool strand_overlapB_,
             const std::string edge_id_ = "*") : edge_id(edge_id_),
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

        typedef std::pair<const Sequence*, const Edge*> GFA_line;

        GFA_Parser();
        GFA_Parser(const std::string& filename);
        GFA_Parser(const std::vector<std::string>& filenames);

        ~GFA_Parser();

        GFA_Parser(GFA_Parser&& o);
        GFA_Parser& operator=(GFA_Parser&& o);

        bool open_write(const size_t version_GFA = 1, const std::string tags_line_header = "", const bool compressed_output = false);
        std::pair<std::string, bool> open_read();

        void close();

        bool write_sequence(const std::string& id, const size_t len, const std::string seq = "*", const std::string tags_line = "");
        bool write_edge(const std::string vertexA_id, const size_t pos_start_overlapA, const size_t pos_end_overlapA, const bool strand_overlapA,
                        const std::string vertexB_id, const size_t pos_start_overlapB, const size_t pos_end_overlapB, const bool strand_overlapB,
                        const std::string edge_id = "*");

        GFA_line read(size_t& file_id);
        GFA_line read(size_t& file_id, bool& new_file_opened, const bool skip_edges = false);

    private:

        std::pair<std::string, bool> open(const size_t idx_filename);

        std::vector<std::string> graph_filenames;

        std::unique_ptr<std::istream> graph_in;
        std::unique_ptr<std::ostream> graph_out;

        size_t v_gfa;
        size_t file_no;

        char buffer_stream[8192];

        bool file_open_write;
        bool file_open_read;

        Sequence s;
        Edge e;
};

#endif
