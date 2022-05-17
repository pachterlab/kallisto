#include "GFA_Parser.hpp"

GFA_Parser::GFA_Parser() : file_open_write(false), file_open_read(false), file_no(0), v_gfa(0), graph_out(nullptr), graph_in(nullptr),
                           graphfile_in(nullptr), graphfile_out(nullptr) {}

GFA_Parser::GFA_Parser(const string& filename) :    file_open_write(false), file_open_read(false), file_no(0),
                                                    v_gfa(0), graph_out(nullptr), graph_in(nullptr),
                                                    graphfile_in(nullptr), graphfile_out(nullptr) {

    graph_filenames.push_back(filename);

    size_t pos_match_point = graph_filenames[0].find_last_of(".");

    if ((pos_match_point == string::npos) || (graph_filenames[0].substr(pos_match_point + 1) != "gfa")) graph_filenames[0].append(".gfa");
}

GFA_Parser::GFA_Parser(const vector<string>& filenames) :   file_open_write(false), file_open_read(false), file_no(0),
                                                            v_gfa(0), graph_out(nullptr), graph_in(nullptr),
                                                            graphfile_in(nullptr), graphfile_out(nullptr) {

    graph_filenames = filenames;

    for (auto& filename : graph_filenames){

        size_t pos_match_point = filename.find_last_of(".");

        if ((pos_match_point == string::npos) || (filename.substr(pos_match_point + 1) != "gfa")) filename.append(".gfa");
    }
}

GFA_Parser::GFA_Parser(GFA_Parser&& o) :    graph_filenames(o.graph_filenames), graphfile_out(move(o.graphfile_out)),
                                            graphfile_in(move(o.graphfile_in)), graph_out(nullptr), graph_in(nullptr),
                                            v_gfa(o.v_gfa), file_no(o.file_no),
                                            file_open_write(o.file_open_write), file_open_read(o.file_open_read) {

    if (file_open_write) graph_out.rdbuf(graphfile_out->rdbuf());
    if (file_open_read) graph_in.rdbuf(graphfile_in->rdbuf());

    o.file_open_write = false;
    o.file_open_read = false;
}

GFA_Parser& GFA_Parser::operator=(GFA_Parser&& o){

    if (this != &o) {

        close();

        graph_filenames = o.graph_filenames;

        graphfile_in = move(o.graphfile_in);
        graphfile_out = move(o.graphfile_out);

        v_gfa = o.v_gfa;
        file_no = o.file_no;

        file_open_write = o.file_open_write;
        file_open_read = o.file_open_read;

        if (file_open_write) graph_out.rdbuf(graphfile_out->rdbuf());
        if (file_open_read) graph_in.rdbuf(graphfile_in->rdbuf());

        o.file_open_write = false;
        o.file_open_read = false;
    }

    return *this;
}

GFA_Parser::~GFA_Parser() {

    close();
}

bool GFA_Parser::open_write(const size_t version_GFA, const string tags_line_header) {

    if (graph_filenames.size() == 0){

        cerr << "GFA_Parser::open_write(): No file specified in input" << endl;
        return false;
    }

    const string& filename = graph_filenames[0];

    FILE* fp = fopen(filename.c_str(), "w");

    if ((file_open_write = (fp != NULL)) == true) {

        fclose(fp);

        if (remove(filename.c_str()) != 0) cerr << "GFA_Parser::open_write(): Could not remove temporary file " << filename << endl;
    }
    else cerr << "GFA_Parser::open_write(): Could not open file " << filename << " for writing" << endl;

    if ((version_GFA != 1) && (version_GFA != 2)) {

        cerr << "GFA_Parser::open_write(): Only supports GFA format version 1 and 2" << endl;
        file_open_write = false;
    }
    else v_gfa = version_GFA;

    if (file_open_write){

        if (graphfile_out == nullptr) graphfile_out = new ofstream();

        graphfile_out->open(filename.c_str(), ios_base::out);
        graphfile_out->rdbuf()->pubsetbuf(buffer_stream, sizeof(buffer_stream));

        graph_out.rdbuf(graphfile_out->rdbuf());
        //graph_out.sync_with_stdio(false);

        graph_out << "H\tVN:Z:" << (v_gfa == 1 ? "1" : "2") << ".0";

        if (!tags_line_header.empty() && (tags_line_header != "")) graph_out << "\t" << tags_line_header;

        graph_out << "\n";
    }

    return file_open_write;
}

bool GFA_Parser::open_read() {

    if (graph_filenames.size() == 0){

        cerr << "GFA_Parser::open_read(): No file specified in input" << endl;
        return false;
    }

    for (const auto& filename : graph_filenames){

        FILE* fp = fopen(filename.c_str(), "r");

        if (fp != NULL) fclose(fp);
        else cerr << "GFA_Parser::open_read(): Could not open file " << filename << " for reading" << endl;
    }

    file_open_read = open(file_no);

    return file_open_read;
}

bool GFA_Parser::open(const size_t idx_filename){

    if (idx_filename < graph_filenames.size()){

        FILE* fp = fopen(graph_filenames[idx_filename].c_str(), "r");

        if ((file_open_read = (fp != NULL)) == true) fclose(fp);
        else cerr << "GFA_Parser::open(): Could not open file " << graph_filenames[idx_filename] << " for reading" << endl;

        if (file_open_read) {

            if (graphfile_in == nullptr) graphfile_in = new ifstream();

            graphfile_in->open(graph_filenames[idx_filename], ios_base::in);
            graphfile_in->rdbuf()->pubsetbuf(buffer_stream, sizeof(buffer_stream));

            graph_in.rdbuf(graphfile_in->rdbuf());
            //graph_in.sync_with_stdio(false);

            string header;

            getline(graph_in, header);

            if (header.empty()) {

                cerr << "GFA_Parser::open(): Empty file: " << graph_filenames[idx_filename] << endl;
                close();
            }
            else if (header[0] != 'H'){

                cerr << "GFA_Parser::open(): Wrong GFA header in " << graph_filenames[idx_filename] << endl;
                close();
            }
            else if (header.substr(0, 10) == "H\tVN:Z:1.0") v_gfa = 1;
            else if (header.substr(0, 10) == "H\tVN:Z:2.0") v_gfa = 2;
            else {

                cerr << "GFA_Parser::open(): Unspecified GFA format version in " << graph_filenames[idx_filename] <<
                ", version 1.0 is assumed by default." << endl;

                v_gfa = 1;
            }
        }
    }
    else file_open_read = false;

    return file_open_read;
}

void GFA_Parser::close(){

    if (file_open_write){

        graphfile_out->close();

        delete graphfile_out;

        graphfile_out = nullptr;
        file_open_write = false;
    }
    else if (file_open_read){

        graphfile_in->close();

        delete graphfile_in;

        graphfile_in = nullptr;
        file_open_read = false;
    }
}

bool GFA_Parser::write_sequence(const string& id, const size_t len, const string seq, const string tags_line){

    if (file_open_write){

        graph_out << "S" << "\t" << id;

        if (v_gfa == 2) graph_out << "\t" << len;

        graph_out << "\t" << seq;

        if (!tags_line.empty() && (tags_line != "")) graph_out << "\t" << tags_line;

        graph_out << "\n";
    }
    else cerr << "GFA_Parser::write_sequence(): Input file is not open in writing mode" << endl;

    return file_open_write;
}

bool GFA_Parser::write_edge(const string vertexA_id, const size_t pos_start_overlapA, const size_t pos_end_overlapA, const bool strand_overlapA,
                            const string vertexB_id, const size_t pos_start_overlapB, const size_t pos_end_overlapB, const bool strand_overlapB,
                            const string edge_id) {

    if (file_open_write){

        if (pos_start_overlapA > pos_end_overlapA){

            cerr << "GFA_Parser::write_edge(): Vertex A overlap start position greater than vertex A overlap end position" << endl;
            close();
            return false;
        }

        if (pos_start_overlapB > pos_end_overlapB){

            cerr << "GFA_Parser::write_edge(): Vertex B overlap start position greater than vertex B overlap end position" << endl;
            close();
            return false;
        }

        if (v_gfa == 1){

            if ((pos_end_overlapB - pos_start_overlapB) != (pos_end_overlapA - pos_start_overlapA)){

                cerr << "GFA_Parser::write_edge(): Overlap lengths must be the same for vertex A and B in GFA format version 1" << endl;
                close();
                return false;
            }

            graph_out << "L" << "\t" <<
            vertexA_id << "\t" << (strand_overlapA ? "+" : "-") << "\t" <<
            vertexB_id << "\t" << (strand_overlapB ? "+" : "-") << "\t" <<
            (pos_end_overlapA - pos_start_overlapA) << "M\n";
        }
        else {

            graph_out << "E" << "\t" << edge_id << "\t" <<
            vertexA_id << (strand_overlapA ? "+" : "-") << "\t" <<
            vertexB_id << (strand_overlapB ? "+" : "-") << "\t" <<
            pos_start_overlapA << "\t" << pos_end_overlapA << "\t" <<
            pos_start_overlapB << "\t" << pos_end_overlapB << "\t" <<
            "*" << "\n";
        }
    }
    else {

        cerr << "GFA_Parser::write_edge(): Input file is not open in writing mode" << endl;
        return false;
    }

    return true;
}

GFA_Parser::GFA_line GFA_Parser::read(size_t& file_id) {

    if (file_open_read){

        vector<string> line_fields;

        string line;

        while (getline(graph_in, line).good()){

            if (line[0] == 'S'){ // Segment line

                const char* buffer = line.c_str() + 2;
                const char* end_buffer = line.c_str() + line.length();
                const char* prev_buffer = strchr(buffer, '\t');

                while (prev_buffer != NULL){

                    line_fields.push_back(string(buffer, prev_buffer - buffer));
                    buffer = prev_buffer + 1;
                    prev_buffer = strchr(buffer, '\t');
                }

                if (end_buffer - buffer != 0) line_fields.push_back(string(buffer, end_buffer - buffer));

                const size_t line_fields_sz = line_fields.size();

                s.clear();

                if (v_gfa == 1){ // GFA format version 1

                    if (line_fields_sz < 2){

                        cerr << "GFA_Parser::read(): Missing fields in Segment line" << endl;
                        close();
                    }

                    s.id = move(line_fields[0]);
                    s.seq = move(line_fields[1]);

                    for (size_t i = 2; i < line_fields_sz; ++i) s.tags.push_back(move(line_fields[i]));
                }
                else {

                    if (line_fields_sz < 3){

                        cerr << "GFA_Parser::read(): Missing fields in Segment line" << endl;
                        close();
                    }

                    s.id = move(line_fields[0]);
                    s.len = sscanf(line_fields[1].c_str(), "%zu", &(s.len));
                    s.seq = move(line_fields[2]);

                    for (size_t i = 3; i < line_fields_sz; ++i) s.tags.push_back(move(line_fields[i]));
                }

                file_id = file_no;

                return make_pair(&s, nullptr);
            }
            else if ((v_gfa == 1) && (line[0] == 'L')){ // Link line, only GFA v1

                const char* buffer = line.c_str() + 2;
                const char* end_buffer = line.c_str() + line.length();
                const char* prev_buffer = strchr(buffer, '\t');

                while (prev_buffer != NULL){

                    line_fields.push_back(string(buffer, prev_buffer - buffer));
                    buffer = prev_buffer + 1;
                    prev_buffer = strchr(buffer, '\t');
                }

                if (end_buffer - buffer != 0) line_fields.push_back(string(buffer, end_buffer - buffer));

                const size_t line_fields_sz = line_fields.size();

                e.clear();

                if (line_fields_sz < 4){

                    cerr << "GFA_Parser::read(): Missing fields in Link line" << endl;
                    close();
                }

                e.vertexA_id = move(line_fields[0]);
                e.vertexB_id = move(line_fields[2]);

                if (line_fields[1] == "+") e.strand_overlapA = true;
                else if (line_fields[1] == "-") e.strand_overlapA = false;
                else {

                    cerr << "GFA_Parser::read(): Orientation of Segment A on Link line is not + or -" << endl;
                    close();
                }

                if (line_fields[3] == "+") e.strand_overlapB = true;
                else if (line_fields[3] == "-") e.strand_overlapB = false;
                else {

                    cerr << "GFA_Parser::read(): Orientation of Segment B on Link line is not + or -" << endl;
                    close();
                }

                file_id = file_no;

                return make_pair(nullptr, &e);
            }
            else if ((v_gfa == 2) && (line[0] == 'E')){ // Edge line, only GFA v2

                const char* buffer = line.c_str() + 2;
                const char* end_buffer = line.c_str() + line.length();
                const char* prev_buffer = strchr(buffer, '\t');

                while (prev_buffer != NULL){

                    line_fields.push_back(string(buffer, prev_buffer - buffer));
                    buffer = prev_buffer + 1;
                    prev_buffer = strchr(buffer, '\t');
                }

                if (end_buffer - buffer != 0) line_fields.push_back(string(buffer, end_buffer - buffer));

                const size_t line_fields_sz = line_fields.size();

                e.clear();

                if (line_fields_sz < 8){

                    cerr << "GFA_Parser::read(): Missing fields in Edge line" << endl;
                    close();
                }

                e.edge_id = line_fields[0];

                const char ca = line_fields[1][line_fields[1].length() - 1]; // Last char. of line_fields[1];

                e.strand_overlapA = (ca != '-');

                if ((ca == '-') || (ca == '-')) e.vertexA_id = line_fields[1].substr(0, line_fields[1].length() - 1);
                else e.vertexA_id = move(line_fields[1]);

                sscanf(line_fields[2].c_str(), "%zu", &(e.pos_start_overlapA));
                sscanf(line_fields[3].c_str(), "%zu", &(e.pos_end_overlapA));

                const char cb = line_fields[4][line_fields[4].length() - 1]; // Last char. of line_fields[4];

                e.strand_overlapB = (cb != '-');

                if ((cb == '-') || (cb == '-')) e.vertexB_id = line_fields[4].substr(0, line_fields[4].length() - 1);
                else e.vertexB_id = move(line_fields[4]);

                sscanf(line_fields[5].c_str(), "%zu", &(e.pos_start_overlapB));
                sscanf(line_fields[6].c_str(), "%zu", &(e.pos_end_overlapB));

                file_id = file_no;

                return make_pair(nullptr, &e);
            }

            line_fields.clear();
        }

        if (getline(graph_in, line).eof()){

            close();

            if ((file_open_read = open(file_no + 1))){

                ++file_no;

                file_id = file_no;

                return read(file_id);
            }
        }
        else if (getline(graph_in, line).fail()){

            cerr << "GFA_Parser::read(): Error while reading" << endl;
            close();
        }
    }
    else cerr << "GFA_Parser::read(): Input file is not open in reading mode" << endl;

    file_id = file_no;

    return make_pair(nullptr, nullptr);
}

GFA_Parser::GFA_line GFA_Parser::read(size_t& file_id, bool& new_file_opened, const bool skip_edges) {

    new_file_opened = false;

    if (file_open_read){

        vector<string> line_fields;

        string line;

        while (getline(graph_in, line).good()){

            if (line[0] == 'S'){ // Segment line

                const char* buffer = line.c_str() + 2;
                const char* end_buffer = line.c_str() + line.length();
                const char* prev_buffer = strchr(buffer, '\t');

                while (prev_buffer != NULL){

                    line_fields.push_back(string(buffer, prev_buffer - buffer));
                    buffer = prev_buffer + 1;
                    prev_buffer = strchr(buffer, '\t');
                }

                if (end_buffer - buffer != 0) line_fields.push_back(string(buffer, end_buffer - buffer));

                const size_t line_fields_sz = line_fields.size();

                s.clear();

                if (v_gfa == 1){ // GFA format version 1

                    if (line_fields_sz < 2){

                        cerr << "GFA_Parser::read(): Missing fields in Segment line: " << endl;
                        cerr << line << endl;
                        close();
                    }

                    s.id = move(line_fields[0]);
                    s.seq = move(line_fields[1]);

                    for (size_t i = 2; i < line_fields_sz; ++i) s.tags.push_back(move(line_fields[i]));
                }
                else {

                    if (line_fields_sz < 3){

                        cerr << "GFA_Parser::read(): Missing fields in Segment line" << endl;
                        close();
                    }

                    s.id = move(line_fields[0]);
                    s.len = sscanf(line_fields[1].c_str(), "%zu", &(s.len));
                    s.seq = move(line_fields[2]);

                    for (size_t i = 3; i < line_fields_sz; ++i) s.tags.push_back(move(line_fields[i]));
                }

                file_id = file_no;

                return make_pair(&s, nullptr);
            }
            else if (!skip_edges){

                if ((v_gfa == 1) && (line[0] == 'L')){ // Link line, only GFA v1

                    const char* buffer = line.c_str() + 2;
                    const char* end_buffer = line.c_str() + line.length();
                    const char* prev_buffer = strchr(buffer, '\t');

                    while (prev_buffer != NULL){

                        line_fields.push_back(string(prev_buffer, prev_buffer - prev_buffer));
                        prev_buffer = prev_buffer + 1;
                        prev_buffer = strchr(buffer, '\t');
                    }

                    if (end_buffer - prev_buffer != 0) line_fields.push_back(string(prev_buffer, end_buffer - prev_buffer));

                    const size_t line_fields_sz = line_fields.size();

                    e.clear();

                    if (line_fields_sz < 4){

                        cerr << "GFA_Parser::read(): Missing fields in Link line" << endl;
                        close();
                    }

                    e.vertexA_id = move(line_fields[0]);
                    e.vertexB_id = move(line_fields[2]);

                    if (line_fields[1] == "+") e.strand_overlapA = true;
                    else if (line_fields[1] == "-") e.strand_overlapA = false;
                    else {

                        cerr << "GFA_Parser::read(): Orientation of Segment A on Link line is not + or -" << endl;
                        close();
                    }

                    if (line_fields[3] == "+") e.strand_overlapB = true;
                    else if (line_fields[3] == "-") e.strand_overlapB = false;
                    else {

                        cerr << "GFA_Parser::read(): Orientation of Segment B on Link line is not + or -" << endl;
                        close();
                    }

                    file_id = file_no;

                    return make_pair(nullptr, &e);
                }
                else if ((v_gfa == 2) && (line[0] == 'E')){ // Edge line, only GFA v2

                    const char* buffer = line.c_str() + 2;
                    const char* end_buffer = line.c_str() + line.length();
                    const char* prev_buffer = strchr(buffer, '\t');

                    while (prev_buffer != NULL){

                        line_fields.push_back(string(buffer, prev_buffer - buffer));
                        buffer = prev_buffer + 1;
                        prev_buffer = strchr(buffer, '\t');
                    }

                    if (end_buffer - buffer != 0) line_fields.push_back(string(buffer, end_buffer - buffer));

                    const size_t line_fields_sz = line_fields.size();

                    e.clear();

                    if (line_fields_sz < 8){

                        cerr << "GFA_Parser::read(): Missing fields in Edge line" << endl;
                        close();
                    }

                    e.edge_id = move(line_fields[0]);

                    const char ca = line_fields[1][line_fields[1].length() - 1]; // Last char. of line_fields[1];

                    e.strand_overlapA = (ca != '-');

                    if ((ca == '-') || (ca == '-')) e.vertexA_id = line_fields[1].substr(0, line_fields[1].length() - 1);
                    else e.vertexA_id = move(line_fields[1]);

                    sscanf(line_fields[2].c_str(), "%zu", &(e.pos_start_overlapA));
                    sscanf(line_fields[3].c_str(), "%zu", &(e.pos_end_overlapA));

                    const char cb = line_fields[4][line_fields[4].length() - 1]; // Last char. of line_fields[4];

                    e.strand_overlapB = (cb != '-');

                    if ((cb == '-') || (cb == '-')) e.vertexB_id = line_fields[4].substr(0, line_fields[4].length() - 1);
                    else e.vertexB_id = move(line_fields[4]);

                    sscanf(line_fields[5].c_str(), "%zu", &(e.pos_start_overlapB));
                    sscanf(line_fields[6].c_str(), "%zu", &(e.pos_end_overlapB));

                    file_id = file_no;

                    return make_pair(nullptr, &e);
                }
            }

            line_fields.clear();
        }

        if (getline(graph_in, line).eof()){

            close();

            if ((file_open_read = open(file_no + 1))){

                ++file_no;

                file_id = file_no;
                new_file_opened = true;
            }
        }
        else if (getline(graph_in, line).fail()){

            cerr << "GFA_Parser::read(): Error while reading" << endl;
            close();
        }
    }
    else cerr << "GFA_Parser::read(): Input file is not open in reading mode" << endl;

    file_id = file_no;

    return make_pair(nullptr, nullptr);
}
