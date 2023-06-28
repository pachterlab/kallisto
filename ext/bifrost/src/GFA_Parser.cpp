#include "GFA_Parser.hpp"


using namespace std;

GFA_Parser::GFA_Parser() : file_open_write(false), file_open_read(false), file_no(0), v_gfa(0), graph_out(nullptr), graph_in(nullptr) {}

GFA_Parser::GFA_Parser(const string& filename) :    file_open_write(false), file_open_read(false), file_no(0),
                                                    v_gfa(0), graph_out(nullptr), graph_in(nullptr) {

    graph_filenames.push_back(filename);
}

GFA_Parser::GFA_Parser(const vector<string>& filenames) :   file_open_write(false), file_open_read(false), file_no(0),
                                                            v_gfa(0), graph_out(nullptr), graph_in(nullptr) {

    graph_filenames = filenames;
}

GFA_Parser::GFA_Parser(GFA_Parser&& o) :    graph_filenames(move(o.graph_filenames)),
                                            graph_out(move(o.graph_out)), graph_in(move(o.graph_in)),
                                            v_gfa(o.v_gfa), file_no(o.file_no),
                                            file_open_write(o.file_open_write), file_open_read(o.file_open_read) {

    o.file_open_write = false;
    o.file_open_read = false;
}

GFA_Parser& GFA_Parser::operator=(GFA_Parser&& o){

    if (this != &o) {

        close();

        graph_filenames = o.graph_filenames;

        graph_in = move(o.graph_in);
        graph_out = move(o.graph_out);

        v_gfa = o.v_gfa;
        file_no = o.file_no;

        file_open_write = o.file_open_write;
        file_open_read = o.file_open_read;

        o.file_open_write = false;
        o.file_open_read = false;
    }

    return *this;
}

GFA_Parser::~GFA_Parser() {

    close();
}

bool GFA_Parser::open_write(const size_t version_GFA, const string tags_line_header, const bool compressed_output) {

    if (graph_filenames.size() == 0){

        cerr << "GFA_Parser::open_write(): No file specified in input" << endl;
        return false;
    }

    string fn = graph_filenames.front();

    {
        const size_t pos_ext = fn.find_last_of(".");

        if (pos_ext == string::npos) fn.append(compressed_output ? ".gfa.gz" : ".gfa");
        else if (!compressed_output && (fn.substr(pos_ext + 1) != "gfa")) fn.append(".gfa");
        else if (compressed_output && (fn.substr(pos_ext + 1) != "gz")) fn.append(".gfa.gz");
    }

    FILE* fp = fopen(fn.c_str(), "w");

    file_open_write = (fp != NULL);

    if (file_open_write) {

        fclose(fp);

        if (remove(fn.c_str()) != 0) cerr << "GFA_Parser::open_write(): Could not remove temporary file " << fn << endl;
    }
    else cerr << "GFA_Parser::open_write(): Could not open file " << fn << " for writing" << endl;

    if ((version_GFA != 1) && (version_GFA != 2)) {

        cerr << "GFA_Parser::open_write(): Only supports GFA format version 1 and 2" << endl;
        file_open_write = false;
    }
    else v_gfa = version_GFA;

    if (file_open_write){

        if (graph_out == nullptr) {

            if (compressed_output) graph_out = unique_ptr<ostream>(new zstr::ofstream(fn, ios_base::out));
            else graph_out = unique_ptr<ostream>(new std::ofstream(fn, ios_base::out));
        }

        {
            ostream& gout = *graph_out;

            gout << "H\tVN:Z:" << (v_gfa == 1 ? "1" : "2") << ".0";

            if (!tags_line_header.empty() && (tags_line_header != "")) gout << "\t" << tags_line_header;

            gout << "\n";
        }
    }

    return file_open_write;
}

pair<string, bool> GFA_Parser::open_read() {

    if (graph_filenames.size() == 0){

        cerr << "GFA_Parser::open_read(): No file specified in input" << endl;
        return {string(), false};
    }

    for (const auto& fn : graph_filenames){

        FILE* fp = fopen(fn.c_str(), "r");

        if (fp != NULL) fclose(fp);
        else {

            cerr << "GFA_Parser::open_read(): Could not open file " << fn << " for reading" << endl;
            return {string(), false};
        }
    }

    return open(file_no);
}

pair<string, bool> GFA_Parser::open(const size_t fn_id){

    string header;

    if (fn_id < graph_filenames.size()){

        const string& fn = graph_filenames[fn_id];

        FILE* fp = fopen(fn.c_str(), "r");

        file_open_read = (fp != NULL);

        if (file_open_read) fclose(fp);
        else cerr << "GFA_Parser::open(): Could not open file " << fn << " for reading" << endl;

        if (file_open_read) {

            if (graph_in == nullptr) graph_in = unique_ptr<istream>(new zstr::ifstream(fn, ios_base::in));

            getline(*graph_in, header);

            if (header.empty()) {

                cerr << "GFA_Parser::open(): Empty file: " << fn << endl;
                close();
            }
            else if (header[0] != 'H'){

                cerr << "GFA_Parser::open(): Wrong GFA header in " << fn << endl;
                close();
            }
            else if (header.substr(0, 10) == "H\tVN:Z:1.0") v_gfa = 1;
            else if (header.substr(0, 10) == "H\tVN:Z:2.0") v_gfa = 2;
            else {

                cerr << "GFA_Parser::open(): Unspecified GFA format version in " << fn <<
                ", version 1.0 is assumed by default." << endl;

                v_gfa = 1;
            }
        }
    }
    else file_open_read = false;

    return {header, file_open_read};
}

void GFA_Parser::close(){

    if (file_open_write){

        delete graph_out.release();

        graph_out = nullptr;
        file_open_write = false;
    }
    else if (file_open_read){

        delete graph_in.release();

        graph_in = nullptr;
        file_open_read = false;
    }
}

bool GFA_Parser::write_sequence(const string& id, const size_t len, const string seq, const string tags_line){

    if (file_open_write){

        ostream& gout = *graph_out;

        gout << "S" << "\t" << id;

        if (v_gfa == 2) gout << "\t" << len;

        gout << "\t" << seq;

        if (!tags_line.empty() && (tags_line != "")) gout << "\t" << tags_line;

        gout << "\n";
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

        ostream& gout = *graph_out;

        if (v_gfa == 1){

            if ((pos_end_overlapB - pos_start_overlapB) != (pos_end_overlapA - pos_start_overlapA)){

                cerr << "GFA_Parser::write_edge(): Overlap lengths must be the same for vertex A and B in GFA format version 1" << endl;
                close();
                return false;
            }

            gout << "L" << "\t" <<
            vertexA_id << "\t" << (strand_overlapA ? "+" : "-") << "\t" <<
            vertexB_id << "\t" << (strand_overlapB ? "+" : "-") << "\t" <<
            (pos_end_overlapA - pos_start_overlapA) << "M\n";
        }
        else {

            gout << "E" << "\t" << edge_id << "\t" <<
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

        istream& gin = *graph_in;

        while (getline(gin, line).good()){

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

        if (getline(gin, line).eof()){

            close();

            file_open_read = open(file_no + 1).second;

            if (file_open_read){

                ++file_no;

                file_id = file_no;

                return read(file_id);
            }
        }
        else if (getline(gin, line).fail()){

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

        istream& gin = *graph_in;

        while (getline(gin, line).good()){

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

        if (getline(gin, line).eof()){

            close();

            file_open_read = open(file_no + 1).second;

            if (file_open_read){

                ++file_no;

                file_id = file_no;
                new_file_opened = true;
            }
        }
        else if (getline(gin, line).fail()){

            cerr << "GFA_Parser::read(): Error while reading" << endl;
            close();
        }
    }
    else cerr << "GFA_Parser::read(): Input file is not open in reading mode" << endl;

    file_id = file_no;

    return make_pair(nullptr, nullptr);
}
