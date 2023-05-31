#include "CompactedDBG.hpp"
#include "ColoredCDBG.hpp"

using namespace std;

void PrintVersion() {

    cout << BFG_VERSION << endl;
}

void PrintUsage() {

    cout << "Bifrost " << BFG_VERSION << endl << endl;

    cout << "Highly parallel construction, indexing and querying of colored and compacted de Bruijn graphs" << endl << endl;

    cout << "Usage: Bifrost [COMMAND] [PARAMETERS]" << endl << endl;

    cout << "[COMMAND]:" << endl << endl;

    cout << "   build                   Build a compacted de Bruijn graph, with or without colors" << endl;
    cout << "   update                  Update a compacted (colored) de Bruijn graph with new sequences" << endl;
    cout << "   query                   Query a compacted (colored) de Bruijn graph" << endl << endl;

    cout << "[PARAMETERS]: build" << endl << endl;

    cout << "   > Mandatory with required argument:" << endl << endl;

    cout << "   -s, --input-seq-file     Input sequence file in fasta/fastq(.gz) format" << endl;
    cout << "                            Multiple files can be provided as a list in a text file (one file per line)" << endl;
    cout << "                            K-mers with exactly 1 occurrence in the input sequence files will be discarded" << endl;
    cout << "   -r, --input-ref-file     Input reference file in fasta/fastq(.gz) or gfa(.gz) format" << endl;
    cout << "                            Multiple files can be provided as a list in a text file (one file per line)" << endl;
    cout << "                            All k-mers of the input reference files are used" << endl;
    cout << "   -o, --output-file        Prefix for output file(s)" << endl << endl;

    cout << "   > Optional with required argument:" << endl << endl;

    cout << "   -t, --threads            Number of threads (default is 1)" << endl;
    cout << "   -k, --kmer-length        Length of k-mers (default is 31)" << endl;
    cout << "   -m, --min-length         Length of minimizers (default is automatically chosen)" << endl;
    cout << "   -B, --bloom-bits         Number of Bloom filter bits per k-mer (default is 14)" << endl;
    cout << "   -l, --load-mbbf          Input Blocked Bloom Filter file, skips filtering step (default is no input)" << endl;
    cout << "   -w, --write-mbbf         Output Blocked Bloom Filter file (default is no output)" << endl;

    cout << "   > Optional with no argument:" << endl << endl;

    cout << "   -c, --colors             Color the compacted de Bruijn graph (default is no coloring)" << endl;
    //cout << "   -y, --keep-mercy         Keep low coverage k-mers connecting tips" << endl;
    cout << "   -i, --clip-tips          Clip tips shorter than k k-mers in length" << endl;
    cout << "   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length" << endl;
    cout << "   -f, --fasta-out          Output file is in fasta format (only sequences) instead of gfa (unless graph is colored)" << endl;
    cout << "   -b, --bfg-out            Output file is in bfg/bfi format (Bifrost graph and index) instead of gfa (unless graph is colored)" << endl;
    cout << "   -n, --no-compress-out    Output files must be uncompressed" << endl << endl;
    cout << "   -N, --no-index-out       No index file is created" << endl << endl;
    cout << "   -v, --verbose            Print information messages during execution" << endl << endl;

    cout << "[PARAMETERS]: update" << endl << endl;

    cout << "  > Mandatory with required argument:" << endl << endl;

    cout << "   -g, --input-graph-file   Input graph file to update in gfa(.gz) or bfg format" << endl;
    cout << "   -s, --input-seq-file     Input sequence file in fasta/fastq(.gz) format" << endl;
    cout << "                            Multiple files can be provided as a list in a text file (one file per line)" << endl;
    cout << "                            K-mers with exactly 1 occurrence in the input sequence files will be discarded" << endl;
    cout << "   -r, --input-ref-file     Input reference file in fasta/fastq(.gz) or gfa(.gz) format" << endl;
    cout << "                            Multiple files can be provided as a list in a text file (one file per line)" << endl;
    cout << "                            All k-mers of the input reference files are used" << endl;
    cout << "   -o, --output-file        Prefix for output file(s)" << endl << endl;

    cout << "   > Optional with required argument:" << endl << endl;

    cout << "   -I, --input-index-file   Input index file associated with graph to update in bfi format" << endl;
    cout << "   -C, --input-color-file   Input color file associated with graph to update in color.bfg format" << endl;
    cout << "   -t, --threads            Number of threads (default is 1)" << endl;
    cout << "   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)" << endl;
    cout << "   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or automatically chosen)" << endl << endl;

    cout << "   > Optional with no argument:" << endl << endl;

    cout << "   -i, --clip-tips          Clip tips shorter than k k-mers in length" << endl;
    cout << "   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length" << endl;
    cout << "   -f, --fasta-out          Output file is in fasta format (only sequences) instead of gfa (unless colors are output)" << endl;
    cout << "   -b, --bfg-out            Output file is in bfg/bfi format (Bifrost graph and index) instead of gfa (unless graph is colored)" << endl;
    cout << "   -n, --no-compress-out    Output files must be uncompressed" << endl << endl;
    cout << "   -N, --no-index-out       No index file is created" << endl << endl;
    cout << "   -v, --verbose            Print information messages during execution" << endl << endl;

    cout << "[PARAMETERS]: query" << endl << endl;

    cout << "  > Mandatory with required argument:" << endl << endl;

    cout << "   -g, --input-graph-file   Input graph file to query in gfa(.gz) or bfg" << endl;
    cout << "   -q, --input-query-file   Input query file in fasta/fastq(.gz)" << endl;
    cout << "                            Multiple files can be provided as a list in a text file (one file per line)" << endl;
    cout << "   -o, --output-file        Prefix for output file" << endl;
    cout << "   -e, --ratio-kmers        Ratio of k-mers from queries that must occur in the graph (default is 0.8)" << endl << endl;

    cout << "   > Optional with required argument:" << endl << endl;

    cout << "   -I, --input-index-file   Input index file associated with graph to query in bfi format" << endl;
    cout << "   -C, --input-color-file   Input color file associated with the graph to query in color.bfg format" << endl;
    cout << "                            Presence/absence of queries will be output for each color" << endl;
    cout << "   -t, --threads            Number of threads (default is 1)" << endl;
    cout << "   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)" << endl;
    cout << "   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or or automatically chosen)" << endl << endl;

    cout << "   > Optional with no argument:" << endl << endl;

    cout << "   -a, --approximate        Graph is searched with exact and inexact k-mers (1 substitution or indel) from queries" << endl;
    cout << "   -v, --verbose            Print information messages during execution" << endl << endl;
}

int parse_ProgramOptions(int argc, char **argv, CCDBG_Build_opt& opt) {

    int option_index = 0, c;

    const char* opt_string = "s:r:q:g:I:C:o:t:k:m:e:B:l:w:aidvcyfbnN";

    static struct option long_options[] = {

        {"input-seq-file",      required_argument,  0, 's'},
        {"input-ref-file",      required_argument,  0, 'r'},
        {"input-query-file",    required_argument,  0, 'q'},
        {"input-graph-file",    required_argument,  0, 'g'},
        {"input-index-file",    required_argument,  0, 'I'},
        {"input-color-file",    required_argument,  0, 'C'},
        {"output-file",         required_argument,  0, 'o'},
        {"threads",             required_argument,  0, 't'},
        {"kmer-length",         required_argument,  0, 'k'},
        {"min-length",          required_argument,  0, 'm'},
        {"ratio-kmers",         required_argument,  0, 'e'},
        {"bloom-bits",          required_argument,  0, 'B'},
        {"load-mbbf",           required_argument,  0, 'l'},
        {"write-mbbf",          required_argument,  0, 'w'},
        {"approximate",         no_argument,        0, 'a'},
        {"clip-tips",           no_argument,        0, 'i'},
        {"del-isolated",        no_argument,        0, 'd'},
        {"verbose",             no_argument,        0, 'v'},
        {"colors",              no_argument,        0, 'c'},
        //{"keep-mercy",          no_argument,        0, 'y'},
        {"fasta-out",           no_argument,        0, 'f'},
        {"bfg-out",             no_argument,        0, 'b'},
        {"no-compress-out",     no_argument,        0, 'n'},
        {"no-index-out",        no_argument,        0, 'N'},
        {0,                     0,                  0,  0 }
    };

    if (strcmp(argv[1], "--version") == 0) return 1; // Print version
    else if (strcmp(argv[1], "--help") == 0) return 2; // print help

    if (strcmp(argv[1], "build") == 0) opt.build = true;
    else if (strcmp(argv[1], "update") == 0) opt.update = true;
    else if (strcmp(argv[1], "query") == 0) opt.query = true;

    if (opt.build || opt.update || opt.query){

        while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

            switch (c) {

                case 's':
                    opt.filename_seq_in.push_back(optarg);
                    break;
                case 'r':
                    opt.filename_ref_in.push_back(optarg);
                    break;
                case 'q':
                    opt.filename_query_in.push_back(optarg);
                    break;
                case 'g':
                    opt.filename_graph_in = optarg;
                    break;
                case 'I':
                    opt.filename_index_in = optarg;
                    break;
                case 'C':
                    opt.filename_colors_in = optarg;
                    break;
                case 'o':
                    opt.prefixFilenameOut = optarg;
                    break;
                case 't':
                    opt.nb_threads = atoi(optarg);
                    break;
                case 'k':
                    opt.k = atoi(optarg);
                    break;
                case 'm':
                    opt.g = atoi(optarg);
                    break;
                case 'e':
                    opt.ratio_kmers = atof(optarg);
                    break;
                case 'B':
                    opt.nb_bits_kmers_bf = atoi(optarg);
                    break;
                case 'w':
                    opt.outFilenameBBF = optarg;
                    break;
                case 'l':
                    opt.inFilenameBBF = optarg;
                    break;
                case 'a':
                    opt.inexact_search = true;
                    break;
                case 'i':
                    opt.clipTips = true;
                    break;
                case 'd':
                    opt.deleteIsolated = true;
                    break;
                case 'v':
                    opt.verbose = true;
                    break;
                case 'c':
                    opt.outputColors = true;
                    break;
                /*case 'y':
                    opt.useMercyKmers = true;
                    break;*/
                case 'f':
                    opt.outputGFA = false;
                    opt.outputFASTA = true;
                    break;
                case 'b':
                    opt.outputGFA = false;
                    opt.outputBFG = true;
                    break;
                case 'n':
                    opt.compressOutput = false;
                    break;
                case 'N':
                    opt.writeIndexFile = false;
                    break;
                default: break;
            }
        }
    }

    return 0;
}

bool check_ProgramOptions(CCDBG_Build_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [&](vector<string>& v_files) {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "Error: File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const int format = FileParser::getFileFormat(file.c_str());

                if (format >= 0) files_tmp.push_back(file); // File is FASTA/FASTQ/GFA/GRAPH.BFG
                else {

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != NULL){

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        size_t i = 0;

                        while (i_file_txt.getline(buffer, 4096).good()){

                            fp = fopen(buffer, "r");

                            if (fp == NULL) {

                                cerr << "Error: Could not open file at line " << i << " in file " << file << " for reading." << endl;
                                ret = false;
                                break;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }

                            ++i;
                        }

                        if (i_file_txt.fail() && (i == 0)) {

                            cerr << "Error: File " << file << " is neither FASTA, FASTQ nor GFA." << endl;
                            cerr << "If it is a list of files, it is either empty or has a line with >4096 characters." << endl;
                            ret = false;
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "Error: Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
            }
        }

        v_files = move(files_tmp);

        delete[] buffer;
    };

    // Check general parameters

    if (!opt.build && !opt.update && !opt.query){

        cerr << "Error: No command selected (can be 'build' or 'update' or 'query')." << endl;
        ret = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if (opt.k <= 2){

        cerr << "Error: Length k of k-mers cannot be less than 3." << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed " << (MAX_KMER_SIZE - 1) << "." << endl;
        cerr << "To enable a larger k, recompile Bifrost with the appropriate MAX_KMER_SIZE variable." << endl;
        ret = false;
    }

    if (opt.g == 0){

        cerr << "Error: Length m of minimizers cannot be equal to 0." << endl;
        ret = false;
    }

    if ((opt.g >= 0) && (opt.g > opt.k - 2)) {

        cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
        ret = false;
    }

    if (opt.query){  // Check param. command build

        if (opt.prefixFilenameOut.length() == 0) {

            cerr << "Error: No output filename prefix given." << endl;
            ret = false;
        }
        else {

            const string out = opt.prefixFilenameOut + ".tsv";

            FILE* fp = fopen(out.c_str(), "w");

            if (fp == NULL) {

                cerr << "Error: Could not open file for writing output of query in TSV format: " << out << "." << endl;
                ret = false;
            }
            else {

                fclose(fp);
                if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
            }
        }

        if (opt.filename_query_in.size() == 0) {

            cerr << "Error: Missing input query files." << endl;
            ret = false;
        }
        else check_files(opt.filename_query_in);

        if ((opt.ratio_kmers < 0.0) || (opt.ratio_kmers > 1.0)) {

            cerr << "Error: Ratio of k-mers from queries that must occur in the graph cannot be less than 0.0 or more than 1.0 (" << opt.ratio_kmers << ")." << endl;
            ret = false;
        }

        if (opt.g > opt.k - 2){

            cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
            ret = false;
        }
    }
    else {

        if (opt.prefixFilenameOut.length() == 0) {

            cerr << "Error: No output filename prefix given." << endl;
            ret = false;
        }
        else {

            const string out = opt.prefixFilenameOut + (opt.outputGFA ? ".gfa" : ".fasta");

            FILE* fp = fopen(out.c_str(), "w");

            if (fp == NULL) {

                cerr << "Error: Could not open file for writing output graph in GFA format: " << out << "." << endl;
                ret = false;
            }
            else {

                fclose(fp);
                if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
            }
        }

        if ((opt.filename_seq_in.size() + opt.filename_ref_in.size()) == 0) {

            cerr << "Error: Missing input files." << endl;
            ret = false;
        }
        else {

            check_files(opt.filename_seq_in);
            check_files(opt.filename_ref_in);
        }
    }

    if (opt.build){ // Check param. command build

        if (opt.outFilenameBBF.length() != 0){

            FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

            if (fp == NULL) {

                cerr << "Error: Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                ret = false;
            }
            else {

                fclose(fp);

                if (remove(opt.outFilenameBBF.c_str()) != 0){

                    cerr << "Error: Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
                }
            }
        }

        if (opt.inFilenameBBF.length() != 0){

            if (check_file_exists(opt.inFilenameBBF)){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "Error: Input Blocked Bloom filter " << opt.inFilenameBBF << " file does not exist." << endl;
                ret = false;
            }
        }
    }

    if (opt.build || opt.update) {

        if (opt.outputFASTA && opt.outputBFG) {

            cerr << "Error: Two output format selected: FASTA and GFA. Rerun with one format selected." << endl;
            ret = false;
        }

        if (opt.outputBFG && !opt.writeIndexFile) {

            cerr << "Error: Bifrost index file output (.bfi) is required to output Bifrost graph file (.bfg). Remove the no index output argument." << endl;
            ret = false;
        }

    }

    if (opt.update || opt.query) {

        if (opt.filename_graph_in.length() == 0){

            cerr << "Error: No graph file was provided in input." << endl;
            ret = false;
        }
        else if (!check_file_exists(opt.filename_graph_in)){

            cerr << "Error: The graph file does not exist or is not a valid input file format." << endl;
            ret = false;
        }
        else {

            FILE* fp = fopen(opt.filename_graph_in.c_str(), "r");

            if (fp == NULL) {

                cerr << "Error: Could not read input graph file " << opt.filename_graph_in << "." << endl;
                ret = false;
            }
            else fclose(fp);
        }

        if (opt.filename_colors_in.length() != 0){

            if (!check_file_exists(opt.filename_colors_in)){

                cerr << "Error: The input color file does not exist or is not a valid input file format." << endl;
                ret = false;
            }
            else {

                FILE* fp = fopen(opt.filename_colors_in.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input color file " << opt.filename_colors_in << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
        }

        if (opt.filename_index_in.length() != 0){

            if (!check_file_exists(opt.filename_index_in)){

                cerr << "Error: The input index file does not exist or is not a valid input file format." << endl;
                ret = false;
            }
            else {

                FILE* fp = fopen(opt.filename_index_in.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read index file " << opt.filename_index_in << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
        }
    }

    return ret;
}

int main(int argc, char **argv){

    if (argc < 2) PrintUsage();
    else {

        CCDBG_Build_opt opt;

        opt.outputColors = false; // We dont know yet if we want colors or not

        const int print = parse_ProgramOptions(argc, argv, opt); // Parse input parameters

        if (print == 1) PrintVersion();
        else if (print == 2) PrintUsage();
        else if (check_ProgramOptions(opt)) {

            bool success = true; // Abort if any operation goes wrong
            
            if (opt.build){ // Build the graph

                if (opt.outputColors){

                    ColoredCDBG<> ccdbg(opt.k, opt.g);

                    success = ccdbg.buildGraph(opt);

                    if (success) success = ccdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                    if (success) success = ccdbg.buildColors(opt);
                    if (success) success = ccdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.writeIndexFile, opt.compressOutput, opt.verbose);
                }
                else {

                    CompactedDBG<> cdbg(opt.k, opt.g);

                    success = cdbg.build(opt);

                    if (success) success = cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);
                    if (success) success = cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.outputGFA, opt.outputFASTA, opt.outputBFG, opt.writeIndexFile, opt.compressOutput, opt.verbose);
                }
            }
            else if (opt.update){

                CCDBG_Build_opt lopt = opt;

                if (lopt.filename_colors_in.size() != 0){ // If colors in or out

                    ColoredCDBG<> ccdbg1(lopt.k, lopt.g);

                    if (lopt.filename_index_in.length() == 0) success = ccdbg1.read(lopt.filename_graph_in, lopt.filename_colors_in, lopt.nb_threads, lopt.verbose);
                    else success = ccdbg1.read(lopt.filename_graph_in, lopt.filename_index_in, lopt.filename_colors_in, lopt.nb_threads, lopt.verbose);

                    if (success) {

                        lopt.k = ccdbg1.getK();
                        lopt.g = ccdbg1.getG();

                        ColoredCDBG<> ccdbg2(lopt.k, lopt.g);

                        if (success) success = ccdbg2.buildGraph(lopt);
                        if (success) success = ccdbg2.buildColors(lopt);

                        if (success) {

                            const size_t ccdbg1_len = ccdbg1.length();
                            const size_t ccdbg2_len = ccdbg2.length();

                            ColoredCDBG<>& ccdbg_a = (ccdbg1_len > ccdbg2_len) ? ccdbg1 : ccdbg2;
                            ColoredCDBG<>& ccdbg_b = (ccdbg1_len > ccdbg2_len) ? ccdbg2 : ccdbg1;

                            if (success) success = ccdbg_a.merge(move(ccdbg_b), lopt.nb_threads, lopt.verbose);

                            if (success) success = ccdbg_a.simplify(lopt.deleteIsolated, lopt.clipTips, lopt.verbose);
                            if (success) success = ccdbg_a.write(lopt.prefixFilenameOut, lopt.nb_threads, lopt.writeIndexFile, lopt.compressOutput, lopt.verbose);
                        }
                    }
                }
                else {

                    CompactedDBG<> cdbg1(lopt.k, lopt.g);

                    if (lopt.filename_index_in.length() == 0) success = cdbg1.read(lopt.filename_graph_in, lopt.nb_threads, lopt.verbose);
                    else success = cdbg1.read(lopt.filename_graph_in, lopt.filename_index_in, lopt.nb_threads, lopt.verbose);

                    if (success) {

                        lopt.k = cdbg1.getK();
                        lopt.g = cdbg1.getG();

                        CompactedDBG<> cdbg2(lopt.k, lopt.g);

                        if (success) success = cdbg2.build(lopt);
                        if (success) {

                            const size_t cdbg1_len = cdbg1.length();
                            const size_t cdbg2_len = cdbg2.length();

                            CompactedDBG<>& cdbg_a = (cdbg1_len > cdbg2_len) ? cdbg1 : cdbg2;
                            CompactedDBG<>& cdbg_b = (cdbg1_len > cdbg2_len) ? cdbg2 : cdbg1;

                            if (success) {

                                success = cdbg_a.merge(cdbg_b, lopt.nb_threads, lopt.verbose);
                                cdbg_b.clear();
                            }

                            if (success) success = cdbg_a.simplify(lopt.deleteIsolated, lopt.clipTips, lopt.verbose);
                            if (success) success = cdbg_a.write(lopt.prefixFilenameOut, lopt.nb_threads, lopt.outputGFA,
                                                                lopt.outputFASTA, lopt.outputBFG, lopt.writeIndexFile, lopt.compressOutput, lopt.verbose);
                        }
                    }
                }
            }
            else if (opt.query){

                if (opt.filename_colors_in.size() != 0){

                    ColoredCDBG<> ccdbg(opt.k, opt.g);

                    if (opt.filename_index_in.length() == 0) success = ccdbg.read(opt.filename_graph_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);
                    else success = ccdbg.read(opt.filename_graph_in, opt.filename_index_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);

                    if (success) success = ccdbg.search(opt.filename_query_in, opt.prefixFilenameOut, opt.ratio_kmers, opt.inexact_search, opt.nb_threads, opt.verbose);
                }
                else {

                    CompactedDBG<> cdbg(opt.k, opt.g);

                    if (opt.filename_index_in.length() == 0) success = cdbg.read(opt.filename_graph_in, opt.nb_threads, opt.verbose);
                    else success = cdbg.read(opt.filename_graph_in, opt.filename_index_in, opt.nb_threads, opt.verbose);

                    if (success) success = cdbg.search(opt.filename_query_in, opt.prefixFilenameOut, opt.ratio_kmers, opt.inexact_search, opt.nb_threads, opt.verbose);
                }
            }

            if (!success) {

                cerr << "Operation aborted." << endl;
                exit(1); 
            }
        }
    }
}
