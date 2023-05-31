# Bifrost

### Parallel construction, indexing and querying of colored and compacted de Bruijn graphs

* **Build**, **index**, **color** and **query** the compacted de Bruijn graph
* **No need to build the uncompacted** de Bruijn graph
* **Reads** or **assembled genomes** as input
* Output **graph in GFA** (can be visualized with [Bandage](https://github.com/rrwick/Bandage))
* **Graph cleaning**: short tip clipping, etc.
* **No disk** usage (adapted for cluster architectures)
* **Multi-threaded** and **SIMD** optimized
* **No parameters to estimate** with other tools
* **Inexact** *k*-mer search of queries
* **C++ API** available:
    * Associate **your data with vertices**
    * **Add** or **remove** (sub-)sequences / *k*-mers / colors
    * **Find unitigs** containing **queried k-mers**

## Table of Contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Binary usage](#binary-usage)
* [API](#api)
* [FAQ](#faq)
* [Troubleshooting](#troubleshooting)
* [Citation](#citation)
* [Contact](#contact)
* [License](#license)

## Requirements

To install Bifrost using Bioconda or Brew, go directly to Section [Installation](#installation). To install from source, you will need:

* C++11 compiler:
    * [GCC](https://gcc.gnu.org/) >= 5.1.0
    * [Clang](http://clang.llvm.org/) >= 3.5
* [Cmake](https://cmake.org/) >= 2.8.12
* [Zlib](https://zlib.net/)

All are probably already installed on your computer as those are installed by default on most operating systems. They can be downloaded and installed by following the instructions on their respective websites. However, it is most likely they are all available via a package manager for your operating system: 

* **Ubuntu/Debian**:
```
sudo apt-get install build-essential cmake zlib1g-dev
```
* **MacOS** (with [Homebrew](https://brew.sh/)):
```
brew install --with-toolchain llvm
brew install cmake zlib
```
* **Windows 10**:

1. Open the Windows Store
2. Search and install the `Ubuntu` app (from `Canonical Group Limited`)
3. Open the Windows main menu and open the `Ubuntu` app (it should open an Ubuntu terminal)
4. Use the following command in the Ubuntu terminal:
```
sudo apt-get install build-essential cmake zlib1g-dev
```
5. Use the opened Ubuntu terminal for compiling, installing and running Bifrost (see next section). See [Troubleshooting](#troubleshooting) if you have any problem during the installation.

## Installation

Compared to the source install, the Conda package do not support *k>31* nor native compilation (including AVX2 instructions). Use the source installation for benchmarking.

* From [Bioconda](https://bioconda.github.io):

  ```
  conda -c bioconda bifrost
  ```

* From source

  ```
  git clone https://github.com/pmelsted/bifrost.git
  cd bifrost && mkdir build && cd build
  cmake ..
  make
  make install
  ```

  `make install` might require `sudo` (`sudo make install`) to proceed. To install Bifrost in the non-default path `/some/path/`, add the option `-DCMAKE_INSTALL_PREFIX=/some/path/` to the `cmake` command.

  By default, the installation creates:
  * a binary (*Bifrost*)
  * a dynamic library (*libbifrost.so* for Unix or *libbifrost.dylib* for MacOS)
  * a static library (*libbifrost.a*)

  **Advanced options**
  * Bifrost compiles by default with `-march=native`: the compiler targets architecture instructions specific to the machine Bifrost is compiled on. Hence, the binary and library produced might not work on a different machine. Native compilation can be disabled by adding the option `-DCOMPILATION_ARCH=OFF` to the `cmake` command (disables all AVX2 optimizations too). Alternatively, you can use this option to specify the architecture you want to target: `x86-64`, `knl`, etc. Default is `-DCOMPILATION_ARCH=native`.
  * Bifrost uses AVX2 instructions during graph construction which can be disabled by adding the option `-DENABLE_AVX2=OFF` to the `cmake` command.

  If you encounter any problem during the installation, see the [Troubleshooting](#troubleshooting) section.

### Large *k*-mers

The default maximum *k*-mer size supported is 31. To work with larger *k* in the binary, you must install Bifrost from source and replace *MAX_KMER_SIZE* with a larger multiple of 32. This can be done in two ways:

* By adding the following option to the `cmake` command:
```
-DMAX_KMER_SIZE=64
```

* By replacing *MAX_KMER_SIZE* in *CMakeLists.txt*:
```
SET(MAX_KMER_SIZE "64" CACHE STRING "MAX_KMER_SIZE")
```

Actual maximum k-mer size is *MAX_KMER_SIZE-1*, e.g maximum *k* is 63 for *MAX_KMER_SIZE=64*. Increasing *MAX_KMER_SIZE* increases Bifrost memory usage (*k*=31 uses 8 bytes of memory per *k*-mer while *k*=63 uses 16 bytes of memory per *k*-mer).

The maximum size of minimizers (*g*-mers) *MAX_GMER_SIZE* can be adjusted the same way as *MAX_KMER_SIZE*. This is especially useful if you want to use a large *k*-mer size but a small *g*-mer size. By default, *MAX_GMER_SIZE* is equal to *MAX_KMER_SIZE*.

To work with larger *k* when using the Bifrost API, the new value *MAX_KMER_SIZE* must be given to the compiler and linker as explained in Section [API](#api)

## Binary usage:

```
Bifrost
```

displays the command line interface:
```
Bifrost x.y

Highly parallel construction, indexing and querying of colored and compacted de Bruijn graphs

Usage: Bifrost [COMMAND] [PARAMETERS]

[COMMAND]:

   build                   Build a compacted de Bruijn graph, with or without colors
   update                  Update a compacted (possible colored) de Bruijn graph with new sequences
   query                   Query a compacted (possible colored) de Bruijn graph

[PARAMETERS]: build

   > Mandatory with required argument:

   -s, --input-seq-file     Input sequence file (FASTA/FASTQ possibly gzipped)
                            Multiple files can be provided as a list in a TXT file (one file per line)
                            K-mers with exactly 1 occurrence in the input sequence files will be discarded
   -r, --input-ref-file     Input reference file (FASTA/FASTQ possibly gzipped and GFA)
                            Multiple files can be provided as a list in a TXT file (one file per line)
                            All k-mers of the input reference files are used
   -o, --output-file        Prefix for output file(s)

   > Optional with required argument:

   -t, --threads            Number of threads (default is 1)
   -k, --kmer-length        Length of k-mers (default is 31)
   -m, --min-length         Length of minimizers (default is 23)
   -b, --bloom-bits         Number of Bloom filter bits per k-mer with 1+ occurrences in the input files (default is 14)
   -B, --bloom-bits2        Number of Bloom filter bits per k-mer with 2+ occurrences in the input files (default is 14)
   -l, --load-mbbf          Input Blocked Bloom Filter file, skips filtering step (default is no input)
   -w, --write-mbbf         Output Blocked Bloom Filter file (default is no output)
   -u, --chunk-size         Read chunk size per thread (default is 64)

   > Optional with no argument:

   -c, --colors             Color the compacted de Bruijn graph (default is no coloring)
   -y, --keep-mercy         Keep low coverage k-mers connecting tips
   -i, --clip-tips          Clip tips shorter than k k-mers in length
   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length
   -a, --fasta              Output file is in FASTA format (only sequences) instead of GFA
   -v, --verbose            Print information messages during execution

[PARAMETERS]: update

  > Mandatory with required argument:

   -g, --input-graph-file   Input graph file to update (GFA format)
   -s, --input-seq-file     Input sequence file (FASTA/FASTQ possibly gzipped)
                            Multiple files can be provided as a list in a TXT file (one file per line)
                            K-mers with exactly 1 occurrence in the input sequence files will be discarded
   -r, --input-ref-file     Input reference file (FASTA/FASTQ possibly gzipped and GFA)
                            Multiple files can be provided as a list in a TXT file (one file per line)
                            All k-mers of the input reference files are used
   -o, --output-file        Prefix for output file(s)

   > Optional with required argument:

   -f, --input-color-file   Input color file associated with the input graph file to update
   -t, --threads            Number of threads (default is 1)
   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)
   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or 23)

   > Optional with no argument:

   -i, --clip-tips          Clip tips shorter than k k-mers in length
   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length
   -v, --verbose            Print information messages during execution

[PARAMETERS]: query

  > Mandatory with required argument:

   -g, --input-graph-file   Input graph file to query (GFA format)
   -q, --input-query-file   Input query file (FASTA/FASTQ possibly gzipped)
                            Multiple files can be provided as a list in a TXT file (one file per line)
   -o, --output-file        Prefix for output file
   -e, --ratio-kmers        Ratio of k-mers from queries that must occur in the graph (default is 0.8)

   > Optional with required argument:

   -f, --input-color-file   Input color file associated with the input graph file to query
                            Presence/absence of queries will be output for each color
   -t, --threads            Number of threads (default is 1)
   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)
   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or 23)

   > Optional with no argument:

   -n, --inexact            Graph is searched with exact and inexact k-mers (1 substitution or indel) from queries         
   -v, --verbose            Print information messages during execution
```

### Examples

- **Build**

  1. **Build a compacted de Bruijn graph from read files and clean the graph**
     ```
     Bifrost build -t 4 -k 31 -i -d -s A.fastq -s B.fastq -o AB_graph 
     ```
     The compacted de Bruijn graph is built (`build`) with 4 threads (`-t 4`) from the 31-mers (`-k 31`) of files *A.fastq* and *B.fastq* (`-s A.fastq -s B.fastq`). By using parameter `-s`, files *A.fastq* and *B.fastq* are filtered: 31-mers occurring exactly once in *A* and *B* are discarded from the construction. Graph simplification steps are performed after building (`-i -d`) and the graph is written to file *AB_graph.gfa* (`-o AB_graph`).

  2. **Build a compacted de Bruijn graph from a reference genome file**
     ```
     Bifrost build -t 4 -k 31 -r C.fasta -o C_graph 
     ```
     The compacted de Bruijn graph is built (`build`) with 4 threads (`-t 4`) from the 31-mers (`-k 31`) of file *C.fasta* (`-r C.fasta`). By using parameter `-r`, file *C.fasta* is NOT filtered: all 31-mers occurring in *C* are used during the construction. The graph is written to file *C_graph.gfa* (`-o C_graph`).

  3. **Build a compacted and colored de Bruijn graph from read files and reference genome files, clean the graph**
     ```
     Bifrost build -t 4 -k 31 -c -i -d -s A.fastq -s B.fastq -r C.fasta -o ABC 
     ```
     Combining the two previous examples, the compacted de Bruijn graph is built (`build`) with 4 threads (`-t 4`) from the 31-mers (`-k 31`) of files *A.fastq*, *B.fastq* (`-s A.fastq -s B.fastq`) and file *C.fasta* (`-r C.fasta`). Graph simplification steps are performed after building (`-i -d`). The graph is colored (`-c`), meaning that each k-mer of the graph unitigs keeps track of whether it occurs in *A*, *B* or *C*. The graph is written to file *ABC.gfa* and the colors are written to file *ABC.bfg_colors* (`-o ABC`).

- **Update**

  1. **Update a compacted de Bruijn graph with a reference genome file**
     ```
     Bifrost update -t 4 -r D.fasta -g C_graph.gfa -o CD_graph 
     ```
     The compacted de Bruijn graph *C* (`-g C_graph.gfa`) is updated (`update`) with 4 threads (`-t 4`) from the *k*-mers of file *D.fasta* (`-r D.fasta`). By using parameter `-r`, file *D.fasta* is NOT filtered: all *k*-mers occurring in *D* are used during the merging. The graph is written to file *CD_graph.gfa* (`-o CD_graph`).

  2. **Update a compacted and colored de Bruijn graph with read files and clean the graph**
     ```
     Bifrost update -t 4 -i -d -s E.fastq -s F.fastq -g ABC.gfa -f ABC.bfg_colors -o ABCEF 
     ```
     The compacted and colored de Bruijn graph *ABC* (`-g ABC.gfa -f ABC.bfg_colors`) is updated (`update`) with 4 threads (`-t 4`) from the *k*-mers of files *E.fastq* and *F.fastq* (`-s E.fastq -s F.fastq`). Graph simplification steps are performed after merging (`-i -d`). The graph is written to file *ABCEF.gfa* and the colors are written to file *ABCEF.bfg_colors* (`-o ABCEF`).

- **Query**

  1. **Query a compacted de Bruijn graph for presence/absence of queries in the graph**
     ```
     Bifrost query -t 4 -e 0.8 -g ABCEF.gfa -q queries.fasta -o presence_queries 
     ```
     The compacted de Bruijn graph *ABCEF* (`-g ABCEF.gfa`) is queried (`query`) with 4 threads (`-t 4`) for the presence/absence of sequences from file *queries.fasta* (`-q queries.fasta`). At least 80% of each query *k*-mers must be found in the graph to have the query reported as present (`-e 0.8`). The results are stored in a binary matrix written to file *presence_queries.tsv* (`-o presence_queries`): rows are the queries, column is presence/absence in graph, intersection of a row and a column is a binary value indicating presence/absence of the query in graph (1 is present, 0 is not present).

  2. **Query a compacted de Bruijn graph for presence/absence of queries in the graph in inexact mode**
     ```
     Bifrost query -t 4 -e 0.8 -n -g ABCEF.gfa -q queries.fasta -o presence_queries 
     ```
     The compacted de Bruijn graph *ABCEF* (`-g ABCEF.gfa`) is queried (`query`) with 4 threads (`-t 4`) for the presence/absence of sequences from file *queries.fasta* (`-q queries.fasta`). At least 80% of each query *k*-mers must be found in the graph to have the query reported as present (`-e 0.8`). Queries are searched for exact and inexact *k*-mers (`-n`): *k*-mers with up to one substitution or indel. The results are stored in a binary matrix written to file *presence_queries.tsv* (`-o presence_queries`): rows are the queries, column is presence/absence in graph, intersection of a row and a column is a binary value indicating presence/absence of the query in graph (1 is present, 0 is not present).

  3. **Query a colored and compacted de Bruijn graph for presence/absence of queries in each color of the graph**
     ```
     Bifrost query -t 4 -e 0.8 -g ABCEF.gfa -f ABCEF.bfg_colors -q queries.fasta -o presence_queries 
     ```
     The compacted and colored de Bruijn graph *ABCEF* (`-g ABCEF.gfa -f ABCEF.bfg_colors`) is queried (`query`) with 4 threads (`-t 4`) for the sequences of file *queries.fasta* (`-q queries.fasta`). At least 80% of each query *k*-mers must be found in a color of the graph to have the query reported as present for that color (`-e 0.8`). The results are stored in a binary matrix written to file *presence_queries.tsv* (`-o presence_queries`): rows are the queries, columns are the colors, intersection of a row and a column is a binary value indicating presence/absence of the query in the color of the graph (1 is present, 0 is not present).

## API

Changes in the API are reported in the [Changelog](https://github.com/pmelsted/bifrost/blob/master/Changelog.md).

### Tutorial

The [API tutorial](doc/tutorial/Intro.md) should help you get started with the C++ API.

### Documentation

Documentation for the Bifrost library is available in the */doc/doxygen* folder (HTML version, open *html/index.html*).

The following command regenerates the documentation:
```
cd <bifrost_directory>
doxygen Doxyfile
```

The documentation contains a description of all the functions and structures of the library.

### Usage

The Bifrost C++ API can be used by adding
```
#include <bifrost/CompactedDBG.hpp>
```
for uncolored compacted de Bruijn graphs and
```
#include <bifrost/ColoredCDBG.hpp>
```
for colored compacted de Bruijn graphs in your C++ headers.

To compile, we recommend using the following compile flags:
```
-O3 -std=c++11
```
Furthermore, Bifrost compiles by default with flag `-march=native` so unless native compilation was disabled when installing Bifrost, use flag `-march=native` too.

Finally, use the following flags for linking:
```
-lbifrost -pthread -lz
```

You can also link to the Bifrost static library (*libbifrost.a*) for better performance:
```
<path_to_lib_folder>/libbifrost.a -pthread -lz
```

The default maximum *k*-mer size supported is 31. To work with larger *k*, the code using the Bifrost C++ API must be compiled and linked with the flag `-DMAX_KMER_SIZE=x` for compiling and linking where `x` is a larger multiple of 32, such as:
```
-DMAX_KMER_SIZE=64
```
Actual maximum k-mer size is *MAX_KMER_SIZE-1*, e.g maximum *k* is 63 for *MAX_KMER_SIZE=64*. Increasing *MAX_KMER_SIZE* increases Bifrost memory usage (*k*=31 uses 8 bytes of memory per *k*-mer while *k*=63 uses 16 bytes of memory per *k*-mer).

## FAQ

**Can I provide in input multiple files?**

Yes, use parameter `-r` or `-s` for each file to input.

**Can I provide in input a file which is a list of files?**

Yes, a text file containing one input filename per line with no empty lines can be given in input.

**What are the accepted input file formats?**

FASTA, FASTQ and GFA. Input FASTA and FASTQ files can be compressed with gzip (extension .gz). If you input a GFA file for the construction, you probably want to use the `-r` parameter for that file.

**Can I mix different file formats in input?**

Yes, as long as they are FASTA, FASTQ and GFA.

**If I input a GFA file for building the de Bruijn graph, does it need to contain an already compacted de Bruijn graph?**

No, it can contain any type of sequence graph (like an uncompacted de Bruijn graph or a sequence graph).

**Can I build a compacted (colored) de Bruijn graph from assembled genomes and reads?**

Yes. Input your assembled genomes with parameter `-r` and your reads with parameter `-s`.

**Can I use the graph file without its color file ?**

Yes. Just do not input the color file and Bifrost will consider it is an **un**colored compacted de Bruijn graph.

**In which order are inserted the colors?**

A color corresponds to an input file the graph was built/updated from. The order in which the colors are inserted is the same as the order of the files given by parameter `-r` and parameter `-s`. However, in case both parameters `-r` and `-s` are used, no assumption can be made on whether the files given by parameter `-s` will be inserted before or after the ones given by parameter `-r`.

**Different runs of Bifrost on the same dataset with the same parameters produces graphs with different unitigs. Which graph is correct?**

All of them. The difference between the graphs resides in circular unitigs (unitigs connecting to themselves) which are their own connected components ("isolated"). These unitigs can have a different sequence from one run to another because the starting position will be different, yet they represent exactly the same sequence. As an example, circular unitig ATAT composed of 3-mers can also be represented with sequence TATA. The number of unitigs will remain the same from one graph to another.

## Troubleshooting

* compilation (`make`) fails because some header files (*.h*) are not found

Assuming the header files (*.h*) are located at the path */usr/local/include/*, the following command set the environment variables *C_INCLUDE_PATH* and *CPLUS_INCLUDE_PATH* correctly for the time of the session:
```
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include/
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/include/
```

* executing the binary of Bifrost fails because *libbifrost.so* or *libbifrost.a* is not found

Assuming that *libbifrost*.(*so*|*dylib*|*a*) is located at the path */usr/local/lib/*, the following command set the environment variables *LD_LIBRARY_PATH*, *LIBRARY_PATH* and *PATH* correctly for the time of the session:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib/
export PATH=$PATH:/usr/local/lib/
```

* Bifrost crashes right at the beginning with error `Illegal instruction`

  You are most likely running Bifrost on a different machine than the one used to compile it. By default, Bifrost is compiled in native mode such to target architecture instructions specific to the machine it is compiled on. Using Bifrost on a different machine with a different architecture might result in this error. To solve this issue, Bifrost must be recompiled with native architecture compilation disabled, as explained in the Advanced options of Section [Installation](#installation).

## Citation

```
@article {holley2019bifrost,
  author = {Holley, Guillaume and Melsted, P{\'a}ll},
  title = "{Bifrost - Highly parallel construction and indexing of colored and compacted de Bruijn graphs}",
  elocation-id = {695338},
  doi = {10.1101/695338},
  journal = {bioRxiv},
  year = {2019}
}
```

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.

## License

* Bifrost is BSD2 licensed (https://github.com/pmelsted/bifrost/blob/master/LICENSE)
* The wyhash library is Unlicense licensed (https://github.com/wangyi-fudan/wyhash)
* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)
* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)
* The kseq library is copyrighted by Heng Li and released under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)
* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)
* The zstr library is MIT licensed (https://github.com/mateidavid/zstr)
* The GetRSS library is Creative Commons Attribution 3.0 licensed
