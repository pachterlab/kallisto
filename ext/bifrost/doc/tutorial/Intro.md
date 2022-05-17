# Bifrost API

Welcome to the introduction of the Bifrost API. In this tutorial, you will learn how to use the Bifrost C++ API to index, modify, query and color compacted de Bruijn graphs. Suggestions are always welcome so do not hesitate to leave some feedbacks.

### Getting starting

This introduction assumes you are familiar with C++ although no advanced knowledge of the language is required. This tutorial is in construction so I will detail and clarify some parts of it over time. It is additionally assumed you are familiar with the general concept of de Bruijn graphs.

I will try to cover as much ground as possible in this tutorial. However, the keyword of this tutorial is 'simplicity' so I will not detail every single function of the API. If you are particularly interested in a function, the documentation (located in `/doc/doxygen/`) describes in detail what it does, what is the input and what is the output. If a function or a type requires special attention, I will point out that having a look at the documentation would be beneficial.

## Table of Contents

* [Compacted de Bruijn graph](#compacted-de-bruijn-graph)
	* [Creating a graph](#creating-a-graph)
	* [Adding sequences](#adding-sequences)
	* [Finding *k*-mers](#finding-k-mers)
	* [Deleting sequences](#deleting-sequences)
	* [Storing unitigs identity](#storing-unitigs-identity)
	* [Building the graph from a data set](#building-the-graph-from-a-data-set)
	* [Cleaning the graph](#cleaning-the-graph)
	* [Reading and writing graphs](#reading-and-writing-graphs)
	* [Traversing the graph](#traversing-the-graph)
	* [Adding data to unitigs](#adding-data-to-unitigs)

## Compacted de Bruijn graph

Graphs in Bifrost are compacted bi-directed de Bruijn graphs. It is important to understand this concept before getting started as the API relies on these notions. Unitigs, compaction and bi-directed edges are all explained in this excellent [short tutorial](https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md) by Paul Medvedev and Rayan Chikhi.
I will try to summarize quickly some of the important concepts for Bifrost:

- Bifrost uses a node-centric representation of the graph. Vertices are represented explicitly while edges are represented implicitly.
- Vertices are unitigs, i.e., sequences of length greater or equal to *k* (the *k*-mer size). Hence, a unitig is composed of at least one *k*-mer.
- *K*-mers in a unitig are not branching (edge in-degree = edge out-degree = 1) in the graph **except** the first and last *k*-mers which **might** be branching.
- A *k*-mer occurs in **at most** one unitig.
- Because the graph is bi-directed, a unitig represents two sequences: itself and its reverse-complement. Hence, a unitig can be traversed/read in two different ways: from left to right (*foward* or *+* direction) or from right to left by complementing the DNA symbols (*reverse* or *-* direction). This is the *strandness* of the unitig.
- Edges are directed: they have a start unitig *A* and an end unitig *B*. Furthermore, edges embed the strandness of the unitigs they connect, which can be: {+,+}, {+,-}, {-,+} and {-,-}. An edge *e* connecting the forward sequence of *A* to the reverse-complement sequence of *B* would be *e={A+,B-}*.

### Creating a graph

Let's start by creating a compacted de Bruijn graph.
```cpp
#include <bifrost/CompactedDBG.hpp>

int main(int argc, char **argv){

	CompactedDBG<> cdbg;
}
```

The first line includes the Bifrost C++ header for compacted de Bruijn graphs and that will be the only `#include` you need. Then, you just created an empty compacted de Bruijn graph named `cdbg` in the main function. That's it. Without parameters, that graph has a default *k*-mer size of `k=31` but you can change that later. Note that if you want to use a *k*-mer size larger than 31, you have to compile the Bifrost library with an extra parameter as explained [here](https://github.com/pmelsted/bifrost#large-k-mers).
The `<>` in `CompactedDBG<>` is a template for the type of data that you want to associate to unitigs. For now, we assume you do not want to associate data to unitigs so it is empty. Later in this tutorial, we will cover how to create a graph of type `CompactedDBG<MyData>` where type `MyData` represents some data you want to associate to unitigs.

In the next code samples, I will skip the include and the main function for simplicity, unless it is required to understand the code snippet.

Let's initialize our de Bruijn graph with a specific *k*-mer size that we will subsequently print.
```cpp
const size_t k = 31;

CompactedDBG<> cdbg(k);

cout << "K-mer size is " << cdbg.getK() << endl;
```

### Adding sequences

Now that we have an empty graph, let's add some data in it:
```cpp
const string seq = "ACGTCGTACGTCCCGTAAACGTTAAACGTAAACGTGTGTGCAAAATGTCTAGTTTTTTTACGCTGATATAGTC";

cdbg.add(seq);
```

Function `CompactedDBG::add()` takes as input a string (a DNA sequence) and inserts it in the graph. The function takes care of splitting the input sequence into unitigs. *K*-mers overlapping non-DNA characters will be automatically discarded.

Since it is the first time we modify a graph, here are two properties of Bifrost to keep in mind:
- Bifrost graphs are always compacted, no matter what. If you edit the graph by adding sequences, removing unitigs or merging graphs, Bifrost will *always* take care of compacting the graph. It is not possible to have an intermediate state where a Bifrost graph is not compacted.
- Most functions in the Bifrost API have an optional parameter *verbose* for printing information messages about the execution of the function. By default, `verbose=false` so no information messages are printed. By setting it to `true`, information message will be printed to the standard output `stdout`.

### Finding ***k***-mers

Searching *k*-mers in the graph is one of the most common task to perform. It works as follows:
```cpp
const string kmer_sequence = "ACGTCGTACGTCCCGTAAACGTTAAACGTAA";
const Kmer km = Kmer(kmer_sequence.c_str());

UnitigMap<> um = cdbg.find(km);
```

An object `Kmer` is first created from the *k*-mer sequence and it is passed to function `CompactedDBG::find()` which returns a `UnitigMap` object. Objects of type `UnitigMap` are fundamental in Bifrost as they are used for nearly everything. A `UnitigMap` object represents the mapping of a sequence on a unitig of the graph. Hence, when searching for a *k*-mer in the graph, the returned `UnitigMap` object is the mapping of the searched *k*-mer sequence on a unitig of the graph or an empty mapping if the *k*-mer is not found:
```cpp
if (um.isEmpty) cout << "Kmer " << kmer_sequence << " was not found" << endl;
else {

	const string unitig = um.referenceUnitigToString();
	const string strandness = um.strand ? "forward" : "reverse-complement";
	const size_t position = um.dist;

	cout << "Kmer " << kmer_sequence << " was found in the " << strandness << " direction of unitig " << unitig << endl;
}
```

The important parameters of `UnitigMap` objects are:
- `UnitigMap::isEmpty`: if `true`, the mapping is empty: the searched sequence was not found in the graph.
- `UnitigMap::dist`: start position of the searched sequence on the unitig
- `UnitigMap::len`: length of the mapping **in *k*-mers** (at least 1)
- `UnitigMap::strand`: strandness of the mapped sequence, `true` for forward and `false` for reverse-complement
- `UnitigMap::size`: length of the unitig **in bp** (at least *k*)

The type `UnitigMap` has tons of useful functions so it is probably a good idea to have a look at its documentation. Here is a quick overview:
- `UnitigMap::referenceUnitigToString()`: outputs the string of the unitig associated with this mapping.
- `UnitigMap::mappedSequenceToString()`: outputs the string of the mapped sequence. Takes into account the strandness.
- `UnitigMap::getUnitigHead()`: returns the first *k*-mer of the unitig associated with this mapping.
- `UnitigMap::getMappedHead()`: returns the first *k*-mer of the unitig in the mapped sequence. takes into account the mapping strandness.
- ...

As you can see above, some functions take into account the strandness and some functions do not. Here is a little bit of terminology to figure out what is what:
- reference unitig: the unitig to which the sequence map to, as stored in the graph (in forward direction)
- mapped sequence: a substring of the reference unitig, in forward or reverse-complement direction

All `UnitigMap` functions containing `referenceUnitig` or just `Unitig` return something related to the reference unitig. Hence, parameter `UnitigMap:strand` is not used here. All `UnitigMap` functions containing `mappedSequence` or just `Mapped` return something related to the mapped sequence (a substring of the reference unitig) and takes into account parameter `UnitigMap::strand`. Here is a figure to explain it all:

### Deleting sequences

Let assume that we want to delete the unitig containing a *k*-mer `km`. `CompactedDBG<>::find()` was used to locate the *k*-mer in the graph and returned a `UnitigMap` object `um`. Removing from the graph the reference unitig used in `um` is:
```cpp
cdbg.remove(um);
```

This function removes the unitig entirely though. What if you want to remove only a substring, say just *k*mer `km`?
```cpp
const string unitig = um.referenceUnitigToString();
const string unitig_pre = unitig.substr(0, um.dist + k - 1); // Unitig substring 'before' k-mer, i.e, prefix
const string unitig_suf = unitig.substr(um.dist + 1, um.size - um.dist - 1); // Unitig substring 'after' k-mer, i.e, suffix

cdbg.remove(um); // remove unitig
cdbg.add(unitig_pre); // Insert prefix
cdbg.add(unitig_suf); // Insert suffix
```

To remove a substring, you must remove the unitig entirely and re-insert the unitig prefix and suffix surrounding the substring to remove.

### Storing unitigs identity

A recurrent question is 'How do I save a unitig identifier or its position in the graph?' which is mostly the same as 'Can I access a unitig in the graph using an identifier?'. With that regard, graphs in Bifrost are like hash tables (type `map`) in C++: You cannot access elements of the data structure with an identifier or a position but you can access them using an iterator. Let assume that for the purpose of your program, you want to "save" a reference to unitig represented by a `UnitigMap` object `um`:
```cpp
vector<UnitigMap<>> v_um;

v_um.push_back(um);
```

That is the closest equivalent of storing an iterator, quick and fast. Beware that same as for `vector` and `map` iterators in C++, the `UnitigMap` object is only valid and functional as long as you do **not** modify the data structure. Using the stored `UnitigMap` object after modifying the graph is undefined behavior and your program might very well crash. The reason is that Bifrost automatically re-compacts the graph after modification: some unitigs might split or merge. Hence, the unitig pointed out by the stored `UnitigMap` object might look different or just not exist anymore in the graph.

There is another way to store unitig identitiers or positions which will be ideal in a number of situations. Remember that each *k*-mer in the graph occur in **at most** one unitig. Which means that a *k*-mer can be used as an identifier for a unitig:
```cpp
vector<Kmer> v_km;

for (auto& um : v_um) {

	const Kmer head_kmer = um.getUnitigHead(); // First k-mer of the unitig

	v_km.push_back(head_kmer); // Push the head k-mer to a vector
}

v_um.clear(); // Vector of UnitigMap is not needed anymore
```
The idea is to use the head *k*-mer of the unitig as its identifier. Hence, instead of having a vector of `UnitigMap`, you have a vector of `Kmer` for which each *k*-mer is the head *k*-mer of a unitig you want to store. To retrieve the unitigs associated to these *k*mers:
```cpp
for (auto& km : v_km) {

	UnitigMap<> um = cdbg.find(km, true);
}
```

The major advantage of using `Kmer` over `UnitigMap` is that a `Kmer` object is a lot less memory consuming than storing a `UnitigMap` object, which is useful for storing lots of unitig identifiers. On the other hand, retrieving the unitig associated to a `Kmer` costs some time. It is a tradeoff you have to decide for yourself. My advise: if you have only a few hundreds/thousands of unitigs to remember, store `UnitigMap` objects. Otherwise, store `Kmer` objects.

Not that in the previous code snippet, I used `cdbg.find(km, true)`. When this last parameter is set to `true` (`false` by default), it indicates that you want to search for this *k*-mer **only** at the extremities of unitigs (head or tail *k*-mers only). Doing this significantly speeds up the search and it comes very handy when you search for *k*-mers that you know are the head *k*-mer of unitigs.

### Building the graph from a data set

Rather than manually adding sequences in the graph, Bifrost enables developers to construct the graph directly from data sets. For this, two modes are proposed:
- *reference* mode: the graph is built from **all** *k*-mers in the input data sets.
- *sequence* mode: the graph is built from *k*-mers occuring **twice or more** in the input data sets.

The reference mode is ideal for building the graph from reference genomes while the sequence mode is better for reads. Note that you can build your graph from assembled genomes and reads by combining both modes! Parameters of the graph construction are set in an object of type `CDBG_Build_opt` that you must configure:
- `CDBG_Build_opt::k`: Length of k-mers.
- `CDBG_Build_opt::filename_seq_in`: Vector of input filenames for the sequence mode.
- `CDBG_Build_opt::filename_ref_in`: Vector of input filenames for the reference mode.
- `CDBG_Build_opt::nb_threads`: Number of threads to use.
- `CDBG_Build_opt::verbose`: Print information messages during the construction.

Input filenames can be provided in FASTA or FASTQ format, gzipped or not, and GFA. Indeed, if you already have a graph in GFA format but not necessarily a de Bruijn graph, Bifrost will take its sequences to build a compacted de Bruijn graph. You can also provide in input a list of filenames as a text file with one filename per line (be careful to not have empty lines). Once you're done with setting these parameters, constructing the graph is as easy as:
```cpp
CDBG_Build_opt opt;

opt.k = 31;
opt.nb_threads = 4;
opt.verbose = true;

opt.filename_ref_in.push_back(string("my_assembled_genome.fasta"));
opt.filename_seq_in.push_back(string("my_short_reads.fastq"));

CompactedDBG<> cdbg(opt.k);

cdbg.build(opt);
```

### Cleaning the graph

### Reading and writing graphs

Bifrost offers to read and write graphs with the GFA and FASTA file formats. The default is GFA output: It is a plain text, tabulation formatted file format that explicitely describes a sequence graph. It has a short and simple [specification](http://gfa-spec.github.io/GFA-spec/GFA1.html), it is easy to visualize with [Bandage](https://rrwick.github.io/Bandage/) and it is now a community standard for sequence graph tools in computational biology. Note that Bifrost outputs specifically GFA v1 and it only uses the Segment (S) and Link (L) fields. Another file format that Bifrost uses to store graphs is FASTA. FASTA was not originally designed to store graphs but as Bifrost has an implicit representation of edges, storing only the unitigs in FASTA is enough. Hence, the advantage of exporting the graph to FASTA over GFA is file size as edges are not stored. Finally, GFA and FASTA are very compressible file formats so you can easily bgzip them.

Writing the graph works as follow:
```cpp
const string output_filename = "/my/output/graph";
const size_t nb_threads = 4;
const bool gfa_output = true;
const bool verbose = true;


cdbg.write(output_filename, nb_threads, gfa_output, verbose);
```

and reading works in a similar fashion:
```cpp
CompactedDBG<> cdbg(k);

const string input_filename = "/my/input/graph.gfa";
const size_t nb_threads = 4;
const bool verbose = true;

cdbg.read(input_filename, nb_threads, verbose);
```

Input file format, GFA or FASTA, is automatically recognized when reading. A detail you might want to have a closer look at is the *k*-mer size. Indeed, if the input graph has been built with Bifrost and is in GFA format, the file already contains the *k*-mer size and its the one that will be used, you don't have to do anything. If you input a compacted de Bruijn graph not built with Bifrost or a graph in FASTA format, Bifrost will try to read the graph with the *k*-mer size given at the declaration of the graph, i.e., `CompactedDBG<> cdbg(k)`. In that case, it is your responsability to make sure the *k*-mer size is correct and the graph is correctly compacted if built with a different tool than Bifrost.

### Traversing the graph

Traversing the graph requires first a starting vertex represented by a `UnitigMap<>` object. Assuming a starting `UnitigMap<> um` object, traversing the neighboring vertices is then as simple as:
```cpp
for (const auto& um_succ : um.getSuccessors()) {

	cout << um_succ.referenceUnitigToString() << endl;
	cout << um_succ.mappedSequenceToString() << endl;
}

for (const auto& um_pred : um.getPredecessors()) {

	cout << um_pred.referenceUnitigToString() << endl;
	cout << um_pred.mappedSequenceToString() << endl;
}
```

Function `UnitigMap::getSuccessors()` and `UnitigMap::getPredecessors()` return the successors and predecessors, respectively, under the form of `UnitigMap` objects (`um_succ` and `um_pred` in the previous example). As the traversal functions work at the unitig level, it does not matter whether the starting `UnitigMap<> um` object maps to a single *k*-mer or the full unitig: the traversal functions will consider the full unitig as the starting vertex. However, `getSuccessors()` and `getPredecessors()` will always return `UnitigMap` objects that are the mapping of a full unitig sequence to itself, i.e, `um_succ.dist = 0` and `um_succ.len = um_succ.size - k + 1`.
**What is really important to consider when traversing is the strandness of the unitigs.** As described in Section [Compacted de Bruijn graph](#compacted-de-bruijn-graph), edges in bi-directed de Bruijn graphs embed the strandness of the unitigs they connect. In the previous example, the strandness of the start unitig is `um.strand`, the strandness of the successors is `um_succ.strand` and the strandness of the predecessors is `um_pred.strand`. What it means in practice is that the successors of `um` with `um.strand = true` are not going to be the same as with `um.strand = false`. To be accurate, the successors of `um` with `um.strand = true` are going to be the predecessors of `um` with `um.strand = false`. This can be a little tricky to understand so let's have an example:
```cpp
cout << "--- Printing successors in forward direction ---" << endl;

um.strand = true;

for (const auto& um_succ : um.getSuccessors()) {

	cout << um_succ.mappedSequenceToString() << end;
}

cout << "--- Printing predecessors in backward direction ---" << endl;

um.strand = false;

for (const auto& um_pred : um.getPredecessors()) {

	cout << um_pred.mappedSequenceToString() << end;
}
```

In this example, we print the unitig sequence of the successors of `um` in forward direction (`um.strand = true`) and the unitig sequence of the predecessors of `um` in backward direction (`um.strand = false`). And to no surprise, they are the same. One last consideration: because unitig strandness changes when traversing the graph, accessing the unitig sequence is done with `UnitigMap::mappedSequenceToString()` as shown in the previous code snippet. Indeed, `mappedSequenceToString()` always takes into account the strandness of the mapped sequence while `referenceUnitigToString()` gives you the unitig sequence the way it is store in the graph, i.e, with an arbitrary strandness.

### Adding data to unitigs
