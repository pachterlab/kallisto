---
layout: page
title: "Manual"
description: ""
group: navigation
---
{% include JB/setup %}


#### Compiling


~~~
git clone https://github.com/pachterlab/kallisto
cd kallisto
mkdir build
cd build
cmake ..
make
make install
~~~

If you don't want to install kallisto system-wide and only want the executable binary, don't run make install -- make will already have generated the binary and have put it at: kallisto/build/src/kallisto

#### kallisto

Typing `kallisto` produces a list of usage options, which are:

~~~
kallisto 0.50.0

Usage: kallisto <CMD> [arguments] ..

Where <CMD> can be one of:

    index         Builds a kallisto index
    quant         Runs the quantification algorithm
    quant-tcc     Runs quantification on transcript-compatibility counts
    bus           Generate BUS files for single-cell data
    h5dump        Converts HDF5-formatted results to plaintext
    inspect       Inspects and gives information about an index
    version       Prints version information
    cite          Prints citation information

Running kallisto <CMD> without arguments prints usage information for <CMD>
~~~
The usage commands are:

#### index

`kallisto index` builds an index from a FASTA formatted file of target sequences. The arguments for the index command are:

~~~
kallisto 0.50.0
Builds a kallisto index

Usage: kallisto index [arguments] FASTA-files

Required argument:
-i, --index=STRING          Filename for the kallisto index to be constructed 

Optional argument:
-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)
-d, --d-list=STRING         Path to a FASTA-file containing sequences to mask from quantification
    --make-unique           Replace repeated target names with unique names
    --aa                    Generate index from a FASTA-file containing amino acid sequences
    --distinguish           Generate index where sequences are distinguished by the sequence name
-t, --threads=INT           Number of threads to use (default: 1)
-m, --min-size=INT          Length of minimizers (default: automatically chosen)
-e, --ec-max-size=INT       Maximum number of targets in an equivalence class (default: no maximum)
~~~

The Fasta file supplied can be either in plaintext or gzipped format. Note: Do not supply the genome Fasta file; the Fasta file must be a transcriptome Fasta.

#### quant

`kallisto quant` runs the quantification algorithm. The arguments for the quant command are:

~~~
kallisto 0.50.0
Computes equivalence classes for reads and quantifies abundances

Usage: kallisto quant [arguments] FASTQ-files

Required arguments:
-i, --index=STRING            Filename for the kallisto index to be used for
                              quantification
-o, --output-dir=STRING       Directory to write output to

Optional arguments:
-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
    --seed=INT                Seed for the bootstrap sampling (default: 42)
    --plaintext               Output plaintext instead of HDF5
    --single                  Quantify single-end reads
    --single-overhang         Include reads where unobserved rest of fragment is
                              predicted to lie outside a transcript
    --fr-stranded             Strand specific reads, first read forward
    --rf-stranded             Strand specific reads, first read reverse
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
                              (default: -l, -s values are estimated from paired
                               end data, but are required when using --single)
-t, --threads=INT             Number of threads to use (default: 1)
~~~

__kallisto__ can process either single-end or paired-end reads. The default running mode is paired-end and requires an even number of FASTQ files represented as pairs, e.g.

~~~
kallisto quant -i index -o output pairA_1.fastq pairA_2.fastq pairB_1.fastq pairB_2.fastq
~~~

For single-end mode you supply the `--single` flag, as well as the `-l` and `-s` options, and list any number of FASTQ files, e.g

~~~
kallisto quant -i index -o output --single -l 200 -s 20 file1.fastq.gz file2.fastq.gz file3.fastq.gz
~~~

FASTQ files can be either plaintext or
gzipped.

**Important note:** only supply one sample at a time to _kallisto_. The
multiple FASTQ (pair) option is for users who have samples that span multiple
FASTQ files.

In the case of single-end reads, the -l option must be used to specify the average fragment length. Typical Illumina libraries produce fragment lengths ranging from 180--200 bp but it's best to determine this from a library quantification with an instrument such as an Agilent Bioanalyzer. For paired-end reads, the average fragment length can be directly estimated from the reads and the program will do so if -l is not used (this is the preferred run mode). For reads that are produced by 3'-end sequencing, the `--single-overhang` option does not discard reads where the expected fragment size goes beyond the transcript start.

The number of bootstrap samples is specified using -b. Note that because of the large amount of data that may be produced when the number of bootstrap samples is high, __kallisto__ outputs bootstrap results in HDF5 format. The `h5dump` command can be used afterwards to convert this output to plaintext, however most convenient is to analyze bootstrap results with __sleuth__.

`kallisto quant` produces three __output__ files by default:

- __abundance.h5__ is a HDF5 binary file containing run info, abundance
  esimates, bootstrap estimates, and transcript length information length. This
  file can be read in by __sleuth__
- __abundance.tsv__ is a plaintext file of the abundance estimates. It does
  not contains bootstrap estimates. Please use the `--plaintext` mode to output
  plaintext abundance estimates. Alternatively, `kallisto h5dump` can be used
  to output an HDF5 file to plaintext. The first line contains a header for
  each column, including [estimated counts, TPM, effective length](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
- __run\_info.json__ is a json file containing information about the run


##### Optional arguments


+ `-t, --threads` specifies the number of threads to be used both for pseudoalignment and running bootstrap. The default value is 1 thread, specifying more than the number of bootstraps or the number of cores on your machine has no additional effect.

+ `--fr-stranded` runs kallisto in strand specific mode, only fragments where the **first** read in the pair pseudoaligns to the **forward** strand of a transcript are processed. If a fragment pseudoaligns to multiple transcripts, only the transcripts that are consistent with the first read are kept.

+ `--rf-stranded` same as `--fr-stranded` but the **first** read maps to the **reverse** strand of a transcript.


##### Legacy options:

These options are only supported in kallisto versions before 0.50.0. To use them install an older version of kallisto (recommended: 0.48.0).

kallisto quant:

+ `--bias` learns parameters for a model of sequences specific bias and corrects the abundances accordlingly.

+ `--fusion` does normal quantification, but additionally looks for reads that do not pseudoalign because they are potentially from fusion genes. All output is written to the file `fusion.txt` in the output folder.

`--pseudobam` outputs all pseudoalignments to a file `pseudoalignments.bam` in the output directory. This BAM file contains the pseudoalignments in BAM format, ordered by reads so that each pseudoalignment of a read is adjacent in the BAM file.

`--genomebam` constructs the pseudoalignments to the transcriptome, but projects the transcript alignments to genome coordinates, resulting in split-read alignments. When the `--genomebam` option is supplied at GTF file must be given with the `--gtf` option. The GTF file, which can be plain text or gzipped, translates transcripts into genomic coordinates. We recommend downloading a the cdna FASTA files and GTF files from the same data source. The `--chromosomes` option can provide a length of the genomic chromosomes, this option is not neccessary, but gives a more consistent BAM header, some programs may require this for downstream analysis. __kallisto__ does not require the genome sequence to do pseudoalignment, but downstream tools such as genome browsers will probably need it.


A detailed description of the SAM output is [here](pseudobam.html).


#### quant-tcc

`kallisto quant` runs the EM algorithm to produce estimated counts from a transcript-compatibility-counts matrix file (which is in a MatrixMarket format where each column is an equivalence class and each row is a sample). The necessary files can all be generated from `bustools count` in bustools.

~~~
kallisto 0.50.0
Generates BUS files for single-cell sequencing

Usage: kallisto bus [arguments] FASTQ-files

Required arguments:
-o, --output-dir=STRING       Directory to write output to

Optional arguments:
-i, --index=STRING            Filename for the kallisto index to be used
                              (required if file with names of transcripts not supplied)
-T, --txnames=STRING          File with names of transcripts
                              (required if index file not supplied)
-e, --ec-file=FILE            File containing equivalence classes
                              (default: equivalence classes are taken from the index)
-f, --fragment-file=FILE      File containing fragment length distribution
                              (default: effective length normalization is not performed)
--long                        Use version of EM for long reads
-p, --platform=STRING         [PacBio or ONT] used for sequencing
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
                              (note: -l, -s values only should be supplied when
                               effective length normalization needs to be performed
                               but --fragment-file is not specified)
-t, --threads=INT             Number of threads to use (default: 1)
-g, --genemap                 File for mapping transcripts to genes
                              (required for obtaining gene-level abundances)
-G, --gtf=FILE                GTF file for transcriptome information
                              (can be used instead of --genemap for obtaining gene-level abundances)
-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
    --matrix-to-files         Reorganize matrix output into abundance tsv files
    --matrix-to-directories   Reorganize matrix output into abundance tsv files across multiple directories
    --seed=INT                Seed for the bootstrap sampling (default: 42)
    --plaintext               Output plaintext only, not HDF5
~~~


#### bus

`kallisto bus` works with raw FASTQ files for single-cell RNA-Seq datasets. For each read the cell barcode and UMI information and the equivalence class resulting from pseudoalignment are stored in a [BUS](https://github.com/BUStools/BUS) file `output.bus` stored in the output directory directory, along with `matrix.ec` and `transcripts.txt` which store information about the equivalence classes and transcript names for downstream processing.

~~~
kallisto 0.50.0
Generates BUS files for single-cell sequencing

Usage: kallisto bus [arguments] FASTQ-files

Required arguments:
-i, --index=STRING            Filename for the kallisto index to be used for
                              pseudoalignment
-o, --output-dir=STRING       Directory to write output to
-x, --technology=STRING       Single-cell technology used 

Optional arguments:
-l, --list                    List all single-cell technologies supported
-B, --batch=FILE              Process files listed in FILE
-t, --threads=INT             Number of threads to use (default: 1)
-b, --bam                     Input file is a BAM file
-n, --num                     Output number of read in flag column (incompatible with --bam)
-N, --numReads                Maximum number of reads to process from supplied input
-T, --tag=STRING              5′ tag sequence to identify UMI reads for certain technologies
    --fr-stranded             Strand specific reads for UMI-tagged reads, first read forward
    --rf-stranded             Strand specific reads for UMI-tagged reads, first read reverse
    --unstranded              Treat all read as non-strand-specific
    --paired                  Treat reads as paired
    --long                    Treat reads as long
    --aa                      Align to index generated from a FASTA-file containing amino acid sequences
    --inleaved                Specifies that input is an interleaved FASTQ file
    --batch-barcodes          Records both batch and extracted barcode in BUS file
    --verbose                 Print out progress information every 1M proccessed reads
~~~

To process the `output.bus` file further use [bustools](https://github.com/BUStools/bustools); examples of downstream processing can be seen in [dataset specific notebooks](https://github.com/BUStools/Notebooks) available from the [bustools repository](https://github.com/BUStools). 

Running `kallisto bus -l` gives a list of currently supported single cell technologies

~~~
List of supported single-cell technologies

short name       description
----------       -----------
10xv1            10x version 1 chemistry
10xv2            10x version 2 chemistry
10xv3            10x version 3 chemistry
Bulk             Bulk RNA-seq
SmartSeq2        Smart-seq2
BDWTA            BD Rhapsody WTA
CELSeq           CEL-Seq
CELSeq2          CEL-Seq version 2
DropSeq          DropSeq
inDropsv1        inDrops version 1 chemistry
inDropsv2        inDrops version 2 chemistry
inDropsv3        inDrops version 3 chemistry
SCRBSeq          SCRB-Seq
SmartSeq3        Smart-seq3
SPLiT-seq        SPLiT-seq
SureCell         SureCell for ddSEQ
Visium           10x Visium Spatial Transcriptomics
~~~

When specifying the input the short name can be used to indicate the technology.

Additionally `kallisto bus` will accept a string specifying a new technology in the format of `bc:umi:seq` where each of `bc`,`umi` and `seq` are a triplet of integers separated by a comma, denoting the file index, start and stop of the sequence used. For example to specify the `10xV2` technology we would use `0,0,16:0,16,26:1,0,0`. The first part `bc` is `0,0,16` indicating it is in the 0-th file (also known as the first file in plain english), the barcode starts at the 0-th bp and ends at the 16-th bp in the sequence (i.e. 16bp barcode), the UMI is similarly in the same file, right after the barcode in position 16-26 (a 10bp UMI), finally the sequence is in a separate file, starts at 0 and ends at 0 (in this case stopping at 0 means there is no limit, we use the entire sequence).

##### Batch files

Using the --batch option in `kallisto bus ` allows you to supply a list of fastq files in a separate text file (each line represents a unique sample or cell).

The format of the batch file is

~~~
#id file1 file 2
cell1 cell1_1.fastq.gz cell1_1.fastq.gz
cell2 cell2_1.fastq.gz cell2_1.fastq.gz
cell3 cell3_1.fastq.gz cell3_1.fastq.gz
...
~~~

#### h5dump

`kallisto h5dump` converts
[HDF5](https://www.hdfgroup.org/HDF5/whatishdf5.html)-formatted results to
plaintext. The arguments for the h5dump command are:

~~~
kallisto 0.46.0
Converts HDF5-formatted results to plaintext

Usage:  kallisto h5dump [arguments] abundance.h5

Required argument:
-o, --output-dir=STRING       Directory to write output to

~~~

#### inspect

`kallisto inspect` can output the Target de Bruijn Graph in the index in two ways, as a file in `GFA` [format](https://github.com/GFA-spec/GFA-spec) or it can map the contigs of the graph and and equivalence classes in a `BED` format that can be visualized using [IGV](http://software.broadinstitute.org/software/igv/)

~~~
kallisto 0.50.0

Usage: kallisto inspect INDEX-file

Optional arguments:
-t, --threads=INT       Number of threads
~~~

#### version

`kallisto version` displays the current version of the software.

#### cite

`kallisto cite` displays the citation for the paper.
