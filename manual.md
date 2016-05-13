---
layout: page
title: "Manual"
description: ""
group: navigation
---
{% include JB/setup %}

Typing `kallisto` produces a list of usage options, which are:

~~~
kallisto 0.42.5

Usage: kallisto <CMD> [arguments] ..

Where <CMD> can be one of:

    index         Builds a kallisto index
    quant         Runs the quantification algorithm
    pseudo        Runs the pseudoalignment step
    h5dump        Converts HDF5-formatted results to plaintext
    version       Prints version information
    cite          Prints citation information

Running kallisto <CMD> without arguments prints usage information for <CMD>

~~~
The usage commands are:

#### index

`kallisto index` builds an index from a FASTA formatted file of target sequences. The arguments for the index command are:

~~~
kallisto 0.42.5
Builds a kallisto index

Usage: kallisto index [arguments] FASTA-files

Required argument:
-i, --index=STRING          Filename for the kallisto index to be constructed

Optional argument:
-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)
    --make-unique           Replace repeated target names with unique names
~~~

The Fasta file supplied can be either in plaintext or gzipped format.

#### quant

`kallisto quant` runs the quantification algorithm. The arguments for the quant command are:

~~~
kallisto 0.42.5
Computes equivalence classes for reads and quantifies abundances

Usage: kallisto quant [arguments] FASTQ-files

Required arguments:
-i, --index=STRING            Filename for the kallisto index to be used for
                              quantification
-o, --output-dir=STRING       Directory to write output to

Optional arguments:
    --bias                    Perform sequence based bias correction
-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
    --seed=INT                Seed for the bootstrap sampling (default: 42)
    --plaintext               Output plaintext instead of HDF5
    --single                  Quantify single-end reads
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
                              (default: value is estimated from the input data)
-t, --threads=INT             Number of threads to use (default: 1)
    --pseudobam               Output pseudoalignments in SAM format to stdout
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

In the case of single-end reads, the -l option must be used to specify the average fragment length. Typical Illumina libraries produce fragment lengths ranging from 180--200 bp but it's best to determine this from a library quantification with an instrument such as an Agilent Bioanalyzer. For paired-end reads, the average fragment length can be directly estimated from the reads and the program will do so if -l is not used (this is the preferred run mode).

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

+ `--bias` learns parameters for a model of sequences specific bias and corrects the abundances accordlingly.

+ `-t, --threads` specifies the number of threads to be used both for pseudoalignment and running bootstrap. The default value is 1 thread, specifying more than the number of bootstraps or the number of cores on your machine has no additional effect.

##### Pseudobam

`--pseudobam` outputs all pseudoalignments in SAM format to the standard output. The stream can either be redirected into a file, or converted to bam using `samtools`.

For example

~~~
kallisto quant -i index -o out --pseudobam r1.fastq r2.fastq > out.sam
~~~

or by piping directly into `samtools`

~~~
kallisto quant -i index -o out --pseudobam r1.fastq r2.fastq | samtools view -Sb - > out.bam
~~~

A detailed description of the SAM output is [here](pseudobam.html).

#### pseudo

`kallisto pseudo` runs only the pseudoalignment step and is meant for usage in single cell RNA-seq. The arguments for the pseudo command are:

~~~
kallisto 0.42.5
Computes equivalence classes for reads and quantifies abundances

Usage: kallisto pseudo [arguments] FASTQ-files

Required arguments:
-i, --index=STRING            Filename for the kallisto index to be used for
                              pseudoalignment
-o, --output-dir=STRING       Directory to write output to

Optional arguments:
-b  --batch=FILE              Process files listed in FILE
    --single                  Quantify single-end reads
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
                              (default: value is estimated from the input data)
-t, --threads=INT             Number of threads to use (default: 1)
    --pseudobam               Output pseudoalignments in SAM format to stdout

~~~

The form of the command and the meaning of the parameters are identical to the quant command. However, pseudo does not run the EM-algorithm to quantify abundances. In addition the pseudo command has an option to specify many cells in a batch file, e.g.

~~~
kallisto pseudo -i index -o output -b batch.txt
~~~

which will read information about each cell in the `batch.txt` file and process all cells simultaneously.

The format of the batch file is

~~~
#id file1 file 2
cell1 cell1_1.fastq.gz cell1_1.fastq.gz
cell2 cell2_1.fastq.gz cell2_1.fastq.gz
cell3 cell3_1.fastq.gz cell3_1.fastq.gz
...
~~~

where the first column is the id of the cell and the next two fields are the corresponding files containing the paired end reads. Any lines starting with `#` are ignored. In the case of single end reads, specified with `--single`, only one file should be specified per cell.


#### h5dump

`kallisto h5dump` converts
[HDF5](https://www.hdfgroup.org/HDF5/whatishdf5.html)-formatted results to
plaintext. The arguments for the h5dump command are:

~~~
kallisto 0.42.5
Converts HDF5-formatted results to plaintext

Usage:  kallisto h5dump [arguments] abundance.h5

Required argument:
-o, --output-dir=STRING       Directory to write output to

~~~

#### version

`kallisto version` displays the current version of the software.

#### cite

`kallisto cite` displays the citation for the paper.
