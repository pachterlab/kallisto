---
layout: page
title: "Manual"
description: ""
group: navigation
---
{% include JB/setup %}

Typing `kallisto` produces a list of usage options, which are:

    Usage: kallisto <CMD> [arguments] ..

    Where <CMD> can be one of:

        index         Builds a kallisto index 
        quant         Runs the quantification algorithm 
        h5dump        Converts HDF5-formatted results to plaintext
        version       Prints version information

    Running kallisto <CMD> without arguments prints usage information for <CMD>

The four usage commands are:

#### index

`kallisto index` builds an index from a Fasta formatted file of target sequences. The arguments for the index command are:

    Usage: kallisto index [arguments] FASTA-file

    Required argument:
    -i, --index=STRING          Filename for the kallisto index to be constructed 
     
    Optional argument:
    -k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)

The Fasta file supplied can be either in plaintext or gzipped format. 

#### quant

`kallisto quant` runs the quantification algorithm. The arguments for the quant command are:

    Usage: kallisto quant [arguments] FASTQ-files
    
    Required arguments:
    -i, --index=STRING            Filename for the kallisto index to be used for
                                  quantification
    -o, --output-dir=STRING       Directory to write output to
    
    Optional arguments:
    -l, --fragment-length=DOUBLE  Estimated average fragment length
                                  (default: value is estimated from the input data)
    -b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)
        --seed=INT                Seed for the bootstrap sampling (default: 42)
        --plaintext               Output plaintext instead of HDF5

__kallisto__ can process either single-end or paired-end reads. The running mode (single-end vs. paired-end) is determined by the number of FASTQ files provided (one vs. two respectively). 

In the case of single-end reads, the -l option must be used to specify the average framgent length. Typical Illumina libraries produce fragment lengths ranging from 180--200 bp but it's best to determine this from a library quantification with an instrument such as an Agilent Bioanalyzer. For paired-end reads, the average fragment length can be directly estimated from the reads and the program will do so if -l is not used (this is the preferred run mode).

The number of bootstrap samples is specified using -b. Note that because of the large amount of data that may be produced when the number of bootstrap samples is high, __kallisto__ outputs bootstrap results in HDF5 format. The `h5dump` command can be used afterwards to convert this output to plaintext, however most convenient is to analyze bootstrap results with __sleuth__.

#### h5dump

`kallisto h5dump` converts [HDF5](https://www.hdfgroup.org/HDF5/whatishdf5.html)-formatted results to plaintext. The arguments for the h5dump command are:

    Usage:  kallisto h5dump [arguments] abundance.h5

    Required argument:
    -o, --output-dir=STRING       Directory to write output to

#### version

`kallisto version` displays the current version of the software.