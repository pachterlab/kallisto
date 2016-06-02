---
layout: page
title: "Getting Started"
description: ""
group: navigation
---
{% include JB/setup %}

The short tutorial below explains how to run __kallisto__ using a small example distributed with the program. If you already know how to use __kallisto__ but just want to get it quickly installed on your Mac, the easiest way to do so is using brew. Type:

`ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"` to install brew, followed by

`brew tap homebrew/science` and then

`brew install kallisto`

to have kallisto universally executable.


#### Download

The example/test distributed with __kallisto__ is included with the binaries we distribute with the program:

- For Mac download the __kallisto__ binary from [here](download.html)

- For Linux, please [download the source code](https://github.com/pachterlab/kallisto/releases) and compile.

#### Install

**\*\*Note\*\*** For users who do not have 'root' access, please follow the [local build
tutorial](local_build.html) or [download](download.html) directly.

If you do not already have __kallisto__ universally executable on your machine, Begin by copying the __kallisto__ executable from the downloaded binary to

`cp kallisto /usr/local/bin/`

Depending on your system configuration, you might need to prepend the command
with `sudo`.

You should be able to type `kallisto` and see:

~~~
kallisto 0.43.0

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

#### Building an index

__kallisto__ quantifies read files directly without the need for read alignment, but it does perform a procedure called pseudoalignment. Pseudoalignment requires processing a transcriptome file to create a "transcriptome index". To begin, first change directories to where the test files distributed with the __kallisto__ executable are located:

`cd kallisto/tests`

Next, build an index type:

`kallisto index -i transcripts.idx transcripts.fasta.gz`

#### Quantification

Now you can quantify abundances of the transcripts using the two read files reads_1.fastq.gz and reads_2.fastq.gz (the .gz suffix means the read files have been gzipped; __kallisto__ can read in either plain-text or gzipped read files). To quantify abundances type:

`kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz`

You can also call __kallisto__ with

`kallisto quant -i transcripts.idx -o output -b 100 <(gzcat reads_1.fastq.gz) <(gzcat reads_2.fastq.gz)`

or with linux, you replace `gzcat` with `zcat` or any other program that writes the FASTQ to stdout. This utilizes an additional core to uncompress the FASTQ files, and speeds up the program by 10--15%.

#### Single end reads

If your reads are single end only you can run kallisto by specifying the `--single` flag,

`kallisto quant -i transcripts.idx -o output -b 100 --single -l 180 -s 20 reads_1.fastq.gz`

however you must supply the length and standard deviation of the fragment length (not the read length).

#### Results

The results of a __kallisto__ run are placed in the specified output directory (the -o option), and therefore the test results should be located in the subdirectory "output". The contents of the directory should look like this:

    total 568
    -rw-r--r--  1 username  staff  282480 May  3 10:10 abundance.h5
    -rw-r--r--  1 username  staff     589 May  3 10:10 abundance.tsv
    -rw-r--r--  1 username  staff     227 May  3 10:10 run_info.json

The results of the main quantification, i.e. the abundance estimate using __kallisto__ on the data is in the `abundance.txt` file. Abundances are reported in "estimated counts" (est_counts) and in [Transcripts Per Million](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/) (TPM). The abundance.txt file you get should look like this:

    target_id       length  eff_length      est_counts      tpm
    NM_001168316    2283    2105.9  164.133 12856.9
    NM_174914       2385    2207.9  1495.6  111741
    NR_031764       1853    1675.9  104.27  10263.4
    NM_004503       1681    1503.9  332.001 36416.5
    NM_006897       1541    1363.9  664     80308.9
    NM_014212       2037    1859.9  55      4878.11
    NM_014620       2300    2122.9  592.584 46046.7
    NM_017409       1959    1781.9  47      4351.04
    NM_017410       2396    2218.9  42      3122.41
    NM_018953       1612    1434.9  227.995 26210.8
    NM_022658       2288    2110.9  4881    381434
    NM_153633       1666    1488.9  359.898 39874.2
    NM_153693       2072    1894.9  72.5147 6312.74
    NM_173860       849     671.903 962     236182
    NR_003084       1640    1462.9  0.00787013      0.887453

The file is tab delimited so that it can easily parsed. The output can also be analyzed with the __sleuth__ tool.

 The `run_info.json` file contains a summary of the run, including data on the number targets used for quantification, the number of bootstraps performed, the version of the program used and how it was called. You should see this:


~~~
{
	"n_targets": 15,
	"n_bootstraps": 0,
	"n_processed": 10000,
	"kallisto_version": "0.42.5",
	"index_version": 10,
	"start_time": "Thu Jun  2 17:45:42 2016",
	"call": "kallisto quant -i transcripts.idx -o output -b 100 reads_1.fastq.gz reads_2.fastq.gz"
}
~~~


 The h5 file contains the main quantification together with the boostraps in [HDF5 format](https://www.hdfgroup.org/HDF5/whatishdf5.html). The reason for this binary format is to compress the large output of runs with many bootstraps. The __h5dump__ command in __kallisto__ can be used to convert the file to plain-text.

That's it.

You can now run kallisto on your dataset of choice. For convenience, we have placed some transcriptome fasta files for human and model organisms [here](http://bio.math.berkeley.edu/kallisto/transcriptomes/). Publicly available RNA-Seq data can be found on the [short read archive](http://www.ncbi.nlm.nih.gov/sra) (a convenient mirror and interface to the SRA is available [here](http://sra.dnanexus.com)). While __kallisto__ cannot process .sra files, such files can be converted to FASTQ with the [fastq-dump](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) tool which is part of the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).
