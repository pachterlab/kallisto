---
layout: page
title: "Getting Started"
description: ""
group: navigation
---
{% include JB/setup %}

The short tutorial below explains how to run __kallisto__ using a small example distributed with the program. 


#### Download and installation

Begin by downloading and installing the program by following instructions on the [download page](https://pachterlab.github.io/kallisto/download). The files needed to confirm that __kallisto__ is working are included with the binaries downloadable from the [download page](https://pachterlab.github.io/kallisto/download).

After downloading and installing __kallisto__ you should be able to type `kallisto` and see:

~~~
kallisto 0.44.0

Usage: kallisto <CMD> [arguments] ..

Where <CMD> can be one of:

    index         Builds a kallisto index 
    quant         Runs the quantification algorithm 
    pseudo        Runs the pseudoalignment step 
    h5dump        Converts HDF5-formatted results to plaintext
    inspect       Inspects and gives information about an index
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

The results of the main quantification, i.e. the abundance estimate using __kallisto__ on the data is in the `abundance.tsv` file. Abundances are reported in "estimated counts" (est_counts) and in [Transcripts Per Million](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/) (TPM). The abundance.tsv file you get should look like this:

~~~
target_id	length	eff_length	est_counts	tpm
ENST00000513300.5	1924	1746.98	102.328	11129.2
ENST00000282507.7	2355	2177.98	1592.02	138884
ENST00000504685.5	1476	1298.98	68.6528	10041.8
ENST00000243108.4	1733	1555.98	343.499	41944.9
ENST00000303450.4	1516	1338.98	664	94221.8
ENST00000243082.4	2039	1861.98	55	5612.36
ENST00000303406.4	1524	1346.98	304.189	42908.2
ENST00000303460.4	1936	1758.98	47	5076.85
ENST00000243056.4	2423	2245.98	42	3553.05
ENST00000312492.2	1805	1627.98	228	26609.9
ENST00000040584.5	1889	1711.98	4295	476675
ENST00000430889.2	1666	1488.98	623.628	79578.2
ENST00000394331.3	2943	2765.98	85.6842	5885.85
ENST00000243103.3	3335	3157.98	962	57879.3
~~~

The file is tab delimited so that it can easily parsed. The output can also be analyzed with the __sleuth__ tool.

 The `run_info.json` file contains a summary of the run, including data on the number targets used for quantification, the number of bootstraps performed, the version of the program used and how it was called. You should see this:

~~~
{
	"n_targets": 14,
	"n_bootstraps": 30,
	"n_processed": 10000,
	"n_pseudoaligned": 9413,
	"n_unique": 7174,
	"p_pseudoaligned": 94.1,
	"p_unique": 71.7,
	"kallisto_version": "0.44.0",
	"index_version": 10,
	"start_time": "Tue Jan 30 09:34:31 2018",
	"call": "kallisto quant -i transcripts.kidx -b 30 -o kallisto_out reads_1.fastq.gz reads_2.fastq.gz"
}
~~~

 The h5 file contains the main quantification together with the boostraps in [HDF5 format](https://www.hdfgroup.org/HDF5/whatishdf5.html). The reason for this binary format is to compress the large output of runs with many bootstraps. The __h5dump__ command in __kallisto__ can be used to convert the file to plain-text.


To visualize the pseudoalignments we need to run __kallisto__ with the `--genomebam` option. To do this we need two additional files, a GTF file, which describes where the transcripts lie in the genome, and a text file containing the length of each chromosome. These files are part of the test directory. To run __kallisto__ we type

~~~
kallisto quant -i transcripts.kidx -b 30 -o kallisto_out --genomebam --gtf transcripts.gtf.gz --chromosomes chrom.txt reads_1.fastq.gz reads_2.fastq.gz
~~~

this is the same run as above, but now we supply `--gtf transcripts.gtf.gz` for the GTF file and the chromoeme file `--chromosomes chrom.txt`. For a larger transcriptome we recommend downloading the GTF file from the same release and data source as the FASTA file used to construct the index. The output now contains two additional files `pseudoalignments.bam` and `pseudoalignments.bam.bai`. The files can be viewed and processed using [Samtools](http://www.htslib.org/) or a genome browser such as [IGV](http://software.broadinstitute.org/software/igv/). There is no need to sort or index the BAM file since kallisto does that directly. For windows users we recommend using the IGV browser, since there are no native Samtools releases (except using Linux Subsystem on Windows 10).

That's it.

You can now run kallisto on your dataset of choice. For convenience, we have placed some transcriptome fasta files for human and model organisms [here](https://github.com/pachterlab/kallisto-transcriptome-indices/releases). Publicly available RNA-Seq data can be found on the [short read archive](http://www.ncbi.nlm.nih.gov/sra) (a convenient mirror and interface to the SRA is available [here](http://sra.dnanexus.com)). While __kallisto__ cannot process .sra files, such files can be converted to FASTQ with the [fastq-dump](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) tool which is part of the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).
