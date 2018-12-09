---
layout: page
title: "Manual"
description: ""
group: navigation
---
{% include JB/setup %}

Typing `kallisto` produces a list of usage options, which are:

~~~
kallisto 0.45.0

Usage: kallisto <CMD> [arguments] ..

Where <CMD> can be one of:

    index         Builds a kallisto index
    quant         Runs the quantification algorithm
    bus           Generate BUS files for single-cell data
    pseudo        Runs the pseudoalignment step
    merge         Merges several batch runs
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
kallisto 0.45.0
Builds a kallisto index

Usage: kallisto index [arguments] FASTA-files

Required argument:
-i, --index=STRING          Filename for the kallisto index to be constructed 

Optional argument:
-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: 31)
    --make-unique           Replace repeated target names with unique names
~~~

The Fasta file supplied can be either in plaintext or gzipped format. Prebuilt indices constructed from [Ensembl reference transcriptomes](https://uswest.ensembl.org/info/data/ftp/index.html) can be download from the [kallisto transcriptome indices](https://github.com/pachterlab/kallisto-transcriptome-indices/releases) site.

#### quant

`kallisto quant` runs the quantification algorithm. The arguments for the quant command are:

~~~
kallisto 0.45.0
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
    --fusion                  Search for fusions for Pizzly
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
    --pseudobam               Save pseudoalignments to transcriptome to BAM file
    --genomebam               Project pseudoalignments to genome sorted BAM file
-g, --gtf                     GTF file for transcriptome information
                              (required for --genomebam)
-c, --chromosomes             Tab separated file with chromosome names and lengths
                              (optional for --genomebam, but recommended)
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

- __abundances.h5__ is a HDF5 binary file containing run info, abundance
  esimates, bootstrap estimates, and transcript length information length. This
  file can be read in by __sleuth__
- __abundances.tsv__ is a plaintext file of the abundance estimates. It does
  not contains bootstrap estimates. Please use the `--plaintext` mode to output
  plaintext abundance estimates. Alternatively, `kallisto h5dump` can be used
  to output an HDF5 file to plaintext. The first line contains a header for
  each column, including [estimated counts, TPM, effective length](https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/).
- __run\_info.json__ is a json file containing information about the run


##### Optional arguments

+ `--bias` learns parameters for a model of sequences specific bias and corrects the abundances accordlingly.

+ `-t, --threads` specifies the number of threads to be used both for pseudoalignment and running bootstrap. The default value is 1 thread, specifying more than the number of bootstraps or the number of cores on your machine has no additional effect.

+ `--fr-stranded` runs kallisto in strand specific mode, only fragments where the **first** read in the pair pseudoaligns to the **forward** strand of a transcript are processed. If a fragment pseudoaligns to multiple transcripts, only the transcripts that are consistent with the first read are kept.

+ `--rf-stranded` same as `--fr-stranded` but the **first** read maps to the **reverse** strand of a transcript.

+ `--fusion` does normal quantification, but additionally looks for reads that do not pseudoalign because they are potentially from fusion genes. All output is written to the file `fusion.txt` in the output folder.

##### Pseudobam

`--pseudobam` outputs all pseudoalignments to a file `pseudoalignments.bam` in the output directory. This BAM file contains the pseudoalignments in BAM format, ordered by reads so that each pseudoalignment of a read is adjacent in the BAM file.


A detailed description of the SAM output is [here](pseudobam.html).

##### GenomeBam

`--genomebam` constructs the pseudoalignments to the transcriptome, but projects the transcript alignments to genome coordinates, resulting in split-read alignments. When the `--genomebam` option is supplied at GTF file must be given with the `--gtf` option. The GTF file, which can be plain text or gzipped, translates transcripts into genomic coordinates. We recommend downloading a the cdna FASTA files and GTF files from the same data source. The `--chromosomes` option can provide a length of the genomic chromosomes, this option is not neccessary, but gives a more consistent BAM header, some programs may require this for downstream analysis. __kallisto__ does not require the genome sequence to do pseudoalignment, but downstream tools such as genome browsers will probably need it.


#### bus

`kallisto bus` works with raw FASTQ files for single-cell RNA-Seq datasets. For each read the cell barcode and UMI information and the equivalence class resulting from pseudoalignment are stored in a [BUS](https://github.com/BUStools/BUS) file `output.bus` stored in the output directory directory, along with `matrix.ec` and `transcripts.txt` which store information about the equivalence classes and transcript names for downstream processing.

~~~
kallisto 0.45.0
Generates BUS files for single-cell sequencing

Usage: kallisto bus [arguments] FASTQ-files

Required arguments:
-i, --index=STRING            Filename for the kallisto index to be used for
                              pseudoalignment
-o, --output-dir=STRING       Directory to write output to
-x, --technology=STRING       Single-cell technology used 

Optional arguments:
-l, --list                    List all single-cell technologies supported
-t, --threads=INT             Number of threads to use (default: 1)
~~~

To process the `output.bus` file further use [bustools](https://github.com/BUStools/bustools); examples of downstream processing can be seen in [dataset specific notebooks](https://github.com/BUStools/Notebooks) available from the [bustools repository](https://github.com/BUStools). 


#### pseudo

`kallisto pseudo` runs only the pseudoalignment step and is meant for usage in single-cell RNA-seq. The arguments for the pseudo command are:

~~~
kallisto 0.45.0
Computes equivalence classes for reads and quantifies abundances

Usage: kallisto pseudo [arguments] FASTQ-files

Required arguments:
-i, --index=STRING            Filename for the kallisto index to be used for
                              pseudoalignment
-o, --output-dir=STRING       Directory to write output to

Optional arguments:
-u  --umi                     First file in pair is a UMI file
-b  --batch=FILE              Process files listed in FILE
    --single                  Quantify single-end reads
-l, --fragment-length=DOUBLE  Estimated average fragment length
-s, --sd=DOUBLE               Estimated standard deviation of fragment length
                              (default: -l, -s values are estimated from paired
                               end data, but are required when using --single)
-t, --threads=INT             Number of threads to use (default: 1)

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

When the `--umi` option is specified the batch file is of the form

~~~
#id umi-file file-1
cell1 cell_1.umi cell_1.fastq.gz
cell2 cell_2.umi cell_2.fastq.gz
cell3 cell_3.umi cell_3.fastq.gz
...
~~~

where the umi-file is a text file of the form

~~~
TTACACTGAC
CCACTCTATG
CAGGAAATCG
...
~~~

listing the Unique Molecular Identifier (UMI) for each read. The order of UMIs and reads in the fastq file must match. Even though the UMI data is single end we do not require or make use of the fragment length.

When run in **UMI** mode kallisto will use the sequenced reads to pseudoalign and find an equivalence class, but rather than count the number of reads for each equivalence class, kallisto counts the number of distinct UMIs that pseudoalign to each equivalence class.

#### h5dump

`kallisto h5dump` converts
[HDF5](https://www.hdfgroup.org/HDF5/whatishdf5.html)-formatted results to
plaintext. The arguments for the h5dump command are:

~~~
kallisto 0.45.0
Converts HDF5-formatted results to plaintext

Usage:  kallisto h5dump [arguments] abundance.h5

Required argument:
-o, --output-dir=STRING       Directory to write output to

~~~

#### merge

`kallisto merge` can merge the results of several batches performed by `pseudo`, this creates a single output as if `kallisto` had ben run on the entire sample.

~~~
kallisto 0.45.0
Computes equivalence classes for reads and quantifies abundances

Usage: kallisto merge [arguments] ouput-directories

Required arguments:
-i, --index=STRING            Filename for the kallisto index to be used for
                              pseudoalignment
-o, --output-dir=STRING       Directory to write output to
~~~

#### inspect

`kallisto inspect` can output the Target de Bruijn Graph in the index in two ways, as a file in `GFA` [format](https://github.com/GFA-spec/GFA-spec) or it can map the contigs of the graph and and equivalence classes in a `BED` format that can be visualized using [IGV](http://software.broadinstitute.org/software/igv/)

~~~
kallisto 0.45.0

Usage: kallisto inspect INDEX-file

Optional arguments:
-G, --gfa=STRING        Filename for GFA output of T-DBG
-g, --gtf=STRING        Filename for GTF file
-b, --bed=STRING        Filename for BED output (default: index + ".bed")
~~~

#### version

`kallisto version` displays the current version of the software.

#### cite

`kallisto cite` displays the citation for the paper.
