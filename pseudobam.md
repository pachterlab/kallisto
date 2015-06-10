---
layout: page
title: "Pseudobam"
description: ""
---
{% include JB/setup %}

The output of `kallisto quant` when run with the `--pseudobam` option is a SAM formatted stream. This stream can be redirected to a SAM file or directly to a BAM file by piping into `samtools view -Sb - > output.bam`. Examples of how to do this are shown in the [Manual](manual.html).

The SAM/BAM format is described in [http://samtools.github.io/hts-specs/SAMv1.pdf](http://samtools.github.io/hts-specs/SAMv1.pdf)

`kallisto` will output a SAM file showing the pseudoalignments of the reads to the transcriptome. Each read will have one SAM entry for each transcript is pseudoaligns to.

For paired end reads the output will first include all entries for the first read, followed by all entries for the second read.

Each SAM line, where a read pseudoaligns to a transcript, follows the normal SAM format with the following interpretation of the fields



+ `QNAME`, read name, any /1 or /2 is discarded from read names and the information is stored in the FLAG field
+ `FLAG`, normal SAM flags
+ `RNAME`, the transcript name
+ `POS`, 1-based position of the leftmost mapping
+ `MAPQ`, 255 for pseudoaligned reads, 0 otherwise
+ `CIGAR`, all matches (probably a lie), e.g. for 75bp reads `75M`
+ `RNEXT`, the transcript of the mate, always equal to `RNAME`
+ `PNEXT`, the leftmost 1-based position of the mate
+ `TLEN`, observed template length
+ `SEQ`, read sequence or the reverse complement if it pseudoaligns to the reverse complement of the transcript
+ `QUAL`, the quality values of the read
+ optional fields
    + `NH`, set to the size of the equivalence class, i.e. the number of transcripts this reads pseudoaligns to.


Example output (skipping header lines)

~~~
ENST00000623973_197_729_1:0:0_2:0:0_8c  99      ENST00000623973 197     255     70M     659     =       532     CCTGGTGACCATCCTTGACTCCCGCCTGCACACGCCCAAGTACTTCTTCCTAGGGAACCTCTCCTTCCTG  2222222222222222222222222222222222222222222222222222222222222222222222  NH:i:3
ENST00000623973_197_729_1:0:0_2:0:0_8c  99      ENST00000341959 197     255     70M     659     =       532     CCTGGTGACCATCCTTGACTCCCGCCTGCACACGCCCAAGTACTTCTTCCTAGGGAACCTCTCCTTCCTG  2222222222222222222222222222222222222222222222222222222222222222222222  NH:i:3
ENST00000623973_197_729_1:0:0_2:0:0_8c  99      ENST00000625188 141     255     70M     603     =       532     CCTGGTGACCATCCTTGACTCCCGCCTGCACACGCCCAAGTACTTCTTCCTAGGGAACCTCTCCTTCCTG  2222222222222222222222222222222222222222222222222222222222222222222222  NH:i:3
ENST00000623973_197_729_1:0:0_2:0:0_8c  147     ENST00000623973 659     255     70M     197     =       -532    GTGACGAATGTGATCTTCCTAGGAGTCCCGGTTCTGTTCTTCTCTTTCTCCTATGTCTTCAGCATCACCA  2222222222222222222222222222222222222222222222222222222222222222222222  NH:i:3
ENST00000623973_197_729_1:0:0_2:0:0_8c  147     ENST00000341959 659     255     70M     197     =       -532    GTGACGAATGTGATCTTCCTAGGAGTCCCGGTTCTGTTCTTCTCTTTCTCCTATGTCTTCAGCATCACCA  2222222222222222222222222222222222222222222222222222222222222222222222  NH:i:3
ENST00000623973_197_729_1:0:0_2:0:0_8c  147     ENST00000625188 603     255     70M     141     =       -532    GTGACGAATGTGATCTTCCTAGGAGTCCCGGTTCTGTTCTTCTCTTTCTCCTATGTCTTCAGCATCACCA  2222222222222222222222222222222222222222222222222222222222222222222222  NH:i:3
~~~

In this case the read `ENST00000623973_197_729_1:0:0_2:0:0_8c`, which was simulated from `ENST00000623973` was aligned to 3 transcripts, one of which is `ENST00000623973` and the mapping positions match the read name. To decode the SAM flags consult [this handy guide](https://broadinstitute.github.io/picard/explain-flags.html).
