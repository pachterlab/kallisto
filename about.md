---
layout: page
title: "About"
group: navigation
---

{% include JB/setup %}

__kallisto__ is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads. It is based on the novel idea of _pseudoalignment_ for rapidly determining the compatibility of reads with targets, without the need
for alignment. On benchmarks with standard RNA-Seq data, __kallisto__ can
    quantify 30 million human reads in less than 3  minutes on a Mac desktop
    computer using only the read sequences and a transcriptome index that
    itself takes less than 10 minutes to build. Pseudoalignment of reads
    preserves the key information needed for quantification, and __kallisto__
    is therefore not only fast, but also as accurate as existing
    quantification tools. In fact, because the pseudoalignment procedure is
    robust to errors in the reads, in many benchmarks __kallisto__
    significantly outperforms existing tools. __kallisto__ is described in detail in:

Nicolas L Bray, Harold Pimentel, Páll Melsted and Lior Pachter, [Near-optimal probabilistic RNA-seq quantification](http://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html), Nature Biotechnology __34__, 525--527 (2016), doi:10.1038/nbt.3519

To use __kallisto__ [download the software](download.html) and visit the
[Getting started](starting.html) page for a quick tutorial. More information about kallisto, including a demonstration of its use, is available in the materials from the [first kallisto-sleuth workshop](https://pachterlab.github.io/kallisto-sleuth-workshop-2016/). Numerous [notebooks](https://github.com/BUStools/bustools-notebooks) provide walk-throughs demonstrating the use of kallisto for single-cell RNA-Seq processing.

