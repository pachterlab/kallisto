---
layout: page
title: "FAQ"
group: navigation
---

{% include JB/setup %}

- I'm having trouble with __kallisto__. Can I get help?

Yes. If you have think you have discovered a bug that needs to be fixed please file a report on the GitHub page. If you have a question about installing or running the program please ask on the [Google user group](https://groups.google.com/forum/#!forum/kallisto-sleuth-users).

- __kallisto__ is fast but is it as accurate as other quantification tools? 

Yes. The source of  __kallisto__'s speed is the novel concept of _pseudoalignment_ of reads. Pseudo alignment is all that is needed for accurate quantification, and in tests on both real and simulated data we find that __kallisto__ is frequently _more_ accurate than existing methods. The increase in accuracy over some methods is due to the robustness of pseudo alignment to read base errors and indels, and the ability to perform many rounds of the Expectation-Maximization algorithm.

- Can __kallisto__  be used for finding differentially expressed genes in RNA-Seq analysis? 

No. __kallisto__ can be used to build indices for quantification, and it can quantify RNA-Seq samples. It is also useful for assessing the reliability of abundance estimates via the bootstrap samples it produces. But it is not a tool for performing differential analysis. However the companion tool __sleuth__ is designed to for that purpose, and __sleuth__  uses the boostrapped samples produced by __kallisto__.

- Is __kallisto__ usable with both single-end and paired-end reads?

Yes.

- Does __kallisto__ require reads to be of the same length?

No.

- Can __kallisto__ be used to quantify single-cell RNA-Seq data?

Yes.

- I've already mapped all my reads. Can I use those mapping with __kallisto__?

No. The mappings are not relevant and not needed for __kallisto__. 


- Are you distributing pre-built indices?

No. Building indices with __kallisto index__ will be faster in practice than downloading index files. For example, the __kallisto__ index for the human transcriptome takes between 5--10 minutes to build on a standard desktop or laptop. We are distributing transcriptome fasta files for model organisms [here](http://bio.math.berkeley.edu/kallisto/transcriptomes/).

- My RNA-Seq was prepared with a stranded library. Is there a special option I need to use with kallisto?

No.

- Do you know how I should pick the k-mer length to use?

Yes. There is a technical condition on the target de Bruijn graph that must be met, namely that the girth should be greater than k (we will soon provide a tool for this calculation with __kallisto__ but it is not yet ready). Subject to satisfying that girth condition, it is best to choose the smallest possible k-mer length as that increases robustness to error. Larger k-mer lengths lead to faster run times but since that is not an issue with kallisto.

- Do I need to unzip my gzipped Fasta or FastQ files before running kallisto?

No.

- Is there a reason you picked the name kallisto for your program?

Yes.



