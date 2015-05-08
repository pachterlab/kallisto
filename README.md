# kallisto

[![Join the chat at https://gitter.im/pachterlab/kallisto](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/pachterlab/kallisto?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

__kallisto__ is a program for quantifying abundances of transcripts from
RNA-Seq data, or more generally of target sequences using high-throughput
sequencing reads. It is based on the novel idea of _pseudoalignment_ for
rapidly determining the compatibility of reads with targets, without the need
for alignment. On benchmarks with standard RNA-Seq data, __kallisto__ can
quantify 30 million human reads in less than 3  minutes on a Mac desktop
computer using only the read sequences and a transcriptome index that
itself takes less than 10 minutes to build. Pseudoalignment of reads
preserves the key information needed for quantification, and __kallisto__
is therefore not only fast, but also as accurate than existing
quantification tools. In fact, because the pseudoalignment procedure is
robust to errors in the reads, in many benchmarks __kallisto__
significantly outperforms existing tools.

## Manual

Please visit http://pachterlab.github.io/kallisto for the manual.

## License

Please read the license before using kallisto. The license is distributed with __kallisto__ in the file `license.txt` also viewable [here](http://pachterlab.github.io/kallisto). 
