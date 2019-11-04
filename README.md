# kallisto

__kallisto__ is a program for quantifying abundances of transcripts from
RNA-Seq data, or more generally of target sequences using high-throughput
sequencing reads. It is based on the novel idea of _pseudoalignment_ for
rapidly determining the compatibility of reads with targets, without the need
for alignment. On benchmarks with standard RNA-Seq data, __kallisto__ can
quantify 30 million human bulk RNA-seq reads in less than 3  minutes on a Mac desktop
computer using only the read sequences and a transcriptome index that
itself takes than 10 minutes to build. Pseudoalignment of reads
preserves the key information needed for quantification, and __kallisto__
is therefore not only fast, but also comparably accurate to other existing
quantification tools. In fact, because the pseudoalignment procedure is
robust to errors in the reads, in many benchmarks __kallisto__
significantly outperforms existing tools. The __kallisto__ algorithms are described in more detail in:

NL Bray, H Pimentel, P Melsted and L Pachter, [Near optimal probabilistic RNA-seq quantification](http://www.nature.com/nbt/journal/v34/n5/abs/nbt.3519.html), Nature Biotechnology __34__, p 525--527 (2016).

Scripts reproducing all the results of the paper are available [here](https://github.com/pachterlab/kallisto_paper_analysis).

__kallisto__ quantified bulk RNA-Seq can be analyzed with [__sleuth__](https://github.com/pachterlab/sleuth/).

__kallisto__ can be used together with [__bustools__](https://bustools.github.io/) to pre-process single-cell RNA-seq data. See the [kallistobus.tools](https://www.kallistobus.tools/) website for instructions.

## Manual

Please visit http://pachterlab.github.io/kallisto/manual.html for the manual.

## License

__kallisto__ is distributed under the BSD-2 license. The license is distributed with __kallisto__ in the file `license.txt`, which is also viewable [here](https://pachterlab.github.io/kallisto/download). Please read the license before using __kallisto__.

## Getting help

For help running __kallisto__, please post to the [kallisto-and-applications Google Group](https://groups.google.com/forum/#!forum/kallisto-and-applications).

## Reporting bugs

Please report bugs to the [Github issues page](https://github.com/pachterlab/kallisto/issues)

## Development and pull requests

We typically develop on separate branches, then merge into devel once features
have been sufficiently tested. `devel` is the latest, stable, development
branch. `master` is used only for official releases and is considered to be
stable. If you submit a pull request (thanks!) please make sure to request to
merge into `devel` and NOT `master`. Merges usually only go into `master`, but
not out.
