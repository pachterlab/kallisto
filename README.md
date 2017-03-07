# kallisto

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

__kallisto__ quantified RNA-Seq can be analyzed with [__sleuth__](https://github.com/pachterlab/sleuth/).

## Manual

Please visit http://pachterlab.github.io/kallisto/manual.html for the manual.

## License

Please read the license before using kallisto. The license is distributed with __kallisto__ in the file `license.txt` also viewable [here](https://pachterlab.github.io/kallisto/download).

## Announcements

There is a low traffic Google Group,
[kallisto-sleuth-announcements](https://groups.google.com/d/forum/kallisto-sleuth-announcements)
where we make announcements about new releases. This is a read-only mailing
list.

## Getting help

For help running __kallisto__, please post to the [kallisto-sleuth-users
Google Group](https://groups.google.com/d/forum/kallisto-sleuth-users).

## Reporting bugs

Please report bugs to the [Github issues
page](https://github.com/pachterlab/kallisto/issues)

## Development and pull requests

We typically develop on separate branches, then merge into devel once features
have been sufficiently tested. `devel` is the latest, stable, development
branch. `master` is used only for official releases and is considered to be
stable. If you submit a pull request (thanks!) please make sure to request to
merge into `devel` and NOT `master`. Merges usually only go into `master`, but
not out.
