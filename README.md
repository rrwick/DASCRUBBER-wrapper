# DASCRUBBER wrapper

Gene Myers produced a nice set of tools for scrubbing long reads: trimming them, breaking chimeras and patching up low quality regions. While raw long reads can contain a fair bit of junk, scrubbed reads are hopefully all contiguous pieces of the underlying sequence, which makes assembly much easier. Read all about it on his blog: [Scrubbing Reads for Better Assembly](https://dazzlerblog.wordpress.com/2017/04/22/1344/)

Unfortunately, there are two issues with DASCRUBBER. First, the method only works for PacBio reads, because PacBio-style FASTA headers are needed to build a Dazzler database. Second, it is not simple to run, involving more than 10 separate tools and commands.

I wrote this wrapper script to solve these issues. It carries out the entire DASCRUBBER process with a single, easy to run command, and it works on any set of long reads (including Oxford Nanopore reads) by faking PacBio read names.

Disclaimer: While nothing about this wrapper is specific to small genomes, I've never tried it on big eukaryotic genomes. If you try it on a large genome and run into problems, please let me know on the [issues page](https://github.com/rrwick/DASCRUBBER-wrapper/issues).



## Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Example commands](#example-commands)
* [Method](#method)
* [Full usage](#full-usage)
* [License](#license)



## Requirements

You'll need a number of Dazzler tools available in your PATH: [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB), [DALIGNER](https://github.com/thegenemyers/DALIGNER), [DAMASKER](https://github.com/thegenemyers/DAMASKER) and [DASCRUBBER](https://github.com/thegenemyers/DASCRUBBER)

This bash loop should clone and build them. Just replace `~/.local/bin` with wherever you want to put the executables:

```
for repo in DAZZ_DB DALIGNER DAMASKER DASCRUBBER; do
    git clone https://github.com/thegenemyers/"$repo"
    cd "$repo" && make -j && cd ..
    find "$repo" -maxdepth 1 -type f -executable -exec cp {} ~/.local/bin \;
done
```



## Installation

This tool is a single Python 3 script with no third-party dependencies. It will run without any installation:
```
git clone https://github.com/rrwick/DASCRUBBER-wrapper
DASCRUBBER-wrapper/dascrubber_wrapper.py --help
```

If you want, you can copy the it to somewhere in your PATH for easy usage:
```
cp DASCRUBBER-wrapper/dascrubber_wrapper.py ~/.local/bin
dascrubber_wrapper.py --help
```



## Example commands

__Default parameters for a 5.5 Mbp bacterial genome:__<br>
`dascrubber_wrapper.py -i reads.fastq.gz -g 5.5M | gzip > scrubbed.fasta.gz`

__Limit daligner memory usage:__<br>
`dascrubber_wrapper.py -i reads.fastq.gz -g 5.5M --daligner_options="-M80" | gzip > scrubbed.fasta.gz`

__Keep Dazzler files after completion:__<br>
`dascrubber_wrapper.py -i reads.fastq.gz -g 5.5M -d working_files --keep | gzip > scrubbed.fasta.gz`



## Method

1. Convert reads to a FASTA file with PacBio-style headers.
    * It is assumed that the input reads are either FASTQ or FASTA (one line per sequence). Gzipped reads are okay.
2. Build a [Dazzler database](https://dazzlerblog.wordpress.com/2016/05/21/dbs-and-dams-whats-the-difference/) of the reads with [`fasta2DB`](https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/).
3. Split the database with [`DBsplit`](https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/).
    * Splitting seems to be necessary or else the `DASedit` command will fail.
4. [Find all read-read overlaps](https://dazzlerblog.wordpress.com/2014/07/10/dalign-fast-and-sensitive-detection-of-all-pairwise-local-alignments/) with [`daligner`](https://dazzlerblog.wordpress.com/command-guides/daligner-command-reference-guide/).
    * This is the slowest and most memory hungry step of the process.
5. [Mask repeats](https://dazzlerblog.wordpress.com/2016/04/01/detecting-and-soft-masking-repeats/) with [`REPmask`](https://dazzlerblog.wordpress.com/command-guides/damasker-commands/).
    * The threshold depth of coverage (`REPmask`'s `-c` option) can be set in two ways. The `--repeat_depth` option will set it to a multiple of the base depth of coverage. E.g. if the base depth (as determined using the genome size) is 50x and `--repeat_depth` is 2 (the default), then regions with 100x or greater depth are considered repeats. Alternatively, you can manually set `REPmask` options: e.g. `--repmask_options="-c40"`
6. Find tandem repeats with [`datander`](https://dazzlerblog.wordpress.com/command-guides/damasker-commands/).
7. Mask tandem repeats with [`TANmask`](https://dazzlerblog.wordpress.com/command-guides/damasker-commands/).
8. Find [intrinsic quality values](https://dazzlerblog.wordpress.com/2015/11/06/intrinsic-quality-values/) with [`DASqv`](https://dazzlerblog.wordpress.com/command-guides/dascrubber-command-guide/).
9. [Trim reads and split chimeras](https://dazzlerblog.wordpress.com/2017/04/22/1344/) with [`DAStrim`](https://dazzlerblog.wordpress.com/command-guides/dascrubber-command-guide/).
10. Patch low-quality regions of reads with [`DASpatch`](https://dazzlerblog.wordpress.com/command-guides/dascrubber-command-guide/).
11. Produce a new database of scrubbed reads with [`DASedit`](https://dazzlerblog.wordpress.com/command-guides/dascrubber-command-guide/).
12. Extract a FASTA of scrubbed reads with [`DB2fasta`](https://dazzlerblog.wordpress.com/command-guides/dazz_db-command-guide/).
13. Restore original read names and output to stdout.
    * A range is appended to the end of the new read names. For example, if the original read was named `read1975`, the scrubbed read might be named `read1975/400_9198`.
    * It is possible for one original read to result in more than one scrubbed read. For example, a chimeric read named `read2392` might result in two scrubbed reads: `read2392/0_12600` and `read2392/12700_25300`.



## Full usage

```
usage: DASCRUBBER_wrapper.py -i INPUT_READS -g GENOME_SIZE [-d TEMPDIR] [-k]
                             [-r REPEAT_DEPTH]
                             [--dbsplit_options DBSPLIT_OPTIONS]
                             [--daligner_options DALIGNER_OPTIONS]
                             [--repmask_options REPMASK_OPTIONS]
                             [--datander_options DATANDER_OPTIONS]
                             [--tanmask_options TANMASK_OPTIONS]
                             [--dasqv_options DASQV_OPTIONS]
                             [--dastrim_options DASTRIM_OPTIONS]
                             [--daspatch_options DASPATCH_OPTIONS]
                             [--dasedit_options DASEDIT_OPTIONS] [-h]

A wrapper tool for the DASCRUBBER pipeline for scrubbing (trimming and chimera
removal) of long read sets (PacBio or ONT reads)

Required arguments:
  -i INPUT_READS, --input_reads INPUT_READS
                        input set of long reads to be scrubbed
  -g GENOME_SIZE, --genome_size GENOME_SIZE
                        approximate genome size (examples: 3G, 5.5M or 800k),
                        used to determine depth of coverage

Optional arguments:
  -d TEMPDIR, --tempdir TEMPDIR
                        path of directory for temporary files (default: use a
                        directory in the current location named
                        dascrubber_temp_PID where PID is the process ID)
  -k, --keep            keep the temporary directory (default: delete the
                        temporary directory after scrubbing is complete)
  -r REPEAT_DEPTH, --repeat_depth REPEAT_DEPTH
                        REPmask will be given a repeat threshold of this
                        depth, relative to the overall depth (e.g. if 2, then
                        regions with twice the base depth are considered
                        repeats) (default: 2)

Command options:
  You can specify additional options for each of the Dazzler commands if you
  do not want to use the defaults (example: --daligner_options="-M80")

  --dbsplit_options DBSPLIT_OPTIONS
  --daligner_options DALIGNER_OPTIONS
  --repmask_options REPMASK_OPTIONS
  --datander_options DATANDER_OPTIONS
  --tanmask_options TANMASK_OPTIONS
  --dasqv_options DASQV_OPTIONS
  --dastrim_options DASTRIM_OPTIONS
  --daspatch_options DASPATCH_OPTIONS
  --dasedit_options DASEDIT_OPTIONS

Help:
  -h, --help            show this help message and exit
```



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
