# MultiMap: Iterative Bayesian Analysis of Ambiguously Mapped Reads in ChIP-seq

This tool is designed to analyze ambiguously mapped reads from paired-end short-read next-generation sequencing by using the alignment scores and the distribution of mapped reads genome-wide to iteratively reweight and refine the weights assigned to each alignment of each read. This tool is primarily intended for and validated on ChIP-seq and simulated ChIP-seq data; however, in principle, this tool can be used for DNA-seq and other selected NGS applications more broadly. Furthermore, this tool was primarily devised to analyze ICeChIP-seq data; however, this tool does not itself perform ICeChIP-specific analyses and can therefore be used more broadly for ChIP-seq.

There are two components to the tool presented here: the MultiMapPrep bash script and the MultiMap compiled binary. The MultiMapPrep script serves to align the fastq files to the genome, filter the output file, and parse it into a BED file that can be used with the MultiMap binary. The MultiMap binary is a compiled program written in C++ to analyze the alignment BED file produced by MultiMapPrep with the iterative Bayesian algorithm and produce a BEDGRAPH file of alignment weights across the genome, which can be treated as the weighted genome coverage BEDGRAPH file.

This tool is designed for use with Bowtie2.

## Documentation for MultiMapPrep

This tool will take as inputs the paired-end FastQ files and the Bowtie2 index for the genome to which the reads will be aligned. The outputs will be a GZipped BED file, containing the alignments with alignment scores; and a directory containing BED files of each alignment with alignment scores, divided by the number of alignments reported fore each read.

In addition to standard Unix tools, including awk, sed, and gzip, the MultiMapPrep script requires the following to be installed and added to the PATH environment variable.

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

We have tested this script on Ubuntu LTS 14.04, 16.04, and 18.04. For the above listed tools, we have tested Bowtie2 v. 2.3.4.1.

### Manual for MultiMapPrep

MultiMapPrep \[options\] -x \[Bowtie2 index\] -o \[output prefix\] -1 \[R1 fastq\] -2 \[R2 fastq\]

Inputs (required):
-x Path to basename of Bowtie2 index for alignment
-o Output prefix prepended to the output files
-1 Fastq file for read mate 1 (can be gzipped)
-2 Fastq file for read mate 2 (can be gzipped)

Options:
-p Number of CPU threads to be used for multithreaded alignment (default: 1)
-I Minimum insert length (default: 100)
-L Maximum insert length (default: 250)
-k Maximum number of alingments reported (default: 51)
-s String to be removed from read names
