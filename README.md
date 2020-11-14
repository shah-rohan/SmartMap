# MultiMap: Iterative Bayesian Analysis of Ambiguously Mapped Reads

This tool is designed to analyze ambiguously mapped reads from paired-end short-read next-generation sequencing by using the alignment scores and the distribution of mapped reads genome-wide to iteratively reweight and refine the weights assigned to each alignment of each read. This tool can be used for both strand-nonspecific and strand-specific applications. In principle, single-end data can also be analyzed with the MultiMap tool, but the MultiMapPrep and MultiMapRNAPrep scripts do not support single-end analysis, though modification of these scripts to use single-end data is possible if desired.

There are three components to the tool presented here: the MultiMapPrep and MultiMapRNAPrep bash scripts and the MultiMap compiled binary. The MultiMapPrep and MultiMapRNAPrep scripts serves to align the FASTQ files to the genome, filter the output file, and parse it into a BED file that can be used with the MultiMap binary, with and without strand information, respectively. The MultiMap binary is a compiled program written in C++ to analyze the alignment BED file produced by MultiMapPrep or MultiMapRNAPrep with the iterative Bayesian algorithm and produce a BEDGRAPH file of alignment weights across the genome, which can be treated as the weighted genome coverage BEDGRAPH file (with strand-specificity, if applicable).

This tool is designed for use with Bowtie2 or Hisat2.

## Documentation for MultiMapPrep

### Description

The MultiMapPrep script is used to align the FastQ files and prepare a BED file output that is ready to be processed using the MultiMap software. This is accomplished by:

1. Aligning the FastQ files to the genome using Bowtie2 using the following settings:
   * End-to-end alignment with very fast preset
   * No discordant alignments
   * No mixed alignments (single-end alignments if paired-end alignment cannot be found)
   * Report up to *k* alignments per read pair (default 51)
   * Insert size between *I* to *L* base pairs (default 100-250 bp)
2. Filtering for reads with one of the following SAM flags:
   * 99 (read paired, mapped in a proper pair, mate reverse strand, first in pair, primary alignment)
   * 163 (read paired, mapped in a proper pair, mate reverse strand, second in pair, primary alignment)
   * 355 (read paired, mapped in a proper pair, mate reverse strand, first in pair, not primary alignment)
   * 419 (read paired, mapped in a proper pair, mate reverse strand, second in pair, not primary alignment)
3. Extracting into an extended BED file:
   * Chromosome
   * Start position
   * End position
   * Read name
   * Read alignment score (AS:i:)
   * Mate alignment score (YS:i:)
4. Split the reads into separate files based on the number of alignments per read.

The output files are a Gzipped file containing all alignments and a directory `splits` with unzipped files containing extended BED files of the alignments split by number of alignments per read, prepared for line counting or use with the MultiMap software to run the iterative reweight algorithm.

### Dependencies

In addition to standard Unix tools, including awk, sed, and gzip, the MultiMapPrep script requires the following to be installed and added to the PATH environment variable.

* [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

We have tested this script on Ubuntu LTS 14.04, 16.04, and 18.04. For the above listed tools, we have tested Bowtie2 v. 2.3.4.1.

### Manual for MultiMapPrep

`MultiMapPrep [options] -x [Bowtie2 index] -o [output prefix] -1 [R1 fastq] -2 [R2 fastq]`

Inputs (required):

    -x Path to basename of Bowtie2 index for alignment
    -o Output prefix prepended to the output files 
    -1 FastQ file for read mate 1 (can be gzipped)  
    -2 FastQ file for read mate 2 (can be gzipped)

Options:

    -p Number of CPU threads to be used for multithreaded alignment (default: 1)\
    -I Minimum insert length (default: 100)\
    -L Maximum insert length (default: 250)\
    -k Maximum number of alingments reported (default: 51)\
    -s String to be removed from read names

    -h Display help message

## Documentation for MultiMapRNAPrep

### Description

The MultiMapRNAPrep script is used to align the FastQ files and prepare a BED file output that is ready to be processed using the MultiMap software. This is accomplished by:

1. Aligning the FastQ files to the genome using Hisat2 using the following settings:
   * No discordant alignments
   * No mixed alignments (single-end alignments if paired-end alignment cannot be found)
   * Report up to *k* alignments per read pair (default 51)
2. Filtering for reads with one of the following SAM flags:
   * 99 (read paired, mapped in a proper pair, mate reverse strand, first in pair, primary alignment)
   * 163 (read paired, mapped in a proper pair, mate reverse strand, second in pair, primary alignment)
   * 355 (read paired, mapped in a proper pair, mate reverse strand, first in pair, not primary alignment)
   * 419 (read paired, mapped in a proper pair, mate reverse strand, second in pair, not primary alignment)
3. Extracting into an extended BED file:
   * Chromosome
   * Start position
   * End position
   * Read name
   * Strand
   * Read alignment score (AS:i:)
   * Mate alignment score (YS:i:)
4. Split the reads into separate files based on the number of alignments per read.

The output files are a Gzipped file containing all alignments and a directory `splits` with unzipped files containing extended BED files of the alignments split by number of alignments per read, prepared for line counting or use with the MultiMap software to run the iterative reweight algorithm.

### Dependencies

In addition to standard Unix tools, including awk, sed, and gzip, the MultiMapPrep script requires the following to be installed and added to the PATH environment variable.

* [Hisat2](https://daehwankimlab.github.io/hisat2/)

We have tested this script on Ubuntu LTS 14.04, 16.04, and 18.04. For the above listed tools, we have tested Hisat2 v. 2.1.0.

### Manual for MultiMapRNAPrep

`MultiMapRNAPrep [options] -x [Hisat2 index] -o [output prefix] -1 [R1 fastq] -2 [R2 fastq]`

Inputs (required):

    -x Path to basename of Hisat2 index for alignment
    -o Output prefix prepended to the output files
    -1 FastQ file for read mate 1 (can be gzipped)
    -2 FastQ file for read mate 2 (can be gzipped)

Options:

    -p Number of CPU threads to be used for multithreaded alignment (default: 1)
    -k Maximum number of alignments reported (default: 51)
    -s String to be removed from read names

    -h Display help message

## Documentation for MultiMap

The MultiMap software serves to assign weights to each mapping of each read by an iterative Bayesian reweighting algorithm using the BED file(s) outputted by the MultiMapPrep script as an input. The output file is a Gzipped BEDGRAPH file of the genome-wide read weights analagous to a genome coverage read depth BEDGRAPH file. If run in strand-specific mode, the output file is a pair of Gzipped BEDGRAPH files of the genome-wide read weights, with one file per strand.

### Dependencies

The MultiMap binary has no dependencies in and of itself; however, building the MultiMap binary from source requires the following to be installed and added to the PATH environment variable.

* g++

Additionally, the [gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/) library is required for building the MultiMap binary; this is included in the appropriate relative directory with the appropriate reference in the makefile. 

We have tested this software on Ubuntu LTS 14.04, 16.04, and 18.04. For the above listed tools, we have tested g++ v. 7.5.0. 

### Building MultiMap binary from source

The source files for the MultiMap software are provided here. To build the MultiMap binary from source, from the main project directory, run the following command:

`make all`

The MultiMap binary will appear in the "Default" directory. It can then be moved to the directory of interest or added to the PATH environment variable for use.

### Manual for MultiMap 

`MultiMap [options] -g [genome length file] -o [output prefix] [BED or BED.gz file input(s)]`

Inputs (required):

    -g : Genome length file listing all chromosomes and lengths. Chromosomes will appear in this order in the output BEDGRAPH file.
    -o : Output prefix used for output BEDGRAPH and log files.
    Input files: Path to BED or BED.gz files, separated by spaces. Must be the last argument(s) passed to the software.

Options:

    -i : Number of iterations. Default 1.
    -x : Fixation percentage threshold. Default 0.
    -m : Maximum number of alignments for a read to be processed. Default 50.
    -s : Minimum score for Bowtie2 display. Default 0 (unscored).
    -v : N for N-fold cross-validation. Default 1 (no cross-validation).
    -c : Flag for continuous output bedgraphs. Default off.
    -S : Flag for strand-specific mode. Default off.

    -h : Display help message.

## Acknowledgments and Contact

These tools were designed and written by Rohan Shah (rohanshah@uchicago.edu).

Contact Rohan Shah (rohanshah@uchicago.edu) or Alex Ruthenburg (aruthenburg@uchicago.edu) with questions, comments, or issues.

This Software is provided under a modified MIT License. Please cite and acknowledge this repository in any works utilizing this software.
