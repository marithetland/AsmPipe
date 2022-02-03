# Asmbl
**Assembly and quality assessment of short-read Illumina data**

Pipeline of tools used for quality assessment, assembly and gene detection from short-read Illumina sequences. Created for easy use at our microbiology research lab at Stavanger University Hopsital.

This script will quality- and adapter-trim your data with trim_galore, create a FastQC report, assemble the trimmed FASTQ reads using Unicycler, create a Quast QC-report, run mlst to identify sequence types and species, and will calculate average read depth (or coverage) of each sample.

The script creates a report summarising for each sample: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and read depth.

## Table of Contents

[Requirements](#Requirements)  
[Basic Usage](#Basic-usage)  
[Usage](#Usage)  
[Output](#Output)  
[Detailed Explanation](#Detailed-explanation)  
[Updates](#Updates)  

## Requirements
These need to be installed and in path for the entire pipeline to work. Other versions of these tools will possibly work too, but these are the ones I have tested.

* Linux or MacOS
* Python 3.9.7
* Pandas (`pip3 install pandas`)
* Paralell (`conda install -c conda-forge parallel`)
* FastQC v0.11.9 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (`conda install -c bioconda fastqc`)
* MultiQC v1.11 (https://multiqc.info/) (`conda install -c bioconda multiqc`)
* CutAdapt v3.5 (for TrimGalore) (`conda install -c bioconda cutadapt`)
* TrimGalore v0.6.7 (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* Unicycler v0.5.0 (https://github.com/rrwick/Unicycler#installation) 
* SPAdes v3.15.3 (http://cab.spbu.ru/software/spades/) (`conda install -c bioconda spades=3.13.0`)
* BLAST+ v2.12.0+ (`conda install -c bioconda blast`)
* Bowtie2 v2.4.5 (`conda install -c bioconda bowtie2`)
* Quast v5.0.2 (http://quast.sourceforge.net/quast) (`conda install -c bioconda quast`)
* mlst v2.19.0 (https://github.com/tseemann/mlst) (`conda install -c conda-forge -c bioconda -c defaults mlst`)
* BWA v0.7.17-r1188 (http://bio-bwa.sourceforge.net/) (`conda install -c bioconda bqa`)
* SAMtools v1.14 (http://www.htslib.org/download/) (`conda install -c bioconda samtools `)
* PicardTools v2.18.29-0 (https://broadinstitute.github.io/picard/) (`conda install -c bioconda picard`)
* Optional: Kleborate v2.20 (https://github.com/katholt/Kleborate) including Kaptive v2.0.0
* Optional but recommended: Install all in a conda environment

## Basic usage

You must be in the directory containing the FASTQ-files to run this pipeline. Output-files will be stored in a specific file-structure in the input-directory.

``` 
cd ~/Directory_with_fastq/ #Enter directory with FASTQ-files
asmbl.py 
```

## Usage

You must be in the directory containing the FASTQ-files to run this pipeline. Output-files will be stored in a specific file-structure in the input-directory. In addition to the default pipeline, you can also run kleborate or abricate.

Usage:

```
ASMBL [-h] [-v] [-t THREADS] [--noex] [--nofqc] [--nomlst]
               [--noquast] [--nocov] [--klebs] [--argannot] [--resfinder]
               [--plasmidfinder] [--card] [--ncbi] [--ecoh] [--abricate_all]

ASMBL

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -t THREADS, --threads THREADS
                        Specify number of threads to use. Default: 4
  --noex                Do not run fastQC, multiQC, Quast, MLST or
                        read depth calculation.
  --nofqc               Do not run fastQC and multiQC
  --nomlst              Do not run MLST
  --noquast             Do not run Quast
  --nocov               Do not calculate read depth (X)
  --klebs               Run Kleborate, with option --all

```

You can add more FASTQ-files to the same output-directories/summary-report, by adding the FASTQ-files to the inital input-directory, leaving the output-folder structure as it was created and re-running the pipeline.

Note: This was initially created for scientist with little/no coding-experience to easily perform assembly, therefore, this script currently only works when you run it from the folder the FASTQ-files are in, and output-files are stored in the same directory. In the future, I will add input and output-options.
 

## Output

The following output-files are created when running AsmPipe:

* Fastq_raw: Any FASTQ-files that have been processed are placed here
* trimmed_reads: The trimmed reads from TrimGalore are placed here
* assembly: The Unicycler-assembled files will be placed here
* assemblies: The FASTA-file from assembly will be copied to this direcory
* QC: Contains reports from FastQC, multiQC, Quast, trimming and overall coverage calculation
* analyses: Contains results from mlst, kleborate and abricate
* sequence_list.txt: List of all samples that have been analysed
* successful_sequences.txt: List of all samples that were successfully assembled
* failed_sequences.txt: List of any samples that failed any stage of the pipeline
* logs: Will contain run-logs from each tool for each sample
* AsmPipe_date_time.csv: Overall summary report of: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and sequence depth.


## Detailed explanation

* FastQC performs quality assessment of raw reads, indicating number of reads, GC%, adapter content, sequence length distribution, and more
* TrimGalore - trims raw reads based on adapter sequences and Phred quality: trims 1 bp off 3' end of every read, removes low-quality (<Phred 20) 3' ends, removes adapter sequences and removes read-pairs if either of the reads' length is <20 bp
* Unicycler functions as a SPAdes optimiser with short-reads only, and pilon polishing attempts to make imporvements on the genome
* Quast quality assessment on assembly outputs the total length, GC%, number of contigs, N50, L50 and more. 
* MLST attempts to identify species and mlst based on the PubMLST schemes. Other tools may be needed for specification, e.g. Kleborate identifies locus variants for Klebsiella samples and separates klebsiella pneumoniae sensu lato into subspecies
* Sequencing depth (X) - maps the reads against their assembled fasta-file to calculate the overall average depth of the genome.

Things to check QC-wise
* That GC% matches the sample species
* That the total length matches the sample species
* That you do not have a high number of contigs (ideally <500)
* That you do not have low average read depth (ideally >40X)
* A low number og long contigs is preferable to a high number of contigs with short contigs

## License
[GNU General Public License, v3](https://www.gnu.org/licenses/gpl-3.0.html)

## Updates
2022-02-03: Updated pipeline to work with newest release of Unicycler (v0.5.0). Unicycler no longer uses pilon for polishing and read error correction is by default turned off, so these options have been removed from the pipeline. Also removed option for abricate as it does not work. Also updated some terms to stop confusion: The folder "assemblies" is now "fasta", "coverage" is now "read depth", and the final report from the pipeline is prefixed with "Asmbl" rather than "AsmPipe".

2021-04-21: Added version checks, fixed a bug with threading and most importantly: Added the flag `--no_correct` to the unicycler command to turn off spades read error correction. This is not needed as the files are QC'd with trim-galore first.

2019-07-18: Added options to find \*fastq-gz files in subdirectories from previuos runs

2019-07-18: Added options to not run parts of the pipeline, and added option to run kleborate (https://github.com/katholt/Kleborate) at the end of the pipeline

2019-10-29: Added options to run ABRICATE as part of the pipeline, and created output-folder "analyses" to put mlst, kleborate and ABRICATE-outputs in. Also added merge_runs.sh which you can use to merge two parent folders with the same structure (from this script). 

2019-10-30: Updated tool versions in README.

