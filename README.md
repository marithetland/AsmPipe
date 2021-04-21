# Asmbl
**Assembly and quality assessment of short-read Illumina data**

Pipeline of tools used for quality assessment, assembly and gene detection from short-read Illumina sequences. Created for easy use at our microbiology research lab at Stavanger University Hopsital.

This script will quality- and adapter-trim your data, create a FastQC report, assemble the trimmed FASTQ-reads using UniCycler, create a Quast QC-report, run mlst to identify sequence types and species, and will also calculate overall coverage (sequence depth) of each sample.

The script creates a report summarising for each sample: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and sequence depth.

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
* Python 3.5.0+
* FastQC v0.11.8 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* MultiQC v1.7 (https://multiqc.info/)
* TrimGalore v0.6.4 (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* SPAdes v3.13.1 (http://cab.spbu.ru/software/spades/)
* Pilon v1.23 (https://www.broadinstitute.org/gaag/pilon)
* Unicycler v0.4.8 (https://github.com/rrwick/Unicycler)
* Quast v4.6.3 (http://quast.sourceforge.net/quast)
* mlst v2.16.3 (https://github.com/tseemann/mlst)
* BWA v0.7.17-r1188 (http://bio-bwa.sourceforge.net/)
* SAMtools v1.9 (http://www.htslib.org/download/)
* PicardTools v2.21.4-SNAPSHOT (https://broadinstitute.github.io/picard/)
* Optional: Kleborate v0.4.0-beta (https://github.com/katholt/Kleborate)
* Optional: Abricate v0.8.10 (https://github.com/tseemann/abricate)

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
                        Specify threads to use. Default: 4
  --noex                Do not run fastQC, multiQC, Quast, MLST or Coverage
                        calculation.
  --nofqc               Do not run fastQC and multiQC
  --nomlst              Do not run MLST
  --noquast             Do not run Quast
  --nocov               Do not calculate X
  --klebs               Run Kleborate, with option --all
  --argannot            Search the ARGannot databse using Abricate
  --resfinder           Search the RESfinder databse using Abricate
  --plasmidfinder       Search the PlasmidFinder databse using Abricate
  --card                Search the CARD databse using Abricate
  --ncbi                Search the NCBI databse using Abricate
  --ecoh                Search the ECOH databse using Abricate
  --abricate_all        Search all the above databases using Abricate


```

You can add more FASTQ-files to the same output-directories/summary-report, by adding the FASTQ-files to the inital input-directory, leaving the output-folder structure as it was created and re-unning the pipeline.

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
* That you do not have a high number of contigs (ideally <700)
* That you do not have low coverage (ideally >30X)
* A low number og long contigs is preferable to a high number of contigs with short contigs

## License
[GNU General Public License, v3](https://www.gnu.org/licenses/gpl-3.0.html)

## Updates
2021-04-21: Added version checks, fixed a bug with threading and most importantly: Added the flag `--no_correct` to the unicycler command to turn off spades read error correction. This is not needed as the files are QC'd with trim-galore first.

2019-07-18: Added options to find \*fastq-gz files in subdirectories from previuos runs

2019-07-18: Added options to not run parts of the pipeline, and added option to run kleborate (https://github.com/katholt/Kleborate) at the end of the pipeline

2019-10-29: Added options to run ABRICATE as part of the pipeline, and created output-folder "analyses" to put mlst, kleborate and ABRICATE-outputs in. Also added merge_runs.sh which you can use to merge two parent folders with the same structure (from this script). 

2019-10-30: Updated tool versions in README

