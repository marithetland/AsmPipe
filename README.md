# AsmPipe

NOTE: This is still under development, there might be bugs when using nondefault options...

Pipeline of tools used for quality assessment, assembly and gene detection from short-read Illumina sequences. Created for easy use at our microbiology research lab at Stavanger University Hopsital.

This script will quality- and adapter-trim your data, create a FastQC report, assemble the trimmed FASTQ-reads using UniCycler, create a Quast QC-report, run mlst to identify sequence types and species, and will also calculate overall coverage (sequence depth) of each sample.

The script creates a report summarising for each sample: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and sequence depth.

2019-07-18: Added options to not run parts of the pipeline, and added option to run kleborate (https://github.com/katholt/Kleborate) at the end of the pipeline


## Table of Contents

[Requirements](#Requirements)  
[Usage](#Usage)  
[ExampleCommand](#Example-command)  
[Output](#Output)  

## Requirements
These need to be installed for the entire pipeline to work. Other versions of these tools will possibly work too, but these are the I have tested with.

* Linux or MacOS
* Python 3.5.0+
* FastQC v0.11.7 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* MultiQC v1.7 (https://multiqc.info/)
* TrimGalore v0.5.0 (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
* SPAdes v3.13.0 (http://cab.spbu.ru/software/spades/)
* Pilon v1.22 (https://www.broadinstitute.org/gaag/pilon)
* Unicycler v0.4.4 (https://github.com/rrwick/Unicycler)
* Quast v4.6.3 (http://quast.sourceforge.net/quast)
* mlst v2.16.2 (https://github.com/tseemann/mlst)
* BWA v0.7.17-r1188 (http://bio-bwa.sourceforge.net/)
* SAMtools v1.7 (http://www.htslib.org/download/)
* PicardTools v2.17.8-SNAPSHOT (https://broadinstitute.github.io/picard/)
* Optional: Kleborate v0.3.0 (https://github.com/katholt/Kleborate)


## Usage

You must be in the directory containing the FASTQ-files to run this pipeline. Output-files will be stored in a specific file-structure in the input-directory.

Usage:

```
AMR-NGS

usage: asmpipe.py [-h] [-v] [-t THREADS] [--noex] [--nofqc] [--nomlst]
                  [--noquast] [--nocov] [--klebs]


Optional arguments:
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
  --klebs               Run Kleborate

```

You can add more FASTQ-files to the same output-directories/summary-report, by adding the FASTQ-files to the inital input-directory, leaving the output-folder structure as it was created and re-unning the pipeline.

Note: This was initially created for scientist with little/no coding-experience to easily perform assembly, therefore, this script currently only works when you run it from the folder the FASTQ-files are in, and output-files are stored in the same directory. In the future, I will add input and output-options.
 

## Example command

``` 
cd ~/Directory_with_fastq/ #Enter directory with FASTQ-files
asmpipe.py #Run pipeline
```

## Output

The following output-files are created when running AsmPipe:

* Fastq_raw: Any FASTQ-files that have been processed are placed here
* trimmed_reads: The trimmed reads from TrimGalore are placed here
* assembly: The Unicycler-assembled files will be placed here
* assemblies: The FASTA-file from assembly will be copied to this direcory
* QC: Contains reports from FastQC, multiQC, Quast, trimming and overall coverage calculation
* sequence_list.txt: List of all samples that have been analysed
* successful_sequences.txt: List of all samples that were successfully assembled
* failed_sequences.txt: List of any samples that failed any stage of the pipeline
* logs: Will contain run-logs from each tool for each sample
* AsmPipe_date_time.csv: Overall summary report of: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and sequence depth.
