# AsmPipe

Pipeline of tools used for quality assessment, assembly and gene detection from short-read Illumina sequences. Created for easy use at our microbiology research lab at Stavanger University Hopsital.

This script will quality- and adapter-trim your data, create a FastQC report, assemble the trimmed FASTQ-reads using UniCycler, create a Quast QC-report, run mlst to identify sequence types and species, and will also calculate overall coverage (sequence depth) of each sample.

The script creates a report summarising for each sample: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and sequence depth.


## Table of Contents

[Requirements](#Requirements)  
[Usage](#Usage)  
[ExampleCommand](#Example-command)  
[Output](#Output)  

## Requirements

* Linux or MacOS
* Python 3.5.0+
* FastQC
* TrimGalore
* SPAdes
* Unicycler
* Quast
* mlst
* kleborate
* BWA
* SAMtools
* PicardTools

## Usage

You must be in the directory containing the FASTQ-files to run this pipeline. Output-files will be stored in a specific file-structure in the input-directory.

Usage:

```
AMR-NGS

Usage: asmpipe.py [-h] [-v] [--noqc]

Input options (required):
  --noqc         Do not run fastQC and multiQC

Optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit


```

You can add more FASTQ-files to the same output-directories/summary-report, by adding the FASTQ-files to the inital input-directory, leaving the output-folder structure as it was created and re-unning the pipeline.

Note: This was initially created for scientist with little/no coding-experience to easily perform assembly, therefore, this script currently only works when you run it from the folder the FASTQ-files are in, and output-files are stored in the same directory. In the future, I will add input and output-options, and will also add options for including/excluding parts of the pipeline.
 

## Example command

``` 
cd ~/Directory_with_fastq/ #Enter directory with FASTQ-files
asmpipe.py #Run pipeline
```

## Output

The following output-files are created when running AsmPipe:

* Fastq_raw: Any FASTQ-files that have been processed are placed here
* trimmed_reds: The trimmed reads from TrimGalore are placed here
* assembly: The Unicycler-assembled files will be placed here
* assemblies: The FASTA-file from assembly will be copied to this direcory
* QC: Contains reports from FastQC, multiQC, Quast, trimming and overall coverage calculation
* sequence_list.txt: List of all samples that have been analysed
* successful_sequences.txt: List of all samples that were successfully assembled
* failed_sequences.txt: List of any samples that failed any stage of the pipeline
* logs: Will contain run-logs from each tool for each sample
* AsmPipe_date_time.csv: Overall summary report of: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and sequence depth.
