# Asmbl-nf
 Short read assembly pipeline (nextflow)

# Asmbl
**Assembly and quality assessment of short-read Illumina data**

Pipeline of tools used for quality assessment, assembly and gene detection from short-read Illumina sequences. Created for easy use at our microbiology research lab at Stavanger University Hopsital.

This script will quality- and adapter-trim your data with Trim Galore, create a FastQC report, assemble the trimmed FASTQ reads using SPAdes (or Unicycler), create a Quast QC-report, run mlst to identify sequence types and species, and will calculate average read depth (or coverage) of each sample.

The script creates a report summarising for each sample: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and read depth.

## Table of Contents

[Installation](#Installation)   
[Usage](#Usage)  
[Output](#Output)  
[Detailed Explanation](#Detailed-explanation)  
[Updates](#Updates)  


## Installation

Install conda environment and install dependencies:

```
mamba create -n asmbl_env -c bioconda -c conda-forge pandas blast fastqc multiqc trim-galore unicycler quast openjdk==17.0.3 biopython perl-moo mlst 
```

Clone the github repo and set up:
```
git clone https://github.com/marithetland/Asmbl.git
```

Optional:

Kleborate: 
```
git clone --recursive https://github.com/katholt/Kleborate.git
cd Kleborate/kaptive
git pull https://github.com/katholt/Kaptive master
```

kmerfinder:
```
mamba install -c bioconda kmerfinder 
```

KmerFinder db:
```
wget https://cge.food.dtu.dk/services/KmerFinder/etc/kmerfinder_db.tar.gz
tar -xvf kmerfinder_db.tar.gz README.md VALIDATE.py bacteria bacteria.md5 config
rm kmerfinder_db.tar.gz
```

## Usage

You must be in the directory containing the FASTQ-files to run this pipeline. Output-files will be stored in a specific file-structure in the input-directory. In addition to the default pipeline, you can also run Kleborate, KmerFinder and rMLST.

``` 
cd ~/Directory_with_fastq/  #Enter directory with FASTQ-files
conda activate asmbl_env    #Activate conda environment
nextflow run Asmbl.nf       #Run the pipeline
```


There are several options that can be set by modifying the nextflow.config file:

| Option                            | Description                                                       | Default                    |
| ----                              | ----                                                              | ----                       |
| `unicycler050`                    | Run Unicycler v0.5.0 for assembly (true/false)                    | false                      |
| `unicycler048`                    | Run Unicycler v0.4.8 for assembly (true/false)                    | false                      |
| `unicycler048_path`               | Specify the path to Unicycler v0.4.8                              |                            |
| `pilon_uni048_path`               | Specify the path to pilon.jar                                     |                            |
| `pilon_version_path`              | Specify the path to pilon                                         |                            |
| `depth_filter`                    | Depth_filter for Unicycler                                        | 0.25                       |
| `trim`                            | Run Trim_galore (true/false)                                      | true                       |
| `kleborate`                       | Run Kleborate (true/false)                                        | false                      |
| `reads_type1`                     | Specify input reads type 1                                        | ./*L001_R{1,2}_001.fastq.gz|
| `reads_type2`                     | Specify input reads type 2                                        | ./*_{1,2}.fastq.gz         |
| `failure_action`                  | Specify the nextflow error strategy (terminate/ignore/finish)     | ignore                     |
| `rmlst`                           | Run rMLST (true/false)                                            | false                      |
| `kmerfinder`                      | Run KmerFinder (true/false)                                       | false                      |
| `kmerfinder_db`                   | Specify the path to kmerfinder_db/bacteria/                       |                            |
| `fast`                            | Run with 72 threads rather than 36 (true/false)                   | false                      |



## Output

The following output-files are created when running AsmPipe:

* fastq: Any FASTQ-files that have been processed are placed here
* trimmed_reads: The trimmed reads from Trim Galore are placed here
* fasta: The FASTA-file from assembly will be placed here
* gfa: The GFA-file from assembly will be placed here
* reports: Contains reports from FastQC, MultiQC, Quast, mlst, fast_count, Kleborate, rMLST and KmerFinder. Also contains an overall summary
  report of: Species, ST, no. reads, GC%, no. contigs, largest contig, total sequence length, N50, L50 and sequence depth.
* logs: Contains run-logs from Trim Galore and SPAdes (or Unicycler).
 

## Detailed explanation

* FastQC performs quality assessment of raw reads, indicating number of reads, GC%, adapter content, sequence length distribution, and more
* Trim Galore - trims raw reads based on adapter sequences and Phred quality: trims 1 bp off 3' end of every read, removes low-quality (<Phred 20) 3' ends, removes adapter sequences and removes read-pairs if either of the reads' length is <20 bp
* SPAdes is a short-read genome assembler. Unicycler functions as a SPAdes optimiser with short-reads only. Unicycler v0.4.8 uses pilon polishing to attempt to make improvements on the genome
* Quast quality assessment on assembly outputs the total length, GC%, number of contigs, N50, L50 and more. 
* mlst attempts to identify species and MLST based on the PubMLST schemes. Other tools may be needed for specification, e.g. Kleborate identifies locus variants for Klebsiella samples and separates klebsiella pneumoniae sensu lato into subspecies
* The overall average depth of the genome is calculated based on total length in bp from fast_count and total_length from quast.


Things to check QC-wise
* That GC% matches the sample species
* That the total length matches the sample species
* That you do not have a high number of contigs (ideally <500)
* That you do not have low average read depth (ideally >40X)
* A low number of long contigs is preferable to a high number of contigs with short contigs

## Updates
2023-10-25: Updated the pipeline to nextflow pipeline. SPAdes is now the default assembler. Added option to run KmerFinder and rMLST. Switched out bwa with fast_count for read depth calculation. 

2022-02-03: Updated pipeline to work with newest release of Unicycler (v0.5.0). Unicycler no longer uses pilon for polishing and read error correction is by default turned off, so these options have been removed from the pipeline. Also removed option for abricate as it does not work. Also updated some terms to stop confusion: The folder "assemblies" is now "fasta", "coverage" is now "read depth", and the final report from the pipeline is prefixed with "Asmbl" rather than "AsmPipe".

2021-04-21: Added version checks, fixed a bug with threading and most importantly: Added the flag `--no_correct` to the Unicycler command to turn off SPAdes read error correction. This is not needed as the files are QC'd with Trim Galore first.

2019-07-18: Added options to find \*fastq-gz files in subdirectories from previuos runs

2019-07-18: Added options to not run parts of the pipeline, and added option to run Kleborate (https://github.com/katholt/Kleborate) at the end of the pipeline

2019-10-29: Added options to run ABRICATE as part of the pipeline, and created output-folder "analyses" to put mlst, Kleborate and ABRICATE-outputs in. Also added merge_runs.sh which you can use to merge two parent folders with the same structure (from this script). 

2019-10-30: Updated tool versions in README.
