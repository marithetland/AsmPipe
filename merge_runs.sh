#!/usr/bin/env bash

### 2019-10-29 ###
# If you want to merge output of ASMBL with 
# another ASMBL-run, use this script.
# Usage: cd into the folder you want to move, and run:
# bash merge_runs.sh /path/to/folder/to/merge/with

destination_folder=$1

##ToDo: Add check for / or not on input-folder

##ToDo: Add check if files already exists

echo "Merging files with Master-directory: " ${destination_folder}
echo "Merging Fastq_raw"
mv Fastq_raw/* ${destination_folder}/Fastq_raw/
echo "Merging Assembly"
mv assembly/* ${destination_folder}/assembly/
echo "Merging Assemblies"
mv assemblies/* ${destination_folder}/assemblies/
echo "Merging Logs"
mv logs/* ${destination_folder}/logs/
echo "Merging FastQC"
mv QC/fastQC/* ${destination_folder}/QC/fastQC/
echo "Merging Coverage"
mv QC/Coverage/* ${destination_folder}/QC/Coverage/
echo "Merging Success"
mv success/* ${destination_folder}/success/
echo "Merging Trimmed_reads"
mv trimmed_reads/* ${destination_folder}/trimmed_reads/
echo "Mergining analyses"
mv analyses/* ${destination_folder}/analyses/
echo "Removing all empty subdirectories"
find . -type d -empty -delete  #Remove any empty folders
