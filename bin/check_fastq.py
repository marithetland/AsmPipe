#!/usr/bin/env python3

import gzip
from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--fastq', type=str, required=True, help='specify the fastq file')
args = parser.parse_args()

#check if empty
with open(args.fastq, 'rb') as f:
    file_content = f.read(1)
    print(len(file_content) > 0)

#Check if fastq format
with gzip.open(args.fastq, "rt") as handle:
    is_fastq = SeqIO.parse(handle, "fastq")
    print(all(is_fastq))
