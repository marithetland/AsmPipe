#!/usr/bin/env python3

from Bio import SeqIO
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--fastq', type=str, required=True, help='specify the fastq file')
args = parser.parse_args()


with gzip.open(args.fastq, "rt") as handle:
    is_fastq = SeqIO.parse(handle, "fastq")
    print(any(fastq))
