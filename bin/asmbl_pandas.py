#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('--mlst', type=str, required=True, help='specify mlst.tsv')
parser.add_argument('--fast_count', type=str, required=True, help='specify fast_count.tsv')
parser.add_argument('--quast', type=str, required=True, help='specify quast_transposed_report.tsv')
parser.add_argument('--multiqc', type=str, required=True, help='specify multiqc_fastqc.txt')
args = parser.parse_args()

#MLST
mlst_file = pd.read_csv(args.mlst, sep='\t', header=None, names=list(['Assembly','species', 'ST', 'al1','al2','al3','al4','al5','al6','al7']))
mlst_df = mlst_file.replace("_assembly.fasta","", regex=True)
mlst_df_sub = mlst_df[['Assembly','species','ST']]

#QUAST
quast_file = pd.read_csv(args.quast, sep='\t')
quast_df = quast_file.replace("_assembly","", regex=True)
quast_df.rename(columns={'# contigs (>= 0 bp)':'#contigs'}, inplace=True)
quast_df.rename(columns={'Total length (>= 0 bp)':'Total_length'}, inplace=True)
quast_df.rename(columns={'Largest contig':'Largest_contig'}, inplace=True)
quast_df.rename(columns={'GC (%)':'GC(%)'}, inplace=True)
quast_df_sub = quast_df[['Assembly', '#contigs','GC(%)','N50', 'L50', 'Total_length', 'Largest_contig']]

#MERGE MLST WITH QUAST
mlst_quast = pd.merge(mlst_df_sub, quast_df_sub, on='Assembly', how='outer')

#FASTCOUNT
fast_count_file = pd.read_csv(args.fast_count, sep='\t')
fast_count_df = fast_count_file.replace("_L001_R[12]_001.fastq.gz","", regex=True)
fast_count_df1 = fast_count_df.replace("_[12].fastq.gz","", regex=True)
fast_count_df1.rename(columns={'filename':'Assembly'}, inplace=True)
fast_count_df1.rename(columns={'total_length':'Read_depth'}, inplace=True)
#SUM bp in read1 and read2
sum_fast_count_df = fast_count_df1.groupby('Assembly')['Read_depth'].sum()

#MERGE FASTCOUNT WITH MLST + QUAST
mlst_quast_fast_count = pd.merge(mlst_quast, sum_fast_count_df, on='Assembly', how='outer')
#CALCULATE READ_DEPTH
mlst_quast_fast_count["Read_depth"] = mlst_quast_fast_count["Read_depth"].div(mlst_quast_fast_count["Total_length"]).round(2)

#MULTIQC-file
multiqc_file = pd.read_csv(args.multiqc, sep='\t')
multiqc_df = multiqc_file.replace("_[12]_val_[12]","", regex=True)
multiqc_df1 = multiqc_df.replace("_L001_R[12]_001_val_[12]","", regex=True)
multiqc_df1.rename(columns={'Sample':'Assembly'}, inplace=True)
multiqc_df1.rename(columns={'Total Sequences':'#Reads'}, inplace=True)
multiqc_df_sub = multiqc_df1[['Assembly', '#Reads']]
multiqc_df_sub = multiqc_df_sub.drop_duplicates() #All pairs should have same number of reads/sequences

#MERGE MULTIQC WITH MLST + QUAST + READ_DEPTH 
mlst_quast_fast_count_multiqc = pd.merge(mlst_quast_fast_count, multiqc_df_sub, on='Assembly', how='outer')

#MAKE FILE
mlst_quast_fast_count_multiqc.to_csv(path_or_buf='final_report.tsv', sep="\t", index=False)