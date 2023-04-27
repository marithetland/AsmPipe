#!/usr/bin/env python3

import pandas as pd
from argparse import ArgumentParser

#input from nextflow will be:
#fast_count.csv (add together bp from read1 and read2)
#quast_transposed_report.tsv (#contigs)


#Final report file:
#assembly(name) species ST  #contigs    GC (%)  N50 L50 Total_length    Largest contig  Avg_readDepth   StDev   #Reads
#whatever   mlst    mlst    quast   quast   quast   quast   quast   quast   (fast_count1+fastcount2)/total_length(quast)    ?   multiqc_fastqc.txt


#species and ST : mlst.tsv

# #contigs, GC (%), N50, L50, Total_length, Largest contig (transposed_report.tsv)

# Avg_readDepth (fast_count1+fastcount2)/total_length(quast)

#read_depth_file = pd.read_csv('fast_count.tsv', sep='\t', header=None, names=list(['Assembly','species', 'ST', 'al1','al2','al3','al4','al5','al6','al7']))


#in nextflow
#process    pandas {
#     input:
#     path(mlst), path(fast_count), path(quast), path(multiqc)

#     script:
#     """
#     python asmbl_pandas.py --mlst $mlst --fast_count $fast_count --quast $quast --multiqc $multiqc
#     """
# }

parser = ArgumentParser()
parser.add_argument('--mlst', type=str, required=True, help='specify mlst.tsv')
parser.add_argument('--fast_count', type=str, required=True, help='specify fast_count.tsv')
parser.add_argument('--quast', type=str, required=True, help='specify quast_transposed_report.tsv')
parser.add_argument('--multiqc', type=str, required=True, help='specify multiqc_fastqc.txt')
args = parser.parse_args()


#MLST-file, reduses the df to only 3 columns: Assembly, species, ST
mlst_file = pd.read_csv(args.mlst, sep='\t', header=None, names=list(['Assembly','species', 'ST', 'al1','al2','al3','al4','al5','al6','al7']))
mlst_df = mlst_file.replace("_assembly.fasta","", regex=True)
mlst_df_sub = mlst_df[['Assembly','species','ST']]

#Make a seq_df to merge it with mlst_file? Is it necessary?
# seq_df_mlst = pd.merge(seq_df, mlst_df_sub, on='Assembly', how='outer')


#QUAST-file
quast_file = pd.read_csv(args.quast, sep='\t')
quast_df = quast_file.replace("_assembly","", regex=True)
#TODO: Edit _ to - in quast
quast_df.rename(columns={'# contigs (>= 0 bp)':'#contigs'}, inplace=True)
quast_df.rename(columns={'Total length (>= 0 bp)':'Total_length'}, inplace=True)
quast_df_sub = quast_df[['Assembly', '#contigs','GC (%)','N50', 'L50', 'Total_length', 'Largest contig']]
#Same question, can I skip the step by making seq_df?
#seq_df_mlst_quast = pd.merge(seq_df_mlst, quast_df_sub, on='Assembly', how='outer')

mlst_quast = pd.merge(mlst_df_sub, quast_df_sub, on='Assembly', how='outer')

#test-file
#mlst_quast.to_csv(path_or_buf='mlst_quast.tsv', sep='\t')

#MULTIQC-file
multiqc_file = pd.read_csv(args.multiqc, sep='\t')
#test-file
#multiqc_file.to_csv(path_or_buf='multiqc_file.tsv', sep='\t')

#Can I write these as "_[12]_val_[12]" or something? Making it into 2 lines instead of 4. Same for fast_count...
#Maybe I can skip the multi_df1= etc. and only write .replace like .rename
multiqc_df1 = multiqc_file.replace("_1_val_1","", regex=True)
multiqc_df2 = multiqc_df1.replace("_2_val_2","", regex=True)
multiqc_df3 = multiqc_df2.replace("_L001_R1_001_val_1","", regex=True)
multiqc_df4 = multiqc_df3.replace("_L001_R2_001_val_2","", regex=True)
#test-file
#multiqc_df4.to_csv(path_or_buf='multiqc_df.tsv', sep='\t')


multiqc_df4.rename(columns={'Sample':'Assembly'}, inplace=True)
multiqc_df4.rename(columns={'Total Sequences':'#Reads'}, inplace=True)
multiqc_df_sub = multiqc_df4[['Assembly', '#Reads']]
multiqc_df_sub = multiqc_df_sub.drop_duplicates() #All pairs should have same number of reads/sequences

#test-file
#multiqc_df_sub.to_csv(path_or_buf='multiqc_sub.tsv', sep='\t')

    
mlst_quast_multiqc = pd.merge(mlst_quast, multiqc_df_sub, on='Assembly', how='outer')
#test-file
mlst_quast_multiqc.to_csv(path_or_buf='mlst_quast_multiqc.tsv', sep='\t')

fast_count_file = pd.read_csv(args.fast_count, sep='\t')

#Can I write these as "_[12]_val_[12]" or something? Making it into 2 lines instead of 4. Same for fast_count...
#Maybe I can skip the multi_df1= etc. and only write .replace like .rename
fast_count_df1 = fast_count_file.replace("_L001_R1_001.fastq.gz","", regex=True)
fast_count_df2 = fast_count_df1.replace("_L001_R2_001.fastq.gz","", regex=True)
fast_count_df3 = fast_count_df2.replace("_1.fastq.gz","", regex=True)
fast_count_df4 = fast_count_df3.replace("_2.fastq.gz","", regex=True)

fast_count_df4.rename(columns={'filename':'Assembly'}, inplace=True)
fast_count_df4.rename(columns={'total_length':'#bp'}, inplace=True)

sum_fast_count_df = fast_count_df4.groupby('Assembly')['#bp'].sum()
#test-file
sum_fast_count_df.to_csv(path_or_buf='sum_fast_count.tsv', sep='\t')

mlst_quast_multiqc_fast_count = pd.merge(mlst_quast_multiqc, sum_fast_count_df, on='Assembly', how='outer')

#test-file
mlst_quast_multiqc_fast_count.to_csv(path_or_buf='final_report.tsv', sep="\t")


#seq_df_mlst_quast_cov_fastqc.to_csv(path_or_buf='Asmbl_'+todays_date+'.csv', sep="\t")
