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
args = parser.parse_args()

mlst_file = pd.read_csv(args.mlst, sep='\t', header=None, names=list(['Assembly','species', 'ST', 'al1','al2','al3','al4','al5','al6','al7']))
#mlst_file.to_csv(path_or_buf='mlst_with_header.tsv', sep='\t')
mlst_df = mlst_file.replace("_assembly.fasta","", regex=True)
mlst_df_sub = mlst_df[['Assembly','species','ST']]
mlst_df_sub.to_csv(path_or_buf='mlst_sub.tsv', sep='\t')
# seq_df_mlst = pd.merge(seq_df, mlst_df_sub, on='Assembly', how='outer')