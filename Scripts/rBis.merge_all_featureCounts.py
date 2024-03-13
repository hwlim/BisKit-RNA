#!/usr/bin/env python3

import argparse
import pandas as pd

options = argparse.ArgumentParser(description="Create merged table of featureCounts results for samples.", usage="python3 rBis.merge_all_featureCounts.py [options] -f sample1.tsv sample2.tsv sample4.tsv ... ")
options.add_argument('-f','--featureCounts', nargs='+',
                        help='Required; list of featureCounts output for both samples.', required=True)
options.add_argument('-o', '--outFile', default='out.tsv',
                        help='Path to output file')
args = options.parse_args()

samples_list = args.featureCounts
outFile = args.outFile

dfs = []

for sample in samples_list:
    tmpDF = pd.read_csv(sample, comment='#', sep="\t", header=0)
    dfs.append(tmpDF)

merged_df = dfs[0]

# Merge the remaining DataFrames using a loop
for df in dfs[1:]:
    merged_df = pd.merge(merged_df, df, on='source:gene', how='inner')

merged_df.to_csv(outFile, sep='\t', index=False)