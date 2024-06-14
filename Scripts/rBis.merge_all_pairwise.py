#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import scipy.stats as stats
# import statistics
# import math
import ast
pd.set_option('display.max_columns', None)

options = argparse.ArgumentParser(description="Merge all pairwise comparison results.", usage="python3 rBis.merge_all_pairwise.py -p pairwise1.tsv pairwise2.tsv ...")
options.add_argument('-p','--pairwises', nargs='+',
                        help='Required; comma separated list of alignStat.txt files for all samples.', required=True)     
options.add_argument('-o', '--outFile', default='~',
                        help='Path to output compare.tsv file')
options.add_argument('-d', '--diffPairNameDict',
                        help='dictionary for diffPairNames')
args = options.parse_args()

args = options.parse_args()
pairwises = args.pairwises
outFile = args.outFile
diffPairNameDict = ast.literal_eval(args.diffPairNameDict)

print(diffPairNameDict)

dfList = []
comparisons = []

## modify format of sample files to merge them easily
for i in range(len(pairwises)):
    dfTmp = pd.read_csv(pairwises[i], sep = '\t')
    pairName = pairwises[i].split("/")[1]
    tmpName = diffPairNameDict[pairName]
    comparisonName="_vs_".join(tmpName)
    
    cols = ['#SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene']
    dfTmp['SeqID'] = dfTmp[cols].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)

    tmpCol = dfTmp.pop('SeqID')
    dfTmp.insert(0, 'SeqID', tmpCol)  
    dfTmp = dfTmp.drop(columns=cols)
    
    ## move identifier col to leftmost position to use as index
    dfTmp = dfTmp.set_index('SeqID')
    dfList.append(dfTmp)
    
    comparisons.append(comparisonName)

## concatenate DFs
masterDF = pd.concat(dfList, axis=1)

## remove duplicates columns like "cov_RG1"
masterDF = masterDF.loc[:,~masterDF.columns.duplicated()].copy()

## formatting index
masterDF.reset_index(inplace=True)
masterDF = masterDF.rename(columns = {'index':'SeqID'})

## splitting identifier to user-readable format
masterDF[['SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene']] = masterDF["SeqID"].str.split(";", expand=True)

## move avg_methRate and methRate to the left
avg_methRate_cols = [col for col in masterDF.columns if col.startswith('avg_methRate')]
methRate_cols = [col for col in masterDF.columns if col.startswith('methRate') and col not in avg_methRate_cols]
other_cols = [col for col in masterDF.columns if col not in avg_methRate_cols and col not in methRate_cols]

# Create the new order of columns
new_order = avg_methRate_cols + methRate_cols + other_cols

# Reorder the DataFrame
masterDF = masterDF[new_order]

## move comparisons to the left
masterDF = masterDF[ comparisons + [col for col in masterDF.columns if col not in comparisons ] ]

## reorder the rest of the columns
newCols = ['Source', 'SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'gene_type', 'gene',]
masterDF = masterDF[ newCols + [col for col in masterDF.columns if col not in newCols ] ]

masterDF[comparisons] = masterDF[comparisons].fillna("No Coverage")

## get sample names from comparisons
sampleNames = []
for comp in comparisons:
    sampleNames = sampleNames + comp.split("_vs_")


## remove duplicates from list and sort
sampleNames = [*set(sampleNames)]
sampleNames.sort()

## rename cols
for comp in comparisons:
    groups = comp.split("_vs_")
    key_found = next((key for key, value in diffPairNameDict.items() if value == groups), None)

    masterDF.rename(columns={comp : key_found}, inplace=True)


## save
masterDF.to_excel(outFile + "/raw_all_sample_comparisons.xlsx", index=False, header=True)
masterDF.to_csv(outFile + "/raw_all_sample_comparisons.tsv", index=False, sep="\t", header=True)

## make simplified version
# List of substrings to search for
substrings = ['cov', 'count', 'uniq', 'CI', 'mState', 'scores', 'delta', 'log2', 'abs', 'numReps', ]

# Create a regular expression pattern from the substrings
pattern = '|'.join(substrings)

filtered_df = masterDF.loc[:, ~masterDF.columns.str.contains(pattern, case=False, regex=True)]
filtered_df.to_excel(outFile + "/all_sample_comparisons.xlsx", index=False, header=True)
filtered_df.to_csv(outFile + "/all_sample_comparisons.tsv", index=False, sep="\t", header=True)