#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import scipy.stats as stats
import statistics
import math
pd.set_option('display.max_columns', None)

options = argparse.ArgumentParser(description="Merge all pairwise comparison results.", usage="python3 rBis.merge_all_pairwise.py -p pairwise1.tsv pairwise2.tsv ...")
options.add_argument('-p','--pairwises', nargs='+',
                        help='Required; comma separated list of alignStat.txt files for all samples.', required=True)     
options.add_argument('-o', '--outFile', default='~',
                        help='Path to output compare.tsv file')
args = options.parse_args()

args = options.parse_args()
pairwises = args.pairwises
outFile = args.outFile

dfList = []
comparisons = []

## modify format of sample files to merge them easily
for i in range(len(pairwises)):
    dfTmp = pd.read_csv(pairwises[i], sep = '\t')
    comparisonName = pairwises[i].split("/")[-3]
    
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
masterDF[['#SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene']] = masterDF["SeqID"].str.split(";", expand=True)

newCols = ['#SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene', ]
masterDF = masterDF[ newCols + [col for col in masterDF.columns if col not in newCols ] ]
masterDF = masterDF[ comparisons + [col for col in masterDF.columns if col not in comparisons ] ]
masterDF = masterDF.drop(columns=['SeqID'])

masterDF[comparisons] = masterDF[comparisons].fillna("No Coverage")

## get sample names from comparisons
sampleNames = []
for comp in comparisons:
    sampleNames = sampleNames + comp.split("_vs_")


## remove duplicates from list and sort
sampleNames = [*set(sampleNames)]
sampleNames.sort()

masterDF.to_excel(outFile + "/all_sample_comparisons.xlsx", index=False, header=True)
masterDF.to_csv(outFile + "/all_sample_comparisons.tsv", index=False, sep="\t", header=True)