#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

options = argparse.ArgumentParser(description="Merge categorizations by source for all pairs.", usage="python rBis.mergeCatStatsBySource (options) -p Categorization_By_Source1.tsv, Categorization_By_Source2.tsv ...")
options.add_argument('-p','--pairwises', nargs='+',
                        help='Required; comma separated list of alignStat.txt files for all samples.', required=True)     
options.add_argument('-o', '--outPrefix', default='~',
                        help='Path to output files')
args = options.parse_args()

args = options.parse_args()
pairwises = args.pairwises
outPrefix = args.outPrefix

dfList = []
comparisons = []

for i in range(len(pairwises)):
    dfTmp = pd.read_csv(pairwises[i], sep = '\t', header=None)
    comparisonName = pairwises[i].split("/")[-3]

    dfList.append(dfTmp)
    
    comparisons.append(comparisonName)

masterDF = pd.concat(dfList, axis=0)

## split first column by space
masterDF[['Comparison', 'Source']] = masterDF[0].str.split(' ', expand=True)
masterDF = masterDF.drop(0, axis=1)
masterDF = masterDF.reindex(columns=['Comparison', 'Source'] + list(masterDF.columns[:-2]))

## rename columns
masterDF = masterDF.rename(columns={
    1: 'UP',
    2: 'DOWN',
    3: 'UNCHANGED',
    4: 'UNIQ1',
    5: 'UNIQ2'
})

## split dataframe by Source
for Source, group in masterDF.groupby('Source'):
    group = group.drop("Source", axis=1)
    group.to_csv(f'{outPrefix}/{Source}_catStats.tsv', sep = "\t", index=False)

    group = group.drop("UNCHANGED", axis=1)
    group = group.drop("UNIQ1", axis=1)
    group = group.drop("UNIQ2", axis=1)
    group.to_csv(f'{outPrefix}/{Source}_catStats_updown.tsv', sep = "\t", index=False)
