#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
pd.options.display.max_colwidth = 200000

options = argparse.ArgumentParser(description="Create merged table for mpileup results.", usage="python rBis.mpileup_mergeTables.py (options) -a mpileup_all.tsv -u mpileup_uniq.tsv ")
options.add_argument('-a','--inTableAll', required=True,
                        help='Required; output files from samtools mpileup.')
options.add_argument('-u','--inTableUniq', required=True,
                        help='Required; output from mpileup.')
options.add_argument('-o', '--outputFile', default='~/mpileup_out_final.tsv',
                        help='Path to output file')
args = options.parse_args()

inTableAll = args.inTableAll
inTableUniq = args.inTableUniq
outFile = args.outputFile

## read information into dataframes
allDF = pd.read_csv(inTableAll, header=0, sep = '\t')
uniqDF = pd.read_csv(inTableUniq, header=0, sep = '\t')

## need to merge two columns to use them as loci for merging
### merge cols in table and move to front
cols = ['ref', 'pos']
allDF['SeqID'] = allDF[cols].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
tmpCol = allDF.pop('SeqID')
allDF.insert(0, 'SeqID', tmpCol)
allDF = allDF.drop(columns=cols)
allDF = allDF.set_index('SeqID')

### merge cols in bed and move to front
uniqDF['SeqID'] = uniqDF[cols].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)
tmpCol = uniqDF.pop('SeqID')
uniqDF.insert(0, 'SeqID', tmpCol)
uniqDF = uniqDF.drop(columns=['ref', 'pos', 'strand'])
uniqDF = uniqDF.set_index('SeqID')

## rename cols
uniqDF = uniqDF.rename(columns={'convertedBaseCount': 'uniqConvertedBaseCount', 
                                'unconvertedBaseCount':'uniqUnconvertedBaseCount'})
masterDF = pd.concat([allDF,uniqDF], axis=1)

masterDF.reset_index(inplace=True)

print(masterDF)

newCols = ['ref', 'pos', 'strand', 'convertedBaseCount', 'unconvertedBaseCount',
            'uniqConvertedBaseCount', 'uniqUnconvertedBaseCount']

## save empty dataframe if empty
if masterDF.empty:
    masterDF = pd.DataFrame(columns=newCols)
    masterDF.to_csv(outFile, index=False, sep="\t", header=True)

else:
    ## splitting identifier to user-readable format
    masterDF[['ref', 'pos']] = masterDF[masterDF.columns[0]].str.split(";", expand=True)
    masterDF = masterDF.drop(columns=[masterDF.columns[0]])

    masterDF = masterDF[newCols + [col for col in masterDF.columns if col not in newCols ] ]
    masterDF = masterDF.dropna(subset=['convertedBaseCount'], axis=0)

    ## save
    masterDF.to_csv(outFile, index=False, sep="\t", header=True)
