#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
import numpy as np
pd.options.display.max_colwidth = 200000
pd.set_option('display.max_rows', 500)

#parse command line arguments
options = argparse.ArgumentParser(description="Add strand info to the mpileup output file.", usage="python mpileup_getStrand.py (options) -c mpileup_out_raw.tsv")
options.add_argument('-c', '--call', required=True,
                        help='mpileup_out_raw.tsv file from run_samtools_mpileup')
options.add_argument('-o', '--outFile', default="mpileup.tsv",
                        help='output file destination.')

args = options.parse_args()
inCall = args.call
outFile = args.outFile

forwardFilterList = [".","A","C","G","T","N"]
reverseFilterList = [",","a","c","g","t","n"]

## get strand, and filter out any alignments that aren't relevant for bis-seq
def getStrand(row):
    ## get reads forward aligned to G or g
    if (row['refBase'] == "G") or (row['refBase'] == "g"):
        if (any(char in row["readBase"] for char in reverseFilterList) == True):
            return "-"
        else:
            return "="
    
    ## get reads reverse aligned to C or c
    elif (row['refBase'] == "C") or (row['refBase'] == "c"):
        if (any(char in row["readBase"] for char in forwardFilterList) == True):
            return "+"
        else:
            return "="
    
    else:
        return "="

    
## get converted Base Count
def getUnconvertedBaseCount(row):
    bases = row['readBase']
    refBase = row['refBase']
    
    if (refBase == "C") or (refBase == "c"):    
        convCount = bases.count('.')
    
    elif (refBase == "G") or (refBase == "g"):    
        convCount = bases.count(',')
    
    return convCount


## get Unconverted Base Count
def getConvertedBaseCount(row):
    bases = row['readBase']
    refBase = row['refBase']
    
    if (refBase == "C") or (refBase == "c"):    
        convCount = bases.count('T')
    
    elif (refBase == "G") or (refBase == "g"):    
        convCount = bases.count('a')

    return convCount

try:
    # Code that may raise pandas.errors.EmptyDataError
    df = pd.read_csv(inCall, sep='\t', header=0)

except pd.errors.EmptyDataError:
    df = pd.DataFrame(columns=["ref", "pos", "strand", "convertedBaseCount", "unconvertedBaseCount"])
    df.to_csv(outFile, index=False, sep = "\t", header=True)
    exit()

filtered = df
filtered['strand'] = filtered.apply(getStrand, axis=1)
newFiltered = filtered.loc[(filtered['strand'] != "=")]
newFiltered['unconvertedBaseCount'] = newFiltered.apply(getUnconvertedBaseCount, axis=1)
newFiltered['convertedBaseCount'] = newFiltered.apply(getConvertedBaseCount, axis=1)

rearrange = ["ref", "pos", "strand", "convertedBaseCount", "unconvertedBaseCount"]
newFiltered=newFiltered[rearrange]

newFiltered.to_csv(outFile, index=False, sep = "\t", header=True)
