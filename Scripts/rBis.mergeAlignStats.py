#!/usr/bin/env python3

import argparse
import os
import re
import sys
import pandas as pd
import numpy as np

#parse command line arguments
options = argparse.ArgumentParser(description="merge align statistics to view as a summary", usage="python3 rBis.mergeAlignStats.py [options] -s sampleName -a alignStats1.tsv alignStats2.tsv...")
options.add_argument('-a','--alignStats', nargs='+',
                        help='Required; comma separated list of alignStat.txt files for all samples.', required=True)      
options.add_argument('-o', '--outputDest', default="mergedAlignStat.tsv",
                        help='Output file destination; default = mergedAlignStat.tsv')
options.add_argument('-s', '--sampleName', required=True,
                        help='Sample name.')

args = options.parse_args()
stats = args.alignStats
outFile = args.outputDest
sampName = args.sampleName

dfs = []
for stat in stats:
    df = pd.read_csv(stat, header=0, sep="\t")
    dfs.append(df)

tmpDF = pd.concat(dfs)

indexed = tmpDF.reset_index()
final = indexed.rename(index={0: "rRNA", 1: "tRNA", 2: "miRNA", 3: "piRNA", 4: "Genome", 5: "circRNA"})
fin = final.drop(['index'], axis=1)
fin.to_csv(f'{outFile}.tsv', sep = "\t", header=True)

