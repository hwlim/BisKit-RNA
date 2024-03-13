#!/usr/bin/env python3

import argparse
import os
import re
import sys
import pandas as pd
import numpy as np

#parse command line arguments
options = argparse.ArgumentParser(description="merge call statistics from all sources for one sample.", usage="python mergeCallStats.py (options) -a call_stats_rRNA.txt call_stats_tRNA.txt... -s sample1")
options.add_argument('-a','--callStats', nargs='+',
                        help='Required; comma separated list of call_stats.txt files for all samples.', required=True)      
options.add_argument('-o', '--outputDest', default="mergedCallStat.tsv",
                        help='Output file destination; default = mergedCallStats.tsv')
options.add_argument('-s', '--sampleName', required=True,
                        help='Sample name.')

args = options.parse_args()
stats = args.callStats
outFile = args.outputDest
sampName = args.sampleName

dfs = []
for stat in stats:
    df = pd.read_csv(stat, header=0, sep="\t")
    dfs.append(df)

tmpDF = pd.concat(dfs)

indexed = tmpDF.reset_index()
final = indexed.rename(index={0: "rRNA", 1: "tRNA", 2: "miRNA", 3: "piRNA", 4: "genome", 5: "circRNA"})
fin = final.drop(['index'], axis=1)
fin.to_csv(f'{outFile}.tsv', sep = "\t", header=True)

indexed = tmpDF.reset_index()
final = indexed.rename(index={0: f"{sampName} rRNA", 1: f"{sampName} tRNA", 2: f"{sampName} miRNA", 3: f"{sampName} piRNA", 4: f"{sampName} genome", 5: f"{sampName} circRNA"})
fin = final.drop(['index'], axis=1)
fin.to_csv(f'{outFile}NoHeader.tsv', sep = "\t", header=False)

