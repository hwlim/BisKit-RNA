#!/usr/bin/env python3

import argparse
import pandas as pd

#parse command line arguments
options = argparse.ArgumentParser(description="merge call statistics", usage="python merge_all_callStats.py [options] -s callStats1.tsv callStats2.tsv callStats3.tsv... -a callStatsAll1.tsv callStatsAll2.tsv ...")
options.add_argument('-s','--sig', nargs='+',
                        help='Required; simplified call stats.tsv file.', required=True)
options.add_argument('-a','--all', nargs='+',
                        help='Required; all call stats.tsv file.', required=True)
options.add_argument('-o', '--outPrefix', default="~",
                        help='Output file destination; default = home directory')

args = options.parse_args()
sigAndNonsig = args.sig
allCandidates = args.all
outPrefix = args.outPrefix


sigDFs = []
for stat in sigAndNonsig:
    df = pd.read_csv(stat, header=0, sep="\t")
    sigDFs.append(df)
    
totalDFs = []
for stat in allCandidates:
    df = pd.read_csv(stat, header=0, sep="\t")
    totalDFs.append(df)


result = sigDFs[0].copy()
for df in sigDFs[1:]:
    result.iloc[:, 1:] += df.iloc[:, 1:]

result.to_csv(f'{outPrefix}/callStats_sigNonsig.tsv', sep = "\t", header=True, index=False)


result = totalDFs[0].copy()
for df in totalDFs[1:]:
    result.iloc[:, 1:] += df.iloc[:, 1:]

result.to_csv(f'{outPrefix}/callStats_all.tsv', sep = "\t", header=True, index=False)
