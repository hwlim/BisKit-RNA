#!/usr/bin/env python3

import argparse
import os
import re
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#parse command line arguments
options = argparse.ArgumentParser(description="draw methRate and coverage distribution plots using call.tsv", usage="python3 rBis.draw_dist_plots.py [options] -c call.tsv")
options.add_argument('-c', '--candidateFile', required=True,
                        help='call.tsv file from call step output')
options.add_argument('-o', '--outFile', default="out",
                        help='output file destination.')
options.add_argument('-n', '--sampName',
                        help='sample name.')
options.add_argument('-s', '--sig', default=0.05,
                        help='Significance threshold; default = 0.05')
options.add_argument('-t', '--sigType', default="pVal",
                        help='Significance type; default = pVal; use pVal or FDR')  


args = options.parse_args()
candidates = args.candidateFile
outFilePrefix = args.outFile
sigThresh = float(args.sig)
sampleName = args.sampName
if args.sigType == "pVal":
    sigType = "p-value_mState"
    sigSave = "pVal"
elif args.sigType == "FDR":
    sigType = "FDR_mState"
    sigSave = "FDR"
else:
    print("need to use pVal or FDR as -t flag.")
    exit()

df = pd.read_csv(candidates, sep='\t', header=0)

## draw distribution plots
# coverage dist
cov = np.array(df["cov"])
cov = np.log2(cov)
covPlt = plt.hist(cov, bins = 100)
covPlt = plt.xlabel("log2(Coverage)")
covPlt = plt.ylabel("Frequency")
covPlt = plt.title("log2(Coverage) Distribution\n" + sampleName)
plt.savefig(f"{outFilePrefix}/log2_cov.png", dpi=300)
plt.close()

# methRate
mr = np.array(df["methRate"])
mr = mr.astype(float)
mrPlt = plt.hist(mr, bins = 100)
mrPlt = plt.xlabel("Methylation Rate")
mrPlt = plt.ylabel("Frequency")
mrPlt = plt.title("Methylation Rate Distribution\n" + sampleName)
plt.savefig(f"{outFilePrefix}/methRate.png", dpi=300)
plt.close()

fdr = df[(df[sigType] < sigThresh )]

# coverage dist w/ FDR threshold
cov = np.array(fdr["cov"])
cov = np.log2(cov)
covPlt = plt.hist(cov, bins = 100)
covPlt = plt.xlabel("log2(Coverage)")
covPlt = plt.ylabel("Frequency")
covPlt = plt.title("log2(Coverage) Distribution after " + sigSave + " filter\n" + sampleName)
plt.savefig(f"{outFilePrefix}/log2_cov_{sigSave}.png", dpi=300)
plt.close()

# methRate w/ FDR threshold
mr = np.array(fdr["methRate"])
mr = mr.astype(float)
mrPlt = plt.hist(mr, bins = 100)
mrPlt = plt.xlabel("Methylation Rate")
mrPlt = plt.ylabel("Frequency")
mrPlt = plt.title("Methylation Rate Distribution after " + sigSave + " filter\n" + sampleName)
plt.savefig(f"{outFilePrefix}/methRate_{sigSave}.png", dpi=300)
plt.close()