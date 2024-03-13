#!/usr/bin/env python3

import argparse
import os
import re
import sys
import pandas as pd
import numpy as np

#parse command line arguments
options = argparse.ArgumentParser(description="Merge m5C candidates among all sources for one sample.", usage="python mergeCandidates.py (options) -c call_rRNA.tsv call_tRNA.tsv ... ")
options.add_argument('-c', '--candidates', nargs='+',
                        help='Required; comma separated list of call.tsv files for all sources of one sample.', required=True)       
options.add_argument('-o', '--outputDest', default="mergedCandidates.tsv",
                        help='Output file destination; default = mergedCandidates.tsv')

args = options.parse_args()
candidates = args.candidates
outFile = args.outputDest

dfs = []
for cand in candidates:
    df = pd.read_csv(cand, header=0, sep="\t")
    dfs.append(df)

## add "Source", "gene_type", and "gene" columns
dfs[0]["Source"] = "rRNA"

dfs[1]["Source"] = "tRNA"

dfs[2]["Source"] = "miRNA"

dfs[3]["Source"] = "piRNA"

dfs[4]["Source"] = "Genome"

dfs[5]["Source"] = "circRNA"

rearrange = ["#SeqID", "refPos", "refStrand", "refBase", "cov", "C_count", "methRate", "uniqCov", "C_count_uniq", "methRate_uniq", "95_CI_lower", "95_CI_upper", "p-value_mState", "FDR_mState", "scores", "seqContext", "genomicCoords", "Source", "gene_type", "gene"]

dfs[5]=dfs[5][rearrange]
tmpDF = pd.concat(dfs)

tmpDF.to_csv(outFile, sep = "\t", header=True, index=False)