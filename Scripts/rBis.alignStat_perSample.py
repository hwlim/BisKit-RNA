#!/usr/bin/env python3

import argparse
import os
import re
import sys
import pandas as pd
import numpy as np
pd.options.display.max_colwidth = 200000

#parse command line arguments
options = argparse.ArgumentParser(description="Organizes the alignment statistics output by hisat-3n to machine-friendly format.", usage="python alignStat_perSample.py [options] -n sampleName -s summary.out ")
options.add_argument('-s', '--hisatSummaryOut', required=True,
                        help='summary.out file from hisat2 alignment')
options.add_argument('-o', '--outputDest', default="alignStat.txt",
                        help='Output file destination; default = alignStat.txt')
options.add_argument('-n', '--sampName', required=True,
                        help='Get sample name.')

args = options.parse_args()
summary = args.hisatSummaryOut
outFile = args.outputDest
sampName = args.sampName

df = pd.read_csv(summary, header=None)
df = df[0].str.split(' ', expand=True)

totalReads = int(df[0][0])
unaligned = int(df[4][2])
uniqAligned = int(df[4][3])
uniqAlignedPerc = float(uniqAligned / totalReads * 100)
multiAligned = int(df[4][4])
multiAlignedPerc = float(multiAligned / totalReads * 100)
overallAligned = uniqAligned + multiAligned
overallAlignedPerc = float(overallAligned / totalReads * 100)

data = {"Sample": [sampName],
        "Uniquely Aligned": [uniqAligned],
        "Multi-Aligned": [multiAligned],
        "Unaligned": [unaligned]}

table = pd.DataFrame(data)
table.to_csv(outFile, index=False, sep = "\t", header=True)