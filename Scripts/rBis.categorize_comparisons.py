#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import scipy.stats as stats
import statistics
import math
import matplotlib.pyplot as plt

options = argparse.ArgumentParser(description="Categorize the comparisons results by added the following labels to each candidate: UP, DOWN, UNCHANGED, NA, UNIQ1 or UNIQ2.", usage="python3 rBis.categorize_comparisons [options] -p comparison.tsv")
options.add_argument('-p','--pairwise',
                        help='Required; comma separated list of alignStat.txt files for all samples.', required=True)     
options.add_argument('-o', '--outFile', default='~',
                        help='Path to output compare.tsv file')
options.add_argument('-s', '--sig', default=0.05,
                        help='Significance threshold; default = 0.05')
options.add_argument('-t', '--sigType', default="pVal",
                        help='Significance type; default = pVal; use pVal or FDR')
options.add_argument('-d', '--diff', default=0.05,
                        help='diff threshold, default=0.05')
options.add_argument('-m', '--minReps', default=2,
                        help='minimum number of replicates a C position should be seen to be considered a candidate, default=2')
options.add_argument('-i', '--sampleFile',
                        help='sample.tsv')
options.add_argument('-n', '--comparisonPair',
                        help='comparison pair')
options.add_argument('-g', '--groupList',
                        help='group list')  

args = options.parse_args()
pairwise = args.pairwise
outFile = args.outFile
sigThresh = float(args.sig)
diffThresh = float(args.diff)
minReps = int(args.minReps)
comparison = args.comparisonPair
groupList = args.groupList.split(",")

## import group as dictionary from command line
sampleFile = args.sampleFile
sampleFile = pd.read_csv(sampleFile, sep="\t", comment="#", na_filter=False)

groups = {}
for group in groupList:
    samplesList = sampleFile.Name[sampleFile.Group == group].tolist()
    groups[group]=samplesList

if args.sigType == "pVal":
    sigType = "p-value_mState"
    sigSave = "pVal"
elif args.sigType == "FDR":
    sigType = "FDR_mState"
    sigSave = "FDR"
else:
    print("need to use pVal or FDR as -t flag.")
    exit()

## are there replicates
n = 1
for key in groups:
    if len(groups[key]) > n:
        n = len(groups[key])

if n > 1:
    reps = True
else:
    reps = False

## get up, down, unchanged, or unique
def differential(row):
    sampleDiff = row[(f'delta_MethRate_{compareSamp}')]
    sampleAbsDiff = row[(f'abs_delta_methRate_{compareSamp}')]
    sampleFDR = row[(f'{sigType}_{compareSamp}')]
    result = ""

    ## if there are replicates for this experiment
    if reps:
        
        numReps = []
        for grou in groups:
            numReps.append(row[(f'numReps_{grou}')])
        
        ## if first sample does not have enough replicates and second sample has enough replicates, uniq to sample 2
        if (numReps[0] < minReps) and (numReps[1] >= minReps):
            result = "uniq2"
        
        ## if second sample does not have enough replicates and first sample has enough replicates, uniq to sample 1
        elif (numReps[1] < minReps) and (numReps[0] >= minReps):
            result = "uniq1"
        
        ## if both samples don't have enough replicates, NA
        elif (numReps[0] < minReps) and (numReps[1] < minReps):
            result = "NA"


    ## if no replicates for this experiment
    elif reps == False:
        covSample1 = row[(f'cov_{sampleNames[0]}')]
        covSample2 = row[(f'cov_{sampleNames[1]}')]
    
        if math.isnan(sampleFDR):
            if math.isnan(float(covSample1)) and math.isnan(float(covSample2)):
                result = "NA"
            elif math.isnan(float(covSample1)) and float(covSample2) > 0:
                result = "uniq2"
            elif math.isnan(float(covSample2)) and float(covSample1) > 0:
                result = "uniq1"
            else:
                result = "UNCHANGED"

    ## if candidate is not unique to a group or is not NA:
    if result == "":
        if (sampleFDR > sigThresh) or (sampleAbsDiff < diffThresh):
            result = "UNCHANGED"
        elif (sampleFDR < sigThresh) and (sampleAbsDiff >= diffThresh) and (sampleDiff > 0):
            result = "UP"
        elif (sampleFDR < sigThresh) and (sampleAbsDiff >= diffThresh) and (sampleDiff < 0):
            result = "DOWN"
        else:
            result = "UNCHANGED"

    ## Error handling
    try:
        return result
    
    except UnboundLocalError:
        print("ERROR")
        print(row)
        print(sampleDiff)
        print(sampleAbsDiff)
        print(sampleFDR)
        print(covSample1)
        print(covSample2)
        print()
        exit()

## Function to move a column to a specific position
def move_column_inplace(df, col, pos):
    if col in df.columns:
        ## Extract the column and delete it from original position
        col_data = df.pop(col)
        ## Reinsert it at the desired position
        df.insert(pos, col, col_data)
    else:
        print(f"Column {col} not found in DataFrame.")

## get sampleNames
sampleNames = groupList

## get comparison
compareSamp = sampleNames[0] + "_vs_" + sampleNames[1]

## make DF
masterDF = pd.read_csv(pairwise, header=0, sep = '\t')

## check if the control only has one replicate
ctrlRep1 = False
if reps:
    ctrlSampleNameCol = masterDF.filter(like='numReps').columns[0]
    numRepsCtrlUnique, numRepsCtrlCounts = np.unique(masterDF[ctrlSampleNameCol], return_counts=True)

    if numRepsCtrlUnique[-1] == 1:
        print("no Reps for control, pretending there are two")
        masterDF[ctrlSampleNameCol] = masterDF[ctrlSampleNameCol].replace(1, 2)
        ctrlRep1 = True

## check if the sample only has one replicate
sampRep1 = False
if reps:
    sampleNameCol = masterDF.filter(like='numReps').columns[-1]
    numRepsSampUnique, numRepsSampCounts = np.unique(masterDF[sampleNameCol], return_counts=True)

    if numRepsSampUnique[-1] == 1:
        print("no Reps for the sample, pretending there are two")
        masterDF[sampleNameCol] = masterDF[sampleNameCol].replace(1, 2)
        sampRep1 = True

## get diff results
masterDF[compareSamp] = masterDF.apply(differential, axis=1)

## organize columns
newCols = [compareSamp, '#SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene' ]
masterDF = masterDF[newCols + [col for col in masterDF.columns if col not in newCols ] ]

if ctrlRep1 == True:
    masterDF[ctrlSampleNameCol] = masterDF[ctrlSampleNameCol].replace(2, 1)

if sampRep1 == True:
    masterDF[sampleNameCol] = masterDF[sampleNameCol].replace(2, 1)


## save categorization
masterDF.to_csv(outFile + "/categorized_pairwise_comparison.tsv", index=False, sep="\t", header=True)

## get stats
unique, counts = np.unique(masterDF[compareSamp], return_counts=True)
unique = list(unique)
counts = list(counts)

## add comparison identifier
unique = ['Comparison'] + unique
counts = [comparison] + counts

## remove NA if there are replicates; they won't be needed for plotting
if reps:
    # Find the index of 'NA' in list1
    if 'NA' in unique:
        index = unique.index('NA')
    
        # Remove the item 'NA' from list1 and the corresponding item from list2
        unique.pop(index)
        counts.pop(index)

## make dataframe and name columns by unique list
df2 = pd.DataFrame(counts).T
df2.columns = unique
# df2.rename(columns={f'Unique to {sampleNames[0]}': 'uniq1'}, inplace=True)
# df2.rename(columns={f'Unique to {sampleNames[1]}': 'uniq2'}, inplace=True)

## organize columns
df2.columns = ['Comparison'] + df2.columns[1:].tolist()
move_column_inplace(df2, 'UP', 1)
move_column_inplace(df2, 'DOWN', 2)
move_column_inplace(df2, 'UNCHANGED', 3)
move_column_inplace(df2, 'uniq1', 4)
move_column_inplace(df2, 'uniq2', 5)

## save stats
df2.to_csv(outFile + "/stats.tsv", index=False, sep = "\t", header=True)

## filter DF and draw density
fdr = masterDF[(masterDF[f'{sigType}_{compareSamp}'] < sigThresh )]

# methRate
mr = np.array(masterDF["log2_methRateFC_" + compareSamp])
mr = mr.astype(float)
mrPlt = plt.hist(mr, bins = 100)
mrPlt = plt.xlabel(f"Methylation Rate Log2({sampleNames[1]}/{sampleNames[0]})")
mrPlt = plt.ylabel("Frequency")
mrPlt = plt.title("Methylation Rate Log2 Fold Change\n" + sampleNames[0] + " vs " + sampleNames[1])
plt.savefig(outFile + '/Plots/methRate_log2foldChange.png', dpi=300)
plt.close()

# fdr methRate
mr = np.array(fdr["log2_methRateFC_" + compareSamp])
mr = mr.astype(float)
mrPlt = plt.hist(mr, bins = 100)
mrPlt = plt.xlabel(f"Methylation Rate Log2({sampleNames[1]}/{sampleNames[0]})")
mrPlt = plt.ylabel("Frequency")
mrPlt = plt.title("Methylation Rate Log2 Fold Change After " + sigSave + " filter\n" + sampleNames[0] + " vs " + sampleNames[1])
plt.savefig(outFile + '/Plots/methRate_log2foldChange_' + sigSave + '.png', dpi=300)
plt.close()

# delta methRate
deltamr = np.array(masterDF["delta_MethRate_" + compareSamp])
deltamr = deltamr.astype(float)
deltamrPlt = plt.hist(deltamr, bins = 100)
deltamrPlt = plt.xlabel(f"Delta Methylation Rate ({sampleNames[1]}−{sampleNames[0]})")
deltamrPlt = plt.ylabel("Frequency")
deltamrPlt = plt.title("Delta Methylation Rate\n" + sampleNames[0] + " vs " + sampleNames[1])
plt.savefig(outFile + '/Plots/delta_methRate.png', dpi=300)
plt.close()

# delta methRate fdr
deltamr = np.array(fdr["delta_MethRate_" + compareSamp])
deltamr = deltamr.astype(float)
deltamrPlt = plt.hist(deltamr, bins = 100)
deltamrPlt = plt.xlabel(f"Delta Methylation Rate ({sampleNames[1]}−{sampleNames[0]})")
deltamrPlt = plt.ylabel("Frequency")
deltamrPlt = plt.title("Delta Methylation Rate After " + sigSave + " filter\n" + sampleNames[0] + " vs " + sampleNames[1])
plt.savefig(outFile + '/Plots/delta_methRate_' + sigSave + '.png', dpi=300)
plt.close()