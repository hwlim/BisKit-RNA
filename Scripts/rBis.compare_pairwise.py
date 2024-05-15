#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import scipy.stats as stats
import statistics
import math
import statsmodels.api as sm

options = argparse.ArgumentParser(description="Create merged table for rMATS results.", usage="python3 rBis.compare_pairwise.py [options] -c alignStat1.txt alignStat2.txt alignStat3.txt ")
options.add_argument('-c','--candidates', nargs='+',
                        help='Required; comma separated list of alignStat.txt files for all samples.', required=True)     
options.add_argument('-o', '--outFile', default='~',
                        help='Path to output compare.tsv file')
options.add_argument('-s', '--sampleFile',
                        help='sample.tsv')
options.add_argument('-p', '--comparisonPair',
                        help='comparison pair')
options.add_argument('-g', '--groupList',
                        help='group list')
options.add_argument('-m', '--minReps', default = 2,
                        help='minimum number of replicates a C location has to appear in to be considered a candidate; default = 1')
args = options.parse_args()

sampleCandidates = args.candidates
outFile = args.outFile
comparison = args.comparisonPair
minReps = args.minReps
groupList = args.groupList.split(",")
compareSamp=f'{groupList[0]}_vs_{groupList[1]}'

## import group as dictionary from command line
sampleFile = args.sampleFile
sampleFile = pd.read_csv(sampleFile, sep="\t", comment="#", na_filter=False)

groups = {}
for group in groupList:
    samplesList = sampleFile.Name[sampleFile.Group == group].tolist()
    groups[group]=samplesList


## are there replicates
n = 1
for key in groups:
    if len(groups[key]) > n:
        n = len(groups[key])

if n > 1:
    reps = True
else:
    reps = False

sampleNames = []

for tsv in sampleCandidates:
    sampName = tsv.split("/")[-3]
    sampleNames.append(sampName)

print(sampleNames)

## get delta methylation rate
def deltaMethRate(row):
    
    methRates_by_group = []
    
    ## get average of groups
    for group in groups:
        methRates = []
        for sample in groups[group]:
            methRates.append(row[f"methRate_{sample}"])
        
        arr = np.array(methRates)
        methRates_by_group.append(np.mean(arr))
        
    return methRates_by_group[1] - methRates_by_group[0]

## methRateFC
def methRateFoldChange(row):
    methRates_by_group = []
    
    ## get average of groups
    for group in groups:
        methRates = []
        for sample in groups[group]:
            methRates.append( math.log2(row[f"methRate_{sample}"] + 0.0001) )
        
        arr = np.array(methRates)
        methRates_by_group.append(np.mean(arr))
        
    return methRates_by_group[1] - methRates_by_group[0]    

## function for fisher's exact test
def fet(row):
    tl = row[('cov_' + sampleNames[1])] - row[('C_count_' + sampleNames[1])]
    tr = row[('C_count_' + sampleNames[1])]
    bl = row[('cov_' + sampleNames[0])] - row[('C_count_' + sampleNames[0])]
    br = row[('C_count_' + sampleNames[0])]
    data = [ [tl, tr],
             [bl, br] ]
    
    if (math.isnan(tr)) or (math.isnan(br)):
        p_value = float('nan')
    else:
        odd_ratio, p_value = stats.fisher_exact(data)
    return p_value

## add combined p-value column for each group
def combinedPval(row):
            
    ## define list of pVals within a group
    listOfPvals = []

    ## loop through all sample groups
    sampleList = groups[group]
    for sample in sampleList:
        pVal = row[('p-value_mState_' + sample)]
        listOfPvals.append(pVal)

    ## calculate combined pVal
    testStatistic, combined_pVal = stats.combine_pvalues(listOfPvals, method='fisher',weights=None)
    return combined_pVal

## add cochran mantel haenszel test 
def cmhTest(row):
    ## need to make a list of lists, where the nested list contains the following:
    ##         [a, b, c, d], where a=sample1_unconvertedBaseCount, b=sample1_coverage, 
    ##                             c=sample2_unconvertedBaseCount, d=sample2_coverage

    ## define matrix for final p-value calculation
    mtx = []
    
    ## get replicates; if only one replicate in control, make it zero coverage.
    for i in range(len(next(iter(groups)))):
        oneCtrl = False
        firstReplicatesList = []
        
        for group in groups:
            
            if i > len(groups[group])-1:
                oneCtrl = True
            else:
                firstReplicate = groups[group][i]
            firstReplicatesList.append(firstReplicate)

        if oneCtrl == False:
            a = row[(f"cov_{firstReplicatesList[1]}")] - row[(f"C_count_{firstReplicatesList[1]}")]
            b = row[(f"C_count_{sampleNames[1]}")]
            c = row[(f"cov_{firstReplicatesList[0]}")] - row[(f"C_count_{firstReplicatesList[0]}")]
            d = row[(f"C_count_{sampleNames[0]}")]
        else:
            a = row[(f"cov_{firstReplicatesList[1]}")] - row[(f"C_count_{firstReplicatesList[1]}")]
            b = row[(f"C_count_{sampleNames[1]}")]
            c = 0
            d = 0        
    
        repComparison = [a,b,c,d]
        if (row[(f"cov_{firstReplicatesList[1]}")] != float("nan")) and (row[(f"cov_{firstReplicatesList[0]}")] != float("nan")):
            mtx.append(repComparison)
        
    mtx = np.asarray(mtx)
    
    tables = [np.reshape(x.tolist(), (2, 2)) for x in mtx]
    # print(tables)

    cmh = sm.stats.StratifiedTable(tables)
    pval = float(cmh.test_null_odds().pvalue)

    return pval

## add column for how many replicates this candidate was found in for each group
def countRepsPair(row):
    nonNANcount = 0
    ## go through first group
    for replicate in groups[group]:
        currentCov = row[(f'cov_{replicate}')]

        if math.isnan(float(currentCov)) == False:
            nonNANcount += 1

    return nonNANcount


## Make a list of DFs
dfList = []

## modify format of sample files to merge them easily
for i in range(len(sampleCandidates)):
    dfTmp = pd.read_csv(sampleCandidates[i], sep = '\t')
    
    cols = ['#SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene']
    dfTmp['SeqID'] = dfTmp[cols].apply(lambda row: ';'.join(row.values.astype(str)), axis=1)

    ## rename columns by adding the sample name as an identifier
    tmpCol = dfTmp.pop('SeqID')
    dfTmp.insert(0, 'SeqID', tmpCol)  
    dfTmp = dfTmp.drop(columns=cols)  
    dfTmp = dfTmp.rename(columns={'cov': 'cov_'+sampleNames[i],
                                 'C_count': 'C_count_'+sampleNames[i],
                                 'methRate': 'methRate_'+sampleNames[i],
                                 'uniqCov': 'uniqCov'+sampleNames[i],
                                 'uniqCovRatio': 'uniqCovRatio'+sampleNames[i],
                                 'C_count_uniq': 'C_count_uniq'+sampleNames[i],
                                 'methRate_uniq': 'methRate_uniq_'+sampleNames[i],
                                 '95_CI_lower': '95_CI_lower_'+sampleNames[i],
                                 '95_CI_upper': '95_CI_upper_'+sampleNames[i],
                                 'p-value_mState': 'p-value_mState_'+sampleNames[i],
                                 'FDR_mState': 'FDR_mState_'+sampleNames[i],
                                 'scores': 'scores_'+sampleNames[i],
                                 'seqContext': 'seqContext_'+sampleNames[i],
                                 'Source': 'Source_'+sampleNames[i],
                                 'gene_type': 'gene_type_'+sampleNames[i],
                                 'gene': 'gene_'+sampleNames[i]})
    
    dfTmp = dfTmp.set_index('SeqID')
    dfList.append(dfTmp)

## merge sample files according to same loci
masterDF = pd.concat(dfList, axis=1)


## for each group, order samples by from high to low total coverage
for group in groups:
    
    samps = []
    totalCovs = []
    
    for samp in groups[group]:
        samps.append(samp)
        totalCovs.append( masterDF[f'cov_{samp}'].sum() )
    
    samplesSortedByHighCov = [x for _,x in sorted(zip(totalCovs,samps))][::-1]    
    groups[group] = samplesSortedByHighCov


## compare
## get delta methylation rate
masterDF[f"delta_MethRate_{compareSamp}"] = masterDF.apply(deltaMethRate, axis=1)

## get MR fold change
masterDF[f"log2_methRateFC_{compareSamp}"] = masterDF.apply(methRateFoldChange, axis=1)

## get absolute values of the following columns
absDeltaMethRate = list(map(abs, masterDF[f"delta_MethRate_{compareSamp}"]))
masterDF[f"abs_delta_methRate_{compareSamp}"] = absDeltaMethRate
absMethRateFC = list(map(abs, masterDF[f"log2_methRateFC_{compareSamp}"]))
masterDF[f"abs_log2_methRateFC_{compareSamp}"] = absMethRateFC

if reps == True:
    # get combined p-value for each group
    for group in groups:
        masterDF[f"numReps_{group}"] = masterDF.apply(countRepsPair, axis=1)
        masterDF[f"combined_p-value_mState_{group}"] = masterDF.apply(combinedPval, axis=1)
    ## get cochran-mantel-haenzsel pval column
    masterDF[f"p-value_mState_{compareSamp}"] = masterDF.apply(cmhTest, axis=1)
else:
    # skip combined pval since no replicates
    for group in groups:
        masterDF[f"combined_p-value_mState_{group}"] = np.nan
    ## get fisher's exact test column
    masterDF[f"p-value_mState_{compareSamp}"] = masterDF.apply(fet, axis=1)

## get FDR values
pVals = np.array(masterDF[f"p-value_mState_{compareSamp}"])
rankedP = stats.rankdata(pVals)
fdr = pVals * len(pVals) / rankedP
fdr[fdr > 1] = 1
## add FDR values to DF
masterDF[f"FDR_mState_{compareSamp}"] = fdr

## get combined FDR
for group in groups:
    pVals = np.array(masterDF[f"combined_p-value_mState_{group}"])
    rankedP = stats.rankdata(pVals)
    fdr = pVals * len(pVals) / rankedP
    fdr[fdr > 1] = 1
    ## add FDR values to DF
    masterDF[f"combined_FDR_mState_{group}"] = fdr


## formatting index
masterDF.reset_index(inplace=True)
masterDF = masterDF.rename(columns = {'index':'SeqID'})

## splitting identifier to user-readable format
masterDF[['#SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene']] = masterDF["SeqID"].str.split(";", expand=True)

newCols = ['#SeqID', 'refPos', 'refStrand', 'refBase', 'seqContext', 'genomicCoords', 'Source', 'gene_type', 'gene', ]
masterDF = masterDF[newCols + [col for col in masterDF.columns if col not in newCols ] ]
masterDF = masterDF.drop(columns=['SeqID'])

masterDF.to_csv(outFile + "/pairwise_comparison.tsv", index=False, sep="\t", header=True)

