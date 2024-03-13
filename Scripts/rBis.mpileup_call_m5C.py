#!/usr/bin/env python

import argparse
import os
import re
import sys
import pybedtools
import pandas as pd
import numpy as np
import scipy.stats as stats
import math

pd.options.display.max_colwidth = 200000
pd.set_option('display.max_rows', 500)
from Bio import SeqIO

#parse command line arguments
options = argparse.ArgumentParser(description="Call m5C candidates using the output from mpileup.", usage="python mpileup_call_m5C.py (options) -c mpileup_final.tsv -m rRNA -r genome.fa -e sample1 -u mm10 ")
options.add_argument('-c', '--call', required=True,
                        help='mpileup_final.tsv file')
options.add_argument('-z', '--z', default=1.96,
                        help='z-val for 0.95 confidence interval; default = 1.96')
options.add_argument('-q', '--q', default=30,
                        help='minimum base quality threshold; default = 30')
options.add_argument('-v', '--minCov', default=10,
                        help='create an additional output file with minimum coverage; default = 10')
options.add_argument('-o', '--outFile', default="call.tsv",
                        help='output file destination.')
options.add_argument('-s', '--sig', default=0.05,
                        help='Significance threshold; default = 0.05')
options.add_argument('-t', '--sigType', default="pVal",
                        help='Significance type; default = pVal; use pVal or FDR')
options.add_argument('-m', '--source', required=True,
                        help='input RNA type; e.g. rRNA, tRNA, miRNA, piRNA, genome, circRNA')
options.add_argument('-g', '--annotations', default="/data/limlab/Resource/exoBS/mm10/miRNA/mmu.gff3",
                        help='path to annotation file; miRNA.gff3 for miRNA and tRNA.fa for tRNA.')
options.add_argument('-l', '--lookupTable', default="/data/limlab/Resource/exoBS/mm10/lookupTable.tsv",
                        help='path to lookup table; default="/data/limlab/Resource/exoBS/mm10/lookupTable.tsv".')
options.add_argument('-w', '--seqContextWindow', default=3,
                        help='sequence context window size; default = 3, 3 on each side of candidate C')
options.add_argument('-r', '--genomeRefFA', required=True,
                        help='path to genome reference fasta file')
options.add_argument('-a', '--rRNArefFA', default="/data/limlab/Resource/exoBS/mm10/rRNA/rRNA.fa",
                        help='path to rRNA reference fasta file')
options.add_argument('-n', '--minCcount', default=3,
                        help='min C count required at a C position for it to be considered a candidate; default = 3')
options.add_argument('-e', '--sampleName', required=True,
                        help='Sample Name')
options.add_argument('-u', '--species', required=True,
                        help='species; hg38 or mm10')
options.add_argument('-y', '--mpileupIn', default="",
                        help='mpileup output file used for genome only')

args = options.parse_args()
inCall = args.call
z = args.z
q = args.q
sigThresh = float(args.sig)
exp_meth = int(10**(-q/10))
outFilePrefix = args.outFile
minCov = int(args.minCov)
source = args.source
rRNAFA = args.rRNArefFA
genomeFA = args.genomeRefFA
gff3in = args.annotations
seqContextWindow = int(args.seqContextWindow)
lookupTableIn = args.lookupTable
minC = args.minCcount
sampleName = args.sampleName
species = args.species
mpileup_in = args.mpileupIn

if args.sigType == "pVal":
    sigType = "p-value_mState"
    sigSave = "pVal"
elif args.sigType == "FDR":
    sigType = "FDR_mState"
    sigSave = "FDR"
else:
    print("need to use pVal or FDR as -t flag.")
    exit()

## define tmpdir
tempdir = f'{os.environ.get("TMPDIR")}/pybedtools/{sampleName}_{source}'
if not os.path.exists(tempdir):
    os.makedirs(tempdir)
pybedtools.helpers.set_tempdir(tempdir)


def getGene(inRow):
    rowList = inRow.split(";")
    if rowList[0] == ".":
        return "no_feature"
    else:
        txt = "gene_name"
    #     print(rowList)
        target = [i for i in rowList if txt in i][0]
        gene = target.replace('"','').split(" ")[-1]
        return gene

def getGeneType(inRow):
    rowList = inRow.split(";")
    if rowList[0] == ".":
        return "no_gene_type"
    else:
        txt = "gene_type"
    #     print(rowList)
        target = [i for i in rowList if txt in i][0]
        gene = target.replace('"','').split(" ")[-1]
        return gene

## Wilson confidence interval
def conf_interval(cCount, cov):
    p = cCount / cov
    denom = 1 + z**2/cov
    cap = p + z*z / (2*cov)
    asd = math.sqrt((p*(1 - p) + z*z / (4*cov)) / cov)
    lower = (cap - z*asd) / denom
    upper = (cap + z*asd) / denom
    
    return [lower, upper]

def pVal_mState(unconvertedBaseCount, cov):
    exp_cov = cov - exp_meth
    data = [[unconvertedBaseCount, cov - unconvertedBaseCount], [exp_meth, exp_cov]]

    odd_ratio, p_value = stats.fisher_exact(data)
    return p_value

def uniqCovRatio(row):
    return row["uniqCov"] / row["cov"]

def getReverseCompliment(context):
    newContext = ""
    for i in range(len(context)):
        if context[i] == "A":
            letter = "T"
        elif context[i] == "a":
            letter = "t"
        elif context[i] == "C":
            letter = "G"
        elif context[i] == "c":
            letter = "g"
        elif context[i] == "G":
            letter = "C"
        elif context[i] == "g":
            letter = "c"
        elif context[i] == "T":
            letter = "A"
        elif context[i] == "t":
            letter = "a"
        elif context[i] == "N":
            letter = "N"
        
        newContext += letter

    return newContext[::-1]

def getCov(row):
    unconvertedBaseCount = row['C_count']
    convertedBaseCount = row['convertedBaseCount']
    return unconvertedBaseCount + convertedBaseCount

def getUniqCov(row):
    unconvertedBaseCount = row['C_count_uniq']
    convertedBaseCount = row['convertedBaseCount_uniq']
    if (math.isnan(unconvertedBaseCount)) or (math.isnan(convertedBaseCount)):
        return 0
    else:
        return int(unconvertedBaseCount + convertedBaseCount)

def getMethRate_all(row):
    unconvertedBaseCount = row['C_count']
    cov = row['cov']
    return unconvertedBaseCount / cov

def getMethRate_uniq(row):
    unconvertedBaseCount = row['C_count_uniq']
    cov = row['uniqCov']
    if np.isnan(unconvertedBaseCount):
        return 0
    elif (np.isnan(cov)) or (cov == 0):
        return 0
    else:
        return unconvertedBaseCount / cov

def getConfsAndScores(row):
    unconvertedBaseCount = row['C_count']
    cov = row['cov']
    confs = conf_interval(unconvertedBaseCount, cov)
    pVal = row['p-value_mState']
    score = unconvertedBaseCount * confs[0] * -1 * math.log10( pVal+0.0001 )
    return f'{confs[0]};{confs[1]};{score}'

def getpVal(row):
    unconvertedBaseCount = row['C_count']
    cov = row['cov']
    return pVal_mState(unconvertedBaseCount, cov)

def getRefBase(row):
    if row['refStrand'] == "+":
        return "C"
    else:
        return "G"

def getSeqContextrRNA(row):
    numNs = 0
    chrom = row['#SeqID']
    seqContextStart = row['refPos'] - seqContextWindow - 1
    seqContextEnd = row['refPos'] + seqContextWindow
    if seqContextStart < 0:
        numNs = abs(seqContextStart)
        seqContextStart = 0
    
    if species == "mm10":
        if chrom == "BK000964.3":
            refLen = refLengths[0]
        else:
            refLen = refLengths[1]
    elif species == "hg38":
        refLen = refLengths[0]
    
    if seqContextEnd > refLen:
        numNs = seqContextEnd - refLen
        seqContextEnd = refLen
        
    a = pybedtools.BedTool(f'{chrom}\t{seqContextStart}\t{seqContextEnd}', from_string=True)
    b = a.sequence(fi=rRNAFA)
    seq = str(open(b.seqfn).read())
    context = seq.split("\n")[1]
    
    addNs = "N" * numNs
    
    if seqContextStart == 0:
        return addNs + context
    else:
        return context + addNs

def getGenomicCoords(row):
    if row.empty:
        return "0"
    
    strand = row['refStrand']
    transcriptName = row['#SeqID']
    relativePos = row['refPos']
    ID = transcriptName.split("|")[0]
    geneName = transcriptName.split("|")[0].split("_")[-1]
    seqID = geneName

    ## get genomic coords
    rowLookup = lookupTable[lookupTable['seqID'] == ID]
    
    ## ignore if no annotation exists
    if len(rowLookup) == 0:
        return "0"

    location = rowLookup['coord'].values[0].split(",")[0]
    chrom = location.split(":")[0]    
    currentStrand = location.split(":")[2]
    
    if currentStrand == "+":
        start = int(location.split(":")[1].split("-")[0]) + relativePos - 2
        end = start + 1
        genomicCoord = f'{chrom}:{start}-{end}:{currentStrand}'
    
    else:
        end = int(location.split(":")[1].split("-")[1]) - relativePos + 1
        start = end - 1
        genomicCoord = f'{chrom}:{start}-{end}:{currentStrand}'        
        
    return genomicCoord

def getSeqContext(row):
    if row.empty:
        return ""
    
    genomicCoord = row['genomicCoords']

    if genomicCoord == "" or genomicCoord == "0":
        return ""
    else:
        chrom = genomicCoord.split(":")[0]
        start = genomicCoord.split(":")[1].split("-")[0]
        end = genomicCoord.split(":")[1].split("-")[1]
        currentStrand = genomicCoord.split(":")[-1]
        relativePos = row['refPos']

        if currentStrand == "+":
            seqContextStart = int(start) - seqContextWindow
            seqContextEnd = int(start) + seqContextWindow + 1
        else:
            seqContextStart = int(start) - seqContextWindow
            seqContextEnd = int(start) + seqContextWindow + 1


        a = pybedtools.BedTool(f'{chrom}\t{seqContextStart}\t{seqContextEnd}\t.\t.\t{currentStrand}', from_string=True)
        b = a.sequence(fi=genomeFA, s=True)
        seq = str(open(b.seqfn).read())
        try:
            context = seq.split("\n")[1]
        except IndexError:
            return ""
        
        return context

def getGeneNameTRNA(row):
    transcriptName = row['#SeqID']
    relativePos = row['refPos']
    ID = transcriptName.split("|")[0]
    geneName = transcriptName.split("|")[0].split("_")[-1]
    return geneName

def getGeneName_piRNA(row):
    transcriptName = row['#SeqID']
    relativePos = row['refPos']
    ID = transcriptName.split("|")[0]
    geneName = transcriptName.split("|")[0].split("_")[-1]
    return geneName

def getGeneName_circRNA(row):
    transcriptName = row['#SeqID']
    return transcriptName

def getGeneType_miRNA(row):
    if row.empty:
        return ""
    
    info = row[0].split("|")
    ID = info[1]
    relativePos = row['refPos']
    seqID = info[0]
    rowLookup = gff3[gff3[8].str.contains(ID)]

    if len(rowLookup) == 0:
        geneType = "-"
    elif len(rowLookup) == 1:
        geneType = "hairpin"
    elif len(rowLookup) > 1:
        ## get genome location from position
        geneName = rowLookup[8].values[0].split(";")[2].split("=")[-1]
        geneStart = rowLookup[3].values[0]
        start = int(geneStart) + relativePos - 1
        geneType = "hairpin"

        i = 1
        while i < len(rowLookup.values):

            ## get start and end of mature miRNA coords
            matureSeqID = rowLookup[8].values[i].split(";")[2].split("=")[-1]
            matureStart = int(rowLookup[3].values[i])
            matureEnd = int(rowLookup[4].values[i])

            ## check if current C falls in mature miRNA region
            if (start > matureStart) and (start <= matureEnd):

                geneName = matureSeqID
                geneType = "mature"
                break

            i += 1
            
    return geneType

def getGeneName_miRNA(row):
    if row.empty:
        return ""
    
    info = row[0].split("|")
    ID = info[1]
    relativePos = row['refPos']
    seqID = info[0]
    rowLookup = gff3[gff3[8].str.contains(ID)]

    if len(rowLookup) == 0:
        geneName = "-"
        geneType = "-"
    elif len(rowLookup) == 1:
        geneName = rowLookup[8].values[0].split(";")[2].split("=")[-1]
        geneType = "hairpin"
    elif len(rowLookup) > 1:
        ## get genome location from position
        geneName = rowLookup[8].values[0].split(";")[2].split("=")[-1]
        geneStart = rowLookup[3].values[0]
        start = int(geneStart) + relativePos - 1
        geneType = "hairpin"

        i = 1
        while i < len(rowLookup.values):

            ## get start and end of mature miRNA coords
            matureSeqID = rowLookup[8].values[i].split(";")[2].split("=")[-1]
            matureStart = int(rowLookup[3].values[i])
            matureEnd = int(rowLookup[4].values[i])

            ## check if current C falls in mature miRNA region
            if (start > matureStart) and (start <= matureEnd):

                geneName = matureSeqID
                geneType = "mature"
                break

            i += 1
            
    return geneName

def getID_miRNA(row):
    info = row[0].split("|")[0]
    return info


## Genome
if source == "genome":
    chrList = []
    posList = []
    strandList = []
    convertedCountListAll = []
    unconvertedCountListAll = []
    convertedCountListUniq = []
    unconvertedCountListUniq = []
    geneList = []
    miniGeneList = []
    geneTypeList = []
    miniGeneTypeList = []
    seqContexts = []
    genomicCoords = []
 
    #read bed file
    annotateBed = pd.read_csv(inCall, sep='\t', header=None)
    annotateBed = annotateBed.drop([1,6,7,8,9,10,11,12,13], axis=1)

    previousChr = ''
    previousPos = 0
    previousStrand = ''
    currentRow = 0

    ## loop through dataframe
    for index, row in annotateBed.iterrows():
        currentChr = row[0]
        currentPos = row[2]
        currentStrand = row[3]
        convertedBaseCount = row[4]
        unconvertedBaseCount = row[5]
        feature = row[14]
        gene = getGene(feature)
        geneType = getGeneType(feature)
        
        if currentStrand == "-":
            genomicCoord = f'{currentChr}:{str(int(currentPos-1))}-{currentPos}:{currentStrand}'
        else:
            genomicCoord = f'{currentChr}:{str(int(currentPos-1))}-{currentPos}:{currentStrand}'


        ## if duplicate row and only gene is different, just add the gene
        if (currentChr == previousChr) and (currentPos == previousPos) and (currentStrand == previousStrand):
            miniGeneList.append(gene)
            miniGeneTypeList.append(geneType)

        ## if not a duplicate row
        else:
            miniGeneList.append(gene)
            miniGeneTypeList.append(geneType)

            chrList.append(currentChr)
            posList.append(currentPos)
            strandList.append(currentStrand)
            convertedCountListAll.append(convertedBaseCount)
            unconvertedCountListAll.append(unconvertedBaseCount)
            geneList.append(",".join(miniGeneList))
            geneTypeList.append(",".join(miniGeneTypeList))
            genomicCoords.append(genomicCoord)

            ## reset mini Gene list
            miniGeneList = []
            miniGeneTypeList = []

        previousChr = currentChr
        previousPos = currentPos
        previousStrand = currentStrand

    ## create dataframe from lists
    data = {'#SeqID': chrList,
            'refPos': posList,
            'refStrand': strandList,
            'convertedBaseCount': convertedCountListAll,
            'C_count': unconvertedCountListAll,
            'gene_type': geneTypeList,
            'gene': geneList,
            'genomicCoords': genomicCoords}
    df = pd.DataFrame(data)

else:
    callBed = pd.read_csv(inCall, sep='\t', header=0)
    if callBed.empty:
        
        rearrange = ["#SeqID", "refPos", "refStrand", "refBase", "cov", "C_count", "methRate", "uniqCov", "C_count_uniq", "methRate_uniq", "95_CI_lower", "95_CI_upper", "p-value_mState", "FDR_mState", "scores", "seqContext", "genomicCoords", "gene_type", "gene"]
        newDF = pd.DataFrame(columns=rearrange)
        newDF.to_csv(outFilePrefix + "/call.tsv", index=False, sep = "\t", header=True)

        cov10DF = newDF[(newDF['cov'] >= int(minCov))]

        numCandidates = 0
        numCandidates_cov10 = 0
        sigFilterDF = 0
        mrFilterDF = 0
        significantDF = 0
        sigRatio = 0

        data = {"Sample": [sampleName],
                "num_Candidates_total": [numCandidates],
                "num_Candidates_cov>=10": [numCandidates_cov10],
                "num_Candidates_MR>=0.1": [mrFilterDF],
                "num_Candidates_" + sigSave + "<" + str(sigThresh): [sigFilterDF],
                "num_Candidates_cov_MR_"+ sigSave +"_filtered": [significantDF]}

        table = pd.DataFrame(data)
        table.to_csv(outFilePrefix + "/call_stats.tsv", index=False, sep = "\t", header=True)

        simplifiedData = {"Sample": [sampleName],
                        "num_sig_candidates": [significantDF],
                        "num_non-sig_candidates": [numCandidates_cov10 - significantDF]}

        simplifiedTable = pd.DataFrame(simplifiedData)
        simplifiedTable.to_csv(outFilePrefix + "/call_stats_simplified.tsv", index=False, sep = "\t", header=True)
        exit()

if source == "genome":
    mpileup = pd.read_csv(mpileup_in, sep = '\t', header=0)
    mpileup = mpileup.rename(columns={"ref": "#SeqID", "pos": "refPos", "strand":"refStrand", "convertedBaseCount":"convertedBaseCount", "unconvertedBaseCount":"unconvertedBaseCount", "uniqConvertedBaseCount":"convertedBaseCount_uniq", "uniqUnconvertedBaseCount":"C_count_uniq" })

    ## merge annotation and mpileup result
    merged = df.merge(mpileup, on=['#SeqID', 'refPos', 'refStrand', 'convertedBaseCount'])
    newDF = merged.drop(['unconvertedBaseCount'], axis=1)

elif source == "rRNA":
    ##get length of reference sequences in rRNA
    fasta_sequences = SeqIO.parse(open(rRNAFA),'fasta')
    refLengths = []
    for fasta in fasta_sequences:
        refLengths.append(len(str(fasta.seq)))
    
    callBed['gene_type'] = "rRNA"
    callBed['gene'] = ""
    callBed['genomicCoords'] = ""
    newDF = callBed.rename(columns={"ref": "#SeqID", "pos": "refPos", "strand":"refStrand", "convertedBaseCount":"convertedBaseCount", "unconvertedBaseCount":"C_count", "uniqConvertedBaseCount":"convertedBaseCount_uniq", "uniqUnconvertedBaseCount":"C_count_uniq", "gene_type":"gene_type", "gene":"gene", "genomicCoords":"genomicCoords"})

elif source == "tRNA":
    lookupTable = pd.read_csv(lookupTableIn, sep = '\t', header=0)
    callBed['gene_type'] = "tRNA"
    newDF = callBed.rename(columns={"ref": "#SeqID", "pos": "refPos", "strand":"refStrand", "convertedBaseCount":"convertedBaseCount", "unconvertedBaseCount":"C_count", "uniqConvertedBaseCount":"convertedBaseCount_uniq", "uniqUnconvertedBaseCount":"C_count_uniq"})
    
elif source == "piRNA":
    lookupTable = pd.read_csv(lookupTableIn, sep = '\t', header=0)
    callBed['gene_type'] = "piRNA"
    newDF = callBed.rename(columns={"ref": "#SeqID", "pos": "refPos", "strand":"refStrand", "convertedBaseCount":"convertedBaseCount", "unconvertedBaseCount":"C_count", "uniqConvertedBaseCount":"convertedBaseCount_uniq", "uniqUnconvertedBaseCount":"C_count_uniq", "gene_type":"gene_type"})
    
elif source == "circRNA":
    lookupTable = pd.read_csv(lookupTableIn, sep = '\t', header=0)
    callBed['gene_type'] = "circRNA"
    newDF = callBed.rename(columns={"ref": "#SeqID", "pos": "refPos", "strand":"refStrand", "convertedBaseCount":"convertedBaseCount", "unconvertedBaseCount":"C_count", "uniqConvertedBaseCount":"convertedBaseCount_uniq", "uniqUnconvertedBaseCount":"C_count_uniq", "gene_type":"gene_type"})
    
elif source == "miRNA":
    gff3 = pd.read_csv(gff3in, sep = '\t', skiprows=13, header=None)
    lookupTable = pd.read_csv(lookupTableIn, sep = '\t', header=0)
    callBed['#SeqID'] = callBed.apply(getID_miRNA, axis=1)
    newDF = callBed.rename(columns={"#SeqID": "#SeqID", "pos": "refPos", "strand":"refStrand", "convertedBaseCount":"convertedBaseCount", "unconvertedBaseCount":"C_count", "uniqConvertedBaseCount":"convertedBaseCount_uniq", "uniqUnconvertedBaseCount":"C_count_uniq", "gene_type":"gene_type"})

## get coverage
newDF['cov'] = newDF.apply(getCov, axis=1)
## get unique coverage
newDF['uniqCov'] = newDF.apply(getUniqCov, axis=1)
originalCovDF = newDF[(newDF['cov'] >= 1)]
newDF = newDF[(newDF['cov'] >= int(minCov))]
newDF = newDF[(newDF['C_count'] >= int(minC))]

if newDF.empty:
    rearrange = ["#SeqID", "refPos", "refStrand", "refBase", "cov", "C_count", "methRate", "uniqCov", "C_count_uniq", "methRate_uniq", "95_CI_lower", "95_CI_upper", "p-value_mState", "FDR_mState", "scores", "seqContext", "genomicCoords", "gene_type", "gene"]
    newDF = pd.DataFrame(columns=rearrange)
    newDF.to_csv(outFilePrefix + "/call.tsv", index=False, sep = "\t", header=True)

    cov10DF = newDF[(newDF['cov'] >= int(minCov))]

    numCandidates = 0
    numCandidates_cov10 = 0
    sigFilterDF = 0
    mrFilterDF = 0
    significantDF = 0
    sigRatio = 0

    data = {"Sample": [sampleName],
            "num_Candidates_total": [numCandidates],
            "num_Candidates_cov>=10": [numCandidates_cov10],
            "num_Candidates_MR>=0.1": [mrFilterDF],
            "num_Candidates_" + sigSave + "<" + str(sigThresh): [sigFilterDF],
            "num_Candidates_cov_MR_"+ sigSave +"_filtered": [significantDF]}

    table = pd.DataFrame(data)
    table.to_csv(outFilePrefix + "/call_stats.tsv", index=False, sep = "\t", header=True)

    simplifiedData = {"Sample": [sampleName],
                    "num_sig_candidates": [significantDF],
                    "num_non-sig_candidates": [numCandidates_cov10 - significantDF]}

    simplifiedTable = pd.DataFrame(simplifiedData)
    simplifiedTable.to_csv(outFilePrefix + "/call_stats_simplified.tsv", index=False, sep = "\t", header=True)
    exit()

if source == "rRNA":
    newDF['seqContext'] = newDF.apply(getSeqContextrRNA, axis=1)

elif source == "tRNA":
    newDF['gene'] = newDF.apply(getGeneNameTRNA, axis=1)
    newDF['#SeqID'] = newDF['gene']
    newDF['genomicCoords'] = newDF.apply(getGenomicCoords, axis=1)
    newDF['seqContext'] = newDF.apply(getSeqContext, axis=1)

elif source == "miRNA":
    newDF['gene_type'] = newDF.apply(getGeneType_miRNA, axis=1)
    newDF['gene'] = newDF.apply(getGeneName_miRNA, axis=1)
    newDF['genomicCoords'] = newDF.apply(getGenomicCoords, axis=1)
    newDF = newDF[(newDF['genomicCoords'] != "0")]
    newDF['seqContext'] = newDF.apply(getSeqContext, axis=1)

elif source == "piRNA":

    newDF['gene'] = newDF.apply(getGeneName_piRNA, axis=1)
    newDF['#SeqID'] = newDF['gene']
    newDF['genomicCoords'] = newDF.apply(getGenomicCoords, axis=1)
    newDF['seqContext'] = newDF.apply(getSeqContext, axis=1)

elif source == "circRNA":
    print("getting genes...")
    newDF['gene'] = newDF.apply(getGeneName_circRNA, axis=1)
    newDF['#SeqID'] = newDF['gene']
    print("getting genomic coords...")
    newDF['genomicCoords'] = newDF.apply(getGenomicCoords, axis=1)
    print("Getting seqContext for circRNA..")
    newDF['seqContext'] = newDF.apply(getSeqContext, axis=1)

elif source == "genome":
    newDF['seqContext'] = newDF.apply(getSeqContext, axis=1)

## get methRate etc
print("Getting methRate...")
newDF['methRate'] = newDF.apply(getMethRate_all, axis=1)
newDF['methRate_uniq'] = newDF.apply(getMethRate_uniq, axis=1)
newDF['p-value_mState'] = newDF.apply(getpVal, axis=1)
newDF['lower_upper_confs_scores'] = newDF.apply(getConfsAndScores, axis=1)

## get fdr and add to DF
pVals = np.array(newDF['p-value_mState'])
rankedP = stats.rankdata(pVals)
fdr = pVals * len(pVals) / rankedP
fdr[fdr > 1] = 1
newDF['FDR_mState'] = fdr

## add reference base
newDF['refBase'] = newDF.apply(getRefBase, axis=1)

## split confidence intervals and scores to separate cols
newDF[['95_CI_lower','95_CI_upper', 'scores']] = newDF['lower_upper_confs_scores'].str.split(';',expand=True)
rearrange = ["#SeqID", "refPos", "refStrand", "refBase", "cov", "C_count", "methRate", "uniqCov", "C_count_uniq", "methRate_uniq", "95_CI_lower", "95_CI_upper", "p-value_mState", "FDR_mState", "scores", "seqContext", "genomicCoords", "gene_type", "gene"]
newDF=newDF[rearrange]
newDF['C_count_uniq'] = newDF['C_count_uniq'].fillna(0)

newDF.to_csv(outFilePrefix + "/call.tsv", index=False, sep = "\t", header=True)

cov10DF = newDF[(newDF['cov'] >= int(minCov))]

numCandidates = len(originalCovDF)
numCandidates_cov10 = len(cov10DF)
sigFilterDF = len(cov10DF[(cov10DF[sigType] < sigThresh )])
mrFilterDF = len(cov10DF[(cov10DF['methRate'] >= 0.1 )])
significantDF = len(cov10DF[(cov10DF[sigType] < sigThresh ) & (cov10DF['methRate'] >= 0.1) ])
sigRatio = significantDF / numCandidates_cov10 * 100

data = {"Sample": [sampleName],
        "num_Candidates_total": [numCandidates],
        "num_Candidates_cov>=10": [numCandidates_cov10],
        "num_Candidates_MR>=0.1": [mrFilterDF],
        "num_Candidates_" + sigSave + "<" + str(sigThresh): [sigFilterDF],
        "num_Candidates_cov_MR_"+ sigSave +"_filtered": [significantDF]}

table = pd.DataFrame(data)
table.to_csv(outFilePrefix + "/call_stats.tsv", index=False, sep = "\t", header=True)

simplifiedData = {"Sample": [sampleName],
                  "num_sig_candidates": [significantDF],
                  "num_non-sig_candidates": [numCandidates_cov10 - significantDF]}

simplifiedTable = pd.DataFrame(simplifiedData)
simplifiedTable.to_csv(outFilePrefix + "/call_stats_simplified.tsv", index=False, sep = "\t", header=True)


## delete all temp files
pybedtools.cleanup(remove_all=True)
