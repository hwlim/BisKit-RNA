#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse

options = argparse.ArgumentParser(description="Create bigwig files for candidate profile for each sample.", usage="python3 rBis.draw_candidate_profile.py [options] -l lookupTable.tsv")
options.add_argument('-l', '--lookupTable', required=True,
                        help='path to lookupTable.tsv')  
options.add_argument('-o', '--outFile', default='~',
                        help='Path to output compare.tsv file')
options.add_argument('-c', '--candidates',
                        help='mergedCandidates.tsv')        
args = options.parse_args()

outFile = args.outFile
mergedCandidatesIn = args.candidates

## read input files
lookupTable = pd.read_csv(args.lookupTable, sep = '\t', header=0)
mergedCandidatesIn = pd.read_csv(mergedCandidatesIn, sep = '\t', header=0)

## get rid of rRNA
filtered = mergedCandidatesIn.loc[(mergedCandidatesIn['Source'] != "rRNA")]

masterList = []

for index, row in filtered.iterrows():
    seqID = row['#SeqID'].split("|")[0]
    pos = row['refPos']
    mr = row['methRate']
    refStrand = row['refStrand']

    if row["Source"] == "tRNA":
        source = 0
        chrom = row['genomicCoords'].split(":")[0]
        strand = row['genomicCoords'].split(":")[2]
        start = row['genomicCoords'].split(":")[1].split("-")[0]
        end = row['genomicCoords'].split(":")[1].split("-")[1]

    elif row["Source"] == "miRNA":
        source = 1
        chrom = row['genomicCoords'].split(":")[0]
        strand = row['genomicCoords'].split(":")[2]
        start = row['genomicCoords'].split(":")[1].split("-")[0]
        end = row['genomicCoords'].split(":")[1].split("-")[1]

    elif row["Source"] == "piRNA":
        source = 2
        chrom = row['genomicCoords'].split(":")[0]
        strand = row['genomicCoords'].split(":")[2]
        start = row['genomicCoords'].split(":")[1].split("-")[0]
        end = row['genomicCoords'].split(":")[1].split("-")[1]

    elif row["Source"] == "Genome":
        source = 3
        chrom = seqID
        start = pos - 1
        end = start + 1
        strand = row['refStrand']

    elif row["Source"] == "circRNA":
        source = 4
        chrom = row['genomicCoords'].split(":")[0]
        strand = row['genomicCoords'].split(":")[2]
        start = row['genomicCoords'].split(":")[1].split("-")[0]
        end = row['genomicCoords'].split(":")[1].split("-")[1]

    else:
        print(row['Source'])
        exit()

    if chrom == "chrMT":
        continue
    elif chrom.startswith("chrJH"):
        continue

    masterList.append([chrom, start, end, mr, strand, source])

plusList = []
minusList = []

## filter by strand
for lst in masterList:
    if (lst[-2] == "+"):
        plusList.append(lst)

    elif (lst[-2] == "-"):
        minusList.append(lst)

        
plusDF = pd.DataFrame(plusList, columns = ['chr', 'start', 'end', 'mr', 'strand', 'source'])
minusDF = pd.DataFrame(minusList, columns = ['chr', 'start', 'end', 'mr', 'strand', 'source'])
plusDF = plusDF.drop('strand', axis=1)
minusDF = minusDF.drop('strand', axis=1)
plusDF['start'] = plusDF['start'].astype(int)
minusDF['start'] = minusDF['start'].astype(int)

## sort by chromosome, start, end, and source (source is ordered by importance from tRNA, miRNA, to Genome)
sortedPlusDF = plusDF.sort_values(by = ['chr', 'start', 'source'], ascending = [True, True, True], na_position = 'first')
sortedMinusDF = minusDF.sort_values(by = ['chr', 'start', 'source'], ascending = [True, True, True], na_position = 'first')

dups_color_and_shape = sortedPlusDF.pivot_table(columns=['chr','start'], aggfunc='size')
print("Total plus candidates:", len(sortedPlusDF))
print("Duplicate plus loci:", list(dups_color_and_shape).count(2))

dups_color_and_shape = sortedMinusDF.pivot_table(columns=['chr','start'], aggfunc='size')
print("Total minus candidates:", len(sortedMinusDF))
print("Duplicate minus loci:", list(dups_color_and_shape).count(2))

## drop duplicates, dropping lower importance source first
## then drop last column what has information on source, since bedGraphToBigWig can't read 5 columns
uniqPlusDF = sortedPlusDF.drop_duplicates(
    subset = ['chr', 'start', 'end'],
    keep = 'first').reset_index(drop = True)

uniqPlusDF = uniqPlusDF.drop('source', axis=1)
uniqPlusDF.drop_duplicates(subset=['chr', 'start'], keep='first', inplace=True)

uniqMinusDF = sortedMinusDF.drop_duplicates(
  subset = ['chr', 'start', 'end'],
  keep = 'first').reset_index(drop = True)

uniqMinusDF = uniqMinusDF.drop('source', axis=1)
uniqMinusDF.drop_duplicates(subset=['chr', 'start'], keep='first', inplace=True)

uniqPlusDF.to_csv(f'{outFile}/plusCandidates.bg', index=False, sep="\t", header=None)
uniqMinusDF.to_csv(f'{outFile}/minusCandidates.bg', index=False, sep="\t", header=None)
