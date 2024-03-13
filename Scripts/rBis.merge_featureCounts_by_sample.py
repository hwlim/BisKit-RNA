#!/usr/bin/env python3

import argparse
import pandas as pd

options = argparse.ArgumentParser(description="Create merged table of featureCounts results for all sources of a sample.", usage="python3 rBis.merge_featureCounts_by_sample.py [options] -s sample1 -f tRNA.tsv miRNA.tsv circRNA.tsv ... ")
options.add_argument('-f','--featureCounts', nargs='+',
                        help='Required; list of featureCounts output for both all and uniq bam files.', required=True)
options.add_argument('-o', '--outFile', default='out.tsv',
                        help='Path to output file')
options.add_argument('-s', '--sampleName', required=True, default='sampleName',
                        help='Sample Name')
args = options.parse_args()

featureCountList = args.featureCounts
outFile = args.outFile
sampleName = args.sampleName

tRNA_all = pd.read_csv(featureCountList[0], comment='#', sep="\t", header=0)
tRNA_all.columns.values[2] = 'Count_all'
tRNA_uniq = pd.read_csv(featureCountList[1], comment='#', sep="\t", header=0)
tRNA_uniq.columns.values[2] = 'Count_uniq'
tRNA = pd.merge(tRNA_uniq, tRNA_all, on=['Gene_id', 'Length'])
tRNA.drop(tRNA.columns[1], axis=1, inplace=True)
tRNA["Gene_id"] = tRNA["Gene_id"].apply(lambda x: 'tRNA:' + x)

miRNA_all = pd.read_csv(featureCountList[2], comment='#', sep="\t", header=0)
miRNA_all.columns.values[2] = 'Count_all'
miRNA_uniq = pd.read_csv(featureCountList[3], comment='#', sep="\t", header=0)
miRNA_uniq.columns.values[2] = 'Count_uniq'
miRNA = pd.merge( miRNA_uniq, miRNA_all, on=['Gene_id', 'Length'])
miRNA.drop(miRNA.columns[1], axis=1, inplace=True)
miRNA["Gene_id"] = miRNA["Gene_id"].apply(lambda x: 'miRNA:' + x)

piRNA_all = pd.read_csv(featureCountList[4], comment='#', sep="\t", header=0)
piRNA_all.columns.values[2] = 'Count_all'
piRNA_uniq = pd.read_csv(featureCountList[5], comment='#', sep="\t", header=0)
piRNA_uniq.columns.values[2] = 'Count_uniq'
piRNA = pd.merge(piRNA_uniq, piRNA_all, on=['Gene_id', 'Length'])
piRNA.drop(piRNA.columns[1], axis=1, inplace=True)
piRNA["Gene_id"] = piRNA["Gene_id"].apply(lambda x: 'piRNA:' + x)

genome_all = pd.read_csv(featureCountList[6], comment='#', sep="\t", header=0)
genome_all.columns.values[2] = 'Count_all'
genome_uniq = pd.read_csv(featureCountList[7], comment='#', sep="\t", header=0)
genome_uniq.columns.values[2] = 'Count_uniq'
genome = pd.merge(genome_uniq, genome_all, on=['Gene_id', 'Length'])
genome.drop(genome.columns[1], axis=1, inplace=True)
genome["Gene_id"] = genome["Gene_id"].apply(lambda x: 'genome:' + x)

circRNA_all = pd.read_csv(featureCountList[8], comment='#', sep="\t", header=0)
circRNA_all.columns.values[2] = 'Count_all'
circRNA_uniq = pd.read_csv(featureCountList[9], comment='#', sep="\t", header=0)
circRNA_uniq.columns.values[2] = 'Count_uniq'
circRNA = pd.merge(circRNA_uniq, circRNA_all, on=['Gene_id', 'Length'])
circRNA.drop(circRNA.columns[1], axis=1, inplace=True)
circRNA["Gene_id"] = circRNA["Gene_id"].apply(lambda x: 'circRNA:' + x)

df = pd.concat([tRNA, miRNA, piRNA, genome, circRNA], axis=0)

df.rename(columns={"Gene_id": "source:gene", "Count_uniq": f"{sampleName}_uniq", "Count_all": f"{sampleName}_all"}, inplace=True)

df.to_csv(outFile, sep='\t', index=False)