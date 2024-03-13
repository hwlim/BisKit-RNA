#!/usr/bin/env python3

import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt

#parse command line arguments
options = argparse.ArgumentParser(description="draw m-bias plot using output bam file.", usage="python3 rBis.draw_mBias_plot.py [options] -r 50 -g 100 -b align.bam")
options.add_argument('-b', '--bamFile', required=True,
                        help='align.bam file from alignment')
options.add_argument('-o', '--outFilePrefix', default="m-bias",
                        help='output file destination and prefix.')
options.add_argument('-r', '--readLenCommon', required=True,
                        help='most common read length in the input bam file.')
options.add_argument('-g', '--readLenGenome', required=True,
                        help='most common read length in the genome_align.bam file.')
options.add_argument('-s', '--minQual', default=35,
                        help='minimum base quality to be considered for m-bias plot; default = 30')

args = options.parse_args()
bam_file = args.bamFile
outFilePrefix = args.outFilePrefix
readLenCommon = int(args.readLenCommon)
readLenGenome = int(args.readLenGenome)
minQual = int(args.minQual)

## define array of most common read length and count the C occurrences
mbias_forward_no_qual_C_common = np.zeros(readLenCommon)
mbias_forward_no_qual_T_common = np.zeros(readLenCommon)
mbias_forward_with_qual_C_common = np.zeros(readLenCommon)
mbias_forward_with_qual_T_common = np.zeros(readLenCommon)

mbias_reverse_no_qual_G_common = np.zeros(readLenCommon)
mbias_reverse_no_qual_A_common = np.zeros(readLenCommon)
mbias_reverse_with_qual_G_common = np.zeros(readLenCommon)
mbias_reverse_with_qual_A_common = np.zeros(readLenCommon)

mbias_forward_no_qual_C_genome = np.zeros(readLenGenome)
mbias_forward_no_qual_T_genome = np.zeros(readLenGenome)
mbias_forward_with_qual_C_genome = np.zeros(readLenGenome)
mbias_forward_with_qual_T_genome = np.zeros(readLenGenome)

mbias_reverse_no_qual_G_genome = np.zeros(readLenGenome)
mbias_reverse_no_qual_A_genome = np.zeros(readLenGenome)
mbias_reverse_with_qual_G_genome = np.zeros(readLenGenome)
mbias_reverse_with_qual_A_genome = np.zeros(readLenGenome)

## open bam file
bam = pysam.AlignmentFile(bam_file, "rb")

## number of reads used to plot mbias
numRelevantReadsForward_common = 0
numRelevantReadsReverse_common = 0
numRelevantReadsForward_genome = 0
numRelevantReadsReverse_genome = 0

## loop through all reads in bam file
for read in bam:

    ## skip invalid reads
    if read.is_unmapped or read.is_secondary:
        continue

    ## skip if read length is not equal to the relevant read length
    if read.query_length == readLenCommon:

        ## get sequence and qualities
        sequence = read.query_sequence
        qualities = read.query_qualities

        ## aggregate number of reads
        if read.is_reverse:
            numRelevantReadsReverse_common += 1
        else:
            numRelevantReadsForward_common += 1

        ## loop through all nucleotides in the read
        for i in range(len(sequence)):
            base = sequence[i]
            quality = qualities[i]

            if read.is_reverse:
                if base == "G":
                    mbias_reverse_no_qual_G_common[i] += 1

                    if quality >= minQual:
                        mbias_reverse_with_qual_G_common[i] += 1

                elif base == "A":
                    mbias_reverse_no_qual_A_common[i] += 1

                    if quality >= minQual:
                        mbias_reverse_with_qual_A_common[i] += 1

            else:
                if base == "C":
                    mbias_forward_no_qual_C_common[i] += 1

                    if quality >= minQual:
                        mbias_forward_with_qual_C_common[i] += 1

                elif base == "T":
                    mbias_forward_no_qual_T_common[i] += 1

                    if quality >= minQual:
                        mbias_forward_with_qual_T_common[i] += 1


    if read.query_length == readLenGenome:

        ## get sequence and qualities
        sequence = read.query_sequence
        qualities = read.query_qualities

        ## aggregate number of reads
        if read.is_reverse:
            numRelevantReadsReverse_genome += 1
        else:
            numRelevantReadsForward_genome += 1

        ## loop through all nucleotides in the read
        for i in range(len(sequence)):
            base = sequence[i]
            quality = qualities[i]

            if read.is_reverse:
                if base == "G":
                    mbias_reverse_no_qual_G_genome[i] += 1

                    if quality >= minQual:
                        mbias_reverse_with_qual_G_genome[i] += 1

                elif base == "A":
                    mbias_reverse_no_qual_A_genome[i] += 1

                    if quality >= minQual:
                        mbias_reverse_with_qual_A_genome[i] += 1

            else:
                if base == "C":
                    mbias_forward_no_qual_C_genome[i] += 1

                    if quality >= minQual:
                        mbias_forward_with_qual_C_genome[i] += 1

                elif base == "T":
                    mbias_forward_no_qual_T_genome[i] += 1

                    if quality >= minQual:
                        mbias_forward_with_qual_T_genome[i] += 1

bam.close()

## draw forward plot for most common read length
mbias_forward_no_qual_plot_common = (mbias_forward_no_qual_C_common / np.add(mbias_forward_no_qual_C_common, mbias_forward_no_qual_T_common) ) * 100
mbias_forward_with_qual_plot_common = (mbias_forward_with_qual_C_common / np.add(mbias_forward_with_qual_C_common, mbias_forward_with_qual_T_common) ) * 100
mbias_reverse_no_qual_plot_common = (mbias_reverse_no_qual_G_common / np.add( mbias_reverse_no_qual_G_common, mbias_reverse_no_qual_A_common) ) * 100
mbias_reverse_with_qual_plot_common = (mbias_reverse_with_qual_G_common / np.add( mbias_reverse_with_qual_G_common, mbias_reverse_with_qual_A_common) ) * 100

# Generate x-axis values based on the length of the arrays
x = np.arange(len(mbias_forward_no_qual_plot_common)) 

## plot the quality array first
plt.plot(x, mbias_forward_no_qual_plot_common, label='% of C', zorder=2)
plt.plot(x, mbias_forward_with_qual_plot_common, label=f'% of C with Q >= {minQual}', zorder=1)

max_value = np.max(mbias_forward_no_qual_plot_common)

# Set the y-axis range
y_range = [0, 50] if max_value <= 50 else [0, 100]

plt.ylim(y_range)

plt.xlabel('Position in Read')
plt.ylabel('Methylation Rate (%)')
plt.title('m-bias plot (forward)')
plt.legend()
plt.grid(True)
plt.savefig(f'{outFilePrefix}_mostCommonReadLength_forward.png', dpi=200)
plt.close()


## draw reverse plot
plt.plot(x, mbias_reverse_no_qual_plot_common, label='% of C', zorder=2)
plt.plot(x, mbias_reverse_with_qual_plot_common, label=f'% of C with Q >= {minQual}', zorder=1)

max_value = np.max(mbias_reverse_no_qual_plot_common)

# Set the y-axis range
y_range = [0, 50] if max_value <= 50 else [0, 100]

plt.ylim(y_range)

plt.xlabel('Position in Read')
plt.ylabel('% Methylation Rate (%)')
plt.title('m-bias plot (reverse)')
plt.legend()
plt.grid(True)
plt.savefig(f'{outFilePrefix}_mostCommonReadLength_reverse.png', dpi=200)
plt.close()


## draw forward plot for most common read length of genome alignment.bed
mbias_forward_no_qual_plot_genome = (mbias_forward_no_qual_C_genome / np.add(mbias_forward_no_qual_C_genome, mbias_forward_no_qual_T_genome) ) * 100
mbias_forward_with_qual_plot_genome = (mbias_forward_with_qual_C_genome / np.add(mbias_forward_with_qual_C_genome, mbias_forward_with_qual_T_genome) ) * 100
mbias_reverse_no_qual_plot_genome = (mbias_reverse_no_qual_G_genome / np.add( mbias_reverse_no_qual_G_genome, mbias_reverse_no_qual_A_genome) ) * 100
mbias_reverse_with_qual_plot_genome = (mbias_reverse_with_qual_G_genome / np.add( mbias_reverse_with_qual_G_genome, mbias_reverse_with_qual_A_genome) ) * 100

# Generate x-axis values based on the length of the arrays
x = np.arange(len(mbias_forward_no_qual_plot_genome)) 

## plot the quality array first
plt.plot(x, mbias_forward_no_qual_plot_genome, label='% of C', zorder=2)
plt.plot(x, mbias_forward_with_qual_plot_genome, label=f'% of C with Q >= {minQual}', zorder=1)

max_value = np.max(mbias_forward_no_qual_plot_genome)

# Set the y-axis range
y_range = [0, 50] if max_value <= 50 else [0, 100]

plt.ylim(y_range)

plt.xlabel('Position in Read')
plt.ylabel('Methylation Rate (%)')
plt.title('m-bias plot (forward)')
plt.legend()
plt.grid(True)
plt.savefig(f'{outFilePrefix}_mostCommonReadLengthOfGenomeAlignment_forward.png', dpi=200)
plt.close()


## draw reverse plot
plt.plot(x, mbias_reverse_no_qual_plot_genome, label='% of C', zorder=2)
plt.plot(x, mbias_reverse_with_qual_plot_genome, label=f'% of C with Q >= {minQual}', zorder=1)

max_value = np.max(mbias_reverse_no_qual_plot_genome)

# Set the y-axis range
y_range = [0, 50] if max_value <= 50 else [0, 100]

plt.ylim(y_range)

plt.xlabel('Position in Read')
plt.ylabel('% Methylation Rate (%)')
plt.title('m-bias plot (reverse)')
plt.legend()
plt.grid(True)
plt.savefig(f'{outFilePrefix}_mostCommonReadLengthOfGenomeAlignment_reverse.png', dpi=200)
plt.close()

## print coords for multiQC
xCoordsList = list(range(1, len(mbias_forward_no_qual_plot_genome)+1))

for i in range(len(xCoordsList)):
    print(f'{xCoordsList[i]}\t{mbias_forward_no_qual_plot_genome[i]}')