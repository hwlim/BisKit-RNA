#!/usr/bin/env python3

import pysam
import argparse

#parse command line arguments
options = argparse.ArgumentParser(description="Organizes the alignment statistics output by hisat-3n to machine-friendly format.", usage="python alignStat_perSample.py [options] -n sampleName -s summary.out ")
options.add_argument('-b', '--bam', required=True,
                        help='align.bam file from hisat-3n')

args = options.parse_args()

def count_chromosomes(bam_file_path):
    # Dictionary to store chromosome counts
    chromosome_counts = {}

    # Open the BAM file
    with pysam.AlignmentFile(bam_file_path, 'rb') as bamfile:
        for read in bamfile:
            chromosome = read.reference_name
            chrom = chromosome.split("|")[0].split("tRNA-")[-1]
            if chrom in chromosome_counts:
                chromosome_counts[chrom] += 1
            else:
                chromosome_counts[chrom] = 1

    return chromosome_counts

if __name__ == "__main__":
    # Path to your BAM file
    bam_file_path = args.bam

    # Count chromosomes
    chromosome_counts = count_chromosomes(bam_file_path)

    print("Gene_id\tLength\tCount")
    
    # Print the chromosome counts
    for chrom, count in chromosome_counts.items():
        print(f'{chrom}\t0\t{count}')
        