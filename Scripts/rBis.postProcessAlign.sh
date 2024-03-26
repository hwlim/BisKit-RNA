#!/usr/bin/env bash

set -e


function printUsage {
	echo -e "Usage: `basename $0` -o <destination directory> <bam file> <original.fq.gz> <unaligned_pre.fq.gz>
Description: Post-processing of output bam file from hisat-3n
Options:
	-o Output file prefix including path, default=<src file name>
		<outPrefix>.align.plus.bam
		<outPrefix>.align.plus.uniq.bam
		<outPrefix>.align.minus.bam
Input:
	align.bam
	original.fq
	unmapped_pre.fq
Output:
	align.plus.bam
    align.plus.uniq.bam
	align.minus.bam
    unmapped.fq" >&2
}

###################################
## option and input file handling
des=""
while getopts ":o:" opt; do
	case $opt in
		o)
			des=$OPTARG
			;;
		\?)
			echo "Invalid options: -$OPTARG" >&2
			printUsage
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			printUsage
			exit 1
			;;
	esac
done


shift $((OPTIND-1))
if [ $# -eq 2 ];then
	printUsage
	exit 1
fi
bam=$1
originalFQ=$2
unmappedPreFQ=$3

###################################
## main code

## split to forward strand
echo "making align.plus.bam."
samtools view -b -F 16 $bam > ${des}/align.plus.bam
samtools index ${des}/align.plus.bam

## split to minus strand
echo "making align.minus.bam."
samtools view -b -f 16 $bam > ${des}/align.minus.bam
samtools index ${des}/align.minus.bam

## get align.plus.uniq.bam
echo "making align.plus.uniq.bam."
samtools view ${des}/align.plus.bam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh -r > ${des}/uniqReads.txt
gawk -F'\t' '$1 ==1 {print $2}' ${des}/uniqReads.txt > ${des}/readNamesOfInterest.txt
samtools view -N ${des}/readNamesOfInterest.txt -o ${des}/align.plus.uniq.bam ${des}/align.plus.bam
samtools index ${des}/align.plus.uniq.bam

## get discarded reads
echo "getting discarded reads."
### get read names in forward and reverse bam
samtools view ${des}/align.plus.bam | cut -f1 | sort -S 1G | uniq > ${des}/plusReads.txt
samtools view ${des}/align.minus.bam | cut -f1 | sort -S 1G | uniq > ${des}/minusReads.txt

### get reads that are only found in align.minus.bam
comm -13 ${des}/plusReads.txt ${des}/minusReads.txt > ${des}/discardedReadNames.txt

### extract discarded read names from 
seqtk subseq $originalFQ ${des}/discardedReadNames.txt > ${des}/tmp_discarded.fq
gzip -f ${des}/tmp_discarded.fq

cat ${des}/tmp_discarded.fq.gz $unmappedPreFQ > ${des}/unaligned.fq.gz

\rm ${des}/tmp_discarded.fq.gz
\rm ${des}/discardedReadNames.txt
\rm ${des}/plusReads.txt
\rm ${des}/minusReads.txt
\rm ${des}/uniqReads.txt
\rm ${des}/readNamesOfInterest.txt
