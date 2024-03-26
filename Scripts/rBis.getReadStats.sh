#!/usr/bin/env bash

set -e


function printUsage {
	echo -e "Usage: `basename $0` -o <destination directory> <original bam file> <plus bam file> <minus bam file> <sample name>
Description: read stat calculation from resulting bam files after post processing
Options:
	-o Output file name including path, default=<src file name>
		<outputFile>
Input:
	align.bam
	align.plus.bam
	align.minus.bam
    sample name
Output:
	readStats.tsv" >&2
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
plusBam=$2
minusBam=$3
sampleName=$4
unalignedFQ=$5

###################################
## main code

totalReads=$(samtools view $bam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | wc -l)
totalUniq=$(samtools view $bam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | cut -f1 | grep -w "1" | wc -l)
totalMulti=$(samtools view $bam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | cut -f1 | grep -vw "1" | wc -l)
totalPlusReads=$(samtools view $plusBam | cut -f1 | sort -S 1G | uniq | wc -l)
plusUniq=$(samtools view $plusBam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | cut -f1 | grep -w "1" | wc -l)
plusMulti=$(samtools view $plusBam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | cut -f1 | grep -vw "1" | wc -l)
totalMinusReads=$(samtools view $minusBam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | wc -l)
minusUniq=$(samtools view $minusBam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | cut -f1 | grep -w "1" | wc -l)
minusMulti=$(samtools view $minusBam | cut -f1 | ${BISKIT_PATH}/Scripts/sortByFreq.sh | cut -f1 | grep -vw "1" | wc -l)
unaligned=$( zcat $unalignedFQ | wc -l)

echo totalReads: $totalReads
echo totalUniq: $totalMulti
echo totalMulti: $totalMulti
echo totalPlusReads: $totalPlusReads
echo plusUniq: $plusUniq
echo plusMulti: $plusMulti
echo totalMinusReads: $totalMinusReads
echo minusUniq: $minusUniq
echo minusMulti: $minusMulti
echo unaligned: $unaligned

## write to output file
echo | awk -vs="$sampleName" -va="totalReadsRegardlessOfStrand" -ve="totalUniqAlignedRegardlessOfStrand" -vf="totalMultiAlignedRegardlessOfStrand" -vb="totalReads" -vc="strandUniq" -vd="strandMulti" -vg="unaligned" '{{ print s"\t"e"\t"f"\t"a"\t"b"\t"c"\t"d"\t"g }}' > ${des}
echo | awk -vs="plus" -va="$totalReads" -ve="$totalUniq" -vf="$totalMulti" -vb="$totalPlusReads" -vc="$plusUniq" -vd="$plusMulti" -vg="$unaligned" '{{ print s"\t"e"\t"f"\t"a"\t"b"\t"c"\t"d"\t"g }}' >> ${des}
echo | awk -vs="minus" -va="$totalReads" -ve="$totalUniq" -vf="$totalMulti" -vb="$totalMinusReads" -vc="$minusUniq" -vd="$minusMulti" -vg="$unaligned" '{{ print s"\t"e"\t"f"\t"a"\t"b"\t"c"\t"d"\t"g }}' >> ${des}
