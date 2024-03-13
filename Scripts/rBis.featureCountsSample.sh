#!/usr/bin/env bash

trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [bam]
Description:
	A wrapper script to run featureCounts for single bam file and produce compact ouptuts
Output:
	- <outPrefix>.txt
	- <outPrefix>.txt.summary
	- <outPrefix>.log
Options:
	-o <outPrefix>: output file prefix. default=featureCount
	-g <gtf/gff/bed>: gene annotation file
	-s <optStr>: other option string for featureCounts. default=\"-O -T 1 -t exon -g gene_id\"
		For example:
		if paired-end: add '-p' option
		for library strand, add:
			-s 0: for unstranded
			-s 1: for same stranded (fr-secondstrand)
			-s 2: for opposite stranded (fr-firststrand)
		for min overlap:
			--fracOverlap [0-1]. default=0"
}

if [ $# -eq 0 ];then
	printUsage
	exit 1
fi


###################################
## option and input file handling
gtf=""
outPrefix=""
optStr="-O -T 1 -t exon -g gene_id"
while getopts ":g:o:s:" opt; do
	case $opt in
		g)
			gtf=$OPTARG
			;;
		o)
			outPrefix=$OPTARG
			;;
		s)
			optStr=${OPTARG}
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

if [ "$gtf" = "" ];then
	echo "Error: gene annotation file must be specified (-g)" >&2
	printUsage
	exit 1
fi

shift $((OPTIND-1))
if [ $# -eq 0 ];then
	printUsage
	exit 1
fi

assertFileExist()
{
	local FILEL=($@)
	local f
	for f in ${FILEL[@]}
	do
		if [ ! -f $f ]; then
			echo -e "Error: File $f does not exist." >&2
			exit 1
		fi
	done
}


bam=$1
assertFileExist $gtf $bam

tmpRaw=${TMPDIR}/__temp__.$$.raw
tmpLog=${TMPDIR}/__temp__.$$.log
tmpSum=${TMPDIR}/__temp__.$$.raw.summary
tmpTxt=${TMPDIR}/__temp__.$$.txt


desDir=`dirname ${outPrefix}`
mkdir -p $desDir

desTxt=${outPrefix}.txt
desSum=${outPrefix}.txt.summary
desLog=${outPrefix}.log

featureCounts -a ${gtf} ${optStr} -o ${tmpRaw} ${bam} 2>&1 | tee ${tmpLog}

set -o pipefail
head -n 1 $tmpRaw > $tmpTxt
echo -e "Gene_id\tLength\tCount" >> $tmpTxt
tail -n +3 $tmpRaw | cut -f 1,6,7 >> $tmpTxt
set +o pipefail

mv -f ${tmpTxt} ${desTxt}
mv -f ${tmpSum} ${desSum}
mv -f ${tmpLog} ${desLog}
