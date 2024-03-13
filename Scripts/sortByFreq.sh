#!/usr/bin/env bash

source $COMMON_LIB_BASE/commonBash.sh
trap 'if [ `ls -1 ${TMPDIR}/__temp__.$$.* 2>/dev/null | wc -l` -gt 0 ];then rm ${TMPDIR}/__temp__.$$.*; fi' EXIT

function printUsage {
	echo "Usage: `basename $0` (options) [input file]
Description:
	For a give list of elements, count frequency and print sorted list by frequency. For example
	Input:
		a
		c
		a
		b
	Output:
		1	b
		1	c
		2	a

Options:
	-r : if set, sort by decreasing order
Input
	- Plaint text file or no input specification for stdin
	  No need to be sorted" >&2

}



###################################
## option and input file handling
decrease=0
while getopts ":rh" opt; do
	case $opt in
		r)
			decrease=1
			;;
		h)
			printUsage
			exit 1
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
if [ $# -eq 0 ];then
	src=/dev/stdin
else
	src=$1
	assertFileExist $1
fi




###################################
## main code

opt="-k1,1n"
if [ $decrease -eq 1 ];then
	opt="-k1,1nr"
fi

cat $src \
	| sort \
	| uniq -c \
	| sed -e 's/^ *//' -e 's/ /\t/' \
	| sort $opt 

