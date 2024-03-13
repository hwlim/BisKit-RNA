#!/usr/bin/env bash

## Written by Hee Woong Lim
##
## Script to submit snakemake job to CCHMC lsf system

source $COMMON_LIB_BASE/commonBash.sh


nJob=50
totalWaitTime="72:00"
#timestamp=$(date +%Y%m%d_%H%M%S)
config=${BISKIT_PATH}/Snakemake/cluster.yml

assertFileExist $config
assertFileExist ${BISKIT_PATH}/Snakemake/Snakefile



if [ ! -f diag.pdf ];then
	module load python3/3.6.3
	module load graphviz/2.40.1
	snakemake -s ${BISKIT_PATH}/Snakemake/Snakefile --dag | dot -Tpdf > diag.pdf
	module purge
fi

usrName=$(whoami)
mkdir -p logs
bsub -M 4000 -W ${totalWaitTime} -eo submit.err -oo submit.out \
	"
	module load python3/3.6.3

	snakemake -s ${BISKIT_PATH}/Snakemake/Snakefile --rerun-incomplete \
		-j $nJob \
		--latency-wait 60 \
		--cluster-config $config \
		--cluster 'bsub -W {cluster.walltime} -n {cluster.cpu} -M {cluster.memory} -J $$.{cluster.name} -R {cluster.resource} -eo {cluster.error} -oo {cluster.output}'
	
	"

