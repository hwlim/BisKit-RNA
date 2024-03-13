#!/usr/bin/env bash

# if [ -z ${COMMON_LIB_BASE+x} ]; then
# 	echo -e "Error: LimLabBase is not properly set up" >&2
# 	exit 1
# fi
# source $COMMON_LIB_BASE/commonBash.sh

if [ -z ${BISKIT_PATH+x} ]; then
	echo -e "Error: Environment variable BISKIT is not defined. Refer to the initial setup section in the github page at https://github.com/hwlim/bisKit for instructions setting up this enrironment variable." >&2
	exit 1
fi



echo -e "Initializing BisKit" >&2
echo -e "BisKit path = ${BISKIT_PATH}" >&2


cp -i -v ${BISKIT_PATH}/Config/sample.tsv .
cp -i -v ${BISKIT_PATH}/Config/config.bis.yml .
