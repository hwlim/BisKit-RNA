#!/usr/bin/env bash
## Script to perform snakemake dry-run

if [ -z ${BISKIT_PATH+x} ]; then
	echo -e "Error: Environment variable BISKIT is not defined. Refer to the initial setup section in the github page at https://github.com/hwlim/bisKit for instructions setting up this enrironment variable." >&2
	exit 1
fi

# echo -e "Purging pre-loaded modules to prevent conflicts..." >&2
module purge
# echo -e "Loading anaconda module..." >&2
module load anaconda3
# echo -e "Loading anaconda virtual environment..." >&2
source activate BisKit

echo -e "Performing dry-run..." >&2
snakemake -np -s ${BISKIT_PATH}/Snakemake/Snakefile


# echo -e "Deactivating anaconda virtual environment..." >&2
conda deactivate
# echo -e "Purging anaconda module..." >&2
module unload anaconda3
