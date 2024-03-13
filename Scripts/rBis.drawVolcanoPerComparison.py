#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#parse command line arguments
options = argparse.ArgumentParser(description="Draw volcano plots for each comparison", usage="python3 drawVolcanoPerComparison.py [options] -c comparison.tsv")
options.add_argument('-c', '--comparison', required=True,
                        help='categorized_pairwise_comparison.tsv file')
options.add_argument('-o', '--outPrefix', default="call.tsv",
                        help='output file destination.')
options.add_argument('-t', '--sigType', default="pVal",
                        help='Significance type; default = pVal; use pVal or FDR')
options.add_argument('-y', '--yThresh', default=30,
                        help='output plot y-axis limit')

args = options.parse_args()
comparisonFile = args.comparison
outPrefix = args.outPrefix
sigType = args.sigType
yThresh = int(args.yThresh)

if args.sigType == "pVal":
    sigType = "p-value_mState"
    sigSave = "pVal"
elif args.sigType == "FDR":
    sigType = "FDR_mState"
    sigSave = "FDR"
else:
    print("need to use pVal or FDR as -t flag.")
    exit()

comparisonFile = pd.read_csv(comparisonFile, sep="\t", header=0)

compName = comparisonFile.columns[0]
samp1 = compName.split('_vs_')[0]
samp2 = compName.split('_vs_')[1]

data = comparisonFile[comparisonFile[compName].isin(['DOWN', 'UP', 'UNCHANGED'])]

data["coords"] = data.apply(lambda row: f"{row['#SeqID']}:{row['refPos']}:{row['refStrand']}", axis=1)

deltaMR = f"delta_MethRate_{compName}"
sig = f"{sigType}_{compName}"


# Convert the 'value' column to numeric type
data[sig] = pd.to_numeric(data[sig])

## get -log10 of significance threshold
data['logSig'] = np.log10(data[sig]) * -1

## cap all values over yThresh to yThresh
data['logSig'] = data['logSig'].clip(upper=yThresh)

columns_to_subset = [compName, "logSig", deltaMR]
df_subset = data[columns_to_subset]

# Filter the DataFrame based on the values in the "RG1_vs_RG2" column
up_data = df_subset[df_subset[compName] == 'UP']
down_data = df_subset[df_subset[compName] == 'DOWN']
unchanged_data = df_subset[df_subset[compName] == 'UNCHANGED']

# Create the volcano plot
plt.scatter(up_data[deltaMR], up_data['logSig'], color='red', label='UP', s=1, zorder=3)
plt.scatter(down_data[deltaMR], down_data['logSig'], color='blue', label='DOWN', s=1, zorder=3)
plt.scatter(unchanged_data[deltaMR], unchanged_data['logSig'], color='gray', label='UNCHANGED', s=1, zorder=2)

# axis limits
plt.xlim(-1, 1)
plt.ylim(0, yThresh+1)

# Add labels and title
plt.xlabel(f'Delta Methylation Rate ({samp2} − {samp1})')
plt.ylabel(f'log10({sigSave})')
plt.title(f'm5C Categorization')
plt.legend()
plt.grid(True, zorder=1, alpha=0.3)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

# Show the plot
plt.savefig(f'{outPrefix}.png', dpi=150)
plt.savefig(f'{outPrefix}.pdf', dpi=150)
