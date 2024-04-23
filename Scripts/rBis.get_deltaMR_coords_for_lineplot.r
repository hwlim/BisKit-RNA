#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))
suppressPackageStartupMessages(library('tidyr', quiet=TRUE))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="delta_MR_all_pairwise", help="Output prefix. default=delta_MR_all_pairwise")
)
parser <- OptionParser(usage = "%prog [options] <tsv>", option_list=option_list,
			description = "Get coords for delta methRate distribution. This will be used to draw a lineplot in MultiQC.")
arguments <- parse_args(parser, positional_arguments = TRUE)

opt=arguments$options
outPrefix=opt$outPrefix

if(length(arguments$args) == 0) {
	print_help(parser)
	q()
} else {
	srcL=arguments$args
}

compNames = c()
distList = c()
n = length(srcL)

for( src in srcL ){
	statsTable <- read.table(file = src, sep = '\t', header = TRUE, comment.char="")

    y_cols = grep("delta_MethRate_", names(statsTable))
    comparison = strsplit(src, "/")[[1]][2]
    deltaMRpos <- subset(statsTable[, y_cols], statsTable[, y_cols] > 0.05 )
    deltaMRneg <- subset(statsTable[, y_cols], statsTable[, y_cols] < -0.05 )
    delta = c(deltaMRpos, deltaMRneg)

    compNames = append(compNames, comparison)
    distList = append(distList, list(delta))
}

# extract X and Y coordinates for each list
xcoords <- list()
ycoords <- list()

for (i in 1:n) {
    density_i <- density(distList[[i]])
    xcoords[[i]] <- density_i$x
    ycoords[[i]] <- density_i$y
}

# return X and Y coordinates
final = list(xcoords = xcoords, ycoords = ycoords)

for (i in 1:length(final$xcoords)) {
    comp = compNames[i]
    xcoord = final$xcoords[[i]]
    ycoord = final$ycoords[[i]]
    table = data.frame( col1=xcoord, col2=ycoord )
    write.table(table, paste0(outPrefix, "/delta_MR_", comp, ".tsv"), sep = "\t", row.names=FALSE, col.names=FALSE)
}
