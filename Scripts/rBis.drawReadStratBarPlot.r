#!/usr/bin/env Rscript

#Create pdf/html report using R
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="Read_Stratification_", help="Output prefix. default=Read_Stratification_"),
	make_option(c("-r","--outPrefixAll"), default="Read_Stratification", help="Output prefix for total read stratification. default=Read_Stratification")
)

parser <- OptionParser(usage = "%prog [options] <tsv>", option_list=option_list,
			description = "Draw bar plot of read stratification.")
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	q()
} else {
	srcL=arguments$args
}

# Option handling
opt=arguments$options
outPrefix=opt$outPrefix
outPrefixAll=opt$outPrefixAll

## split into six groups
newL <- split(srcL, cut(seq_along(srcL), 7, labels = FALSE))

alignStats <- newL[[1]]
rRNAreadStats <- newL[[2]]
tRNAreadStats <- newL[[3]]
miRNAreadStats <- newL[[4]]
piRNAreadStats <- newL[[5]]
genomeReadStats <- newL[[6]]
circRNAreadStats <- newL[[7]]

## x-axis
samples = c()

## legends
sources = c()

## y-axis
counts = c()
tableCounts = c()

for (x in 1:length(alignStats)) {
	alignStatsTable <- read.table(file = alignStats[x], sep = '\t', header = TRUE)
	rRNAtable <- read.table(file = rRNAreadStats[x], sep = '\t', header = TRUE)
	tRNAtable <- read.table(file = tRNAreadStats[x], sep = '\t', header = TRUE)
	miRNAtable <- read.table(file = miRNAreadStats[x], sep = '\t', header = TRUE)
	piRNAtable <- read.table(file = piRNAreadStats[x], sep = '\t', header = TRUE)
	genomeTable <- read.table(file = genomeReadStats[x], sep = '\t', header = TRUE)
	circRNAtable <- read.table(file = circRNAreadStats[x], sep = '\t', header = TRUE)
	sampleName <- colnames(tRNAtable)[[1]]

	rRNA_multi <- rRNAtable$totalReads[1]
	tRNA_multi <- tRNAtable$totalReads[1]
	miRNA_multi <- miRNAtable$totalReads[1]
	piRNA_multi <- piRNAtable$totalReads[1]
	genome_multi <- genomeTable$totalReadsRegardlessOfStrand[1]
	circRNA_multi <- circRNAtable$totalReads[1]
	unaligned <- (circRNAtable$unaligned[1])/4

	samples = append(samples, c(sampleName, sampleName, sampleName, sampleName, sampleName, sampleName, sampleName))
	sources = append(sources, c("rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA", "Unaligned"))
	counts = append(counts, c(rRNA_multi, tRNA_multi, miRNA_multi, piRNA_multi, genome_multi, circRNA_multi, unaligned))
	tableCounts = append(tableCounts, c( sampleName, rRNA_multi, tRNA_multi, miRNA_multi, piRNA_multi, genome_multi, circRNA_multi, unaligned))

}

## save table with all info
sampleNames = unique(samples)
columns = unique(sources)
columns = append(columns, "Sample", 0)
tab = matrix( tableCounts, ncol=8, byrow=TRUE )
colnames(tab) = columns
tab = as.table(tab)
write.table(tab, file=paste0(outPrefixAll, "_All.tsv"), quote=FALSE, sep='\t', row.names=FALSE)

## save table without the unaligned info
sampleNames = unique(samples)
columns = unique(sources)
columns = append(columns, "Sample", 0)
tab = matrix( tableCounts, ncol=8, byrow=TRUE )
colnames(tab) = columns
newTab = subset(tab, select = -Unaligned)
newTab = as.table(newTab)
write.table(newTab, file=paste0(outPrefixAll, "_Aligned.tsv"), quote=FALSE, sep='\t', row.names=FALSE)

tmpList <- list(samples, sources, counts)
testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Sample','Source','Proportion')

testDF$Source <- factor(testDF$Source, levels=c("rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA", "Unaligned"))
testDF$Sample <- factor(testDF$Sample, levels=(unique(testDF$Sample)))

bars <- ggplot(testDF,
    aes(x = Sample, y = Proportion, fill = Source)) +
    geom_bar(position = "stack", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#A3A500", "#FF61CC", "#00B0F6", "#E76BF3", "#00CC33", "#FF9900")) +
	ggtitle("Read Stratification including Unaligned") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_Count.png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Count.pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

barsProportion <- ggplot(testDF,
	aes(x = Sample, y = Proportion, fill = Source)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#A3A500", "#FF61CC", "#00B0F6", "#E76BF3", "#00CC33", "#FF9900")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("Read Stratification including Unaligned (%)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_Percentage.png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Percentage.pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")
