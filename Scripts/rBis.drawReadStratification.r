#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="m5cSignificant", help="Output prefix. default=allReadsStratification"),
    make_option(c("-a","--alignStats"), help="merged alignment stats for all sources"),
    make_option(c("-r","--rRNAreadStats"), help="read stats for filtered rRNA"),	
    make_option(c("-t","--tRNAreadStats"), help="read stats for filtered tRNA"),
    make_option(c("-m","--miRNAreadStats"), help="read stats for filtered miRNA"),
    make_option(c("-p","--piRNAreadStats"), help="read stats for filtered piRNA"),
    make_option(c("-c","--circRNAreadStats"), help="read stats for filtered circRNA")
)

parser <- OptionParser(usage = "%prog [options] <tsv>", option_list=option_list,
			description = "Draw bar plot of read stratification for all sources")
arguments <- parse_args(parser, positional_arguments = TRUE)

# Option handling
opt=arguments$options
outPrefix=opt$outPrefix
alignStats=opt$alignStats
rRNAreadStats=opt$rRNAreadStats
tRNAreadStats=opt$tRNAreadStats
miRNAreadStats=opt$miRNAreadStats
piRNAreadStats=opt$piRNAreadStats
circRNAreadStats=opt$circRNAreadStats

## x-axis
samp = c()

## legends
source = c()

## y-axis
counts = c()

alignStatsTable <- read.table(file = alignStats, sep = '\t', header = TRUE)
rRNAtable <- read.table(file = rRNAreadStats, sep = '\t', header = TRUE)
tRNAtable <- read.table(file = tRNAreadStats, sep = '\t', header = TRUE)
miRNAtable <- read.table(file = miRNAreadStats, sep = '\t', header = TRUE)
piRNAtable <- read.table(file = piRNAreadStats, sep = '\t', header = TRUE)
circRNAtable <- read.table(file = circRNAreadStats, sep = '\t', header = TRUE)
sampleName <- colnames(tRNAtable)[[1]]

## Multi align stats
rRNA_multi <- rRNAtable$totalReads[1] + alignStatsTable$Multi.Aligned[1]
tRNA_multi <- tRNAtable$totalReads[1] + alignStatsTable$Multi.Aligned[2]
miRNA_multi <- miRNAtable$totalReads[1] + alignStatsTable$Multi.Aligned[3]
piRNA_multi <- piRNAtable$totalReads[1] + alignStatsTable$Multi.Aligned[4]
genome <- alignStatsTable$Uniquely.Aligned[5]
circRNA_multi <- circRNAtable$totalReads[1] + alignStatsTable$Multi.Aligned[6]
unaligned <- alignStatsTable$Unaligned[6] + alignStatsTable$Multi.Aligned[6]

samples = append(samp, c(sampleName, sampleName, sampleName, sampleName, sampleName, sampleName, sampleName))
sources = append(source, c("rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA", "Unaligned"))
counts = append(counts, c(rRNA_multi, tRNA_multi, miRNA_multi, piRNA_multi, genome, circRNA_multi, unaligned))

tmpList <- list(samples, sources, counts)
testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Sample','Source','Proportion')

testDF$Source <- factor(testDF$Source, levels=c("rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA", "Unaligned"))
testDF$Sample <- factor(testDF$Sample, levels=(unique(testDF$Sample)))

barsProportionUnaligned <- ggplot(testDF,
               aes(x = Sample,
                   y = Proportion,
                   fill = Source)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#A3A500", "#FF61CC", "#00B0F6", "#E76BF3", "#00CC33", "#FF9900")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("Read Stratification including Unaligned (%)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_With_Unaligned.png"), barsProportionUnaligned, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_With_Unaligned.pdf"), barsProportionUnaligned, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

### without unaligned
## x-axis
samp = c()
## legends
source = c()
## y-axis
counts = c()

samples = append(samp, c(sampleName, sampleName, sampleName, sampleName, sampleName, sampleName))
sources = append(source, c("rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA"))
counts = append(counts, c(rRNA_multi, tRNA_multi, miRNA_multi, piRNA_multi, genome, circRNA_multi ))

tmpList <- list(samples, sources, counts)
testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Sample','Source','Proportion')

testDF$Source <- factor(testDF$Source, levels=c("rRNA", "tRNA", "miRNA", "piRNA", "Genome", "circRNA"))
testDF$Sample <- factor(testDF$Sample, levels=(unique(testDF$Sample)))

barsProportion <- ggplot(testDF,
               aes(x = Sample,
                   y = Proportion,
                   fill = Source)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#A3A500", "#FF61CC", "#00B0F6", "#E76BF3", "#00CC33")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("Read Stratification without Unaligned (%)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_Without_Unaligned.png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Without_Unaligned.pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")
