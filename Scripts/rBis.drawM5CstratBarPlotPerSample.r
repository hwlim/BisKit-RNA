#!/usr/bin/env Rscript

#Create pdf/html report using R
suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))
suppressPackageStartupMessages(library('dplyr', quiet=TRUE))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="m5cStratification", help="Output prefix. default=m5cStratification"),
	make_option(c("-s","--sigThresh"), default=0.05, help="Output prefix. default=0.05"),
	make_option(c("-t","--sigType"), default="pVal", help="type of stat; either 'pVal' or 'FDR'. default=pVal"),
	make_option(c("-r","--mrThresh"), default=0.1, help="methylation rate threshold. default=0.1"),
	make_option(c("-c","--covThresh"), default=10, help="coverage threshold. default=10")
)
parser <- OptionParser(usage = "%prog [options]", option_list=option_list,
			description = "Draw read stratification bar plots using allReadStats.tsv")

arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	stop("Error: Requires mergedCandidates.tsv file(s) for all samples")
} else {
	src <- arguments$args
}

# Option handling
opt=arguments$options
outPrefix=opt$outPrefix
covThresh=as.numeric(opt$covThresh)
sigThresh=as.numeric(opt$sigThresh)
mrThresh=as.numeric(opt$mrThresh)

if (opt$sigType == "pVal") {
	sigSave = "pVal"
	sigType = "p.value_mState"
} else {	
	sigSave = "FDR"
	sigType = "FDR_mState"
}


## x-axis
samples = c()
samplesAnnotation = c()

## legends
sources = c()
annotations = c()

## y-axis
countsAll = c()
countsSig = c()
countsGenome = c()


## og df
statsTable <- read.table(file = src, sep = '\t', header = TRUE, comment.char = "")

## get sample name
samplename = tail(strsplit(src, "/")[[1]], n=3)[1]

## Get total number of C locations per source
rRNAall = length(which(statsTable$Source=="rRNA"))
tRNAall = length(which(statsTable$Source=="tRNA"))
miRNAall = length(which(statsTable$Source=="miRNA"))
genomeall = length(which(statsTable$Source=="Genome"))
    
samples = append(samples, c(samplename, samplename, samplename, samplename))
sources = append(sources, c("rRNA", "tRNA", "miRNA", "Genome"))
countsAll = append(countsAll, c(rRNAall, tRNAall, miRNAall, genomeall))

## filtered df
sigTable = statsTable %>% filter(.[[sigType]]<sigThresh & cov>=covThresh & methRate>=mrThresh)

## Get significant number of C locations per source
rRNAsig = length(which(sigTable$Source=="rRNA"))
tRNAsig = length(which(sigTable$Source=="tRNA"))
miRNAsig = length(which(sigTable$Source=="miRNA"))
genomesig = length(which(sigTable$Source=="Genome"))
    
countsSig = append(countsSig, c(rRNAsig, tRNAsig, miRNAsig, genomesig))

## genome annotated
annotated = length(which(sigTable$Source=="Genome" & sigTable$gene_type!="no_gene_type"))
notAnnotated = length(which(sigTable$Source=="Genome" & sigTable$gene_type=="no_gene_type"))

samplesAnnotation = append(samplesAnnotation, c(samplename, samplename))
annotations = append(annotations, c("Annotated", "Not Annotated"))
countsGenome = append(countsGenome, c(annotated, notAnnotated))


## all
tmpList <- list(samples, sources, countsAll)
testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Sample','Source','Counts')

testDF$Source <- factor(testDF$Source, levels=c("rRNA", "tRNA", "miRNA", "Genome"))
testDF$Sample <- factor(testDF$Sample, levels=(unique(testDF$Sample)))

barsProportion <- ggplot(testDF,
               aes(x = Sample,
                   y = Counts,
                   fill = Source)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#A3A500", "#FF61CC", "#00B0F6", "#E76BF3")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("% of C locations with at least 1 coverage") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Samples") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_(All).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_(All).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")


## significant only
tmpList <- list(samples, sources, countsSig)
testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Sample','Source','Counts')

testDF$Source <- factor(testDF$Source, levels=c("rRNA", "tRNA", "miRNA", "Genome"))
testDF$Sample <- factor(testDF$Sample, levels=(unique(testDF$Sample)))

barsProportion <- ggplot(testDF,
               aes(x = Sample,
                   y = Counts,
                   fill = Source)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#A3A500", "#FF61CC", "#00B0F6", "#E76BF3")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("% of significant m5C candidates") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Samples") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_(", sigSave,").png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_(", sigSave,").pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")


## draw bar plot for genome annotated and not annotated
tmpList <- list(samplesAnnotation, annotations, countsGenome)
testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Sample','Annotations','Counts')

testDF$Source <- factor(testDF$Annotations, levels=c("Annotated", "Not Annotated"))
testDF$Sample <- factor(testDF$Sample, levels=(unique(testDF$Sample)))

barsProportion <- ggplot(testDF,
               aes(x = Sample,
                   y = Counts,
                   fill = Annotations)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#00B0F6")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("m5C Stratification (% Based on annotation)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Samples") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_(", sigSave, ")_(Annotated).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_(", sigSave, ")_(Annotated).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")
