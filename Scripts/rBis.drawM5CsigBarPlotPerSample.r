#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="m5cSignificant", help="Output prefix. default=allReadsStratification"),
	make_option(c("-t","--sigType"), help="Significant type; pVal or FDR")
)
parser <- OptionParser(usage = "%prog [options] <tsv>", option_list=option_list,
			description = "Draw bar plot of proportion of significant m5C candidates per Source")
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	q()
} else {
	src=arguments$args
}

# Option handling
opt=arguments$options
outPrefix=opt$outPrefix
sigType=opt$sigType

statsTable <- read.table(file = src, sep = '\t', header = TRUE)

rRNAsig = statsTable[[paste0("num_Candidates_cov_MR_", sigType, "_filtered")]][1]
tRNAsig = statsTable[[paste0("num_Candidates_cov_MR_", sigType, "_filtered")]][2]
miRNAsig = statsTable[[paste0("num_Candidates_cov_MR_", sigType, "_filtered")]][3]
genomeSig = statsTable[[paste0("num_Candidates_cov_MR_", sigType, "_filtered")]][4]

rRNAnotSig = statsTable$num_Candidates_cov..10[1] - rRNAsig
tRNAnotSig = statsTable$num_Candidates_cov..10[2] - tRNAsig
miRNAnotSig = statsTable$num_Candidates_cov..10[3] - miRNAsig
genomeNotSig = statsTable$num_Candidates_cov..10[4] - genomeSig

tmpList <- list( 
			source <- c("rRNA", "rRNA", "tRNA", "tRNA", "miRNA", "miRNA", "Genome", "Genome"),
			sig <- c("Significant", "Not Significant", "Significant", "Not Significant", "Significant", "Not Significant", "Significant", "Not Significant"),
			candidates <- c(rRNAsig, rRNAnotSig, tRNAsig, tRNAnotSig, miRNAsig, miRNAnotSig, genomeSig, genomeNotSig)
			)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Source','Significance','cand')

testDF$Significance <- factor(testDF$Significance, levels=c('Significant', 'Not Significant'))
testDF$Source <- factor(testDF$Source, levels=(unique(testDF$Source)))

bars <- ggplot(testDF,
               aes(x = Source,
                   y = cand,
                   fill = Significance)) +
    geom_bar(position = "stack", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#00B0F6")) +
	ggtitle("Significant m5C Candidates") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Source") +
	ylab("Counts")

ggsave(paste0(outPrefix, "_Count.png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Count.pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

barsProportion <- ggplot(testDF,
               aes(x = Source,
                   y = cand,
                   fill = Significance)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("#F8766D", "#00B0F6")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("Significant m5C Candidates (%)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Source") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_Percentage.png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Percentage.pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")
