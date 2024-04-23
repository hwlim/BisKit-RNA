#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="m5cSignificant", help="Output prefix. default=m5cSignificant"),
	make_option(c("-g","--groupList"), help="a string of groups separated by comma")
)
parser <- OptionParser(usage = "%prog [options] <tsv>", option_list=option_list,
			description = "Draw bar plot of proportion of significant m5C candidates per Source and output table.
			Input:
				- categorized_pairwise_comparison.tsv
				- list of groups involved
			Output:
				- static plots for read categorization
				- table containing read categorization
			")
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
samp1=strsplit( opt$groupList, "," )[[1]][1]
samp2=strsplit( opt$groupList, "," )[[1]][2]

## x-axis
sources = c()

## legend
trend = c()

## y-axis
counts = c()

statsTable <- read.table(file = src, sep = '\t', header = TRUE, comment.char = "")
# compName = tail(strsplit(src, "/")[[1]], n=3)[1]
compName = paste0(samp1, "_vs_", samp2)

## select row based on column values
uprRNA = nrow(subset(statsTable, statsTable[[compName]] == "UP" & Source == "rRNA"))
downrRNA = nrow(subset(statsTable, statsTable[[compName]] == "DOWN" & Source == "rRNA"))
unchangedrRNA = nrow(subset(statsTable, statsTable[[compName]] == "UNCHANGED" & Source == "rRNA"))
uniq1rRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp1) & Source == "rRNA"))
uniq2rRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp2) & Source == "rRNA"))

uptRNA = nrow(subset(statsTable, statsTable[[compName]] == "UP" & Source == "tRNA"))
downtRNA = nrow(subset(statsTable, statsTable[[compName]] == "DOWN" & Source == "tRNA"))
unchangedtRNA = nrow(subset(statsTable, statsTable[[compName]] == "UNCHANGED" & Source == "tRNA"))
uniq1tRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp1) & Source == "tRNA"))
uniq2tRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp2) & Source == "tRNA"))

upmiRNA = nrow(subset(statsTable, statsTable[[compName]] == "UP" & Source == "miRNA"))
downmiRNA = nrow(subset(statsTable, statsTable[[compName]] == "DOWN" & Source == "miRNA"))
unchangedmiRNA = nrow(subset(statsTable, statsTable[[compName]] == "UNCHANGED" & Source == "miRNA"))
uniq1miRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp1) & Source == "miRNA"))
uniq2miRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp2) & Source == "miRNA"))

uppiRNA = nrow(subset(statsTable, statsTable[[compName]] == "UP" & Source == "piRNA"))
downpiRNA = nrow(subset(statsTable, statsTable[[compName]] == "DOWN" & Source == "piRNA"))
unchangedpiRNA = nrow(subset(statsTable, statsTable[[compName]] == "UNCHANGED" & Source == "piRNA"))
uniq1piRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp1) & Source == "piRNA"))
uniq2piRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp2) & Source == "piRNA"))

upGenome = nrow(subset(statsTable, statsTable[[compName]] == "UP" & Source == "Genome"))
downGenome = nrow(subset(statsTable, statsTable[[compName]] == "DOWN" & Source == "Genome"))
unchangedGenome = nrow(subset(statsTable, statsTable[[compName]] == "UNCHANGED" & Source == "Genome"))
uniq1Genome = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp1) & Source == "Genome"))
uniq2Genome = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp2) & Source == "Genome"))

upcircRNA = nrow(subset(statsTable, statsTable[[compName]] == "UP" & Source == "circRNA"))
downcircRNA = nrow(subset(statsTable, statsTable[[compName]] == "DOWN" & Source == "circRNA"))
unchangedcircRNA = nrow(subset(statsTable, statsTable[[compName]] == "UNCHANGED" & Source == "circRNA"))
uniq1circRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp1) & Source == "circRNA"))
uniq2circRNA = nrow(subset(statsTable, statsTable[[compName]] == paste0('Unique to ', samp2) & Source == "circRNA"))

sumUp = uprRNA + uptRNA + upmiRNA + uppiRNA + upGenome + upcircRNA
sumDown = downrRNA + downtRNA + downmiRNA + downpiRNA + downGenome + downcircRNA
sumUniq1 = uniq1rRNA + uniq1tRNA + uniq1miRNA + uniq1piRNA + uniq1Genome + uniq1circRNA
sumUniq2 = uniq2rRNA + uniq2tRNA + uniq2miRNA + uniq2piRNA + uniq2Genome + uniq2circRNA

## UP and DOWN
sources = append(sources, c("rRNA", "rRNA", "tRNA", "tRNA", "miRNA", "miRNA", "piRNA", "piRNA", "Genome", "Genome", "circRNA", "circRNA"))
trend = append(trend, c("Up", "Down", "Up", "Down", "Up", "Down", "Up", "Down", "Up", "Down", "Up", "Down"))
counts = append(counts, c(uprRNA, downrRNA, uptRNA, downtRNA, upmiRNA, downmiRNA, uppiRNA, downpiRNA, upGenome, downGenome, upcircRNA, downcircRNA))

tmpList <- list(sources, trend, counts)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Source','Trend','Counts')

testDF$Trend <- factor(testDF$Trend, levels=c('Up', 'Down'))
testDF$Source <- factor(testDF$Source, levels=(unique(testDF$Source)))


# bars <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "stack", stat = "identity") +
# 	scale_fill_manual(values=c("Red","Blue")) +
# 	ggtitle("m5C Categorization Between Comparisons") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Source") +
# 	ylab("Counts")

# ggsave(paste0(outPrefix, "_Count_(UPvDOWN).png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Count_(UPvDOWN).pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

# barsProportion <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "fill", stat = "identity") +
# 	scale_fill_manual(values=c("Red","Blue")) +
# 	scale_y_continuous(labels = scales::percent_format()) +
# 	ggtitle("m5C Categorization Between Comparisons (%)") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Source") +
# 	ylab("Proportion")

# ggsave(paste0(outPrefix, "_Percentage_(UPvDOWN).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Percentage_(UPvDOWN).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")


## UNIQ1 and UNIQ2
## x-axis
sources = c()

## legend
trend = c()

## y-axis
counts = c()

sources = append(sources, c("rRNA", "rRNA", "tRNA", "tRNA", "miRNA", "miRNA", "piRNA", "piRNA", "Genome", "Genome", "circRNA", "circRNA"))
trend = append(trend, c(paste0("Unique to ", samp1), paste0("Unique to ", samp2), paste0("Unique to ", samp1), paste0("Unique to ", samp2), paste0("Unique to ", samp1), paste0("Unique to ", samp2), paste0("Unique to ", samp1), paste0("Unique to ", samp2), paste0("Unique to ", samp1), paste0("Unique to ", samp2), paste0("Unique to ", samp1), paste0("Unique to ", samp2)))
counts = append(counts, c(uniq1rRNA, uniq2rRNA, uniq1tRNA, uniq2tRNA, uniq1miRNA, uniq2miRNA, uniq1piRNA, uniq2piRNA, uniq1Genome, uniq2Genome, uniq1circRNA, uniq2circRNA))

tmpList <- list(sources, trend, counts)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Source','Trend','Counts')

testDF$Trend <- factor(testDF$Trend, levels=c(paste0("Unique to ", samp1), paste0("Unique to ", samp2)))
testDF$Source <- factor(testDF$Source, levels=(unique(testDF$Source)))

# bars <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "stack", stat = "identity") +
# 	scale_fill_manual(values=c("Green","Yellow")) +
# 	ggtitle("m5C Categorization Between Comparisons") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Source") +
# 	ylab("Counts")

# ggsave(paste0(outPrefix, "_Count_(UNIQ1v2).png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Count_(UNIQ1v2).pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

# barsProportion <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "fill", stat = "identity") +
# 	scale_fill_manual(values=c("Green","Yellow")) +
# 	scale_y_continuous(labels = scales::percent_format()) +
# 	ggtitle("m5C Categorization Between Comparisons (%)") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Source") +
# 	ylab("Proportion")

# ggsave(paste0(outPrefix, "_Percentage_(UNIQ1v2).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Percentage_(UNIQ1v2).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

## All
## x-axis
sources = c()

## legend
trend = c()

## y-axis
counts = c()

sources = append(sources, c("rRNA", "rRNA", "rRNA", "rRNA", "rRNA", "tRNA", "tRNA", "tRNA", "tRNA", "tRNA", "miRNA", "miRNA", "miRNA", "miRNA", "miRNA", "piRNA", "piRNA", "piRNA", "piRNA", "piRNA", "Genome", "Genome", "Genome", "Genome", "Genome", "circRNA", "circRNA", "circRNA", "circRNA", "circRNA"))
trend = append(trend, c("Up", "Down", "Unchanged", paste0("Unique to ", samp1), paste0("Unique to ", samp2), "Up", "Down", "Unchanged", paste0("Unique to ", samp1), paste0("Unique to ", samp2), "Up", "Down", "Unchanged", paste0("Unique to ", samp1), paste0("Unique to ", samp2), "Up", "Down", "Unchanged", paste0("Unique to ", samp1), paste0("Unique to ", samp2), "Up", "Down", "Unchanged", paste0("Unique to ", samp1), paste0("Unique to ", samp2), "Up", "Down", "Unchanged", paste0("Unique to ", samp1), paste0("Unique to ", samp2)))
counts = append(counts, c(uprRNA, downrRNA, unchangedrRNA, uniq1rRNA, uniq2rRNA, uptRNA, downtRNA, unchangedtRNA, uniq1tRNA, uniq2tRNA, upmiRNA, downmiRNA, unchangedmiRNA, uniq1miRNA, uniq2miRNA, uppiRNA, downpiRNA, unchangedpiRNA, uniq1piRNA, uniq2piRNA, upGenome, downGenome, unchangedGenome, uniq1Genome, uniq2Genome, upcircRNA, downcircRNA, unchangedcircRNA, uniq1circRNA, uniq2circRNA ))

## save table with all info
tab = matrix(c( uprRNA, downrRNA, unchangedrRNA, uniq1rRNA, uniq2rRNA, uptRNA, downtRNA, unchangedtRNA, uniq1tRNA, uniq2tRNA, upmiRNA, downmiRNA, unchangedmiRNA, uniq1miRNA, uniq2miRNA, uppiRNA, downpiRNA, unchangedpiRNA, uniq1piRNA, uniq2piRNA, upGenome, downGenome, unchangedGenome, uniq1Genome, uniq2Genome, upcircRNA, downcircRNA, unchangedcircRNA, uniq1circRNA, uniq2circRNA ), ncol=5, byrow=TRUE)
rownames(tab) = c( paste(compName, "rRNA"), paste(compName, "tRNA"), paste(compName, "miRNA"), paste(compName, "piRNA"), paste(compName, "Genome"), paste(compName, "circRNA") )
colnames(tab) = c( "UP", "DOWN", "UNCHANGED", "uniq1", "uniq2")
tab = as.table(tab)
write.table(tab, file=paste0(outPrefix, ".tsv"), quote=FALSE, sep='\t', col.names = FALSE)

tmpList <- list(sources, trend, counts)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Source','Trend','Counts')

testDF$Trend <- factor(testDF$Trend, levels=c("Up", "Down", paste0("Unique to ", samp1), paste0("Unique to ", samp2)))
testDF$Source <- factor(testDF$Source, levels=(unique(testDF$Source)))

# bars <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "stack", stat = "identity") +
# 	scale_fill_manual(values=c("Red", "Blue", "Green","Yellow")) +
# 	ggtitle("m5C Categorization Between Comparisons") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Source") +
# 	ylab("Counts")

# ggsave(paste0(outPrefix, "_Count_(All).png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Count_(All).pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

# barsProportion <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "fill", stat = "identity") +
# 	scale_fill_manual(values=c("Red", "Blue", "Green","Yellow")) +
# 	scale_y_continuous(labels = scales::percent_format()) +
# 	ggtitle("m5C Categorization Between Comparisons (%)") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Source") +
# 	ylab("Proportion")

# ggsave(paste0(outPrefix, "_Percentage_(All).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Percentage_(All).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")


## All Sources in one bar
## x-axis
sources = c()

## legend
trend = c()

## y-axis
counts = c()

sources = append(sources, c(compName, compName, compName, compName))
trend = append(trend, c("Up", "Down", paste0("Unique to ", samp1), paste0("Unique to ", samp2)))
counts = append(counts, c(sumUp, sumDown, sumUniq1, sumUniq2))

tmpList <- list(sources, trend, counts)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Source','Trend','Counts')

testDF$Trend <- factor(testDF$Trend, levels=c("Up", "Down", paste0("Unique to ", samp1), paste0("Unique to ", samp2)))
testDF$Source <- factor(testDF$Source, levels=(unique(testDF$Source)))

# bars <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "stack", stat = "identity") +
# 	scale_fill_manual(values=c("Red", "Blue", "Green","Yellow")) +
# 	ggtitle("m5C Categorization Between Comparisons") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Comparison") +
# 	ylab("Counts")

# ggsave(paste0(outPrefix, "_Count_(Single).png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Count_(Single).pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

# barsProportion <- ggplot(testDF,
#                aes(x = Source,
#                    y = Counts,
#                    fill = Trend)) +
#     geom_bar(position = "fill", stat = "identity") +
# 	scale_fill_manual(values=c("Red", "Blue", "Green","Yellow")) +
# 	scale_y_continuous(labels = scales::percent_format()) +
# 	ggtitle("m5C Categorization Between Comparisons (%)") +
# 	theme(plot.title = element_text(hjust = 0.5)) +
# 	xlab("Comparison") +
# 	ylab("Proportion")

# ggsave(paste0(outPrefix, "_Percentage_(Single).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
# ggsave(paste0(outPrefix, "_Percentage_(Single).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")
