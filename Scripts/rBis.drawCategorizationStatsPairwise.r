#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('optparse', quiet=TRUE))
suppressPackageStartupMessages(library('ggplot2', quiet=TRUE))
suppressPackageStartupMessages(library('data.table', quiet=TRUE))

# command line option handling
option_list <- list(
	make_option(c("-o","--outPrefix"), default="Differential_Analysis_Categorization", help="Output prefix. default=Differential_Analysis_Categorization"),
    make_option(c("-m","--mqcOut"), default="Differential_Analysis_Categorization", help="Output prefix. default=Differential_Analysis_Categorization")
)
parser <- OptionParser(usage = "%prog [options] <tsv>", option_list=option_list,
			description = "Draw Categorization stats for each pairwise comparison.")
arguments <- parse_args(parser, positional_arguments = TRUE)
if(length(arguments$args) == 0) {
	print_help(parser)
	q()
} else {
	srcL=arguments$args
}

print(srcL)

# Option handling
opt=arguments$options
outPrefix=opt$outPrefix
mqcOut=opt$mqcOut

## UP v DOWN
comps = c()
trend = c()
counts = c()
allTables = c()

for( src in srcL ){
	statsTable <- read.table(file = src, sep = '\t', header = TRUE)

    allTables = rbind(allTables, statsTable)

    comp = statsTable$Comparison
    up = statsTable$UP
    down = statsTable$DOWN
    unchanged = statsTable$UNCHANGED
    
    comps = append(comps, comp)
    comps = append(comps, comp)
    trend = append(trend, "Up")
    trend = append(trend, "Down")
    counts = append(counts, up)
    counts = append(counts, down)
}

write.table(allTables, file=paste0(mqcOut, ".tsv"), quote=FALSE, sep='\t', row.names = FALSE)

tmpList <- list(comps, trend, counts)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Comparison','Trend','Counts')

testDF$Trend <- factor(testDF$Trend, levels=c('Up', 'Down'))
testDF$Comparison <- factor(testDF$Comparison, levels=(unique(testDF$Comparison)))

bars <- ggplot(testDF,
               aes(x = Comparison,
                   y = Counts,
                   fill = Trend)) +
    geom_bar(position = "stack", stat = "identity") +
	scale_fill_manual(values=c("Red","Blue")) +
	ggtitle("m5C Categorization Between Comparisons") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Counts")

ggsave(paste0(outPrefix, "_Count_(UPvDOWN).png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Count_(UPvDOWN).pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

barsProportion <- ggplot(testDF,
               aes(x = Comparison,
                   y = Counts,
                   fill = Trend)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("Red","Blue")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("m5C Categorization Between Comparisons (%)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_Percentage_(UPvDOWN).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Percentage_(UPvDOWN).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")


## UNIQ1 v UNIQ2
comps = c()
trend = c()
counts = c()

for( src in srcL ){
	statsTable <- read.table(file = src, sep = '\t', header = TRUE)

    comp = statsTable$Comparison
    uniq1 = statsTable$uniq1
    uniq2 = statsTable$uniq2
    
    comps = append(comps, comp)
    comps = append(comps, comp)
    trend = append(trend, "Uniq1")
    trend = append(trend, "Uniq2")
    counts = append(counts, uniq1)
    counts = append(counts, uniq2)
}

tmpList <- list(comps, trend, counts)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Comparison','Trend','Counts')

testDF$Trend <- factor(testDF$Trend, levels=c('Uniq1', 'Uniq2'))
testDF$Comparison <- factor(testDF$Comparison, levels=(unique(testDF$Comparison)))

bars <- ggplot(testDF,
               aes(x = Comparison,
                   y = Counts,
                   fill = Trend)) +
    geom_bar(position = "stack", stat = "identity") +
	scale_fill_manual(values=c("Green","Yellow")) +
	ggtitle("m5C Categorization Between Comparisons") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Counts")

ggsave(paste0(outPrefix, "_Count_(UNIQ1v2).png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Count_(UNIQ1v2).pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

barsProportion <- ggplot(testDF,
               aes(x = Comparison,
                   y = Counts,
                   fill = Trend)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("Green","Yellow")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("m5C Categorization Between Comparisons (%)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_Percentage_(UNIQ1v2).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Percentage_(UNIQ1v2).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")


## all
comps = c()
trend = c()
counts = c()

for( src in srcL ){
	statsTable <- read.table(file = src, sep = '\t', header = TRUE)

    comp = statsTable$Comparison
    up = statsTable$UP
    down = statsTable$DOWN
    uniq1 = statsTable$uniq1
    uniq2 = statsTable$uniq2
    
    comps = append(comps, comp)
    comps = append(comps, comp)
    comps = append(comps, comp)
    comps = append(comps, comp)
    trend = append(trend, "Up")
    trend = append(trend, "Down")
    trend = append(trend, "Uniq1")
    trend = append(trend, "Uniq2")
    counts = append(counts, up)
    counts = append(counts, down)
    counts = append(counts, uniq1)
    counts = append(counts, uniq2)
}

tmpList <- list(comps, trend, counts)

testDF = data.table::copy(tmpList)
setDT(testDF)

colnames(testDF) <- c('Comparison','Trend','Counts')

testDF$Trend <- factor(testDF$Trend, levels=c('Up', 'Down', 'Uniq1', 'Uniq2'))
testDF$Comparison <- factor(testDF$Comparison, levels=(unique(testDF$Comparison)))

bars <- ggplot(testDF,
               aes(x = Comparison,
                   y = Counts,
                   fill = Trend)) +
    geom_bar(position = "stack", stat = "identity") +
	scale_fill_manual(values=c("Red", "Blue", "Green","Yellow")) +
	ggtitle("m5C Categorization Between Comparisons") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Counts")

ggsave(paste0(outPrefix, "_Count_(All).png"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Count_(All).pdf"), bars, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")

barsProportion <- ggplot(testDF,
               aes(x = Comparison,
                   y = Counts,
                   fill = Trend)) +
    geom_bar(position = "fill", stat = "identity") +
	scale_fill_manual(values=c("Red", "Blue", "Green","Yellow")) +
	scale_y_continuous(labels = scales::percent_format()) +
	ggtitle("m5C Categorization Between Comparisons (%)") +
	theme(plot.title = element_text(hjust = 0.5)) +
	xlab("Comparison") +
	ylab("Proportion")

ggsave(paste0(outPrefix, "_Percentage_(All).png"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "png")
ggsave(paste0(outPrefix, "_Percentage_(All).pdf"), barsProportion, width = 5, height = 4, dpi = 150, units = "in", device = "pdf")