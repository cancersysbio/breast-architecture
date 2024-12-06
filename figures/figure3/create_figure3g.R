### GENERATE FIGURE 3G ############################################################################
# create boxplot of ER-induced r-loop enrichment in ecDNA in primary tumors

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### CALCULATE R LOOPS #############################################################################
# calculate overlap 
calculate_rloops_overlap <- function(rloops, ecdnas) {
	rloops <- rloops[!is.na(rloops$start) & !is.na(rloops$end),]
	# create genomic ranges for rloops 
	rloopsgr <- makeGRangesFromDataFrame(
		rloops,
		ignore.strand = TRUE
		)
	# create genomic ranges for ecdna 
	ecdnagr <- makeGRangesFromDataFrame(
		ecdnas[,c('sample','ampliconID','chr','start','end')],
		keep.extra.columns=TRUE,
		ignore.strand = TRUE
		)
	# calculate overlaps 
	ecdnarloops <- findOverlaps(ecdnagr, rloopsgr, ignore.strand = TRUE)
	# find number of overlaps per segment
	overlapcounts <- as.data.frame(table(factor(queryHits(ecdnarloops), levels = 1:length(ecdnagr))))
	overlapcounts$sample <- as.data.frame(mcols(ecdnagr))$sample
	overlapcounts$ampliconID <- as.data.frame(mcols(ecdnagr))$ampliconID
	overlapcounts$width <- as.data.frame(ranges(ecdnagr))$width
	# calculate number of ecdnas that overlap 
	overlapcounts_peramplicon <- aggregate(overlapcounts[,c('Freq','width')], overlapcounts[,c('sample','ampliconID')], sum)
	overlapcounts_peramplicon$FreqNorm <- (overlapcounts_peramplicon$Freq/overlapcounts_peramplicon$width)*1e6
	# return
	return(overlapcounts_peramplicon)
}

### MAIN ##########################################################################################
# read in differential R loops from GSE81851
rloops_file <- file.path(main_repo_path, 'data', 'GSE81851_DESeq_differential_calls.tsv.gz')
if (!file.exists(rloops_file)) {
	stop("Please download DRIPSeq DESeq differential calls from GSE81851 and put in data directory ...")
}
# read in differential R loops 
rloops <- read.delim(
	rloops_file, 
	as.is = TRUE
	)
rloops$chr <- gsub('chr','',rloops$chr)

# find sig regions 
t24vs0 <- rloops[rloops$padj.T0_T24_DRIP < 0.05 & rloops$log2FoldChange.T0_T24_DRIP > 2,]

# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))
# read in megatable 
met_megatable <- read.delim(
	file.path(main_repo_path, 'data','metastatic_megatable.txt'),
	as.is = TRUE
	)
colnames(met_megatable) <- gsub('Sample','sample', colnames(met_megatable))

# read in amplicons 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
} 
amplicons <- read.delim(amplicons_files, as.is = TRUE)
pamplicons <- amplicons[which(amplicons$sample %in% megatable$sample),]
# split into categories
pecdna <- pamplicons[which(pamplicons$ecDNA == 'Positive'),] 
pnoncyclic <- pamplicons[which(pamplicons$amplicon_type == 'Complex non-cyclic'),] 
# extract amplicons in metastatic tumours
mamplicons <- amplicons[which(amplicons$sample %in% met_megatable$sample),]
# split into categories
mecdna <- mamplicons[which(mamplicons$ecDNA == 'Positive'),] 
mnoncyclic <- mamplicons[which(mamplicons$amplicon_type == 'Complex non-cyclic'),] 

# calculate overlap 
pt24ecdna <- calculate_rloops_overlap(t24vs0, pecdna)
pt24ecdna$comparison <- 'ecdna'
pt24noncyclic <- calculate_rloops_overlap(t24vs0, pnoncyclic)
pt24noncyclic$comparison <- 'noncyclic'

pplot_data <- rbind(pt24ecdna, pt24noncyclic)
pplot_data$stage <- 'a'

mt24ecdna <- calculate_rloops_overlap(t24vs0, mecdna)
mt24ecdna$comparison <- 'ecdna'
mt24noncyclic <- calculate_rloops_overlap(t24vs0, mnoncyclic)
mt24noncyclic$comparison <- 'noncyclic'

mplot_data <- rbind(mt24ecdna, mt24noncyclic)
mplot_data$stage <- 'b'

plot_data <- rbind(pplot_data, mplot_data)
plot_data$index <- paste(plot_data$stage, plot_data$comparison)
# add subtype 
plot_data <- merge(
	plot_data, 
	rbind(megatable[,c('sample','group')], met_megatable[,c('sample', 'group')]),
	by = 'sample'
	)

res <- list()
for (i in c('a','b')) {
	tmp <- plot_data[plot_data$stage == i,]
	stats <- wilcox.test(
		tmp[tmp$comparison == 'ecdna','FreqNorm'],
		tmp[tmp$comparison == 'noncyclic','FreqNorm']
		)
	res[[i]] <- data.frame(
		stage = i,
		es = median(tmp[tmp$comparison == 'ecdna','FreqNorm'])-median(tmp[tmp$comparison == 'noncyclic','FreqNorm']),
		p = stats$p.value
		)
}
res <- do.call(rbind, res)

# plot
create.boxplot(
	FreqNorm ~ index,
	data = plot_data,
	add.stripplot = TRUE,
	xaxis.cex = 1,
	xaxis.lab = c('Cyclic','Non-cyclic','Cyclic','Non-cyclic'),
	ylimits = c(0,2),
	yat = seq(0,2,0.5),
	ylab.label = 'Number R-loops (per Mbp)',
	xlab.label = 'Amplicon Type',
	add.text = TRUE,
	text.labels = c(
		'Primary','Metastatic',
		expression('ES=0.10'), expression('P=1.7x10'^'-6'),
		expression('ES=0.08'), expression('P=1.6x10'^'-7')
		),
	add.rectangle = TRUE,
	xleft.rectangle = 2.5,
	ybottom.rectangle = 0,
	xright.rectangle = 6,
	ytop.rectangle = 35,
	#ytop.rectangle = 10,
	col.rectangle = 'gold',
	alpha.rectangle = 0.25,
	text.x = c(1.5,3.5, 1.5, 1.5, 3.5, 3.5),
	text.y = c(1.95,1.95,1.87, 1.8, 1.87, 1.8),
	text.anchor = 'centre',
	text.col = 'black',
	text.cex = 1,
	text.fontface = 'bold',
	filename = 'figure3g.pdf',
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'Figure3g_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

