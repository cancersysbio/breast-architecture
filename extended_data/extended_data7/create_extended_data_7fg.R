### CREATE EXTENDED DATA 7FG ######################################################################
# create boxplot of R-loop density in cyclic vs non-cyclic across subgroups 

### MAIN ##########################################################################################
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

### CALCULATE OVERLAP IN PRIMARY AND METASTATIC ###################################################
calculate_overlap_primary_metastatic <- function(pamplicons, mamplicons, rloops) {
	# split into categories
	pecdna <- pamplicons[which(pamplicons$ecDNA == 'Positive'),] 
	pnoncyclic <- pamplicons[which(pamplicons$amplicon_type == 'Complex non-cyclic'),] 
	mecdna <- mamplicons[which(mamplicons$ecDNA == 'Positive'),] 
	mnoncyclic <- mamplicons[which(mamplicons$amplicon_type == 'Complex non-cyclic'),] 

	# calculate overlap primary 
	ptecdna <- calculate_rloops_overlap(rloops, pecdna)
	ptecdna$comparison <- 'ecdna'
	ptnoncyclic <- calculate_rloops_overlap(rloops, pnoncyclic)
	ptnoncyclic$comparison <- 'noncyclic'

	pplot_data <- rbind(ptecdna, ptnoncyclic)
	pplot_data$stage <- 'primary'

	# calculate overlap primary 
	mtecdna <- calculate_rloops_overlap(rloops, mecdna)
	mtecdna$comparison <- 'ecdna'
	mtnoncyclic <- calculate_rloops_overlap(rloops, mnoncyclic)
	mtnoncyclic$comparison <- 'noncyclic'

	mplot_data <- rbind(mtecdna, mtnoncyclic)
	mplot_data$stage <- 'metastatic'
	# create plot data 
	plot_data <- rbind(pplot_data, mplot_data)
	return(plot_data)
}


### CALCULATE RLOOP OVERLAP DIFFERENCE ############################################################
calculate_rloop_overlap_difference <- function(plot_data) {
	res <- list()
	for (i in c('primary','metastatic')) {
		for (j in c('IC10','ER+ High','ER+ Typical','HER2+','all')) {
			if (j == 'all') {
				tmp <- plot_data[which(plot_data$stage == i),]
			} else {
				tmp <- plot_data[which(plot_data$stage == i & plot_data$group == j),]
			} 
		stats <- wilcox.test(
			tmp[tmp$comparison == 'ecdna','FreqNorm'],
			tmp[tmp$comparison == 'noncyclic','FreqNorm']
			)
		res[[paste0(i, j)]] <- data.frame(
			stage = i,
			subtype = j,
			es = median(tmp[tmp$comparison == 'ecdna','FreqNorm'])-median(tmp[tmp$comparison == 'noncyclic','FreqNorm']),
			p = stats$p.value
			)
		}
	}
	res <- do.call(rbind, res)
	return(res)
}

### PLOT R LOOP OVERLAP ###########################################################################
plot_rloop_overlap <- function(plot_data, xaxis.lab, text.labels, text.x, text.y,
	filename, ylimits, yat, xright.rectangle, xleft.rectangle, width = 6, xaxis.rot = 0) {
		create.boxplot(
			FreqNorm ~ index,
			data = plot_data,
			add.stripplot = TRUE,
			xaxis.cex = 1,
			xaxis.lab = xaxis.lab,
			ylimits = ylimits,
			yat = yat,
			ylab.label = 'Number R-loops (per Mbp)',
			xlab.label = 'Amplicon Type',
			xaxis.rot = xaxis.rot,
			add.text = TRUE,
			text.labels = text.labels,
			add.rectangle = TRUE,
			xleft.rectangle = xleft.rectangle,
			ybottom.rectangle = 0,
			xright.rectangle = xright.rectangle,
			ytop.rectangle = ylimits[2]+10,
			col.rectangle = 'gold',
			alpha.rectangle = 0.25,
			text.x = text.x,
			text.y = text.y,
			text.anchor = 'centre',
			text.col = 'black',
			text.cex = 1.1,
			text.fontface = 'bold',
			filename = filename,
			width = width,
			resolution = 300
			)

}

### RLOOPS ########################################################################################
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
t2vsinput <- rloops[rloops$padj.T2_DRIP_over_Input < 0.05 & rloops$log2FoldChange.T2_DRIP_over_Input > 2,]
t24vsinput <- rloops[rloops$padj.T24_DRIP_over_Input < 0.05 & rloops$log2FoldChange.T24_DRIP_over_Input > 2,]
t2vs0 <- rloops[rloops$padj.T0_T2_DRIP < 0.05 & rloops$log2FoldChange.T0_T2_DRIP > 2,]
t24vs0 <- rloops[rloops$padj.T0_T24_DRIP < 0.05 & rloops$log2FoldChange.T0_T24_DRIP > 2,]

### MEGATABLE #####################################################################################
# read in megatable 
megatable <- read.delim(
		file.path(main_repo_path, 'data', 'primary_megatable.txt'),
		as.is = TRUE
		)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))
# read in megatable 
met_megatable <- read.delim(
		file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
		as.is = TRUE
		)
colnames(met_megatable) <- gsub('Sample','sample', colnames(met_megatable))

### MAIN ##########################################################################################
# read in amplicons 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
}
amplicons <- read.delim(amplicons_files, as.is = TRUE)
pamplicons <- amplicons[which(amplicons$sample %in% megatable$sample),]
mamplicons <- amplicons[which(amplicons$sample %in% met_megatable$sample),]

# calculate overlap at 10
plot_data_t24 <- calculate_overlap_primary_metastatic(
	pamplicons,
	mamplicons,
	t24vs0
	)
# add subtype 
plot_data_t24 <- merge(
	plot_data_t24, 
	rbind(megatable[,c('sample','group','iC10','ER')], met_megatable[,c('sample', 'group','iC10','ER')]),
	by = 'sample'
	)
plot_data_t24$index <- paste0(plot_data_t24$stage, plot_data_t24$group, plot_data_t24$comparison)

# create primary and metastatic plot data
plot_data_t24_primary <- plot_data_t24[which(plot_data_t24$stage == 'primary' & plot_data_t24$group %in% c('ER+ High','ER+ Typical','HER2+','IC10')),]
# read in metastatic source data
plot_data_t24_metastatic <- plot_data_t24[which(plot_data_t24$stage == 'metastatic' & plot_data_t24$group %in% c('ER+ High','ER+ Typical','HER2+','IC10')),]

# calculate stats and create plots
plot_data_t24_res <- calculate_rloop_overlap_difference(rbind(plot_data_t24_primary, plot_data_t24_metastatic))
plot_rloop_overlap(
	plot_data_t24_primary,
	xaxis.lab = rep(c('Cyclic','Non-cyclic'), 5),
	text.labels = c('ER+ High','  ER+\nTypical','HER2+','IC10',
		expression('ES=0.19'),expression('P=3.0x10'^'-3'), expression('ES=0.11'),
		expression('P=3.2x10'^'-4')),
	text.x = c(seq(1.5,8,2),1.5,1.5,3.5,3.5),
	text.y = c(rep(2.1,4),rep(c(1.98, 1.9),2)),
	filename = 'extended_data7f.pdf', 
	ylimits = c(0,2.3),
	yat = seq(0, 2, 0.5),
	xright.rectangle = c(0,4.5), 
	xleft.rectangle = c(2.5,6.5),
	#width = 7,
	xaxis.rot = 45
	)
plot_rloop_overlap(
	plot_data_t24_metastatic,
	xaxis.lab = rep(c('Cyclic','Non-cyclic'), 5),
	text.labels = c('ER+ High','  ER+\nTypical','HER2+','IC10',
		expression('ES=0.18'),expression('P=1.9x10'^'-6'), expression('ES=0.08'),
		expression('P=0.02')),
	text.x = c(seq(1.5,8,2),1.5,1.5,5.5,5.5),
	text.y = c(rep(2,4),rep(c(1.88, 1.8),2)),
	filename = 'extended_data7g.pdf', 
	ylimits = c(0,2.2),
	yat = seq(0, 2, 0.5),
	xright.rectangle = c(0,4.5), 
	xleft.rectangle = c(2.5,6.5),
	#width = 7,
	xaxis.rot = 45
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data_t24_primary,
	file = file.path(main_repo_path, 'data', 'ExtendedData7f_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
write.table(
	plot_data_t24_metastatic,
	file = file.path(main_repo_path, 'data', 'ExtendedData7g_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
