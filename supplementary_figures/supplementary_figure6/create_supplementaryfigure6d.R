### CREATE SUPPLEMENTARY FIGURE 6B ################################################################
# create boxplot of R-loop delta with and without A3B knockdown

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)
library(rtracklayer)

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
	# import chain file
	chainfile <- file.path(main_repo_path, 'data', 'hg19ToHg38.over.chain')
	if (!file.exists(chainfile)) {
		stop("Error: please download hg19ToHg38.over.chain from UCSC and put it in data directory ...")
	}
	chain <- import.chain(chainfile)
	rloopsgr_hg19 <- liftOver(rloopsgr, chain)
	# create genomic ranges for ecdna 
	ecdnas$chr <- paste0('chr', ecdnas$chr)
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

compare_wt_ko_baseline_treatment <- function(ko_treatment, wt_treatment, ko_baseline, wt_baseline, ecdna, noncyclic) {
	# calculate rloop overlap 
	ko_dtf <- calculate_rloops_overlap(ko_treatment, ecdna)
	wt_dtf <- calculate_rloops_overlap(wt_treatment, ecdna)
	ko_baseline_dtf <- calculate_rloops_overlap(ko_baseline, ecdna)
	wt_baseline_dtf <- calculate_rloops_overlap(wt_baseline, ecdna)

	ko_complex_dtf <- calculate_rloops_overlap(ko_treatment, noncyclic)
	wt_complex_dtf <- calculate_rloops_overlap(wt_treatment, noncyclic)
	ko_complex_baseline_dtf <- calculate_rloops_overlap(ko_baseline, noncyclic)
	wt_complex_baseline_dtf <- calculate_rloops_overlap(wt_baseline, noncyclic)

	# create boxplot 
	plot_data <- merge(ko_dtf, wt_dtf, by = c('sample','ampliconID'))
	colnames(plot_data) <- gsub('.x','_ko', colnames(plot_data), fixed = TRUE)
	colnames(plot_data) <- gsub('.y','_wt', colnames(plot_data), fixed = TRUE)

	plot_data_baseline <- merge(ko_baseline_dtf, wt_baseline_dtf, by = c('sample','ampliconID'))
	colnames(plot_data_baseline) <- gsub('.x','_ko', colnames(plot_data_baseline), fixed = TRUE)
	colnames(plot_data_baseline) <- gsub('.y','_wt', colnames(plot_data_baseline), fixed = TRUE)


	# create boxplot 
	plot_data_complex <- merge(ko_complex_dtf, wt_complex_dtf, by = c('sample','ampliconID'))
	colnames(plot_data_complex) <- gsub('.x','_ko', colnames(plot_data_complex), fixed = TRUE)
	colnames(plot_data_complex) <- gsub('.y','_wt', colnames(plot_data_complex), fixed = TRUE)

	plot_data_complex_baseline <- merge(ko_complex_baseline_dtf, wt_complex_baseline_dtf, by = c('sample','ampliconID'))
	colnames(plot_data_complex_baseline) <- gsub('.x','_ko', colnames(plot_data_complex_baseline), fixed = TRUE)
	colnames(plot_data_complex_baseline) <- gsub('.y','_wt', colnames(plot_data_complex_baseline), fixed = TRUE)

	plot_data_both <- data.frame(
		Freq = c(
			plot_data$Freq_ko-plot_data$Freq_wt,
			plot_data_baseline$Freq_ko-plot_data_baseline$Freq_wt,
			plot_data_complex$Freq_ko-plot_data_complex$Freq_wt,
			plot_data_complex_baseline$Freq_ko-plot_data_complex_baseline$Freq_wt
			),
		Experiment = c(rep('ecdna',nrow(plot_data)), rep('ecdna_baseline', nrow(plot_data_baseline)), 
			rep('complex', nrow(plot_data_complex)), rep('complex_baseline', nrow(plot_data_complex_baseline)))
		)
	ecdna_pvalue <- wilcox.test(plot_data$Freq_ko-plot_data$Freq_wt, plot_data_baseline$Freq_ko-plot_data_baseline$Freq_wt)$p.value
	ecdna_es <- median(plot_data$Freq_ko-plot_data$Freq_wt)-median(plot_data_baseline$Freq_ko-plot_data_baseline$Freq_w)
	complex_pvalue <- wilcox.test(plot_data_complex$Freq_ko-plot_data_complex$Freq_wt, plot_data_complex_baseline$Freq_ko-plot_data_complex_baseline$Freq_wt)$p.value
	complex_es <- median(plot_data_complex$Freq_ko-plot_data_complex$Freq_wt)-median(plot_data_complex_baseline$Freq_ko-plot_data_complex_baseline$Freq_wt)
	out <- list()
	out$plot_data <- plot_data_both
	out$ecdna_stats <- c(ecdna_es, ecdna_pvalue)
	out$complex_stats <- c(complex_es, complex_pvalue)
	return(out)
}

### MAIN ##########################################################################################
# set file names
dripseq_files <- list()
dripseq_files[['ko']] <- file.path(main_repo_path, 'ecDNA/MCF10A_DRIPSeq_McCann/', 'GSE148581_DRIP_KO_PMA_6h_Called_Peaks.txt.gz')
dripseq_files[['wt']] <- file.path(main_repo_path, 'ecDNA/MCF10A_DRIPSeq_McCann/', 'GSE148581_DRIP_WT_PMA_6h_Called_Peaks.txt.gz')
dripseq_files[['ko_baseline']] <- file.path(main_repo_path, 'ecDNA/MCF10A_DRIPSeq_McCann/', 'GSE148581_DRIP_KO_DMSO_Called_Peaks.txt.gz')
dripseq_files[['wt_baseline']] <- file.path(main_repo_path, 'ecDNA/MCF10A_DRIPSeq_McCann/', 'GSE148581_DRIP_WT_DMSO_Called_Peaks.txt.gz')

# check if files exist
if (!all(file.exists(unlist(dripseq_files)))) {
	stop("Please download DRIP-Seq 'Called Peaks' files from GSE148581 and put them in data directory ...")
}
# read in drip seq 
ko <- read.delim(
	dripseq_files$ko,
	as.is = TRUE
	)
wt <- read.delim(
	dripseq_files$wt,
	as.is = TRUE
	)
ko_baseline <- read.delim(
	dripseq_files$ko_baseline,
	as.is = TRUE
	)
wt_baseline <- read.delim(
	dripseq_files$wt_baseline,
	as.is = TRUE
	)
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))
megatable <- megatable[megatable$group %in% c('HER2+','ER+ High'),]

# read in amplicons 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
}
amplicons <- read.delim(amplicons_files, as.is = TRUE)
pamplicons <- amplicons[which(amplicons$sample %in% megatable$sample),]
# find ecdna
pecdna <- pamplicons[which(pamplicons$ecDNA == 'Positive'),] 
pnoncyclic <- pamplicons[which(pamplicons$amplicon_type == 'Complex non-cyclic'),] 

primary_out <- compare_wt_ko_baseline_treatment(ko_treatment = ko, wt_treatment = wt, ko_baseline = ko_baseline,
	wt_baseline = wt_baseline, ecdna = pecdna, noncyclic = pnoncyclic)
plot_data_both <- primary_out$plot_data
pecdna_stats <- primary_out$ecdna_stats
pcomplex_stats <- primary_out$complex_stats

plot_data_both$code <- NA
plot_data_both[plot_data_both$Experiment == 'ecdna','code'] <- 'a'
plot_data_both[plot_data_both$Experiment == 'ecdna_baseline','code'] <- 'b'
plot_data_both[plot_data_both$Experiment == 'complex','code'] <- 'c'
plot_data_both[plot_data_both$Experiment == 'complex_baseline','code'] <- 'd'

create.boxplot(
	Freq ~ code,
	data = plot_data_both,
	add.stripplot = TRUE,
	filename = 'supplementary_figure6d.pdf',
	xaxis.lab = c('Cyclic\nPMA 6h','Cyclic\nBaseline','Non-cylic\nPMA 6h','Non-cylic\nBaseline'),
	ylab.label = 'Difference in Number of R-loops\n(AB3 KO - AB3 WT)',
	xaxis.cex = 1,
	ylimits = c(-350,350),
	yat = seq(-300,300,100),
	add.text = TRUE,
	text.labels = c(
		expression('ES=31.0'), expression('P=1.93x10'^'-7'),
		expression('ES=0.0'), expression('P=0.50')
		),
	text.x = c(1.5,1.5,3.5,3.5),
	text.y = c(335,310,335,310),
	text.anchor = 'centre',
	text.col = 'black',
	text.cex = 1,
	text.fontface = 'bold',
	add.rectangle = TRUE,
	xleft.rectangle = 0,
	ybottom.rectangle = -500,
	xright.rectangle = 2.5,
	ytop.rectangle = 500,
	col.rectangle = 'gold',
	alpha.rectangle = 0.25,
	#xaxis.rot = 45,
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data_both,
	file = file.path(main_repo_path, 'data', 'SupplementaryFigure6d_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
