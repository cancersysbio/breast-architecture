### CREATE FIGURE 3H #############################################################################
# create figure 3h, percent BP overlap with ER-induced r-loops

### PREAMBLE #####################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### CALCULATE STD ERR #############################################################################
std.err <- function(x) sd(x)/sqrt(length(x))

### CALCUALTE BP OVERLAP WITH R-LOOPS #############################################################
calculate_bp_overlap_with_rloops <- function(bpfile, rloopsgr) {
	# read in data
	tmp <- read.delim(
		bpfile,
		as.is = TRUE,
		col.names = c('chromosome','start','end','info','x','strand'),
		header = FALSE
		)
	tmp$start <- tmp$start-5000
	tmp$end <- tmp$end+5000
	# create genomic ranges for ecdna 
	svgr <- makeGRangesFromDataFrame(
			tmp
			)
	# calculate overlaps 
	svrloops <- subsetByOverlaps(svgr, rloopsgr)
	# format output data 
	parsefile <- strsplit(bpfile, '_')[[1]]
	out <- data.frame(
			cell = parsefile[3],
			condition = parsefile[4],
			location = parsefile[5],
			replicate = gsub('.bed.gz', '', parsefile[6]),
			overlap = length(svrloops),
			bps = nrow(tmp)
			)
	return(out)
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

# find sig regions 
t24vs0 <- rloops[rloops$padj.T0_T24_DRIP < 0.05 & rloops$log2FoldChange.T0_T24_DRIP > 2,]
# create genomic ranges for rloops 
rloopsgr <- makeGRangesFromDataFrame(
		t24vs0,
		strand.field="str"
		)

# iterate over breakpoint files to calculate overlap with rloops
if (!dir.exists(main_repo_path, 'data','GSE227369_RAW')) {
	stop("Please download HTGTS data from GSE227369. Download GSE227369_RAW folder and put folder in data directory ...")
}
files <- list.files(
	path = file.path(main_repo_path, 'data','GSE227369_RAW'), 
	pattern = 'bed.gz'
	)
out <- do.call(rbind, sapply(
	files, 
	calculate_bp_overlap_with_rloops,
	rloopsgr = rloopsgr,
	simplify = FALSE
	))
# keep only MCF7
out <- out[out$cell == 'MCF7',]
# calculte percent overlap
out$prop <- out$overlap/out$bps*100
# calculate median and standard deviation
outagg_median <- aggregate(out$prop, out[,c('cell','condition','location')], median)
outagg_sd <- aggregate(out$prop, out[,c('cell','condition','location')], sd)
outagg <- merge(outagg_median, outagg_sd, by = c('cell','condition','location'))
	colnames(outagg) <- c('cell','condition','location','median','se')
# calculate stats
p <- wilcox.test(out[out$condition == 'withoutE2','prop'], out[out$condition == 'withE2','prop'], paired = TRUE)$p.value
es <- median(out[out$condition == 'withE2','prop'])/median(out[out$condition == 'withoutE2','prop'])

# plot barplot
outagg$index <- paste(outagg$cell, outagg$location, outagg$condition, sep = '_')
create.barplot(
		median ~ index,
		data = outagg,
		ylimits = c(0,0.22),
		yat = seq(0,0.2,0.05),
		xaxis.lab = rep(c('+', '-'),2),
		xlab.label = 'E2 Treatment',
		ylab.label = 'Percent BP Overlap',
		y.error.up = outagg$se,
		col = '#e2b2b2',
		y.error.bar.col = 'black',
		abline.v = 2.5,
		abline.lty = 2,
		error.whisker.width = 0.1,
		add.text = TRUE,
		text.labels = c('RARA','SHANK2'),
		text.x = c(1.5,3.5),
		text.y = c(0.21,0.21),
		main = paste0('FC=', round(es, digits = 2),'; P=', round(p, digits = 3)),
		main.cex = 1.5,
		filename = 'figure3h.pdf',
		resolution = 300
		)

### SAVE ##########################################################################################
# save source data
write.table(
	outagg,
	file = file.path(main_repo_path, 'data', 'Figure3h_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)