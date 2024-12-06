### GENERATE EXTENDED DATA 7D ############################################################################
# create boxplot of ER and APOBEC3B overlap with ecDNA

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### CALCULATE OVERLAP #############################################################################
# calculate overlap 
calculate_overlap <- function(regions, ecdnas) {
	regions <- regions[!is.na(regions$start) & !is.na(regions$end),]
	# create genomic ranges for regions 
	regionsgr <- makeGRangesFromDataFrame(
		regions,
		ignore.strand = TRUE
		)
	# create genomic ranges for ecdna 
	ecdnagr <- makeGRangesFromDataFrame(
		ecdnas[,c('sample','ampliconID','chr','start','end')],
		keep.extra.columns=TRUE,
		ignore.strand = TRUE
		)
	# calculate overlaps 
	ecdnarloops <- findOverlaps(ecdnagr, regionsgr, ignore.strand = TRUE)
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
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))
megatable <- megatable[which(megatable$group %in% c('ER+ High','HER2+')),]

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

# set up list of files
chipseq <- list()
chipseq[['APOBEC3B']] <- file.path(main_repo_path, 'data', 'APOBEC3B_merged.bed')
chipseq[['ESR1']] <- file.path(main_repo_path, 'data', 'ESR1_merged.bed')

res <- list()
stats <- list()
for (tf in c('APOBEC3B','ESR1')) {
	bed <- read.delim(
		chipseq[[tf]],
		header = FALSE,
		col.names = c('chr','start','end','id','width','x'),
		as.is = TRUE
		)
	bed$chr <- gsub('chr','', bed$chr)
	cyclic <- calculate_overlap(bed,  pecdna)
	cyclic$tf <- tf
	cyclic$type <- 'cyclic'
	noncyclic <- calculate_overlap(bed,  pnoncyclic)
	noncyclic$tf <- tf
	noncyclic$type <- 'noncyclic'
	
	res[[tf]] <- rbind(cyclic, noncyclic)
	# calculate stats
	statstmp <- wilcox.test(
		cyclic$FreqNorm,
		noncyclic$FreqNorm
		)
	stats[[tf]] <- data.frame(
		tf = tf,
		es = median(cyclic$FreqNorm)-median(noncyclic$FreqNorm),
		p = statstmp$p.value
		)
}
res <- do.call(rbind, res)
stats <- do.call(rbind, stats)

res$group <- paste(res$tf, res$type, sep = '_')
create.boxplot(
	FreqNorm ~ group,
	data = res,
	add.stripplot = TRUE,
	xaxis.cex = 1,
	ylimits = c(0,30),
	yat = seq(0,30,10),
	ylab.label = 'Number peaks (per Mbp)',
	xlab.label = 'Amplicon Type',
	#xaxis.rot = 45,
	xaxis.lab = c('Cyclic','Non-cyclic','Cyclic','Non-cyclic'),
	add.text = TRUE,
	text.labels = c(
		'APOBEC3B','ESR1',
		expression('ES=1.9'), expression('P=0.003'),
		expression('ES=1.1'), expression('P=0.02')
		),
	add.rectangle = TRUE,
	xleft.rectangle = c(0,2.5),
	ybottom.rectangle = 0,
	xright.rectangle = c(2.5,6),
	ytop.rectangle = 35,
	#ytop.rectangle = 10,
	col.rectangle = c('darkseagreen4','mediumpurple2'),
	alpha.rectangle = 0.25,
	text.x = c(1.5,3.5, 1.5, 1.5, 3.5, 3.5),
	#text.y = c(1.95,1.95,1.87, 1.8, 1.87, 1.8),
	text.y = c(29,29,28, 27, 28, 27),
	text.anchor = 'centre',
	text.col = 'black',
	text.cex = 1,
	text.fontface = 'plain',
	filename = 'extended_data7d.pdf',
	resolution = 300
	)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  res,
  file = file.path(main_repo_path, 'data', 'ExtendedData7d_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )