### SUPPLEMENTARY FIGURE 5B #######################################################################
# find overlap between SV clusters and ecdna 

### PREAMBLE ######################################################################################
library(GenomicRanges)

## set data directory
data_path <- ""
if ((!exists("data_path")) | data_path == "") {
	stop("Error: Path for SV data direcotry not set. Please set path to SV data directory ...")
}

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}
### FIND OVERLAP WITH SV ########################################################################## 
find_overlap_with_sv <- function(sample, ecdna, data_path, threshold = 2) {
	print(sample)
	filename <- file.path(data_path, paste0(sample, '.sv_clusters_and_footprints.tsv'))
	if (file.exists(filename)) {
		svcluster <- read.delim(
			filename,
			as.is = TRUE,
			header = FALSE
			)
		# only keep clusters with three or more SVs 
		svcluster <- svcluster[svcluster$V12 > threshold,]
		if (nrow(svcluster) > 0) {
			# create granges
			clustgr <- c(
				makeGRangesFromDataFrame(
					svcluster,
					seqnames.field = 'V1',
					start.field = 'V2',
					end.field = 'V3'
					),
				makeGRangesFromDataFrame(
					svcluster,
					seqnames.field = 'V4',
					start.field = 'V5',
					end.field = 'V6'
					)
				)
			# create ecdna granges 
			ecdnagr <- makeGRangesFromDataFrame(
					ecdna[ecdna$sample == sample,c('sample','chr','start','end','ampliconID')],
					seqnames.field = 'chr',
					start.field = 'start',
					end.field = 'end',
					keep.extra.columns=TRUE
					)
			# find overlaps 
			overlaps <- subsetByOverlaps(ecdnagr, clustgr)
			# return ecdna that overlap 
			return(unique(as.data.frame(mcols(overlaps))))
		}
	} else {
		return(data.frame(sample = sample, ampliconID = NA))
		}
}

### TEST THRESHOLDS ###############################################################################
test_thresolds <- function(ecdna, threshold, data_path) {
	# for each sample read in sv clusters and get for overlap 
	eoverlap <- do.call(rbind, sapply(
		unique(ecdna$sample),
		find_overlap_with_sv,
		ecdna = ecdna,
		threshold = threshold,
		data_path = data_path,
		simplify = FALSE
		))
	eoverlap$sv <- TRUE
	# find how many ecdna overlap sv
	ecdna_amp <- unique(ecdna[,c('sample','ampliconID')])
	ecdna_amp <- merge(ecdna_amp, eoverlap, by = c('sample','ampliconID'), all.x = TRUE)
	eprop <- sum(ecdna_amp$sv == TRUE, na.rm = TRUE)/nrow(ecdna_amp)
	res <- data.frame(
		threshold = threshold, 
		proportion = eprop
		)
	return(res)
}
### MAIN ##########################################################################################
# read in amplicons 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
}
amplicons <- read.delim(
	amplicons_file,
	as.is = TRUE
	)
ecdna <- amplicons[amplicons$ecDNA == 'Positive',]
complex <- amplicons[amplicons$amplicon_type == 'Complex non-cyclic',]

# for each sample read in sv clusters and get for overlap 
eoverlap <- do.call(rbind, sapply(c(2,4,9), test_thresolds, ecdna = ecdna, data_path = data_path, simplify = FALSE))
coverlap <- do.call(rbind, sapply(c(2,4,9), test_thresolds, ecdna = complex, data_path = data_path, simplify = FALSE))

eoverlap$type <- 'ecdna'
coverlap$type <- 'complex'

# create plot 
plot_data <- rbind(eoverlap, coverlap)
plot_data$index <- 1:nrow(plot_data)
create.barplot(
	proportion ~ index,
	data = plot_data,
	ylimits = c(0,1),
	filename = 'supplementary_figure5b.pdf',
	ylab.label = 'Proportion',
	xlab.label = 'Min # SVs in Cluster',
	xaxis.lab = c('3','5','10','3','5','10'),
	add.rectangle = TRUE,
	xleft.rectangle = 3.5,
	ybottom.rectangle = 0,
	xright.rectangle = 7,
	ytop.rectangle = 1,
	add.text = TRUE,
	yat = seq(0,1,0.2),
	text.labels = c('Cyclic','Non-Cyclic'),
	text.cex = 0.8,
	text.y = 0.97,
	text.x = c(2,5),
	col.rectangle = 'grey50',
	alpha.rectangle = 0.5,
	resolution = 300
	)

### WRITE TO FILE #################################################################################
# write to file 
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data','SupplementaryFigure5b_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)



