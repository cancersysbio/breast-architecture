### CREATE EXTENDED DATA 7B #######################################################################
# calculate overlap of translocations in cyclic and complex non-cyclic amplifications
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)

## set data directory
data_path <- ""
if ((!exists("data_path")) | data_path == "") {
	stop("Error: Path for SV data direcotry not set. Please set path to SV data directory ...")
}

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}
### CALCULATE OVERLAP #############################################################################
calculate_translocation_overlap <- function(sample, ecdna, datadir) {
	print(sample)
	tmp <- ecdna[ecdna$sample == sample,]
	# read in svs
	svfile <- file.path(datadir, paste0(sample, '_input.bedpe'))
	if (file.exists(svfile)) {
		# read in sv 
		sv <- read.delim(
			svfile,
			as.is = TRUE,
			header = FALSE
			)
		if (nrow(sv) > 0) {
			# turn both breakpoints into granges 
			svgr1 <- makeGRangesFromDataFrame(sv[,c('V1','V2','V3','V7','V12')], seqnames.field = 'V1', 
				start.field = 'V2', end.field = 'V3', keep.extra.columns=TRUE)
			svgr2 <- makeGRangesFromDataFrame(sv[,c('V4','V5','V6','V7','V12')], seqnames.field = 'V4', 
				start.field = 'V5', end.field = 'V6', keep.extra.columns=TRUE)
			# iterate over amplicon 
			output <- list()
			for (amplicon in unique(tmp$ampliconID)) {
				tmp2 <- tmp[tmp$ampliconID == amplicon,]
				# reformat as gr 
				egr <- makeGRangesFromDataFrame(tmp2, keep.extra.columns=TRUE)
				# find overlaps 
				toverlap1 <- suppressWarnings(subsetByOverlaps(svgr1, egr))
				toverlap2 <- suppressWarnings(subsetByOverlaps(svgr2, egr))
				# find unique translocation ids
				tra <- unique(c(
					mcols(toverlap1[mcols(toverlap1)$V12 == 'TRA',])$V7, 
					mcols(toverlap2[mcols(toverlap2)$V12 == 'TRA',])$V7
					))
				# find unique sv ids
				allsvs <- unique(c(
					mcols(toverlap1)$V7, 
					mcols(toverlap2)$V7
					))
				# summarize overlap 
				if (length(toverlap1) > 0 | length(toverlap2) > 0) {
					output[[amplicon]] <- data.frame(
						sample = sample,
						ampliconID = amplicon,
						#intervals = length(unique(c(mcols(toverlap1)$intervalID, mcols(toverlap2)$intervalID))),
						total_intervals = length(unique(tmp2$intervalID)),
						number_SV = length(allsvs),
						number_TRA = length(tra)
						)
				} else {
					output[[amplicon]] <- data.frame(
						sample = sample,
						ampliconID = amplicon,
						#intervals = 0,
						total_intervals = length(unique(tmp2$intervalID)),
						number_SV = 0,
						number_TRA = 0
						)
					}
				}
				output <- do.call(rbind, output)
			} else {
				output <- data.frame(
					sample = sample,
					ampliconID = tmp$ampliconID,
					#intervals = 0,
					total_intervals = as.data.frame(table(tmp$ampliconID))$Freq,
					number_SV = 0,
					number_TRA = 0
					)
			}
		} else {
			output <- data.frame(
				sample = sample,
				ampliconID = NA,
				#intervals = NA,
				total_intervals = NA,
				number_SV = NA,
				number_TRA = NA
				)
		}
	return(output)
}

### CALCULATE OVERLAP WRAPPER #####################################################################
calculate_translocation_overlap_wrapper <- function(megatable, disease) {
	# read in amplicons 
	amplicons <- do.call(rbind, sapply(
		c(7,17),
		function(i) {
			tmp <- read.delim(
				file.path(main_repo_path, 'data', paste0('amplicon_segments_project', i, '_all.tsv')), 
				as.is = TRUE
				)
			return(tmp[which(tmp$sample %in% megatable$sample),])
			},
		simplify = FALSE
		))
	ecdna <- amplicons[which(amplicons$ecDNA == 'Positive'),1:8]
	complex <- amplicons[which(amplicons$amplicon_type == 'Complex non-cyclic'),1:8]
	# find overlapping translocations in cyclic amplifications
	ecdna_df <- do.call(rbind, sapply(
		unique(ecdna$sample),
		calculate_translocation_overlap,
		ecdna = ecdna,
		datadir = file.path(data_path, disease),
		simplify = FALSE
		))
	ecdna_df$prop_TRA <- ecdna_df$number_TRA/ecdna_df$number_SV
	ecdna_df$type <- 'cyclic'
	ecdna_df$type_code <- 'a'
	# find overlapping translocations in complex non cyclic amplifications
	complex_df <- do.call(rbind, sapply(
		unique(complex$sample),
		calculate_translocation_overlap,
		ecdna = complex,
		datadir = file.path(data_path, disease),
		simplify = FALSE
		))
	complex_df$prop_TRA <- complex_df$number_TRA/complex_df$number_SV
	complex_df$type <- 'complex'
	complex_df$type_code <- 'b'
	# merge data
	plot_data <- rbind(ecdna_df, complex_df)
	plot_data <- merge(plot_data, megatable[,c('sample','group','ENiClust')], by = 'sample')
	plot_data$lognumber_TRA <- log10(plot_data$number_TRA+1)
	plot_data <- merge(plot_data, megatable[,c('sample','ER')], by = 'sample')
	plot_data[which(plot_data$group == 'HER2+' & plot_data$ER == 1),'group'] <- 'HER2+/ERa'
	plot_data[which(plot_data$group == 'HER2+' & plot_data$ER == 0),'group'] <- 'HER2+/ERb'
	return(plot_data)
}
### MAIN ##########################################################################################
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

# create primary plot data 
primary_plot_data <- calculate_translocation_overlap_wrapper(
	megatable = megatable,
	disease = 'primary'
	)
metastatic_plot_data <- calculate_translocation_overlap_wrapper(
	megatable = met_megatable,
	disease = 'metastatic'
	)

### WRITE TO FILE #################################################################################
# write primary data to file
write.table(
	primary_plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData7b_sourcetable_primary.txt'),
	sep = '\t',
	row.names = FALSE
	)
# write metastatic data to file
write.table(
	metastatic_plot_data,
	file = file.path(main_repo_path, 'data', 'SupplementaryFigure6b_sourcetable_metastatic.txt'),
	sep = '\t',
	row.names = FALSE
	)