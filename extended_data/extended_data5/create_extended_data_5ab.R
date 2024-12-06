### CREATE EXTENDED DATA 5AB ######################################################################
# create barplot of recurrent oncogenes in primary and metastatic 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)
library(biovizBase)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### ADD CYTOBAND ANNOTATIONS ######################################################################
add_cytoband_annotations <- function(ecdna) {
	# find all genes 
	allgenes <- unique(unlist(sapply(
				ecdna$oncogenes,
				function(y) {
					unlist(strsplit(y, ','))
				})))
	data(hg19IdeogramCyto)
	data(genesymbol)
	cytobands <- do.call(rbind, sapply(
		allgenes[allgenes != ''],
		function(x) {
			print(x)
			if (any(mcols(genesymbol)$symbol == x)) {
				tmp <- subsetByOverlaps(hg19IdeogramCyto, genesymbol[x,])
				return(data.frame(gene = x, cytoband = paste0(as.character(seqnames(tmp)), mcols(tmp)$name[[1]])))
			} else {
				return(data.frame(gene = x, cytoband = NA))
				}
			},
		simplify = FALSE
		))
	cytobands <- unique(cytobands)
	cytobands$cytoband_subset <- sapply(cytobands$cytoband, function(x) strsplit(as.character(x), '\\.')[[1]][1])
	cytobands$chromosome <- sapply(cytobands$cytoband, function(x) strsplit(as.character(x), 'p|q')[[1]][1])
	cytobands_interest <- c('chr8q24','chr17q12','chr17q23','chr8p12','chr8p11','chr11q13')
	genes_in_cytobands <- cytobands[which(cytobands$cytoband_subset %in% cytobands_interest),]
	genes_on_chromosome <- cytobands[which(cytobands$chromosome %in% c('chr8','chr11','chr17')),]
	# create output 
	out <- list()
	out[['genes_cyto']] <- genes_in_cytobands
	out[['genes_chr']] <- genes_on_chromosome
	out[['cytobands']] <- cytobands
	return(out)
}

### ADD SUBTYPE ###################################################################################
add_subtype <- function(ecdna, megatable) {
	ecdna <- merge(ecdna, megatable[,c('sample','group','ENiClust','ER')], by ='sample')
	ecdna$subtype <- ecdna$ENiClust
	ecdna[which(ecdna$ENiClust %in% c('Other','ic8') | (ecdna$ENiClust == 'ic4' & ecdna$ER == 1)),'subtype'] <- 'ER+ Typical'
	ecdna[which((ecdna$ENiClust == 'ic4' & ecdna$ER == 0)),'subtype'] <- 'ic4ER-'
	ecdna[which((ecdna$ENiClust == 'ic5')),'subtype'] <- 'HER2+'
	return(ecdna)
}

### FIND RECURRENT GENES ##########################################################################
find_recurrent_genes <- function(ecdna) {
	recurrence <- do.call(rbind, sapply(
		unique(ecdna$subtype)[!is.na(unique(ecdna$subtype))],
		function(x) {
			tmp <- unique(ecdna[which(ecdna$subtype == x),c('sample','oncogenes')])
			tmpgenes <- unlist(sapply(
				tmp$oncogenes,
				function(y) {
					unlist(strsplit(y, ','))
				}))
			tmpdf <- as.data.frame(table(tmpgenes))
			colnames(tmpdf) <- c('gene','freq')
			tmpdf$prop <- tmpdf$freq/nrow(unique(ecdna[which(ecdna$subtype == x),c('sample','ampliconID')]))
			tmpdf$group <- x
			return(tmpdf[tmpdf$gene != '',])
			},
		simplify = FALSE
		))
	return(recurrence)
}

### CREATE PLOT DATA ##############################################################################
create_plot_data <- function(ecdna, recurrence, cosmic) {
	# find oncogenes
	recurrence_cosmic <- recurrence[which(recurrence$gene %in% cosmic$GENE_SYMBOL),]
	recurrence_cosmic10 <- recurrence_cosmic[which(recurrence_cosmic$prop > 0.2),]
	recurrence_cosmic10$group <- factor(recurrence_cosmic10$group, levels = c('ER+ Typical','HER2+','ic4ER-','ic10','ic1','ic2','ic6','ic9'))
	plot_data <- as.data.frame(table(recurrence_cosmic10$group))
	colnames(plot_data) <- c('group','count')

	# calculate avg length per subtype 
	avg_length <- aggregate(ecdna$amplicon_length, list(ecdna$subtype), median)
	colnames(avg_length) <- c('group','length')
	plot_data <- merge(plot_data, avg_length, by = 'group')
	plot_data$norm_count <- (plot_data$count/plot_data$length)*1e6

	plot_data <- plot_data[order(-plot_data$norm_count),]
	plot_data$index <- 1:nrow(plot_data)
	return(plot_data)
}

### MAIN ##########################################################################################
# read in segments 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
}
amplicons <- read.delim(amplicons_files, as.is = TRUE)
ecdna_all <- amplicons[amplicons$ecDNA == 'Positive',]

# add cytoband annotations 
cyout <- add_cytoband_annotations(ecdna_all)
cytobands <- cyout$cytobands
genes_in_cytobands <- cyout$genes_cyto
genes_on_chromosome <- cyout$genes_chr

# find oncogenes
cosmic <- read.delim(
		file.path(main_repo_path, 'data', 'Cosmic_CancerGeneCensus_v98_GRCh37.tsv.gz'),
		as.is = TRUE
		)

### PRIMARY 5A ####################################################################################
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data','primary_megatable.txt'),
	as.is = TRUE
	)
#megatable <- add_four_main_groups(megatable)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))

# add subtype 
ecdna <- add_subtype(ecdna_all, megatable)

# find recurrence in each subtype 
recurrence <- find_recurrent_genes(ecdna)

# create plot data
plot_data <- create_plot_data(ecdna = ecdna, recurrence = recurrence, cosmic = cosmic)
# only two IC4ER- samples with ecdna so removing
plot_data <- plot_data[plot_data$group != 'ic4ER-',]

create.barplot(
		norm_count ~ index,
		data = plot_data,
		xaxis.lab = gsub('^ic','IC', plot_data$group),
		xaxis.rot = 45,
		xaxis.cex = 1.2,
		ylimits = c(0,1),
		yat = seq(0,1,0.25),
		#col = default.colours(10)[c(9:10,5)],
		add.rectangle = TRUE,
		xleft.rectangle = 1.5,
		ybottom.rectangle = 0,
		xright.rectangle = 5.5,
		ytop.rectangle = 10,
		col.rectangle = 'grey50',
		alpha.rectangle = 0.5,
		main = 'Primary',
		main.cex = 2,
		ylab.label = 'Number of Recurrent\nOncogenes per Mbp',
		xlab.label = 'Subtype',
		add.text = TRUE,
		text.labels = 'ER+ High',
		text.x = 3.5,
		text.y = 0.95,
		filename = 'extended_data5a.pdf',
		resolution = 300
		)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData5a_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### METASTATIC 5B #################################################################################
# read in megatable 
metmegatable <- read.delim(
	file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
	as.is = TRUE
	)
colnames(metmegatable) <- gsub('Sample','sample', colnames(metmegatable))

# add subtype 
metecdna <- add_subtype(ecdna_all, metmegatable)

# find recurrence in each subtype 
metrecurrence <- find_recurrent_genes(metecdna)

# find oncogenes
plot_data <- create_plot_data(metecdna, metrecurrence, cosmic)
# all have recurrence rate = 0, so reordering to make labeling easier
plot_data[plot_data$group == 'ic9','index'] <- 5
plot_data[plot_data$group == 'ER+ Typical','index'] <- 6
plot_data[plot_data$group == 'ic10','index'] <- 7
plot_data <- plot_data[order(plot_data$index),]

create.barplot(
		norm_count ~ index,
		data = plot_data,
		xaxis.lab = gsub('^ic','IC', plot_data$group),
		xaxis.rot = 45,
		xaxis.cex = 1.2,
		ylimits = c(0,1),
		yat = seq(0,1,0.25),
		#col = default.colours(10)[c(9:10,5)],
		add.rectangle = TRUE,
		xleft.rectangle = c(0,3.5),
		ybottom.rectangle = 0,
		xright.rectangle = c(2.5,5.5),
		ytop.rectangle = 10,
		main = 'Metastatic',
		main.cex = 2,
		col.rectangle = 'grey50',
		alpha.rectangle = 0.5,
		ylab.label = 'Number of Recurrent\nOncogenes per Mbp',
		xlab.label = 'Subtype',
		add.text = TRUE,
		text.labels = c('ER+ High','ER+ High'),
		text.x = c(1.5,4.5),
		text.y = rep(0.95,2),
		filename ='extended_data5b.pdf',
		resolution = 300
		)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData5b_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)


