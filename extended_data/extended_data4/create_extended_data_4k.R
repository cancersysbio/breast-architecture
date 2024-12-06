### EXTENDED DATA 4K ##############################################################################
# create barplot of oncogenes recurrently included in ER+ Typical risk tumors

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
	ecdna <- merge(ecdna, megatable[,c('sample','group','eniclust','ER')], by ='sample')
	ecdna$subtype <- ecdna$eniclust
	ecdna[which(ecdna$eniclust %in% c('Other','ic8') | (ecdna$eniclust == 'ic4' & ecdna$ER == 1)),'subtype'] <- 'ER+ Low'
	ecdna[which((ecdna$eniclust == 'ic4' & ecdna$ER == 0)),'subtype'] <- 'ic4ER-'
	ecdna[which((ecdna$eniclust == 'ic5')),'subtype'] <- 'HER2+'
	return(ecdna)
}


### FIND COSMIC GENES IN ECDNA ####################################################################
find_cosmic_genes_in_ecdna <- function(ecdna, genes_in_cytobands, cosmic) {
	# find cosmic genes
	ecdna$oncogene <- apply(
		ecdna,
		1,
		function(x) {
			tmp <- strsplit(x['oncogenes'], ',')[[1]]
			onco <- tmp[which(tmp %in% cosmic$GENE_SYMBOL)]
			if (length(onco) > 0) {
				if (any(onco %in% genes_in_cytobands$gene)) {
					return('IC oncogene')
				} else {
					return('Non-IC oncogene')
				}
			} else {
					return('No oncogene')
				}
			})
	ecdna$ic_oncogene <- apply(
		ecdna,
		1,
		function(x) {
			tmp <- strsplit(x['oncogenes'], ',')[[1]]
			onco <- tmp[which(tmp %in% cosmic$GENE_SYMBOL)]
			if (length(onco) > 0) {
				if (any(onco %in% genes_in_cytobands$gene)) {
					return(paste(onco, collapse = ','))
				} else {
					return('Non-IC oncogene')
				}
			} else {
					return('No oncogene')
				}
			})
	return(ecdna)
}


### MAIN ##########################################################################################
# read in amplicons 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
}
amplicons <- read.delim(amplicons_files, as.is = TRUE)
ecdna_all <- amplicons[amplicons$ecDNA == 'Positive',]

### PRIMARY #######################################################################################
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))

# find oncogenes
cosmic <- read.delim(
		file.path(main_repo_path, 'data', 'Cosmic_CancerGeneCensus_v98_GRCh37.tsv.gz'),
		as.is = TRUE
		)

# add cytoband annotations 
cyout <- add_cytoband_annotations(ecdna_all)
cytobands <- cyout$cytobands
genes_in_cytobands <- cyout$genes_cyto
genes_on_chromosome <- cyout$genes_chr

# add subtype
ecdna <- merge(ecdna_all, megatable[,c('sample','group')], by = 'sample')

# find cosmic genes
ecdna <- find_cosmic_genes_in_ecdna(ecdna, genes_in_cytobands, cosmic = cosmic)
# keep er typical
erlow_genes <- as.data.frame(table(unlist(sapply(
	unique(ecdna[ecdna$group == 'ER+ Typical' & ecdna$oncogene == 'IC oncogene',c('sample','ampliconID','ic_oncogene')])$ic_oncogene,
	function(x) {strsplit(x, ',')[[1]]}
	))))
erlow_genes <- erlow_genes[order(-erlow_genes$Freq),]
erlow_genes <- erlow_genes[erlow_genes$Freq > 2,]
erlow_genes$proportion <- erlow_genes$Freq/nrow( unique(ecdna[ecdna$group == 'ER+ Typical',c('sample','ampliconID')]))
erlow_genes$index <- 1:nrow(erlow_genes)
create.barplot(
	proportion ~ index,
	data = erlow_genes,
	xaxis.lab = erlow_genes$Var1,
	xaxis.rot = 45,
	xaxis.cex = 1,
	xlab.label = 'Genes',
	ylab.label = 'Proportion of ER+ Typical\necDNA',
	ylimits = c(0,0.3),
	yat = seq(0,0.3,0.1),
	filename = 'extended_data4k.pdf',
	#width = 10,
	resolution = 300 
	)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  erlow_genes,
  file = file.path(main_repo_path, 'data', 'ExtendedData4k_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )
