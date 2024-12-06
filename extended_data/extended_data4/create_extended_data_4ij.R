### CREATE FIGURE 4JK ##############################################################################
# create barplot of number of oncogenes per Mbp across subgroups

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

### FIND COSMIC GENES IN ECDNA ####################################################################
find_cosmic_genes_in_ecdna <- function(ecdna, cosmic) {
	# find cosmic genes
	ecdna$oncogene <- apply(
		ecdna,
		1,
		function(x) {
			tmp <- strsplit(x['oncogenes'], ',')[[1]]
			onco <- tmp[which(tmp %in% cosmic$GENE_SYMBOL)]
			if (length(onco) > 0) {
				return(paste(onco, collapse = '|'))
			} else {
					return('No oncogene')
				}
			})
	has_onco <- unique(ecdna[,c('sample','ampliconID','subtype','oncogene','amplicon_length')])
	return(has_onco)
}

### ADD SUBTYPE ###################################################################################
add_subtype <- function(ecdna, megatable) {
	ecdna <- merge(ecdna, megatable[,c('sample','group','ENiClust','ER')], by ='sample')
	ecdna$subtype <- ecdna$ENiClust
	ecdna[which(ecdna$ENiClust %in% c('Other','ic8') | (ecdna$ENiClust == 'ic4' & ecdna$ER == 1)),'subtype'] <- 'ER+ Low'
	ecdna[which((ecdna$ENiClust == 'ic4' & ecdna$ER == 0)),'subtype'] <- 'ic4ER-'
	ecdna[which((ecdna$ENiClust == 'ic5' & ecdna$ER == 1)),'subtype'] <- 'HER2+/ER+'
	ecdna[which((ecdna$ENiClust == 'ic5' & ecdna$ER == 0)),'subtype'] <- 'HER2+/ER-'
	return(ecdna)
}

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
		file.path(main_repo_path, 'data','Cosmic_CancerGeneCensus_v98_GRCh37.tsv.gz'),
		as.is = TRUE
		)

# add subtype 
ecdna <- add_subtype(ecdna_all, megatable)

# find cosmic genes
has_onco <- find_cosmic_genes_in_ecdna(ecdna, cosmic = cosmic)
has_onco$count <- sapply(
	has_onco$oncogene,
	function(x) {
		tmp <- unlist(strsplit(x, '\\|'))
		return(ifelse(any(tmp == 'No oncogene'), 0, length(tmp)))
		})
has_onco$count_mbp <- has_onco$count/(has_onco$amplicon_length/1e6)

cyto_map <- c('ic6','ic1','ic2','HER2+','ic9','ic6')
names(cyto_map) <- c('chr8p11','chr17q23','chr11q13','chr17q12','chr8q24','chr8p12')

# find co-amp IC2/IC6
ic26 <- has_onco[has_onco$subtype %in% c('ic2','ic6'),]
ic26_coamp <- ic26[grepl('FGFR1|ZNF703|IKBKB', ic26$oncogene) & grepl('CCND1|PAK1|FADD',ic26$oncogene),]


# add cytoband annotations 
cyout <- add_cytoband_annotations(ecdna_all)
cytobands <- cyout$cytobands
genes_in_cytobands <- cyout$genes_cyto
genes_on_chromosome <- cyout$genes_chr
# find oncogenes in cytobands
genes_in_cytobands$cosmic <- sapply(genes_in_cytobands$gene, function(x) x %in% cosmic$GENE_SYMBOL)
genes_in_cytobands$subtype <- cyto_map[as.character(genes_in_cytobands$cytoband_subset)]
oncogenes_cyto <- as.data.frame(table(genes_in_cytobands[,c('subtype','cosmic')]))
oncogenes_cyto <- oncogenes_cyto[oncogenes_cyto$cosmic == TRUE,]

oncogenes_cyto$subtype <- as.character(oncogenes_cyto$subtype)
oncogenes_cyto[oncogenes_cyto$subtype == 'HER2+','subtype'] <- 'HER2+/ER+'
oncogenes_cyto <- rbind(oncogenes_cyto, oncogenes_cyto[oncogenes_cyto$subtype == 'HER2+/ER+',])
oncogenes_cyto[6,'subtype'] <- 'HER2+/ER-'

has_onco <- has_onco[which(has_onco$subtype != 'IC4ER-'),]

# set order
ic_order <- c('HER2+/ER-','HER2+/ER+','ic1','ic9','ic6','ic2','ER+ Low','ic10')
has_onco$index <- NA
for (i in 1:length(ic_order)) { 
	has_onco[has_onco$subtype == ic_order[i],'index'] <- i
}
has_onco <- has_onco[order(has_onco$index),]
has_onco$index <- factor(has_onco$index)

col <- rep(NA, nrow(has_onco))
col[which(has_onco$subtype == 'ic1')] <- '#ff5500ff'
col[which(has_onco$subtype == 'ic2')] <- '#00ee77ff'
col[which(has_onco$subtype == 'ic6')] <- '#fffe41ff'
col[which(has_onco$subtype == 'ic9')] <- '#ee82edff'
col[which(has_onco$subtype == 'ER+ Low')] <- '#a1c1d4ff'
col[which(has_onco$subtype == 'HER2+/ER-')] <- '#8b0100ff'
col[which(has_onco$subtype == 'HER2+/ER+')] <- '#8b0100b3'
col[which(has_onco$subtype == 'ic10')] <- '#7c26ccff'
col[which(has_onco$subtype == 'ic4ER-')] <- '#c2b7dbff'

# create boxplot
create.boxplot(
	count_mbp ~ index,
	data = has_onco,
	xaxis.lab = c('HER2+/ER-','HER2+/ER+','IC1','IC9','IC6','IC2','ER+ Typical','IC10','IC4ER-'),
	xaxis.rot = 45,
	xaxis.cex = 1,
	points.col = col,
	ylab.label = 'Number Oncogenes (per Mbp)',
	ylab.cex = 1.8,
	xlab.label = 'Subtype',
	ylimits = c(0,2),
	yat = seq(0,2,0.5),
	yaxis.cex = 1,
	add.stripplot = TRUE,
	filename = 'extended_data4j.pdf',
	resolution = 300
	)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  has_onco,
  file = file.path(main_repo_path, 'data', 'ExtendedData4j_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )

### EXTENDED DATA 4K ##############################################################################
# only ER+ high and HER2+
has_onco$subtype <- as.character(has_onco$subtype)
has_onco <- merge(has_onco, oncogenes_cyto, by = 'subtype')
has_onco$ratio <- ((has_onco$count/has_onco$Freq)/has_onco$amplicon_length)*1e6

# set order
ic_order <- c('HER2+/ER-','HER2+/ER+','ic1','ic9','ic6','ic2')
has_onco$index <- NA
for (i in 1:length(ic_order)) { 
	has_onco[has_onco$subtype == ic_order[i],'index'] <- i
}
has_onco <- has_onco[order(has_onco$index),]
has_onco$index <- factor(has_onco$index)

col <- rep(NA, nrow(has_onco))
col[which(has_onco$subtype == 'ic1')] <- '#ff5500ff'
col[which(has_onco$subtype == 'ic2')] <- '#00ee77ff'
col[which(has_onco$subtype == 'ic6')] <- '#fffe41ff'
col[which(has_onco$subtype == 'ic9')] <- '#ee82edff'
col[which(has_onco$subtype == 'ER+ Low')] <- '#a1c1d4ff'
col[which(has_onco$subtype == 'HER2+/ER-')] <- '#8b0100ff'
col[which(has_onco$subtype == 'HER2+/ER+')] <- '#8b0100b3'
col[which(has_onco$subtype == 'ic10')] <- '#7c26ccff'
col[which(has_onco$subtype == 'ic4ER-')] <- '#c2b7dbff'

# create boxplot
create.boxplot(
	ratio ~ index,
	data = has_onco,
	xaxis.lab = c('HER2+/ER-\n(17q22)','HER2+/ER+\n(17q22)','IC1  \n(17q23)','IC9 \n(8q24)','IC6 \n(8p11)','IC2 \n(11q13)'),
	xaxis.rot = 45,
	xaxis.cex = 1,
	points.col = col,
	ylab.label = 'Cyclic/Cytoband Oncogenes\n(per Mbp)',
	ylab.cex = 1.8,
	xlab.label = 'Subtype (Cytoband)',
	ylimits = c(0,0.6),
	yat = seq(0,0.6,0.2),
	yaxis.cex = 1,
	add.stripplot = TRUE,
	filename = 'extended_data4k.pdf',
	resolution = 300
	)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  has_onco,
  file = file.path(main_repo_path, 'data', 'ExtendedData4k_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )

