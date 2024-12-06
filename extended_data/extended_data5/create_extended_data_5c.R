### CREATE EXTENDED DATA 5C #######################################################################
# create oncogene plot for metastatic tumors

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
	ecdna <- merge(ecdna, megatable[,c('sample','group','iC10','ER')], by ='sample')
	ecdna$subtype <- ecdna$iC10
	ecdna[which(ecdna$eniclust %in% c('Other','ic8') | (ecdna$eniclust == 'ic4' & ecdna$ER == 1)),'subtype'] <- 'ER+ Low'
	ecdna[which((ecdna$eniclust == 'ic4' & ecdna$ER == 0)),'subtype'] <- 'ic4ER-'
	ecdna[which((ecdna$eniclust == 'ic5')),'subtype'] <- 'HER2+'
	return(ecdna)
}

### FIND COSMIC GENES IN ECDNA ####################################################################
find_cosmic_genes_in_ecdna <- function(ecdna, genes_in_cytobands) {
	# find oncogenes
	cosmic <- read.delim(
		file.path(main_repo_path, 'data', 'Cosmic_CancerGeneCensus_v98_GRCh37.tsv.gz'),
		as.is = TRUE
		)
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
	has_onco <- unique(ecdna[,c('sample','ampliconID','group','oncogene')])
	return(has_onco)
}


### MAIN ##########################################################################################
# read in amplicons 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
}
amplicons <- read.delim(amplicons_files, as.is = TRUE)
ecdna_all <- amplicons[amplicons$ecDNA == 'Positive',]

### METASTATIC ####################################################################################
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))

# add cytoband annotations 
cyout <- add_cytoband_annotations(ecdna_all)
cytobands <- cyout$cytobands
genes_in_cytobands <- cyout$genes_cyto
genes_on_chromosome <- cyout$genes_chr

# add subtype 
ecdna <- merge(ecdna_all, megatable[,c("sample","group")], by = "sample")

# find cosmic genes
has_onco <- find_cosmic_genes_in_ecdna(ecdna, genes_in_cytobands)
oncogene_dist <- as.data.frame(table(has_onco[,c('oncogene','group')]))
oncogene_dist$Prop <- NA 
for (i in unique(oncogene_dist$group)) {
	oncogene_dist[which(oncogene_dist$group == i),'Prop'] <- oncogene_dist[which(oncogene_dist$group == i),'Freq']/sum(oncogene_dist[which(oncogene_dist$group == i),'Freq'])
}

# reformat groups so plot in wanted order
oncogene_dist$oncocode <- NA
oncogene_dist[oncogene_dist$oncogene == 'IC oncogene','oncocode'] <- 'a'
oncogene_dist[oncogene_dist$oncogene == 'Non-IC oncogene','oncocode'] <- 'b'
oncogene_dist[oncogene_dist$oncogene == 'No oncogene','oncocode'] <- 'c'
oncogene_dist <- oncogene_dist[order(oncogene_dist$group, oncogene_dist$oncocode),]

oncogene_dist$oncocode <- paste(oncogene_dist$group, oncogene_dist$oncocode, sep = '_')

# create plot
create.barplot(
	Prop ~ group,
	groups = oncocode,
	stack = TRUE,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	xaxis.cex = 1.1,
	#xaxis.rot = 45,
	main = 'Metastatic',
	xaxis.lab = c('ER+\nHigh','ER+\nTypical','HER2+','IC10','IC4ER-'),
	main.cex = 2,
	legend = list(
             right = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('white', 'darkgrey','black')
                             ),
                         text = list(
                             lab = c('No\nOncogene','Alternative\nOncogene','IC-specific\nOncogene')
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     )
                 )
             ),
	add.text = TRUE,
	text.labels = oncogene_dist$Freq[c(3,2,1,6,5,4,9,8,7,12,11,10,15,14,13)],
	text.x = rep(1:5, each = 3),
	text.y = c(0.97, 0.79, 0.38, 0.98, 0.7, 0.2, 0.99,
		0.9, 0.5, 0.9090909, 0.45, .07, 0.9, 0.52, 0.1),
	text.col = c(rep('black',8),'white','black','black','white',rep('black',3)),
	#text.col = rep(c('black','black','white'), 5),
	data = oncogene_dist,
	col = c(
		'#fa954eff','#fccaa7','white',
		'#a1c1d4ff','#d0e0e9','white',
		'#8b0100ff','#c58087','white',
		'#7c26ccff','#7c26cc6c','white',
		'#c2b7dbff','#e6e2ff','white'
		),
	#col = rev(c('white','darkgrey','black')),
	ylab.label = 'Proportion',
	xlab.label = 'Subtype',
	filename = 'extended_data5c.pdf',
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	oncogene_dist,
	file = file.path(main_repo_path, 'data', 'ExtendedData5c_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)


