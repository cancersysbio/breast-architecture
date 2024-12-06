### CREATE FIGURE 3I #############################################################################
# create figure 3i
# calculate distance of ER-induced R-loops to IC amplicon genes

### PREAMBLE #####################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)
library(biovizBase)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### FIND COSMIC GENES IN ECDNA ####################################################################
find_cosmic_genes_in_ecdna <- function(ecdna, cosmic) {
	reuda <- c('RPS6KB1','ZNF703','FGF3','PAK1','RSF1','NARS2')
	# find cosmic genes
	ecdna$oncogene <- apply(
		ecdna,
		1,
		function(x) {
			tmp <- strsplit(x['oncogenes'], ',')[[1]]
			onco <- tmp[which(tmp %in% c(cosmic$GENE_SYMBOL,reuda))]
			if (length(onco) > 0) {
				return(paste(onco, collapse = '|'))
			} else {
					return('No oncogene')
				}
			})
	has_onco <- unique(ecdna[,c('sample','ampliconID','oncogene','ENiClust','amplicon_length')])
	return(has_onco)
}

### CALCULATE RECURRENCE IN ECDNA #################################################################
calculate_gene_recurrence_in_ecdna <- function(cosmic, megatable) {
	# read in amplicons 
	amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
	if (!file.exists(amplicons_file)) {
		stop("Please put amplicon segments as outputted by AA in data directory ...")
	}
	amplicons <- read.delim(amplicons_files, as.is = TRUE)
	amplicons <- amplicons[which(amplicons$sample %in% megatable$sample),]
	# find ecdna
	ecdna <- amplicons[which(amplicons$ecDNA == 'Positive'),] 
	ecdna <- merge(ecdna, megatable[,c('sample','ENiClust')], by = 'sample', all.x = TRUE)
	has_onco <- find_cosmic_genes_in_ecdna(ecdna = ecdna, cosmic = cosmic)
	# calculate recurrence
	recurrence <- as.data.frame(table(unlist(sapply(
		has_onco$oncogene,
		function(x) {
			tmp <- unlist(strsplit(x, '\\|')[[1]])
			return(tmp)
			}))))
	colnames(recurrence) <- c('gene','recurrence')
	out[['recurrence']] <- recurrence
	out[['has_onco']] <- has_onco
	return(out)
	}

### CALCULATE DISTANCE TO ELEMENT #################################################################
calculate_distance_to_element <- function(genesymbol_cosmic, element_gr, recurrence, has_onco) {
	# calculate distance to element
	edist_tmp <- distanceToNearest(genesymbol_cosmic, element_gr)
	# create plot data
	plot_data <- data.frame(
		gene = mcols(genesymbol_cosmic)$symbol[queryHits(edist_tmp)],
		distance = mcols(edist_tmp)$distance
		)
	# calculate min per gene coordinates
	plot_data <- aggregate(plot_data$distance, list(plot_data$gene), min)
	colnames(plot_data) <- c('gene','distance')
	plot_data$index <- 'all'
	plot_data <- merge(plot_data, recurrence, by = 'gene', all.x = TRUE)
	ecdnanum <- as.data.frame(table(has_onco$ENiClust))
	# estimate recurrence proportion
	plot_data$ecdna <- NA
	for (i in 1:length(iconcogenes)) {
		plot_data[which(plot_data$gene %in% iconcogenes[[i]]),'ecdna'] <- ecdnanum[ecdnanum$Var1 == tolower(names(iconcogenes[i])),'Freq']
	}
	plot_data$prop <- plot_data$recurrence/plot_data$ecdna
	return(plot_data)
}

### SET UP COLS, CEX AND ALPHA ####################################################################
set_up_cols_cex_alpha <- function(plot_data, recurrence) {
	# set up colours and alpha
	col <- rep('darkgrey', nrow(plot_data))
	col[which(plot_data$gene %in% iconcogenes[['IC1']] & plot_data$prop >= 0.5)] <- '#ff5500ff'
	col[which(plot_data$gene %in% iconcogenes[['IC2']] & plot_data$prop >= 0.5)] <- '#00ee77ff'
	col[which(plot_data$gene %in% iconcogenes[['IC5']] & plot_data$prop >= 0.5)] <- '#8b0100ff'
	col[which(plot_data$gene %in% iconcogenes[['IC6']] & plot_data$prop >= 0.5)] <- '#fffe41ff'
	col[which(plot_data$gene %in% iconcogenes[['IC9']] & plot_data$prop >= 0.5)] <- '#ee82edff'

	cex <- rep(0.5, nrow(plot_data))
	cex[which(plot_data$prop >= 0.5 & plot_data$gene %in% unlist(iconcogenes))] <- 1

	alpha <- rep(0.5, nrow(plot_data))
	alpha[which(plot_data$prop >= 0.5 & plot_data$gene %in% unlist(iconcogenes))] <- 1
	# return plotting parameters
	out <- list()
	out[['col']] <- col
	out[['cex']] <- cex
	out[['alpha']] <- alpha
	return(out)
}

### MAIN ##########################################################################################
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))

# find oncogenes
cosmic_file <- file.path(main_repo_path, 'data', 'Cosmic_CancerGeneCensus_v98_GRCh37.tsv.gz')
if (!file.exists(rloops_file)) {
	stop("Please download Cosmic Cancer Genes v98 GRCh37 from COSMIC and put in data directory ...")
}
cosmic <- read.delim(
		cosmic_file,
		as.is = TRUE
		)
reuda <- c('RPS6KB1','ZNF703','FGF3','PAK1','RSF1','NARS2')
# set IC specific
iconcogenes <- list(
	IC5 = c('ERBB2','CDK12','IKZF3'),
	IC1 = c('DDX5','CLTC','PPM1D','BRIP1','CD79B','RPS6KB1'),
	IC2 = c('FADD','CCND1','NUMA1','MEN1','PAK1','RSF1','FGF3','NARS2','FAT3','TAF15'),
	IC6 = c('FGFR1','IKBKB','HOOK3','KAT6A','NRG1','WRN','ZNF703'),
	IC9 = c('MYC','EXT1','RECQL4','NDRG1','FAM135B','RAD21')
	)

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
t24vs0 <- rloops[rloops$padj.T0_T24_DRIP < 0.05 & rloops$log2FoldChange.T0_T24_DRIP > 2,]
# convert to granges
t24vs0gr <- makeGRangesFromDataFrame(
		t24vs0,
		ignore.strand = TRUE
		)

# calculate recurrence
out <- calculate_gene_recurrence_in_ecdna(cosmic = cosmic, megatable = megatable)
recurrence <- out$recurrence
has_onco <- out$has_onco

# get gene coordinates 
data(genesymbol)
genesymbol_cosmic <- genesymbol[mcols(genesymbol)$symbol %in% c(cosmic$GENE_SYMBOL,reuda),]

## MEASURE DISTANCE TO RLOOPS #####################################################################
# calculate distance to element
plot_data <- calculate_distance_to_element(
	genesymbol_cosmic = genesymbol_cosmic,
	element_gr = t24vs0gr,
	recurrence = recurrence, 
	has_onco = has_onco
	)

# set up colours and alpha
plot_param <- set_up_cols_cex_alpha(plot_data, recurrence)

ic9 <- 1-(sum(plot_data[plot_data$gene == 'MYC','distance'] < plot_data$distance)/nrow(plot_data))
ic2 <- 1-(sum(plot_data[plot_data$gene == 'PAK1','distance'] < plot_data$distance)/nrow(plot_data))
ic6 <- 1-(sum(plot_data[plot_data$gene == 'ZNF703','distance'] < plot_data$distance)/nrow(plot_data))
ic1 <- 1-(sum(plot_data[plot_data$gene == 'RPS6KB1','distance'] < plot_data$distance)/nrow(plot_data))
ic5 <- 1-(sum(plot_data[plot_data$gene == 'ERBB2','distance'] < plot_data$distance)/nrow(plot_data))

plot_data$logdistance <- log10(plot_data$distance+1)
create.boxplot(
	logdistance ~ index,
	data = plot_data,
	add.stripplot = TRUE,
	points.col = plot_param$col,
	jitter.factor = 7,
	ylimits = c(4,9),
	yat = seq(4,9,2),
	yaxis.lab = c(
		expression(10^4),
		expression(10^6),
		expression(10^8)
		),
	ylab.label = 'Distance to R-loop (Mbp)',
	xlab.label = '',
	xaxis.lab = 'Genes',
	add.text = TRUE,
	text.labels = c('MYC','ZNF703','PAK1','ERBB2','RPS6KB1',
		expression('P'['IC9']*'=0.07'),
		expression('P'['IC6']*'=0.02'),
		expression('P'['IC2']*'=0.08'),
		expression('P'['HER2+']*'=0.40'),
		expression('P'['IC1']*'=0.66')
		),
	text.col = c(
		'#ee82edff','gold', '#00ee77ff','#8b0100ff','#ff5500ff',
		rep('black',5)
		),
	text.y = c(
		plot_data[plot_data$gene == 'MYC','logdistance'],
		plot_data[plot_data$gene == 'ZNF703','logdistance'],
		plot_data[plot_data$gene == 'PAK1','logdistance'],
		plot_data[plot_data$gene == 'ERBB2','logdistance'],
		plot_data[plot_data$gene == 'RPS6KB1','logdistance'],
		8.8,8.5,8.2,7.9,7.6
		),
	text.x = c(1.25, 1.25, 1.25, 1.25, 1.25, rep(0.55, 3),0.575, 0.55),
	text.cex = 1,
	text.fontface = 4,
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 1.5,
                             fill = c('#ff5500ff','#00ee77ff','#8b0100ff','#fffe41ff','#ee82edff')
                             ),
                         text = list(
                             lab = c('IC1','IC2','IC5','IC6','IC9')
                             ),
                         cex = 2
                         )
                     ),
                 x = 0.99,
                 y = 0.99,
                 corner = c(1,1),
                 draw = FALSE
                 )
             ),
	points.cex = plot_param$cex,
	points.alpha = plot_param$alpha,
	filename = 'figure3i.pdf',
	height = 12,
	width = 5,
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'Figure3i_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
