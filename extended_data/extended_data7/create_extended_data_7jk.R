### CREATE EXTENDED DATA 7JK ######################################################################

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

### MAIN ##########################################################################################
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))

### PRIMARY 7J ####################################################################################
# read in overlapping translocations
primary_ecdna_plot <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7j_sourcetable.txt'),
	as.is = TRUE
	)

create.boxplot(
	logdistance ~ index,
	data = primary_ecdna_plot,
	add.stripplot = TRUE,
	ylab.label = 'Median Distance',
	ylimits = c(4.7,8.5),
	yat = 5:8,
	yaxis.lab = c(expression(10^5), expression(10^6), expression(10^7), expression(10^8)),
	xlab.label = 'Location Translocations to ecDNA',
	xaxis.lab = rep(c('Within','Outside'),4),
	#xaxis.rot = 45,
	xaxis.cex = 0.9,
	xlab.cex = 1.5,
	ylab.cex = 1.5,
	yaxis.cex = 1.2,
	add.rectangle = TRUE,
	xleft.rectangle = c(0,2.5,4.5,6.5),
	ybottom.rectangle = 0,
	xright.rectangle = c(2.5,4.5,6.5,10),
	ytop.rectangle = 10,
	col.rectangle = c('darkorange2','dodgerblue4','darkorange2', 'dodgerblue4'),
	alpha.rectangle = 0.25,
	add.text = TRUE,
	text.labels = c(
		'ER+ High','ER+ Typical', 'HER2+','IC10',
		expression('FC=0.44'), expression('FDR=1.7x10'^'-6'),
		expression('FC=0.46'), expression('FDR=0.08'),
		expression('FC=0.66'), expression('FDR=5.2x10'^'-3'),
		expression('FC=0.60'), expression('FDR=0.30')
		),
	text.x = c(1.5,3.5, 5.5,7.5, 1.5, 1.5, 3.5, 3.5, 5.5, 5.5, 7.5, 7.5),
	text.y = c(rep(8.3,4),rep(c(5,4.85),4)),
	text.anchor = 'centre',
	text.col = 'black',
	text.cex = 1,
	text.fontface = 'bold',
	#main = 'Primary',
	#main.cex = 2,
	filename = 'extended_data7j.pdf',
	resolution = 300 
	)

res <- list()
for (i in unique(primary_ecdna_plot$group)[!is.na(unique(primary_ecdna_plot$group))]) {
	tmp <- primary_ecdna_plot[which(primary_ecdna_plot$group == i),]
	stats <- wilcox.test(
		tmp[tmp$type == 'median_ecdna_rloops','distance'],
		tmp[tmp$type == 'median_tra_rloops','distance']
		)
	es <- median(tmp[tmp$type == 'median_ecdna_rloops','distance'], na.rm = TRUE)/median(tmp[tmp$type == 'median_tra_rloops','distance'], na.rm = TRUE)
	res[[i]] <- data.frame(
		group = i,
		es = es,
		p = stats$p.value
		)
}
res <- do.call(rbind, res)
res$fdr <- p.adjust(res$p, method = 'fdr')

#### METASTATIC 7K ###################################################################################
# read in megatable 
met_megatable <- read.delim(
	file.path(main_repo_path, 'data','metastatic_megatable.txt'),
	as.is = TRUE
	)
colnames(met_megatable) <- gsub('Sample','sample', colnames(met_megatable))

# read in overlapping translocations
met_ecdna_plot <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7k_sourcetable.txt'),
	as.is = TRUE
	)

create.boxplot(
	logdistance ~ index,
	data = met_ecdna_plot,
	add.stripplot = TRUE,
	ylab.label = 'Median Distance',
	yat = 5:8,
	yaxis.lab = c(expression(10^5), expression(10^6), expression(10^7), expression(10^8)),
	ylimits = c(4.5,8.8),
	xlab.label = 'Location Translocations to ecDNA',
	xaxis.lab = rep(c('Within','Outside'),4),
	#xaxis.rot = 45,
	xaxis.cex = 0.9,
	xlab.cex = 1.5,
	ylab.cex = 1.5,
	yaxis.cex = 1.2,
	add.rectangle = TRUE,
	xleft.rectangle = c(0,2.5,4.5,6.5),
	ybottom.rectangle = 0,
	xright.rectangle = c(2.5,4.5,6.5,10),
	ytop.rectangle = 10,
	col.rectangle = c('darkorange2','dodgerblue4','darkorange2', 'dodgerblue4'),
	alpha.rectangle = 0.25,
	add.text = TRUE,
	text.labels = c(
		'ER+ High','ER+ Typical', 'HER2+','IC10',
		expression('FC=0.50'), expression('FDR=1.83x10'^'-8'),
		expression('FC=0.76'), expression('FDR=0.30'),
		expression('FC=0.52'), expression('FDR=8.30x10'^'-4'),
		expression('FC=0.83'), expression('FDR=0.86')
		),
	text.x = c(1.5,3.5, 5.5,7.5,1.5, 1.5, 3.5, 3.5, 5.5, 5.5, 7.5, 7.5),
	text.y = c(rep(8.6,4), rep(c(4.85,4.7),4)),
	text.anchor = 'centre',
	text.col = 'black',
	text.cex = 1,
	text.fontface = 'bold',
	#main = 'Primary',
	#main.cex = 2,
	filename = 'extended_data7k.pdf',
	resolution = 300 
	)


res <- list()
for (i in unique(met_ecdna_plot$group)[!is.na(unique(met_ecdna_plot$group))]) {
	tmp <- met_ecdna_plot[which(met_ecdna_plot$group == i),]
	stats <- wilcox.test(
		tmp[tmp$type == 'median_ecdna_rloops','distance'],
		tmp[tmp$type == 'median_tra_rloops','distance']
		)
	es <- median(tmp[tmp$type == 'median_ecdna_rloops','distance'], na.rm = TRUE)/median(tmp[tmp$type == 'median_tra_rloops','distance'], na.rm = TRUE)
	res[[i]] <- data.frame(
		group = i,
		es = es,
		p = stats$p.value
		)
}
res <- do.call(rbind, res)
res$fdr <- p.adjust(res$p, method = 'fdr')