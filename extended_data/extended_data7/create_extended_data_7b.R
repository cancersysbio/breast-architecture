### GENERATE EXTENDED DATA 7B #####################################################################
# create boxplot of translocations overlapping with ecDNA

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)

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

### OVERLAPS ##########################################################################################
# read in primary plot data
primary_plot_data <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7b_sourcetable_primary.txt'),
	as.is = TRUE
	)
primary_plot_data <- primary_plot_data[which(primary_plot_data$group %in% c('HER2+/ERb','HER2+/ERa','ER+ High','ER+ Typical','IC10')),]
primary_plot_data <- primary_plot_data[primary_plot_data$number_SV >= 3,]
primary_plot_data$index <- paste0(primary_plot_data$group, primary_plot_data$type_code)

res <- list()
for (i in unique(primary_plot_data$group)) {
	tmp <- primary_plot_data[which(primary_plot_data$group == i),]
	stats <- wilcox.test(
		tmp[tmp$type == 'cyclic','number_TRA'],
		tmp[tmp$type == 'complex','number_TRA']
		)
	es <- median(tmp[tmp$type == 'cyclic','number_TRA'], na.rm = TRUE)/median(tmp[tmp$type == 'complex','number_TRA'], na.rm = TRUE)
	res[[i]] <- data.frame(
		group = i,
		es = es,
		p = stats$p.value
		)
}
res <- do.call(rbind, res)
res$fdr <- p.adjust(res$p, method = 'fdr')

create.boxplot(
	lognumber_TRA ~ index,
	data = primary_plot_data,
	add.stripplot = TRUE,
	xlab.label = 'Amplicon Type',
	ylimits = c(0,3),
	yat = 0:3,
	yaxis.lab = c(expression(0), expression(10^1), expression(10^2), expression(10^3)),
	ylab.label = 'Number of Translocations',
	xaxis.lab = rep(c('Cyclic','Non-\nCyclic'),5),
	xaxis.cex = 1.2,
	yaxis.cex = 1.2,
	add.rectangle = TRUE,
	xleft.rectangle = c(0,2.5,4.5,6.5,8.5),
	ybottom.rectangle = 0,
	xright.rectangle = c(2.5,4.5,6.5,8.5,11),
	ytop.rectangle = 10,
	col.rectangle = c('darkorange2', 'dodgerblue4','darkorange2','dodgerblue4','darkorange2'),
	alpha.rectangle = 0.25,
	add.text = TRUE,
	text.labels = c(
		'ER+ High','ER+ Typical', 'HER2+/ER+','HER2+/ER-','IC10',
		expression('FC=3.4'), expression('FDR=1.1x10'^'-3'),
		expression('FC=8.5'), expression('FDR=3.7x10'^'-6'),
		expression('FC=3.0'), expression('FDR=7.4x10'^'-3'),
		expression('FC=2.7'), expression('FDR=0.46'),
		expression('FC=1.0'), expression('FDR=0.81')
		),
	text.x = c(1.5,3.5, 5.5,7.5, 9.5, 1.5, 1.5, 3.5, 3.5, 5.5, 5.5, 7.5, 7.5,9.5,9.5),
	text.y = c(rep(2.9,5),rep(c(2.78,2.68),5)),
	text.anchor = 'centre',
	text.col = 'black',
	text.cex = 1.2,
	width = 8.7,
	height = 7,
	text.fontface = 'bold',
	filename = 'extended_data7b.pdf',
	resolution = 300 
	)
