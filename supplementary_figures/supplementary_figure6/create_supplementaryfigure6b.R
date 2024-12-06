### CREATE EXTENDED DATA SF6B #######################################################################
# create boxplot of translocation proportions between cyclic and non-cyclic amplifications

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}
### CALCULATE STATS ###############################################################################
calculate_stats <- function(plot_data) {
	res <- list()
	for (i in c('IC10','ER+ High')) {
		tmp <- plot_data[which(plot_data$group == i),]
		if (i == 'ER+ High') {
			stats <- wilcox.test(
				tmp[tmp$type == 'cyclic','prop_TRA'],
				tmp[tmp$type == 'complex','prop_TRA'],
				alternative = 'greater'
				)
		} else {
			stats <- wilcox.test(
				tmp[tmp$type == 'cyclic','prop_TRA'],
				tmp[tmp$type == 'complex','prop_TRA'],
				alternative = 'less'
				)
			}
		es <- median(tmp[tmp$type == 'cyclic','prop_TRA'], na.rm = TRUE)/median(tmp[tmp$type == 'complex','prop_TRA'], na.rm = TRUE)
		res[[i]] <- data.frame(
			group = i,
			es = es,
			p = stats$p.value
			)
	}
	res <- do.call(rbind, res)
	res$fdr <- p.adjust(res$p, method = 'fdr')
	return(res)
}

### MAIN ##########################################################################################
# source data generated using create_translocation_overlaps.R found in analysis
# read in primary plot data
primary_plot_data <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7b_sourcetable_primary.txt'),
	as.is = TRUE
	)

# read in primary plot data
met_plot_data <- read.delim(
	file.path(main_repo_path, 'data', 'SupplementaryFigure6b_sourcetable_metastatic.txt'),
	as.is = TRUE
	)

# reformat for plotting
pplot_data_subset <- primary_plot_data[primary_plot_data$group %in% c('ER+ High','IC10'),]
pplot_data_subset$stage <- 'aprim'
mplot_data_subset <- met_plot_data[met_plot_data$group %in% c('ER+ High','IC10'),]
mplot_data_subset$stage <- 'bmet'
plot_data_subset <- rbind(
	pplot_data_subset,
	mplot_data_subset
	)
plot_data_subset$index <- paste0(plot_data_subset$stage, plot_data_subset$group, plot_data_subset$type_code)

create.boxplot(
	prop_TRA ~ index,
	data = plot_data_subset,
	add.stripplot = TRUE,
	xlab.label = 'Amplicon Type',
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	ylab.label = 'Translocation Proportion',
	xlab.cex = 1.5,
	ylab.cex = 1.5,
	xaxis.lab = rep(c('Cyclic','Non-\ncyclic'),5),
	xaxis.cex = 1,
	add.rectangle = TRUE,
	xleft.rectangle = c(0,2.5,4.5,6.5),
	ybottom.rectangle = 0,
	xright.rectangle = c(2.5,4.5,6.5,9),
	ytop.rectangle = 10,
	col.rectangle = c('darkorange2','dodgerblue4','darkorange2','dodgerblue4'),
	alpha.rectangle = 0.25,
	add.text = TRUE,
	text.labels = c(
		'Primary','Metastatic',
		'ER+ High','IC10','ER+ High','IC10',
		expression('FC=1.22'), expression('P=0.12'),
		expression('FC=0.75'), expression('P=0.06'),
		expression('FC=1.30'), expression('P=0.15'),
		expression('FC=0.81'), expression('P=0.22')
		),
	abline.v = 4.5,
	abline.lty = 2,
	text.x = c(2.5,6.5,1.5,3.5, 5.5,7.5,1.5, 1.5, 3.5, 3.5, 5.5, 5.5, 7.5, 7.5),
	text.y = c(rep(0.975, 2), rep(0.935,4),rep(c(0.90,0.865),4)),
	text.anchor = 'centre',
	text.col = 'black',
	text.cex = 1,
	text.fontface = 'bold',
	filename = 'supplementary_figure6b.pdf',
	resolution = 300 
	)

### REPORT STATS ##################################################################################
# calculate stats
met_res <- calculate_stats(plot_data = met_plot_data)
prim_res <- calculate_stats(plot_data = primary_plot_data)
