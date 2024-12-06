### CREATE EXTENDED DATA FIGURE 2B ################################################################
# create boxplot of fraction of genome altered across subgroups

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in plot data 
plot_data <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData2b_sourcetable.txt'),
	as.is = TRUE
	)

# create boxplot
create.boxplot(
	fraction_cna ~ cohort | group,
	data = plot_data,
	add.stripplot = TRUE,
	xaxis.rot = 45,
	#ylimits = c(0,1000),
	xaxis.lab = c('DCIS sWGS','Primary sWGS','Primary dWGS','Metastatic dWGS'),
	xlab.label  = 'Disease Stage',
	#ylab.label = 'Number of Segments',
	ylab.label = 'Fraction Genome Altered',
	xaxis.cex = 1,
	yaxis.cex= 1,
	filename = 'extended_data2b.pdf',
	resolution = 300
	)







