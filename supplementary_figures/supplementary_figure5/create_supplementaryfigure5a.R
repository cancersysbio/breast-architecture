### SUPPLEMENTARY FIGURE 5A: COMPLEX AMPLICON DISTRIBUTION #########################################
# create stacked barplot comparing distribution of amplicons to 
# Kim et al. 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in megatable
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
# reformat
plot_data <- data.frame(
	type = c('cyclic','complex','linear','none'),
	count = c(
		sum(megatable$AA_ecdna > 0, na.rm = TRUE),
		sum(megatable$AA_ecdna == 0 & megatable$AA_complex > 0, na.rm = TRUE),
		sum(megatable$AA_ecdna == 0 & megatable$AA_complex == 0 & megatable$AA_amplicon > 0, na.rm = TRUE),
		sum(megatable$AA_ecdna == 0 & megatable$AA_complex == 0 & megatable$AA_amplicon == 0, na.rm = TRUE)
		))
plot_data$prop <- plot_data$count/sum(plot_data$count)

# read in kim et al 
# this is supplementary table 3 from Kim et al.
kim_etal_file <- file.path(main_repo_path, 'data', 'kim_etal_supplementary_table3.csv')
if (!file.exists(kim_etal_file)) {
	stop("Please download Supplementary Table 3 from Kim et al. and place it in data directory named kim_etal_supplementary_table3.csv ...")
}
kim <- read.csv(
	kim_etal_file,
	as.is = TRUE
	)
kim <- kim[kim$lineage == 'Breast',]
kim[kim$sample_classification == 'BFB','sample_classification'] <- 'cyclic'
kim[kim$sample_classification == 'Circular','sample_classification'] <- 'cyclic'
kim[kim$sample_classification == 'Heavily-rearranged','sample_classification'] <- 'complex'
kim[kim$sample_classification == 'Linear','sample_classification'] <- 'linear'
kim[kim$sample_classification == 'No-fSCNA','sample_classification'] <- 'none'

kim_plot_data <- as.data.frame(table(kim$sample_classification))
colnames(kim_plot_data) <- c('type','count')
kim_plot_data$prop <- kim_plot_data$count/sum(kim_plot_data$count)
# set cohort
plot_data$cohort <- 'ours'
kim_plot_data$cohort <- 'kim'
# create plot data
plot_data <- rbind(plot_data, kim_plot_data)

# create barplot
create.barplot(
	prop ~ cohort,
	groups = plot_data$type,
	data = plot_data,
	stack = TRUE,
	ylab.label = 'Proportion of Samples',
	xlab.label = 'Cohort',
	xaxis.lab = c('Kim et al.', 'Current'),
	legend = list(
             right = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 2,
                             fill = default.colours(5)[-4]
                             ),
                         text = list(
                             lab = c('Cyclic','Complex Non-cyclic','Linear','None')
                             ),
                         padding.text = 3,
                         cex = 1
                         )
                     )
                 )
             ),
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	col = default.colours(5)[-4],
	filename = 'supplementary_figure5a.pdf',
	resolution = 300
	)

### WRITE TO TABLE ################################################################################
# write to file 
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'SupplementaryFigure5a_sourcetable.txt'),
	sep = '\t', 
	row.names = FALSE
	)
