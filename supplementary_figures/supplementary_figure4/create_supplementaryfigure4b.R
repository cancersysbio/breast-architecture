### CREATE SUPPLEMENTARY FIGURE 4B ################################################################
# create supplementary figure 4b
# correlation heatmap comparing primary to metastatic signatures
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in primary signatures 
primary <- read.delim(
	file.path(main_repo_path, 'data', 'sv_signature_composition_primary.txt'),
	as.is = TRUE
	)
colnames(primary) <- gsub('SV32','Primary_SV32', colnames(primary))
# read in metastatic signatures 
metastatic <- read.delim(
	file.path(main_repo_path, 'data', 'sv_signature_composition_metastatic.txt'),
	as.is = TRUE
	)
colnames(metastatic) <- gsub('SV32','Metastatic_SV32', colnames(metastatic))
# merge 
plot_data <- merge(primary, metastatic, by = 'MutationType')
plot_data <- cor(plot_data[,-1], method = 'spearman')

col <- rep(NA, ncol(plot_data))
col[grep('Primary', colnames(plot_data))] <- '#9ebff1ff'
col[grep('Metastatic', colnames(plot_data))] <- '#fe85cfff'

covariate <- list(
	rect = list(
		col = 'white',
		fill = col,
		lwd = 1.5
		)
	)

covariate.legend <- list(
	legend = list(
		colours = c('#9ebff1ff','#fe85cfff'),
		labels = c('Primary','Metastatic'),
		border = 'white',
		title = 'Disease'
		)
	)

plot_data <- plot_data[grep('Metastatic', rownames(plot_data)), grep('Primary', colnames(plot_data))]
rownames(plot_data) <- gsub('32', '', rownames(plot_data))
colnames(plot_data) <- gsub('32', '', colnames(plot_data))


# create heatmap
create.heatmap(
	plot_data,
	colour.scheme = c('dodgerblue3','white','firebrick3'),
	#ylab.label = 'Primary',
	#xlab.label = 'Metastasis',
	xaxis.lab = gsub('Primary_|Metastatic_', '', colnames(plot_data)),
	yaxis.lab = gsub('Primary_|Metastatic_', '', colnames(plot_data)),
	xaxis.cex = 1,
	yaxis.cex = 1,
	covariates.top = covariate,
	covariates = covariate,
	covariate.legends = covariate.legend,
	covariates.top.grid.border = list(col = 'white', lwd = 1.5),
	colourkey.cex = 1.5,
	filename = 'supplementary_figure4b.pdf',
	resolution = 300,
	height = 5.5,
	width = 6
	)
