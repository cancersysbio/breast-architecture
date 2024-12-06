### CREATE SUPPLEMENTARY FIGURE 2CD ################################################################
# create supplementary figure 2cd 
# create heatmap of SV and CNA signatures correlations in primary cohort
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
# read in megatable
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)

### SUPPLEMENTARY FIGURE 2D #######################################################################
# find correlations 
sigs <- megatable[,grep('^CN|^SV32', colnames(megatable))]
i <- apply(sigs, 1, function(x) any(is.na(x)))
plot_data <- cor(sigs[!i,], method = 'spearman')
# create dotmap
plot_data <- plot_data[grep('^CN', rownames(plot_data)),grep('^SV', colnames(plot_data))]
colnames(plot_data) <- gsub('32', '', colnames(plot_data))

cov_cols <- c(
	rep('darkseagreen',3),
	rep('darkorange',5),
	rep('deepskyblue',4),
	rep('deeppink',4),
	'darkolivegreen',
	rep('gold',4)
	)
covariate <- list(
         rect = list(
             col = 'transparent',
             fill = cov_cols,
             lwd = 1.5
             )
         )
cov.legend <- list(
         legend = list(
             colours = rev(c('darkseagreen', 'darkorange','deepskyblue','deeppink','darkolivegreen','gold')),
             labels = rev(c('Ploidy','Chromothripsis','fLOH','cLOH','HRD','Unknown')),
             title = 'CNA Signatures'
             )
         )

# create heatmap
create.heatmap(
	t(plot_data[-c(22:24),]),
	colour.scheme = c('dodgerblue3','white','firebrick3'),
	xaxis.lab = NA,
	yaxis.lab = NA,
	cluster.dimensions = 'columns',
	xaxis.cex = 1,
	yaxis.cex = 1,
	at = seq(-0.5,0.5,0.01),
	covariates = covariate,
	covariate.legends = cov.legend,
	colourkey.cex = 1.5,
	filename = 'supplementary_figure2d.pdf',
	resolution = 300
	)

### SUPPLEMENTARY FIGURE 2C #######################################################################
# set file name
nk_file <- file.path(main_repo_path, 'data', 'NikZainal_supplementary_table21_rearrangement_sig_composition.csv')
if (!file.exists(nk_file)) {
	stop("Please download supplementary table 21 from Nik-Zainal et al. and save in data directory as NikZainal_supplementary_table21_rearrangement_sig_composition.csv ...")
}
# read in nk composition 
nk <- read.csv(
	nk_file,
	skip = 1,
	as.is = TRUE
	)
for (i in 3:8) {
	nk[,i] <- as.numeric(gsub('%','', nk[,i]))
}
sig6 <- read.delim(
	file.path(main_repo_path, 'data', 'sv_signature_composition_primary.txt'),
	as.is = TRUE
	)

plot_data2 <- cor(
	cbind(nk[,-c(1,2)], sig6[,-1]),
	method = 'spearman'
	)
rownames(plot_data2) <- gsub('Rearrangement.Signature.','RS',rownames(plot_data2))
plot_data2 <- plot_data2[grep('RS', rownames(plot_data2)),grep('SV', colnames(plot_data2))]
colnames(plot_data2) <- gsub('32','', colnames(plot_data2))
plot_data2 <- plot_data2[c('RS6','RS4','RS2','RS1','RS3','RS5'),c('SVA','SVD','SVB','SVF','SVC','SVE')]

# create heatmap
create.heatmap(
	t(plot_data2),
	colour.scheme = c('dodgerblue3','white','firebrick3'),
	xaxis.lab = NA,
	yaxis.lab = NA,
	clustering.method = 'none',
	xaxis.cex = 1,
	yaxis.cex = 1,
	at = seq(-0.5,1,0.01),
	right.padding = 3,
	ylab.label = 'Nik-Zainal RS',
	xlab.label = 'De novo RS',
	colourkey.cex = 1.5,
	filename = 'supplementary_figure2c.pdf',
	resolution = 300
	)