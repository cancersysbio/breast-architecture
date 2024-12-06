### CREATE FIGURE 2C ##############################################################################
# plot archetype correlations
# also generates supplementary figure 3g

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### TEST FEATURES AGAINST ARCHETYPES ##############################################################
test_features_against_archtypes <- function(arch, features) {
	# test correlation with archetypes
	res <- do.call(rbind, sapply(
		features,
		function(x) {
			print(x)
			arch1 <- cor.test(arch$Arc1, arch[,x], method = 'spearman')
			arch2 <- cor.test(arch$Arc2, arch[,x], method = 'spearman')
			arch3 <- cor.test(arch$Arc3, arch[,x], method = 'spearman')
			data.frame(
				rho = c(arch1$estimate, arch2$estimate, arch3$estimate),
				p = c(arch1$p.value, arch2$p.value, arch3$p.value),
				archetype = c('Arc1','Arc2','Arc3'),
				feature = x
				)
			},
		simplify = FALSE
		))
	res$fdr <- p.adjust(res$p, method = 'fdr')
	return(res)
}

### SET FEATURE ORDER #############################################################################
set_feature_order <- function(plot_data) {
	tmp1 <- plot_data[plot_data$archetype == 'Arc1',]
	feat_order1 <- tmp1[order(tmp1$rho),'feature']
	tmp3 <- plot_data[plot_data$archetype == 'Arc3',]
	feat_order3 <- tmp3[order(tmp3$rho),'feature']
	tmp2 <- plot_data[plot_data$archetype == 'Arc2',]
	feat_order2 <- tmp2[order(tmp2$rho),'feature']
	# minor edits to feature order
	feat_order <- c(
		grep(
			'SBS',
			unique(c(as.character(feat_order1), as.character(feat_order3),as.character(feat_order2))),
			value = TRUE
			),
		grep(
			'SBS',
			unique(c(as.character(feat_order1), as.character(feat_order3),as.character(feat_order2))),
			invert = TRUE,
			value = TRUE
			)
		)
	# swapping dm and inv to match categorization 
	feat_order <- feat_order[c(1:10,17,12,11,13:16,18:27)]
	return(feat_order)
}
 
### MAIN ##########################################################################################
# read in source data 
arch <- read.delim(
	file.path(main_repo_path, 'data', 'Figure2c_sourcetable_primary.txt'),
	as.is = TRUE
	)
# reformat genome doubling
arch$genome_doubled <- (arch$genome_doubled == 'True')*1
# add dm and cmplxdm together
arch$dm <- arch$dm+arch$cpxdm

# test correlation with archetypes
features <- colnames(arch)
features <- features[!features %in% c('sample','Arc1','Arc2','Arc3','group','cpxdm')]
res <- test_features_against_archtypes(arch = arch, features = features)

### METASTATIC ####################################################################################
# read in source data 
metarch <- read.delim(
	file.path(main_repo_path, 'data', 'Figure2c_sourcetable_metastatic.txt'),
	as.is = TRUE
	)
# reformat genome doubled
metarch$genome_doubled <- (metarch$genome_doubled == 'True')*1
# add dm and cmplxdm together
metarch$dm <- metarch$dm+metarch$cpxdm
# test features
metres <- test_features_against_archtypes(arch = metarch, features = features)

### PLOT ################################################################################
# find sig features and create plot data
sig_features <- unique(res[res$fdr < 0.05,'feature'])
primplot_data <- res[res$feature %in% sig_features,]
metplot_data <- metres[metres$feature %in% sig_features,]
metplot_data$stage <- 'metastatic'
primplot_data$stage <- 'primary'
plot_data <- rbind(primplot_data, metplot_data)

plot_data$code <- NA
plot_data[plot_data$archetype == 'Arc2' & plot_data$stage == 'metastatic','code'] <- '2. ER+ Typical -enriched'
plot_data[plot_data$archetype == 'Arc3' & plot_data$stage == 'metastatic','code'] <- '1. HER2+/ER+ High -enriched'
plot_data[plot_data$archetype == 'Arc1' & plot_data$stage == 'metastatic','code'] <- '3. TNBC -enriched'
plot_data[plot_data$archetype == 'Arc1' & plot_data$stage == 'primary','code'] <- '3. TNBC -enriched'
plot_data[plot_data$archetype == 'Arc2' & plot_data$stage == 'primary','code'] <- '2. ER+ Typical -enriched'
plot_data[plot_data$archetype == 'Arc3' & plot_data$stage == 'primary','code'] <- '1. HER2+/ER+ High -enriched'

# find top 10 most correlated features per archetype
feat_imp <- list()
for (i in c('Arc1','Arc2','Arc3')) {
	tmp <- plot_data[plot_data$stage == 'primary' & plot_data$archetype == i,]
	feat_imp[[i]] <- tmp[order(-abs(tmp$rho)),'feature'][1:10]
}
feat_imp <- unique(unlist(feat_imp))
# don't have both ecDNA and dm in main figure 
feat_imp <- feat_imp[!feat_imp %in% c('chromoplexy','dm')]
feat_imp <- c(as.character(feat_imp), 'AA_complex')

# set feature order
feat_order <- set_feature_order(primplot_data)
main <- feat_order[feat_order %in% feat_imp]
# edit feature name
feat_xlab <- main
feat_xlab <- gsub('fraction_cna','fga', feat_xlab)
feat_xlab <- gsub('genome_doubled','wgd', feat_xlab)
feat_xlab <- gsub('AA_ecdna','ecdna', feat_xlab)
feat_xlab <- gsub('AA_complex','complex amp', feat_xlab)
feat_xlab <- gsub('hrd_loh','HRD LOH', feat_xlab)
feat_xlab <- gsub('hrd_snv','HRD SBS', feat_xlab)
feat_xlab <- gsub('apobec','APOBEC SBS', feat_xlab)

# keep features in main plot
main_plot_data <- plot_data[plot_data$feature %in% main,]

main_plot_data$index <- NA
for (i in 1:length(main)) {
	feattmp <- main[i]
	main_plot_data[main_plot_data$feature == feattmp,'index'] <- i
}

cov.legend <- list(
         legend = list(
             colours = c('#9ebff1ff','#fe85cfff'),
             labels = c('Primary', 'Metastatic'),
             border = 'transparent'
             )
         )
cov.legend.grob <- legend.grob(legends = cov.legend)

create.scatterplot(
	rho ~ index | code,
	data = main_plot_data,
	groups = main_plot_data$stage,
	type = c('h','p'),
	abline.h = 0,
	abline.v = c(4.5,8.5),
	abline.lty = 2,
	xaxis.rot = 45,
	xaxis.cex = 1.25,
	ylimits = c(-1,1),
	xlimits = c(0,16),
	#xaxis.rot = 45,
	layout = c(1,3),
	yat = seq(-1,1,0.5),
	yaxis.cex = 1,
	col = c('#fe85cfff','#9ebff1ff'),
	ylab.label = 'Spearman\'s Rho',
	xlab.label = 'Genomic Features',
	add.text = TRUE,
	text.labels = c('Simple SVs','Complex\nAmplifications','SNV\nSignatures'),
	text.x = c(12.5,6.5,2.5),
	text.y = c(0.87,0.82,0.82),
	#layout = c(3,1),
	legend = list(inside = list(fun = cov.legend.grob, x = 0, y = 0.71)),
	text.cex = 1,
	xaxis.lab = feat_xlab,
	xat = 1:length(main),
	width = 13,
	height = 8,
	filename = 'figure2c.pdf',
	resolution = 300
	)

### CREATE SUPPLEMENTARY FIGURE 3G ################################################################
# set feature order
feat_order <- set_feature_order(primplot_data)
supp <- feat_order[!feat_order %in% feat_imp]
# edit feature name
feat_xlab <- supp
feat_xlab <- gsub('fraction_cna','fga', feat_xlab)
feat_xlab <- gsub('genome_doubled','wgd', feat_xlab)
feat_xlab <- gsub('AA_ecdna','ecdna', feat_xlab)
feat_xlab <- gsub('AA_complex','complex amp', feat_xlab)
feat_xlab <- gsub('hrd_loh','HRD LOH', feat_xlab)
feat_xlab <- gsub('hrd_snv','HRD SBS', feat_xlab)
feat_xlab <- gsub('dm','ecDNA (JaBbA)', feat_xlab)

# keep features in main plot
supp_plot_data <- plot_data[plot_data$feature %in% supp,]

supp_plot_data$index <- NA
for (i in 1:length(supp)) {
	feattmp <- supp[i]
	supp_plot_data[supp_plot_data$feature == feattmp,'index'] <- i
}

# plot figure
create.scatterplot(
	rho ~ index | code,
	data = supp_plot_data,
	#horizontal = TRUE,
	groups = supp_plot_data$stage,
	type = c('h','p'),
	abline.h = 0,
	abline.v = 5.5,
	abline.lty = 2,
	xaxis.rot = 45,
	xaxis.cex = 1,
	ylimits = c(-1,1),
	xlimits = c(0,13),
	#xaxis.rot = 45,
	yat = seq(-1,1,0.5),
	yaxis.cex = 1.2,
	col = c('#fe85cfff','#9ebff1ff'),
	ylab.label = 'Spearman\'s Rho',
	xlab.label = 'Genomic Features',
	add.text = TRUE,
	text.labels = c('Structural\nVariants','SNV\nSignatures'),
	text.x = c(9.5,3),
	text.y = c(0.82,0.82),
	layout = c(3,1),
	legend = list(inside = list(fun = cov.legend.grob, x = 0, y = 0.12)),
	text.cex = 1,
	xaxis.lab = feat_xlab,
	xat = 1:length(supp),
	width = 17,
	height = 5,
	filename = 'supplementary_figure3g.png',
	resolution = 300
	)

### WRITE TO FILE #################################################################################
# write to file
write.table(
	plot_data,
	file = 'Figure2c_sourcetable_statistics.txt',
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
write.table(
	supp_plot_data,
	file = 'SupplementaryFigure3g_sourcetable_statistics.txt',
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

