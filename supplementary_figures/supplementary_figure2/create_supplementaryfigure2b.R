### CREATE SUPPLEMENTARY FIGURE 2B ##############################################################
# create supplementary figure 2b 
# plot NMF diagnostic measures for choice of number of SV signatures in primary samples 
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### PRIMARY 2B ####################################################################################
# read in diagnostic measures
diagnostic <- read.csv(
	file.path(main_repo_path, 'data', 'sv_signature_diagnostics_primary.csv'),
	as.is = TRUE
	)
colnames(diagnostic) <- gsub('Stability..Avg.Silhouette.','Silhouette', colnames(diagnostic))
colnames(diagnostic) <- gsub('Mean.Cosine.Distance','Cosine', colnames(diagnostic))
diagnostic$Signatures <- gsub('*','',diagnostic$Signatures, fixed = TRUE)
# create plot data
plot_data <- gather(
	diagnostic[diagnostic$Signatures %in% 1:10,c('Signatures','Silhouette','Cosine')],
	key = "key",
	value = "value",
	-Signatures
	)
plot_data$Signatures <- as.numeric(plot_data$Signatures)

# create scatterplot 
p1 <- create.scatterplot(
	value ~ Signatures,
	data = plot_data[plot_data$key == 'Cosine',],
	ylimits = c(0,0.4),
	yaxis.tck = 0,
	xaxis.tck = 0,
	xat = 1:10,
	add.rectangle = TRUE,
	xleft.rectangle = 5.5,
	ybottom.rectangle = 0,
	xright.rectangle = 6.5,
	ytop.rectangle = 0.4,
	col.rectangle = 'grey50',
	alpha.rectangle = 0.5,
	yat = seq(0,0.4,0.1),
	ylab.label = 'Cosine Distance',
	xlab.label = 'Signature',
	resolution = 300
	)
p2 <- create.scatterplot(
	value ~ Signatures,
	data = plot_data[plot_data$key == 'Silhouette',],
	ylimits = c(0.6,1.05),
	yaxis.tck = 0,
	xaxis.tck = 0,
	xat = 1:10,
	add.rectangle = TRUE,
	xleft.rectangle = 5.5,
	ybottom.rectangle = 0.6,
	xright.rectangle = 6.5,
	ytop.rectangle = 1.05,
	col.rectangle = 'grey50',
	alpha.rectangle = 0.5,
	yat = seq(0.6,1,0.1),
	ylab.label = 'Silhouette',
	xlab.label = '',
	resolution = 300
	)
create.multipanelplot(
	list(p2, p1),
	filename = file.path(main_repo_path, 'supplementary_figure2b.pdf'),
	resolution = 300
	)