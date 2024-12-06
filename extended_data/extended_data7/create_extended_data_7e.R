### CREATE EXTENDED DATA 7E #######################################################################
# plot ESR1 levels in DCIS vs IBC in ER+ typical and ER+ high

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(DESeq2)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in source data 
paired <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7e_sourcetable.txt'),
	as.is = TRUE
	)
# create plot data
plot_data <- gather(paired[,-1], value = 'esr1', key = 'type', -group)
plot_data$index <- paste(plot_data$group, plot_data$type, sep = '_')

# calculate stats
erhigh_stats <- wilcox.test(paired[paired$group == 'ER+ High','DCIS_Epithelial'], paired[paired$group == 'ER+ High','IBC_Epithelial'], paired = TRUE, conf.int =TRUE)
ertypical_stats <- wilcox.test(paired[paired$group == 'ER+ Typical','DCIS_Epithelial'], paired[paired$group == 'ER+ Typical','IBC_Epithelial'], paired = TRUE, conf.int =TRUE)

create.boxplot(
	esr1 ~ index,
	data = plot_data,
	add.stripplot = TRUE,
	ylimits = c(-0.5,0.5),
	yat = seq(-0.4,0.4,0.2),
	add.rectangle = TRUE,
	xleft.rectangle = 0,
	ybottom.rectangle = -0.5,
	xright.rectangle = 2.5,
	ytop.rectangle = 12,
	col.rectangle = 'mediumpurple2',
	alpha.rectangle = 0.25,
	add.text = TRUE,
	text.x = c(1.5,3.5, 1.5, 3.5, 1.5,3.5),
	text.y = c(0.47, 0.47, 0.43, 0.43, 0.39, 0.39),
	#text.y = c(11.7,11.7,11.1,11.1,10.7,10.7),
	text.labels = c('ER+ High','ER+ Typical', paste('ES=', round(erhigh_stats$estimate[[1]], digits = 2)),
		paste('ES=', round(ertypical_stats$estimate[[1]], digits = 2)), paste('P=', round(erhigh_stats$p.value, digits = 2)),
		paste('P=', round(ertypical_stats$p.value, digits = 2))),
	text.fontface = c('bold','bold','plain','plain','plain','plain'),
	xaxis.lab = c('DCIS','Invasive','DCIS','Invasive'),
	xlab.label = 'Stage of Progression',
	ylab.label = 'ER Early Signaling',
	filename = 'extended_data7e.pdf',
	resolution = 300
	)