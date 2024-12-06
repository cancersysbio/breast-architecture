### CREATE EXTENDED DATA 4L #####################################################################
# create boxplot of archetypes comparing ER+ typical with ecDNA to high and typical tumors

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
	file.path(main_repo_path, 'data','primary_megatable.txt'),
	as.is = TRUE
	)
megatable <- megatable[which(megatable$group != ''),]
# find er typical with ecDNA 
erlow_ecdna <- megatable[which(megatable$AA_ecdna > 0 & megatable$group == 'ER+ Typical'),'Sample']

# read archetypes 
arch <- read.delim(
	file.path(main_repo_path, 'data', 'Figure2btop_sourcetable.txt',
	as.is = TRUE
	)
arch <- merge(megatable[,c('Sample','group')], arch, by = 'Sample')
arch <- arch[arch$group %in% c('ER+ High','ER+ Typical'),]
arch[arch$Sample %in% erlow_ecdna,'group'] <- 'ER+ Typical\necDNA'

plot_data <- gather(arch[,1:5], key = 'archetype', value = 'value', -Sample, -group)
plot_data[plot_data$archetype == 'Arc1','archetype'] <- 'TNBC-enriched'
plot_data[plot_data$archetype == 'Arc2','archetype'] <- 'HER2+/ER+ High-enriched'
plot_data[plot_data$archetype == 'Arc3','archetype'] <- 'ER+ Typical-enriched'

plot_data$value <- as.numeric(as.character(plot_data$value))

create.boxplot(
	value ~ group | archetype,
	data = plot_data,
	add.stripplot = TRUE,
	xlab.label = 'Subtype',
	ylab.label = 'Proportion',
	xaxis.lab = c('ER+\nHigh','ER+\nTypical','ER+\nTypical\necDNA'),
	xaxis.cex = 0.8,
	filename = 'extended_data4l.pdf',
	resolution = 300
	)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  plot_data,
  file = file.path(main_repo_path, 'data', 'ExtendedData4l_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )
