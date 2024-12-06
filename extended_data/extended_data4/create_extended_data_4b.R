### EXTENDED DATA 4B: COMPLEX AMPLICON DISTRIBUTION ACROSS STAGES HER2+ ###########################
# create barplot of ecDNA in HER2 stratified by ER status in primary and metastatic

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
megatable$any_complex <- NA
megatable[which(megatable$AA_ecdna > 0),'any_complex'] <- 'cyclic'
megatable[which(megatable$AA_ecdna == 0 & megatable$AA_complex > 0),'any_complex'] <- 'complex'
megatable[which(megatable$AA_ecdna == 0 & megatable$AA_complex == 0),'any_complex'] <- 'none'
megatable[which(megatable$group == 'HER2+' & megatable$ER == 1),'group'] <- 'HER2+/ER+'
megatable[which(megatable$group == 'HER2+' & megatable$ER == 0),'group'] <- 'HER2+/ER-'
plot_data <- as.data.frame(table(megatable[,c('any_complex','group')]))
plot_data$prop <- NA
for (i in unique(plot_data$group)) {
	plot_data[plot_data$group == i,'prop'] <- plot_data[plot_data$group == i,'Freq']/sum(plot_data[plot_data$group == i,'Freq'])

}
plot_data$stage <- 'primary'

# read in megatable
metmegatable <- read.delim(
	file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
	as.is = TRUE
	)
metmegatable[which(metmegatable$AA_ecdna > 0),'any_complex'] <- 'cyclic'
metmegatable[which(metmegatable$AA_ecdna == 0 & metmegatable$AA_complex > 0),'any_complex'] <- 'complex'
metmegatable[which(metmegatable$AA_ecdna == 0 & metmegatable$AA_complex == 0),'any_complex'] <- 'none'
metmegatable[which(metmegatable$group == 'HER2+' & metmegatable$ER == 1),'group'] <- 'HER2+/ER+'
metmegatable[which(metmegatable$group == 'HER2+' & metmegatable$ER == 0),'group'] <- 'HER2+/ER-'
met_plot_data <- as.data.frame(table(metmegatable[,c('any_complex','group')]))
met_plot_data$prop <- NA
for (i in unique(met_plot_data$group)) {
	met_plot_data[met_plot_data$group == i,'prop'] <- met_plot_data[met_plot_data$group == i,'Freq']/sum(met_plot_data[met_plot_data$group == i,'Freq'])

}
met_plot_data$stage <- 'metastatic'

plot_data <- rbind(plot_data[plot_data$group %in% c('HER2+/ER+', 'HER2+/ER-'),], met_plot_data[met_plot_data$group %in% c('HER2+/ER+', 'HER2+/ER-'),])
plot_data$stage_code <- NA
plot_data[plot_data$stage == 'primary','stage_code'] <- 'a'
plot_data[plot_data$stage == 'metastatic','stage_code'] <- 'b'

plot_data$amp_code <- NA
plot_data[plot_data$any_complex == 'cyclic','amp_code'] <- 'd'
plot_data[plot_data$any_complex == 'complex','amp_code'] <- 'e'
plot_data[plot_data$any_complex == 'none','amp_code'] <- 'f'


plot_data$index <- paste(plot_data$group, plot_data$stage_code, sep ='_')
plot_data$gindex <- paste(plot_data$stage_code, plot_data$amp_code, sep ='_')

# create barplot
p1 <- create.barplot(
	prop ~ index,
	groups = gindex,
	col = c('#3d85c6ff','#9ebff1ff','white','#b55e93ff','#fe85cfff','white'),
	stack = TRUE,
	xaxis.lab = rep('', 100),
	data = plot_data,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	height = 8,
	width = 8,
	abline.v = 2.5,
	abline.lty = 2,
	ylab.label = 'Proportion of Samples',
	xlab.label = '',
	resolution = 300
	)
p2 <- create.barplot(
	Freq ~ index,
	groups = gindex,
	col = c('#3d85c6ff','#9ebff1ff','white','#b55e93ff','#fe85cfff','white'),
	data = plot_data,
	stack = TRUE,
	xat = seq(1.5, nrow(plot_data), 2),
	xaxis.lab = unique(plot_data$group),
	ylimits = c(0,300),
	yat = seq(0,300,50),
	height = 8,
	width = 8,
	abline.v = 2.5,
	abline.lty = 2,
	ylab.label = 'Number of Samples',
	xlab.label = 'Subtype',
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('#3d85c6ff','#9ebff1ff')
                             ),
                         text = list(
                             lab = c('Cyclic','Complex non-cyclic')
                             ),
                         title = 'Primary',
                         padding.text = 5,
                         cex = 1.2
                         )
                     ),
                     # Positioning legend on plot
                     x = 0.01,
                     y = 0.99,
                     corner = c(0,1)
                 ),
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('#b55e93ff','#fe85cfff')
                             ),
                         text = list(
                             lab = c('Cyclic','Complex non-cyclic')
                             ),
                         title = 'Metastatic',
                         padding.text = 5,
                         cex = 1.2
                         )
                     ),
                     # Positioning legend on plot
                     x = 0.01,
                     y = 0.6,
                     corner = c(0,1)
                 )
             ),
	resolution = 300
	)
create.multipanelplot(
	list(p1,p2),
	filename = 'extended_data4b.pdf',
	height = 9,
	width = 7,
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData4b_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)