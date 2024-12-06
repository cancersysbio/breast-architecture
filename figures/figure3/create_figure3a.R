### FIGURE 3A: COMPLEX AMPLICON DISTRIBUTION ACROSS STAGES ########################################
# create boxplot of complex amplicon distribution across subgroups and stages
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in primary megatable
megatable <- read.delim(
	file.path(main_repo_path, 'data','primary_megatable.txt'),
	as.is = TRUE
	)
megatable <- megatable[which(megatable$group != ''),]
megatable$any_complex <- NA
megatable[which(megatable$AA_ecdna > 0),'any_complex'] <- 'cyclic'
megatable[which(megatable$AA_ecdna == 0 & megatable$AA_complex > 0),'any_complex'] <- 'complex'
megatable[which(megatable$AA_ecdna == 0 & megatable$AA_complex == 0),'any_complex'] <- 'none'
megatable$group <- gsub('Low','Typical',megatable$group)
plot_data <- as.data.frame(table(megatable[,c('any_complex','group')]))
plot_data$prop <- NA
for (i in unique(plot_data$group)) {
	plot_data[plot_data$group == i,'prop'] <- plot_data[plot_data$group == i,'Freq']/sum(plot_data[plot_data$group == i,'Freq'])

}
plot_data$stage <- 'primary'

# read in metastatic megatable
metmegatable <- read.delim(
	file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
	as.is = TRUE
	)
metmegatable <- metmegatable[which(megatable$group != ''),]
metmegatable[which(metmegatable$AA_ecdna > 0),'any_complex'] <- 'cyclic'
metmegatable[which(metmegatable$AA_ecdna == 0 & metmegatable$AA_complex > 0),'any_complex'] <- 'complex'
metmegatable[which(metmegatable$AA_ecdna == 0 & metmegatable$AA_complex == 0),'any_complex'] <- 'none'
metmegatable$group <- gsub('Low','Typical',metmegatable$group)
met_plot_data <- as.data.frame(table(metmegatable[,c('any_complex','group')]))
met_plot_data$prop <- NA
for (i in unique(met_plot_data$group)) {
	met_plot_data[met_plot_data$group == i,'prop'] <- met_plot_data[met_plot_data$group == i,'Freq']/sum(met_plot_data[met_plot_data$group == i,'Freq'])

}
met_plot_data$stage <- 'metastatic'

# create plot data with dummy code to ensure correct ordering in plot
plot_data <- rbind(plot_data, met_plot_data)
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
	abline.v = c(2.5,4.5,6.5,8.5),
	abline.lty = 2,
	ylab.label = 'Proportion of Samples',
	xlab.label = '',
	#stack= TRUE,
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
	ylimits = c(0,450),
	yat = seq(0,450,50),
	height = 8,
	width = 8,
	abline.v = c(2.5,4.5,6.5,8.5),
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
                     x = 0.99,
                     y = 0.99,
                     corner = c(1,1)
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
                     x = 0.99,
                     y = 0.6,
                     corner = c(1,1)
                 )
             ),
	resolution = 300
	)
create.multipanelplot(
	list(p1,p2),
	filename = 'figure3a.pdf',
	resolution = 300
	)

# create stats for text
tmp <- plot_data[plot_data$stage == 'primary',]
ct <- data.frame(
		subtype = c(
			sum(tmp[tmp$group %in% c('ER+ High','HER2+') & tmp$any_complex != 'none','Freq']),
			sum(tmp[tmp$group %in% c('ER+ High','HER2+') & tmp$any_complex == 'none','Freq'])
			),
		notsubtype = c(
			sum(tmp[!tmp$group %in% c('ER+ High','HER2+') & tmp$any_complex != 'none','Freq']),
			sum(tmp[!tmp$group %in% c('ER+ High','HER2+') & tmp$any_complex == 'none','Freq'])
			)
		)
print("Enrichment of any complex amplifications in ER+ High and HER2+ tumours")
fisher.test(ct)

res <- list()
for (i in c('ER+ High','HER2+','IC10','ER+ Typical')) {
	tmp <- plot_data[plot_data$group == i,]
	ct <- data.frame(
		primary = c(
			sum(tmp[tmp$stage == 'primary' & tmp$any_complex != 'none','Freq']),
			sum(tmp[tmp$stage == 'primary' & tmp$any_complex == 'none','Freq'])
			),
		metastatic = c(
			sum(tmp[tmp$stage == 'metastatic' & tmp$any_complex != 'none','Freq']),
			sum(tmp[tmp$stage == 'metastatic' & tmp$any_complex == 'none','Freq'])
			)
		)
	stats <- fisher.test(ct)
	res[[i]] <- data.frame(
		subtype = i,
		or = stats$estimate[[1]],
		p = stats$p.value
		)
}
res <- do.call(rbind, res)
print("Enrichment of any complex amplification in metastatic vs primary tumours by subtype")
print(res)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  plot_data,
  file = file.path(main_repo_path, 'data', 'Figure3a_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )

