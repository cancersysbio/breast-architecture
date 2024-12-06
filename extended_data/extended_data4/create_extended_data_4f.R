### EXTENDED DATA 4F ##############################################################################
# create barplot of ecDNA in DCIS
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
seq_plot_data <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData4f_sourcetable.txt'),
	as.is = TRUE
	)

# create barplot
p1 <- create.barplot(
	prop ~ group,
	groups = amp_code,
	col = c('#519a71ff','#8fdab0','white'),
	stack = TRUE,
	xaxis.lab = rep('', 100),
	data = dcis_plot_data,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	height = 8,
	width = 8,
	ylab.label = 'Proportion of Samples',
	xlab.label = '',
	resolution = 300
	)
p2 <- create.barplot(
	prop ~ seq,
	groups = amp_code,
	col = c('#519a71ff','#8fdab0','white'),
	stack = TRUE,
	xaxis.lab = rep('', 100),
	data = seq_plot_data,
	yaxis.lab =rep('', 100),
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	height = 8,
	width = 8,
	ylab.label = '',
	xlab.label = '',
	resolution = 300
	)
p3 <- create.barplot(
	Freq ~ group,
	groups = amp_code,
	col = c('#519a71ff','#8fdab0','white'),
	data = dcis_plot_data,
	stack = TRUE,
	xaxis.lab = c('ER+\nHigh','ER+\nTypical','HER2+','IC10','IC4ER-'),
	ylimits = c(0,410),
	yat = seq(0,400,50),
	height = 8,
	width = 8,
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
                             fill = c('#519a71ff','#8fdab0')
                             ),
                         text = list(
                             lab = c('Cyclic','Complex non-cyclic')
                             ),
                         title = 'DCIS',
                         padding.text = 5,
                         cex = 1.2
                         )
                     ),
                     # Positioning legend on plot
                     x = 0.99,
                     y = 0.99,
                     corner = c(1,1)
                 )
             ),
	#stack= TRUE,
	resolution = 300
	)
p4 <- create.barplot(
	x ~ seq,
	groups = amp_code,
	col = c('#519a71ff','#8fdab0','white'),
	data = seq_plot_data,
	stack = TRUE,
	xaxis.lab = c('sWGS','WGS'),
	ylimits = c(0,410),
	yat = seq(0,400,50),
	yaxis.lab =rep('', 100),
	height = 8,
	width = 8,
	ylab.label = '',
	xlab.label = 'Coverage',
	resolution = 300
	)
create.multipanelplot(
	list(p1,p2,p3,p4),
	plot.objects.heights = c(1,1),
	plot.objects.widths = c(0.75,0.25),
	layout.width = 2,
	layout.height = 2,
	filename = 'extended_data4f.pdf',
	resolution = 300
	)



