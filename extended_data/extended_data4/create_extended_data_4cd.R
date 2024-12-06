### CREATE EXTENDED DATA 4CD ######################################################################
### CLASSIFY SVS ##################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### EXTENDED DATA 4C ##############################################################################
# check if jabba results have been generated
# file.path(basedir,'svs/jabba/primary/2023-05-08_jabba_events.txt')
jabba_file <- file.path(main_repo_path, 'data','jabba_events.txt')
if (!file.exists(jabba_file)) {
	stop("Please provide summarized matrix of jabba events as jabba_events.txt in the data direcotry...")
}
# read in jabba 
jabba <- read.delim(
	jabba_file,
	as.is = TRUE
	)

# read in primary megatable
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
colnames(megatable) <- gsub('Sample','sample',colnames(megatable))
megatable <- megatable[which(megatable$group != ''),]

eventsdf <- merge(jabba[,1:15], megatable[,c('sample','ENiClust','group','AA_ecdna','AA_complex')], by = 'sample')
eventsdf$AA <- (eventsdf$AA_ecdna > 0)*1
eventsdf$jabba <- 'none'
eventsdf[which(eventsdf$dm > 0 | eventsdf$cpxdm > 0),'jabba'] <- 'dm'

plot_data <- as.data.frame(table(eventsdf[,c('jabba','group')]))
plot_data$prop <- NA
for (i in unique(plot_data$group)) {
	plot_data[plot_data$group == i,'prop'] <- plot_data[plot_data$group == i,'Freq']/sum(plot_data[plot_data$group == i,'Freq'])
}

plot_data$code <- NA
plot_data[plot_data$jabba == 'none','code'] <- 'c'
plot_data[plot_data$jabba == 'dm','code'] <- 'b'

# create barplot
p1 <- create.barplot(
	prop ~ group,
	groups = code,
	col = c('#3d85c6ff','white'),
	stack = TRUE,
	xaxis.lab = rep('', 100),
	#col = c('#8fdab0','#9ebff1ff','#fe85cfff'),
	data = plot_data,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	height = 8,
	width = 8,
	ylab.label = 'Proportion of Samples',
	xlab.label = '',
	#stack= TRUE,
	resolution = 300
	)
p2 <- create.barplot(
	Freq ~ group,
	groups = code,
	col = c('#3d85c6ff','white'),
	data = plot_data,
	stack = TRUE,
	xat = 1:5,
	xaxis.lab = unique(plot_data$group),
	ylimits = c(0,300),
	yat = seq(0,300,50),
	height = 8,
	width = 8,
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('#3d85c6ff')
                             ),
                         text = list(
                             lab = c('Cyclic')
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
                 )
             ),
	ylab.label = 'Number of Samples',
	xlab.label = 'Subtype',
	#stack= TRUE,
	resolution = 300
	)
create.multipanelplot(
	list(p1,p2),
	filename = 'extended_data4c.pdf',
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData4c_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### EXTENDED DATA 4D ##############################################################################
# create extended data 4e
missing <- eventsdf[eventsdf$jabba == 'none' & eventsdf$AA == 1,]
missing$category <- apply(
	missing,
	1,
	function(x) {
		if (x['tyfonas'] > 0 & x['bfb'] == 0) {
			return('tyfona')
		} else if (x['tyfonas'] == 0 & x['bfb'] > 0) {
			return('bfb')
		} else if (x['tyfonas'] > 0 & x['bfb'] > 0) {
			return('both')
		} else {
			return('neither')
		}
	})

plot_data <- as.data.frame(table(missing[,c('category','group')]))
plot_data$prop <- NA
for (i in unique(plot_data$group)) {
	plot_data[plot_data$group == i,'prop'] <- plot_data[plot_data$group == i,'Freq']/sum(plot_data[plot_data$group == i,'Freq'])
}

create.barplot(
	prop ~ group,
	groups = plot_data$category,
	stack = TRUE,
	col = default.colours(4),
	ylimits = c(0,1),
	yat = seq(0,1, 0.2),
	data = plot_data,
	ylab.label = 'Proportion',
	xlab.label = '',
	xaxis.lab = c('ER+\nHigh','ER+\nTypical','HER2+','IC10','IC4ER-'),
	filename = 'extended_data4d.pdf',
	add.text = TRUE,
	text.labels = plot_data$Freq[c(4,3,2,1,8,7,6,5,12,10,9,15,14,13,17)],
	text.x =c(rep(1:2,each = 4), rep(3:4, each = 3), 5),
	text.y = c(
		0.93, 0.78, 0.62, 0.3,
		0.92, 0.45, 0.2, 0.11, 
		0.75, 0.47, 0.2, 
		0.8, 0.45, .2, 
		0.5),
	text.col = c('black', rep('white',3),'black',rep('white', 3), 'black',rep('white',6)),
	legend = list(
             right = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = default.colours(4)
                             ),
                         text = list(
                             lab = c('bfb','both','neither','tyfona')
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     )
                 )
             ),
	width = 8,
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData4d_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

