### CREATE EXTENDED DATA 2GH ######################################################################
# create extended data 2gh
# boxplot of SV signatures by subtype
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
megatable <- megatable[which(megatable$group != ''),]
# create plot data 
plot_data <- gather(
	megatable[,c('Sample','SV32A','SV32B','SV32C','SV32D','SV32E','SV32F','group')],
	key = 'key',
	value = 'value',
	-Sample,
	-group
	)
plot_data[plot_data$key == 'SV32A','key'] <- 'RS6'
plot_data[plot_data$key == 'SV32B','key'] <- 'RS2'
plot_data[plot_data$key == 'SV32C','key'] <- 'RS3'
plot_data[plot_data$key == 'SV32D','key'] <- 'RS4'
plot_data[plot_data$key == 'SV32E','key'] <- 'RS5'
plot_data[plot_data$key == 'SV32F','key'] <- 'RS1'
plot_data$logvalue <- log10(plot_data$value+1)


col <- rep(NA, nrow(plot_data))
col[which(plot_data$group == 'ER+ High')] <- '#fa954eff'
col[which(plot_data$group == 'ER+ Typical')] <- '#a1c1d4ff'
col[which(plot_data$group == 'HER2+')] <- '#8b0100ff'
col[which(plot_data$group == 'IC10')] <- '#7c26ccff'
col[which(plot_data$group == 'IC4ER-')] <- '#c2b7dbff'

res <- do.call(rbind, sapply(
	unique(plot_data$key),
	function(i) {
		tmp <- plot_data[plot_data$key ==  i,]
		stats <- kruskal.test(tmp$value, tmp$group)
		data.frame(
			svsig = i,
			p = stats$p.value
			)
		},
	simplify = FALSE
	))

# create boxplot 
create.boxplot(
	logvalue ~ group | key,
	data = plot_data,
	xaxis.rot = 45,
	ylab.label = 'Activity',
	yat = 0:3,
	ylimits = c(0,3.2),
	points.col = col,
	yaxis.lab = c(
		expression(0),
		expression(10^1),
		expression(10^2),
		expression(10^3)
		),
	xlab.label = 'Primary',
	filename = 'extended_data2g.pdf',
	add.stripplot = TRUE,
	xaxis.cex = 1,
	resolution = 300
	)

# test stats
res <- do.call(rbind, sapply(
	unique(plot_data$key),
	function(i) {
		tmp <- plot_data[plot_data$key ==  i,]
		stats <- wilcox.test(tmp[which(tmp$group == 'HER2+'),'value'], tmp[which(tmp$group == 'ER+ High'),'value'],)
		data.frame(
			svsig = i,
			es = median(tmp[which(tmp$group == 'HER2+'),'value'], na.rm = TRUE)-median(tmp[which(tmp$group == 'ER+ High'),'value'], na.rm = TRUE),
			p = stats$p.value
			)
		},
	simplify = FALSE
	))
res$fdr <- p.adjust(res$p, method = 'fdr')


res <- do.call(rbind, sapply(
	unique(plot_data$key),
	function(i) {
		tmp <- plot_data[plot_data$key ==  i,]
		stats <- wilcox.test(tmp[which(tmp$group == 'IC10'),'value'], tmp[which(tmp$group != 'IC10'),'value'],)
		data.frame(
			svsig = i,
			es = median(tmp[which(tmp$group == 'IC10'),'value'], na.rm = TRUE)-median(tmp[which(tmp$group != 'IC10'),'value'], na.rm = TRUE),
			p = stats$p.value
			)
		},
	simplify = FALSE
	))
res$fdr <- p.adjust(res$p, method = 'fdr')

res <- do.call(rbind, sapply(
	unique(plot_data$key),
	function(i) {
		tmp <- plot_data[plot_data$key ==  i,]
		stats <- wilcox.test(tmp[which(tmp$group %in% c('HER2+','ER+ High')),'value'], tmp[which(!tmp$group %in% c('HER2+','ER+ High')),'value'],)
		data.frame(
			svsig = i,
			es = median(tmp[which(tmp$group %in% c('HER2+','ER+ High')),'value'], na.rm = TRUE)-median(tmp[which(!tmp$group %in% c('HER2+','ER+ High')),'value'], na.rm = TRUE),
			p = stats$p.value
			)
		},
	simplify = FALSE
	))
res$fdr <- p.adjust(res$p, method = 'fdr')

### WRITE TO FILE #################################################################################
# write plot data to file
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData2g_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### EXTENDED DATA 2H ##############################################################################
high_erlow <- megatable[which(megatable$group %in% c('ER+ High','ER+ Typical','HER2+')),]
high_erlow[high_erlow$ENiClust %in% c('Other','ic4','ic8'),'ENiClust'] <- 'ER+ Typical'

plot_data <- gather(
	high_erlow[,c('Sample','SV32A','SV32B','SV32C','SV32D','SV32E','SV32F','ENiClust')],
	key = 'key',
	value = 'value',
	-Sample,
	-ENiClust
	)
plot_data[plot_data$key == 'SV32A','key'] <- 'RS6'
plot_data[plot_data$key == 'SV32B','key'] <- 'RS2'
plot_data[plot_data$key == 'SV32C','key'] <- 'RS3'
plot_data[plot_data$key == 'SV32D','key'] <- 'RS4'
plot_data[plot_data$key == 'SV32E','key'] <- 'RS5'
plot_data[plot_data$key == 'SV32F','key'] <- 'RS1'
plot_data$logvalue <- log10(plot_data$value+1)

plot_data$ENiClust <- gsub('^ic','IC', plot_data$ENiClust)
plot_data[which(plot_data$ENiClust == 'IC5'),'ENiClust'] <- 'HER2+'

col <- rep(NA, nrow(plot_data))
col[which(plot_data$ENiClust == 'IC1')] <- '#ff5500ff'
col[which(plot_data$ENiClust == 'IC2')] <- '#00ee77ff'
col[which(plot_data$ENiClust == 'IC6')] <- '#fffe41ff'
col[which(plot_data$ENiClust == 'IC9')] <- '#ee82edff'
col[which(plot_data$ENiClust == 'ER+ Typical')] <- '#a1c1d4ff'
col[which(plot_data$ENiClust == 'HER2+')] <- '#8b0100ff'

# create boxplot 
create.boxplot(
	logvalue ~ ENiClust | key,
	data = plot_data,
	xaxis.rot = 45,
	ylab.label = 'Activity',
	yat = 0:3,
	points.col = col,
	ylimits = c(0,3.2),
	yaxis.lab = c(
		expression(0),
		expression(10^1),
		expression(10^2),
		expression(10^3)
		),
	xlab.label = 'Primary',
	filename = 'extended_data2h.pdf',
	add.stripplot = TRUE,
	xaxis.cex = 1.2,
	resolution = 300
	)

plot_data_er <- plot_data[which(plot_data$ENiClust %in% c('IC1','IC2','IC6','IC9')),]
res <- list()
for (j in c('IC1','IC2','IC6','IC9')) {
	res[[j]] <- do.call(rbind, sapply(
		unique(plot_data$key),
		function(i) {
			tmp <- plot_data_er[plot_data_er$key ==  i,]
			stats <- wilcox.test(tmp[which(tmp$ENiClust == j),'value'], tmp[which(!tmp$ENiClust == j),'value'],)
			data.frame(
				subtype = j,
				svsig = i,
				es = median(tmp[which(tmp$ENiClust == j),'value'], na.rm = TRUE)-median(tmp[which(!tmp$ENiClust == j),'value'], na.rm = TRUE),
				p = stats$p.value
				)
			},
		simplify = FALSE
		))
}
res <- do.call(rbind, res)
res$fdr <- p.adjust(res$p, method = 'fdr')

### WRITE TO FILE #################################################################################
# write plot data to file
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData2h_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)