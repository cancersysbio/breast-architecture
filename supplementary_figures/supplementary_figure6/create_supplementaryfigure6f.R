### SUPPLEMENTARY FIGURE 6f #####################################################################
# create supplementary figure 6f
# plot replication timing in ecDNA vs non-cyclic

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyverse)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}


date <- Sys.Date()
### MAIN ##########################################################################################
# read in WA 
plot_data <- read.csv(
	file.path(main_repo_path, 'data', 'SupplementaryFigure6f_sourcetable.csv'),
	as.is = TRUE,
	header = TRUE
	)
plot_data[plot_data$Cycle.Type == 'Complex non-cyclic','Cycle.Type'] <- 'Non-cyclic' 
plot_data$WA <- as.numeric(as.character(plot_data$WA))
plot_data$index <- paste0(plot_data$Subtype, plot_data$Cycle.Type)

res <- list()
for (i in unique(plot_data$Subtype)) {
	tmp <- plot_data[which(plot_data$Subtype == i),]
	stats <- wilcox.test(
		tmp[tmp$Cycle.Type == 'Cyclic','WA'],
		tmp[tmp$Cycle.Type == 'Non-cyclic','WA']
		)
	es <- median(tmp[tmp$Cycle.Type == 'Cyclic','WA'])-median(tmp[tmp$Cycle.Type == 'Non-cyclic','WA'])
	res[[i]] <- data.frame(
		subtype = i,
		es = es,
		p = stats$p.value
		)
}
res <- do.call(rbind, res)

# create plot
create.boxplot(
			WA ~ index,
			data = plot_data,
			add.stripplot = TRUE,
			xaxis.cex = 1,
			xaxis.lab = rep(c('Cyclic','Non-cyclic'), 5),
			ylimits = c(0,43),
			yat = seq(0,40,10),
			ylab.label = 'Replication Timing\nWeighted Average',
			xlab.label = 'Amplicon Type',
			xaxis.rot = 45,
			add.text = TRUE,
			text.labels = c('ER+ High','  ER+\nTypical','HER2+','IC10/IC4ER-',
				paste0('ES=', round(res[res$subtype == 'ER+ High Risk','es'], 2)), paste0('P=', round(res[res$subtype == 'ER+ High Risk','p'], 2)),
				paste0('ES=', round(res[res$subtype == 'ER+ Typical Risk','es'], 2)), paste0('P=', round(res[res$subtype == 'ER+ Typical Risk','p'], 2)),
				paste0('ES=', round(res[res$subtype == 'HER2+','es'], 2)), paste0('P=', round(res[res$subtype == 'HER2+','p'], 2)),
				paste0('ES=', round(res[res$subtype == 'IC4ER-/IC10','es'], 2)), paste0('P=', round(res[res$subtype == 'IC4ER-/IC10','p'], 2))
				),
			add.rectangle = TRUE,
			xleft.rectangle = c(2.5,6.5),
			ybottom.rectangle = 0,
			xright.rectangle = c(0,4.5),
			ytop.rectangle = 1000+10,
			col.rectangle = 'grey75',
			alpha.rectangle = 0.25,
			text.x = c(seq(1.5,8,2),rep(seq(1.5,8,2), each = 2)),
			text.y = c(rep(40,4), rep(c(36.5,34.5), 4)),
			text.anchor = 'centre',
			text.col = 'black',
			text.cex = 1.1,
			text.fontface = 'bold',
			filename = 'supplementary_figure6f.pdf',
			resolution = 300
			)
