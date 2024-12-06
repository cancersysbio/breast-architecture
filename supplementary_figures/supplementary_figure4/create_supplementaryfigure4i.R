### SUPPLEMENTARY FIGURE 4I #######################################################################
# create supplementary figure 4i
# create violin plot of CN signatures in metastatic disease
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
        file.path(basedir, 'data', 'metastatic_megatable.txt'),
        as.is = TRUE
        )

# create plot data 
plot_data <- gather(
        megatable[,c('Sample','CN1','CN2','CN3','CN7','CN9','CN11','CN17','CN20','CN21','group')],
        key = 'key',
        value = 'value',
        -Sample,
        -group
        )
plot_data$logvalue <- log10(plot_data$value+1)
plot_data$group <- gsub('Low','Typical', plot_data$group)

col <- rep(NA, nrow(plot_data))
col[which(plot_data$group == 'ER+ High')] <- '#fa954eff'
col[which(plot_data$group == 'ER+ Typical')] <- '#a1c1d4ff'
col[which(plot_data$group == 'HER2+')] <- '#8b0100ff'
col[which(plot_data$group == 'IC10')] <- '#7c26ccff'
col[which(plot_data$group == 'IC4ER-')] <- '#c2b7dbff'

# create boxplot 
create.violinplot(
        logvalue ~ group | key,
        data = plot_data,
        xaxis.rot = 45,
        ylab.label = 'Activity',
        yat = 0:3,
        col = c('#fa954eff','#a1c1d4ff','#8b0100ff','#7c26ccff','#c2b7dbff'),
        ylimits = c(0,3.2),
        yaxis.lab = c(
                expression(0),
                expression(10^1),
                expression(10^2),
                expression(10^3)
                ),
        xlab.label = 'Metastatic',
        filename = 'supplementary_figure4i.pdf',
        xaxis.cex = 1,
        resolution = 300
        )