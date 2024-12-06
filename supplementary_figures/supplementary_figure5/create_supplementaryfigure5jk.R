### SUPPLEMENTARY FIGURE 5JK #####################################################################
# create scatterplot of age and clock-like mutations

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyverse)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### SUPPLEMENTARY FIGURE 5J #######################################################################
# read in source data
dtf <- read.delim(
    file.path(main_repo_path, 'data', 'SupplementaryFigure5j_sourcetable.txt'),
    as.is = TRUE
    )

# read in eniclust 
eniclust <- read.delim(
        file.path(main_repo_path, 'data', 'tcga_predictions.tsv'),
        as.is = TRUE
        )
dtf <- merge(dtf, eniclust[,c('Sample','voting')], by ='Sample')
dtf$group <- (dtf$voting == 'ic9')*1

dtf <- dtf[order(-dtf$group),]

col <- rep('black', nrow(dtf))
col[which(dtf$group == 1)] <- '#ee82edff'

alpha <- rep(0.25, nrow(dtf))
alpha[which(dtf$group == 1)] <- 1

create.scatterplot(
        log10(myc_rna) ~ myc_cn,
        data = dtf,
        col = col,
        alpha = alpha,
        ylab.label = expression('log'['10']*'MYC mRNA'),
        xlab.label = expression('MYC CN'),
        filename = 'supplementary_figure4j.pdf',
        add.points = TRUE,
        points.x = dtf[which(dtf$group == 1),'myc_cn'],
        points.y = log(dtf[which(dtf$group == 1),'myc_rna']),
        points.col = '#ee82edff',
        points.col.border = '#ee82edff',
        #points.cex = 0.75,
        key = list(
             text = list(
                 lab = c('IC9','Other'),
                 cex = 1.5,
                 col = 'black'
                 ),
             points = list(
                 pch = 19,
                 col = c('#ee82edff','black'),
                 cex = 1.5
                 ),
             x = 0.99,
             y = 0.01,
             corner = c(1,0),
             padding.text = 2
             ),
        legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = dtf$myc_rna,
                         y = dtf$myc_cn,
                         label.items = c('spearman'),
                         alpha.background = 0,
                         key.cex = 1.5
                         )
                     ),
                 x = 0.99,
                 y = 0.99,
                 corner = c(1,1)
                 ),
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = dtf[which(dtf$group == 1),'myc_rna'],
                         y = dtf[which(dtf$group == 1),'myc_cn'],
                         label.items = c('spearman'),
                         alpha.background = 0,
                         key.cex = 1.5
                         )
                     ),
                 x = 0.99,
                 y = 0.90,
                 corner = c(1,1)
                 )
             ),
        resolution = 300
        )

### SUPPLEMENTARY FIGURE 5K #######################################################################
# read in source data
dtf <- read.delim(
    file.path(main_repo_path, 'data', 'SupplementaryFigure5k_sourcetable.txt'),
    as.is = TRUE
    )

# read in megatable 
mtable <- read.delim(
        file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
        as.is = TRUE
        )
colnames(mtable) <- gsub('Sample','sample', colnames(mtable))
dtf <- merge(dtf, mtable[,c('sample','ENiClust')], by = 'sample')
# identify IC9
dtf$group <- (dtf$ENiClust == 'ic9')*1

dtf <- dtf[order(-dtf$group),]

col <- rep('black', nrow(dtf))
col[which(dtf$group == 1)] <- '#ee82edff'

alpha <- rep(0.25, nrow(dtf))
alpha[which(dtf$group == 1)] <- 1

dtf$rna <- as.numeric(dtf$rna)
create.scatterplot(
        log10(myc_rna) ~ myc_cn,
        data = dtf,
        col = col,
        alpha = alpha,
        xlimits = c(0,55),
        xat = seq(0,50,10),
        ylab.label = expression('log'['10']*'MYC mRNA'),
        xlab.label = expression('MYC CN'),
        filename = 'supplementary_figure5k.pdf',
        add.points = TRUE,
        points.x = dtf[which(dtf$group == 1),'myc_cn'],
        points.y = log(dtf[which(dtf$group == 1),'myc_rna']),
        points.col = '#ee82edff',
        points.col.border = '#ee82edff',
        #points.cex = 0.75,
        key = list(
             text = list(
                 lab = c('IC9','Other'),
                 cex = 1.5,
                 col = 'black'
                 ),
             points = list(
                 pch = 19,
                 col = c('#ee82edff','black'),
                 cex = 1.5
                 ),
             x = 0.99,
             y = 0.01,
             corner = c(1,0),
             padding.text = 2
             ),
        legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = dtf$myc_cn,
                         y = dtf$myc_rna,
                         label.items = c('spearman'),
                         alpha.background = 0,
                         key.cex = 1.5
                         )
                     ),
                 x = 0.99,
                 y = 0.99,
                 corner = c(1,1)
                 ),
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = dtf[which(dtf$group == 1),'myc_rna'],
                         y = dtf[which(dtf$group == 1),'myc_cn'],
                         label.items = c('spearman'),
                         alpha.background = 0,
                         key.cex = 1.5
                         )
                     ),
                 x = 0.99,
                 y = 0.90,
                 corner = c(1,1)
                 )
             ),
        resolution = 300
        )
