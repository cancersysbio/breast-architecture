### CREATE EXTENDED DATA 3F ###############################################################
# create extended data 3c
# boxplot of SV signatures by subtype in metastatic disease 
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
        file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
        as.is = TRUE
        )

# rename SV signatures to match correlated sig in primary 
colnames(megatable) <- gsub('SV32A','RS6', colnames(megatable))
colnames(megatable) <- gsub('SV32B','RS2', colnames(megatable))
colnames(megatable) <- gsub('SV32C','RS4', colnames(megatable))
colnames(megatable) <- gsub('SV32D','RS5', colnames(megatable))
colnames(megatable) <- gsub('SV32E','RS1', colnames(megatable))
colnames(megatable) <- gsub('SV32F','RS3', colnames(megatable))

# create plot data 
plot_data <- gather(
        megatable[,c('Sample','RS1','RS2','RS3','RS4','RS5','RS6','group')],
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
create.boxplot(
        logvalue ~ group | key,
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
        xlab.label = 'Metastatic',
        filename = 'extended_data3f.pdf',
        add.stripplot = TRUE,
        xaxis.cex = 1,
        resolution = 300
        )

res <- list()
for (j in unique(plot_data$group[!is.na(plot_data$group)])) {
        tmpres <- do.call(rbind, sapply(
                unique(plot_data$key),
                function(i) {
                        tmp <- plot_data[plot_data$key ==  i,]
                        stats <- wilcox.test(tmp[which(tmp$group == j),'value'], tmp[which(tmp$group != j),'value'],)
                        data.frame(
                                svsig = i,
                                es = median(tmp[which(tmp$group == j),'value'], na.rm = TRUE)-median(tmp[which(tmp$group != j),'value'], na.rm = TRUE),
                                p = stats$p.value
                                )
                        },
                simplify = FALSE
                ))
        tmpres$subtype <- j
        res[[j]] <- tmpres
}
res <- do.call(rbind, res)
res$fdr <- p.adjust(res$p, method = 'fdr')



tmpres <- do.call(rbind, sapply(
                unique(plot_data$key),
                function(i) {
                        tmp <- plot_data[plot_data$key ==  i,]
                        stats <- wilcox.test(tmp[which(tmp$group %in% c('ER+ High','HER2+')),'value'], tmp[which(!tmp$group %in% c('ER+ High','HER2+')),'value'],)
                        data.frame(
                                svsig = i,
                                es = median(tmp[which(tmp$group %in% c('ER+ High','HER2+')),'value'], na.rm = TRUE)-median(tmp[which(!tmp$group %in% c('ER+ High','HER2+')),'value'], na.rm = TRUE),
                                p = stats$p.value
                                )
                        },
                simplify = FALSE
                ))
tmpres$fdr <- p.adjust(tmpres$p, method = 'fdr')

### SAVE ##########################################################################################
# write source data to file 
write.table(
  plot_data,
  file = file.path(main_repo_path, 'data', 'ExtendedData3f_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )
