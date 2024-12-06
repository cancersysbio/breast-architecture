### GENERATE EXTENDED DATA 1HI ####################################################################
# create boxplot of mrna in er+ high vs typical

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### EXTENDED DATA 1H ###############################################################################
# read in ER signature 
plot_data <- read.delim(
        file.path(main_repo_path, 'data','na_megatable.txt'),
        as.is = TRUE
        )

plot_data <- plot_data[plot_data$group %in% c('ER+ High','ER+ Typical','HER2+') & grepl('TCGA', plot_data$Sample),]
plot_data$index <- NA
plot_data[which(plot_data$group == 'HER2+' & plot_data$ER == 0),'index'] <- 'a'
plot_data[which(plot_data$group == 'HER2+' & plot_data$ER == 1),'index'] <- 'b'
plot_data[which(plot_data$group == 'ER+ High'),'index'] <- 'c'
plot_data[which(plot_data$group == 'ER+ Typical'),'index'] <- 'd'

stats_highlow <- wilcox.test(
        plot_data[plot_data$group == 'ER+ High','ER_early_sig'],
        plot_data[plot_data$group == 'ER+ Typical','ER_early_sig']
        )
es_highlow <- median(plot_data[plot_data$group == 'ER+ High','ER_early_sig'], na.rm = TRUE)-median(plot_data[plot_data$group == 'ER+ Typical','ER_early_sig'], na.rm = TRUE)

stats_highher2 <- wilcox.test(
        plot_data[plot_data$group == 'ER+ High','ER_early_sig'],
        plot_data[plot_data$group == 'HER2+' & plot_data$ER == 1,'ER_early_sig']
        )
es_highher2 <- median(plot_data[plot_data$group == 'ER+ High' ,'ER_early_sig'], na.rm = TRUE)-median(plot_data[plot_data$group == 'HER2+' & plot_data$ER == 1,'ER_early_sig'], na.rm = TRUE)


stats_highher2neg <- wilcox.test(
        plot_data[plot_data$group == 'ER+ High','ER_early_sig'],
        plot_data[plot_data$group == 'HER2+' & plot_data$ER == 0,'ER_early_sig']
        )
es_highher2neg <- median(plot_data[plot_data$group == 'ER+ High' ,'ER_early_sig'], na.rm = TRUE)-median(plot_data[plot_data$group == 'HER2+' & plot_data$ER == 0,'ER_early_sig'], na.rm = TRUE)


create.boxplot(
        ER_early_sig ~ index,
        data = plot_data,
        add.stripplot = TRUE,
        ylab.label = 'ER Early Signaling',
        xlab.label = 'Subtype',
        ylimits = c(-0.7,0.7),
        yat = seq(-0.6,0.6,0.3),
        add.text = TRUE,
        text.labels = c(
                expression('ES=0.31'), expression('P=3.1x10'^'-10'),
                expression('ES=0.09'), expression('P=0.06'),
                expression('ES=-0.08'), expression('P=7.0x10'^'-5')
                ),
        col.rectangle =  'mediumpurple2',
        add.rectangle = TRUE,
        xleft.rectangle = c(-0.5,2.5),
        ybottom.rectangle = -2,
        xright.rectangle = c(1.5,3.5),
        ytop.rectangle = 4.5,
        alpha.rectangle = 0.25,
        text.x = c(1,1,2, 2, 4,4),
        text.y = c(0.67,0.59,0.67, 0.59, -0.59,-0.67),
        xaxis.lab = c('HER2+\nER-','HER2+\nER+','ER+\nHigh','ER+\nTypical'),
        filename = 'extended_data1h.pdf',
        resolution = 300
        )

### WRITE TO FILE #################################################################################
# write source data to file
write.table(
        plot_data[,c('ER_early_sig','index','group','ER')],
        file = file.path(main_repo_path, 'data', 'ExtendedData1h_sourcetable.txt'),
        sep = '\t',
        row.names = FALSE
        )

#### EXTENDED DATA 1I #############################################################################
# read in megatable
megatable <- read.delim(
        file.path(main_repo_path, 'data', 'primary_megatable.txt'),
        as.is = TRUE
        )
metmegatable <- read.delim(
        file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
        as.is = TRUE
        )
# read in ER signature 
plot_data <- read.delim(
        file.path(main_repo_path, 'data','rna_megatable.txt'),
        as.is = TRUE
        )
plot_data <- plot_data[which(gsub('WT01','WG01',plot_data$Sample) %in% c(megatable$Sample, metmegatable$Sample)),]
plot_data <- plot_data[plot_data$Cohort == 'Hartwig' & plot_data$group %in% c('ER+ High','ER+ Typical','HER2+'),]
plot_data$stage_code <- 'b'
plot_data[which(plot_data$Stage == 'Primary'),'stage_code'] <- 'a'

plot_data$index <- paste0(plot_data$group, plot_data$stage_code)
col <- rep('#9ebff1ff', nrow(plot_data))
col[which(plot_data$Stage == 'Metastatic')] <- '#fe85cfff'

res <- list()
for (i in c('ER+ High','ER+ Typical','HER2+')) {
        stats <- wilcox.test(
                plot_data[plot_data$group == i & plot_data$Stage == 'Primary','ER_early_sig'],
                plot_data[plot_data$group == i & plot_data$Stage == 'Metastatic','ER_early_sig']
                )
        es <- median(plot_data[plot_data$group == i & plot_data$Stage == 'Primary','ER_early_sig'])-median(plot_data[plot_data$group == i & plot_data$Stage == 'Metastatic','ER_early_sig'])
        res[[i]] <- data.frame(
                subtype = i,
                es = es,
                p = stats$p.value
                )
}
res <- do.call(rbind,res)


create.boxplot(
        ER_early_sig ~ index,
        data = plot_data,
        points.col = col,
        ylimits = c(-0.7, 0.7),
        xaxis.lab = rep(c("Primary","Metastatic"), 3),
        add.text = TRUE,
        text.labels = c(
                'ER+ High','ER+ Typical','HER2+',
                'ES=0.14', 'P=0.46',
                'ES=0.07', 'P=0.64',
                'ES=0.20', 'P=0.92'
                ),
        text.x = c(1.5, 3.5, 5.5, 1.5, 1.5, 3.5, 3.5, 5.5, 5.5),
        text.y = c(rep(0.66, 3), rep(c(-0.6, -0.65),3)), 
        abline.v = c(2.5, 4.5),
        abline.lty = 2,
        ylab.label = 'ER Early Signaling',
        xlab.label = 'Stage of Progression',
        xaxis.rot = 45,
        xaxis.cex = 1,
        yaxis.cex = 1,
        filename = paste0(date, '_ER_signaling_primary_vs_met_boxplot.pdf'),
        add.stripplot = TRUE,
        resolution = 300
        )

### WRITE TO FILE #################################################################################
# write source data to file
write.table(
        plot_data[,c('ER_early_sig','index','stage_code','Stage','group')],
        file = file.path(main_repo_path, 'data', 'ExtendedData1i_sourcetable.txt'),
        sep = '\t',
        row.names = FALSE
        )

