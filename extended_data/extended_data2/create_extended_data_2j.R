<<<<<<< HEAD:extended_data/extended_data2/create_extended_data_2j.R
### CREATE EXTENDED DATA 2J #######################################################################
# create barplot of jabba events across subgroups 

### PREAMEBLE #####################################################################################
library(gGnome)
library(fishHook)
library(tidyr)

=======
### CREATE EXTENDED DATA FIGURE 2i #############################################################################
# creates extended data figure 2i
# extended data figure 2i provides the Copy number and SV profiles of primary and metastatic samples, each representative of either a TNBC -enriched, ER+ Typical -enriched, ER+ High/HER2+ -enriched, or of mixed profile in the center of the Pareto front.

### PREAMBLE #####################################################################################
library(yaml)
library(tidyverse)
library(readxl)
library(ggpubr)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)
library(circlize)
library(RColorBrewer)
library(dinamic, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")

>>>>>>> c398e855a1f91e3baab879e682684528cdf0bfff:extended_data/extended_data2/create_extended_data_2i.R
# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

<<<<<<< HEAD:extended_data/extended_data2/create_extended_data_2j.R
date <- Sys.Date()
### MAIN ##########################################################################################
# read in jabba
jabba_file <- file.path(main_repo_path, 'data', 'jabba_events.txt')
if (!file.exists(jabba_file)) {
	stop("Please put JaBbA events as outputted by JaBbA in data directory ...")
}
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

eventsdf <- merge(jabba[,1:15], megatable[,c('sample','ENiClust','group')], by = 'sample')

# not considering small, simple svs that dominate (i.e. del, inv, dup)
plot_data <- aggregate(eventsdf[,c(2:5,7,10:15)], list(eventsdf$group), sum)
#plot_data[,2:13] <- plot_data[,2:13]/rowSums(plot_data[,2:13])
plot_data <- gather(plot_data, key = 'event', value = 'count', -Group.1)

# calculate proportion 
plot_data$prop <- NA
for (i in unique(plot_data$Group.1)) {
	plot_data[which(plot_data$Group.1 == i),'prop'] <- plot_data[which(plot_data$Group.1 == i),'count']/sum(plot_data[which(plot_data$Group.1 == i),'count'])
}

# create plot
create.barplot(
	prop ~ Group.1,
	groups = plot_data$event,
	data = plot_data,
	filename = 'extended_data2j.pdf',
	stack = TRUE,
	xaxis.cex = 1.2,
	xaxis.rot = 45,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	xlab.label = 'Subtype',
	ylab.label = 'Proportion of Events',
	legend = list(
             right = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = default.colours(11)
                             ),
                         text = list(
                             lab = unique(plot_data$event)[order(unique(plot_data$event))]
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     )
                 )
             ),
	col = default.colours(11), 
	resolution = 300
	)

### WRITE TABLE ###################################################################################
# write to file
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData2j_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
=======
config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

colCNVs = c( colorRampPalette(c(rgb(0.8,0.2,0.1), rgb(1,1,1)))(100) , colorRampPalette(c(rgb(1,1,1), rgb(0,0.4,0.7,1) ))(401) )

# Plotting function
plot_crc1 <- function(samples){
  for(sample in samples){
    print(sample)
    
    circos.clear()
    col_text <- rgb(0.2,0.2,0.2)
    circos.par("track.height" = 0.20, start.degree = 90, circle.margin = 0.35 , 
               gap.degree = 3,
               cell.padding = c(0, 0, 0, 0))
    circos.initialize(sectors = rownames(hg19_bed), xlim = hg19_bed)
    
    data_cytoband = data.frame(chr = c("chr8","chr8","chr11","chr17","chr17","chr20"),
                               start = lim_start,
                               end = lim_end,
                               labels = c("8p12 (IC6)","8q24 (IC9)","11q13 (IC2)","17q12 (IC5)","17q23 (IC1)","20q13 (IC1/IC9)"), stringsAsFactors = F)
    
    circos.genomicLabels(data_cytoband, labels.column=4, cex=0.8, col="black", line_lwd=1, line_col="black", #"grey80",
                         side="outside", connection_height=0.07, labels_height=0.04)
    
    circos.track(ylim=c(0,1),panel.fun=function(x,y) {
      chr=str_remove(CELL_META$sector.index,"chr")
      xlim=CELL_META$xlim
      ylim=CELL_META$ylim
      circos.text(mean(xlim),mean(ylim),chr,cex=0.8,col="#333333",
                  facing="bending.inside",niceFacing=TRUE)
    },bg.col=rep(c(rgb(0.85,0.85,0.85),rgb(0.7,0.7,0.7)),12),bg.border=F,track.height=0.06)
    
    # Add CN levels
    CNVs_org.tmp = dataset.cn %>% dplyr::filter(ID == sample) %>% 
      dplyr::mutate(chr=paste0("chr",chrom),start=loc.start,end=loc.end,value=1,value1=tcn.em-lcn.em,value2=lcn.em)
    
    circos.genomicTrack(CNVs_org.tmp %>% dplyr::select(chr,start,end,value,value1,value2), ylim=c(0,4), bg.border=NA, 
                        panel.fun = function(region, value, ...) {
                          region1 = region[value$value1!=1,]
                          region2 = region[value$value2!=1,]
                          if(nrow(region1)>0) circos.genomicRect(region1,value[value$value1!=1,1],
                                                                 col=colCNVs[round(value$value1[value$value1!=1]*100+1)],
                                                                 border = NA, ytop = 2,ybottom = 0 )
                          if(nrow(region2)>0) circos.genomicRect(region2,value[value$value2!=1,1],
                                                                 col=colCNVs[round(value$value2[value$value2!=1]*100+1)],
                                                                 border = NA, ytop = 4,ybottom = 2 )
                        })
    
    # Add SVs
    Link1 = dataset.sv[as.data.frame(dataset.sv)[,"sample"]==sample,] %>% dplyr::mutate(chrom=paste0("chr",chrom1),start=pos1,end=pos1_end) %>% dplyr::select(chrom,start,end) %>% as.data.frame()
    Link2 = dataset.sv[as.data.frame(dataset.sv)[,"sample"]==sample,] %>% dplyr::mutate(chrom=paste0("chr",chrom2),start=pos2,end=pos2_end) %>% dplyr::select(chrom,start,end) %>% as.data.frame()
    circos.genomicLink(Link1,Link2,col=c(alpha("#a3a3a3",0.7),"#000000")[as.numeric(dataset.sv[as.data.frame(dataset.sv)[,"sample"]==sample,]$chrom_type=="inter")+1])
    
    #title(paste0(names(samples)[which(samples == sample)], " sample ", paste(unlist(strsplit(sample,"_"))[2:3],collapse = "_")," (", toupper(mega$eniclust[which(mega$Sample == sample)]), ")"))
    title(paste0(names(samples)[which(samples == sample)], " sample ", paste(unlist(strsplit(sample,"_"))[2:3],collapse = "_")))
  }
}

### MAIN ##########################################################################################
basedir <- "/oak/stanford/groups/ccurtis2/users/"
outdir <- paste0(basedir,"lisem/bc_landscape/submission/figures/")

## Load genomic data
chr_sizes = read.table(paste0(basedir,"srinivap/IC_classifiers/Data/General/hg19.chrom.sizes.txt"), header = F, sep = "\t")
hg19_bed = data.frame(start=0,end=chr_sizes$V2[which(chr_sizes$V1 %in% paste0("chr",c(1:22,"X","Y","M")))],row.names = paste0("chr",c(1:22,"X","Y","M")) )

data(annot.file)
lim_start = sapply(c("8p12","8q24","11q13","17q12","17q23","20q13"), function(x) min(annot.file$start[which(startsWith(paste0(annot.file$chr,annot.file$cytoband), x))]))
lim_end = sapply(c("8p12","8q24","11q13","17q12","17q23","20q13"), function(x) max(annot.file$end[which(startsWith(paste0(annot.file$chr,annot.file$cytoband), x))]))

## Load metadata
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])
mega$ENiClust[which(mega$ENiClust == "Other")] = "IC3/IC7"

## Load data
dataset.cn_prim = read.delim(paste0(basedir,"lisem/ICGC/02_segment.tsv"))
dataset.cn_mets = read.delim(paste0(basedir,"lisem/Hartwig/eniclust/02_segment.txt"))
dataset.sv_prim = read.delim(paste0(basedir,"khoulaha/BreastLandscape/data/sv_pairs.bedpe"))
dataset.sv_mets = read.delim(paste0(basedir,"khoulaha/BreastLandscape/data/sv_pairs_project17.tsv"))

dataset.cn <- rbind(dataset.cn_prim, dataset.cn_mets)
dataset.cn$ID = gsub("_hisens","",dataset.cn$ID)
dataset.sv <- rbind(dataset.sv_prim, dataset.sv_mets)

dataset.sv$chrom_type = ifelse(dataset.sv$chrom1 == dataset.sv$chrom2,"intra","inter")

### SAVE ##########################################################################################
## Figure----
# Primary
samples_prim = c("TNBC -enriched" = "CURTIS_H002087_T01_01_WG01",
            "ER+ Typical -enriched" = "CURTIS_H000491_T01_01_WG01",
            "ER+ High/HER2+ -enriched" = "CURTIS_H001567_T01_01_WG01",
            "Mixed" = "CURTIS_H004066_T01_01_WG01")

png(filename = file.path(outdir, 'Circos_prim.png'), res = 300, width = 8.85, height = 6.59, units = 'in')
par(mfrow=c(2,2),mar = c(3, 3, 3, 3)*0.5)
plot_crc1(samples=samples_prim)
dev.off()

# Metastatic
samples_met = c("TNBC -enriched" = "CURTIS_H004112_T01_01_WG01",
            "ER+ Typical -enriched" = "CURTIS_H004163_T01_01_WG01",
            "ER+ High/HER2+ -enriched" = "CURTIS_H004681_T01_01_WG01",
            "Mixed" = "CURTIS_H004010_T01_01_WG01")

png(filename = file.path(outdir, 'Circos_mets.png'), res = 300, width = 8.85, height = 6.59, units = 'in')
par(mfrow=c(2,2),mar = c(3, 3, 3, 3)*0.5)
plot_crc1(samples=samples_met)
dev.off()

## SourceData---- 
data_cn <- data.frame()
data_sv <- data.frame()
samples <- list(Primary = samples_prim, Metastatic = samples_met)
for(stage in names(samples)){
  for(sample in samples[[stage]]){
    print(sample)
    
    # Add CN levels
    CNVs_org.tmp = dataset.cn %>% dplyr::filter(ID == sample) %>% 
      dplyr::mutate(chr=paste0("chr",chrom),start=loc.start,end=loc.end,value=1,value1=tcn.em-lcn.em,value2=lcn.em)
    
    data_cn_tmp <- CNVs_org.tmp %>% dplyr::select(chr,start,end,value,value1,value2) %>% dplyr::mutate(Sample = sample, Stage = stage)
    data_cn <- rbind(data_cn, data_cn_tmp)
    
    
    # Add SVs
    Link1 = dataset.sv[as.data.frame(dataset.sv)[,"sample"]==sample,] %>% dplyr::mutate(chrom=paste0("chr",chrom1),start=pos1,end=pos1_end) %>% dplyr::select(chrom,start,end) %>% as.data.frame()
    Link2 = dataset.sv[as.data.frame(dataset.sv)[,"sample"]==sample,] %>% dplyr::mutate(chrom=paste0("chr",chrom2),start=pos2,end=pos2_end) %>% dplyr::select(chrom,start,end) %>% as.data.frame()
    
    data_sv_tmp <- cbind(Link1 %>% dplyr::mutate(Sample = sample, breakpoint = "B1", chrom1 = chrom, start1 = start, end1 = end) %>% dplyr::select(Sample,breakpoint,chrom1,start1,end1),
                         Link2 %>% dplyr::mutate(breakpoint = "B2", chrom2 = chrom, start2 = start, end2 = end, Stage = stage) %>% dplyr::select(breakpoint,chrom2,start2,end2))
    data_sv <- rbind(data_sv, data_sv_tmp)
  }
}
write.table(data_cn, paste0(basedir,"lisem/bc_landscape/github/ExtendedData2j_CNV_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
write.table(data_sv, paste0(basedir,"lisem/bc_landscape/github/ExtendedData2j_SV_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

>>>>>>> c398e855a1f91e3baab879e682684528cdf0bfff:extended_data/extended_data2/create_extended_data_2i.R
