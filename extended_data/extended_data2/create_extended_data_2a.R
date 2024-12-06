### CREATE EXTENDED DATA FIGURE 2a #############################################################################
# creates extended data figure 2a
# extended data figure 2a provides the IC group-level copy number profile with SV burden as overlay in DCIS samples.

### PREAMBLE #####################################################################################
library(yaml)
library(data.table)
require(GenomicRanges)
require(dplyr)

library(ggplot2)
library(ggpubr)
library(patchwork)
require(ggtext)
require(cowplot)

require(GenVisR, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

col_groups = list("ER+ High" = c("#ff6a00", "#ffb07c"), 
                  "ER+ Typical" = c("#7aa6c2",alpha("#7aa6c2",0.7)), 
                  "HER2+" = c("#8b0100ff","#dd6867ff"), 
                  "TNBC" = c("#6300a1","#bc8ad9"))

### Customized cnFreq() R function from GenVisR
cnFreq_qual <- function(x)
{
  # Check that x is a data frame
  if(!is.data.frame(x)){
    memo <- paste0("Did not detect a data frame in argument supplied",
                   " to x... attempting to coerce")
    warning(memo)
    x <- as.data.frame(x)
    x <- droplevels(x)
  }
  
  # Check that x has at least 1 row
  if(nrow(x) < 1){
    memo <- paste0("x needs at least one row")
    stop(memo)
  }
  
  # remove any NA values in the data
  if(any(is.na(x))){
    na_rows_removed <- nrow(x) - nrow(na.omit(x))
    memo <- paste0("Removing ", na_rows_removed, " rows containing NA values")
    message(memo)
    x <- na.omit(x)
  }
  
  if(all(c('chromosome', 'start','end', 'segmean', 'sample') %in% colnames(x))){
    
    # make sure columns are of the correct type
    x$chromosome <- as.factor(x$chromosome)
    x$start <- as.integer(as.character(x$start))
    x$end <- as.integer(as.character(x$end))
    x$segmean <- as.numeric(as.character(x$segmean))
    x$sample <- as.factor(x$sample)
    
    # make sure windows are consistent if not disjoin them
    tmp <- split(x, x$sample)
    tmp_vec <- tmp[[1]]$end
    if(any(!unlist(sapply(tmp, function(x) x[,"end"] %in% tmp_vec), use.names=F))){
      memo <- paste0("Did not detect identical genomic segments for all samples",
                     " ...Performing disjoin operation")
      message(memo) 
      
      # here we split the DF up in an attempt to avoid complaints that lists are to large
      x <- split(x, f=x$chromosome)
      x <- lapply(x, cnFreq_disjoin)
      
      
      x <- do.call(rbind, x)
    }
    rm(tmp)
    rm(tmp_vec)
  } else {
    memo <- paste0("Did not detect correct columns in argument supplied",
                   " to x!")
    stop(memo)
  }
  
  # Check chromosome column in x
  if(!all(grepl("^chr", x$chromosome))){
    memo <- paste0("Did not detect the prefix \"chr\" in the chromosome",
                   " column of x... adding prefix")
    message(memo)
    x$chromosome <- paste0("chr", x$chromosome)
    x$chromosome <- as.factor(x$chromosome)
  } else if(all(grepl("^chr", x$chromosome))) {
    memo <- paste0("Detected \"chr\" in the chromosome column of x...",
                   " proceeding")
    message(memo)
  } else {
    memo <- paste0("Detected unknown or mixed prefixes in the chromosome ",
                   " column of x, should either have a chr prefix or ",
                   "none at all!")
    stop(memo)
  }
  
  return(x)
}

cnFreq_disjoin <- function(x){
  
  # create the Granges object for the data
  x <- GenomicRanges::GRanges(seqnames=x$chromosome,
                              ranges=IRanges::IRanges(start=x$start, end=x$end),
                              "sample"=x$sample, "segmean"=x$segmean)
  
  # disjoin with grange, get a mapping of meta columns and expand it
  disJoint_x <- GenomicRanges::disjoin(x, with.revmap=TRUE)
  revmap <- GenomicRanges::mcols(disJoint_x)$revmap
  disJoint_x <- rep(disJoint_x, lengths(revmap))
  
  
  # exract the meta columns and map them back to the disJoint GRanges object
  sample <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$sample, revmap))
  segmean <- unlist(IRanges::extractList(GenomicRanges::mcols(x)$segmean, revmap))
  GenomicRanges::mcols(disJoint_x)$sample <- sample
  GenomicRanges::mcols(disJoint_x)$segmean <- segmean
  
  # convert the GRanges Object back to a data frame
  disJoint_x <- as.data.frame(disJoint_x)[,c("seqnames", "start", "end", "width",
                                             "sample", "segmean")]
  colnames(disJoint_x) <- c("chromosome", "start", "end", "width", "sample", "segmean")
  return(disJoint_x)
}

multi_chrBound <- function(x)
{
  # Check that input has size
  if(nrow(x) < 1)
  {
    memo <- paste0("input has 0 rows, it is possible that the UCSC",
                   " MySQL query has failed")
    stop(memo)
  }
  
  # Extract the columns needed
  data <- x[,c('chrom' ,'chromStart' , 'chromEnd')]
  
  # Obtain max for each chromosome
  maxChrom <- stats::aggregate(chromEnd ~ chrom, data=data, max)
  maxChrom <- cbind(maxChrom, maxChrom[,2])
  colnames(maxChrom) <- c('chromosome', 'start', 'end')
  
  # Obtain min for each chromosome
  minChrom <- stats::aggregate(chromStart ~ chrom, data=data, min)
  minChrom <- cbind(minChrom, minChrom[,2])
  colnames(minChrom) <- c('chromosome', 'start', 'end')
  
  # bind all the data together
  data <- rbind(maxChrom, minChrom)
  
  return(data)
}

cnFreq_LM <- function(x, CN_low_cutoff=1.5, CN_high_cutoff=2.5, plot_title=NULL,
                      CN_Loss_colour='#002EB8', CN_Gain_colour='#A30000',
                      x_title_size=12, y_title_size=12, facet_lab_size=10,
                      plotLayer=NULL, plotType="proportion", genome="hg19",
                      plotChr=NULL, out="plot")
{
  # Perform quality check on input data
  x <- cnFreq_qual(x)
  samples <- unique(x$sample)
  
  # Calculate a columns of Observed CN gains/losses/and obs samples in the
  # cohort for each segment
  gainFreq <- function(x){length(x[x >= CN_high_cutoff])}
  gainFrequency <- aggregate(segmean ~ chromosome + start + end, data=x, gainFreq)$segmean
  
  lossFreq <- function(x){length(x[x <= CN_low_cutoff])}
  lossFrequency <- aggregate(segmean ~ chromosome + start + end, data=x, lossFreq)$segmean
  
  x <- aggregate(segmean ~ chromosome + start + end, data=x, length)
  colnames(x)[which(colnames(x) %in% "segmean")] <- "sampleFrequency"
  x$gainFrequency <- gainFrequency
  x$lossFrequency <- lossFrequency
  
  # check for coordinate space, if any widths are 1 it might indicate a problem
  if(max(x$sampleFrequency) > length(samples)){
    memo <- paste0("Detected additional sample rows after disjoin operation",
                   " typically this indicates coordinates are 0-based, please convert",
                   " coordinates to 1-base for accurate results")
    warning(memo)
  }
  
  # Calculate the proportion
  x$gainProportion <- x$gainFrequency/length(samples)
  x$lossProportion <- x$lossFrequency/length(samples)
  
  # get the dummy data for plot boundaries
  preloaded <- c("hg38", "hg19", "mm10", "mm9", "rn5")
  if(any(genome == preloaded)){
    message("genome specified is preloaded, retrieving data...")
    UCSC_Chr_pos <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == genome,]
    UCSC_Chr_pos <- multi_chrBound(UCSC_Chr_pos)
  } else {
    # Obtain data for UCSC genome and extract relevant columns
    memo <- paste0("attempting to query UCSC mySQL database for chromosome",
                   " positions")
    message(memo)
    cyto_data <- suppressWarnings(multi_cytobandRet(genome))
    UCSC_Chr_pos <- multi_chrBound(cyto_data)        
  }     
  
  # check that the dummy data has a size
  if(nrow(UCSC_Chr_pos) < 1)
  {
    memo <- paste0("did not recognize genome ", genome,
                   ", plotting provided data and ignoring chromosome ",
                   "boundaries! Output could be decieving!")
    warning(memo)
  }
  
  dummy_data <- lapply(unique(x$sample),
                       function(sample, chr_pos) cbind(chr_pos, sample),
                       UCSC_Chr_pos)
  dummy_data <- do.call("rbind", dummy_data)
  chr_order <- gtools::mixedsort(unique(dummy_data$chromosome))
  dummy_data$chromosome <- factor(dummy_data$chromosome, levels=chr_order)
  
  # select chromosomes to plot
  if(!is.null(plotChr)){
    if(any(!plotChr %in% dummy_data$chromosome)) {
      missingChr <- plotChr[!plotChr %in% dummy_data$chromosome]
      plotChr <- plotChr[!plotChr %in% missingChr]
      memo <- paste0("The following chromosomes: ", toString(missingChr),
                     ", could not be found! Valid chromosomes are: ",
                     toString(unique(dummy_data$chromosome)))
      warning(memo)
    }
    
    dummy_data <- dummy_data[dummy_data$chromosome %in% plotChr,]
    dummy_data$chromosome <- factor(dummy_data$chromosome, levels=plotChr)
    x <- x[x$chromosome %in% plotChr,]
    x$chromosome <- factor(x$chromosome, levels=plotChr)
  }
  res = list(x, dummy_data)
  return(res)
}

### MAIN ##########################################################################################
basedir <- "/oak/stanford/groups/ccurtis2/users/"
datadir <- paste0(basedir,"lisem/bc_landscape/github/")
outdir <- paste0(basedir,"lisem/bc_landscape/github/")

### Load chromosomes size
chr_sizes = read.table("/oak/stanford/groups/ccurtis2/users/srinivap/IC_classifiers/Data/General/hg38.chrom.sizes.txt", header = F, sep = "\t")

### Load DCIS data
## Load IC calls
labels_dcis = read.table(paste0(basedir,"khoulaha/BreastLandscape/clinical/DCIS_clinical.csv"), header = T, sep=",")
labels_dcis$ic10 = labels_dcis$IC10_RNA
labels_dcis$ic10[which(!is.na(labels_dcis$ic10))] = paste0("ic", labels_dcis$ic10[which(!is.na(labels_dcis$ic10))])

labels_dcis$group = NA
labels_dcis$group[which(labels_dcis$ic10 == "ic5")] = "HER2+"
labels_dcis$group[which(labels_dcis$ic10 %in% c("ic1", "ic2", "ic6", "ic9"))] = "ER+ High"
labels_dcis$group[which(labels_dcis$ic10 %in% c("ic3", "ic7", "ic8") | (labels_dcis$ic10 == "ic4" & labels_dcis$ER_RNA == "+"))] = "ER+ Typical"
labels_dcis$group[which(labels_dcis$ic10 == "ic10" | (labels_dcis$ic10 == "ic4" & labels_dcis$ER_RNA == "-"))] = "TNBC"

## Load CN
cn_segment <- data.frame()
for(f in list.files("/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/QDNAseq/hg38_50kb/")){
  print(f)
  if(paste0(f,".recal.seg") %in% list.files(paste0("/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/QDNAseq/hg38_50kb/",f,"/"))){
    f = paste0("/oak/stanford/groups/ccurtis2/users/azizk/data/htan-pca/QDNAseq/hg38_50kb/",f,"/",f,".recal.seg")
    cn_segment_tmp <- read.delim(f)
    cn_segment <- rbind(cn_segment_tmp, cn_segment)
    rm(cn_segment_tmp)
  }
}
cn_segment$SAMPLE_NAME <- gsub("-","_", cn_segment$SAMPLE_NAME)

## Matching IDs
labels_dcis = labels_dcis[which(!is.na(labels_dcis$ic10)),]
labels_dcis$Sample_ID[which(grepl("DCIS",labels_dcis$Sample_ID))] <- paste0(sapply(strsplit(labels_dcis$Sample_ID[which(grepl("DCIS",labels_dcis$Sample_ID))],"_"),"[[",1),"_DCIS")
cn_segment$SAMPLE_NAME[which(grepl("DCIS",cn_segment$SAMPLE_NAME) & !startsWith(cn_segment$SAMPLE_NAME,"H_SL"))] <- paste0(sapply(strsplit(cn_segment$SAMPLE_NAME[which(grepl("DCIS",cn_segment$SAMPLE_NAME) & !startsWith(cn_segment$SAMPLE_NAME,"H_SL"))],"_"),"[[",1),"_DCIS")

table(labels_dcis$Sample_ID %in% gsub(".recal","",cn_segment$SAMPLE_NAME))

labels_dcis$Sample_ID[which(!labels_dcis$Sample_ID %in% gsub(".recal","",cn_segment$SAMPLE_NAME))]
unique(cn_segment$SAMPLE_NAME[which(!gsub(".recal","",cn_segment$SAMPLE_NAME) %in% labels_dcis$Sample_ID)])

ids <- labels_dcis$Sample_ID[which(!labels_dcis$Sample_ID %in% gsub(".recal","",cn_segment$SAMPLE_NAME))][which(sapply(sapply(strsplit(labels_dcis$Sample_ID[which(!labels_dcis$Sample_ID %in% gsub(".recal","",cn_segment$SAMPLE_NAME))],"_"),"[[",1), function(x) any(grepl(x,cn_segment$SAMPLE_NAME))))]

## Add IC calls to CN segments
cn_segment$group <- sapply(gsub(".recal","",cn_segment$SAMPLE_NAME), function(x) if(x %in% labels_dcis$Sample_ID){labels_dcis$group[which(labels_dcis$Sample_ID == x)]}else{NA})
table(cn_segment$group, useNA = "always")

table_cna <- cn_segment
length(unique(table_cna$SAMPLE_NAME[which(table_cna$SAMPLE_NAME %in% labels_dcis$Sample_ID[which(labels_dcis$IC10_RNA == 2)])]))

## Generate plot datasets
groups = list("ER+ High" = unique(table_cna$SAMPLE_NAME[which(table_cna$group == "ER+ High")]), 
              "ER+ Typical" = unique(table_cna$SAMPLE_NAME[which(table_cna$group == "ER+ Typical")]), 
              "HER2+" = unique(table_cna$SAMPLE_NAME[which(table_cna$group == "HER2+")]), 
              "TNBC" = unique(table_cna$SAMPLE_NAME[which(table_cna$group == "TNBC")]))

for(g in 1:length(groups)){
  group = names(groups)[g]
  
  ids <- groups[[group]]
  
  table_cna.cncf <- table_cna %>% dplyr::filter(SAMPLE_NAME %in% ids) %>% dplyr::select(chromosome=CHROMOSOME, start=START, end=STOP, segmean=LOG2_RATIO_MEAN, sample=SAMPLE_NAME)
  table_cna.cncf$segmean <- (2**(table_cna.cncf$segmean))*2 
  res <- cnFreq_LM(table_cna.cncf, CN_low_cutoff = 2*(2^-0.3), CN_high_cutoff = 2*(2^0.3), genome = 'hg38')
  
  data.acnv <- res[[1]][c("chromosome","start","end","lossProportion")]
  colnames(data.acnv)[1] <- "V2"
  data.acnv$V1 <- "<2"
  data.acnv = data.acnv[which(data.acnv$lossProportion > 0),]
  data.acnv$lossProportion = -data.acnv$lossProportion
  colnames(data.acnv)[4] <- "V5"
  data.acnv$V3 <- apply(data.acnv[,c("start","end")],1,mean)
  data.acnv$V4 <- sapply(1:nrow(data.acnv), function(i) data.acnv$end[i] - data.acnv$start[i])
  data.acnv = data.acnv[,which(colnames(data.acnv) %in% paste0("V",1:5))]
  
  data.acnv_bis <- res[[1]][c("chromosome","start","end","gainProportion")]
  colnames(data.acnv_bis)[1] <- "V2"
  data.acnv_bis$V1 <- ">2"
  data.acnv_bis = data.acnv_bis[which(data.acnv_bis$gainProportion > 0),]
  colnames(data.acnv_bis)[4] <- "V5"
  data.acnv_bis$V3 <- apply(data.acnv_bis[,c("start","end")],1,mean)
  data.acnv_bis$V4 <- sapply(1:nrow(data.acnv_bis), function(i) data.acnv_bis$end[i] - data.acnv_bis$start[i])
  data.acnv_bis = data.acnv_bis[,which(colnames(data.acnv_bis) %in% paste0("V",1:5))]
  
  data.acnv <- rbind(data.acnv, data.acnv_bis)
  rm(data.acnv_bis)
  
  data.acnv = data.acnv[which(data.acnv$V2 %in% paste0("chr",1:22)),]
  data.acnv$V2 <- factor(data.acnv$V2, levels = paste0("chr",1:22))
  
  subdat1 = subset(data.acnv, V1==">2")
  subdat2 = subset(data.acnv, V1=="<2")
  subdat2$V1 <- factor(subdat2$V1, levels = levels(factor(subdat2$V1)))
  subdat1$V1 <- factor(subdat1$V1, levels = rev(levels(factor(subdat1$V1))))
  centro = data.frame(pos=c(123400000, 93900000, 90900000, 50000000, 48750000, 60550000, 60100000, 45200000, 43850000, 39800000, 53400000, 35500000, 17700000, 17150000, 19000000, 36850000, 25050000, 18450000, 26150000, 28050000, 11950000, 15550000),
                      V2=paste0("chr",c(1:22)))
  
  cols <- c(">2" = "darkorange1", "<2" = "steelblue4")
  if ("expand" %in% names(formals(coord_cartesian))) {
    yLim <- coord_cartesian(ylim=c(-100,100), expand = FALSE)
  } else{
    yLim <- coord_cartesian(ylim=c(-100,100))
  }
  
  centro$V2 <- factor(centro$V2, levels = levels(data.acnv$V2))
  
  save(data.acnv, subdat1, subdat2, file = paste0(datadir, 'CN_SVB_', group, '_DCIS.RData'))
}

### SAVE ##########################################################################################
## Figure----
for(g in 1:length(groups)){
  group = names(groups)[g]
  
  rm(data.acnv, subdat1, subdat2)
  load(file = paste0(datadir, 'CN_SVB_', group, '_DCIS.RData'))
  
  chr.labels <- gsub('chr','',levels(data.acnv$V2))
  names(chr.labels) <- levels(data.acnv$V2)
  
  tot = length(unique(groups[[group]]))
  print(paste0("Proportion in DCIS (n = ",tot, ")"))
  
  p <- ggplot(data=data.acnv) + geom_hline(yintercept=-25, colour="white", size=0.5) +
    facet_grid(~V2, space="free_x", scales="free_x", labeller=labeller(V2 = chr.labels)) +
    theme_LM +
    theme(axis.text.x=element_blank(), axis.ticks=element_blank(),
          axis.title.x=element_blank(),text = element_text(size=15),
          axis.text.y = element_text(size=15), legend.text = element_text(size=10),
          legend.position = "bottom", legend.box = "horizontal",
          panel.background = element_rect(fill = "white", colour = "grey", color = "grey"),
          panel.grid.major = element_line(color="grey"),
          panel.grid.minor = element_line(color="grey"),
          strip.background = element_rect(fill = "#8fdab0", colour = "transparent")) +
    scale_fill_manual(name = "", values = c(">2" = col_groups[which(names(col_groups) == group)][[1]][1], "<2" = col_groups[which(names(col_groups) == group)][[1]][2]), labels = c(">2" = "AMP", "<2" = "DEL")) +
    geom_bar(data=subdat1 ,aes(x=V3, y=V5, fill=factor(V1), width=V4), stat="identity", alpha = 0.5) +
    geom_bar(data=subdat2 ,aes(x=V3, y=V5, fill=factor(V1), width=V4), stat="identity", alpha = 0.5) +
    ylim(min(subdat2$V5),max(subdat1$V5)) +
    ylab("") +
    guides(fill=NULL) + 
    geom_point(aes(x=pos, y=0), centro , size=1.5) +
    geom_hline(yintercept=0, colour="black", size=0.5) + 
    scale_colour_manual(values = c(col_groups[which(names(col_groups) == group)][[1]][1]), labels=c(""), name="SV burden")
  
  png(filename = file.path(outdir, paste0('CN_SVB_', group, '_DCIS.png')), res = 300, width = 16.13, height = 3.71/1.5, units = 'in')
  print(p)
  dev.off()
}

## SourceData---- 
sourcetable <- data.frame()
for(g in 1:length(groups)){
  group = names(groups)[g]
  rm(data.acnv)
  
  load(file = paste0(datadir, 'CN_SVB_', group, '_DCIS.RData'))
  colnames(data.acnv) <- c("Chromosome","Proportion","CN","Segment_location","Segment_length")
  data.acnv[,"Subgroup"] <- group
  data.acnv[,"Stage"] <- "DCIS"
  sourcetable <- rbind(sourcetable, data.acnv)
}
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/github/ExtendedData2a_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
