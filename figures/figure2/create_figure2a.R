### CREATE FIGURE 2a CN/SV profiles #############################################################################
# creates figure 2a CN/SV profiles
# figure 2a provides the IC group-level copy number profile with structural variant burden as overlay. 

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

col_groups = list("ER+ High" = c("#ff6a00", "#ffb07c"), "ER+ Typical" = c("#7aa6c2",alpha("#7aa6c2",0.7)), 
                  "HER2+" = c("#8b0100ff","#dd6867ff"), 
                  "IC10" = c("#6300a1","#bc8ad9"), "IC4 ER-" = c("#8e7abcff","#d1c9e4ff"))

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
      #x <- sapply(x, function(x_chr) cnFreq_disjoin(x_chr))
      
      #for(x_chr in x){
      #  print(unique(as.character(x_chr$chromosome)))
      #  cnFreq_disjoin(x_chr)
      #}
      
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
  #disJoint_x <- rep(GRangesList(disJoint_x), lengths(GRangesList(revmap))) 
  #disJoint_x <- rep(as(disJoint_x, "SimpleList"), lengths(revmap))
  #disJoint_x <- rep(as(disJoint_x, "SimpleList"), lengths(as(revmap, "SimpleList")))
  
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
} # with too large dataset in x this error occurs: https://github.com/griffithlab/GenVisR/issues/312

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
datadir <- paste0(basedir,"lisem/bc_landscape/submission/data/")
outdir <- paste0(basedir,"lisem/bc_landscape/submission/figures/")

### Load megatables
## Primary/Metastasis
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])
mega$Sample_Type = NA
mega$Sample_Type[which(mega$Sample %in% mega_1$Sample)] = "Primary"
mega$Sample_Type[which(mega$Sample %in% mega_2$Sample)] = "Metastatic"

### Load chromosomes size
chr_sizes = read.table(paste0(basedir,"srinivap/IC_classifiers/Data/General/hg19.chrom.sizes.txt"), header = F, sep = "\t")

### CN profile
## Load Gene-level CN
table_cna = read.table(paste(basedir, "lisem/ICGC/02_segment.tsv", sep=""), sep="\t", header = T)
table_cna_bis = read.table(paste(basedir, "lisem/Hartwig/eniclust/02_segment.txt", sep=""), sep="\t", header = T) %>%
  dplyr::mutate(ID = gsub('_hisens','',ID))

table_cna = dplyr::bind_rows(table_cna, table_cna_bis)
rm(table_cna_bis)

table_cna = table_cna[which(table_cna$ID %in% mega$Sample),]
table_cna$group <- sapply(table_cna$ID, function(x) mega$group[which(mega$Sample == x)])
table(table_cna$group, useNA = "always")

### SVB profile
## Load SV tables
# data_pBC is a dataframe of SV breakpoints position instead of SV event positions (dataset.sv)
load(file = paste0(basedir,"lisem/bc_landscape/sv/data/data_pBC.RData")) # see /oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/sv/scripts/00_calculate_sv_burden.R
load(file = paste0(basedir,"lisem/bc_landscape/sv/data/data_mBC.RData")) # see /oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/sv/scripts/00_calculate_sv_burden.R
data_BC = rbind(data_pBC, data_mBC)
rm(data_pBC, data_mBC)

## Quantify SVB : count is the number of samples with at least one SV in the given region
resolution = 10

# Primary
groups = list("ER+ High" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Primary" & mega$group == "ER+ High")])], 
              "ER+ Typical" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Primary" & mega$group == "ER+ Typical")])], 
              "HER2+" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Primary" & mega$group == "HER2+")])], 
              "IC10" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Primary" & mega$group == "IC10")])],
              "IC4 ER-" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Primary" & mega$group == "IC4ER-")])])

df_merge = data.frame()
for(g in 1:length(groups)){
  
  group = groups[[g]]
  
  data = data_BC[which(data_BC$sample %in% group),]
  
  start_seq = sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) seq(from=1, to=chr_sizes$V2[which(chr_sizes$V1 == x)], by = 100000*resolution)[-length(seq(from=1, to=chr_sizes$V2[which(chr_sizes$V1 == x)], by = 100000*resolution))])
  end_seq = sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) (seq(from=1, to=chr_sizes$V2[which(chr_sizes$V1 == x)], by = 100000*resolution)-1)[-1])
  
  query.gr <- GRanges(
    seqnames = paste0("chr", data$chr),
    ranges = IRanges(start = data$start, end = data$end+1),
    strand = data$strand
  )
  
  subject.obj = as.data.table(data.frame(chr = unlist(sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) rep(x,length(start_seq[[x]])))),
                                         bp = unlist(end_seq)+1))
  
  subject.gr <- GRanges(
    seqnames = unlist(sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) rep(x,length(start_seq[[x]])))),
    ranges = IRanges(start = unlist(start_seq), end = unlist(end_seq))
  )
  
  hits <- findOverlaps(query.gr, subject.gr,
                       ignore.strand=F,
                       type="any",
                       select="all")
  
  hits <- hits[!duplicated(hits@from)]
  
  df = data.frame(chr = sapply(unique(hits@to), function(x) subject.obj$chr[x]),
                  count = sapply(unique(hits@to), function(x) length(unique(data$sample[unique(hits@from[which(hits@to == x)])]))),
                  bp = sapply(unique(hits@to), function(x) subject.obj$bp[x]), stringsAsFactors = F)
  
  df$chr = sapply(strsplit(df$chr, "chr"),"[[", 2)
  df$chr[which(df$chr == "X")] = "23"
  table(df$chr, useNA = "always")
  df$chr = as.numeric(df$chr)
  
  df$group = names(groups)[g] 
  df_merge = rbind(df_merge, df)
  
}
df_cum <- df_merge %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarise(max_bp = max(bp)) %>% 
  dplyr::mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  dplyr::select(chr, bp_add)

df_merge <- df_merge %>% 
  inner_join(df_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

# Metastatic
groups_mets = list("ER+ High" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Metastatic" & mega$group == "ER+ High")])], 
                   "ER+ Typical" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Metastatic" & mega$group == "ER+ Typical")])], 
                   "HER2+" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Metastatic" & mega$group == "HER2+")])], 
                   "IC10" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Metastatic" & mega$group == "IC10")])],
                   "IC4 ER-" = data_BC$sample[which(data_BC$sample %in% mega$Sample[which(mega$Sample_Type == "Metastatic" & mega$group == "IC4ER-")])])

resolution = 10

df_merge_mets = data.frame()
for(g in 1:length(groups_mets)){
  
  group = groups_mets[[g]]
  
  data = data_BC[which(data_BC$sample %in% group),]
  
  start_seq = sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) seq(from=1, to=chr_sizes$V2[which(chr_sizes$V1 == x)], by = 100000*resolution)[-length(seq(from=1, to=chr_sizes$V2[which(chr_sizes$V1 == x)], by = 100000*resolution))])
  end_seq = sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) (seq(from=1, to=chr_sizes$V2[which(chr_sizes$V1 == x)], by = 100000*resolution)-1)[-1])
  
  query.gr <- GRanges(
    seqnames = paste0("chr", data$chr),
    ranges = IRanges(start = data$start, end = data$end+1),
    strand = data$strand
  )
  
  subject.obj = as.data.table(data.frame(chr = unlist(sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) rep(x,length(start_seq[[x]])))),
                                         bp = unlist(end_seq)+1))
  
  subject.gr <- GRanges(
    seqnames = unlist(sapply(c(paste0("chr",unique(sort(as.numeric(data$chr)))), "chrX"), function(x) rep(x,length(start_seq[[x]])))),
    ranges = IRanges(start = unlist(start_seq), end = unlist(end_seq))
  )
  
  hits <- findOverlaps(query.gr, subject.gr,
                       ignore.strand=F,
                       type="any",
                       select="all")
  
  hits <- hits[!duplicated(hits@from)]
  
  df = data.frame(chr = sapply(unique(hits@to), function(x) subject.obj$chr[x]),
                  count = sapply(unique(hits@to), function(x) length(unique(data$sample[unique(hits@from[which(hits@to == x)])]))),
                  bp = sapply(unique(hits@to), function(x) subject.obj$bp[x]), stringsAsFactors = F)
  
  df$chr = sapply(strsplit(df$chr, "chr"),"[[", 2)
  df$chr[which(df$chr == "X")] = "23"
  table(df$chr, useNA = "always")
  df$chr = as.numeric(df$chr)
  
  df$group = names(groups_mets)[g] 
  df_merge_mets = rbind(df_merge_mets, df)
  
}
df_cum <- df_merge_mets %>% 
  dplyr::group_by(chr) %>% 
  dplyr::summarise(max_bp = max(bp)) %>% 
  dplyr::mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  dplyr::select(chr, bp_add)

df_merge_mets <- df_merge_mets %>% 
  inner_join(df_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

### Merged plot
# generate/save data to plot
for(g in 1:length(groups)){
  group = names(groups)[g]
  mygroup = ifelse(group == "IC4 ER-", "IC4ER-", group)
  
  for(stage in c("Primary","Metastatic")){
    
    ids <- mega$Sample[which(mega$group == mygroup & mega$Sample_Type == stage)]
    
    #table_cna.cncf <- table_cna %>% dplyr::filter(ID %in% ids) %>% dplyr::select(chromosome=chrom, start=loc.start, end=loc.end, segmean=tcn.em, sample=ID)
    table_cna.cncf <- table_cna %>% dplyr::filter(ID %in% ids, (tcn.em != 2 | lcn.em != 1)) %>% 
      dplyr::select(chromosome=chrom, start=loc.start, end=loc.end, segmean=tcn.em, sample=ID)
    res <- cnFreq_LM(table_cna.cncf, CN_low_cutoff = 2*(2^-0.3), CN_high_cutoff = 2*(2^0.3), genome = 'hg19')
    
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
    
    # add SVB
    
    my_df_merge <- list("Primary" = df_merge, "Metastatic" = df_merge_mets)
    my_df_merge = my_df_merge[[stage]]
    
    my_df_merge$V2 <- paste0("chr",my_df_merge$chr)
    my_df_merge$V2 <- factor(my_df_merge$V2, levels = levels(data.acnv$V2))
    
    chr.labels <- gsub('chr','',levels(data.acnv$V2))
    names(chr.labels) <- levels(data.acnv$V2)
    
    save(data.acnv, subdat1, subdat2, my_df_merge, file = paste0(datadir, 'CN_SVB_', group, '_', stage,'.RData'))
  }
}

### SAVE ##########################################################################################
## Figure----
for(g in 1:length(groups)){
  group = names(groups)[g]
  mygroup = ifelse(group == "IC4 ER-", "IC4ER-", group)
  
  for(stage in c("Primary","Metastatic")){
    rm(data.acnv, subdat1, subdat2, my_df_merge)
    load(file = paste0(datadir, 'CN_SVB_', group, '_', stage,'.RData'))
    
    centro = data.frame(pos=c(123400000, 93900000, 90900000, 50000000, 48750000, 60550000, 60100000, 45200000, 43850000, 39800000, 53400000, 35500000, 17700000, 17150000, 19000000, 36850000, 25050000, 18450000, 26150000, 28050000, 11950000, 15550000),
                        V2=paste0("chr",c(1:22)))
    centro$V2 <- factor(centro$V2, levels = levels(data.acnv$V2))
    
    chr.labels <- gsub('chr','',levels(data.acnv$V2))
    names(chr.labels) <- levels(data.acnv$V2)
    
    tot = list("Primary" = length(unique(groups[[group]])), "Metastatic" = length(unique(groups_mets[[group]])))
    tot = tot[[stage]]
    print(paste0("Proportion in ", stage," (n = ",tot, ")"))
    
    p <- ggplot(data=data.acnv) + geom_hline(yintercept=-25, colour="white", linewidth=0.5) +
      facet_grid(~V2, space="free_x", scales="free_x", labeller=labeller(V2 = chr.labels)) +
      theme_LM +
      theme(axis.text.x=element_blank(), axis.ticks=element_blank(),
            axis.title.x=element_blank(),text = element_text(size=15),
            axis.text.y = element_text(size=15), legend.text = element_text(size=10),
            #legend.position = "bottom", legend.box = "horizontal",
            legend.position = "none",
            panel.background = element_rect(fill = "white", colour = "grey", color = "grey"),
            panel.grid.major = element_line(color="grey"),
            panel.grid.minor = element_line(color="grey"),
            strip.background = element_rect(fill = ifelse(stage == "Primary","#9ebff1ff","#fe85cfff"), colour = "transparent")) +
      scale_fill_manual(name = "", values = c(">2" = col_groups[which(names(col_groups) == group)][[1]][1], "<2" = col_groups[which(names(col_groups) == group)][[1]][2]), labels = c(">2" = "AMP", "<2" = "DEL")) +
      geom_bar(data=subdat1 ,aes(x=V3, y=V5, fill=factor(V1), width=V4), stat="identity", alpha = 0.5) +
      geom_bar(data=subdat2 ,aes(x=V3, y=V5, fill=factor(V1), width=V4), stat="identity", alpha = 0.5) +
      ylim(min(subdat2$V5),ifelse(group %in% c("HER2+","IC10", "IC4 ER-"),1,max(subdat1$V5))) +
      ylab("") +
      guides(fill=NULL) + 
      geom_point(aes(x=pos, y=0), centro , size=1.5) +
      geom_hline(yintercept=0, colour="black", linewidth=0.5) + 
      geom_line(data=my_df_merge[which(my_df_merge$group == group & my_df_merge$chr %in% 1:22),] , aes(x=bp, y=count/tot, color = "black"), stat = "identity", linewidth = 0.6) + 
      scale_colour_manual(values = c(col_groups[which(names(col_groups) == group)][[1]][1]), labels=c(""), name="SV burden")
    
    png(filename = file.path(outdir, paste0('CN_SVB_', group, '_', stage,'.png')), res = 300, width = 16.13, height = 2.47 - (2.47*0.25), units = 'in')
    print(p)
    dev.off()
  }
}

## SourceData---- 
for(stage in c("Primary","Metastatic")){
  sourcetable <- data.frame()
  
  for(g in 1:length(groups)){
    group = names(groups)[g]
    mygroup = ifelse(group == "IC4 ER-", "IC4ER-", group)
    rm(data.acnv)
    
    load(file = paste0(datadir, 'CN_SVB_', group, '_', stage,'.RData'))
    colnames(data.acnv) <- c("Chromosome","Proportion","CN","Segment_location","Segment_length")
    data.acnv[,"Subgroup"] <- mygroup
    data.acnv[,"Stage"] <- stage
    sourcetable <- rbind(sourcetable, data.acnv)
  }
  write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/github/Figure2a_",stage,"_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
}



