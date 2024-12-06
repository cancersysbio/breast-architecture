### CREATE FIGURE 3e #############################################################################
# creates figure 3e
# figure 3e provides the proportion of tumors within each ER+ High-risk and HER2+ IC subtype that harbor cyclic, complex non-cyclic or linear amplification in IC-specific oncogenes.

### PREAMBLE #####################################################################################
library(yaml)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(tidyr)
require(patchwork)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
col.ic = c("ic1" = "#ff5500ff", "ic2" = "#00ee77ff", "ic3" = "#cd3279ff",
           "ic4" = "#bbe0e2ff", "ic4+" = "#00c4cdff", "ic4-" = "#b8d0d5ff",
           "ic5" = "#8b0100ff", "ic6" = "#fffe41ff", "ic7" = "#0000cdff",
           "ic8" = "#feaa00ff","ic9" = "#ee82edff", "ic10" = "#7c26ccff",
           "ic19" = "#ff192c")

### MAIN ##########################################################################################
basedir <- "/oak/stanford/groups/ccurtis2/users/"
outdir <- paste0(basedir,"lisem/bc_landscape/publication/figures/")

### Load GENCODE annotations
ref.39 = rtracklayer::import("/oak/stanford/groups/ccurtis2/isabl/assemblies/GRCh37/STARIndex/gencode.v39lift37.annotation.gtf", format="gtf")
ref.39 = as.data.table(as.data.frame(ref.39))
ref.39.gene = ref.39[which(ref.39$type == "gene"),]
rm(ref.39)

### Load COSMIC genes (from: https://cancer.sanger.ac.uk/census)
cosmic = read.delim(paste0(basedir,"lisem/general/Census_allWed Aug  2 00_50_58 2023.tsv"))

### Load Rueda et al. genes
rueda_ic1 = c("HSF5", "PPM1E", "PRR11", "DHX40", "TUBD1", "RPS6KB1", "CA4", "C17orf64", "BCAS3", "TBX2", "BRIP1", "TBC1D3P2")
rueda_ic2 = c("FGF3", "CCND1", "CTTN", "CLPB", "P2RY2", "UCP2", "CHRDL2", "MAP6", "OMP", "PAK1", "RSF1", "NARS2")
rueda_ic6 = c("ZNF703", "EIF4EBP1", "LETM2", "STAR", "FGFR1")
rueda_ic9 = c("FBXO32", "SQLE", "LINC00861", "PCAT1", "MYC", "LINC00977", "MIR5194", "ADCY8")

rueda_ic = list(ic1 = rueda_ic1, ic2=rueda_ic2, ic6=rueda_ic6, ic9=rueda_ic9, ic5=NA, ic19=NA)

### Load megatables: without donor-duplicates
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")

### Load considered genes
load(file = paste0(basedir,"lisem/bc_landscape/github/genes.ic.RData"))

mytitles = c(ic1="IC1 - 17q23 AMP",ic2="IC2 - 11q13 AMP",ic6="IC6 - 8p12 AMP",ic9="IC9 - 8q24 AMP",ic5="IC5 - 17q12 AMP",ic19="IC1/IC9 - 20q13 AMP")

sourcetable = data.frame()
freq = c()

plot_list = list()
for(g in 1:length(genes.ic)){
  
  ic = names(genes.ic)[g]
  mytitle = mytitles[which(names(mytitles) == ic)][[1]]
  print(mytitle)
  genes = genes.ic[[g]]
  chr = as.character(ref.39.gene$seqnames[which(ref.39.gene$gene_name == genes[1])])
  rueda = rueda_ic[which(names(rueda_ic) == ic)][[1]]
  
  ### Load amplicon annotations
  ## Load Amplicon segments
  amplicon = read.table(paste0(basedir, "khoulaha/BreastLandscape/ecDNA/amplicon_segments_project7_all_v2.tsv"), sep="\t", header = T) %>%
    dplyr::bind_rows(read.table(paste0(basedir, "khoulaha/BreastLandscape/ecDNA/amplicon_segments_project17_all_v2.tsv"), sep="\t", header = T)) %>%
    dplyr::filter(sample %in% mega_1$Sample) # Primary only
  
  any(grepl(paste0(c("CURTIS_H000535","CURTIS_H004681"), collapse = "|"),amplicon$sample))
  
  ## Map GENCODE genes with Amplicon segments
  query.gr <- GRanges(
    seqnames = ref.39.gene$seqnames[which(ref.39.gene$gene_name %in% genes)],
    ranges = IRanges(start = ref.39.gene$start[which(ref.39.gene$gene_name %in% genes)], end = ref.39.gene$end[which(ref.39.gene$gene_name %in% genes)]),
    strand = ref.39.gene$strand[which(ref.39.gene$gene_name %in% genes)]
  )
  
  subject.gr <- GRanges(
    seqnames = paste0("chr", amplicon$chr),
    ranges = IRanges(start = amplicon$start, end = amplicon$end)
  )
  
  hits <- findOverlaps(query.gr, subject.gr,
                       ignore.strand=F,
                       type="any",
                       select="all")
  
  amplicon$genes_ic = NA
  for(i in unique(hits@to)){
    amplicon$genes_ic[i] = paste0(ref.39.gene$gene_name[which(ref.39.gene$gene_name %in% genes)][hits@from[which(hits@to == i)]], collapse = ",")
  }
  
  ### Load copy number annotations
  ## Load Gene-level CN
  gene_cn = read.table(paste0(basedir, "lisem/ICGC/01_gene_level.tsv"), sep="\t", header = T) %>%
    dplyr::bind_rows(read.table(paste0(basedir, "lisem/Hartwig/eniclust/01_gene_level.txt"), sep="\t", header = T)) %>%
    dplyr::mutate(Sample = sapply(strsplit(Sample, "_vs_"),"[[",1)) %>%
    dplyr::filter(Sample %in% mega_1$Sample) # Primary only        
  
  ## Load CN segments
  segment_cn = read.table(paste0(basedir, "lisem/ICGC/02_segment.tsv"), sep="\t", header = T) %>%
    dplyr::bind_rows(read.table(paste0(basedir, "lisem/Hartwig/eniclust/02_segment.txt"), sep="\t", header = T)) %>%
    dplyr::filter(ID %in% mega_1$Sample) # Primary only   
  
  ## Map GENCODE genes with Amplicon segments
  subject.gr <- GRanges(
    seqnames = paste0("chr", segment_cn$chrom[which(segment_cn$chrom == gsub("chr","",chr))]),
    ranges = IRanges(start = segment_cn$loc.start[which(segment_cn$chrom == gsub("chr","",chr))], 
                     end = segment_cn$loc.end[which(segment_cn$chrom == gsub("chr","",chr))])
  )
  
  hits <- findOverlaps(query.gr, subject.gr,
                       ignore.strand=F,
                       type="any",
                       select="all")
  
  segment_cn$genes_ic = NA
  for(i in unique(hits@to)){
    segment_cn$genes_ic[i] = paste0(ref.39.gene$gene_name[which(ref.39.gene$gene_name %in% genes)][hits@from[which(hits@to == i)]], collapse = ",")
  }
  
  ### Plotting
  if(ic == "ic19"){
    mysamples = mega_1$Sample[which(mega_1$ENiClust %in% c("ic1","ic9"))]
  }else{
    mysamples = mega_1$Sample[which(mega_1$ENiClust == ic)]
  }
  
  stopifnot(all(mysamples %in% gene_cn$Sample))
  data <- gene_cn %>% dplyr::filter(Sample %in% mysamples) %>%
    dplyr::select(c(Sample,genes[genes %in% colnames(gene_cn)])) %>%
    dplyr::mutate_at(vars(-c(Sample)), ~replace(., .>2, "YES")) %>%
    pivot_longer(!Sample, names_to = "gene", values_to = "amp") %>%
    dplyr::rowwise() %>% dplyr::mutate(amp = ifelse(!is.na(amp) & amp != "YES","NO",amp)) %>%
    dplyr::rowwise() %>% dplyr::mutate(segment = ifelse(gene %in% unlist(strsplit(segment_cn$genes_ic[which(segment_cn$tcn.em > 2 & segment_cn$sample == Sample & segment_cn$chrom == chr)],",")),"YES","NO")) %>%
    dplyr::rowwise() %>% dplyr::mutate(ecDNA = ifelse(gene %in% unlist(strsplit(amplicon$genes_ic[which(amplicon$amplicon_type == "Cyclic" & amplicon$ecDNA == "Positive" & amplicon$sample == Sample & paste0("chr",amplicon$chr) == chr)],",")),"YES","NO")) %>%
    dplyr::rowwise() %>% dplyr::mutate(complex = ifelse(gene %in% unlist(strsplit(amplicon$genes_ic[which(amplicon$amplicon_type == "Complex non-cyclic" & amplicon$sample == Sample & paste0("chr",amplicon$chr) == chr)],",")),"YES","NO")) %>%
    dplyr::rowwise() %>% dplyr::mutate(linear = ifelse(gene %in% unlist(strsplit(amplicon$genes_ic[which(amplicon$amplicon_type == "Linear amplification" & amplicon$sample == Sample & paste0("chr",amplicon$chr) == chr)],",")),"YES","NO")) %>%
    dplyr::rowwise() %>% dplyr::mutate(type = ifelse(amp == "NO", "NO", ifelse(amp == "YES" & ecDNA == "YES", "ecDNA-AMP", ifelse(amp == "YES" & ecDNA == "NO" & complex == "YES", "Complex-AMP", ifelse(amp == "YES" & ecDNA == "NO" & complex == "NO", "AMP", NA))))) %>%
    dplyr::rowwise() %>% dplyr::mutate(label = ifelse(gene %in% cosmic$Gene.Symbol[which(cosmic$Role.in.Cancer %in% c("oncogene","oncogene, fusion"))],paste0("<span style = 'color: red;'>",gene,"</span>"), paste0("<span style = 'color: black;'>",gene,"</span>")))
  
  stopifnot(all(is.na(data$segment[which((is.na(data$type) | data$type == "NO") & !data$Sample %in% data$Sample[which(data$type == "YES")])]) | data$segment[which((is.na(data$type) | data$type == "NO") & !data$Sample %in% data$Sample[which(data$type == "YES")])] == "NO"))
  
  gene_levels = data %>% 
    dplyr::filter(type == "NO") %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::arrange(n) %>% 
    dplyr::select(gene) %>% pull()
  
  gene_levels = c(unique(as.character(data$gene)[which(!as.character(data$gene) %in% gene_levels)]), as.character(gene_levels))
  
  data$type = factor(data$type, levels = c("NO","AMP","Complex-AMP","ecDNA-AMP"))
  data$gene = factor(data$gene, levels = c(rueda, genes[which(!genes %in% rueda)]))
  data$label = factor(data$label, levels = sapply(gene_levels, function(x) unique(data$label[which(data$gene == x)])))
  
  p <- data %>%
    dplyr::filter(!is.na(type)) %>% 
    ggplot(aes(x=label, fill=type)) + 
    geom_bar(position="fill", col = "black", linewidth = 0.5) +
    scale_fill_manual(name = "", values = c("ecDNA-AMP"=as.character(col.ic[which(names(col.ic) == ic)]), "Complex-AMP"=as.character(alpha(col.ic[which(names(col.ic) == ic)],0.8)),
                                            "AMP"=as.character(alpha(col.ic[which(names(col.ic) == ic)],0.6)), "NO"=as.character(alpha(col.ic[which(names(col.ic) == ic)],0.2))), guide = "none") +
    xlab(mytitle) +
    ylab(ifelse(ic == "ic1","Frequency","")) +
    theme_bw() + ggtitle(paste0("Primary " ,ifelse(ic == "ic19", "IC1/IC9", gsub("ic","IC",ic)), " (n=",length(unique(data$Sample[which(!is.na(data$type))])),")")) +
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1, size = 5*4),
          axis.text.y = element_text(size = 5*4),
          axis.title = element_text(size = 6*4),
          title = element_text(size = 7*4))
  
  plot_list[[length(plot_list) + 1]] <- p
  
  freq_ecDNA <- max(sapply(unique(data$gene), function(g) length(unique(data$Sample[which(data$gene == g & data$type == "ecDNA-AMP")]))))
  print(freq_ecDNA/length(unique(data$Sample[which(!is.na(data$type))])))
  freq = c(freq, freq_ecDNA/length(unique(data$Sample[which(!is.na(data$type))])))
  
  data$ic = ic
  sourcetable = rbind(sourcetable, data)
}

### SAVE ##########################################################################################
## Figure----
mysubtitles = c(ic1 = "Primary IC1 (n=40)", ic2 = "Primary IC2 (n=22)", ic6 = "Primary IC6 (n=25)", 
                ic9 = "Primary IC9 (n=50)", ic5 = "Primary IC5 (n=70)", "ic19" = "Primary IC1/IC9 (n=90)")

ind = 3
p <- sourcetable %>%
  dplyr::rowwise() %>% dplyr::mutate(mytitle = mytitles[which(names(mytitles) == ic)]) %>%
  dplyr::rowwise() %>% dplyr::mutate(mysubtitle = mysubtitles[which(names(mysubtitles) == ic)]) %>%
  dplyr::filter(!is.na(type)) %>% 
  ggplot(aes(x=label, alpha=factor(type, levels = rev(c("ecDNA-AMP","Complex-AMP","AMP","NO"))), fill=ic)) + 
  geom_bar(position="fill", col = "black", linewidth = 0.5) +
  scale_fill_manual(name = "", values = col.ic, guide = "none") +
  scale_alpha_manual(name = "", values = c("ecDNA-AMP"=1, "Complex-AMP"=0.8,"AMP"=0.6, "NO"=0.2), guide = "none") +
  xlab("") +
  ylab("Frequency") +
  #facet_grid(. ~ factor(mytitle, levels = c("IC1 - 17q23 AMP","IC2 - 11q13 AMP","IC6 - 8p12 AMP","IC9 - 8q24 AMP","IC5 - 17q12 AMP","IC1/IC9 - 20q13 AMP")), scales = "free_x") +
  facet_grid(. ~ factor(mysubtitle, levels = c("Primary IC1 (n=40)", "Primary IC2 (n=22)", "Primary IC6 (n=25)", "Primary IC9 (n=50)", "Primary IC5 (n=70)", "Primary IC1/IC9 (n=90)")), scales = "free_x", space = "free_x") +
  theme_bw() + #ggtitle(paste0("Primary " ,ifelse(ic == "ic19", "IC1/IC9", gsub("ic","IC",ic)), " (n=",length(unique(data$Sample[which(!is.na(data$type))])),")")) +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1, size = 5*ind),
        axis.text.y = element_text(size = 5*ind),
        axis.title = element_text(size = 6*ind),
        strip.text = element_text(size = 6*ind), 
        strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"))

outdir = "/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/github/"
pdf(paste0(outdir, 'Figure3e.pdf'), width = 11.20*1.7, height = 5.78/1.5)
p
dev.off()

## SourceData---- 
sourcetable = sourcetable[,c(1,10,2:4,7,6,5,8)]
colnames(sourcetable) = c("Sample","ENiClust", "Gene", "AMP","Segment", "linear","complex","ecDNA","Summary")
sourcetable$ENiClust[which(sourcetable$ENiClust == "ic19")] = sapply(sourcetable$Sample[which(sourcetable$ENiClust == "ic19")], function(x) paste0("IC1 or IC9: ",mega_1$ENiClust[which(mega_1$Sample == x)]))
write.table(sourcetable, paste0(basedir, "lisem/bc_landscape/github/Figure3e_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

## For main text----
amplicon = read.table(paste0(basedir, "khoulaha/BreastLandscape/ecDNA/amplicon_segments_project7_all_v2.tsv"), sep="\t", header = T) %>%
  dplyr::bind_rows(read.table(paste0(basedir, "khoulaha/BreastLandscape/ecDNA/amplicon_segments_project17_all_v2.tsv"), sep="\t", header = T)) %>%
  dplyr::filter(sample %in% mega_1$Sample)

mega_1 %>%
  dplyr::filter(ENiClust == "ic9") %>%
  dplyr::rowwise() %>% dplyr::mutate(ecDNA = Sample %in% amplicon$sample[which(amplicon$amplicon_type == "Cyclic" & amplicon$ecDNA == "Positive")]) %>%
  dplyr::group_by(ecDNA) %>%
  dplyr::summarize(n = n_distinct(Sample))
21/(21+29)
