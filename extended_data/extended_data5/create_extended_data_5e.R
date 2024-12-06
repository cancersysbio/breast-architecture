### CREATE EXTENDED DATA FIGURE 5e #############################################################################
# creates extended data figure 5e
# extended data figure 5e provides the proportion of metastatic tumors within each ER+ High-risk and HER2+ IC subtype that harbor cyclic, complex non-cyclic or linear amplification in the IC-specific oncogenes. Genes indicated in red are reported as oncogenes in COSMIC. The number of tumors within each subtype are indicated at the top of each subpanel.

### PREAMBLE #####################################################################################
library(yaml)
library(dplyr)
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
           "ic8" = "#feaa00ff","ic9" = "#ee82edff", "ic10" = "#7c26ccff","ic19" = "#ff192c")

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

### Load megatables
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")

### Load considered genes: see: /oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/scripts/Figure4C.R
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
    dplyr::filter(sample %in% mega_2$Sample) # Metastatic only
  
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
    dplyr::filter(Sample %in% mega_2$Sample) # Metastatic only        
  
  ## Load CN segments
  segment_cn = read.table(paste0(basedir, "lisem/ICGC/02_segment.tsv"), sep="\t", header = T) %>%
    dplyr::bind_rows(read.table(paste0(basedir, "lisem/Hartwig/eniclust/02_segment.txt"), sep="\t", header = T)) %>%
    dplyr::filter(ID %in% mega_2$Sample) # Metastatic only   
  
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
    mysamples = mega_2$Sample[which(mega_2$ENiClust %in% c("ic1","ic9"))]
  }else{
    mysamples = mega_2$Sample[which(mega_2$ENiClust == ic)]
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
    ylab("Frequency") +
    theme_bw() + ggtitle(paste0("Metastatic " ,ifelse(ic == "ic19", "IC1/IC9", gsub("ic","IC",ic)), " (n=",length(unique(data$Sample[which(!is.na(data$type))])),")")) +
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1))
  
  plot_list[[length(plot_list) + 1]] <- p
  
  freq_ecDNA <- max(sapply(unique(data$gene), function(g) length(unique(data$Sample[which(data$gene == g & data$type == "ecDNA-AMP")]))))
  print(freq_ecDNA/length(unique(data$Sample[which(!is.na(data$type))])))
  freq = c(freq, freq_ecDNA/length(unique(data$Sample[which(!is.na(data$type))])))
  
  data$ic = ic
  sourcetable = rbind(sourcetable, data)
}

### SAVE ##########################################################################################
## Figure----
outdir = "/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/github/"
#png(filename = file.path(outdir, 'ExtendedData5e.png'), res = 300, width = 11.20, height = 5.78, units = 'in')
pdf(paste0(outdir, 'ExtendedData5e.pdf'), width = 11.20, height = 5.78)
(plot_list[[1]] + plot_list[[2]] + plot_list[[3]]) / (plot_list[[4]] + plot_list[[5]] + plot_list[[6]])
dev.off()

sourcetable = sourcetable[,c(1,10,2:4,7,6,5,8)]
colnames(sourcetable) = c("Sample","ENiClust", "Gene", "AMP","Segment", "linear","complex","ecDNA","Summary")
sourcetable$ENiClust[which(sourcetable$ENiClust == "ic19")] = sapply(sourcetable$Sample[which(sourcetable$ENiClust == "ic19")], function(x) paste0("IC1 or IC9: ",mega_2$ENiClust[which(mega_2$Sample == x)]))
write.table(sourcetable, paste0(basedir, "lisem/bc_landscape/github/ExtendedData5e_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
