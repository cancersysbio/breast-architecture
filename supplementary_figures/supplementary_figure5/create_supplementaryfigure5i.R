### CREATE SUPPLEMENTARY FIGURE 5I #############################################################################
# Co-amplification of MYC and PVT1

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(GenomicRanges)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

# Data 
seg_primary_path <- "/path/to/merged/primary/CN/segments.txt"
seg_metastatic_path <- "/path/to/merged/metastatic/CN/segments.txt"

seg_primary <- read.delim(seg_primary_path)
seg_metastatic <- read.delim(seg_metastatic_path)

seg <- rbind(seg_primary, seg_metastatic) %>% mutate(ID = substr(ID, 1, 26)) %>% dplyr::rename(Sample = ID)

eniclust_wgs <- rbind(
  read.delim(file.path(main_repo_path, 'data','primary_megatable.txt')) %>% mutate(Cohort = "Primary"),
  read.delim(file.path(main_repo_path, 'data','metastatic_megatable.txt')) %>% mutate(Cohort = "Metastatic")
) %>% select(Sample, eniclust = ENiClust, genome_doubled, Cohort, ER)

seg <- seg %>% inner_join(eniclust_wgs)

gr <- regioneR::toGRanges(seg[, c(2:4, 1, 5:22)])

genes_df <- data.frame(
  chr = c(8, 8), 
  start = c(128747680, 128806772), 
  end = c(128755197, 129199347), 
  gene  = c("MYC", "PVT1")
)

genes <- regioneR::toGRanges(genes_df)

overlapMYC <- findOverlaps(genes[1], gr)
overlapPVT1 <- findOverlaps(genes[2], gr)

# Add gene column GRanges 
mcols(gr)$gene <- NA

grMYC <- gr
mcols(grMYC)$gene[subjectHits(overlapMYC)] <- "MYC"
grMYC <- plyranges::filter(grMYC, !is.na(gene))

grPVT1 <- gr
mcols(grPVT1)$gene[subjectHits(overlapPVT1)] <- "PVT1"
grPVT1 <- plyranges::filter(grPVT1, !is.na(gene))

# If gene is split across segments, take the minimum
MYC <- data.frame(mcols(grMYC)) %>% arrange(tcn.em) %>% filter(!duplicated(Sample))
PVT1 <- data.frame(mcols(grPVT1)) %>% arrange(tcn.em) %>% filter(!duplicated(Sample))

# Define AMP based on CN and WGD as usual 
MYC <- MYC %>% mutate(AMP = ifelse(genome_doubled == "True", ifelse(tcn.em > 5, 1, 0), 
                                   ifelse(tcn.em > 4, 1, 0)), 
                      AMP = factor(AMP))
PVT1 <- PVT1 %>% mutate(AMP = ifelse(genome_doubled == "True", ifelse(tcn.em > 5, 1, 0), 
                                     ifelse(tcn.em > 4, 1, 0)), 
                        AMP = factor(AMP))

# Combine data 
amp <- inner_join(MYC %>% dplyr::rename(MYC = AMP) %>% select(-gene), 
                  PVT1 %>% select(Sample, PVT1 = AMP)) %>% 
  mutate(IC9 = ifelse(eniclust == "ic9", "IC9", "Other"), 
         Amplification = ifelse(MYC == 0, ifelse(PVT1 == 0, "Neither", "PVT1 only"), 
                                ifelse(PVT1 == 0, "MYC only", "MYC and PVT1")), 
         Cohort = factor(Cohort, levels = c("Primary", "Metastatic")), 
         Amplification = factor(Amplification, levels = c("Neither", "MYC only", "PVT1 only", "MYC and PVT1")))

source_data <- amp %>% group_by(Cohort, IC9, Amplification) %>% summarize(n = n()) %>% 
  group_by(Cohort, IC9) %>% mutate(total = sum(n), Proportion = n/total)

## Save plot ----
pdf(file = file.path(main_repo_path, "plots", "SupplementaryFigure5J.pdf"), 
    width = 7, height = 5)
ggbarplot(source_data, "IC9", "Proportion", fill = "IC9", alpha = "Amplification", 
          facet.by = c("Cohort"), xlab = "", palette =c("#EE82EE", "gray"), 
          title = "Co-amplification of MYC and PVT1") +
  theme_LM
dev.off()


### SAVE ##########################################################################################
source_data %>% 
  select(Stage = Cohort, IC9, Amplification, Proportion, Total_Samples = total) %>% 
  write.table(file.path(main_repo_path, "data", "SupplementaryFigure5i_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
