### CREATE EXTENDED DATA 9B #############################################################################
# Coamplifications of IDO1 and high risk oncogenes

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

# Data ----
# WGS CN 
primary_cn_path <- "/path/to/FACETS/gene_level/merged/primary.tsv"
metastatic_cn_path <- "/path/to/FACETS/gene_level/merged/metastatic.tsv"

cols <- rep("NULL", 25)
cols[c(1, 2, 16:17, 24:25)] <- "character"
primary_cn <- read.delim(primary_cn_path, colClasses = cols)
metastatic_cn <- read.delim(metastatic_cn_path, colClasses = cols)

filtered <- rbind(primary_cn %>% filter(gene %in% c("FGFR1", "ZNF703", "IDO1", "KAT6A")), 
                  metastatic_cn %>% filter(gene %in% c("FGFR1", "ZNF703", "IDO1", "KAT6A"))) %>% dplyr::rename(Sample = sample)

primary_megatable <- read.delim(file.path(main_repo_path, 'data','primary_megatable.txt'))
metastatic_megatable <- read.delim(file.path(main_repo_path, 'data','metastatic_megatable.txt'))

samples <- rbind(primary_megatable %>% mutate(Stage = "Primary"), 
                 metastatic_megatable %>% mutate(Stage = "Metastatic")) %>% select(Sample, ENiClust, group, Stage) %>% 
  mutate(eniclust = gsub("ic", "IC", ENiClust), 
         Stage = factor(Stage, levels = c("Primary", "Metastatic")))

filtered <- filtered %>% inner_join(samples) 

coamp <- function(data, gene_interest, other_genes) {
  # Return list of samples with an amplification in gene_interest and in any of other_genes
  amp_gene_int <- data %>% filter(grepl("AMP", cn_state), gene == gene_interest) %>% pull(Sample)
  amp_other <- data %>% filter(grepl("AMP", cn_state), gene %in% other_genes) %>% pull(Sample)
  
  return(amp_gene_int[amp_gene_int %in% amp_other])
}

# IDO1 and FGFR1 or ZNF703 
samples_er_high <- samples %>% mutate(IDO1_coamp = Sample %in% coamp(filtered, "IDO1", c("FGFR1", "ZNF703"))) %>% 
  filter(group == "ER+ High")

source_data <- samples_er_high %>% group_by(Stage, eniclust, IDO1_coamp) %>% summarize(n = n()) %>% group_by(Stage, eniclust) %>% 
  mutate(total = sum(n), proportion = n/total) %>% filter(IDO1_coamp)

## Save plot ----
pdf(file.path(main_repo_path, "plots", "ExtendedData9B.pdf"), height = 4, width = 5)
ggbarplot(source_data, "eniclust", "proportion", facet.by = "Stage", color = "eniclust", fill = "eniclust", 
          palette = ic_colors, xlab = "IC subtype", 
          ylab = "Proportion of samples with IDO1 and\nFGFR1 or ZNF703 coamplification") +
  theme_LM +
  rremove("legend")
dev.off()

# Numbers 
filtered %>% filter(grepl("AMP", cn_state), gene == "IDO1") %>% select(eniclust, Stage) %>% table()
samples %>% select(eniclust, Stage) %>% table()

samples_er_high %>% filter(Stage == "Primary", eniclust == "IC6") %>% select(IDO1_coamp) %>% table() # WGS 



### SAVE ##########################################################################################
## Save source data ----
source_data %>% 
  select(Stage, ENiClust = eniclust, IDO1_coamplification = IDO1_coamp, Proportion = proportion, Total_Samples = total) %>% 
  write.table(file.path(main_repo_path, "data", "Extended_Data9b_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
