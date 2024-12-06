### CREATE EXTENDED DATA 9C #############################################################################
# Coamplifications of IDO1 and high risk oncogenes vs TME 

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
# TCGA data
cn_tcga_path <- "/path/to/TCGA/WES/gene_level/CN/data.txt"
wgd_tcga_path <- "/path/to/TCGA/WES/WGD/data.txt" 

cn_tcga <- read.delim(cn_tcga_path, check.names = F) %>% mutate(Sample = substr(sample, 1, 12))

tcga_tme <- read.delim(file.path(main_repo_path, "data", "rna_megatable.txt")) %>%
  filter(Cohort == "TCGA") %>% select(Sample, TME_subtype, eniclust)

facets_data <- read.delim(wgd_tcga_path) %>% 
  mutate(Sample = substr(sample_x, 1, 12)) %>% select(Sample, genome_doubled)

tcga_data <- cn_tcga %>% select(Sample, FGFR1, ZNF703, IDO1, KAT6A) %>% inner_join(tcga_tme) %>% 
  inner_join(facets_data)

# METABRIC data
cn_metabric_path <- "/path/to/METABRIC/gene_level/CN/data.txt"
wgd_metabric_path <- "/path/to/METABRIC/WGD/data.txt"

metabric_tme <- read.delim(file.path(main_repo_path, "data", "rna_megatable.txt")) %>%
  filter(Cohort == "METABRIC") %>% select(Sample, TME_subtype, eniclust)

cn_data_metabric <- read.delim(cn_metabric_path) %>% dplyr::rename(Sample = sample)
wgd_metabric <- read.delim(wgd_metabric_path) %>% select(Sample, genome_doubled)

metabric_data <- cn_data_metabric %>% select(Sample, FGFR1, ZNF703, IDO1, KAT6A) %>% inner_join(metabric_tme) %>% 
  inner_join(wgd_metabric) 

# TCGA
tcga_data <- tcga_data %>% mutate(
  FGFR1_amp = ((FGFR1 > 5)&(genome_doubled == "True"))|((FGFR1 > 4)&(genome_doubled == "False")), 
  ZNF703_amp = ((ZNF703 > 5)&(genome_doubled == "True"))|((ZNF703 > 4)&(genome_doubled == "False")), 
  IDO1_amp = ((IDO1 > 5)&(genome_doubled == "True"))|((IDO1 > 4)&(genome_doubled == "False")), 
  KAT6A_amp = ((KAT6A > 5)&(genome_doubled == "True"))|((KAT6A > 4)&(genome_doubled == "False")), 
  IDO1_coamp = (FGFR1_amp | ZNF703_amp)&IDO1_amp
)

model_data <- tcga_data %>% filter(eniclust == "ic6") %>% mutate(IE = TME_subtype %in% c("IE", "IE/F"))

model <- glm(IE ~ IDO1_coamp + FGFR1, data = model_data, family = "binomial")

# METABRIC
metabric_data <- metabric_data %>% mutate(
  FGFR1_amp = ((FGFR1 > 5)&(genome_doubled == 1))|((FGFR1 > 4)&(genome_doubled == 0)), 
  ZNF703_amp = ((ZNF703 > 5)&(genome_doubled == 1))|((ZNF703 > 4)&(genome_doubled == 0)), 
  IDO1_amp = ((IDO1 > 5)&(genome_doubled == 1))|((IDO1 > 4)&(genome_doubled == 0)), 
  KAT6A_amp = ((KAT6A > 5)&(genome_doubled == 1))|((KAT6A > 4)&(genome_doubled == 0)), 
  IDO1_coamp = (FGFR1_amp | ZNF703_amp)&IDO1_amp
)

model_data_metabric <- metabric_data %>% filter(eniclust == "ic6") %>% mutate(IE = TME_subtype %in% c("IE", "IE/F"))

model_metabric <- glm(IE ~ IDO1_coamp + FGFR1, data = model_data_metabric, family = "binomial")

# Plot proportions
ic6 <- rbind(
  model_data_metabric %>% mutate(Cohort = "METABRIC"), 
  model_data %>% mutate(Cohort = "TCGA")
) %>% mutate(`Immune enriched` = ifelse(IE, "Yes", "No"), `Immune enriched` = factor(`Immune enriched`, levels = c("Yes", "No")), 
             Coamplification = ifelse(IDO1_coamp, "Yes", "No"), Coamplification = factor(Coamplification, levels = c("Yes", "No")))

test_tcga <- ic6 %>% filter(Cohort == "TCGA") %>% select(`Immune enriched`, IDO1_coamp) %>% table() %>% fisher.test(alternative = "g")
test_metabric <- ic6 %>% filter(Cohort == "METABRIC") %>% select(`Immune enriched`, IDO1_coamp) %>% table() %>% fisher.test(alternative = "g")

p_or <- data.frame(
  Cohort = c("TCGA", "METABRIC"), 
  OR = c(test_tcga$estimate[["odds ratio"]], test_metabric$estimate[["odds ratio"]]), 
  p = c(test_tcga$p.value, test_metabric$p.value)
)

ie <- brewer.pal(2, "RdBu")
ie <- ie[c(1, 3)]
names(ie) <- c("Yes", "No")

source_data <- ic6 %>% group_by(`Immune enriched`, Coamplification, Cohort) %>% summarize(n = n()) %>% group_by(Coamplification, Cohort) %>% 
  mutate(total = sum(n), proportion = n/total) %>% left_join(p_or) 

## Save plot ----
pdf(file.path(main_repo_path, "plots", "ExtendedData9C.pdf"), height = 4, width = 5)
ggbarplot(source_data, "Coamplification", "proportion", facet.by = "Cohort", color = "Immune enriched", fill = "Immune enriched", 
          palette = ie, ylab = "Proportion of IC6 samples", xlab = "Coamplification of IDO1 and\nFGFR1 or ZNF703") + 
  geom_text(aes(label = paste0("OR = ", round(OR, 3), "\np = ", round(p, 3)), 
                x = 1.5, y = 1.075), check_overlap = T) +
  scale_y_continuous(breaks = (0:5)/5, limits = c(0, 1.09)) +
  theme_LM
dev.off()

# Numbers 
tcga_data %>% filter(eniclust == "ic6") %>% select(IDO1_coamp) %>% table() # TCGA
metabric_data %>% filter(eniclust == "ic6") %>% select(IDO1_coamp) %>% table() # METABRIC


### SAVE ##########################################################################################
source_data %>% 
  select(`Immune enriched`, Coamplification, Cohort, Proportion = proportion, Total_Samples = total, OR, p) %>% 
  write.table(file.path(main_repo_path, "data", "Extended_Data9c_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
