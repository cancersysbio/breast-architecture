### CREATE SUPPLEMENTARY FIGURE 7A #############################################################################
# Comparison of CIBERSORTx estimates and TME subtypes 

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

# Read data -----
data <- read.delim(file.path(main_repo_path, 'data','rna_megatable.txt'), check.names = F)
data <- data %>% mutate(Stage = factor(Stage, levels = c("Primary", "Metastatic")))

cell_types_interest <- c("Cancer", "B-cells", "T-cells", "Myeloid", "CAFs", "Plasmablasts")
tme_pairs <- list(c("F", "IE/F"), c("IE", "D"), c("F", "D"), c("IE", "IE/F"))

source_data <- data %>% 
  filter(Cohort %in% c("TCGA", "Hartwig", "METABRIC"), !is.na(TME_subtype)) %>%
  dplyr::rename(Cancer = `Cancer Epithelial`) %>% 
  pivot_longer(cols = all_of(cell_types_interest), names_to = "Cell type", values_to = "Proportion") %>%  
  mutate(TME_subtype = factor(TME_subtype, levels = sort(unique(data$TME_subtype))), 
         `Cell type` = factor(`Cell type`, levels = cell_types_interest), 
         Cohort = factor(Cohort, levels = c("TCGA", "METABRIC", "Hartwig")))

## Save plot -----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure7A.pdf"), height = 6, width = 8)
ggboxplot(source_data, "TME_subtype", "Proportion", color = "TME_subtype", 
          palette = tme_colors, facet.by = c("Cohort", "Cell type"), xlab = "TME subtype", 
          ylab = "Estimated proportion of cell type in sample", scales = "free_y",
          title = "") +
  stat_compare_means(aes(label = ..p.signif..), comparisons = tme_pairs, 
                     label.y = c(0.7, 0.8, 0.9, 1)) +
  guides(color = guide_legend(title = "TME subtype")) +
  coord_cartesian(ylim = c(0,1.1)) +
  theme_LM +
  rremove("legend")
dev.off()

### SAVE ##########################################################################################
source_data %>% 
  select(Sample, TME_subtype, Cohort, `Cell type`, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "SupplementaryFigure7a_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
