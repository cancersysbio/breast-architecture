### CREATE EXTENDED DATA 8G #############################################################################
# TME subtype distribution for high risk and HER2+ primary and liver metastatic samples

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

# Read data 
data <- read.delim(file.path(main_repo_path, "data", "rna_megatable.txt"))
data <- data %>% mutate(Stage = factor(Stage, levels = c("Primary", "Metastatic")))

source_data <- data %>% filter(((Cohort == "TCGA")&(Stage == "Primary"))|((Cohort == "Hartwig")&(Stage == "Metastatic")), 
                               !is.na(group), !is.na(TME_subtype))  %>% filter(biopsySite %in% c("Primary", "Liver")) %>% 
  mutate(Stage = as.character(Stage), 
         Stage = ifelse(Stage== "Metastatic", "Liver metastases", "Primary"),
         Stage = factor(Stage, levels = c("Primary", "Liver metastases"))) %>% 
  group_by(TME_subtype, group, Stage) %>% 
  summarize(n = n()) %>% group_by(group, Stage) %>% mutate(total = sum(n), proportion = n/total) 

## Save plot -----
pdf(file.path(main_repo_path, "plots", "ExtendedData8G.pdf"), height = 5, width = 7)
ggbarplot(source_data, "Stage", "proportion", facet.by = "group", color = "TME_subtype", 
          fill = "TME_subtype", nrow = 1, palette = tme_colors, 
          ylab = "Proportion of samples by TME subtype", xlab = "",
          title = "Primary (TCGA) and Metastatic (Hartwig)") +
  geom_text(aes(label = total, y = 1.03, size = 13), check_overlap = T) +
  labs(fill = "TME subtype", color = "TME subtype") +
  theme_LM +
  rotate_x_text(35) +
  rremove("legend")
dev.off()

### SAVE ##########################################################################################
source_data %>% 
  select(TME_subtype, Subtype = group, Stage, Proportion = proportion, Total_Samples = total) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData8g_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
