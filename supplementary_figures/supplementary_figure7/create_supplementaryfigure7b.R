### CREATE SUPPLEMENTARY FIGURE 7B #############################################################################
# TME subtype distribution by subtype in METABRIC 

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

source_data <- data %>% filter(Cohort == "METABRIC", !is.na(group)) %>% 
  group_by(TME_subtype, group) %>% summarize(n = n()) %>% group_by(group) %>% 
  mutate(total = sum(n), Proportion = n/total)

## Save plot -----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure7B.pdf"), height = 4, width = 6)
ggbarplot(source_data, "group", "Proportion", color = "TME_subtype", 
          fill = "TME_subtype", palette = tme_colors, 
          ylab = "Proportion of samples by TME subtype", xlab = "",
          title = "Primary (METABRIC)") +
  geom_text(aes(label = total, y = 1.03, size = 13), check_overlap = T) +
  labs(fill = "TME subtype", color = "TME subtype") +
  theme_LM +
  rremove("legend")
dev.off()


### SAVE ##########################################################################################
source_data %>% 
  select(TME_subtype, Subtype = group, Total_Samples = total, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "SupplementaryFigure7b_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
