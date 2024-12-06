### CREATE EXTENDED DATA 8F #############################################################################
# TME subtype distribution for ER+ vs ER- samples

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

data_ext8f <- data %>% filter(((Cohort == "TCGA")&(Stage == "Primary"))|((Cohort == "Hartwig")&(Stage == "Metastatic")), 
                              !is.na(group), !is.na(TME_subtype)) %>%
  mutate(ER_groups = ifelse(
    (group %in% c("ER+ High", "ER+ Typical"))|((eniclust == "ic5")&(ER == 1)), "ER+", 
    ifelse((group %in% c("IC10", "IC4ER-"))|((eniclust == "ic5")&(ER == 0)), "ER-", NA)
  ), ER_groups = factor(ER_groups, levels = c("ER+", "ER-"))) %>% filter(!is.na(ER_groups))

source_data <- data_ext8f %>% group_by(TME_subtype, ER_groups, Stage) %>% 
  summarize(n = n()) %>% group_by(ER_groups, Stage) %>% mutate(total = sum(n), proportion = n/total)

## Save plot -----
pdf(file.path(main_repo_path, "plots", "ExtendedData8F.pdf"), height = 5, width = 4)
ggbarplot(source_data, "Stage", "proportion", facet.by = "ER_groups", color = "TME_subtype", 
          fill = "TME_subtype", nrow = 1, palette = tme_colors, 
          ylab = "Proportion of samples by TME subtype", xlab = "",
          title = "Primary (TCGA) and\nMetastatic (Hartwig)") +
  geom_text(aes(label = total, y = 1.03, size = 13), check_overlap = T) +
  labs(fill = "TME subtype", color = "TME subtype") +
  theme_LM +
  rotate_x_text(35) +
  rremove("legend")
dev.off()


### SAVE ##########################################################################################
source_data %>% 
  select(TME_subtype, ER = ER_groups, Stage, Proportion = proportion, Total_Samples = total) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData8f_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
