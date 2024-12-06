### CREATE EXTENDED DATA 8C #############################################################################
# TME subtype distribution for high risk and HER2+ primary and metastatic samples

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

data_ext8c <- data %>% filter(((Cohort == "TCGA")&(Stage == "Primary"))|((Cohort == "Hartwig")&(Stage == "Metastatic")), 
                              !is.na(group), !is.na(TME_subtype)) %>% 
  filter(eniclust %in% c("ic1", "ic2", "ic5", "ic6", "ic9")) %>% 
  mutate(eniclust = gsub("ic", "IC", eniclust))

data_ext8c <- data_ext8c %>% filter(!((eniclust == "IC5")&(is.na(ER)))) %>% 
  mutate(eniclust = ifelse((eniclust == "IC5")&(ER == 1), "IC5/ER+", 
                           ifelse((eniclust == "IC5")&(ER == 0), "IC5/ER-", eniclust)), 
         eniclust = factor(eniclust, levels = c("IC1", "IC2", "IC6", "IC9", "IC5/ER+", "IC5/ER-")))

source_data <- data_ext8c %>% 
  group_by(TME_subtype, eniclust, Stage) %>% 
  summarize(n = n()) %>% group_by(eniclust, Stage) %>% mutate(total = sum(n), proportion = n/total)  

## Save plot -----
pdf(file.path(main_repo_path, "plots", "ExtendedData8C.pdf"), height = 5, width = 7)
ggbarplot(source_data, "Stage", "proportion", facet.by = "eniclust", color = "TME_subtype", 
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
  select(TME_subtype, ENiClust = eniclust, Stage, Proportion = proportion, Total_Samples = total) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData8c_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
