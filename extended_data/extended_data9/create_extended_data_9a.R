### CREATE EXTENDED DATA 9A #############################################################################
# Alterations in immune pathways for high risk subtypes

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(patchwork)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

# Read data -- generated in immune_pathway_alterations.R
immune_alterations_data_path <- file.path(main_repo_path, 'data','immune_pathway_alterations.txt')
if(!file.exists(immune_alterations_data_path)) {
  stop("Error: Immune alteration data not found, run immune_pathway_alterations.R analysis script. ")
} 

data <- read.delim(immune_alterations_data_path)

primary_megatable <- read.delim(file.path(main_repo_path, 'data','primary_megatable.txt'))
metastatic_megatable <- read.delim(file.path(main_repo_path, 'data','metastatic_megatable.txt'))

samples <- rbind(primary_megatable, metastatic_megatable)

data <- data %>% inner_join(samples) %>% filter(!is.na(group))

# Extended Data 9A ----
total_high_risk <- rbind(
  data.frame(table(primary_megatable$ENiClust)) %>% rename(ENiClust = Var1, total = Freq) %>% mutate(Stage = "Primary"),
  data.frame(table(metastatic_megatable$ENiClust)) %>% rename(ENiClust = Var1, total = Freq) %>% mutate(Stage = "Metastatic")
) %>% filter(ENiClust %in% c("ic1", "ic2", "ic6", "ic9")) %>% mutate(ENiClust = gsub("ic", "IC", ENiClust))

data_ext9a <- data %>% filter(group == "ER+ High") %>% mutate(ENiClust = gsub("ic", "IC", ENiClust))

alterations_by_pathway <- data_ext9a %>% select(Pathway, Sample, ENiClust, Stage) %>% unique() %>% 
  group_by(ENiClust, Pathway, Stage) %>% summarize(n = n()) %>% left_join(total_high_risk) %>% 
  mutate(proportion = n/total, 
         Stage = factor(Stage, levels = c("Primary", "Metastatic"))) 

alteration_count <- data_ext9a %>% select(Pathway, Sample, ENiClust, Stage) %>% unique() %>% 
  group_by(ENiClust, Sample, Stage) %>% summarize(pathways_altered = n()) %>% 
  group_by(ENiClust, pathways_altered, Stage) %>% summarise(n = n()) %>% 
  left_join(total_high_risk) %>% mutate(proportion = n/total, 
                                        Stage = factor(Stage, levels = c("Primary", "Metastatic")), 
                                        pathways_altered = ifelse(pathways_altered > 2, 3, pathways_altered), 
                                        ENiClust = factor(ENiClust, levels = names(ic_colors))) 

a <- ggbarplot(alteration_count, "ENiClust", "proportion", fill = "ENiClust", color = "ENiClust", 
               facet.by = "Stage", nrow = 1,
               palette = ic_colors, xlab = "", ylab = "Proportion of samples with altered GIE genes") +
  theme_LM +
  rotate_x_text(35) +
  geom_bar(stat="identity", fill="white", aes(alpha=pathways_altered), width = 0.7) +
  guides(alpha = guide_legend(override.aes = list(fill = "darkgray"),
                              title = "Pathways altered", 
                              ncol = 1), 
         fill = "none", color = "none") +
  scale_alpha(breaks = 1:3, labels = c("3+", "2", "1")) +
  theme(legend.position = c(0.22, 0.85), 
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10)) +
  coord_cartesian(ylim = c(0, 0.7))

alterations_by_pathway$Pathway[alterations_by_pathway$Pathway == "Antigen Presentation "] <- "Antigen\nPresentation"
pathways <- c("HLA-II", "HLA-I", "Epigenetic", "IFN-γ pathway", "Antigen\nPresentation", "PD-L1")
alterations_by_pathway$Pathway <- factor(alterations_by_pathway$Pathway, levels = pathways)

b <- ggbarplot(alterations_by_pathway, "ENiClust", "proportion", fill = "ENiClust", color = "ENiClust", 
               facet.by = c("Stage", "Pathway"),
               palette = ic_colors, xlab = "", ylab = "") +
  theme_LM +
  rotate_x_text(35) +
  rremove("legend") 

## Save plot ----
pdf(file.path(main_repo_path, "plots", "ExtendedData9A.pdf"), height = 5, width = 14)
a+b +
  plot_layout(widths = c(3, 7))
dev.off()

### SAVE ##########################################################################################
# Left
alteration_count %>% 
  select(ENiClust, Pathways_Altered = pathways_altered, Stage, Proportion = proportion) %>% 
  write.table(file.path(main_repo_path, "data", "Extended_Data9a_left_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")

# Right 
alterations_by_pathway %>% 
  select(ENiClust, Pathway, Stage, Proportion = proportion) %>% 
  write.table(file.path(main_repo_path, "data", "Extended_Data9a_right_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
