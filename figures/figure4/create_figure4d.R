### CREATE FIGURE 4D #############################################################################
# Alterations in immune pathways by type of alteration

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

# Figure 4D ----
data_4d <- data %>% mutate(Pathway = gsub("Antigen Presentation ", "Antigen\nPresentation", Pathway), 
                           Pathway = gsub("IFN-γ pathway", "IFN-g\npathway", Pathway))
pathways <- c("HLA-II", "HLA-I", "Epigenetic", "IFN-g\npathway", "Antigen\nPresentation", "PD-L1", "CD58")
data_4d$Pathway <- factor(data_4d$Pathway, levels = pathways)

total_groups <- rbind(
  data.frame(table(primary_megatable$group)) %>% rename(group = Var1, total = Freq) %>% mutate(Stage = "Primary"),
  data.frame(table(metastatic_megatable$group)) %>% rename(group = Var1, total = Freq) %>% mutate(Stage = "Metastatic")
)

data_n_alt <- data_4d %>% select(Sample, Alteration, Pathway, Stage) %>% distinct() %>% group_by(Sample, Pathway, Stage) %>% 
  mutate(n_alt = n())

data_n_alt <- rbind(
  data_n_alt %>% filter(n_alt == 1) %>% select(-n_alt), 
  data_n_alt %>% filter(n_alt > 1) %>% select(-n_alt) %>% mutate(Alteration = "Multiple") %>% distinct()
)

total_samples <- data.frame(Stage = c("Primary", "Metastatic"), total = c(nrow(primary_megatable), nrow(metastatic_megatable)))

source_data <- data_n_alt %>% group_by(Alteration, Pathway, Stage) %>% summarize(n = n()) %>% 
  group_by(Pathway, Stage) %>% left_join(total_samples) %>% 
  mutate(Proportion = n/total)
source_data$Stage <- factor(source_data$Stage, levels = c("Primary", "Metastatic"))

## Save plot ---- 
pdf(file.path(main_repo_path, "plots", "Figure4d.pdf"), height = 1.77, width = 7)
ggbarplot(source_data, "Stage", "Proportion", color = "Alteration", fill = "Alteration", 
          facet.by = "Pathway", nrow = 1, palette = immune_alt_colors, xlab = "", 
          ylab = "Proportion of samples\nwith alterations") +
  theme_LM +
  rotate_x_text(30) +
  theme(legend.position = "right", 
        legend.margin=margin(c(0,0,0,0)),
        legend.key.size = unit(0.3, 'cm'),
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        axis.title.y=element_text(size=6), 
        strip.text.x=element_text(size=6), 
        axis.title.x=element_blank(), 
        axis.text.y=element_text(size=6), 
        axis.text.x=element_text(size=6))
dev.off()


### SAVE ##########################################################################################
source_data %>% 
  select(Alteration, Pathway, Stage, Total_Samples = total, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "Figure4d_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
