### CREATE FIGURE 4C #############################################################################
# Alterations in immune pathways by subtype

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

# Figure 4C ----
total_samples <- data.frame(Stage = c("Primary", "Metastatic"), total = c(nrow(primary_megatable), nrow(metastatic_megatable)))

total_groups <- rbind(
  data.frame(table(primary_megatable$group)) %>% rename(group = Var1, total = Freq) %>% mutate(Stage = "Primary"),
  data.frame(table(metastatic_megatable$group)) %>% rename(group = Var1, total = Freq) %>% mutate(Stage = "Metastatic")
)

alterations_by_pathway <- data %>% select(Pathway, Sample, group, Stage) %>% unique() %>% 
  group_by(group, Pathway, Stage) %>% summarize(n = n()) %>% left_join(total_groups) %>% 
  mutate(proportion = n/total, 
         Stage = factor(Stage, levels = c("Primary", "Metastatic"))) 

alterations_by_pathway_no_group <- data %>% select(Pathway, Sample, Stage) %>% unique() %>% 
  group_by(Pathway, Stage) %>% summarize(n = n()) %>% left_join(total_samples) %>% 
  mutate(proportion = n/total, 
         Stage = factor(Stage, levels = c("Primary", "Metastatic"))) 

alteration_count <- data %>% select(Pathway, Sample, group, Stage) %>% unique() %>% 
  group_by(group, Sample, Stage) %>% summarize(pathways_altered = n()) %>% 
  group_by(group, pathways_altered, Stage) %>% summarise(n = n()) %>% 
  left_join(total_groups) %>% mutate(proportion = n/total, 
                                     Stage = factor(Stage, levels = c("Primary", "Metastatic")), 
                                     pathways_altered = ifelse(pathways_altered > 2, 3, pathways_altered)) 

a <- ggbarplot(alteration_count, "group", "proportion", fill = "group", color = "group", 
               facet.by = "Stage", nrow = 1,
               palette = group_colors, xlab = "", ylab = "Proportion of samples\nwith altered GIE genes") +
  theme_LM +
  geom_bar(stat="identity", fill="white", aes(alpha=pathways_altered), width = 0.7) +
  guides(alpha = guide_legend(override.aes = list(fill = "darkgray"),
                              title = "Pathways\naltered", 
                              ncol = 3), 
         fill = "none", color = "none") +
  scale_alpha(breaks = 1:3, labels = c("3+", "2", "1")) +
  theme(legend.position = "top", 
        legend.margin=margin(c(0,0,0,0)),
        legend.key.size = unit(0.3, 'cm'),
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        axis.title.y=element_text(size=6), 
        strip.text.x=element_text(size=6), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.y=element_text(size=6))

alterations_by_pathway$Pathway[alterations_by_pathway$Pathway == "Antigen Presentation "] <- "Antigen\nPresentation"
alterations_by_pathway$Pathway[alterations_by_pathway$Pathway == "IFN-γ pathway"] <- "IFN-g\npathway"
pathways <- alterations_by_pathway_no_group %>% filter(Stage == "Primary") %>% mutate(Pathway = ifelse(Pathway == "Antigen Presentation ", "Antigen\nPresentation", Pathway), 
                                                                                      Pathway = ifelse(Pathway == "IFN-γ pathway", "IFN-g\npathway", Pathway)) %>%
  arrange(desc(proportion)) %>% pull(Pathway)
alterations_by_pathway$Pathway <- factor(alterations_by_pathway$Pathway, levels = pathways)

b <- ggbarplot(alterations_by_pathway, "group", "proportion", fill = "group", color = "group", 
               facet.by = c("Stage", "Pathway"),
               palette = group_colors, xlab = "", ylab = "") +
  theme_LM +
  rotate_x_text(35) +
  guides(fill = guide_legend(title = "IC subgroup", ncol = 5), 
         color = guide_legend(title = "IC subgroup", ncol = 5)) +
  theme(legend.position = "top", 
        legend.margin=margin(c(0,0,0,0)),
        legend.key.size = unit(0.3, 'cm'),
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 6), 
        strip.text.x=element_text(size=6),
        strip.text.y=element_text(size=6),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.text.y=element_text(size=6))

## Save plot ----
pdf(file.path(main_repo_path, "plots", "Figure4c.pdf"), height = 2, width = 7)
a+b +
  plot_layout(widths = c(3, 10))
dev.off()


# Numbers
## Number of altered samples
data %>% select(Sample, Stage) %>% distinct() %>% pull(Stage) %>% table()

## Number of samples by pathways altered
data %>% select(Sample, Pathway, Stage) %>% distinct() %>% group_by(Sample, Stage) %>% summarize(alt_path = n()) %>% 
  group_by(Stage, alt_path) %>% summarize(n = n()) 

## Number of altered samples by group
data %>% select(Sample, group, Stage) %>% distinct() %>% select(group, Stage) %>% table()

## Alterations in epigenetic pathway
data %>% filter(Pathway == "Epigenetic", Stage == "Primary") %>% pull(Sample) %>% unique() %>% length()

## SV/complex alteration prevalence 
data %>% filter(Alteration %in% c("Complex non-cyclic amplification", "Cyclic amplification", "Structural variant")) %>%
  select(Sample, Stage) %>% distinct() %>% pull(Stage) %>% table()

data %>% filter(Alteration %in% c("Structural variant")) %>%
  select(Sample, Stage) %>% distinct() %>% pull(Stage) %>% table()

## SV/complex # alterations by pathway
data %>% filter(Alteration %in% c("Complex non-cyclic amplification", "Cyclic amplification", "Structural variant")) %>%
  select(Pathway, Stage) %>% table()
data %>% filter(Alteration %in% c("Complex non-cyclic amplification", "Cyclic amplification", "Structural variant")) %>%
  select(Stage) %>% table()



### SAVE ##########################################################################################
# Left
alteration_count %>% 
  select(Subtype = group, Pathways_Altered = pathways_altered, Stage, Proportion = proportion) %>% 
  write.table(file.path(main_repo_path, "data", "Figure4c_sourcetable_left.txt"), 
              quote = F, row.names = F, sep = "\t")

# Right 
alterations_by_pathway %>% 
  select(Subtype = group, Pathway, Stage, Proportion = proportion) %>% 
  write.table(file.path(main_repo_path, "data", "Figure4c_sourcetable_right.txt"), 
              quote = F, row.names = F, sep = "\t")
