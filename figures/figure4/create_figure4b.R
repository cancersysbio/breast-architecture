### CREATE FIGURE 4B #############################################################################
# TME subtype distribution by groups for primary and metastatic samples

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(yaml)

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

data_fig4b <- data %>% filter(((Cohort == "TCGA")&(Stage == "Primary"))|((Cohort == "Hartwig")&(Stage == "Metastatic")), 
                              !is.na(group), !is.na(TME_subtype)) 

source_data <- data_fig4b %>% group_by(TME_subtype, group, Stage) %>% 
  summarize(n = n()) %>% group_by(group, Stage) %>% mutate(total = sum(n), proportion = n/total) %>% 
  mutate(Stage = recode(Stage, "Primary" = "Primary (TCGA)", "Metastatic" = "Metastatic (Hartwig)"))

# Save plot ----
pdf(file.path(main_repo_path, "plots", "Figure4b.pdf"), height = 1.77, width = 3.35)
ggbarplot(source_data, "group", "proportion", facet.by = "Stage", color = "TME_subtype", 
          fill = "TME_subtype", ncol = 2, palette = tme_colors, 
          ylab = "Proportion of samples\nby TME subtype", 
          #title = "Primary (TCGA) and Metastatic (Hartwig)", 
          xlab = "") +
  geom_text(aes(label = total, y = 1.05), size = 7, size.unit = "pt", check_overlap = T) +
  labs(fill = "TME subtype", color = "TME subtype") +
  theme_LM +
  rotate_x_text(30) + 
  rremove("legend") +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=7), 
        strip.text.x=element_text(size=8), 
        axis.title.x=element_blank())
dev.off()


# Statistical tests --------
tme_data <- data_fig4b %>% mutate(Stage = factor(Stage, levels = c("Primary", "Metastatic")), 
                                  IE = ifelse(TME_subtype %in% c("IE", "IE/F"), "Yes", "No"),
                                  IE = factor(IE, levels = c("Yes", "No")))

## Compare enrichment of IE and IE/F in TNBC ------
tme_data %>% mutate(group2 = ifelse(group %in% c("ER+ High", "ER+ Typical", "HER2+"), "ER+/HER2+", "TNBC")) %>% 
  filter(biopsySite == "Primary") %>% mutate(group2 = factor(group2, levels = c("TNBC", "ER+/HER2+"))) %>% 
  select(group2, IE) %>% table() %>% fisher.test(alternative = "g") 

tme_data %>% mutate(group2 = ifelse(group %in% c("ER+ High", "ER+ Typical", "HER2+"), "ER+/HER2+", "TNBC")) %>% 
  filter(biopsySite != "Primary") %>% mutate(group2 = factor(group2, levels = c("TNBC", "ER+/HER2+"))) %>% 
  select(group2, IE) %>% table() %>% fisher.test(alternative = "g") 


## Compare enrichment of D in ER+ High and HER2+ ---------
tme_data %>% mutate(group2 = ifelse(group %in% c("ER+ High", "HER2+"), "ER+ High/HER2+", "ER+ Typical/TNBC"), 
                    D = factor(ifelse(TME_subtype == "D", 1, 0), levels = c(1, 0))) %>% 
  filter(biopsySite == "Primary") %>%
  select(group2, D) %>% table() %>% fisher.test(alternative = "g") 

tme_data %>% mutate(group2 = ifelse(group %in% c("ER+ High", "HER2+"), "ER+ High/HER2+", "ER+ Typical/TNBC"), 
                    D = factor(ifelse(TME_subtype == "D", 1, 0), levels = c(1, 0))) %>% 
  filter(biopsySite != "Primary") %>%
  select(group2, D) %>% table() %>% fisher.test(alternative = "g") 

## Compare enrichment of F and IE/F in ER+ Typical and IC4ER-  ------
tme_data %>% mutate(group2 = ifelse(group %in% c("ER+ Typical", "IC4ER-"), "ER+ Typical/IC4ER-", "ER+ High/HER2+/IC10"),
                    group2 = factor(group2, levels = c("ER+ Typical/IC4ER-", "ER+ High/HER2+/IC10")),
                    `F` = ifelse(TME_subtype %in% c("F", "IE/F"), "Yes", "No"), 
                    `F` = factor(`F`, levels = c("Yes", "No"))) %>% 
  filter(biopsySite == "Primary") %>%
  select(group2, `F`) %>% table() %>% fisher.test(alternative = "g") 

## IE in primary to met by ER ----
tme_data %>% filter(ER == 0) %>% select(Stage, IE) %>% 
  table() %>% fisher.test(alternative = "g")

tme_data %>% filter(ER == 1) %>% select(Stage, IE) %>% 
  table() %>% fisher.test(alternative = "g")


## ER only ----
tme_data %>% filter(grepl("ER+", group, fixed = T)) %>%
  mutate(D = factor(ifelse(TME_subtype == "D", 1, 0), levels = c(1, 0))) %>% 
  filter(biopsySite == "Primary") %>%
  select(group, D) %>% table() %>% fisher.test(alternative = "g") 

tme_data %>% filter(grepl("ER+", group, fixed = T)) %>%
  mutate(D = factor(ifelse(TME_subtype == "D", 1, 0), levels = c(1, 0))) %>% 
  filter(biopsySite != "Primary") %>%
  select(group, D) %>% table() %>% fisher.test(alternative = "g") 



### SAVE ##########################################################################################
source_data %>% 
  select(TME_subtype, Subtype = group, Stage, Total_Samples = total, Proportion = proportion) %>% 
  write.table(file.path(main_repo_path, "data", "Figure4b_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
