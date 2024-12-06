### CREATE EXTENDED DATA FIGURE 9D #############################################################################
# Barplot of GIE in ER+ Typical IDC vs ILC (Metabric)

### PREAMBLE #####################################################################################
library(rstudioapi)
library(dplyr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(yaml)
library(grid)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### FUNCTIONS ##########################################################################################
theme_LM = theme_classic() + grids()

plot_gie_prop <- function(df) {
  df %>%
    select(Sample, group, HistoType, gie) %>%
    ggplot(aes(x=HistoType, fill=HistoType, alpha=gie)) +
    geom_bar(position="fill", stat="count", color="black") + 
    stat_count(geom = "text", 
               position=position_fill(vjust = 0.5),
               aes(label = ..count..),
               color="black", 
               size=4) +
    labs(y="Proportion of samples \n with altered GIE genes", x="") + 
    facet_wrap(~group, labeller = labeller(name = c(group = "ER+ Typical"))) +
    scale_fill_manual(values = c('#a3cfffff','#335e8dff')) +
    scale_alpha_discrete(range = c(0.4, 1), 
                         guide = guide_legend(override.aes = list(fill = "darkgray"),
                                              title = "Altered")) +
    guides(alpha = guide_legend(override.aes = list(fill = "darkgray"),
                                title = "GIE",
                                ncol = 1),
           fill = "none", color = "none") +
    theme_LM +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14)) +
    ylim(0,1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

plot_gie_type_prop <- function(df) {
  df %>%
    mutate(Pathway=factor(Pathway, levels=c("HLA-II", "HLA-I", "Epigenetic",
                                            "IFN-γ pathway", "Antigen Presentation", "CD58"))) %>%
    ggplot(aes(x=HistoType, fill=Pathway)) +
    geom_bar(position="fill", stat="count") + 
    stat_count(geom = "text", 
               position=position_fill(vjust = 0.5),
               aes(label = ..count..),
               color="black", 
               size=3) +
    labs(y="", x="") + 
    facet_wrap(~group, labeller = labeller(name = c(group = "ER+ Typical"))) +
    scale_fill_brewer(palette = "Set2", name="Pathway") +
    theme_LM +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14)) +
    ylim(0,1)
  dev.off()
}

### MAIN ##########################################################################################
gie_wgs <- read.delim(file.path(main_repo_path, "data", "source_data", "extended_data","ExtendedData8d_sourcetable.txt"), sep = ",")
gie_wgs_count <- gie_wgs %>% 
  count(Sample) %>%
  mutate(Sample=substr(Sample, 1, 21))
all_wgs_count <- all_wgs %>%
  filter(HistoType %in% c("IDC", "ILC")) %>%
  filter(group == "ER+ Typical") %>%
  left_join(gie_wgs_count, by = "Sample") %>%
  mutate(n=case_when(is.na(n) ~ 0 ,
                     T ~ n)) %>%
  mutate(gie=case_when(n==0 ~ F,
                       T ~ T))

gie_wgs_types <- gie_wgs %>%
  select(Sample, Pathway, Alteration, HistoType, group) %>%
  filter(HistoType %in% c("IDC", "ILC") &
           group == "ER+ Typical")

df <- data.frame(
  "ILC" = c(all_wgs_count %>%
              count(group, HistoType, gie) %>% 
              filter(HistoType == "ILC" & gie == TRUE) %>%
              pull(n), 
            all_wgs_count %>%
              count(group, HistoType, gie) %>% 
              filter(HistoType == "ILC" & gie == FALSE) %>%
              pull(n)),
  "IDC" = c(all_wgs_count %>%
              count(group, HistoType, gie) %>% 
              filter(HistoType == "IDC" & gie == TRUE) %>%
              pull(n), 
            all_wgs_count %>%
              count(group, HistoType, gie) %>% 
              filter(HistoType == "IDC" & gie == FALSE) %>%
              pull(n)),
  row.names = c("TRUE", "FALSE"),
  stringsAsFactors = FALSE
)
ss <- fisher.test(df, alternative = 'lower')
ss_df <- data.frame(p.value=ss$p.value,
                    or=ss$estimate)


### SAVE ##########################################################################################
pdf(file.path(main_repo_path, "plots", "ExtendedData9d_left.pdf"), height=4, width=3.5)
plot_gie_prop(all_wgs_count)
dev.off()

cairo_pdf(file.path(main_repo_path, "plots", "ExtendedData9d_right.pdf"), height=4, width=4)
plot_gie_type_prop(gie_wgs_types)
dev.off()

ss_df %>%
  write.table(file.path(main_repo_path, "data", "ExtendedData9d_stats.txt"), 
              quote = F, row.names = F, sep = "\t")