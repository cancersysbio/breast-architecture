# CREATE SUPPLEMENTARY FIGURE 5c ----
# Barplot of ecDNA+ samples split by HRD

# PREAMBLE ----
library(rstudioapi)
library(tidyr)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(yaml)
library(ggplot2)
library(ggpubr)
library(patchwork)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# FUNCTIONS ----
theme_LM = theme_classic() + grids()

plot_prop_samples <- function(df, label_hrd, label_group, alpha_col, pivot_col, legend_name) {
  df %>%
    ggplot(aes(x=!!sym(label_hrd), fill = !!sym(label_group), alpha=!!sym(alpha_col))) +
    geom_bar(stat='count', position='fill', colour='black', show.legend = c(fill = FALSE, alpha = TRUE)) +
    facet_grid(as.formula(paste0(pivot_col, "~", label_group)),
               scales='free_x') +
    scale_fill_manual(values = subgroup_colors) +
    scale_alpha_discrete(range = c(0.1, 1),
                         labels = c('FALSE'='No', 'TRUE'='Yes'),
                         name = legend_name) +
    labs(y='Proportion of Samples') +
    theme_LM + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
                     axis.text.y = element_text(size = 18),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 18),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 12),
                     title = element_text(size = 20),
                     legend.position = 'right',
                     strip.text.x = element_text(size = 20),
                     strip.text.y = element_text(size = 18))
  
}

fisher_test_single_category <- function(contingency_table, hrd_labels, category1, category2) {
  fisher_matrix <- as.matrix(contingency_table[,-c(1,2)])
  rownames(fisher_matrix) <- hrd_labels
  
  fisher_test <- fisher.test(fisher_matrix)
  
  return(list(
    category1 = category1, 
    category2 = category2, 
    p_value = fisher_test$p.value, 
    odds_ratio = fisher_test$estimate
  ))
}

subgroup_colors <- c("ER+ High" = "#fa954eff", 
                     "ER+ Typical"="#a1c1d4ff", 
                     "HER2+" = "#bd854dff", 
                     "IC10" = "#904dbdff", 
                     "IC4ER-" = "#c2b7dbff",
                     "TNBC" = "#5B308B")

# MAIN ----
pri_samples <- fread(file.path(main_repo_path, "data", "SupplementaryFigure5d_sourcetable.txt"))
t1 <- pri_samples %>%
  filter(Type == 'AmpliconArchitect' & label_group == 'ER+ High') %>%
  group_by(label_group, label_chord_only, ecDNA_detected) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = ecDNA_detected,
              values_from = n) %>%
  fisher_test_single_category(., c('BRCA2-like', 'non-HRD-like'), 'TRUE', 'FALSE')
t2 <- pri_samples %>%
  filter(Type == 'JaBbA' & label_group == 'ER+ High') %>%
  group_by(label_group, label_chord_only, ecDNA_detected) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = ecDNA_detected,
              values_from = n) %>%
  fisher_test_single_category(., c('BRCA2-like', 'non-HRD-like'), 'TRUE', 'FALSE')
stat_results <- bind_rows(t1, t2)

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure5d.pdf"), height=6, width=7.5)
plot_prop_samples(df=pri_samples, 
                  label_hrd='label_chord_only', 
                  label_group='label_group', 
                  alpha_col='ecDNA_detected', 
                  pivot_col='Type', 
                  legend_name='ecDNA')
dev.off()

fwrite(stat_results, file.path(main_repo_path, "plots", "SupplementaryFigure5d_stats.tsv"))