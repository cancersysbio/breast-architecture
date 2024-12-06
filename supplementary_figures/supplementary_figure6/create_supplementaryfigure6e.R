# CREATE SUPPLEMENTARY FIGURE 5K ----
# Create barplot of large dup+ samples split by HRD

# PREAMBLE ----
library(rstudioapi)
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(yaml)
library(ggplot2)
library(patchwork)
library(ggpubr)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# FUNCTIONS ----
theme_LM = theme_classic() + grids()

plot_prop_samples <- function(df, label_group, alpha_col, pivot_col, legend_name) {
  df %>%
    ggplot(aes(x=ecDNA_status, fill=!!sym(label_group), alpha=!!sym(alpha_col))) +
    geom_bar(stat='count', position='fill', colour='black', show.legend = c(fill = FALSE, alpha = TRUE)) +
    facet_grid(as.formula(paste0(pivot_col, "~", label_group))) +
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
pri_samples <- fread(file.path(main_repo_path, "data", "SupplementaryFigure5k_sourcetable.txt"))

test_prop <- function(df, type2test, group2test) {
  results <- df %>%
    filter(Type == type2test & label_group == group2test) %>%
    group_by(label_group, ecDNA_status, large_dup_positive) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = large_dup_positive,
                values_from = n) %>%
    replace_na(list(`FALSE` = 0)) %>%
    fisher_test_single_category(., c('ecDNA+', 'ecDNA-'), 'TRUE', 'FALSE') 
  
  return(data.frame(group = group2test,
                    tool = type2test,
                    p_value = results$p_value,
                    or = results$odds_ratio,
                    stringsAsFactors = FALSE))
}

stat_results <- data.frame()
for (group2test in c('ER+ High', 'ER+ Typical', 'TNBC')) {
  for (type2test in c('AmpliconArchitect', 'JaBbA')) {
    result <- test_prop(df = pri_samples,
                        type2test = type2test,
                        group2test = group2test)
    stat_results <- bind_rows(stat_results, result)
  }
}

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure5k.pdf"),height=6, width=8)
plot_prop_samples(df=pri_samples, 
                  label_group='label_group', 
                  alpha_col='large_dup_positive', 
                  pivot_col='Type', 
                  legend_name='DUP (>100Kb)')
dev.off()

fwrite(stat_results, file.path(main_repo_path, "plots", "SupplementaryFigure5k_stats.tsv"))