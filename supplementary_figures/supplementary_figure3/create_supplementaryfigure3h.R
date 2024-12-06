# CREATE SUPPLEMENTARY FIGURE 3H ----
# Create boxplot of ternary archetype proportion in HRD-like tumors

# PREAMBLE ----
library(rstudioapi)
library(tidyr)
library(data.table)
library(dplyr)
library(yaml)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(broom)
library(scales)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))


# FUNCTIONS ----
theme_LM = theme_classic() + grids()

plot_arc_group_arc_nosplit <- function(df, grouping, hrd_label_type) {
  df_arc <- df %>%
    drop_na(!!sym(hrd_label_type)) %>%
    filter(!!sym(grouping) != "Remove") %>% # plot only TNBC and ER+
    filter(!is.na(Arc1))  %>%
    distinct(Sample, Arc1, .keep_all = TRUE)
  
  comps <- df_arc %>% distinct(!!sym(hrd_label_type)) %>% pull(!!sym(hrd_label_type)) %>% sort()
  stat_test <- as.list(combn(comps, 2, simplify = FALSE))
  
  df_arc %>%
    ggplot(aes(x = !!sym(hrd_label_type), y = Arc1)) +
    geom_boxplot(width = 0.9, outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.15), size = 2,  alpha = 0.5) +
    stat_compare_means(comparisons = stat_test,
                       method = "wilcox.test",
                       aes(label = paste0("p=", format(..p.val, digits = 3))),
                       size = 4) +
    facet_wrap(as.formula(paste0("~", grouping)), ncol = 2) +
    labs(y="TNBC-Enriched Archetype") +
    scale_y_continuous(expand = c(0, .1)) +
    theme_LM + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
                     axis.text.y = element_text(size = 18),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 18),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 12),
                     title = element_text(size = 20),
                     legend.position = "top",
                     legend.key.width = unit(2, "cm"),
                     strip.text.x = element_text(size = 22),
                     strip.text.y = element_text(size = 18))
}


# MAIN ----
sig_df <- fread(file.path(main_repo_path, "data", "SupplementaryFigure3hi_sourcetable.txt"))

pdf(file.path(main_repo_path, "plots", "SupplementaryFigure3h.pdf"), height=5, width=4)
plot_arc_group_arc_nosplit(df = sig_df,
                           grouping = 'group_label_type1', 
                           hrd_label_type = 'chord_label')
dev.off()
