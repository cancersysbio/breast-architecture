# CREATE SUPPLEMENTARY FIGURE 1I ----
# Boxplot of purity IC subgroups

# PREAMBLE ----
library(rstudioapi)
library(dplyr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(yaml)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# FUNCTIONS ----
theme_LM = theme_classic() + grids()
add_alpha <- function(ggplt) {
  idx <- which(sapply(ggplt$layers, function(l) "PositionJitter" %in% class(l$position)))
  ggplt$layers[[idx]]$aes_params$alpha <- 0.3
}

group_colors <- c("ER+ High" = "#fa954eff", 
                  "ER+ Typical"="#a1c1d4ff", 
                  "HER2+" = "#bd854dff", 
                  "IC10" = "#904dbdff", 
                  "IC4ER-" = "#c2b7dbff")

plot_purity_boxplot <- function(df) {
  df %>%
    ggboxplot(., "group", "value", add = "jitter", xlab = "", ylab = "",
              color = "group", palette = group_colors,
              facet.by = "name", scales = "free", nrow=1, 
              panel.labs = list(name = c("Purity"))) +
    theme_LM + 
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "none",
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14)) +
    rotate_x_text(35) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
}

# MAIN ----
df <- fread(file.path(main_repo_path, "data", "SupplementaryFigure1i_sourcetable.txt"))
p1 <- df %>%
  pivot_longer(cols=c("purity")) %>%
  plot_purity_boxplot(.)
add_alpha(p1)
ss_purity_all_adj <- compare_means(purity ~ group, data=df) %>%
  select(-c(p.format, p.signif))

# SAVE ----
pdf(file.path(main_repo_path, "data", "SupplementaryFigure1i.pdf"), height=4, width=3)
p1
dev.off()

fwrite(ss_purity_all_adj, file.path(main_repo_path, "data", "SupplementaryFigure1i_stats.tsv"))
