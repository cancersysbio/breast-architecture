# CREATE SUPPLEMENTARY FIGURE 4H ----
# Barplot of WGD in IC subgroups, pri and met

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
library(scales)

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

plot_barplot_wgd <- function(df) {
  df %>%
    ggplot(aes(x=Stage, alpha=genome_doubled, fill=group)) + 
    geom_bar(position="fill", stat="count", color="black") + 
    stat_count(geom = "text", 
               position=position_fill(vjust = 0.5),
               aes(label = ..count..),
               color="black", 
               size=4) +
    facet_wrap(~group, ncol = 5) +
    scale_alpha_discrete(range = c(0.4, 1), 
                         guide = guide_legend(override.aes = list(fill = "darkgray"),
                                              title = "WGD")) +
    scale_fill_manual(values = group_colors) +
    guides(alpha = guide_legend(override.aes = list(fill = "darkgray"),
                                title = "WGD",
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
    labs(y="Proportion of samples", x="") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}

# GENERATE SUPPLEMENTARY FIGURE 4H ----
all_mt <- fread(file.path(main_repo_path, "data", "SupplementaryFigure4h_sourcetable.txt")) %>%
  mutate(Stage=factor(Stage, levels=c("Primary", "Metastatic")))

plt <- all_mt %>%
  drop_na(genome_doubled, group) %>%
  select(Stage, group, genome_doubled) %>%
  plot_barplot_wgd(.)

count_wgd <- all_mt %>%
  drop_na(genome_doubled, group) %>%
  count(Stage, group, genome_doubled) %>%
  group_by(group) %>%
  pivot_wider(names_from = genome_doubled, values_from = n)

groups <- c("ER+ High", "ER+ Typical", "HER2+", "IC10", "IC4ER-")

ss_wgd <- list()
for (i in 1:length(groups)) {
  c_group <- groups[i]
  
  count_wgd_temp <- count_wgd %>%
    filter(group==c_group)
  
  xtab <- as.table(rbind(
    c(count_wgd_temp[[1,3]], count_wgd_temp[[1,4]]),
    c(count_wgd_temp[[2,3]], count_wgd_temp[[2,4]])
  ))
  
  dimnames(xtab) <- list(
    Stage = c("Primary", "Metastatic"),
    WGD = c("No", "Yes")
  )
  
  ss_wgd[[i]] <- pairwise_fisher_test(xtab)
}

ss_wgd_all <- rbind(ss_wgd[[1]], ss_wgd[[2]], ss_wgd[[3]],
                    ss_wgd[[4]], ss_wgd[[5]]) %>%
  mutate(group = groups)

p_adj <- ss_wgd_all %>%
  pull(p) %>%
  p.adjust(method = "fdr")

ss_wgd_all_adj <- ss_wgd_all %>%
  mutate(p.adj = p_adj) %>%
  select(-c(p.adj.signif)) %>%
  mutate(method="fisher.test")

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure4h.pdf"), height=4, width=8)
plt + 
  geom_text(data = ss_wgd_all_adj, aes(x = 1, y = 1.0, label = paste0('FDR=', scientific(p.adj, digits = 3))), hjust = .2, vjust = -1, size = 4, inherit.aes = FALSE)
dev.off()

fwrite(ss_wgd_all_adj, file.path(main_repo_path, "plots", "SupplementaryFigure4h_stats.tsv"))
