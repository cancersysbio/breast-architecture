### CREATE EXTENDED DATA FIGURE 2DE #############################################################################
# Boxplot of genomic alterations in IC10 vs IC4ER- (ED2d) &
# Barplot of pathways altered in IC10 vs IC4ER- (ED2e)

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(yaml)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### FUNCTIONS #####################################################################################
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

plot_boxplot_tnbc <- function(df) {
  df %>%
  ggboxplot(., "group", "value", add = "jitter", xlab = "", ylab = "",
            color = "group", palette = group_colors,
            facet.by = "name", scales = "free", nrow=1, 
            panel.labs = list(name = c("Fraction \n genome altered", "Fraction LOH", "No. of \n damaging SVs"))) +
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

plot_barplot_tnbc <- function(prop_df) {
  prop_df %>%
    ggplot(aes(x=group, fill=group, alpha=value)) +
    geom_bar(position="fill", stat="count", color="black") + 
    stat_count(geom = "text", 
               position=position_fill(vjust = 0.5),
               aes(label = ..count..),
               color="black", 
               size=4) +
    facet_wrap(~name, labeller = labeller(name = c(cell_cycle = "Cell cycle", ddr = "DDR", ubi = "Ubiquitination"))) +
    labs(y="Proportion of samples \n with alteration in pathway", x="") + 
    scale_fill_manual(values = c('#7c26ccff','#c2b7dbff')) +
    scale_alpha_discrete(range = c(0.4, 1.5), 
                         guide = guide_legend(override.aes = list(fill = "darkgray"),
                                              title = "Altered")) +
    guides(alpha = guide_legend(override.aes = list(fill = "darkgray"),
                                title = "Pathway\naltered",
                                ncol = 1),
           fill = "none", color = "none") +
    theme_LM +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 11), 
          strip.text.y = element_text(size = 11)) +
    ylim(0,1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
}

### MAIN ##########################################################################################
## GENERATE EXTENDED DATA 2D ----
tnbc_genomic_df <- read.delim(file.path(main_repo_path, "data", "ExtendedData2de_sourcetable_primary.txt"), sep = ",")

p2 <- plot_boxplot_tnbc(tnbc_genomic_df)
add_alpha(p2)

stat_test <- tnbc_genomic_df %>%
  group_by(name) %>%
  wilcox_test(value ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "group")

ss_tnbc <- stat_test %>%
  select(name, group1, group2, p.adj)

eff_sizes <- tnbc_genomic_df %>%
  group_by(name) %>%
  wilcox_effsize(value ~ group) %>%
  select(-c(n1,n2,magnitude))

ss_tnbc_save <- ss_tnbc %>%
  left_join(eff_sizes, by=c("name", "group1", "group2")) %>%
  select(-c(".y."))

## GENERATE EXTENDED DATA 2E ----
primary_w_mutations_long_tnbc <- read.delim(file.path(main_repo_path, "data", "ExtendedData2de_sourcetable_primary_pathways.txt"), sep = ",")

primary_w_mutations_long_tnbc_prop <- primary_w_mutations_long_tnbc %>%
  group_by(group, name, value) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  mutate(value=case_when(value == TRUE ~ "Altered",
                         value == FALSE ~ "Not altered")) %>%
  mutate(value=factor(value, levels=c("Altered", "Not altered")))

pathways <- c('ddr', 'ubi', 'cell_cycle')
ss_tnbc_pathways <- list()
for (i in 1:length(pathways)) {
  pathway <- pathways[i]
  
  count_pathway_temp <- primary_w_mutations_long %>% 
    filter(name == pathway) %>% 
    filter(group %in% c('IC10', 'IC4ER-')) %>% 
    count(group, value)
  
  df <- data.frame(
    "IC10" = c(count_pathway_temp %>% filter(group == "IC10" & value == TRUE) %>% pull(n), count_pathway_temp %>% filter(group == "IC10" & value == FALSE) %>% pull(n)), 
    "IC4ER-" = c(count_pathway_temp %>% filter(group == "IC4ER-" & value == TRUE) %>% pull(n), count_pathway_temp %>% filter(group == "IC4ER-" & value == FALSE) %>% pull(n)), 
    row.names = c("TRUE", "FALSE"))
  
  ss_tnbc_pathways[[i]] <- data.frame(
    p=fisher.test(df)$p,
    or=fisher.test(df)$estimate[[1]],
    pathway=pathway
  )
}

ss_tnbc_pathways_all <- rbind(ss_tnbc_pathways[[1]], ss_tnbc_pathways[[2]], ss_tnbc_pathways[[3]])

p_adj_tnbc <- ss_tnbc_pathways_all %>%
  pull(p) %>%
  p.adjust(method = "fdr")

ss_tnbc_all_adj <- ss_tnbc_pathways_all %>%
  mutate(p.adj = p_adj_tnbc) %>%
  mutate(method="fisher.test")

pdf(file.path(main_repo_path, "plots", "ExtendedData2d.pdf"), height=4, width=6)
p2
dev.off()
ss_tnbc_save %>%
  write.table(file.path(main_repo_path, "plots", "ExtendedData2d_stats.txt"), 
              quote = F, row.names = F, sep = "\t")

pdf(file.path(main_repo_path, "plots", "ExtendedData2e.pdf"), height=4, width=5)
plot_barplot_tnbc(primary_w_mutations_long %>% filter(group %in% c("IC10", "IC4ER-")))
dev.off()
ss_tnbc_all_adj %>%
  write.table(file.path(main_repo_path, "plots", "ExtendedData2e_stats.txt"), 
              quote = F, row.names = F, sep = "\t")
