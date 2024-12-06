# CREATE SUPPLEMENTARY FIGURE 2f ----
# Create scatterplot of RS4 and RS6 by SV type

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

# FUNCTIONS ----
theme_LM = theme_classic() + grids()

colour_palette <- c(
  "BFB" = "#f9a41bff",
  "Chromoplexy" = "#428c40ff",
  "Chromothripsis" = "#682c8bff",
  "CPXDM" = "#fbd605ff",
  "DEL" = "#aa1a91ff",
  "DM" = "#4a86c5ff",
  "DUP" = "#63aaafff",
  "INV" = "#000000ff",
  "INVDUP" = "#cf2528ff",
  "Pyrgo" = "#9aca3bff",
  "Rigma" = "#f47f21ff",
  "TIC" = "#463c8dff",
  "TRA" = "#55bc7cff",
  "Tyfona" = "#ce2f79ff"
)

group_shapes <- c(
  "ER+ High" = 17,
  "HER2+" = 18,
  "ER+ Typical" = 15,
  'IC10' = 3
) 

plot_scatter_group <- function(df) {
  df %>%
    pivot_wider(names_from = 'RS', values_from = 'correlation') %>%
    ggscatter(., x = "RS4", y = "RS6", color = "Type", shape = 'group',
              legend.title = "name",
              ylab = "Correlation with RS6",
              xlab = "Correlation with RS4",
              title = "Per subgroup") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    # stat_cor(method = "spearman") +
    stat_cor() +
    xlim(0, 1) + 
    ylim(0, 1) +
    guides(color = guide_legend(ncol = 2), 
           shape = guide_legend(nrow = 5)) +
    theme_LM +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14),
          legend.position = "bottom",
          legend.title = element_blank()) +
    guides(color = guide_legend(ncol = 2)) +
    scale_color_manual(values = colour_palette) +
    scale_shape_manual(values = group_shapes)
}

plot_scatter_all <- function(df) {
  df %>%
    pivot_wider(names_from = 'RS', values_from = 'correlation') %>%
    ggscatter(., x = "RS4", y = "RS6", color = "Type",
              legend.title = "name",
              ylab = "Correlation with RS6",
              xlab = "Correlation with RS4",
              title = "All subgroups") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    # stat_cor(method = "spearman") +
    stat_cor() +
    xlim(0, 1) + 
    ylim(0, 1) +
    theme_LM +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14),
          legend.position = "bottom",
          legend.title = element_blank()) +
    guides(color = guide_legend(ncol = 2)) +
    scale_color_manual(values = colour_palette) 
}

# MAIN ----
input_df <- fread(file.path(main_repo_path, "data", "SupplementaryFigure2fg_sourcetable.txt"))
jabba_types <- c(
  'BFB', 'Chromoplexy', 'Chromothripsis', 'CPXDM',
  'DEL', 'DM', 'DUP', 'INV', 'INVDUP', 'Pyrgo', 
  'Rigma', 'TIC', 'TRA', 'Tyfona'
)

rs_cor <- list()
for (rs in c('RS4', 'RS6')) {
  for (jabba_type in jabba_types) {
    correlation <- input_df[!is.na(group), .(
      correlation = abs(cor(get(rs), get(jabba_type), use = 'complete.obs'))
    ), by = group]
    
    correlation[, RS := rs]
    correlation[, Type := jabba_type]
    
    rs_cor[[rs]][[jabba_type]] <- correlation
  }
}

df <- rbindlist(unlist(rs_cor, recursive = FALSE), use.names = TRUE, fill = TRUE)
p_bygroup <- plot_scatter_group(df)

rs_cor_grouped <- list()
for (rs in c('RS4', 'RS6')) {
  for (jabba_type in jabba_types) {
    correlation <- input_df[, .(
      correlation = abs(cor(get(rs), get(jabba_type), use = 'complete.obs'))
    )]
    
    correlation[, RS := rs]
    correlation[, Type := jabba_type]
    
    rs_cor_grouped[[rs]][[jabba_type]] <- correlation
  }
}

df_grouped <- rbindlist(unlist(rs_cor_grouped, recursive = FALSE), use.names = TRUE, fill = TRUE)
p_all <- plot_scatter_all(df_grouped)

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure2f.pdf"), width=7, height=8)
(p_all / p_bygroup) + plot_layout(guides="collect")
dev.off()
