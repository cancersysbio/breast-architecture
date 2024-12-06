# CREATE SUPPLEMENTARY FIGURE 3I ----
# Create boxplot of SV signatures in HRD-like tumors

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

plot_sig_group_nosplit_per_signature <- function(df, grouping, hrd_label_type) {
  
  df_arc <- df %>%
    drop_na(!!sym(hrd_label_type)) %>%
    filter(!!sym(grouping) != "Remove") %>%
    select(Sample, !!sym(grouping), !!sym(hrd_label_type), all_of(signatures)) %>%
    pivot_longer(cols=4:ncol(.), names_to = 'Signature', values_to = 'Contribution') %>%
    filter(Signature == 'SBS3' | Signature == 'SBS8') %>%
    distinct(Sample, Signature, .keep_all = TRUE)
  
  sbs_sig <- df_arc %>%
    ggplot(aes(x = !!sym(grouping), y = Contribution)) +
    geom_boxplot(width = 0.9, outlier.shape = NA) + 
    geom_jitter(position = position_jitter(width = 0.15), size = 2,  alpha = 0.5) +
    facet_grid(as.formula(paste0("Signature", "~", hrd_label_type)), scales = 'free_y', space = 'free') +
    stat_compare_means(comparisons = list(c('ER+', 'TNBC')),
                       method = "wilcox.test",
                       label = "p.format",
                       size=4) +
    labs(y="Signature Activity") +
    scale_y_continuous(labels = function(x) paste(x, "%", sep = ""), 
                       expand = c(.1, .3)) +
    theme_LM + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
                     axis.text.y = element_text(size = 18),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 18),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 12),
                     title = element_text(size = 20),
                     legend.position = "top",
                     legend.key.width = unit(2, "cm"),
                     strip.text.x = element_text(size = 18),
                     strip.text.y = element_text(size = 18))
  
  # sv signature
  df_arc_sv <- df %>%
    drop_na(!!sym(hrd_label_type)) %>%
    filter(!!sym(grouping) != "Remove") %>%
    select(Sample, !!sym(grouping), !!sym(hrd_label_type), all_of(signatures)) %>%
    dplyr::rename('RS1' = 'SV32F',
                  'RS5' = 'SV32E',
                  'RS3' = 'SV32C') %>%
    pivot_longer(cols=4:ncol(.), names_to = 'Signature', values_to = 'Contribution') %>%
    filter(Signature %in% c('RS1', 'RS3', 'RS5')) %>%
    distinct(Sample, Signature, .keep_all = TRUE)
  
  sv_sig <- df_arc_sv %>%
    ggplot(aes(x = !!sym(grouping), y = Contribution)) +
    geom_boxplot(width = 0.9, outlier.shape = NA) + 
    geom_jitter(position = position_jitter(width = 0.15), size = 2,  alpha = 0.5) +
    facet_grid(as.formula(paste0("Signature", "~", hrd_label_type)), scales = 'free_y', space = 'free') +
    stat_compare_means(comparisons = list(c('ER+', 'TNBC')),
                       method = "wilcox.test",
                       label = "p.format",
                       size=4) +
    labs(y="Signature Activity", color = 'TNBC Archetype') +
    theme_LM + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
                     axis.text.y = element_text(size = 18),
                     axis.title.x = element_blank(),
                     axis.title.y = element_text(size = 18),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size = 12),
                     title = element_text(size = 20),
                     legend.position = "top",
                     legend.key.width = unit(2, "cm"),
                     strip.text.x = element_text(size = 18),
                     strip.text.y = element_text(size = 18)) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^floor(x)),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = expansion(mult = c(0.2, 0.2))
    )
  (sbs_sig | sv_sig)
  
}


# MAIN ----
sig_df <- fread(file.path(main_repo_path, "data", "SupplementaryFigure3hi_sourcetable.txt")) %>%
  filter(chord_label %in% c("non-HRD-like", "BRCA1-like", "BRCA2-like"))

pdf(file.path(main_repo_path, "plots", "SupplementaryFigure3h"),height=10, width=8)
plot_sig_group_nosplit_per_signature(df = sig_df,
                                     grouping = 'group_label_type1',
                                     hrd_label_type = 'chord_label')
dev.off()
