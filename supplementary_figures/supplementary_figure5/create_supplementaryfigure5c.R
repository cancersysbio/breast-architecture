# CREATE SUPPLEMENTARY FIGURE 5c ----
# Barplot of amplicon overlap (Metastatic)

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
group_colours <- c("ER+ High" = "#fa954eff",
                   "ER+ Typical" = "#a1c1d4ff",
                   "HER2+" = "#8b0100ff",
                   "IC10" = "#7c26ccff",
                   "IC4ER-" = "#c2b7dbff")
ic_colours <- c("11q13" = "#00ee77ff",
                "17q23" = "#ff5500ff",
                "8q24" = "#ee82edff")

sv_type_colour = c(
  "DEL" = "#912F40",
  "DUP" = "#BAD29F",
  "INV" = "#6F73D2",
  "TRA" = "#06BEE1"
)


plot_prop_her2_samples <- function(df) {
  df %>%
    ggplot(., aes(x = Cytoband_hg19, y = Proportion, alpha = Type, fill=Cytoband_hg19)) +
    geom_bar(stat = "identity", position = "stack", color="black") +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4, color = "black") +
    scale_fill_manual(values = c("Total" = "gray", "Proportion" = "skyblue")) +
    labs(title = "",
         x = "",
         y = "Proportion of HER2+ Samples") +
    scale_fill_manual(name = "", values=ic_colours) +
    scale_alpha_discrete(range = c(0.4, 1), 
                         guide = guide_legend(override.aes = list(fill = "darkgray"),
                                              title = "Altered")) +
    guides(alpha = guide_legend(override.aes = list(fill = "darkgray"),
                                title = "",
                                ncol = 1),
           fill = "none", color = "none") +
    theme_LM +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14)) +
    ylim(0,1) 
}

plot_sv_cytoband <- function(df) {
  df %>%
    ggplot(., aes(x=cytoband, y=n, fill=sv_type)) +
    geom_bar(stat='identity', position='fill') +
    geom_text(aes(label = n), position = position_fill(vjust = 0.5), size = 4, color = "black") +
    labs(title = "",
         x = "",
         y = "Proportion of SV type") +
    theme_LM +
    scale_fill_manual(name = "", values=sv_type_colour) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14))
}

theme_LM = theme_classic() + grids()

# MAIN ----
coamp_long <- fread(file.path(main_repo_path, "data", "SupplementaryFigure5c_sourcetable1.txt"))

her2_amplicon_df <- coamp_long %>%
  distinct(Sample, .keep_all = TRUE) %>%
  group_by(Cytoband_hg19) %>%
  count() %>%
  data_frame()

# SV 
samples_sv_type_in_amplicon <- fread(file.path(main_repo_path, "data", "SupplementaryFigure5c_sourcetable2.txt"))


ecdna <- fread(config$aa$gene_proj17) %>%
  filter(sample_name %in% coamp_matrix_clean$Sample) %>%
  left_join(ic_amplicon_genes, by=c("gene"="Gene_symbol")) %>%
  filter(!is.na(Cytoband_hg19))
ecdna_annot <- fread(config$aa$proj17)

her2_ecdna_df <- fread(file.path(main_repo_path, "data", "SupplementaryFigure5c_sourcetable3.txt")) %>%
  group_by(Cytoband_hg19) %>%
  distinct(sample_name) %>%
  count()

df <- her2_amplicon_df %>%
  dplyr::rename(nTRUE=n) %>%
  mutate(total = coamp_matrix_clean$Sample %>% length()) %>%
  mutate(ProportionTRUE = nTRUE/total) %>%
  mutate(nFALSE=total-nTRUE) %>%
  mutate(ProportionFALSE = nFALSE/total) %>%
  pivot_longer(cols = c(starts_with("n"), starts_with("Proportion")),
               names_to = c(".value", "Type"),
               names_pattern = "(n|Proportion)(TRUE|FALSE)") %>%
  mutate(Cytoband_hg19=factor(Cytoband_hg19, levels=c("17q23", "11q13", "8q24")))
plot_prop_her2_samples(df)

p2 <- samples_sv_type_in_amplicon %>%
  mutate(cytoband=factor(cytoband, levels=c("17q23", "11q13", "8q24"))) %>%
  group_by(cytoband) %>%
  count(sv_type) %>%
  plot_sv_cytoband(.)
p2

p1_df <- df %>%
  mutate(Type=case_when(Type == TRUE ~ 'Amplified',
                        Type == FALSE ~ 'Not Amplified'))
p3_df <- her2_ecdna_df %>%
  dplyr::rename(nTRUE=n) %>%
  mutate(total = coamp_matrix_clean$Sample %>% length()) %>%
  mutate(ProportionTRUE = nTRUE/total) %>%
  mutate(nFALSE=total-nTRUE) %>%
  mutate(ProportionFALSE = nFALSE/total) %>%
  pivot_longer(cols = c(starts_with("n"), starts_with("Proportion")),
               names_to = c(".value", "Type"),
               names_pattern = "(n|Proportion)(TRUE|FALSE)") %>%
  mutate(Cytoband_hg19=factor(Cytoband_hg19, levels=c("17q23", "11q13", "8q24"))) %>%
  mutate(Type=case_when(Type == TRUE ~ 'Amplified in ecDNA',
                        Type == FALSE ~ 'Not amplified in ecDNA')) %>%
  filter(Type != "Not amplified in ecDNA") %>%
  mutate(Cytoband_hg19=as.character(Cytoband_hg19)) %>%
  dplyr::left_join(., p1_df %>% filter(Type == 'Amplified') %>% select(Cytoband_hg19, Proportion, n),
                   by = 'Cytoband_hg19') %>%
  mutate(Proportion = Proportion.y - Proportion.x,
         n = n.y - n.x) %>%
  select(Cytoband_hg19, total, Type, n, Proportion) %>%
  mutate(Type = 'Amplified')

p3_df2 <- her2_ecdna_df %>%
  dplyr::rename(nTRUE=n) %>%
  mutate(total = coamp_matrix_clean$Sample %>% length()) %>%
  mutate(ProportionTRUE = nTRUE/total) %>%
  mutate(nFALSE=total-nTRUE) %>%
  mutate(ProportionFALSE = nFALSE/total) %>%
  pivot_longer(cols = c(starts_with("n"), starts_with("Proportion")),
               names_to = c(".value", "Type"),
               names_pattern = "(n|Proportion)(TRUE|FALSE)") %>%
  mutate(Cytoband_hg19=factor(Cytoband_hg19, levels=c("17q23", "11q13", "8q24"))) %>%
  mutate(Type=case_when(Type == TRUE ~ 'Amplified in ecDNA',
                        Type == FALSE ~ 'Not amplified in ecDNA')) %>%
  filter(Type != "Not amplified in ecDNA")

p1_p3_df <- rbind(p1_df %>% filter(Type != 'Amplified'), # samples without SV
                  p3_df) %>% # samples with SV but not ecDNA
  rbind(., p3_df2) %>%
  mutate(Type = factor(Type, levels = c('Not Amplified', 'Amplified', 'Amplified in ecDNA'))) %>%
  mutate(Cytoband_hg19=factor(Cytoband_hg19, levels=c("17q23", "11q13", "8q24")))

p_concat <- plot_prop_her2_samples(p1_p3_df)

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure5c.pdf"), height=4, width=7)
(p_concat | p2) + plot_layout(guides="collect")
dev.off()
