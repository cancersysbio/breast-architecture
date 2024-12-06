# CREATE SUPPLEMENTARY FIGURE 5J ----
# Boxplot of Type I IFN in HER2+ and ER+ High Risk, separated by ecDNA+/-

# PREAMBLE ----
library(rstudioapi)
library(tidyr)
library(data.table)
library(dplyr)
library(yaml)
library(ggplot2)
library(ggpubr)
library(stringr)
library(patchwork)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# FUNCTIONS ----
ecdna_colours <- c("HER2+"="#8b0100ff", 
                   "IC1"="#ff5500ff",
                   "IC6"="#fffe41ff",
                   "IC9"="#ee82edff",
                   "IC10"="#7c26ccff")
ic_subgroup_colours <- c("ER+ High"="#fa954eff", "HER2+"="#8b0100ff", "IC10"="#7c26ccff")
theme_LM = theme_classic() + grids()

# MAIN ----
wgs_mtable_rep <- fread(file.path(main_repo_path, "data", "SupplementaryFigure5j_sourcetable.txt"))

p1 <- wgs_mtable_rep %>%
  filter(Cohort != "Hartwig") %>%
  filter(label != "ic2") %>%
  mutate(ecdna_label = case_when(ecdna==TRUE ~ "ecDNA+",
                                 ecdna==FALSE ~ "ecDNA-"),
         label = toupper(label)) %>%
  mutate(label=factor(label, levels=c("HER2+", "IC1", "IC6", "IC9", "IC10"))) %>%
  ggboxplot(x = "ecdna_label", 
            y = "TypeI_IFN", 
            add = "jitter", 
            color = "label",
            facet.by = "label", 
            nrow = 1) +
  labs(y="Type I IFN Signature") +
  theme_LM + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
                   axis.text.y = element_text(size = 18),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 18),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 20),
                   legend.position = "none",
                   strip.text.x = element_text(size = 22)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_color_manual(values = ecdna_colours)


her2 <- wgs_mtable_rep %>%
  filter(label=="HER2+") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~  TypeI_IFN + ERBB2 + Cohort, data=her2, family=binomial())
her2_df <- data.frame(group="her2+",
                      sig_estimate=summary(model)$coefficient[2,1],
                      pval=summary(model)$coefficient[2,4])

ic1 <- wgs_mtable_rep %>%
  filter(label=="ic1") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ TypeI_IFN + PPM1D + MSI2 + Cohort, data=ic1, family=binomial())
ic1_df <- data.frame(group="ic1",
                     sig_estimate=summary(model)$coefficient[2,1],
                     pval=summary(model)$coefficient[2,4])

ic6 <- wgs_mtable_rep %>%
  filter(label=="ic6") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ TypeI_IFN + FGFR1 + Cohort, data=ic6, family = binomial())
ic6_df <- data.frame(group="ic6",
                     sig_estimate=summary(model)$coefficient[2,1],
                     pval=summary(model)$coefficient[2,4])

ic9 <- wgs_mtable_rep %>%
  filter(label=="ic9") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ TypeI_IFN + MYC + Cohort, data=ic9, family = binomial())
ic9_df <- data.frame(group="ic9",
                     sig_estimate=summary(model)$coefficient[2,1],
                     pval=summary(model)$coefficient[2,4])

ic10 <- wgs_mtable_rep %>%
  filter(label=="IC10") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ TypeI_IFN + ERBB2 + Cohort, data=ic10, family=binomial())
ic10_df <- data.frame(group="ic10",
                      sig_estimate=summary(model)$coefficient[2,1],
                      pval=summary(model)$coefficient[2,4])

ss_ecdna_ifn <- rbind(her2_df, ic1_df, ic6_df, ic9_df, ic10_df)
annotations_ifn <- data.frame(
  label = factor(c("HER2+", "IC1", "IC6", "IC9", "IC10"), levels = c("HER2+", "IC1", "IC6", "IC9", "IC10")),
  x = rep(0, 5),
  y = rep(0.55, 5),
  text = c(paste0('ES=', signif(ss_ecdna_ifn$sig_estimate[1], 3), '\n', "P=", signif(ss_ecdna_ifn$pval[1], 3)), 
           paste0('ES=', signif(ss_ecdna_ifn$sig_estimate[2], 3), '\n', "P=", signif(ss_ecdna_ifn$pval[2], 3)), 
           paste0('ES=', signif(ss_ecdna_ifn$sig_estimate[3], 3), '\n', "P=", signif(ss_ecdna_ifn$pval[3], 3)), 
           paste0('ES=', signif(ss_ecdna_ifn$sig_estimate[4], 3), '\n', "P=", signif(ss_ecdna_ifn$pval[4], 3)), 
           paste0('ES=', signif(ss_ecdna_ifn$sig_estimate[5], 3), '\n', "P=", signif(ss_ecdna_ifn$pval[5], 3))
  )
)

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure5j.pdf"), height=6, width=9)
p1 + geom_text(data = annotations_ifn, aes(x = x + 0.5, y = y, label = text),
               color = "black", size = 5, hjust = 0, vjust = 0)
dev.off()

fwrite(ss_ecdna_ifn, paste0(config$output_dir, "SupplementaryFigure5j_stats.tsv"))