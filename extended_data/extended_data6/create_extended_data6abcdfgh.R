### CREATE EXTENDED DATA FIGURE 6 #############################################################################
# Create boxplot of replication stress signature in TCGA/Metabric, by subgroup or by subtype, and
# replication stress in HER2+ and ER+ High ecDNA+/-

### PREAMBLE #####################################################################################
library(data.table)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(yaml)
library(grid)
library(metafor)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### FUNCTIONS ##########################################################################################
theme_LM = theme_classic() + grids()
colours = c("#fa954eff", "#a1c1d4ff", "#8b0100ff", "#8b0100ff", "#7c26ccff", "#c2b7dbff", "#a3cfffff", "#335e8dff")
ic_colours <- c("IC1"="#ff5500ff", "IC2"="#00ee77ff", "IC6"="#fffe41ff",
                "IC9"="#ee82edff", "IC5"="#8b0100ff", "IC3"="#cd32caff", "IC3/IC7"="#cd32caff",
                "IC8"="#feaa00ff", "IC10"="#7c26ccff", "IC4ER+"="#00c4cdff", "IC4ER-"="#c2b7dbff")

ecdna_colours <- c("HER2+"="#8b0100ff", 
                   "IC1"="#ff5500ff",
                   "IC6"="#fffe41ff",
                   "IC9"="#ee82edff",
                   "IC10"="#7c26ccff")
ic_subgroup_colours <- c("ER+ High"="#fa954eff", "HER2+"="#8b0100ff", "IC10"="#7c26ccff")

add_alpha <- function(ggplt) {
  idx <- which(sapply(ggplt$layers, function(l) "PositionJitter" %in% class(l$position)))
  ggplt$layers[[idx]]$aes_params$alpha <- 0.3
}

### MAIN ##########################################################################################
## EXTENDED DATA 6A  ----
group1 <- c("ER+ High", "ER+ Typical", "HER2+/ER-", "HER2+/ER+","IC10", "IC4ER-")
group2 <- c("IDC", "ILC")
group3 <- c("IC1", "IC2", "IC6", "IC9", 
            "IC3", "IC7", "Other", "IC4ER+","IC8",
            "IC5",
            "IC10", "IC4ER-")

rna_sigs <- read.delim(file.path(main_repo_path, 'data', 'ExtendedData6_sourcetable1.txt'))

tcga_eniclust <- rna_sigs %>%
  mutate(group=na_if(group, "")) %>%
  filter(Cohort=="TCGA") %>%
  drop_na(Replication_stress_sig, group) %>%
  mutate(label=group) %>%
  mutate(label=case_when(group == "HER2+" & ER == 1 ~ "HER2+/ER+",
                         group == "HER2+" & ER == 0 ~ "HER2+/ER-",
                         T ~ group))

tcga_histo <- rna_sigs %>%
  filter(Cohort=="TCGA") %>%
  drop_na(Replication_stress_sig, group) %>%
  mutate(label=case_when(group == "ER+ Typical" & HistoType == "ILC" ~ "ILC",
                         group == "ER+ Typical" & HistoType == "IDC" ~ "IDC",
                         T ~ group)) %>%
  filter(label %in% group2)

tcga <- rbind(tcga_eniclust, tcga_histo) %>%
  mutate(label=factor(label, levels=c("ER+ High", "ER+ Typical", "HER2+/ER-", "HER2+/ER+", 
                                      "IC10", "IC4ER-",
                                      "IDC", "ILC"))) %>%
  mutate(plotting_group=case_when(label %in% group1 ~ "IC subgroup",
                                  label %in% group2 ~ "ER+ Typical")) %>%
  mutate(plotting_group=factor(plotting_group, levels=c("IC subgroup", "ER+ Typical")))

p1 <- ggboxplot(tcga, x = "label", y = "Replication_stress_sig", add = "jitter", color = "label", palette=colors) +
  labs(title="Primary (TCGA)", y="Replication stress signature") +
  theme_LM + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
                   axis.text.y = element_text(size = 12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 14),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 16),
                   legend.position = "none")
add_alpha(p1)
plt1a <- facet(p1, facet.by = "plotting_group", scales = "free", space = "free") +
  theme(strip.text.x = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))

## EXTENDED DATA 6C  ----
tcga_sig1_ic <- tcga %>%
  select(Replication_stress_sig, eniclust, ER, plotting_group) %>%
  mutate(ic10ER=case_when(eniclust == "IC4" & ER == 1 ~ "IC4ER+",
                          eniclust == "IC4" & ER == 0 ~ "IC4ER-",
                          T ~ eniclust)) %>%
  mutate(ic10ER=toupper(ic10ER)) %>%
  mutate(label=ic10ER) %>%
  mutate(label=case_when(label=="OTHER" ~ "IC3/IC7",
                         T ~ label)) %>%
  mutate(label=factor(label, levels=c("IC1", "IC2", "IC6", "IC9", 
                                      "IC3", "IC7", "IC8", "IC3/IC7", "IC4ER+",
                                      "IC5",
                                      "IC10", "IC4ER-")))

ss_tcga_ic <- compare_means(Replication_stress_sig ~ label, group.by = "plotting_group", data=tcga_sig1_ic) %>%
  adjust_pvalue(method = "fdr") %>%
  mutate(p.adj.format=signif(p.adj, 3)) %>%
  mutate(p.adj.true = case_when(p.adj < 0.05 ~ T, 
                                T ~ F))

p2 <- ggboxplot(tcga_sig1_ic, x = "label", y = "Replication_stress_sig", add = "jitter", color="label", palette=ic_colours) +
  labs(title="Primary (TCGA)", y="Replication stress signature") +
  theme_LM + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
                   axis.text.y = element_text(size = 12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 14),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 16),
                   legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4)))

add_alpha(p2)
plt1b <- p2
plt1b

## EXTENDED DATA 6B  ----
rep_stress_mb <- rna_sigs %>%
  filter(Cohort=="METABRIC")

mb_eniclust <- rep_stress_mb %>%
  drop_na(Replication_stress_sig, group) %>%
  mutate(label=group) %>%
  mutate(label=case_when(group == "HER2+" & ER == 1 ~ "HER2+/ER+",
                         group == "HER2+" & ER == 0 ~ "HER2+/ER-",
                         T ~ group))

mb_histo <- rep_stress_mb %>%
  drop_na(Replication_stress_sig, group) %>%
  mutate(label=case_when(group == "ER+ Typical" & HistoType == "ILC" ~ "ILC",
                         group == "ER+ Typical" & HistoType == "IDC" ~ "IDC",
                         T ~ group)) %>%
  filter(label %in% c("ILC", "IDC"))

mb <- rbind(mb_eniclust, mb_histo) %>%
  mutate(label=factor(label, levels=c("ER+ High", "ER+ Typical", "HER2+/ER-", "HER2+/ER+", 
                                      "IC10", "IC4ER-",
                                      "IDC", "ILC"))) %>%
  mutate(plotting_group=case_when(label %in% group1 ~ "IC subgroup",
                                  label %in% group2 ~ "ER+ Typical")) %>%
  mutate(plotting_group=factor(plotting_group, levels=c("IC subgroup", "ER+ Typical"))) %>%
  drop_na(plotting_group)

ss_mb_eniclust <- mb %>% 
  group_by(plotting_group) %>% 
  compare_means(Replication_stress_sig ~ label, data=.) %>%
  adjust_pvalue(method = "fdr") %>%
  add_xy_position(data = tcga, formula = Replication_stress_sig ~ label, y.trans = NULL, scales = "free") %>%
  filter(group1 == "ER+ High" & group2 == "ER+ Typical" |
           group1 == "ER+ Typical" & group2 == "IC10" |
           group1 == "ER+ Typical" & group2 == "IC4ER-" |
           group1 == "ER+ Typical" & group2 == "HER2+/ER+" |
           group1 == "ER+ Typical" & group2 == "HER2+/ER-" |
           group1 == "ER+ High" & group2 == "HER2+/ER+" |
           group1 == "ER+ High" & group2 == "HER2+/ER-" |
           group1 == "IC10" & group2 == "IC4ER-" |
           group1 == "IDC" & group2 == "ILC") %>% # select the comparisons we want to highlight 
  mutate(y.position = c(1.1, 2.1, 2.3, 1.3, 1.7, 1.5, 1.9, 1.1, 1.1)) %>% # manually adjust position of labels
  mutate(p.adj.format=paste0("q = ",signif(p.adj, 3))) %>%
  mutate(xmin = ifelse(row_number() == 9, 1, xmin),
         xmax = ifelse(row_number() == 9, 2, xmax)) %>%
  mutate(plotting_group = c(rep("IC subgroup", 8), "ER+ Typical")) %>%
  mutate(plotting_group = factor(plotting_group, levels=c("IC subgroup", "ER+ Typical")))

p4 <- mb %>%
  ggboxplot(x = "label", y = "Replication_stress_sig", add = "jitter", color = "label", palette=colours) +
  labs(title="Primary (METABRIC)", y="Replication stress signature") +
  theme_LM + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
                   axis.text.y = element_text(size = 12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 14),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 16),
                   legend.position = "none")
add_alpha(p4)

plt2a <- facet(p4, facet.by = "plotting_group", scales = "free") +
  theme(strip.text.x = element_text(size = 14)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4)))

## EXTENDED DATA 6D  ----
mb_sig1_ic <- mb %>%
  select(Replication_stress_sig, eniclust, ER, plotting_group) %>%
  mutate(ic10ER=case_when(eniclust == "IC4" & ER == 1 ~ "IC4ER+",
                          eniclust == "IC4" & ER == 0 ~ "IC4ER-",
                          T ~ eniclust)) %>%
  mutate(ic10ER=toupper(ic10ER)) %>%
  mutate(label=ic10ER) %>%
  mutate(label=case_when(label=="OTHER" ~ "IC3/IC7",
                         T ~ label)) %>%
  mutate(label=factor(label, levels=c("IC1", "IC2", "IC6", "IC9", 
                                      "IC3", "IC7", "IC8", "IC3/IC7", "IC4ER+",
                                      "IC5",
                                      "IC10", "IC4ER-")))

ss_mb_ic <- compare_means(Replication_stress_sig ~ eniclust, group.by = "plotting_group", data=mb) %>%
  adjust_pvalue(method = "fdr") %>%
  mutate(p.adj.format=signif(p.adj, 3)) %>%
  mutate(p.adj.true = case_when(p.adj < 0.05 ~ T, 
                                T ~ F))

p5 <- mb_sig1_ic %>%
  ggboxplot(., x = "label", y = "Replication_stress_sig", add = "jitter", color="label", palette=ic_colours) +
  labs(title="Primary (Metabric)", y="Replication stress signature") +
  theme_LM + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
                   axis.text.y = element_text(size = 12),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 14),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 16),
                   legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
  stat_pvalue_manual(ss_mb_eniclust, tip.length = 0.01, step.increase = 0,
                     label = "p.adj.format", size = 3.3)

add_alpha(p5)
plt2b <- p5

## EXTENDED DATA 6G ----
tcga_eniclust_cgas <- rna_sigs %>%
  mutate(group=na_if(group, "")) %>%
  filter(Cohort=="TCGA") %>%
  drop_na(cGAS_STING_sig, group) %>%
  mutate(label=group) %>%
  mutate(label=case_when(group == "HER2+" & ER == 1 ~ "HER2+/ER+",
                         group == "HER2+" & ER == 0 ~ "HER2+/ER-",
                         T ~ group))

tcga_histo_cgas <- rna_sigs %>%
  filter(Cohort=="TCGA") %>%
  drop_na(cGAS_STING_sig, group) %>%
  mutate(label=case_when(group == "ER+ Typical" & HistoType == "ILC" ~ "ILC",
                         group == "ER+ Typical" & HistoType == "IDC" ~ "IDC",
                         T ~ group)) %>%
  filter(label %in% group2)

tcga_cgas <- rbind(tcga_eniclust_cgas, tcga_histo_cgas) %>%
  mutate(label=factor(label, levels=c("ER+ High", "ER+ Typical", "HER2+/ER-", "HER2+/ER+", 
                                      "IC10", "IC4ER-",
                                      "IDC", "ILC"))) %>%
  mutate(plotting_group=case_when(label %in% group1 ~ "IC subgroup",
                                  label %in% group2 ~ "ER+ Typical")) %>%
  mutate(plotting_group=factor(plotting_group, levels=c("IC subgroup", "ER+ Typical")))

p1_cgas <-
  tcga_cgas %>%
  ggboxplot(x = "label", y = "cGAS_STING_sig", add = "jitter", color = "label", palette=colours) +
  labs(title="Primary (TCGA)", y="cGAS-STING signature") +
  theme_LM + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
                   axis.text.y = element_text(size = 18),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 18),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 20),
                   legend.position = "none")
add_alpha(p1_cgas)

plt1a_cgas <- facet(p1_cgas, facet.by = "plotting_group", scales = "free", space = "free") +
  theme(strip.text.x = element_text(size = 26)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)))

## EXTENDED DATA 6H ----
cgas_mb <- rna_sigs %>%
  filter(Cohort=="METABRIC")

mb_eniclust_cgas <- cgas_mb %>%
  drop_na(cGAS_STING_sig, group) %>%
  mutate(label=group) %>%
  mutate(label=case_when(group == "HER2+" & ER == 1 ~ "HER2+/ER+",
                         group == "HER2+" & ER == 0 ~ "HER2+/ER-",
                         T ~ group))

mb_histo_cgas <- cgas_mb %>%
  drop_na(cGAS_STING_sig, group) %>%
  mutate(label=case_when(group == "ER+ Typical" & HistoType == "ILC" ~ "ILC",
                         group == "ER+ Typical" & HistoType == "IDC" ~ "IDC",
                         T ~ group)) %>%
  filter(label %in% c("ILC", "IDC"))

mb_cgas <- rbind(mb_eniclust_cgas, mb_histo_cgas) %>%
  mutate(label=factor(label, levels=c("ER+ High", "ER+ Typical", "HER2+/ER-", "HER2+/ER+", 
                                      "IC10", "IC4ER-",
                                      "IDC", "ILC"))) %>%
  mutate(plotting_group=case_when(label %in% group1 ~ "IC subgroup",
                                  label %in% group2 ~ "ER+ Typical")) %>%
  mutate(plotting_group=factor(plotting_group, levels=c("IC subgroup", "ER+ Typical"))) %>%
  drop_na(plotting_group)

p4_cgas <- mb_cgas %>%
  ggboxplot(x = "label", y = "cGAS_STING_sig", add = "jitter", color = "label", palette=colours) +
  labs(title="Primary (METABRIC)", y="cGAS-STING signature") +
  theme_LM + theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
                   axis.text.y = element_text(size = 18),
                   axis.title.x = element_blank(),
                   axis.title.y = element_text(size = 18),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 20),
                   legend.position = "none")
add_alpha(p4_cgas)

plt2a_cgas <- facet(p4_cgas, facet.by = "plotting_group", scales = "free") +
  theme(strip.text.x = element_text(size = 26)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)))


## EXTENDED DATA 6F ----
wgs_mtable_cn <- read.delim(file.path(main_repo_path, "data", "ExtendedData6_sourcetable2.txt"), sep = ",")

rep_stress <- rna_sigs %>%
  select(Sample, IsablID, Replication_stress_sig, Cohort) %>%
  mutate(IsablID=na_if(IsablID, "")) %>%
  drop_na(IsablID, Replication_stress_sig)

wgs_mtable_rep <- wgs_mtable_cn %>%
  dplyr::rename("IsablID"="Sample") %>%
  left_join(rep_stress, by = "IsablID") %>%
  drop_na(Replication_stress_sig) %>%
  mutate(label=case_when(group=="HER2+" ~ "HER2+",
                         ENiClust=="ic1" ~ "ic1",
                         ENiClust=="ic2" ~ "ic2",
                         ENiClust=="ic9" ~ "ic9",
                         ENiClust=="ic6" ~ "ic6",
                         ENiClust=="ic10" ~ "IC10",
                         T ~ NA)) %>%
  drop_na(label)

her2 <- wgs_mtable_rep %>%
  filter(label=="HER2+") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~  Replication_stress_sig + ERBB2 + Cohort, data=her2, family=binomial())
her2_df <- data.frame(group="her2+",
                      sig_estimate=summary(model)$coefficient[2,1],
                      pval=summary(model)$coefficient[2,4])

ic1 <- wgs_mtable_rep %>%
  filter(label=="ic1") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ Replication_stress_sig + PPM1D + MSI2 + Cohort, data=ic1, family=binomial())
ic1_df <- data.frame(group="ic1",
                     sig_estimate=summary(model)$coefficient[2,1],
                     pval=summary(model)$coefficient[2,4])

ic6 <- wgs_mtable_rep %>%
  filter(label=="ic6") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ Replication_stress_sig + FGFR1 + Cohort, data=ic6, family = binomial())
ic6_df <- data.frame(group="ic6",
                     sig_estimate=summary(model)$coefficient[2,1],
                     pval=summary(model)$coefficient[2,4])

ic9 <- wgs_mtable_rep %>%
  filter(label=="ic9") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ Replication_stress_sig + MYC + Cohort, data=ic9, family = binomial())
ic9_df <- data.frame(group="ic9",
                     sig_estimate=summary(model)$coefficient[2,1],
                     pval=summary(model)$coefficient[2,4])

ic10 <- wgs_mtable_rep %>%
  filter(label=="IC10") %>%
  filter(Cohort != "Hartwig")
model <- glm(ecdna ~ Replication_stress_sig + ERBB2 + Cohort, data=ic10, family=binomial())
ic10_df <- data.frame(group="ic10",
                      sig_estimate=summary(model)$coefficient[2,1],
                      pval=summary(model)$coefficient[2,4])

ss_ecdna <- rbind(her2_df, ic1_df, ic6_df, ic9_df, ic10_df)

annotations <- data.frame(
  label = factor(c("HER2+", "IC1", "IC6", "IC9", "IC10"), levels = c("HER2+", "IC1", "IC6", "IC9", "IC10")),
  x = rep(-0.1, 5),
  y = rep(1, 5),
  text = c(paste0('ES=', signif(ss_ecdna$sig_estimate[1], 3), '\n', "P=", signif(ss_ecdna$pval[1], 3)), 
           paste0('ES=', signif(ss_ecdna$sig_estimate[2], 3), '\n', "P=", signif(ss_ecdna$pval[2], 3)), 
           paste0('ES=', signif(ss_ecdna$sig_estimate[3], 3), '\n', "P=", signif(ss_ecdna$pval[3], 3)), 
           paste0('ES=', signif(ss_ecdna$sig_estimate[4], 3), '\n', "P=", signif(ss_ecdna$pval[4], 3)), 
           paste0('ES=', signif(ss_ecdna$sig_estimate[5], 3), '\n', "P=", signif(ss_ecdna$pval[5], 3))
  )
)

rep_stress_plot <- wgs_mtable_rep %>%
  filter(Cohort != "Hartwig") %>%
  filter(label != "ic2") %>%
  mutate(ecdna_label = case_when(ecdna==TRUE ~ "ecDNA+",
                                 ecdna==FALSE ~ "ecDNA-"),
         label = toupper(label)) %>%
  mutate(label=factor(label, levels=c("HER2+", "IC1", "IC6", "IC9", "IC10"))) %>%
  ggboxplot(x = "ecdna_label", 
            y = "Replication_stress_sig", 
            add = "jitter", 
            color = "label",
            facet.by = "label", 
            nrow = 1) +
  labs(y="Replication Stress") +
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
rep_stress_plot <- rep_stress_plot + geom_text(data = annotations, aes(x = x, y = y, label = text),
                                               color = "black", size = 5, hjust = 0, vjust = 0)

## Meta-analysis
dtf <- data.frame(coef = ss_ecdna$sig_estimate, se = ss_ecdna$se)
meta_data <- escalc(measure="OR", yi = coef, sei = se, 
                    data = dtf)
res <- rma(yi, vi, data=meta_data)

ss_ecdna_her_er <- ss_ecdna %>%
  filter(group %in% c('her2+', 'ic1', 'ic6'))
dtf_her_er <- data.frame(coef = ss_ecdna_her_er$sig_estimate, se = ss_ecdna_her_er$se)
meta_data_her_er <- escalc(measure="OR", yi = coef, sei = se, 
                           data = dtf_her_er)
res_her_er <- rma(yi, vi, data=meta_data_her_er)

# SAVE ----
pdf(file.path(main_repo_path, "plots", "ExtendedData6a.pdf"), height=6, width=11)
gt = ggplot_gtable(ggplot_build(plt1a))
gt$layout$l[grep('panel-1-1', gt$layout$name)]
gt$widths[5] = 2.5*gt$widths[5]
grid.draw(gt)
dev.off()

rbind(
  tcga %>% 
    filter(plotting_group == "IC subgroup") %>% 
    compare_means(sig1 ~ label, data=.) %>%
    adjust_pvalue(method = "fdr") %>%
    filter(group1 == "ER+ High" & group2 == "ER+ Typical" |
             group1 == "ER+ Typical" & group2 == "IC10" |
             group1 == "ER+ Typical" & group2 == "IC4ER-" |
             group1 == "ER+ Typical" & group2 == "HER2+/ER+" |
             group1 == "ER+ Typical" & group2 == "HER2+/ER-" |
             group1 == "IC10" & group2 == "IC4ER-") %>%
    mutate(p.adj.format=paste0("FDR = ",signif(p.adj, 3))) %>%
    select(group1, group2, p.adj, p.adj.format),
  tcga %>% 
    filter(plotting_group == "ER+ Typical") %>% 
    compare_means(sig1 ~ label, data=.) %>%
    adjust_pvalue(method = "fdr") %>%
    mutate(p.adj.format=paste0("FDR = ",signif(p.adj, 3))) %>%
    select(group1, group2, p.adj, p.adj.format)) %>%
  write.table(., file.path(main_repo_path, "plots", 'ExtendedData6a_stats.txt'),
              quote = F, row.names = F, sep = "\t")

pdf(file.path(main_repo_path, "plots", "ExtendedData6c.pdf"), height=6, width=11)
plt1b
dev.off()

ss_tcga_ic %>%
  write.table(., file.path(main_repo_path, "plots", 'ExtendedData6c_stats.txt'),
              quote = F, row.names = F, sep = "\t")

pdf(file.path(main_repo_path, "plots", "ExtendedData6b.pdf"), height=6, width=11)
gt2 = ggplot_gtable(ggplot_build(plt2a))
gt2$layout$l[grep('panel-1-1', gt2$layout$name)]
gt2$widths[5] = 2.5*gt2$widths[5]
grid.draw(gt2)
dev.off()

ss_mb_eniclust %>%
  write.table(., file.path(main_repo_path, "plots", 'ExtendedData6b_stats.txt'),
              quote = F, row.names = F, sep = "\t")

pdf(file.path(main_repo_path, "plots", "ExtendedData6d.pdf"), height=6, width=12)
plt2b
dev.off()

ss_mb_ic %>%
  write.table(., file.path(main_repo_path, "plots", 'ExtendedData6d_stats.txt'),
              quote = F, row.names = F, sep = "\t")

pdf(file.path(main_repo_path, "plots", "ExtendedData6g.pdf"), height=6, width=9)
gt = ggplot_gtable(ggplot_build(plt1a_cgas))
gt$layout$l[grep('panel-1-1', gt$layout$name)]
gt$widths[5] = 2.5*gt$widths[5]
grid.draw(gt)
dev.off()


rbind(
  tcga_cgas %>% 
    filter(plotting_group == "IC subgroup") %>%
    compare_means(cGAS_STING_sig ~ label, data=.) %>%
    adjust_pvalue(method = "fdr") %>%
    mutate(p.adj.format=paste0("FDR = ",signif(p.adj, 3))) %>%
    select(group1, group2, p.adj, p.adj.format),
  tcga_cgas %>% 
    filter(plotting_group == "ER+ Typical") %>%
    compare_means(cGAS_STING_sig ~ label, data=.) %>%
    adjust_pvalue(method = "fdr") %>%
    mutate(p.adj.format=paste0("FDR = ",signif(p.adj, 3))) %>%
    select(group1, group2, p.adj, p.adj.format)) %>%
  write.table(., file.path(main_repo_path, "plots", 'ExtendedData6g_stats.txt'),
              quote = F, row.names = F, sep = "\t")

pdf(file.path(main_repo_path, "data", "ExtendedData6h.pdf"), height=6, width=9)
gt2 = ggplot_gtable(ggplot_build(plt2a_cgas))
gt2$layout$l[grep('panel-1-1', gt2$layout$name)]
gt2$widths[5] = 2.5*gt2$widths[5]
grid.draw(gt2)
dev.off()

rbind(
  mb_cgas %>% 
    filter(plotting_group == "IC subgroup") %>%
    compare_means(cGAS_STING_sig ~ label, data=.) %>%
    adjust_pvalue(method = "fdr") %>%
    mutate(p.adj.format=paste0("FDR = ",signif(p.adj, 3))) %>%
    select(group1, group2, p.adj, p.adj.format),
  mb_cgas %>% 
    filter(plotting_group == "ER+ Typical") %>%
    compare_means(cGAS_STING_sig ~ label, data=.) %>%
    adjust_pvalue(method = "fdr") %>%
    mutate(p.adj.format=paste0("FDR = ",signif(p.adj, 3))) %>%
    select(group1, group2, p.adj, p.adj.format))%>%
  write.table(., file.path(main_repo_path, "plots", 'ExtendedData6h_stats.txt'),
              quote = F, row.names = F, sep = "\t")

pdf(file.path(main_repo_path, "data", "ExtendedData6f.pdf"), height=6, width=9)
rep_stress_plot
dev.off()

ss_ecdna %>%
  write.table(., file.path(main_repo_path, "plots", 'ExtendedData6f_stats.txt'),
              quote = F, row.names = F, sep = "\t")
