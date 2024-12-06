# CREATE SUPPLEMENTARY FIGURE 6F ----
# Boxplot of downsampling to met

# PREAMBLE ----
library(rstudioapi)
library(data.table)
library(magrittr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(yaml)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# FUNCTIONS ----
theme_LM = theme_classic() + grids()
colours <- c("primary"="#9ebff1ff","met"="#fe85cfff")
boot_prop <- function(data, n, rep) {
  resamples <- lapply(1:rep, function(i) sample(t(data), replace=T, size=n))
  prop_all <- lapply(resamples, function(x) prop.table(table(x)))
  prop_d <- sapply(prop_all, function(x) x[1] + x[2]) #D and F are 1 and 2
  # std.err <- sqrt(var(prop_d))
  boot_prop <- list(resamples=resamples, prop_d=prop_d) 
  return(boot_prop)
}

plot_boxplot_icsubgroup <- function(df) {
  df %>%
    ggplot(aes(x=group, y=prop)) +
    geom_jitter(aes(colour=Stage, alpha=Stage)) +
    geom_boxplot(alpha=0) +
    geom_text(data=n_all_df, 
              aes(x=rownames(n_all_df), y=rep(1.1, 5), label=paste0("n=",n_all)), size=5) +
    theme_pubr(base_size = 12) +
    labs(x="", y="Proportion of immune \n depleted samples (D or F)", title = "", color="") +
    scale_color_manual(values = colours, labels=c("Primary", "Metastic")) +
    scale_alpha_manual(values=c(0.1,1), guide = 'none') +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="bottom") +
    scale_y_continuous(breaks=seq(0, 1.2, 0.2)) +
    theme_LM + 
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.position = "top",
          title = element_text(size = 16))
}

# MAIN ---
tme_calls <- fread(file.path(main_repo_path, "data", "SupplementaryFigure6f_sourcetable.txt")) %>%
  drop_na(Cohort)

# Check liver only ----
subtypes <- c("ER+ High", "ER+ Typical", "HER2+", "IC10", "IC4ER-")
n_liver <- tme_calls %>%
  filter(group %in% subtypes & Stage == "Metastatic" & biopsySite == "Liver") %>%
  count(group) %>%
  pull(n) %>%
  c()
names(n_liver) <- subtypes

set.seed(24601)
boot_liver <- list()
for (i in 1:length(subtypes)) {
  subtype <- subtypes[i]
  n_boot <- n_liver[i]
  
  primary <- tme_calls %>% 
    filter(Stage == "Primary" & group == subtype) %>%
    as.data.frame() %>%
    select(TME_subtype)
  boot_liver[[i]] <- boot_prop(primary, n = n_boot, rep = 1000)
}

boot_liver_sd <- list()
for (i in 1:length(subtypes)) {
  boot_liver_sd[[i]] <- boot_liver[[i]]$prop_d
}

boot_liver_sd <- data.frame(sapply(boot_liver_sd, function(x) x[1:max(lengths(boot_liver_sd))]))
names(boot_liver_sd) <- subtypes
boot_liver_sd_pri <- boot_liver_sd %>%
  pivot_longer(names_to = "group", cols = `ER+ High`:`IC4ER-`, values_to = "prop") %>%
  mutate(Cohort="primary")

# All ----
subtypes <- c("ER+ High", "ER+ Typical", "HER2+", "IC10", "IC4ER-")
n_all <- tme_calls %>%
  filter(group %in% subtypes & Stage == "Metastatic" & Cohort == "Hartwig") %>%
  drop_na(TME_subtype) %>%
  count(group) %>%
  pull(n) %>%
  c()
names(n_all) <- subtypes

set.seed(24601)
boot_all <- list()
for (i in 1:length(subtypes)) {
  subtype <- subtypes[i]
  n_boot <- n_all[i]
  
  primary <- tme_calls %>% 
    filter(Stage == "Primary" & group == subtype & Cohort == "TCGA") %>%
    drop_na(TME_subtype) %>%
    as.data.frame() %>%
    select(TME_subtype)
  boot_all[[i]] <- boot_prop(primary, n = n_boot, rep = 1000)
}

boot_all_sd <- list()
for (i in 1:length(subtypes)) {
  boot_all_sd[[i]] <- boot_all[[i]]$prop_d
}

boot_all_sd <- data.frame(sapply(boot_all_sd, function(x) x[1:max(lengths(boot_all_sd))]))
names(boot_all_sd) <- subtypes

liver_mets <- tme_calls %>%
  drop_na(TME_subtype) %>%
  filter(group %in% subtypes & Stage == "Metastatic" & Cohort == "Hartwig" & biopsySite == 'Liver') %>%
  group_by(group, TME_subtype) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(prop=n/sum(n)) %>%
  filter(TME_subtype %in% c("D", "F")) %>%
  summarise(prop=sum(prop)) %>%
  mutate(Cohort="met")

boot_liver_all <- boot_liver_sd %>%
  pivot_longer(names_to = "group", cols = `ER+ High`:`IC4ER-`, values_to = "prop") %>%
  mutate(Cohort="primary") %>%
  bind_rows(liver_mets)

all_mets <- tme_calls %>%
  drop_na(TME_subtype) %>%
  filter(group %in% subtypes & Stage == "Metastatic" & Cohort == "Hartwig") %>%
  group_by(group, TME_subtype) %>%
  summarise(n=n()) %>%
  group_by(group) %>%
  mutate(prop=n/sum(n)) %>%
  filter(TME_subtype %in% c("D", "F")) %>%
  summarise(prop=sum(prop)) %>%
  mutate(Cohort="met")

boot_all_all <- boot_all_sd %>%
  pivot_longer(names_to = "group", cols = `ER+ High`:`IC4ER-`, values_to = "prop") %>%
  mutate(Cohort="primary") %>%
  bind_rows(all_mets)

boot_all_data <- boot_all_all %>%
  drop_na()

n_all_df <- data.frame(n_all)

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure6f.pdf"), width=4, height=6)
boot_all_data %>%
  mutate(Stage=factor(Cohort, levels=c('primary', 'met'))) %>%
  plot_boxplot_icsubgroup(.)
dev.off()