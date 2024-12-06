### CREATE EXTENDED DATA FIGURE 8D #############################################################################
# Barplot of TME in ER+ Typical IDC vs ILC (Metabric)

### PREAMBLE #####################################################################################
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
tme_colors <- c("IE/F"="#5CBDB9", 
                "IE"="#5036C4", 
                "F"="#AA1946", 
                "D"="#B48913")

plot_tme_prop <- function(df) {
  df %>%
    ggplot(aes(x=histo_type, y=n, fill=TME_subtype)) +
    geom_bar(stat="identity", position="fill") +
    geom_text(aes(y=c(1), label=label), vjust=-.2, size=5)+
    facet_wrap(~group, labeller = labeller(name = c(group = "ER+ Typical"))) +
    scale_fill_manual(values = tme_colors) +
    theme_LM +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 14),
          strip.text.x = element_text(size = 18),
          strip.text.y = element_text(size = 14),
          legend.position = "none") +
    labs(y="Proportion of samples by TME subtype", x="", title="Primary (METABRIC)") +
    coord_cartesian(ylim=c(0,1.01))
}

### MAIN ##########################################################################################
tme <- read.delim(file.path(main_repo_path, "data", "ExtendedData8d_sourcetable.txt"), sep = ",")

tme_count <- tme %>%
  group_by(histo_type, TME_subtype) %>%
  summarise(n=n()) %>%
  mutate(label=ifelse(TME_subtype=="D", paste0(sum(n)),""))

df <- data.frame("ILC" = c(tme %>%
                             count(histo_type, TME_subtype) %>%
                             filter(histo_type == "ILC" & TME_subtype == "IE/F") %>%
                             pull(n), 
                           tme %>%
                             count(histo_type, TME_subtype) %>%
                             filter(histo_type == "ILC" & TME_subtype != "IE/F") %>%
                             summarise(sum=sum(n)) %>%
                             pull(sum)), 
                 "IDC" = c(tme %>%
                             count(histo_type, TME_subtype) %>%
                             filter(histo_type == "IDC" & TME_subtype == "IE/F") %>%
                             pull(n), 
                           tme %>%
                             count(histo_type, TME_subtype) %>%
                             filter(histo_type == "IDC" & TME_subtype != "IE/F") %>%
                             summarise(sum=sum(n)) %>%
                             pull(sum)), 
                 row.names = c("TRUE", "FALSE"))

ss_tme_idc_ilc <- data.frame(
  p=fisher.test(df)$p,
  or=fisher.test(df)$estimate[[1]]
)

# check immune signature
tme %>%
  mutate(group=factor(group, levels=c("ER+ Typical", "ER+ High", "HER2+", "IC10", "IC4ER-"))) %>%
  drop_na(group, histo_type) %>%
  pivot_longer(cols = c('Immune', 'Fibrotic')) %>%
  group_by(group, name) %>%
  wilcox_test(value ~ histo_type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

# SAVE ----
pdf(file.path(main_repo_path, "plots", "ExtendedData8d.pdf"), height=5, width=3)
plot_tme_prop(tme_count %>% mutate(group = "ER+ Typical"))
dev.off()

ss_tme_idc_ilc %>%
  write.table(file.path(main_repo_path, "data", "ExtendedData8d_stats.txt"), 
              quote = F, row.names = F, sep = "\t")
