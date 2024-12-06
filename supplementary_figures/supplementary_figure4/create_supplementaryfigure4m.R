# CREATE SUPPLEMENTARY FIGURE 4m ----
# Barplot of TME in ER+ Typical IDC vs ILC (Metabric)

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
group_colours <- c("ER+ High" = "#fa954eff",
                   "ER+ Typical" = "#a1c1d4ff",
                   "HER2+" = "#8b0100ff",
                   "IC10" = "#7c26ccff",
                   "IC4ER-" = "#c2b7dbff")

plot_barplot_histo_ic <- function(df) {
  df %>%
    ggplot(aes(x=histo_type, y=n, fill=group)) +
    geom_bar(stat="identity", position="fill") +
    geom_text(aes(y=c(1), label=label), vjust=-.2, size=4)+
    labs(x="", y="Proportion of samples by subgroups", title = "METABRIC") +
    scale_fill_manual(name = "", values=group_colours) +
    theme_LM + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
                     axis.text.y = element_text(size = 12),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     legend.text = element_text(size = 14),
                     title = element_text(size = 16))
}

# MAIN ----
metabric <- fread(file.path(main_repo_path, "data", "ExtendedData3l_sourcetable.txt")) %>%
  filter(Cohort=="METABRIC") %>%
  distinct(Sample, .keep_all = TRUE) %>%
  select(Sample, Cohort, group, histo_type)

metabric_df <- metabric %>%
  filter(histo_type %in% c("IDC", "ILC")) %>%
  group_by(histo_type, group) %>%
  summarise(n=n()) %>%
  mutate(label=ifelse(group=="IC4ER-", paste0("n=", sum(n)),""))

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure4m.pdf"), height=4, width=4)  
plot_barplot_histo_ic(metabric_df)
dev.off()

dat <- data.frame(
  "ILC" = c(metabric_df %>%
              filter(histo_type == 'ILC' & group == "ER+ Typical") %>%
              pull(n), 
            metabric_df %>%
              filter(histo_type == 'ILC' & group == "ER+ High") %>%
              pull(n)),
  "IDC" = c(metabric_df %>%
              filter(histo_type == 'IDC' & group == "ER+ Typical") %>%
              pull(n), 
            metabric_df %>%
              filter(histo_type == 'IDC' & group == "ER+ High") %>%
              pull(n)),
  row.names = c("Low", "High"),
  stringsAsFactors = FALSE
)
ilc_enrichment_ss <- data.frame(
  p=fisher.test(dat)$p,
  or=fisher.test(dat)$estimate[[1]]
)

fwrite(ilc_enrichment_ss, file.path(main_repo_path, "plots", "SupplementaryFigure4m_stats.tsv"))

