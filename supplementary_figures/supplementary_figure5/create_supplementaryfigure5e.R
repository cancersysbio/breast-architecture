# CREATE SUPPLEMENTARY FIGURE 5E ----
# Barplot of AA ecDNA calls

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
theme_LM = theme_classic() + grids()
cat_colour = c(
  'Both' = '#6AE484',
  'None' = '#EF2D56',
  'Full' = '#6A8EAE',
  'Downsampled' = '#BBAB8B'
)
ic_colours <- c("11q13" = "#00ee77ff",
                "17q23" = "#ff5500ff",
                "8q24" = "#ee82edff")

plot_one_barplot <- function(df) {
  df %>%
    mutate(category = factor(category, levels = c('Both', 'Full', 'Downsampled', 'None'))) %>%
    ggplot(., aes(x = "", y = count, fill = category)) +
    geom_bar(stat = "identity", position = 'fill') +
    geom_text(aes(label = count), size = 4, position = position_fill(vjust = 0.5)) +
    labs(title = "",
         x = NULL, y = "Proportion of Samples") +
    scale_fill_manual(name = "", values=cat_colour) +
    theme_LM +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_text(size = 14), 
          strip.text.y = element_text(size = 14))
}

# MAIN ----
ecdna <- fread(file.path(main_repo_path, "data", "SupplementaryFigure5e_sourcetable.txt"))

count_categories <- ecdna[, .(
  Both = uniqueN(sample_name[`ecDNA+.x` == "Positive" & `ecDNA+.y` == "Positive"]),
  Full = uniqueN(sample_name[`ecDNA+.x` == "Positive" & `ecDNA+.y` != "Positive"]),
  Downsampled = uniqueN(sample_name[`ecDNA+.x` != "Positive" & `ecDNA+.y` == "Positive"]),
  None = uniqueN(sample_name[`ecDNA+.x` != "Positive" & `ecDNA+.y` != "Positive"])
)]

plot_data <- melt(count_categories, variable.name = "category", value.name = "count")

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure5e.pdf"), height=4, width=3)
plot_one_barplot(plot_data)
dev.off()