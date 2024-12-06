### CREATE EXTENDED DATA 1E #############################################################################
# Primary to metastatic subtype stability 
### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggsankey)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

# Read data
data <- read.delim(file.path(main_repo_path, 'data', 'paired_sample_metadata.tsv'))

# IC subtype ------
mets <- data %>% filter(tumor_site == "metastasis") %>% select(Individual.ID, Sample.ID, 
                                                               Metastatic = ENiClust)
primaries <- data %>% filter(tumor_site == "primary") %>% select(Individual.ID, 
                                                                 Primary = ENiClust)

order_eniclust <- c(paste0("IC", c(1, 2, 6, 9, "3/IC7", "4ER+", 8, 4, 5, 10, "4ER-" )))

plot_data <- primaries %>% right_join(mets) %>% 
  make_long(Primary, Metastatic) %>% mutate(
    node = factor(node, levels = rev(order_eniclust)), 
    next_node = factor(next_node, levels = rev(order_eniclust))) 

pdf(file.path(main_repo_path, "plots", "ExtendedData1E_top.pdf"), height = 6, width = 6)
ggplot(plot_data, aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = node, 
                      label = node)) +
  scale_fill_manual(values = ic_colors, guide="legend") +
  theme_sankey(base_size = 16) +
  geom_sankey(flow.alpha=0.75) +
  geom_sankey_label() +
  labs(x="") + 
  rremove("legend")
dev.off()

# Subgroup ------
mets_subgroup <- data %>% filter(tumor_site == "metastasis", !is.na(Subgroup)) %>% 
  select(Individual.ID, Sample.ID, Metastatic = Subgroup)
primaries_subgroup <- data %>% filter(tumor_site == "primary", !is.na(Subgroup)) %>% 
  select(Individual.ID, Primary = Subgroup)

plot_data <- primaries_subgroup %>% right_join(mets_subgroup) %>% 
  filter(!is.na(Primary)) %>%
  make_long(Primary, Metastatic) %>% mutate(
    node = factor(node, levels = rev(names(group_colors))), 
    next_node = factor(next_node, levels = rev(names(group_colors)))) 

pdf(file.path(main_repo_path, "plots", "ExtendedData1E_bottom.pdf"), height = 6, width = 6)
ggplot(plot_data, aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = node, 
                      label = node)) +
  scale_fill_manual(values = group_colors, guide="legend") +
  theme_sankey(base_size = 16) +
  geom_sankey(flow.alpha=0.75) +
  geom_sankey_label() +
  labs(x="") + 
  rremove("legend")
dev.off()

### SAVE ##########################################################################################
primaries %>% right_join(mets) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData1e_sourcetable_top.txt"), 
              quote = F, row.names = F, sep = "\t")

primaries_subgroup %>% right_join(mets_subgroup) %>% 
  filter(!is.na(Primary)) %>%
  write.table(file.path(main_repo_path, "data", "ExtendedData1e_sourcetable_bottom.txt"), 
              quote = F, row.names = F, sep = "\t")

