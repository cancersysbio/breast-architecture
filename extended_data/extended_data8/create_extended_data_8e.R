### CREATE EXTENDED DATA 8E #############################################################################
# Bootstaping IMC data for TNBC 

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

imc_data <- file.path(main_repo_path, 'data','single_cell_data.csv')
if(!file.exists(imc_data)) {
  stop("Error: METABRIC IMC data path not found, download here: ")
} 

# Read data
imc <- read.delim(imc_data, sep = ",")
clinical_metabric <- read.delim(file.path(main_repo_path, "data", "rna_megatable.txt")) %>% filter(Cohort == "METABRIC") %>% select(Sample, group)

imc_prop <- imc %>% select(Sample = metabricId, description) %>% mutate(`Cell type` = ifelse(description %in% c("Fibroblasts", "Fibroblasts CD68+", "Vascular SMA+", "Myofibroblasts"), "Fibroblasts", 
                                                                                             ifelse(grepl("Macrophages", description), "Macrophages", 
                                                                                                    ifelse(grepl("HR|Basal|HER2+", description), "Cancer", description)))) %>% 
  group_by(Sample, `Cell type`) %>% 
  summarize(n = n()) %>% group_by(Sample) %>% mutate(total = sum(n), proportion = n/total) 

cell_types <- unique(imc_prop$`Cell type`)

imc_prop <- imc_prop %>% left_join(clinical_metabric) %>% filter(!is.na(group)) %>% data.frame(check.names = F)

samples_by_group <- imc_prop %>% select(Sample, group) %>% unique() %>% pull(group) %>% table() %>% data.frame() %>% 
  dplyr::rename(group = 1, n = 2)

bootstrap <- lapply(sort(samples_by_group$group), function(x) {
  imc_subset <- imc_prop %>% filter(group == x)
  
  boot <- sapply(cell_types, function(subtype) {
    imc_subset_cell_type <- imc_subset %>% filter(`Cell type` == subtype)
    n <- nrow(imc_subset_cell_type)
    meds <- sapply(1:1000, function(a) {
      s <- sample(1:n, n, replace = T)
      return(subtype = mean(imc_subset_cell_type[s, "proportion"]))
    })
  })
})

data <- lapply(1:5, function(i) {
  bootstrap[[i]] %>% data.frame() %>% mutate(index = 1:1000, group = sort(samples_by_group$group)[i]) %>% pivot_longer(values_to = "Mean proportion", names_to = "Cell type", cols = 1:8)
}) %>% do.call(rbind, .)

source_data <- data %>% filter(group %in% c("IC10", "IC4ER-"), 
                               `Cell type` %in% c("Fibroblasts", "T.cells")) %>% 
  mutate(`Cell type` = gsub("T.cells", "T cells", `Cell type`, fixed = T))

## Save plot ----
pdf(file.path(main_repo_path, "plots", "ExtendedData8E.pdf"), height = 5, width = 5)
ggboxplot(source_data, "group", "Mean proportion", color = "group", 
          palette = group_colors, facet.by = "Cell type", scales = "free", 
          xlab = "", main = "Mean proportions from bootstrapping\nIMC data (n = 49)") +
  theme_LM +
  rremove("legend") + 
  rotate_x_text(35) +
  stat_compare_means()
dev.off()

### SAVE ##########################################################################################
source_data %>% 
  select(Bootstrap_Index = index, Subtype = group, `Cell type`, `Mean proportion`) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData8e_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")