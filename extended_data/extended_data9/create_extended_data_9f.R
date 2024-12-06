### CREATE EXTENDED DATA 9F #############################################################################
# Alterations in immune pathways in metastatic vs primary tumors 

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(patchwork)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

# Read data -- generated in immune_pathway_alterations.R
immune_alterations_data_path <- file.path(main_repo_path, 'data','immune_pathway_alterations.txt')
if(!file.exists(immune_alterations_data_path)) {
  stop("Error: Immune alteration data not found, run immune_pathway_alterations.R analysis script. ")
} 

data <- read.delim(immune_alterations_data_path)

primary_megatable <- read.delim(file.path(main_repo_path, 'data','primary_megatable.txt'))
metastatic_megatable <- read.delim(file.path(main_repo_path, 'data','metastatic_megatable.txt'))

samples <- rbind(primary_megatable, metastatic_megatable)

data <- data %>% inner_join(samples) %>% filter(!is.na(group))

# Extended Data 9F
total_groups <- rbind(
  data.frame(table(primary_megatable$group)) %>% rename(group = Var1, total = Freq) %>% mutate(Stage = "Primary"),
  data.frame(table(metastatic_megatable$group)) %>% rename(group = Var1, total = Freq) %>% mutate(Stage = "Metastatic")
)

alterations_by_pathway <- data %>% select(Pathway, Sample, group, Stage) %>% unique() %>% 
  group_by(group, Pathway, Stage) %>% summarize(n = n()) %>% left_join(total_groups) %>% 
  mutate(proportion = n/total, 
         Stage = factor(Stage, levels = c("Primary", "Metastatic"))) 

data_tests <- alterations_by_pathway %>% mutate(Yes = n, No = total - n) %>% 
  group_by(group, Pathway) %>% filter(n() == 2) 

test_results <- data_tests %>% select(group, Pathway) %>% unique()
tests <- data_tests %>% select(Yes, No) %>% group_map(~fisher.test(.x))

test_results$p <- sapply(tests, function(x) x$p.value)
test_results$or <- sapply(tests, function(x) x$estimate) 

test_results <- test_results %>% left_join(total_groups %>% group_by(group) %>% summarize(Samples = sum(total)))

# For any pathway 
any_pathway_tests <- data %>% select(Sample, group, Stage) %>% unique() %>% 
  group_by(group, Stage) %>% summarize(Yes = n()) %>% left_join(total_groups) %>% 
  mutate(No = total - Yes) %>% select(Yes, No) %>%  group_map(~fisher.test(.x))

any_pathway_results <- data.frame(
  group = sort(unique(data$group)), 
  Pathway = "Any",
  p = sapply(any_pathway_tests, function(x) x$p.value), 
  or = sapply(any_pathway_tests, function(x) x$estimate)
) %>% left_join(total_groups %>% group_by(group) %>% summarize(Samples = sum(total)))


# For all samples by pathway
total_samples <- data.frame(Stage = c("Primary", "Metastatic"), total = c(nrow(primary_megatable), nrow(metastatic_megatable)))

all_samples_tests <- data %>% select(Sample, Pathway, Stage) %>% unique() %>% 
  group_by(Pathway, Stage) %>% summarize(Yes = n()) %>% left_join(total_samples) %>% 
  mutate(No = total - Yes) %>% select(Yes, No) %>%  group_map(~fisher.test(.x))

all_samples_results <- data.frame(
  Pathway = sort(unique(data$Pathway)), 
  group = "All Samples",
  p = sapply(all_samples_tests, function(x) x$p.value), 
  or = sapply(all_samples_tests, function(x) x$estimate),
  Samples = sum(total_groups$total)
) 

# For all samples and any pathway 
all_samples_any <- data %>% select(Sample, Stage) %>% unique() %>% 
  group_by(Stage) %>% summarise(Yes = n()) %>% left_join(total_samples) %>% 
  mutate(No = total - Yes) %>% select(Yes, No) %>% fisher.test()

all_samples_any_results <- data.frame(
  Pathway = "Any", 
  group = "All Samples", 
  p = all_samples_any$p.value, 
  or = all_samples_any$estimate, 
  Samples = sum(total_groups$total)
)


data_ext9f <- rbind(test_results, any_pathway_results, all_samples_results, all_samples_any_results)
data_ext9f$padj <- p.adjust(data_ext9f$p, method = "fdr")
data_ext9f$Pathway <- factor(data_ext9f$Pathway, levels = sort(unique(data_ext9f$Pathway))[c(2, 1, 3:8)])
limit <- max(abs(log(data_ext9f$or)))

## Save plot ----
pdf(file.path(main_repo_path, "plots", "ExtendedData9F.pdf"), height = 7, width = 6)
ggplot(data_ext9f, aes(x = group, y = Pathway, color = log(or))) +
  geom_tile(aes(fill = -log(padj)), alpha = 0.5) +
  geom_point(aes(size = Samples)) +
  theme_LM +
  scale_color_distiller(palette = "RdBu", limit = c(-limit, limit)) +
  scale_fill_gradient(low = "white", high = "black") +
  scale_size_continuous(range = c(6, 12), breaks = c(50, 100, 200, 500, 1000)) +
  scale_y_discrete(limits = rev(levels(data_ext9f$Pathway))) +
  labs(x = "") +
  rotate_x_text(35)
dev.off()


### SAVE ##########################################################################################
data_ext9f %>% 
  select(Subtype = group, Pathway, p, OR = or, Samples, padj) %>% 
  write.table(file.path(main_repo_path, "data", "Extended_Data9f_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
