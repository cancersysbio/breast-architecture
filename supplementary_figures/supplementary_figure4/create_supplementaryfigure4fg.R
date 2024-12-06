### CREATE SUPPLEMENTARY FIGURE 4f-g #############################################################################
# creates supplementary figure 4f and 4g
# supplementary figure 4f provides the concordance of archetypes defined by Pareto and NMF in primary, metastatic and primary+metastatic samples combined.
# supplementary figure 4g provides the association between the three genomic archetypes and the breast cancer subgroups with archetypes defined by Pareto and NMF.

### PREAMBLE #####################################################################################
library(yaml)
library(dplyr)
library(insight, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
library(modelsummary, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
library(performance, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
library(ggplot2)
library(ggpubr)
require(patchwork)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

col.group = c("ER+ High" = "#fa954eff","ER+ Typical" = "#a1c1d4ff","HER2+" = "#8b0100ff",
              "IC10" = "#7c26ccff", "IC4ER-" = "#c2b7dbff")

### MAIN ##########################################################################################
basedir <- "/oak/stanford/groups/ccurtis2/users/"
outdir = paste0(basedir,"lisem/bc_landscape/revision/figures/")

## Load megatable----
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")
mega_1$Sample_Type = "Primary"
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")
mega_2$Sample_Type = "Metastatic"
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])

mega_1[which(is.na(mega_1$group)),] # 2 IC4 samples with missing hormone receptor status
mega_2[which(is.na(mega_2$group)),] # 7 IC4 samples with missing hormone receptor status

## Load Pareto----
tern_arc_prim <- read.table(paste0(basedir,"lisem/bc_landscape/submission/data/Figure2btop_sourcetable.txt"), header = T, sep="\t")
tern_arc_met <- read.table(paste0(basedir,"lisem/bc_landscape/revision/data/Figure2bbottom_sourcetable.txt"), header = T, sep="\t")
tern_arc_all <- read.table(paste0(basedir,"lisem/bc_landscape/revision/data/ExtendedData3g_sourcetable.txt"), header = T, sep="\t")

## Load NMF----
nmf_rank_prim <- read.table(paste0(basedir,"lisem/bc_landscape/revision/data/nmf_rank_prim.txt"), header = T, sep="")
nmf_rank_met <- read.table(paste0(basedir,"lisem/bc_landscape/revision/data/nmf_rank_met.txt"), header = T, sep="")
nmf_rank_all <- read.table(paste0(basedir,"lisem/bc_landscape/revision/data/nmf_rank_all.txt"), header = T, sep="")

panelA <- as.data.frame(cor(nmf_rank_prim, tern_arc_prim[,paste0("Arc",1:3)])) %>%
  tibble::rownames_to_column(var = "LF") %>%
  tidyr::pivot_longer(!LF, names_to = "Arc", values_to = "Correlation") %>%
  dplyr::mutate(dataset = "Primary") %>%
  bind_rows(as.data.frame(cor(nmf_rank_met, tern_arc_met[,paste0("Arc",1:3)])) %>%
              tibble::rownames_to_column(var = "LF") %>%
              tidyr::pivot_longer(!LF, names_to = "Arc", values_to = "Correlation") %>%
              dplyr::mutate(dataset = "Metastatic")) %>%
  bind_rows(as.data.frame(cor(nmf_rank_all, tern_arc_all[,paste0("Arc",1:3)])) %>%
              tibble::rownames_to_column(var = "LF") %>%
              tidyr::pivot_longer(!LF, names_to = "Arc", values_to = "Correlation") %>%
              dplyr::mutate(dataset = "Primary+Metastatic")) %>%
  ggplot(aes(x=LF,y=Arc,fill=Correlation)) +
  geom_tile() + geom_text(aes(label = round(Correlation,2))) +
  facet_grid(. ~ factor(dataset, levels = c("Primary","Metastatic","Primary+Metastatic"))) +
  scale_fill_gradient2(name = "Correlation", low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-1,1)) +
  xlab("") + ylab("") +
  theme_LM + ggtitle("Pareto archetypes vs. NMF latent factors") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        title = element_text(size=16))

## Stratification IC subgroups----
data <- tern_arc_prim %>% bind_cols(nmf_rank_prim)
results_prim <- lapply(c("Arc1", "Arc2", "Arc3", "LF1", "LF2", "LF3"), function(x) lm(data[,x] ~ data$Subgroup)) |>
  modelplot(draw = FALSE) |>
  transform("Function" = "lm()") %>% 
  dplyr::mutate(term = ifelse(grepl("data",term),gsub("data[$]Subgroup","",term),"ER+ High"), 
                model = c("Arc1", "Arc2", "Arc3", "LF1", "LF2", "LF3")[as.numeric(gsub("[)]","",gsub("[(]","",model)))],
                dataset = "Primary")

data <- tern_arc_met %>% bind_cols(nmf_rank_met)
results_met <- lapply(c("Arc1", "Arc2", "Arc3", "LF1", "LF2", "LF3"), function(x) lm(data[,x] ~ data$Subgroup)) |>
  modelplot(draw = FALSE) |>
  transform("Function" = "lm()") %>% 
  dplyr::mutate(term = ifelse(grepl("data",term),gsub("data[$]Subgroup","",term),"ER+ High"), 
                model = c("Arc1", "Arc2", "Arc3", "LF1", "LF2", "LF3")[as.numeric(gsub("[)]","",gsub("[(]","",model)))],
                dataset = "Metastatic")

data <- tern_arc_all %>% bind_cols(nmf_rank_all)
results_all <- lapply(c("Arc1", "Arc2", "Arc3", "LF1", "LF2", "LF3"), function(x) lm(data[,x] ~ data$Subgroup)) |>
  modelplot(draw = FALSE) |>
  transform("Function" = "lm()") %>% 
  dplyr::mutate(term = ifelse(grepl("data",term),gsub("data[$]Subgroup","",term),"ER+ High"), 
                model = c("Arc1", "Arc2", "Arc3", "LF1", "LF2", "LF3")[as.numeric(gsub("[)]","",gsub("[(]","",model)))],
                dataset = "Primary+Metastatic")

results <- rbind(rbind(results_prim, results_met), results_all)
head(results)

panelB <- results %>%
  dplyr::rowwise() %>% dplyr::mutate(axis = ifelse(grepl("Arc",model),gsub("Arc","",model),gsub("LF","",model))) %>%
  ggplot(aes(y = factor(term, levels = rev(c("IC10","IC4ER-","ER+ Typical","ER+ High","HER2+"))), x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(color = ifelse(grepl("Arc", model),"Pareto","NMF")), position = position_dodge(width = .5)) + 
  facet_grid(axis ~ factor(dataset, levels = c("Primary","Metastatic","Primary+Metastatic"))) + 
  scale_color_manual(name = "", values = c("#ffa600","#6d0047")) +
  ylab("") + xlab("Estimate") +
  theme_LM + ggtitle("Association with IC subgroups (linear regression)") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        title = element_text(size=16))

### SAVE ##########################################################################################
## Figure----
pdf(paste0(outdir, 'SupplementaryFigure4f.pdf'), width = 8.04, height = 7.96)
panelA
dev.off()

pdf(paste0(outdir, 'SupplementaryFigure4g.pdf'), width = 8.04, height = 7.96)
panelB
dev.off()

## SourceData---- 
sourcedata <- as.data.frame(cor(nmf_rank_prim, tern_arc_prim[,paste0("Arc",1:3)])) %>%
  tibble::rownames_to_column(var = "LF") %>%
  tidyr::pivot_longer(!LF, names_to = "Arc", values_to = "Correlation") %>%
  dplyr::mutate(dataset = "Primary") %>%
  bind_rows(as.data.frame(cor(nmf_rank_met, tern_arc_met[,paste0("Arc",1:3)])) %>%
              tibble::rownames_to_column(var = "LF") %>%
              tidyr::pivot_longer(!LF, names_to = "Arc", values_to = "Correlation") %>%
              dplyr::mutate(dataset = "Metastatic")) %>%
  bind_rows(as.data.frame(cor(nmf_rank_all, tern_arc_all[,paste0("Arc",1:3)])) %>%
              tibble::rownames_to_column(var = "LF") %>%
              tidyr::pivot_longer(!LF, names_to = "Arc", values_to = "Correlation") %>%
              dplyr::mutate(dataset = "Primary+Metastatic"))
write.table(sourcedata, paste0(basedir,"lisem/bc_landscape/github/SupplementaryFigure4f_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

sourcedata <- results %>% dplyr::rowwise() %>% dplyr::mutate(axis = ifelse(grepl("Arc",model),gsub("Arc","",model),gsub("LF","",model)))
write.table(sourcedata, paste0(basedir,"lisem/bc_landscape/github/SupplementaryFigure4g_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
