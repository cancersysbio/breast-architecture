### CREATE SUPPLEMENTARY FIGURE 3e-f #############################################################################
# creates supplementary figure 3e and 3f
# supplementary figure 3e provides the Pareto front model fits robustness statistics for each number of PCs/number of Archetypes combination for each architectural PCAs.
# supplementary figure 3f provides the Pareto front model fits robustness statistics across the architectural PCAs with k=3 archetypes within PC1-PC2 spaces.

### PREAMBLE #####################################################################################
library(yaml)
library(dplyr)
library(ggplot2)
library(ggpubr)
require(patchwork)
library(ggh4x)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

### MAIN ##########################################################################################
basedir = "/oak/stanford/groups/ccurtis2/users/"
outdir = paste0(basedir,"lisem/bc_landscape/revision/figures/")

## Load statistics
# see: /oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/revision/scripts/pareto_qc.R
load(file = paste0(basedir,"lisem/bc_landscape/revision/data/arc_ks.PCs_prim.RData"))
arc_ks.PCs_prim <- arc_ks.PCs
load(file = paste0(basedir,"lisem/bc_landscape/revision/data/arc_ks.PCs_met.RData"))
arc_ks.PCs_met <- arc_ks.PCs
load(file = paste0(basedir,"lisem/bc_landscape/revision/data/arc_ks.PCs_all.RData"))
arc_ks.PCs_all <- arc_ks.PCs

## Build table
arcs_ks.PCs.stats = bind_rows(lapply(1:length(arc_ks.PCs_prim), function(i) arc_ks.PCs_prim[[i]]$summary %>% dplyr::mutate(nb.PCs=i+1))) %>% 
  dplyr::mutate(nb.archetypes=factor(k), dataset = "Primary") %>%
  bind_rows(bind_rows(lapply(1:length(arc_ks.PCs_met), function(i) arc_ks.PCs_met[[i]]$summary %>% dplyr::mutate(nb.PCs=i+1))) %>% 
              dplyr::mutate(nb.archetypes=factor(k), dataset = "Metastatic")) %>%
  bind_rows(bind_rows(lapply(1:length(arc_ks.PCs_all), function(i) arc_ks.PCs_all[[i]]$summary %>% dplyr::mutate(nb.PCs=i+1))) %>% 
              dplyr::mutate(nb.archetypes=factor(k), dataset = "Discovery"))

# Looking at each fit
colfunc <- colorRampPalette(c("red","yellow","springgreen","royalblue"))
colfunc(9)

panelA <- as.data.frame(arcs_ks.PCs.stats) %>%
  dplyr::filter(dataset != "METABRIC") %>%
  dplyr::mutate(dataset = ifelse(dataset == "Discovery","Primary+Metastatic",dataset)) %>%
  ggplot(aes(x=nb.archetypes,y=varexpl,col=factor(nb.PCs))) + ylab("t-ratio") +
  geom_line(aes(group=factor(nb.PCs)), linewidth = 1) + geom_point(size = 3, alpha = 0.7) + xlab("Number of Archetypes (k)") + ylab("explained variance") +
  scale_color_manual(name = "Number of PCs (n)", values = colfunc(9)) + 
  facet_grid(. ~ factor(dataset, levels = c("Primary","Metastatic","Primary+Metastatic","METABRIC"))) +
  theme_LM + ggtitle("ParetoTI explained variance statistic of each fit") +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        title = element_text(size=16))

panelB <- arcs_ks.PCs.stats %>%
  dplyr::filter(nb.archetypes != 2 & dataset != "METABRIC") %>%
  dplyr::mutate(dataset = ifelse(dataset == "Discovery","Primary+Metastatic",dataset)) %>%
  ggplot(aes(x=nb.PCs,y=t_ratio,col=nb.archetypes)) + 
  geom_line(linewidth = 1) + geom_point(size = 3, alpha = 0.7) + xlab("Number of PCs (n)") + ylab("t-ratio") +
  scale_color_discrete(name = "Number of Archetypes (k)") + 
  facet_grid(. ~ factor(dataset, levels = c("Primary","Metastatic","Primary+Metastatic")), scales = "free_x", space = "free_x",
             labeller=label_wrap_gen(width = 10, multi_line = TRUE)) +
  scale_x_continuous(breaks = 2:10) +
  force_panelsizes(cols = c(1,1,1,0.5)) +
  theme_LM + ggtitle("ParetoTI t-ratio statistic of each fit") +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        title = element_text(size=16))

# Looking at 3 Arc 2 PCs only
p <- arcs_ks.PCs.stats %>%
  dplyr::filter(nb.PCs == 2 & nb.archetypes == 3 & dataset != "METABRIC") %>%
  dplyr::mutate(dataset = ifelse(dataset == "Discovery","Primary+Metastatic",dataset)) %>%
  dplyr::select(varexpl, t_ratio, dataset) %>%
  tidyr::pivot_longer(!dataset, names_to = "statistic", values_to = "value") %>%
  ggplot() +
  geom_col(aes(x=factor(dataset, levels = c("Primary","Metastatic","Primary+Metastatic")), y=value, fill=factor(statistic, levels = c("varexpl","t_ratio"))), position = "dodge") +
  xlab("") + ylab("Value") + scale_fill_manual(name = "", values = c("varexpl"="#58508d","t_ratio"="#ff6361"), labels = c("varexpl"="explained variance","t_ratio"="t-ratio")) +
  theme_LM + ggtitle("ParetoTI statistics value for 3 archetypes in PC1-PC2") +
  theme(axis.text.x = element_text(size=12),
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        title = element_text(size=16))

### SAVE ##########################################################################################
## Figure----
png(filename = file.path(outdir, 'SupplementaryFigure3e.png'), res = 300, width = 9.46, height = 6.84, units = 'in')
panelA / panelB
dev.off()

png(filename = file.path(outdir, 'SupplementaryFigure3f.png'), res = 300, width = 8.94, height = 3.88, units = 'in')
p
dev.off()

# SourceData----
sourcetable = as.data.frame(arcs_ks.PCs.stats) %>%
  dplyr::mutate(dataset = ifelse(dataset == "Discovery","Primary+Metastatic",dataset))
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/revision/data/SupplementaryFigure3e_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

sourcetable = arcs_ks.PCs.stats %>%
  dplyr::filter(nb.PCs == 2 & nb.archetypes == 3 & dataset != "METABRIC") %>%
  dplyr::mutate(dataset = ifelse(dataset == "Discovery","Primary+Metastatic",dataset)) %>%
  dplyr::select(varexpl, t_ratio, dataset) %>%
  tidyr::pivot_longer(!dataset, names_to = "statistic", values_to = "value")
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/revision/data/SupplementaryFigure3f_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")



