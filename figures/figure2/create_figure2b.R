### CREATE FIGURE 2b #############################################################################
# creates figure 2b
# figure 2b provides the Pareto front projections on ternary plot of copy number and SV signature profiles from primary and metastatic tumors independently.

### PREAMBLE #####################################################################################
library(yaml)
require(tibble)
library(tidyverse)
library(ParetoTI, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
library(Ternary, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
col.ic.group = c("ER+ High" = "#fa954eff", "ER+ Typical" = "#a1c1d4ff", "HER2+" = "#8b0100ff",
                 "IC4ER-" = "#c2b7dbff", "IC10" = "#7c26ccff")

### MAIN ##########################################################################################
basedir <- "/oak/stanford/groups/ccurtis2/users/"
outdir <- paste0(basedir,"lisem/bc_landscape/revision/figures/")

## PCA Primary-only----
## Load PCs
pca_prim = read.table(paste0(basedir, "khoulaha/BreastLandscape/data/primary_pcas_cna_sv_signatures_proportions.txt"), header = T, sep = "\t")
rownames(pca_prim) = pca_prim$Sample

## Load Megatable Primary
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")

## Run Pareto
# WARNING: each k_fit_pch() run will give slightly different approximations for the Archetypes' positions leading to slightly different samples' coordinates on the Ternary plot
# if you want to use the results at the running time, please refer to the related Source Data table
arc_ks.LFs.noboot = k_fit_pch(t(pca_prim[,1:2]), ks = 3:5, check_installed = T, bootstrap = F, volume_ratio = "t_ratio") 

tern_arc = as_tibble(t(arc_ks.LFs.noboot$pch_fits$S[[1]]))
tern_arc$Sample = rownames(t(arc_ks.LFs.noboot$pch_fits$S[[1]]))

## Extract Arc positions
arc_noboot = fit_pch(t(pca_prim[,1:2]), noc = 3, delta = 0) # would need to redo it with bootstrapping
arc = as.data.frame(arc_noboot$XC)

## Projection
tern_arc$ENiClust = sapply(tern_arc$Sample, function(x) mega_1$ENiClust[which(mega_1$Sample == x)])
tern_arc$Subgroup = sapply(tern_arc$Sample, function(x) mega_1$group[which(mega_1$Sample == x)])

## Identify the Arcs
colnames(tern_arc)[1:3] = c("Arc1", "Arc2", "Arc3")

ggplot(tern_arc, aes(x=Arc1,y=Arc2,color=Subgroup)) +
  geom_point() + scale_color_manual(values=col.ic.group)

colnames(tern_arc)[1:3] = c("Arc1", "Arc3", "Arc2")
colnames(arc) = c("Arc1","Arc3","Arc2")

tern_arc_prim = tern_arc
arc_prim = arc

## PCA Metastatic-only----
## Load PCs
pca_mets = read.table(paste0(basedir, "khoulaha/BreastLandscape/data/metastatic_pcas_cna_sv_signatures_proportions.txt"), header = T, sep = "\t")
rownames(pca_mets) = pca_mets$Sample

## Load Megatable Primary
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")

## Run Pareto
# WARNING: each k_fit_pch() run will give slightly different approximations for the Archetypes' positions leading to slightly different samples' coordinates on the Ternary plot
# if you want to use the results at the running time, please refer to the related Source Data table
arc_ks.LFs.noboot = k_fit_pch(t(pca_mets[,1:2]), ks = 3:5, check_installed = T, bootstrap = F, volume_ratio = "t_ratio")

tern_arc = as_tibble(t(arc_ks.LFs.noboot$pch_fits$S[[1]]))
tern_arc$Sample = rownames(t(arc_ks.LFs.noboot$pch_fits$S[[1]]))

## Extract Arc positions
arc_noboot = fit_pch(t(pca_mets[,1:2]), noc = 3, delta = 0) 
arc = as.data.frame(arc_noboot$XC)

## Projection
tern_arc$ENiClust = sapply(tern_arc$Sample, function(x) mega_2$ENiClust[which(mega_2$Sample == x)])
tern_arc$Subgroup = sapply(tern_arc$Sample, function(x) mega_2$group[which(mega_2$Sample == x)])

## Identify the Arcs
colnames(tern_arc)[1:3] = c("Arc1", "Arc2", "Arc3")

ggplot(tern_arc, aes(x=Arc1,y=Arc2,color=Subgroup)) +
  geom_point() + scale_color_manual(values=col.ic.group)

colnames(tern_arc)[1:3] = c("Arc1", "Arc3", "Arc2")
colnames(arc) = c("Arc1","Arc3","Arc2")

tern_arc_met = tern_arc
arc_met = arc

rm(tern_arc, arc)

### SAVE ##########################################################################################
## Figure----
# Primary
svg(paste0(outdir,"Figure2b_top.svg"),h=3*2,w=3*2)
TernaryPlot(point = "up", axis.labels = F, axis.tick = F,
            lab.cex = 0.8, grid.minor.lines = 0,
            grid.lines = 5,
            lab.offset = 0.18,
            grid.lty = 'solid', col = rgb(0.9, 0.9, 0.9), grid.col = 'white', 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            padding = 0.12)

TernaryLines(list(c(0, 255, 0), rep(50,3)), col = 'grey')
TernaryLines(list(c(0, 0, 255), rep(50,3)), col = 'grey')
TernaryLines(list(c(255, 0, 0), rep(50,3)), col = 'grey')

for(g in c("ER+ High","ER+ Typical","HER2+","IC10","IC4ER-")){
  AddToTernary(points, tern_arc_prim[which(tern_arc_prim$Subgroup == g),rev(c("Arc1", "Arc3", "Arc2"))], pch = 21, cex = 2, bg = alpha(col.ic.group[names(col.ic.group) == g], 0.7), col = "white")
}

legend('bottomright', 
       legend = c("ER+ High","ER+ Typical","HER2+","IC10","IC4ER-"),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8, y.intersp = 0.2,
       pt.bg = sapply(c("ER+ High","ER+ Typical","HER2+","IC10","IC4ER-"), function(g) alpha(col.ic.group[names(col.ic.group) == g], 0.7)), col="white")

data_points <- list(
  "1" = c(0, 0, 255),
  "2" = c(255, 0, 0), 
  "3" = c(0, 255, 0))

xy <- CoordinatesToXY(data_points)
points(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), pch = 19, cex = 2.8, col = alpha(c("purple", "blue", "brown"),0.7))
text(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), names(data_points), cex = 0.8, font = 2)

arc_name = c("TNBC -enriched", "ER+ Typical -enriched", "ER+ High/HER2+ -enriched")
text(xy[1, ] + c(-sqrt(0.1**2),0,sqrt(0.1**2)), xy[2, ]+c(-sqrt(0.1**2),0.2,-sqrt(0.1**2)), arc_name, cex = 0.8, font = 2)
dev.off()

# Metastatic
svg(paste0(outdir,"Figure2b_bottom.svg"),h=3*2,w=3*2)
TernaryPlot(point = "up", axis.labels = F, axis.tick = F,
            lab.cex = 0.8, grid.minor.lines = 0,
            grid.lines = 5,
            lab.offset = 0.18,
            grid.lty = 'solid', col = rgb(0.9, 0.9, 0.9), grid.col = 'white', 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            padding = 0.12)

TernaryLines(list(c(0, 255, 0), rep(50,3)), col = 'grey')
TernaryLines(list(c(0, 0, 255), rep(50,3)), col = 'grey')
TernaryLines(list(c(255, 0, 0), rep(50,3)), col = 'grey')

for(g in c("ER+ High","ER+ Typical","HER2+","IC10","IC4ER-")){
  AddToTernary(points, tern_arc_met[which(tern_arc_met$Subgroup == g),rev(c("Arc1", "Arc3", "Arc2"))], pch = 21, cex = 2, bg = alpha(col.ic.group[names(col.ic.group) == g], 0.7), col = "white")
}

legend('bottomright', 
       legend = c("ER+ High","ER+ Typical","HER2+","IC10","IC4ER-"),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8, y.intersp = 0.2,
       pt.bg = sapply(c("ER+ High","ER+ Typical","HER2+","IC10","IC4ER-"), function(g) alpha(col.ic.group[names(col.ic.group) == g], 0.7)), col="white")

data_points <- list(
  "1" = c(0, 0, 255),
  "2" = c(255, 0, 0), 
  "3" = c(0, 255, 0))

xy <- CoordinatesToXY(data_points)
points(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), pch = 19, cex = 2.8, col = alpha(c("purple", "blue", "brown"),0.7))
text(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), names(data_points), cex = 0.8, font = 2)

arc_name = c("TNBC -enriched", "ER+ Typical -enriched", "ER+ High/HER2+ -enriched")
text(xy[1, ] + c(-sqrt(0.1**2),0,sqrt(0.1**2)), xy[2, ]+c(-sqrt(0.1**2),0.2,-sqrt(0.1**2)), arc_name, cex = 0.8, font = 2)
dev.off()

## SourceData---- 
# Primary
tern_arc_prim$ENiClust[which(tern_arc_prim$ENiClust == "Other")] = "ic3/ic7"
tern_arc_prim = tern_arc_prim[,c(4,1,3,2,5,6)]
write.table(tern_arc_prim, paste0(basedir,"lisem/bc_landscape/submission/data/Figure2btop_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
arc_prim$PC = rownames(arc_prim)
arc_prim = arc_prim[,c(4,1,3,2)]
write.table(arc_prim, paste0(basedir,"lisem/bc_landscape/submission/data/Figure2btop_ArcPositions.txt"),row.names = T,col.names = T,quote = F, sep="\t")

# Metastatic
tern_arc_met$ENiClust[which(tern_arc_met$ENiClust == "Other")] = "ic3/ic7"
tern_arc_met = tern_arc_met[,c(4,1,3,2,5,6)]
write.table(tern_arc_met, paste0(basedir,"lisem/bc_landscape/revision/data/Figure2bbottom_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
arc_met$PC = rownames(arc_met)
arc_met = arc_met[,c(4,1,3,2)]
write.table(arc_met, paste0(basedir,"lisem/bc_landscape/revision/data/Figure2bbottom_ArcPositions.txt"),row.names = T,col.names = T,quote = F, sep="\t")

