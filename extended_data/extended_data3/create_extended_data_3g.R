### CREATE EXTENDED DATA FIGURE 3g #############################################################################
# creates extended data figure 3g
# extended data figure 3g provides the distribution of primary and metastatic tumors on Pareto front.

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

### PCA All (Primary/Metastatic)
## Load PCs
pca_all = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/primary_metastatic_pcas_cna_sv_signatures_proportions.txt"), header = T, sep = "\t")
rownames(pca_all) = pca_all$Sample

## Load Megatables
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])
mega$Sample_Type = NA
mega$Sample_Type[which(mega$Sample %in% mega_1$Sample)] = "Primary"
mega$Sample_Type[which(mega$Sample %in% mega_2$Sample)] = "Metastatic"

## Run Pareto
# WARNING: each k_fit_pch() run will give slightly different approximations for the Archetypes' positions leading to slightly different samples' coordinates on the Ternary plot
# if you want to use the results at the running time, please refer to the related Source Data table
arc_ks.LFs.noboot = k_fit_pch(t(pca_all[,1:2]), ks = 3:5, check_installed = T, bootstrap = F, volume_ratio = "t_ratio" )

tern_arc = as_tibble(t(arc_ks.LFs.noboot$pch_fits$S[[1]]))
tern_arc$Sample = rownames(t(arc_ks.LFs.noboot$pch_fits$S[[1]]))

## Extract Arc positions
arc_noboot = fit_pch(t(pca_all[,1:2]), noc = 3, delta = 0) # would need to redo it with bootstrapping
arc = as.data.frame(arc_noboot$XC)
colnames(arc) = c("Arc1","Arc2","Arc3")

## Projection
all(tern_arc$Sample %in% mega$Sample)
tern_arc$ENiClust = sapply(tern_arc$Sample, function(x) mega$ENiClust[which(mega$Sample == x)])
tern_arc$Subgroup = sapply(tern_arc$Sample, function(x) mega$group[which(mega$Sample == x)])
tern_arc$Stage = sapply(tern_arc$Sample, function(x) mega$Sample_Type[which(mega$Sample == x)])

## Identify the Arcs
colnames(tern_arc)[1:3] = c("Arc1", "Arc2", "Arc3")

ggplot(tern_arc, aes(x=Arc1,y=Arc2,color=Subgroup)) +
  geom_point() + scale_color_manual(values=col.ic.group)

colnames(tern_arc)[1:3] = c("Arc1", "Arc3", "Arc2")
colnames(arc) = c("Arc1", "Arc3", "Arc2")

### SAVE ##########################################################################################
## Figure----
svg(paste0(basedir,"lisem/bc_landscape/github/ExtendedData3g_Stage.svg"),h=3*2,w=3*2)
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

AddToTernary(points, tern_arc[which(tern_arc$Stage == "Primary"),rev(c("Arc1", "Arc3", "Arc2"))], pch = 21, cex = 2, bg = alpha("#77a4ecff", 0.7), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$Stage == "Metastatic"),rev(c("Arc1", "Arc3", "Arc2"))], pch = 21, cex = 2, bg = alpha("#fe53bcff", 0.7), col = "white")

legend('bottomright', 
       legend = c('Primary', 'Metastasis'),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,y.intersp = 0.2,
       pt.bg = c(alpha("#77a4ecff", 0.7), alpha("#fe53bcff", 0.7)), col="white")

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
tern_arc$ENiClust[which(tern_arc$ENiClust == "Other")] = "ic3/ic7"
tern_arc = tern_arc[,c(4,3,1,2,5,6,7)]
write.table(tern_arc, paste0(basedir,"lisem/bc_landscape/github/ExtendedData3g_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
arc$PC = rownames(arc)
arc = arc[,c(4,1,3,2)]
write.table(arc, paste0(basedir,"lisem/bc_landscape/github/ExtendedData3g_ArcPositions.txt"),row.names = T,col.names = T,quote = F, sep="\t")
