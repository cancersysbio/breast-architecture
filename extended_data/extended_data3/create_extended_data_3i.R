### CREATE EXTENDED DATA FIGURE 3i #############################################################################
# creates extended data figure 3i
# extended data figure 3i provides the difference in position on the Pareto fronts between the centroid of primary samples and the centroid of metastatic samples in each IC group.

### PREAMBLE #####################################################################################
library(yaml)
library(Ternary, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
library(dplyr)
library(scales)
library(sda, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.2.2-Seurat/")

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

## Load Arc proportions
tern_arc = read.table(paste0(basedir,"lisem/bc_landscape/github/ExtendedData3g_sourcetable.txt"), header = T, sep="\t") #hasn't been saved in submission

cent_prim <- centroids(as.matrix(tern_arc[which(!is.na(tern_arc$Subgroup) & tern_arc$Stage == "Primary"),c("Arc1","Arc2","Arc3")]), tern_arc$Subgroup[which(!is.na(tern_arc$Subgroup) & tern_arc$Stage == "Primary")])
cent_mets <- centroids(as.matrix(tern_arc[which(!is.na(tern_arc$Subgroup) & tern_arc$Stage == "Metastatic"),c("Arc1","Arc2","Arc3")]), tern_arc$Subgroup[which(!is.na(tern_arc$Subgroup) & tern_arc$Stage == "Metastatic")])

### SAVE ##########################################################################################
## Figure----
svg(paste0(outdir,"ExtendedData3i.svg"),h=3*2,w=3*2)
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
  TernaryArrows(cent_prim$means[,g][c(2,3,1)], cent_mets$means[,g][c(2,3,1)], length = 0.2, col = col.ic.group[names(col.ic.group) == g])
}
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
sourcetable = as.data.frame(cent_prim$means)
sourcetable$Arc = rownames(sourcetable)
sourcetable = sourcetable[,c(7,1:6)]
sourcetable = sourcetable %>% pivot_longer(!Arc,names_to="Subgroup", values_to="Centroid_means")
sourcetable$Stage = "Primary"

sourcetable_bis = as.data.frame(cent_mets$means)
sourcetable_bis$Arc = rownames(sourcetable_bis)
sourcetable_bis = sourcetable_bis[,c(7,1:6)]
sourcetable_bis = sourcetable_bis %>% pivot_longer(!Arc,names_to="Subgroup", values_to="Centroid_means")
sourcetable_bis$Stage = "Metastatic"

sourcetable = rbind(sourcetable, sourcetable_bis)
write.table(as.data.frame(sourcetable), paste0(basedir,"lisem/bc_landscape/github/ExtendedData3i_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
