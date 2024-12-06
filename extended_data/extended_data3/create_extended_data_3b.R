### CREATE EXTENDED DATA 3B #######################################################################
# create ternary plot with BRCA-like tumors labeled

### PREAMEBLE #####################################################################################
require(tibble)
library(tidyverse)
library(Ternary)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
# read in plot data 
tern_arc <- read.delim(
  file.path(main_repo_path, 'data', 'ExtendedData3b_sourcetable.txt'),
  as.is = TRUE
  )

res <- list()
for (i in c('ER+ High','HER2+','ER+ Typical','IC10')) {
  stats <- fisher.test(table(tern_arc[tern_arc$group == i,c('brcalike','ecdna')]))
  res[[i]] <- data.frame(
    subgroup = i,
    or = stats$estimate[[1]],
    l95 = stats$conf.int[[1]],i
    u95 = stats$conf.int[[2]],
    p = stats$p.value
    )
}
res <- do.call(rbind, res)


## Ternary
pdf('extended_data3b.pdf',height=3*2,width=3*2)
TernaryPlot(point = "up", alab = "Non amp-associated SVs",
            blab = "Diploidy - Amp-associated SV", clab = "WGD",
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

AddToTernary(points, tern_arc[which(tern_arc$Subgroup == "ER+ High"),c(3,4,2)], pch = 21, cex = 2, bg = alpha("grey", 0.7), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$Subgroup == "ER+ Typical"),c(3,4,2)], pch = 21, cex = 2, bg = alpha("grey", 0.7), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$Subgroup == "HER2+"),c(3,4,2)], pch = 21, cex = 2, bg = alpha("grey", 0.7), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$Subgroup == "IC10"),c(3,4,2)], pch = 21, cex = 2, bg = alpha("grey", 0.7), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$Subgroup == "IC4ER-"),c(3,4,2)], pch = 21, cex = 2, bg = alpha("grey", 0.7), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$brcalike == 1),c(3,4,2)], pch = 21, cex = 2, bg = alpha("black", 0.7), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$brcalike == 1 & tern_arc$ER == 1),c(3,4,2)], pch = 21, cex = 2, bg = alpha("pink", 0.7), col = "white")

legend('topright', 
       legend = c('None', 'ER- HRD-like','ER+ HRD-like'),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,
       pt.bg = c(alpha("grey", 0.7), alpha("black", 0.7), alpha("pink", 0.7)), col="white")

data_points <- list(
  "1" = c(0, 0, 255),
  "2" = c(255, 0, 0), 
  "3" = c(0, 255, 0))
xy <- CoordinatesToXY(data_points)
points(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), pch = 19, cex = 2.8, col = alpha(c("purple", "blue", "brown"),0.7))
text(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), names(data_points), cex = 0.8, font = 2)

arc_name = c("TNBC -enriched", "ER+ Typical -enriched", "ER+ High/HER2+ -enriched")
text(xy[1, ] + c(-sqrt(0.01**2),0,sqrt(0.1**2)), xy[2, ]+c(-sqrt(0.15**2),0.2,-sqrt(0.15**2)), arc_name, cex = 0.8, font = 2)
dev.off()

tern_arc$tnbc <- (tern_arc$dominant == 'Arc1')*1
fisher.test(table(tern_arc[tern_arc$ER == 1,c('tnbc','brcalike')]))
fisher.test(table(tern_arc[which(tern_arc$Subgroup == 'ER+ High'),c('tnbc','brcalike')]))
fisher.test(table(tern_arc[which(tern_arc$Subgroup == 'ER+ Typical'),c('tnbc','brcalike')]))

### SAVE ##########################################################################################
# save source data
write.table(
  tern_arc,
  file = file.path(main_repo_path, 'data', 'ExtendedData3b_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )
