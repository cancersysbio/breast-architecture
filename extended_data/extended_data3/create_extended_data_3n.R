### CREATE EXTENDED DATA FIGURE 3O #############################################################################
# Ternary of IDC/ILC (Primary WGS)

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(yaml)
library(Ternary)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
tern <- read.delim(file.path(main_repo_path, "data", "ExtendedData3o_sourcetable.txt"), sep = ",")
tern_arc <- tern %>% 
  filter(colorby == "ER+ Typical") %>%
  filter(histo_type %in% c("ILC", "IDC"))

pdf(file.path(main_repo_path, "plots", "Figure2f.pdf"),h=3*2, w=3*2)
TernaryPlot(point = "up", 
            # alab = "Non amp-associated SVs",
            # blab = "Diploidy - Amp-associated SV", 
            # clab = "WGD",
            lab.cex = 0.8, grid.minor.lines = 0,
            grid.lines = 5,
            lab.offset = 0.18,
            grid.lty = 'solid', col = rgb(0.9, 0.9, 0.9), grid.col = 'white', 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            axis.labels = F, 
            axis.tick = F,
            padding = 0.12)

TernaryLines(list(c(0, 255, 0), rep(50,3)), col = 'grey')
TernaryLines(list(c(0, 0, 255), rep(50,3)), col = 'grey')
TernaryLines(list(c(255, 0, 0), rep(50,3)), col = 'grey')

AddToTernary(points, tern_arc[which(tern_arc$histo_type == "IDC"),c(3,2,1)], pch = 21, cex = 2, bg = alpha("#a3cfffff", 0.4), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$histo_type == "ILC"),c(3,2,1)], pch = 21, cex = 2, bg = alpha("#335e8dff", 0.7), col = "white")

legend('topright', 
       legend = c('IDC', 'ILC'),
       col="white",
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,
       pt.bg = c(alpha("#a3cfffff", 0.7), alpha("#335e8dff", 0.7)))

data_points <- list(
  "1" = c(0, 0, 255),
  "2" = c(255, 0, 0),
  "3" = c(0, 255, 0))
xy <- CoordinatesToXY(data_points)
points(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), pch = 19, cex = 2.8, col = alpha(c("purple", "blue", "brown"),0.7))
text(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), names(data_points), cex = 0.8, font = 2)
dev.off()
