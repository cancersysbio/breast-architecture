### CREATE EXTENDED DATA 3D #############################################################################
# BRCA1 vs BRCA2 HRD in primary tumors
### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

plot_data <- read.delim(file.path(main_repo_path, "data", "CHORD_HRD_predictions.txt"))
primary_megatable <- read.delim(file.path(main_repo_path, 'data','primary_megatable.txt'))

# Plot by groups-------
plot_data <- plot_data %>% filter(Stage == "Primary", CHORD != "cannot_be_determined", !is.na(group)) %>%
  mutate(Class = ifelse(CHORD == "HR_proficient", "HR proficient", ifelse(p_BRCA1 > p_BRCA2, "BRCA1-like", "BRCA2-like"))) %>% 
  group_by(group, Class) %>% summarize(n = n()) %>% group_by(group) %>% mutate(total = sum(n), Proportion = n/total) %>% data.frame()

plot_data[nrow(plot_data) + 1, ] <- c("HER2+", "BRCA2-like", 0, 55, 0.0)

plot_data <- plot_data %>% arrange(group, Class)

plot_data$amp_code <- NA
plot_data[plot_data$Class == 'BRCA1-like','amp_code'] <- 'd'
plot_data[plot_data$Class == 'BRCA2-like','amp_code'] <- 'e'
plot_data[plot_data$Class == 'HR proficient','amp_code'] <- 'f'

plot_data$gindex <- paste(plot_data$group, plot_data$amp_code, sep ='_')
plot_data$Proportion <- as.numeric(plot_data$Proportion)
plot_data$n <- as.numeric(plot_data$n)

plot_data <- tibble(plot_data)

# create barplot
p1 <- create.barplot(
  Proportion ~ group,
  groups = gindex,
  col = c(
    '#fa954eff','#fccaa7','white',
    '#a1c1d4ff','#d0e0e9','white',
    '#8b0100ff','#c58087','white',
    '#7c26ccff','#7c26cc6c','white',
    '#c2b7dbff','#e6e2ff','white'
  ),
  stack = TRUE,
  xaxis.lab = rep('', 100),
  yaxis.lab = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
  #col = c('#8fdab0','#9ebff1ff','#fe85cfff'),
  data = plot_data,
  ylimits = c(0,1),
  height = 8,
  width = 8,
  ylab.label = 'Proportion of Samples',
  xlab.label = '',
  #stack= TRUE,
  resolution = 300
)

p2 <- create.barplot(
  n ~ group,
  groups = gindex,
  col = c(
    '#fa954eff','#fccaa7','white',
    '#a1c1d4ff','#d0e0e9','white',
    '#8b0100ff','#c58087','white',
    '#7c26ccff','#7c26cc6c','white',
    '#c2b7dbff','#e6e2ff','white'
  ),
  stack = TRUE,
  #xaxis.lab = rep('', 100),
  #yaxis.lab = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
  #col = c('#8fdab0','#9ebff1ff','#fe85cfff'),
  data = plot_data,
  ylimits = c(0,275),
  height = 8,
  width = 8,
  ylab.label = 'Number of Samples',
  xlab.label = '',
  #stack= TRUE,
  resolution = 300, 
  xaxis.rot = 45, 
  legend = list(
    inside = list(
      fun = draw.key,
      args = list(
        key = list(
          points = list(
            col = 'black',
            pch = 22,
            cex = 3,
            fill = c('black','darkgrey')
          ),
          text = list(
            lab = c('BRCA1-like','BRCA2-like')
          ),
          #title = 'Metastatic',
          padding.text = 5,
          cex = 1.2
        )
      ),
      # Positioning legend on plot
      x = 0.98,
      y = 0.96,
      corner = c(1,1)
    )))

## Save plot ----
create.multipanelplot(
  list(p1,p2),
  filename = file.path(main_repo_path, "plots", "ExtendedData3D.pdf"),
  resolution = 300)

### SAVE ##########################################################################################
plot_data %>% 
  select(Subtype = group, Class, n, Total = total, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData3d_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
