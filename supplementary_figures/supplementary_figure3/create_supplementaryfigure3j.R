### CREATE SUPPLEMENTARY FIGURE 3J #############################################################################
# BRCA1 vs BRCA2 HRD in primary high risk tumors

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

# Plot by subtype-------
plot_data <- plot_data %>% inner_join(primary_megatable %>% select(Sample, ENiClust)) %>% filter(group == "ER+ High", CHORD != "cannot_be_determined") %>%
  mutate(Class = ifelse(CHORD == "HR_proficient", "HR proficient", ifelse(p_BRCA1 > p_BRCA2, "BRCA1-like", "BRCA2-like")), 
         Class = factor(Class), ENiClust = factor(gsub("ic", "IC", ENiClust))) %>% 
  group_by(ENiClust, Class) %>% summarize(n = n()) %>% group_by(ENiClust) %>% mutate(total = sum(n), Proportion = n/total) %>% data.frame()

plot_data[nrow(plot_data) + 1, ] <- c("IC1", "BRCA1-like", 0, 32, 0.0)
plot_data[nrow(plot_data) + 1, ] <- c("IC2", "BRCA1-like", 0, 19, 0.0)
plot_data[nrow(plot_data) + 1, ] <- c("IC2", "BRCA2-like", 0, 19, 0.0)
plot_data[nrow(plot_data) + 1, ] <- c("IC6", "BRCA1-like", 0, 23, 0.0)

plot_data <- plot_data %>% arrange(ENiClust, Class)

plot_data$amp_code <- NA
plot_data[plot_data$Class == 'BRCA1-like','amp_code'] <- 'd'
plot_data[plot_data$Class == 'BRCA2-like','amp_code'] <- 'e'
plot_data[plot_data$Class == 'HR proficient','amp_code'] <- 'f'

plot_data$gindex <- paste(plot_data$ENiClust, plot_data$amp_code, sep ='_')
plot_data$Proportion <- as.numeric(plot_data$Proportion)
plot_data$n <- as.numeric(plot_data$n)

plot_data <- tibble(plot_data)

# create barplot
p1 <- create.barplot(
  Proportion ~ ENiClust,
  groups = gindex,
  col = c(
    '#FF5500','#FF550088','white',
    '#00EE76','#00EE7688','white',
    '#FFFF40','#FFFF4088','white',
    '#EE82EE','#EE82EE88','white'
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
  n ~ ENiClust,
  groups = gindex,
  col = c(
    '#FF5500','#FF550088','white',
    '#00EE76','#00EE7688','white',
    '#FFFF40','#FFFF4088','white',
    '#EE82EE','#EE82EE88','white'
  ),
  stack = TRUE,
  #xaxis.lab = rep('', 100),
  #yaxis.lab = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
  #col = c('#8fdab0','#9ebff1ff','#fe85cfff'),
  data = plot_data,
  ylimits = c(0,45),
  height = 8,
  width = 8,
  ylab.label = 'Number of Samples',
  xlab.label = '',
  #stack= TRUE,
  resolution = 300, 
  xaxis.rot = 0, 
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
      x = 0.6,
      y = 0.96,
      corner = c(1,1)
    )))

## Save plot ----
create.multipanelplot(
  list(p1,p2),
  filename = file.path(main_repo_path, "plots", "SupplementaryFigure3J.pdf"),
  resolution = 300)


### SAVE ##########################################################################################
plot_data %>% 
  select(ENiClust, Class, n, Total = total, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "SupplementaryFigure3j_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
