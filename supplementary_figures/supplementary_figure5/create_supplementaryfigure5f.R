### SUPPLEMENTARY FIGURE 5F #####################################################################
# create scatterplot of age and clock-like mutations

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyverse)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# set data path
age_dt_path <- file.path(main_repo_path, 'data', "SupplementaryFigure5f_sourcetable.txt")

# read in data
age_dt <- read.delim(age_dt_path, as.is = TRUE)
age_dt <- age_dt[which(age_dt$group != 'IC4ER-'),]
age_dt$total <- age_dt$SBS1+age_dt$SBS5+age_dt$SBS40

col <- rep(NA, nrow(age_dt))
col[which(age_dt$group == 'ER+ High')] <- '#fa954eff'
col[which(age_dt$group == 'ER+ Typical')] <- '#a1c1d4ff'
col[which(age_dt$group == 'HER2+')] <- '#8b0100ff'
col[which(age_dt$group == 'IC10')] <- '#7c26ccff'

create.scatterplot(
  total ~ age,
  data = age_dt,
  filename = 'supplementary_figure5f.pdf',
  ylab.label = 'Number of clock-like\nmutations',
  col = col,
  ylimits = c(0,15000),
  yat = seq(0,15e3,5e3),
  key = list(
             text = list(
                 lab = c('ER+ High','ER+ Typical','HER2+','IC10'),
                 cex = 1,
                 col = 'black'
                 ),
             points = list(
                 pch = 19,
                 col = c('#fa954eff','#a1c1d4ff','#8b0100ff','#7c26ccff'),
                 cex = 1
                 ),
             x = 0.01,
             y = 0.99,
             corner = c(0,1),
             padding.text = 2
             ),
  yaxis.lab = c(expression(0),expression(5^3), expression(10^3), expression(15^3)),
  xlab.label = 'Age',
  resolution = 300
  )
