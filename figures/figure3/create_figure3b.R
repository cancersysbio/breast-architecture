### CREATE FIGURE 3B ##############################################################################
# create boxplot of mutation timing estimates across subgroups
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
#library(tidyverse)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
#data path
pri_amp_path <- file.path(main_repo_path, "data", "pri_amp_dt.rds")
meta_amp_path <- file.path(main_repo_path, "data", "meta_amp_dt.rds")

#No. of pre-amplified mutations
pri_amp_dt <- as.data.frame(readRDS(pri_amp_path))
meta_amp_dt <- as.data.frame(readRDS(meta_amp_path))

# read in megatable 
megatable <- read.delim(
  file.path(main_repo_path, 'data','primary_megatable.txt'),
  as.is = TRUE
  )
megatable <- megatable[which(megatable$group != ''),]
pri_amp_dt <- merge(pri_amp_dt[,-which(colnames(pri_amp_dt) == 'group')], megatable[,c('Sample','group')], by = 'Sample')
# remove samples where AA failed 
failed_samples <- megatable[which(is.na(megatable$AA_ecdna)),'Sample']
pri_amp_dt <- pri_amp_dt[which(!pri_amp_dt$Sample %in% failed_samples),]

# read in megatable 
meta_megatable <- read.delim(
  file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
  as.is = TRUE
  )
meta_megatable <- meta_megatable[which(meta_megatable$group != ''),]
meta_amp_dt <- merge(meta_amp_dt[,-which(colnames(meta_amp_dt) == 'group')], meta_megatable[,c('Sample','group')], by = 'Sample')

pplot_data <- pri_amp_dt[which(!is.na(pri_amp_dt$group) & pri_amp_dt$group != 'IC4ER-'),]
mplot_data <- meta_amp_dt[which(!is.na(meta_amp_dt$group) & meta_amp_dt$group != 'IC4ER-'),]
pcol <- rep(NA, nrow(pplot_data))
pcol[which(pplot_data$group == 'ER+ High')] <- '#fa954eff'
pcol[which(pplot_data$group == 'ER+ Typical')] <- '#a1c1d4ff'
pcol[which(pplot_data$group == 'HER2+')] <- '#8b0100ff'
pcol[which(pplot_data$group == 'IC10')] <- '#7c26ccff'
#col[which(pri_amp_dt$group == 'IC4ER-')] <- '#c2b7dbff'

mcol <- rep(NA, nrow(mplot_data))
mcol[which(mplot_data$group == 'ER+ High')] <- '#fa954eff'
mcol[which(mplot_data$group == 'ER+ Typical')] <- '#a1c1d4ff'
mcol[which(mplot_data$group == 'HER2+')] <- '#8b0100ff'
mcol[which(mplot_data$group == 'IC10')] <- '#7c26ccff'

# test difference in pre mutation no in ER+ High vs IC10 in primary tumours
perhigh <- wilcox.test(
  pplot_data[pplot_data$group == 'ER+ High','d_pre_nokat_clock'],
  pplot_data[pplot_data$group == 'IC10','d_pre_nokat_clock']
  )
print("Wilcox test of pre mutation number between ER+ High vs IC10 in primary tumours")
print(perhigh)
# test difference in pre mutation no in HER2+ vs IC10 in primary tumours
pher2 <- wilcox.test(
  pplot_data[pplot_data$group == 'HER2+','d_pre_nokat_clock'],
  pplot_data[pplot_data$group == 'IC10','d_pre_nokat_clock']
  )
print("Wilcox test of pre mutation number between HER2+ vs IC10 in primary tumours")
print(pher2)

# test difference in pre mutation no in ER+ High vs IC10 in metastatic tumours
merhigh <- wilcox.test(
  mplot_data[mplot_data$group == 'ER+ High','d_pre_nokat_clock'],
  mplot_data[mplot_data$group == 'IC10','d_pre_nokat_clock']
  )
print("Wilcox test of pre mutation number between ER+ High vs IC10 in metastatic tumours")
print(merhigh)
# test difference in pre mutation no in HER2+ vs IC10 in metastatic tumours
mher2 <- wilcox.test(
  mplot_data[mplot_data$group == 'HER2+','d_pre_nokat_clock'],
  mplot_data[mplot_data$group == 'IC10','d_pre_nokat_clock']
  )
print("Wilcox test of pre mutation number between HER2+ vs IC10 in metastatic tumours")
print(mher2)


#### PLOT #########################################################################################
primary <- create.boxplot(
  d_pre_nokat_clock ~ group,
  data = pplot_data,
  add.stripplot = TRUE,
  points.col = pcol,
  ylab.label = 'Density of SNVs\nbefore amplification',
  xaxis.lab = rep('', 4),
  xlab.label = '',
  add.text = TRUE,
  main = 'Primary',
  main.cex = 1.5,
  text.labels = c( 
    expression('P'['ER+ High vs IC10']*'=1.3x10'^'-6'),
    expression('P'['HER2+ vs IC10']*'=7.8x10'^'-7')
    ),
  text.x = c(1.1,1.05),
  text.y = c(1.8,1.65),
  text.cex = c(1.2,1.2),
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  ylimits = c(0,2),
  yat = seq(0,3,1),
  resolution = 300
  )
metastatic <- create.boxplot(
  d_pre_nokat_clock ~ group,
  data = mplot_data,
  add.stripplot = TRUE,
  points.col = mcol,
  ylab.label = 'Density of SNVs\nbefore amplification',
  xaxis.lab = c('ER+ High','ER+ Typical','HER2+','IC10'),
  xlab.label = 'Subtype',
  add.text = TRUE,
  main = 'Metastatic',
  main.cex = 1.5,
  text.labels = c(
    expression('P'['ER+ High vs IC10']*'=4.4x10'^'-5'),
    expression('P'['HER2+ vs IC10']*'=3.4x10'^'-2')
    ),
  text.x = c(1.1,1.05),
  text.y = c(1.8,1.65),
  text.cex = c(1.2,1.2),
  xaxis.cex = 1.5,
  yaxis.cex = 1.5,
  ylimits = c(0,2),
  yat = seq(0,3,1),
  resolution = 300
  )
create.multipanelplot(
  list(primary, metastatic),
  filename = 'figure3b.pdf',
  resolution = 300
  )

### SAVE ##########################################################################################
# write source data to file 
write.table(
  pplot_data[,c('Sample','d_pre_nokat_clock','group')],
  file = file.path(main_repo_path, 'data', 'Figure3b_sourcetable_primary.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )
write.table(
  mplot_data[,c('Sample','d_pre_nokat_clock','group')],
  file = file.path(main_repo_path, 'data', 'Figure3b_sourcetable_metastatic.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )