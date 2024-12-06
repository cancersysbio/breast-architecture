### CREATE SUPPLEMENTARY FIGURE 1g #############################################################################
# creates supplementary figure 1g
# supplementary figure 1g provides the proportion and count of IC subtype (left) and subgroup (right) tumors in sensitive, intermediate and resistant cases from the clinical trial (NCT00651976) as described in Giltnane et al. Science Translational Medicine 2017.

### PREAMBLE #####################################################################################
library(yaml)
library(ggplot2)
library(ggpubr)
require(patchwork)
library(dplyr)
library(openxlsx)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

col.group = c("ER+ High" = "#fa954eff","ER+ Typical" = "#a1c1d4ff","HER2+" = "#8b0100ff","TNBC" = "#6300a1",
              "IC10" = "#7c26ccff", "IC4ER-" = "#c2b7dbff")

col.ic = c("ic1" = "#ff5500ff", "ic2" = "#00ee77ff", "ic3" = "#cd3279ff",
           "ic4" = "#bbe0e2ff", "ic4ER+" = "#00c4cdff", "ic4ER-" = "#c2b7dbff",
           "ic5" = "#8b0100ff", "ic6" = "#fffe41ff", "ic7" = "#0000cdff", "ic3/ic7" = "#cd32caff",
           "ic8" = "#feaa00ff","ic9" = "#ee82edff", "ic10" = "#7c26ccff",
           "Other" = "grey", "concordant" = "#d0d7d9")

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'
outdir <- paste0(basedir,"lisem/bc_landscape/revision/figures/")

## Load Clinical data of Giltnane with iC10 RNA-only calls
data <- read.table(paste0(basedir,"lisem/bc_landscape/github/giltnane_clinical.txt"), header = TRUE)

## Statistical tests
data$pCR = data[,"Response_Category.Updated"]
data$histo = data[,"Histotype_NST_Lobular_Other"]

data$ic_group = NA
data$ic_group[which(gsub("ic","",data$iC10) %in% c(1,2,6,9))] = "High"
data$ic_group[which(gsub("ic","",data$iC10) %in% c(3,4,7,8))] = "Typical"
data$ic_group[which(gsub("ic","",data$iC10) == 5)] = "HER2+"
data$ic_group[which(gsub("ic","",data$iC10) == 10)] = "TNBC"

## Preparing some plots----
panelA <- data %>%
  dplyr::filter(!is.na(pCR)) %>%
  dplyr::mutate(iC10 = ifelse(iC10 == "ic4", "ic4ER+", iC10)) %>%
  dplyr::group_by(pCR, iC10) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  ggplot(aes(x=factor(pCR, levels = c("Sensitive","Intermediate","Resistant")), y=freq, fill = factor(iC10,levels=paste0("ic",c(1,2,6,9,3,"4ER+",7,8,5,10))))) +
  geom_col(position="fill") +
  scale_fill_manual(name = "", values = col.ic) +
  theme_LM + xlab("Endocrine therapy response") + ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12))

panelB <- data %>%
  dplyr::filter(!is.na(pCR)) %>%
  dplyr::mutate(iC10 = ifelse(iC10 == "ic4", "ic4ER+", iC10)) %>%
  ggplot(aes(x=factor(pCR, levels = c("Sensitive","Intermediate","Resistant")), fill = factor(iC10,levels=paste0("ic",c(1,2,6,9,3,"4ER+",7,8,5,10))))) +
  geom_bar() +
  scale_fill_manual(name = "", values = col.ic) +
  theme_LM + xlab("Endocrine therapy response") + ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12))

panelC <- data %>%
  dplyr::filter(!is.na(pCR)) %>%
  dplyr::group_by(pCR, ic_group) %>% dplyr::mutate(ic_group = ifelse(ic_group == "TNBC",ic_group,paste("ER+",ic_group))) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  ggplot(aes(x=factor(pCR, levels = c("Sensitive","Intermediate","Resistant")), y=freq, fill = ic_group)) +
  geom_col(position="fill") +
  scale_fill_manual(name = "", values = col.group) +
  theme_LM + xlab("Endocrine therapy response") + ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12))

panelD <- data %>%
  dplyr::filter(!is.na(pCR)) %>% dplyr::mutate(ic_group = ifelse(ic_group == "TNBC",ic_group,paste("ER+",ic_group))) %>%
  ggplot(aes(x=factor(pCR, levels = c("Sensitive","Intermediate","Resistant")), fill = ic_group)) +
  geom_bar() +
  scale_fill_manual(name = "", values = col.group) +
  theme_LM + xlab("Endocrine therapy response") + ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=12),
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), 
        legend.text = element_text(size=12))

### SAVE ##########################################################################################
## Figures----
png(filename = file.path(outdir, 'giltane_barplots.png'), res = 300, width = 8.88, height = 7.27, units = 'in')
(panelA + panelB) / (panelC + panelD)
dev.off()

## SourceData----
sourcetable = data[,c(1:4,7)]
colnames(sourcetable)[4:5] = c("Subtype","Subgroup")
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/github/SupplementaryFigure1g_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

## Enrichement tests----
# Without Intermediate: High vs. Typical
x = length(which(data$pCR == "Resistant" & data$ic_group == "High"))
y = length(which(data$pCR == "Sensitive" & data$ic_group == "High"))
z = length(which(data$pCR == "Resistant" & data$ic_group == "Typical"))
w = length(which(data$pCR == "Sensitive" & data$ic_group == "Typical"))
mat = matrix(c(x,y,z,w), byrow = T, nrow = 2)
p1 <- fisher.test(mat)$p.value
e1 <- fisher.test(mat)$estimate

# Without Intermediate: High vs. Rest
x = length(which(data$pCR == "Resistant" & data$ic_group == "High"))
y = length(which(data$pCR == "Sensitive" & data$ic_group == "High"))
z = length(which(data$pCR == "Resistant" & data$ic_group != "High"))
w = length(which(data$pCR == "Sensitive" & data$ic_group != "High"))
mat = matrix(c(x,y,z,w), byrow = T, nrow = 2)
p2 <- fisher.test(mat)$p.value
e2 <- fisher.test(mat)$estimate

# With Resistant-Intermediate: High vs. Typical
x = length(which(data$pCR != "Sensitive" & data$ic_group == "High"))
y = length(which(data$pCR == "Sensitive" & data$ic_group == "High"))
z = length(which(data$pCR != "Sensitive" & data$ic_group == "Typical"))
w = length(which(data$pCR == "Sensitive" & data$ic_group == "Typical"))
mat = matrix(c(x,y,z,w), byrow = T, nrow = 2)
p3 <- fisher.test(mat)$p.value
e3 <- fisher.test(mat)$estimate

# With Resistant-Intermediate: High vs. Typical
x = length(which(data$pCR != "Sensitive" & data$ic_group == "High"))
y = length(which(data$pCR == "Sensitive" & data$ic_group == "High"))
z = length(which(data$pCR != "Sensitive" & data$ic_group != "High"))
w = length(which(data$pCR == "Sensitive" & data$ic_group != "High"))
mat = matrix(c(x,y,z,w), byrow = T, nrow = 2)
p4 <- fisher.test(mat)$p.value
e4 <- fisher.test(mat)$estimate

# With Sensitive-Intermediate: High vs. Typical
x = length(which(data$pCR == "Resistant" & data$ic_group == "High"))
y = length(which(data$pCR != "Resistant" & data$ic_group == "High"))
z = length(which(data$pCR == "Resistant" & data$ic_group == "Typical"))
w = length(which(data$pCR != "Resistant" & data$ic_group == "Typical"))
mat = matrix(c(x,y,z,w), byrow = T, nrow = 2)
p5 <- fisher.test(mat)$p.value
e5 <- fisher.test(mat)$estimate

# With Sensitive-Intermediate: High vs. Typical
x = length(which(data$pCR == "Resistant" & data$ic_group == "High"))
y = length(which(data$pCR != "Resistant" & data$ic_group == "High"))
z = length(which(data$pCR == "Resistant" & data$ic_group != "High"))
w = length(which(data$pCR != "Resistant" & data$ic_group != "High"))
mat = matrix(c(x,y,z,w), byrow = T, nrow = 2)
p6 <- fisher.test(mat)$p.value
e6 <- fisher.test(mat)$estimate

# Summary
max(c(p1,p2,p3,p4,p5,p6))
min(c(e1,e2,e3,e4,e5,e6))
