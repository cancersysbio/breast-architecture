### CREATE FIGURE 1f AND EXTENDED DATA FIGURE 1a #############################################################################
# creates figure 1f and extended data figure 1a
# figure 1f provides the proportion of the IC subgroups and IC subtypes across stages in all breast cancer patients 
# extended data figure 1a provides the proportion of the IC subgroups and IC subtypes across stages in ER+ breast cancer patients only

### PREAMBLE #####################################################################################
library(yaml)
library(ggplot2)
library(ggpubr)
require(patchwork)
library(grid)
library(gridExtra) 
library(ggnewscale)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

col.ic = c("ic1" = "#ff5500ff", "ic2" = "#00ee77ff", 
           "ic4ER+" = "#00c4cdff", "ic4ER-" = "#c2b7dbff",
           "ic5" = "#8b0100ff", "ic6" = "#fffe41ff", "ic3/ic7" = "#cd32caff",
           "ic8" = "#feaa00ff","ic9" = "#ee82edff", "ic10" = "#7c26ccff")

col.group = c("ER+ High" = "#fa954eff","ER+ Typical" = "#a1c1d4ff","HER2+" = "#8b0100ff","TNBC" = "#6300a1", "IC10" = "#7c26ccff", "IC4ER-" = "#c2b7dbff")

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'
outdir <- '/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/publication/figures/'

### Load ENiClust predictions----
## DCIS
# Using the iC10-RNA-only
labels_dcis = read.table(paste0(basedir,"khoulaha/BreastLandscape/clinical/DCIS_clinical.csv"), header = T, sep=",")
colnames(labels_dcis)[which(colnames(labels_dcis) == "IC10_RNA")] = "ic10_class"
labels_dcis = labels_dcis[which(labels_dcis$Clustering_RAHBT  == 1 | labels_dcis$Clustering_TBCRC == 1),] #stick true to the original paper
any(duplicated(labels_dcis$Patient_ID))

labels_dcis = labels_dcis[which(!is.na(labels_dcis$ic10_class)),]
labels_dcis$ic10_class = paste0("ic", labels_dcis$ic10_class)
labels_dcis$ic10_class[which(labels_dcis$ic10_class %in% c("ic3", "ic7"))] = "ic3/ic7"
labels_dcis$ic10_class[which(labels_dcis$ic10_class == "ic4" & labels_dcis$ER_RNA == "+")] = "ic4ER+"
labels_dcis$ic10_class[which(labels_dcis$ic10_class == "ic4" & labels_dcis$ER_RNA == "-")] = "ic4ER-"
table(labels_dcis$ic10_class)

# Making groups                           
labels_dcis$group = NA
labels_dcis$group[which(labels_dcis$ic10_class == "ic5")] = "HER2+"
labels_dcis$group[which(labels_dcis$ic10_class %in% c("ic1", "ic2", "ic6", "ic9"))] = "ER+ High"
labels_dcis$group[which(labels_dcis$ic10_class %in% c("ic3/ic7","ic4ER+","ic8"))] = "ER+ Typical"
labels_dcis$group[which(labels_dcis$ic10_class == "ic10")] = "IC10"
labels_dcis$group[which(labels_dcis$ic10_class == "ic4ER-")] = "IC4ER-"
table(labels_dcis$group, useNA = "always")

## Primary/Metastasis
# Load megatable
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-04-24_primary_megatable.txt"), header = T, sep = "\t")
mega_1$Sample_Type = "Primary"
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-04-24_metastatic_megatable.txt"), header = T, sep = "\t")
mega_2$Sample_Type = "Metastatic"
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])
table(mega$group, useNA = "always") # 2 IC4 primary samples with missing hormone receptor status + 7 IC4 metastatic samples with missing hormone receptor status

# Load predictions for sanity checks
labels_icgc = read.table(paste0(basedir,"lisem/bc_landscape/ic_subtypes/data/eniclust/icgc_predictions.txt"), header = T, sep="\t")
all(labels_icgc$Sample %in% mega$Sample)

hmf_predictions <- read.table(paste0(basedir,"lisem/bc_landscape/github/hmf_predictions.txt"), sep="\t", header = T)
all(hmf_predictions$Sample %in% mega$Sample)

## Merging
labels_prim = mega[which(mega$Sample_Type == "Primary" & !is.na(mega$ENiClust)),]
labels_mets = mega[which(mega$Sample_Type == "Metastatic" & !is.na(mega$ENiClust)),]

labels_prim$ENiClust[which(labels_prim$ENiClust == "Other")] = "ic3/ic7"
labels_mets$ENiClust[which(labels_mets$ENiClust == "Other")] = "ic3/ic7"

labels_prim$ENiClust[which(labels_prim$ENiClust == "ic4" & labels_prim$Sample %in% mega$Sample[which(mega$group == "ER+ Typical")])] = "ic4ER+"
labels_prim$ENiClust[which(labels_prim$ENiClust == "ic4" & labels_prim$Sample %in% mega$Sample[which(mega$group == "IC4ER-")])] = "ic4ER-"

labels_mets$ENiClust[which(labels_mets$ENiClust == "ic4" & labels_mets$Sample %in% mega$Sample[which(mega$group == "ER+ Typical")])] = "ic4ER+"
labels_mets$ENiClust[which(labels_mets$ENiClust == "ic4" & labels_mets$Sample %in% mega$Sample[which(mega$group == "IC4ER-")])] = "ic4ER-"

table(labels_prim$ENiClust)
table(labels_mets$ENiClust)
table(labels_prim$group, useNA = "always")
table(labels_mets$group, useNA = "always")

## Plot
df = data.frame(cohort = c(rep("DCIS", nrow(labels_dcis)),
                           rep("Primary", nrow(labels_prim)),
                           rep("Metastatic", nrow(labels_mets))),
                group = c(labels_dcis$group, 
                          labels_prim$group,
                          labels_mets$group),
                ic10 = c(labels_dcis$ic10_class, 
                         labels_prim$ENiClust,
                         labels_mets$ENiClust), stringsAsFactors = F)

# All
df_calculated <- df[which(!is.na(df$group)),] %>% 
  dplyr::count(cohort, group) %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$cohort = factor(df_calculated$cohort, levels = c("DCIS", "Primary", "Metastatic"))
df_calculated$group = factor(df_calculated$group, levels = c("ER+ High","ER+ Typical","HER2+","IC10","IC4ER-"))

panelA = ggplot(df_calculated, aes(x = cohort, y = Percent, fill = group))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.group) +
  ylab("Proportion of samples") +
  xlab("") + scale_x_discrete(labels = c("DCIS" = "DCIS\n(RNA)")) +
  geom_label(data = df_calculated %>% dplyr::group_by(cohort) %>% dplyr::summarise(Sum = sum(n)), aes(label = paste0("n=",Sum), y = 1.05, x = cohort),inherit.aes=FALSE,label.size=NA,alpha=0,size=5) +
  theme_LM + ggtitle("Discovery") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14),
        axis.text.y = element_text(size=14), 
        axis.title = element_text(size=20), 
        title = element_text(size=20), 
        legend.text = element_text(size=14))

df_calculated <- df[which(df$ic10 != "ic4"),] %>% dplyr::count(cohort, ic10) %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$cohort = factor(df_calculated$cohort, levels = c("DCIS", "Primary", "Metastatic"))
df_calculated$ic10 = factor(df_calculated$ic10, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8","ic5","ic10","ic4ER-"))

panelB = ggplot(df_calculated, aes(x = cohort, y = Percent, fill = ic10))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.ic[which(names(col.ic) %in% df_calculated$ic10)],
                    labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-")) +
  ylab("") +
  xlab("") + scale_x_discrete(labels = c("DCIS" = "DCIS\n(RNA)")) +
  geom_label(data = df_calculated %>% dplyr::group_by(cohort) %>% dplyr::summarise(Sum = sum(n)), aes(label = paste0("n=", Sum), y = 1.05, x = cohort),inherit.aes = F,label.size=NA,alpha=0,size=5) +
  theme_LM + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size=20), 
        legend.text = element_text(size=14))

mylegend = cowplot::get_legend(ggplot(labels_prim[which(labels_prim$ENiClust != "ic4"),],aes(x=CN1,y=CN2)) +
                                 geom_col(aes(fill=group)) +
                                 scale_fill_manual(name = "IC subgroups", values = col.group) +
                                 new_scale_fill() +
                                 geom_col(aes(fill=factor(ENiClust, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8","ic5","ic10","ic4ER-")))) +
                                 scale_fill_manual(name = "IC subtypes", values = col.ic[which(names(col.ic) %in% df$ic10)],
                                                   labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-"))+
                                 theme(legend.text = element_text(size=14), legend.title = element_text(size=16)) )

### Look at ER+ BC patients----
# ER+ only
df_calculated <- df[which(!is.na(df$group) & grepl("ER[+]",df$group)),] %>% dplyr::count(cohort, group) %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$cohort = factor(df_calculated$cohort, levels = c("DCIS", "Primary", "Metastatic"))
df_calculated$group = factor(df_calculated$group, levels = c("ER+ High","ER+ Typical"))

panelAsub = ggplot(df_calculated, aes(x = cohort, y = Percent, fill = group))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.group) +
  ylab("Proportion of ER+ samples") +
  xlab("") + scale_x_discrete(labels = c("DCIS" = "DCIS\n(RNA)")) +
  geom_label(data = df_calculated %>% dplyr::group_by(cohort) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = cohort), inherit.aes = F,label.size=NA,alpha=0,size=5) +
  theme_LM +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14), 
        axis.text.y = element_text(size=14), 
        axis.title = element_text(size=20), 
        legend.text = element_text(size=14))

df_calculated <- df[which(df$ic10 %in% c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8")),] %>% 
  dplyr::count(cohort, ic10) %>%
  dplyr::group_by(cohort) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$cohort = factor(df_calculated$cohort, levels = c("DCIS", "Primary", "Metastatic"))
df_calculated$ic10 = factor(df_calculated$ic10, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8"))

panelBsub = ggplot(df_calculated, aes(x = cohort, y = Percent, fill = ic10))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.ic[which(names(col.ic) %in% df_calculated$ic10)],
                    labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-")) +
  ylab("") +
  xlab("") + scale_x_discrete(labels = c("DCIS" = "DCIS\n(RNA)")) +
  geom_label(data = df_calculated %>% dplyr::group_by(cohort) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = cohort), inherit.aes = F, label.size=NA, alpha=0, size=5) +
  theme_LM +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_text(size=20), 
        legend.text = element_text(size=14))

mylegend_ER = cowplot::get_legend(ggplot(labels_prim[which(labels_prim$group %in% c("ER+ Typical","ER+ High")),],aes(x=CN1,y=CN2)) +
                                    geom_col(aes(fill=group)) +
                                    scale_fill_manual(name = "IC subgroups", values = c("ER+ High" = "#ff6a00","ER+ Typical" = "#7aa6c2","HER2+" = "#a15300","IC10" = "#6300a1","IC4ER-" = "#c2b7dbff")) +
                                    new_scale_fill() +
                                    geom_col(aes(fill=factor(ENiClust, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8","ic5","ic10","ic4ER-")))) +
                                    scale_fill_manual(name = "IC subtypes", values = col.ic[which(names(col.ic) %in% df$ic10)],
                                                      labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-"))+
                                    theme(legend.text = element_text(size=14), legend.title = element_text(size=16)) )

### SAVE ##########################################################################################
## Figures----
png(filename = file.path(outdir, 'Figure1f.png'), res = 300, width = 8.59, height = 5.40, units = 'in')
(panelA + rremove("legend")) + (panelB + rremove("legend")) + mylegend
dev.off()

png(filename = file.path(outdir, 'ExtendedData1a.png'), res = 300, width = (11.5/1.4), height = (7/1.4), units = 'in')
(panelAsub + rremove("legend")) + (panelBsub + rremove("legend")) + mylegend_ER + plot_layout(widths = c(0.4,0.4,0.2))
dev.off()

## SourceData----
sourcetable = df[,c(1,3,2)]
colnames(sourcetable) = c("Stage","Subtype","Subgroup")
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/submission/data/Figure1f_ExtendedData1a_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

## Enrichement tests----
# Enrichment for ER+ High in Metastatic? Yes.
x = length(which(!is.na(df$group) & grepl("ER[+]",df$group) & df$cohort == "Metastatic" & df$group == "ER+ High"))
y = length(which(!is.na(df$group) & grepl("ER[+]",df$group) & df$cohort == "Metastatic" & df$group != "ER+ High"))
z = length(which(!is.na(df$group) & grepl("ER[+]",df$group) & df$cohort == "Primary" & df$group == "ER+ High"))
w = length(which(!is.na(df$group) & grepl("ER[+]",df$group) & df$cohort == "Primary" & df$group != "ER+ High"))

mat = matrix(c(x,y,z,w), nrow = 2)
fisher.test(mat) # 0.019

# Enrichment for ER+ Typical ICs in DCIS? No.
x = length(which(!is.na(df$group) & df$cohort == "DCIS" & df$ic10 == "ic3/ic7"))
y = length(which(!is.na(df$group) & df$cohort == "DCIS" & df$ic10 != "ic3/ic7"))
z = length(which(!is.na(df$group) & df$cohort != "DCIS" & df$ic10 == "ic3/ic7"))
w = length(which(!is.na(df$group) & df$cohort != "DCIS" & df$ic10 != "ic3/ic7"))

mat = matrix(c(x,y,z,w), nrow = 2)
fisher.test(mat) # 0.42

# Enrichment for HER2+ in DCIS? Yes.
x = length(which(!is.na(df$group) & df$cohort == "DCIS" & df$ic10 == "ic5"))
y = length(which(!is.na(df$group) & df$cohort == "DCIS" & df$ic10 != "ic5"))
z = length(which(!is.na(df$group) & df$cohort != "DCIS" & df$ic10 == "ic5"))
w = length(which(!is.na(df$group) & df$cohort != "DCIS" & df$ic10 != "ic5"))

mat = matrix(c(x,y,z,w), nrow = 2)
fisher.test(mat) # 2.98x10-6




