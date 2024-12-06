### CREATE EXTENDED DATA FIGURE 1f and 1g #############################################################################
# creates extended data figure 1f and 1g
# extended data figure 1f provides the subtyping of each samples within each matched pairs in the Hartwig dataset
# extended data figure 1g provides the concordance from the IC subtypes/subgroups to the PAM50 subtypes at each stage

### PREAMBLE #####################################################################################
library(yaml)
library(ggplot2)
library(ggpubr)
require(ggraph)
library(ggsankey, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.2.2-Seurat/")
library(tidygraph)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

col.ic = c("ic1" = "#ff5500ff", "ic2" = "#00ee77ff", "ic3" = "#cd3279ff", 
           "ic4ER+" = "#00c4cdff", "ic4ER-" = "#c2b7dbff",
           "ic5" = "#8b0100ff", "ic6" = "#fffe41ff", "ic7" = "#0000cdff", "ic3/ic7" = "#cd32caff",
           "ic8" = "#feaa00ff","ic9" = "#ee82edff", "ic10" = "#7c26ccff")
names(col.ic) <- toupper(names(col.ic))

col.group = c("ER+ High" = "#fa954eff","ER+ Typical" = "#a1c1d4ff","HER2+" = "#8b0100ff","TNBC" = "#6300a1")

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'
outdir <- '/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/publication/figures/'

### ITH Hartwig----
## Load metadata with donor-duplicates
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-03-01_primary_megatable_with_duplicates.txt"), header = T, sep = "\t")
mega_1$Sample_Type = "Primary"
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-03-01_metastatic_megatable_with_duplicates.txt"), header = T, sep = "\t")
mega_2$Sample_Type = "Metastatic"
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])
table(mega$group, useNA = "always")

mega$ENiClust[which(mega$ENiClust == "Other")] = "IC3/IC7"
mega$ENiClust[which(mega$ENiClust == "ic4" & mega$ER == 1)] = "IC4ER+"
mega$ENiClust[which(mega$ENiClust == "ic4" & mega$ER == 0)] = "IC4ER-"
mega$ENiClust[which(startsWith(mega$ENiClust,"ic"))] = gsub('ic',"IC",mega$ENiClust[which(startsWith(mega$ENiClust,"ic"))])
table(mega$ENiClust, useNA = "always")

mega$group[which(mega$group %in% c("IC10", "IC4ER-"))] = "TNBC"
table(mega$group, useNA = "always")
table(mega$ENiClust[which(is.na(mega$group))], useNA = "always")

hmf_clinic = read.table(paste0(basedir,"lisem/bc_landscape/github/Hartwig_clinical.txt"), header = T, sep = "\t")
all(hmf_clinic$Individual.System.ID %in% mega$Sample)

mega$biopsy <- sapply(mega$Sample, function(x) if(x %in% hmf_clinic$Individual.System.ID){hmf_clinic$biopsySite_LM[which(hmf_clinic$Individual.System.ID == x)]}else{NA})

## Load pam50 calls: see apply_pam50_hartwig.R
pam50_dup <- read.table(paste0(basedir,'lisem/bc_landscape/rna/2022-09-14_hartwig_vst_pam50.txt'), header = T)

# Add PAM50 to metadata
mega$pam50 <- sapply(mega$Sample, function(x) if(substr(x,1,nchar("CURTIS_H004175_T01_01")) %in% substr(pam50_dup$sample,1,nchar("CURTIS_H004175_T01_01"))){pam50_dup$pam50[which(substr(pam50_dup$sample,1,nchar("CURTIS_H004175_T01_01")) == substr(x,1,nchar("CURTIS_H004175_T01_01")))]}else{NA})
table(mega$pam50)

## Determine ITH samples
table(mega$Sample[which(!is.na(mega$pam50))] %in% mega$Sample[which(mega$cohort == "Hartwig" & substr(mega$Sample,1,nchar("CURTIS_H004175_T")) %in% substr(mega$Sample[which(duplicated(substr(mega$Sample,1,nchar("CURTIS_H004175_T"))))],1,nchar("CURTIS_H004175_T")))])
ith_samp <- mega$Sample[which(mega$cohort == "Hartwig" & substr(mega$Sample,1,nchar("CURTIS_H004175_T")) %in% substr(mega$Sample[which(duplicated(substr(mega$Sample,1,nchar("CURTIS_H004175_T"))))],1,nchar("CURTIS_H004175_T")))]
ith_samp <- unique(substr(ith_samp, 1, nchar("CURTIS_H004343")))

# changes to expect: CURTIS_H004082, CURTIS_H004092, CURTIS_H004159, CURTIS_H004188, CURTIS_H004423, CURTIS_H004627
sourcetable_ed_1f = mega %>% dplyr::filter(substr(Sample, 1, nchar("CURTIS_H004343")) %in% ith_samp) %>% dplyr::select(Sample,biopsy,ENiClust,group,pam50) %>% dplyr::arrange(Sample)
colnames(sourcetable_ed_1f) = c("Sample","Stage/Site","Subtype","Subgroup","PAM50")

## Network
mydata = data.frame()
for(x in ith_samp){
  i = length(which(grepl(x, mega$Sample)))
  mydata_tmp = data.frame(from = mega$Sample[which(grepl(x, mega$Sample))][1:i],
                          to = mega$Sample[which(grepl(x, mega$Sample))][c(2:i,1)],
                          stringsAsFactors = F)
  mydata = rbind(mydata, mydata_tmp)
} 
mydata_ic <- mydata[which(mydata$from %in% mega$Sample[which(!is.na(mega$ENiClust))] & mydata$to %in% mega$Sample[which(!is.na(mega$ENiClust))]),]
mydata_ic$ic = sapply(1:nrow(mydata_ic), function(i) if(as.character(mega$ENiClust[which(mega$Sample == mydata_ic$to[i])]) == as.character(mega$ENiClust[which(mega$Sample == mydata_ic$from[i])])){"Same"}else{"Different"})
mydata_group <- mydata[which(mydata$from %in% mega$Sample[which(!is.na(mega$group))] & mydata$to %in% mega$Sample[which(!is.na(mega$group))]),]
mydata_group$group = sapply(1:nrow(mydata_group), function(i) if(as.character(mega$group[which(mega$Sample == mydata_ic$to[i])]) == as.character(mega$group[which(mega$Sample == mydata_ic$from[i])])){"Same"}else{"Different"})

md_graph <- as_tbl_graph(mydata_group, directed = FALSE)

md_graph <- md_graph %>% activate(nodes) %>% 
  dplyr::mutate(group = sapply(name, function(x) as.character(mega$group[which(mega$Sample == x)])))

md_graph <- md_graph %>% activate(nodes) %>% 
  dplyr::mutate(site = sapply(name, function(x) mega$biopsy[which(mega$Sample == x)]))

md_graph <- md_graph %>% activate(nodes) %>% 
  dplyr::mutate(type = sapply(name, function(x) mega$pam50[which(mega$Sample == x)]))

ed_1f = ggraph(md_graph, layout = "fr") +
  geom_node_voronoi(aes(fill = type), max.radius = 0.5, colour = 'white', alpha = 0.7) +
  geom_edge_link(aes(colour = factor(group))) + 
  geom_node_point(aes(colour = group), size = 4) +
  geom_node_label(aes(label = site), repel = TRUE) +
  scale_fill_manual(name = "PAM50", values = c("Basal" = "#7a5195", "Her2" = "#ef5675", "LumA" = "#003f5c", "LumB" = "#ffa600")) +
  scale_colour_manual(name = "IC subgroups", values = col.group) +
  scale_edge_colour_manual(name = "", values = c("Different" = "red", "Same" = "black")) +
  theme_bw() + ggtitle(paste0("Pairs in Hartwig (n=",length(ith_samp),")")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(), title = element_text(size=24), legend.title = element_text(size=16), legend.text = element_text(size=14))

### Sankey----
# Load metadata without donor-duplicates
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")
mega_1$Sample_Type = "Primary"
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")
mega_2$Sample_Type = "Metastatic"
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])

mega$ENiClust[which(mega$ENiClust == "Other")] = "IC3/IC7"
mega$ENiClust[which(mega$ENiClust == "ic4" & mega$ER == 1)] = "IC4ER+"
mega$ENiClust[which(mega$ENiClust == "ic4" & mega$ER == 0)] = "IC4ER-"
mega$ENiClust[which(startsWith(mega$ENiClust,"ic"))] = gsub('ic',"IC",mega$ENiClust[which(startsWith(mega$ENiClust,"ic"))])
mega$group[which(mega$group %in% c("IC10", "IC4ER-"))] = "TNBC"

hmf_clinic = read.table(paste0(basedir,"khoulaha/BreastLandscape/clinical/Hartwig_clinical.txt"), header = T, sep = "\t")
mega$biopsy <- sapply(mega$Sample, function(x) if(x %in% hmf_clinic$Individual.System.ID){hmf_clinic$biopsySite_LM[which(hmf_clinic$Individual.System.ID == x)]}else{NA})

# Load pam50 calls: see apply_pam50_hartwig.R
pam50_dup <- read.table(paste0(basedir,'lisem/bc_landscape/rna/2022-09-14_hartwig_vst_pam50.txt'), header = T)

# Add PAM50 to metadata
table(substr(mega$Sample,1,nchar("CURTIS_H004175_T01_01")) %in% substr(pam50_dup$sample,1,nchar("CURTIS_H004175_T01_01"))) # Hartwig samples with RNA
mega$pam50 <- sapply(mega$Sample, function(x) if(substr(x,1,nchar("CURTIS_H004175_T01_01")) %in% substr(pam50_dup$sample,1,nchar("CURTIS_H004175_T01_01"))){pam50_dup$pam50[which(substr(pam50_dup$sample,1,nchar("CURTIS_H004175_T01_01")) == substr(x,1,nchar("CURTIS_H004175_T01_01")))]}else{NA})
table(mega$pam50)

df_mets <- mega %>% dplyr::filter(mega$Sample %in% mega_2$Sample & !is.na(ENiClust) & !is.na(group) & !is.na(pam50))
length(which(df_mets$group == "ER+ High" & df_mets$pam50 == "LumB"))/length(which(df_mets$group == "ER+ High"))
length(which(df_mets$group == "ER+ Typical" & df_mets$pam50 == "LumB"))/length(which(df_mets$group == "ER+ Typical"))

# in Metastatic (Hartwig+2PCAWG)
ed_1g_right <- mega %>% 
  dplyr::filter(mega$Sample %in% mega_2$Sample & !is.na(ENiClust) & !is.na(group) & !is.na(pam50)) %>%
  make_long(ENiClust, group, pam50) %>%
  dplyr::mutate(node = factor(node, levels = rev(c(paste0("IC",c(1,2,6,9)),"IC3/IC7",paste0("IC",c("4ER+",8,5,10,"4ER-")),
                                                   "ER+ High","ER+ Typical","HER2+","TNBC",
                                                   "LumB","LumA","Normal","Her2","Basal")))) %>%
  dplyr::mutate(next_node = factor(next_node, levels = rev(c(paste0("IC",c(1,2,6,9)),"IC3/IC7",paste0("IC",c("4ER+",8,5,10,"4ER-")),
                                                             "ER+ High","ER+ Typical","HER2+","TNBC",
                                                             "LumB","LumA","Normal","Her2","Basal")))) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = factor(node),
             label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = 1) +
  geom_sankey_label(size = 2) +
  #geom_sankey_label(aes(color = ifelse(node %in% c("IC5","IC7","IC10","HER2+","TNBC","LumA","Basal"), "1", "0"))) + # not working
  #scale_color_manual(values = c("1"="white", "0"="black")) + 
  scale_x_discrete(labels = c("ENiClust" = "ENiClust IC", "group" = "IC subgroup", "pam50" = "PAM50")) +
  scale_fill_manual(name = "", values = c(col.ic[which(names(col.ic) %in% mega$ENiClust)],
                                          col.group, "Basal" = "#7a5195", "Her2" = "#ef5675", "LumA" = "#003f5c", "LumB" = "#ffa600","Normal" = "#ef94d5")) +
  xlab('') +
  theme_sankey(base_size = 5) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 7), 
        axis.text.x = element_text(size=6)) + 
  labs(title = paste0("Metastatic (n=",nrow(mega %>% filter(mega$Sample %in% mega_2$Sample & !is.na(ENiClust) & !is.na(group) & !is.na(pam50))),")"))

# in Primary (use TCGA-WES)
tcga_ic10 = read.table(paste0(basedir,"lisem/bc_landscape/ic_subtypes/data/eniclust/tcga_predictions.txt"), header = T, sep = "")
any(duplicated(tcga_ic10$Sample))
er_status <- read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-02-09_rna_megatable_primary_met_eniclust_v2_er_ihc.txt"),sep="\t",header=T) %>%
  dplyr::select(Sample, ER) %>%
  dplyr::filter(!is.na(ER))
er_status <- read.table(paste0(basedir,"lisem/bc_landscape/ic_subtypes/data/eniclust/tcga_transformed_matrix.txt"),sep="\t",header=T) %>%
  dplyr::filter(!Sample %in% er_status$Sample) %>% dplyr::select(Sample, ER) %>% 
  dplyr::rowwise() %>% dplyr::mutate(ER = ER/4) %>%
  dplyr::bind_rows(er_status)

table(tcga_ic10$Sample %in% er_status$Sample)

tcga_ic10$voting[which(tcga_ic10$voting == "ic4" & tcga_ic10$Sample %in% er_status$Sample[which(er_status$ER == 1)])] = "ic4ER+"
tcga_ic10$voting[which(tcga_ic10$voting == "ic4" & tcga_ic10$Sample %in% er_status$Sample[which(er_status$ER == 0)])] = "ic4ER-"
tcga_ic10$voting[which(tcga_ic10$voting == "Other")] = "ic3/ic7"
table(tcga_ic10$voting, useNA = "always")

tcga_ic10$group = NA
tcga_ic10$group[which(tcga_ic10$voting == "ic5")] = "HER2+"
tcga_ic10$group[which(tcga_ic10$voting %in% c("ic1", "ic2", "ic6", "ic9"))] = "ER+ High"
tcga_ic10$group[which(tcga_ic10$voting %in% c("ic3/ic7", "ic4ER+", "ic8"))] = "ER+ Typical"
tcga_ic10$group[which(tcga_ic10$voting %in% c("ic10", "ic4ER-"))] = "TNBC"
table(tcga_ic10$group, useNA = "always")

tcga_pam50 = read.table(paste0(basedir,"lisem/TCGA/tcga_vst_pam50.txt"), header = T, sep = "")
table(substr(tcga_pam50$sample,1,nchar('TCGA-PE-A5DC')) %in% tcga_ic10$Sample)
tcga_pam50$ENiClust <- toupper(sapply(substr(tcga_pam50$sample,1,nchar('TCGA-PE-A5DC')), function(x) if(x %in% tcga_ic10$Sample){tcga_ic10$voting[which(tcga_ic10$Sample == x)]}else{NA}))
tcga_pam50$group <- sapply(substr(tcga_pam50$sample,1,nchar('TCGA-PE-A5DC')), function(x) if(x %in% tcga_ic10$Sample){tcga_ic10$group[which(tcga_ic10$Sample == x)]}else{NA})

ed_1g_middle <- tcga_pam50 %>% 
  dplyr::filter(!is.na(ENiClust) & !is.na(group) & !is.na(pam50)) %>%
  make_long(ENiClust, group, pam50) %>%
  dplyr::mutate(node = factor(node, levels = rev(c(paste0("IC",c(1,2,6,9)),"IC3/IC7",paste0("IC",c("4ER+",8,5,10,"4ER-")),
                                                   "ER+ High","ER+ Typical","HER2+","TNBC",
                                                   "LumB","LumA","Normal","Her2","Basal")))) %>%
  dplyr::mutate(next_node = factor(next_node, levels = rev(c(paste0("IC",c(1,2,6,9)),"IC3/IC7",paste0("IC",c("4ER+",8,5,10,"4ER-")),
                                                             "ER+ High","ER+ Typical","HER2+","TNBC",
                                                             "LumB","LumA","Normal","Her2","Basal")))) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = factor(node),
             label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = 1) +
  geom_sankey_label(size = 2) +
  scale_x_discrete(labels = c("ENiClust" = "ENiClust IC", "group" = "IC subgroup", "pam50" = "PAM50")) +
  scale_fill_manual(name = "", values = c(col.ic[which(names(col.ic) %in% tcga_pam50$ENiClust)],
                                          col.group, "Basal" = "#7a5195", "Her2" = "#ef5675", "LumA" = "#003f5c", "LumB" = "#ffa600","Normal" = "#ef94d5")) +
  xlab('') +
  theme_sankey(base_size = 5) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 7), 
        axis.text.x = element_text(size=6)) + 
  labs(title = paste0("Primary (n=",nrow(tcga_pam50 %>% dplyr::filter(!is.na(ENiClust) & !is.na(group) & !is.na(pam50))),")"))

# in DCIS
labels_dcis = read.table(paste0(basedir,"khoulaha/BreastLandscape/clinical/DCIS_clinical.csv"), header = T, sep=",")
colnames(labels_dcis)[which(colnames(labels_dcis) == "IC10_RNA")] = "ic10"
labels_dcis = labels_dcis[which(labels_dcis$Clustering_RAHBT  == 1 | labels_dcis$Clustering_TBCRC == 1),] #stick true to the original paper
any(duplicated(labels_dcis$Patient_ID))

labels_dcis$ic10[which(!is.na(labels_dcis$ic10))] = paste0("ic", labels_dcis$ic10[which(!is.na(labels_dcis$ic10))])

labels_dcis$ic10[which(labels_dcis$ic10 == "ic4" & labels_dcis$ER_RNA == "+")] = "IC4ER+"
labels_dcis$ic10[which(labels_dcis$ic10 == "ic4" & labels_dcis$ER_RNA == "-")] = "IC4ER-"
labels_dcis$ic10[which(startsWith(labels_dcis$ic10,"ic"))] = gsub('ic',"IC",labels_dcis$ic10[which(startsWith(labels_dcis$ic10,"ic"))])
table(labels_dcis$ic10, useNA = "always")

labels_dcis = labels_dcis[which(!is.na(labels_dcis$ic10)),]
table(duplicated(labels_dcis$Patient_ID))
sample_dup = labels_dcis$Patient_ID[which(duplicated(labels_dcis$Patient_ID))]

labels_dcis %>% dplyr::filter(Patient_ID %in% sample_dup) %>% dplyr::select(Patient_ID, Sample_ID, Cohort, Diagnostic_Group) %>% dplyr::arrange(Patient_ID)

table(labels_dcis$Diagnostic_Group[which(!labels_dcis$Patient_ID %in% sample_dup)])
table(labels_dcis$Diagnostic_Group[which(duplicated(labels_dcis$Patient_ID))])

labels_dcis$group = NA
labels_dcis$group[which(labels_dcis$ic10 == "IC5")] = "HER2+"
labels_dcis$group[which(labels_dcis$ic10 %in% c("IC1", "IC2", "IC6", "IC9"))] = "ER+ High"
labels_dcis$group[which(labels_dcis$ic10 %in% c("IC3", "IC4ER+","IC7", "IC8"))] = "ER+ Typical"
labels_dcis$group[which(labels_dcis$ic10 %in% c("IC10","IC4ER-"))] = "TNBC"
table(labels_dcis$group, useNA = "always")

ed_1g_left <- labels_dcis %>% 
  dplyr::filter(!is.na(ic10) & !is.na(group) & !is.na(PAM50)) %>%
  make_long(ic10, group, PAM50) %>%
  dplyr::mutate(node = factor(node, levels = rev(c(paste0("IC",c(1,2,6,9,3,7,"4ER+",8,5,10,"4ER-")),
                                                   "ER+ High","ER+ Typical","HER2+","TNBC",
                                                   "LumB","LumA","Normal","Her2","Basal")))) %>%
  dplyr::mutate(next_node = factor(next_node, levels = rev(c(paste0("IC",c(1,2,6,9,3,7,"4ER+",8,5,10,"4ER-")),
                                                             "ER+ High","ER+ Typical","HER2+","TNBC",
                                                             "LumB","LumA","Normal","Her2","Basal")))) %>%
  ggplot(aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = factor(node),
             label = node)) +
  geom_sankey(flow.alpha = 0.6, node.color = 1) +
  geom_sankey_label(size = 2) +
  scale_x_discrete(labels = c("ic10" = "iC10 IC", "group" = "IC subgroup", "PAM50" = "PAM50")) +
  scale_fill_manual(name = "", values = c(col.ic[which(names(col.ic) %in% labels_dcis$ic10)],
                                          col.group,"Basal" = "#7a5195", "Her2" = "#ef5675", "LumA" = "#003f5c", "LumB" = "#ffa600","Normal" = "#ef94d5")) +
  xlab('') +
  theme_sankey(base_size = 5) + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 7), 
        axis.text.x = element_text(size=6)) + 
  labs(title = paste0("DCIS-RNA (n=",nrow(labels_dcis %>% dplyr::filter(!is.na(ic10) & !is.na(group) & !is.na(PAM50))),")"))

### SAVE ##########################################################################################
## Figures----
png(filename = file.path(outdir, 'ExtendedData1f.png'), res = 300, width = 7.81, height = 5.38, units = 'in')
ed_1f
dev.off()

pdf(paste0(outdir, 'ExtendedData1g_left.pdf'), width = (6.97/1.5)*0.5, height = (4.65/1.2)*0.5)
ed_1g_left
dev.off()

pdf(paste0(outdir, 'ExtendedData1g_middle.pdf'), width = (6.97/1.5)*0.5, height = (4.65/1.3)*0.5)
ed_1g_middle
dev.off()

pdf(paste0(outdir, 'ExtendedData1g_right.pdf'), width = (6.97/1.5)*0.5, height = (4.65/1.3)*0.5)
ed_1g_right
dev.off()

## SourceData----
write.table(sourcetable_ed_1f, paste0(basedir,"lisem/bc_landscape/submission/data/ExtendedData1f_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

df_mets$Stage = "Metastatic"
tcga_pam50$Stage = "Primary"
labels_dcis$Stage = "DCIS"

colnames(tcga_pam50)[1] = "Sample"
labels_dcis = labels_dcis[,c(1,23,11,22,10)]
colnames(labels_dcis) = c("Sample","Stage","ENiClust", "group", "pam50")

sourcetable = rbind(labels_dcis[,c("Sample","Stage","ENiClust", "group", "pam50")],rbind(tcga_pam50[,c("Sample","Stage","ENiClust", "group", "pam50")], df_mets[,c("Sample","Stage","ENiClust", "group", "pam50")]))
colnames(sourcetable) = c("Sample","Stage","Subtype","Subgroup","PAM50")
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/submission/data/ExtendedData1g_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

## Calculating the delta----
df_dcis <- labels_dcis %>% dplyr::filter(!is.na(ENiClust) & !is.na(group) & !is.na(pam50))
df_prim <- tcga_pam50 %>% dplyr::filter(!is.na(ENiClust) & !is.na(group) & !is.na(pam50))

length(which(df_prim$pam50 == "LumB"))/length(which(df_prim$pam50 %in% c("LumA","LumB"))) - length(which(df_dcis$pam50 == "LumB"))/length(which(df_dcis$pam50 %in% c("LumA","LumB"))) #11%
length(which(df_mets$pam50 == "LumB"))/length(which(df_mets$pam50 %in% c("LumA","LumB"))) - length(which(df_prim$pam50 == "LumB"))/length(which(df_prim$pam50 %in% c("LumA","LumB"))) #29%

length(which(df_prim$group == "ER+ High"))/length(which(df_prim$group %in% c("ER+ Typical","ER+ High"))) - length(which(df_dcis$group == "ER+ High"))/length(which(df_dcis$group %in% c("ER+ Typical","ER+ High")))
length(which(df_mets$group == "ER+ High"))/length(which(df_mets$group %in% c("ER+ Typical","ER+ High"))) - length(which(df_prim$group == "ER+ High"))/length(which(df_prim$group %in% c("ER+ Typical","ER+ High"))) #9%
length(which(df_mets$group == "ER+ High"))/length(which(df_mets$group %in% c("ER+ Typical","ER+ High"))) - length(which(df_dcis$group == "ER+ High"))/length(which(df_dcis$group %in% c("ER+ Typical","ER+ High")))

# Looking at proportions
length(which(df_mets$pam50 == "LumB" & df_mets$group == "ER+ High"))/length(which(df_mets$pam50 == "LumB"))
length(which(df_mets$pam50 == "LumB" & df_mets$group == "ER+ Typical"))/length(which(df_mets$pam50 == "LumB"))

length(which(df_mets$group == "ER+ High" & df_mets$pam50 == "LumB"))/length(which(df_mets$group == "ER+ High"))
length(which(df_mets$group == "ER+ Typical" & df_mets$pam50 == "LumB"))/length(which(df_mets$group == "ER+ Typical"))
