### CREATE EXTENDED DATA FIGURE 1b, 1c and 1d #############################################################################
# creates extended data figure 1b, 1c and 1d
# extended data figure 1b provides the number of cases for each IC subgroup across ancestries at the primary and metastatic stages
# extended data figure 1c provides the proportion of the IC subgroups and IC subtypes across ancestries in all breast cancer patients only at the primary and metastatic stages
# extended data figure 1d provides the proportion of the IC subgroups and IC subtypes across ancestries in ER+ breast cancer patients only at the primary and metastatic stages

### PREAMBLE #####################################################################################
library(yaml)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggrepel)

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

col.group = c("ER+ High" = "#fa954eff","ER+ Typical" = "#a1c1d4ff","HER2+" = "#8b0100ff","TNBC" = "#6300a1")

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'
outdir <- '/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/publication/figures/'

### In primary: TCGA is the only one having the ancestry annotations----
## Loading data
tcga_ic10 = read.table(paste0(basedir,"lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v14/Pipeline_2/tcga_predictions.txt"), header = T, sep = "")
er_status <- read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-02-09_rna_megatable_primary_met_eniclust_v2_er_ihc.txt"),sep="\t",header=T) %>%
  dplyr::select(Sample, ER) %>%
  dplyr::filter(!is.na(ER))
er_status <- read.table(paste0(basedir,"lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v14/Pipeline_2/tcga_transformed_matrix.txt"),sep="\t",header=T) %>%
  dplyr::filter(!Sample %in% er_status$Sample) %>% dplyr::select(Sample, ER) %>% 
  dplyr::rowwise() %>% dplyr::mutate(ER = ER/4) %>%
  dplyr::bind_rows(er_status)

table(tcga_ic10$Sample %in% er_status$Sample)

tcga_ic10$voting[which(tcga_ic10$voting == "ic4" & tcga_ic10$Sample %in% er_status$Sample[which(er_status$ER == 1)])] = "ic4ER+"
tcga_ic10$voting[which(tcga_ic10$voting == "ic4" & tcga_ic10$Sample %in% er_status$Sample[which(er_status$ER == 0)])] = "ic4ER-"
tcga_ic10$voting[which(tcga_ic10$voting == "Other")] = "ic3/ic7"
table(tcga_ic10$voting)

tcga_ic10$group = NA
tcga_ic10$group[which(tcga_ic10$voting == "ic5")] = "HER2+"
tcga_ic10$group[which(tcga_ic10$voting %in% c("ic1", "ic2", "ic6", "ic9"))] = "ER+ High"
tcga_ic10$group[which(tcga_ic10$voting %in% c("ic3/ic7", "ic4ER+", "ic8"))] = "ER+ Typical"
tcga_ic10$group[which(tcga_ic10$voting %in% c("ic10", "ic4ER-"))] = "TNBC"

clinic = read.table(paste0(basedir,"srinivap/Work/Data/TCGA/Clinical/TCGA-CDR-SupplementalTableS1.csv"),header = T,sep="\t",stringsAsFactors = F,quote = "",fill=T, )
clinic = clinic[which(clinic$type == "BRCA"),]

ancestry = read.delim(paste0(basedir,"khoulaha/BreastLandscape/data/yuan_tcga_breast_ancestry.csv"), sep = ",", header = F)

all(tcga_ic10$Sample %in% ancestry$V1)

tcga_ic10$ancestry <- sapply(tcga_ic10$Sample, function(x) ancestry$V5[which(ancestry$V1 == x)])

tcga_ic10$ancestry[which(tcga_ic10$ancestry %in% c("OA"))] <- NA
tcga_ic10$ancestry[which(tcga_ic10$ancestry == "AA")] <- "AFRICAN"
tcga_ic10$ancestry[which(tcga_ic10$ancestry == "EA")] <- "EUROPEAN"
tcga_ic10$ancestry[which(tcga_ic10$ancestry == "EAA")] <- "EAST ASIAN"
table(tcga_ic10$ancestry)

## Plotting
df_calculated <- tcga_ic10[which(!is.na(tcga_ic10$group) & !is.na(tcga_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, group) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","EAST ASIAN"))
df_calculated$group = factor(df_calculated$group, levels = c("ER+ High","ER+ Typical","HER2+","TNBC"))

ed_1c_top_gp = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = group)) +
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.group) +
  ylab("Proportion of samples") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry),inherit.aes=FALSE,label.size=NA,alpha=0,size=2) +
  theme_LM + ggtitle("Primary (TCGA)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_text(size=5), 
        axis.title = element_text(size=7), 
        title = element_text(size=5), 
        legend.text = element_text(size=6))

df_calculated <- tcga_ic10[which(!is.na(tcga_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, voting) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","EAST ASIAN"))
df_calculated$voting = factor(df_calculated$voting, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8","ic5","ic10","ic4ER-"))

ed_1c_top_ic = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = voting))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.ic[which(names(col.ic) %in% df_calculated$voting)],
                    labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-")) +
  ylab("") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry),inherit.aes = F,label.size=NA,alpha=0,size=2) +
  theme_LM + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size=7), 
        legend.text = element_text(size=6))

mylegend = cowplot::get_legend(ggplot(tcga_ic10,aes(x=voting_proba_ic1,y=voting_proba_ic10)) +
                                 geom_col(aes(fill=group)) +
                                 scale_fill_manual(name = "IC subgroups", values = col.group) +
                                 new_scale_fill() +
                                 geom_col(aes(fill=factor(voting, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8","ic5","ic10","ic4ER-")))) +
                                 scale_fill_manual(name = "IC subtypes", values = col.ic,
                                                   labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-"))+
                                 theme(legend.text = element_text(size=14), legend.title = element_text(size=16)) )

# ER+ only
df_calculated <- tcga_ic10[which(grepl("ER[+]",tcga_ic10$group) & !is.na(tcga_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, group) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","EAST ASIAN"))
df_calculated$group = factor(df_calculated$group, levels = c("ER+ High","ER+ Typical"))

ed_1d_top_gp = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = group))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.group) +
  ylab("Proportion of ER+ samples") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry), inherit.aes = F,label.size=NA,alpha=0,size=2) +
  theme_LM + ggtitle("Primary (TCGA)") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_text(size=5), 
        axis.title = element_text(size=7), 
        title = element_text(size=7), 
        legend.text = element_text(size=6))

df_calculated <- tcga_ic10[which(grepl("ER[+]",tcga_ic10$group) & !is.na(tcga_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, voting) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","EAST ASIAN"))
df_calculated$voting = factor(df_calculated$voting, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8"))

ed_1d_top_ic = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = voting))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.ic[which(names(col.ic) %in% df_calculated$voting)],
                    labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-")) +
  ylab("") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry), inherit.aes = F,label.size=NA,alpha=0,size=2) +
  theme_LM + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size=7), 
        legend.text = element_text(size=6))

### In metastatic----
hmf_ancestry_pca <- read.table(paste0(basedir,"khoulaha/GEB/hartwig/germline/VariantCalling/AIM/2021-09-24_hartwig_PCs.txt"), sep = "\t", header = T)
head(hmf_ancestry_pca)

ggplot(hmf_ancestry_pca, aes(x=PC1,y=PC2)) +
  geom_point()

data = hmf_ancestry_pca[which(!is.na(hmf_ancestry_pca$PC1) & !is.na(hmf_ancestry_pca$PC2)),2:3]
rownames(data) = hmf_ancestry_pca$sample[which(!is.na(hmf_ancestry_pca$PC1) & !is.na(hmf_ancestry_pca$PC2))]

set.seed(42)
k = 5
clustering <- kmeans(data,k)

df = rbind(data,clustering$centers)
df$col <- "point"
df$col[(nrow(df)-(k-1)):nrow(df)] <- "center"
df$label <- NA
df$label[which(df$col == "center")] <- 1:5
ggplot(df, aes(x=PC1,y=PC2, color = col)) +
  geom_point() +
  scale_color_manual(name="",values=c('point'='black','center'='red')) +
  geom_text_repel(aes(label = label)) +
  theme_LM + theme(legend.position = "none")

# ancestry
hmf_map_id = read.table(paste0(basedir, "khoulaha/BreastLandscape/clinical/Hartwig_clinical.txt"), header = T, sep = "\t")
table(rownames(data) %in% hmf_map_id$X.patientId)
white <- unlist(sapply(rownames(data)[which(clustering$cluster %in% c(1,4))], function(x) if(x %in% hmf_map_id$X.patientId){hmf_map_id$Individual.System.ID[which(hmf_map_id$X.patientId == x)]}))
asian <- unlist(sapply(rownames(data)[which(clustering$cluster == 2)], function(x) if(x %in% hmf_map_id$X.patientId){hmf_map_id$Individual.System.ID[which(hmf_map_id$X.patientId == x)]}))
black <- unlist(sapply(rownames(data)[which(clustering$cluster == 3)], function(x) if(x %in% hmf_map_id$X.patientId){hmf_map_id$Individual.System.ID[which(hmf_map_id$X.patientId == x)]}))

# load data
hmf_ic10 = read.table(paste0(basedir,"lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v14/Pipeline_2/hmf_predictions.txt"), header = T, sep = "")
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")
hmf_ic10 = hmf_ic10[which(hmf_ic10$Sample %in% mega_2$Sample),] #looking only at Metastatic stage
rm(mega_2)
er_status <- read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"),sep="\t",header=T) %>%
  dplyr::select(Sample, ER) %>%
  dplyr::filter(!is.na(ER))
er_status <- read.table(paste0(basedir,"lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v14/Pipeline_2/hmf_transformed_matrix.txt"),sep="\t",header=T) %>%
  dplyr::filter(!Sample %in% er_status$Sample) %>% dplyr::select(Sample, ER) %>% 
  dplyr::rowwise() %>% dplyr::mutate(ER = ER/4) %>%
  dplyr::bind_rows(er_status)

table(hmf_ic10$Sample %in% er_status$Sample)

hmf_ic10$voting[which(hmf_ic10$voting == "ic4" & hmf_ic10$Sample %in% er_status$Sample[which(er_status$ER == 1)])] = "ic4ER+"
hmf_ic10$voting[which(hmf_ic10$voting == "ic4" & hmf_ic10$Sample %in% er_status$Sample[which(er_status$ER == 0)])] = "ic4ER-"
hmf_ic10$voting[which(hmf_ic10$voting == "Other")] = "ic3/ic7"
table(hmf_ic10$voting)

hmf_ic10$group = NA
hmf_ic10$group[which(hmf_ic10$voting == "ic5")] = "HER2+"
hmf_ic10$group[which(hmf_ic10$voting %in% c("ic1", "ic2", "ic6", "ic9"))] = "ER+ High"
hmf_ic10$group[which(hmf_ic10$voting %in% c("ic3/ic7", "ic4ER+", "ic8"))] = "ER+ Typical"
hmf_ic10$group[which(hmf_ic10$voting %in% c("ic10", "ic4ER-"))] = "TNBC"

hmf_ic10$ancestry <- NA
hmf_ic10$ancestry[which(hmf_ic10$Sample %in% white)] <- "EUROPEAN"
hmf_ic10$ancestry[which(hmf_ic10$Sample %in% black)] <- "AFRICAN"
hmf_ic10$ancestry[which(hmf_ic10$Sample %in% asian)] <- "ASIAN\nAMERINDIAN\nOCEANIAN"
table(hmf_ic10$ancestry, useNA = "always")

## Plotting
df_calculated <- hmf_ic10[which(!is.na(hmf_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, group) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","ASIAN\nAMERINDIAN\nOCEANIAN"))
df_calculated$group = factor(df_calculated$group, levels = c("ER+ High","ER+ Typical","HER2+","TNBC"))

ed_1c_bottom_gp = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = group))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.group) +
  ylab("Proportion of samples") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry),inherit.aes=FALSE,label.size=NA,alpha=0,size=2) +
  theme_LM + ggtitle("Metastatic") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_text(size=5), 
        axis.title = element_text(size=7), 
        title = element_text(size=7), 
        legend.text = element_text(size=6))

df_calculated <- hmf_ic10[which(!is.na(hmf_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, voting) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","ASIAN\nAMERINDIAN\nOCEANIAN"))
df_calculated$voting = factor(df_calculated$voting, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8","ic5","ic10","ic4ER-"))

ed_1c_bottom_ic = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = voting))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.ic[which(names(col.ic) %in% df_calculated$voting)],
                    labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-")) +
  ylab("") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry),inherit.aes = F,label.size=NA,alpha=0,size=2) +
  theme_LM + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size=7), 
        legend.text = element_text(size=6))

# ER+ only
df_calculated <- hmf_ic10[which(grepl("ER[+]",hmf_ic10$group) & !is.na(hmf_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, group) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","ASIAN\nAMERINDIAN\nOCEANIAN"))
df_calculated$group = factor(df_calculated$group, levels = c("ER+ High","ER+ Typical"))

ed_1d_bottom_gp = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = group))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.group) +
  ylab("Proportion of ER+ samples") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry), inherit.aes = F,label.size=NA,alpha=0,size=2) +
  theme_LM + ggtitle("Metastatic") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_text(size=5), 
        axis.title = element_text(size=7), 
        title = element_text(size=7), 
        legend.text = element_text(size=6))

df_calculated <- hmf_ic10[which(grepl("ER[+]",hmf_ic10$group) & !is.na(hmf_ic10$ancestry)),] %>% 
  dplyr::count(ancestry, voting) %>%
  dplyr::group_by(ancestry) %>%
  dplyr::mutate(Percent = n/sum(n))

df_calculated$ancestry = factor(df_calculated$ancestry, levels = c("EUROPEAN","AFRICAN","ASIAN\nAMERINDIAN\nOCEANIAN"))
df_calculated$voting = factor(df_calculated$voting, levels = c("ic1","ic2","ic6","ic9","ic3/ic7","ic4ER+","ic8"))

ed_1d_bottom_ic = ggplot(df_calculated, aes(x = ancestry, y = Percent, fill = voting))+
  geom_col(position="fill") +
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(name = "", values = col.ic[which(names(col.ic) %in% df_calculated$voting)],
                    labels = c("ic1"="IC1", "ic2"="IC2", "ic6"="IC6", "ic9"="IC9","ic8"="IC8", "ic3/ic7"="IC3/IC7", "ic4ER+"="IC4ER+", "ic5"="IC5", "ic10"="IC10", "ic4ER-"="IC4ER-")) +
  ylab("") +
  xlab("Inferred ancestry") +
  geom_label(data = df_calculated %>% dplyr::group_by(ancestry) %>% dplyr::summarise(Sum = sum(n)), 
             aes(label = Sum, y = 1.05, x = ancestry), inherit.aes = F,label.size=NA,alpha=0,size=2) +
  theme_LM + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_text(size=7), 
        legend.text = element_text(size=6))

### ExtendedData 1b----
tcga_ic10$ancestry = factor(tcga_ic10$ancestry, levels = c("EUROPEAN","AFRICAN","EAST ASIAN"))
p1 <- tcga_ic10 %>% dplyr::filter(!is.na(ancestry)) %>%
  ggplot(aes(x = group, fill = group)) +
  scale_fill_manual(name = "", values = col.group) +
  facet_wrap(.~ancestry, scale = "free_y") +
  xlab('') + ylab('Count') + ggtitle('Primary') +
  geom_bar(position="dodge") + theme_LM + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=16),
        strip.background = element_blank(),
        strip.text = element_text(size=12, vjust = 0, hjust=0),
        title = element_text(size=18),
        legend.position = "none")

hmf_ic10$ancestry = factor(hmf_ic10$ancestry, levels = c("EUROPEAN","AFRICAN","ASIAN\nAMERINDIAN\nOCEANIAN"))
p2 <- hmf_ic10 %>% dplyr::filter(!is.na(ancestry)) %>%
  ggplot(aes(x = group, fill = group)) +
  scale_fill_manual(name = "", values = col.group) +
  facet_wrap(.~ancestry, scale = "free_y") +
  xlab('') + ylab('Count') + ggtitle('Metastatic') +
  geom_bar(position="dodge") + theme_LM + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=14),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=16),
        strip.background = element_blank(),
        strip.text = element_text(size=12, vjust = 0, hjust=0),
        title = element_text(size=18),
        legend.position = "none")

### SAVE ##########################################################################################
## Figures----
pdf(paste0(outdir, 'ExtendedData1b.pdf'), width=10.00, height=5)
p1 + p2
dev.off()

pdf(paste0(outdir, 'ExtendedData1c_top.pdf'), width = (11.5/2.8)*0.55, height = (7/3)*0.7)
(ed_1c_top_gp + rremove("legend")) + (ed_1c_top_ic + rremove("legend")) + mylegend
dev.off()

pdf(paste0(outdir, 'ExtendedData1c_bottom.pdf'), width = (11.5/2.8)*0.55, height = (7/2.8)*0.7)
(ed_1c_bottom_gp + rremove("legend")) + (ed_1c_bottom_ic + rremove("legend")) + mylegend
dev.off()

pdf(paste0(outdir, 'ExtendedData1d_top.pdf'), width = (11.5/2.8)*0.55, height = (7/3)*0.7)
(ed_1d_top_gp + rremove("legend")) + (ed_1d_top_ic + rremove("legend")) + mylegend
dev.off()

pdf(paste0(outdir, 'ExtendedData1d_bottom.pdf'), width = (11.5/2.8)*0.55, height = (7/2.8)*0.7)
(ed_1d_bottom_gp + rremove("legend")) + (ed_1d_bottom_ic + rremove("legend")) + mylegend
dev.off()

## SourceData----
hmf_ic10$Stage = "Metastatic"
tcga_ic10$Stage = "Primary"
sourcetable = rbind(hmf_ic10[,c("Sample","Stage","voting", "group", "ancestry")],tcga_ic10[,c("Sample","Stage","voting", "group", "ancestry")])
colnames(sourcetable) = c("Sample","Stage","Subtype","Subgroup","Ancestry")
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/github/ExtendedData1bcd_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

## Enrichement tests----
# Enrichment for ER+ Typical in WHITE? Yes.
x = length(which(tcga_ic10$ancestry == "EUROPEAN" & tcga_ic10$group == "ER+ Typical"))
y = length(which(tcga_ic10$ancestry != "EUROPEAN" & tcga_ic10$group == "ER+ Typical"))
w = length(which(tcga_ic10$ancestry == "EUROPEAN" & tcga_ic10$group != "ER+ Typical"))
z = length(which(tcga_ic10$ancestry != "EUROPEAN" & tcga_ic10$group != "ER+ Typical"))
mat <- matrix(c(x,w,y,z), nrow = 2, byrow = T)
fisher.test(mat)$p.value

# Enrichment for black between TCGA and Hartwig? Yes.
x = length(which(hmf_ic10$ancestry == "AFRICAN"))
y = length(which(hmf_ic10$ancestry != "AFRICAN"))
w = length(which(tcga_ic10$ancestry == "AFRICAN"))
z = length(which(tcga_ic10$ancestry != "AFRICAN"))
mat <- matrix(c(x,w,y,z), nrow = 2, byrow = T)
fisher.test(mat)$p.value


