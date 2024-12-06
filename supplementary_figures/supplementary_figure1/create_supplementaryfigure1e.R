### CREATE SUPPLEMENTARY FIGURE 1e #############################################################################
# creates supplementary figure 1e
# supplementary figure 1e provides the Kaplan-Meier curves of distant relapse free (DRF) survival of the ER+ Typical-risk reclassified as ER+ High-risk and the ER+ High-risk reclassified as ER+ Typical-risk by ENiClust (purple) vs. iC10 DNA-only (blue) or compared with the reference (iC10 DNA+RNA, red).

### PREAMBLE #####################################################################################
library(yaml)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(dplyr)
library(tibble)
library(tidyverse)
library(caret)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'
outdir <- "/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/submission/figures/"

## Outcome in METABRIC
Typical = c("3", "4ER+", "7", "8", "3/7")
High = c("1", "2", "6", "9")

clinic = read.table(paste0(basedir,"lisem/bc_landscape/github/clinic_METABRIC.txt"), header = T, sep = "\t") # see: clinic_METABRIC.R

sourcetable = clinic[,c("METABRIC.ID","TDR_20", "DR_20","AGE","GRADE","SIZE","ER.Status_LM","LN","eniclust", "iC10", "iC10_DNA", "iC10_RNA")]
colnames(sourcetable) = c("METABRIC.ID","TDR_20", "DR_20","AGE","GRADE","SIZE","ER.Status","LN","ENiClust", "iC10 DNA+RNA", "iC10 DNA-only", "iC10 RNA-only")
sourcetable$ENiClust[which(sourcetable$ENiClust == "Other")] = "ic3/ic7"
sourcetable[,"iC10 DNA+RNA"] = paste0("ic", sourcetable[,"iC10 DNA+RNA"])

clinic$iC10_DNA = gsub("ic","", clinic$iC10_DNA)
clinic$eniclust[which(clinic$eniclust != "Other")] = sapply(strsplit(clinic$eniclust[which(clinic$eniclust != "Other")],"ic"),"[[",2)
clinic$eniclust[which(clinic$eniclust == "4")] = paste0("4ER",ifelse(is.na(clinic$ER.Status[which(clinic$eniclust == "4")]),clinic$ER.Expr[which(clinic$eniclust == "4")],ifelse(clinic$ER.Status[which(clinic$eniclust == "4")] == "pos","+","-")) )
clinic$eniclust[which(clinic$eniclust == "Other")] = "3/7"

# Mis-classified ER+ Typical vs. ER+ High
# Typical -> High
# iC10
clinic$type = NA
clinic$type[which(clinic$iC10 %in% Typical & clinic$eniclust %in% Typical)] = "Typical_as_Typical"
clinic$type[which(clinic$iC10 %in% Typical & clinic$eniclust %in% High)] = "Typical_as_High"
clinic$type = factor(clinic$type, levels = c("Typical_as_High","Typical_as_Typical"))

fit_drf <- survfit(Surv(TDR_20, DR_20) ~ type, data = clinic) # Distant relapse free survival
cox_fit = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", "type+AGE+GRADE+SIZE+ER.Status_LM+LN")), data = clinic)

dfs_low <- ggsurvplot(fit_drf,conf.int = F,pval = F,risk.table = F,
                      palette = c("#825CA6","red"), legend.title="", legend.labs=c("ER+ Typical as High by ENiClust","ER+ Typical by iC10")) + 
  ggtitle("ER+ Typical re-classified as High by ENiClust") + xlab("Time (years)") + ylab("DRF Survival")

p_cox = summary(cox_fit)$coefficients["typeTypical_as_Typical","Pr(>|z|)"]
hr_cox = summary(cox_fit)$coefficients["typeTypical_as_Typical","exp(coef)"]
dfs_low <- dfs_low$plot + 
  ggplot2::annotate(geom="text", x=0, y=0.2, label=paste0("P = " ,round(p_cox,digits=2)), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=0, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

# iC10 DNA-only
clinic$type = NA
clinic$type[which(clinic$iC10 %in% Typical & clinic$eniclust %in% High)] = "ENiClust"
clinic_eniclust = clinic[which(!is.na(clinic$type)),]
clinic$type = NA
clinic$type[which(clinic$iC10 %in% Typical & clinic$iC10_DNA %in% High)] = "iC10 DNA-only"
clinic_iC10_DNA = clinic[which(!is.na(clinic$type)),]

data = rbind(clinic_eniclust, clinic_iC10_DNA)
data = data[which(!data$METABRIC.ID %in% data$METABRIC.ID[which(duplicated(data$METABRIC.ID))]),]
data$type = factor(data$type, levels = c("ENiClust","iC10 DNA-only"))

mydata = data
fit_drf <- survfit(Surv(TDR_20, DR_20) ~ type, data = mydata) # Distant relapse free survival
cox_fit = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", ifelse(length(unique(data$ER.Status_LM)) == 1, "type+AGE+GRADE+SIZE+LN","type+AGE+GRADE+SIZE+ER.Status_LM+LN"))), data = data)

dfs_1 <- ggsurvplot(fit_drf,conf.int = F,pval = F,risk.table = F,
                    palette = c("#825CA6","#317EC2"), legend.title="", legend.labs=c("ENiClust","iC10 DNA-only")) + 
  ggtitle("ER+ Typical re-classified as High") + xlab("Time (years)") + ylab("DRF Survival")

p_cox = summary(cox_fit)$coefficients["typeiC10 DNA-only","Pr(>|z|)"]
hr_cox = summary(cox_fit)$coefficients["typeiC10 DNA-only","exp(coef)"]
dfs_1 <- dfs_1$plot + 
  ggplot2::annotate(geom="text", x=0, y=0.2, label=paste0("P = " ,round(p_cox,digits=2)), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=0, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

# High -> Typical
# iC10
clinic$type = NA
clinic$type[which(clinic$iC10 %in% High & clinic$eniclust %in% High)] = "High_as_High"
clinic$type[which(clinic$iC10 %in% High & clinic$eniclust %in% Typical)] = "High_as_Typical"
clinic$type = factor(clinic$type, levels = c("High_as_Typical","High_as_High"))

fit_drf <- survfit(Surv(TDR_20, DR_20) ~ type, data = clinic) # Distant relapse free survival
cox_fit = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", "type+AGE+GRADE+SIZE+ER.Status_LM+LN")), data = clinic)

dfs_high <- ggsurvplot(fit_drf,conf.int = F,pval = F,risk.table = F,
                       palette = c("#825CA6","red"), legend.title="", legend.labs=c("ER+ High as Typical by ENiClust","ER+ High by iC10")) + 
  ggtitle("ER+ High re-classified as Typical by ENiClust") + xlab("Time (years)") + ylab("DRF Survival")

p_cox = summary(cox_fit)$coefficients["typeHigh_as_High","Pr(>|z|)"]
hr_cox = summary(cox_fit)$coefficients["typeHigh_as_High","exp(coef)"]
dfs_high <- dfs_high$plot + 
  ggplot2::annotate(geom="text", x=0, y=0.2, label=paste0("P = " ,round(p_cox,digits=2)), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=0, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

# iC10 DNA-only
clinic$type = NA
clinic$type[which(clinic$iC10 %in% High & clinic$eniclust %in% Typical)] = "ENiClust"
clinic_eniclust = clinic[which(!is.na(clinic$type)),]
clinic$type = NA
clinic$type[which(clinic$iC10 %in% High & clinic$iC10_DNA %in% Typical)] = "iC10 DNA-only"
clinic_iC10_DNA = clinic[which(!is.na(clinic$type)),]

data = rbind(clinic_eniclust, clinic_iC10_DNA)
data = data[which(!data$METABRIC.ID %in% data$METABRIC.ID[which(duplicated(data$METABRIC.ID))]),]
data$type = factor(data$type, levels = c("iC10 DNA-only","ENiClust"))

mydata = data
fit_drf <- survfit(Surv(TDR_20, DR_20) ~ type, data = mydata) # Distant relapse free survival
cox_fit = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", "type+AGE+GRADE+SIZE+ER.Status_LM+LN")), data = data)

dfs_2 <- ggsurvplot(fit_drf,conf.int = F,pval = F,risk.table = F,
                    palette = c("#825CA6","#317EC2"), legend.title="", legend.labs=c("ENiClust","iC10 DNA-only")) + 
  ggtitle("ER+ High re-classified as Typical") + xlab("Time (years)") + ylab("DRF Survival")

p_cox = summary(cox_fit)$coefficients["typeENiClust","Pr(>|z|)"]
hr_cox = summary(cox_fit)$coefficients["typeENiClust","exp(coef)"]
dfs_2 <- dfs_2$plot + 
  ggplot2::annotate(geom="text", x=0, y=0.2, label=paste0("P = " ,round(p_cox,digits=2)), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=0, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

### SAVE ##########################################################################################
## Figures----
svg(filename = paste0(outdir,"SuppFigure1e.svg"), w=7,h=6.44)
(dfs_1 + dfs_low) / (dfs_2 + dfs_high)
dev.off()

## SourceData----
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/github/SupplementaryFigure1e_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

