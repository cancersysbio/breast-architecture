### CREATE FIGURE 1c #############################################################################
# creates figure 1c
# figure 1c provides the Kaplan-Meier curves of distant relapse free (DRF) survival of the ER+ Typical-risk and ER+ High-risk classes detected by the four IC subtypes classifiers (ENiClust, iC10 DNA+RNA, iC10 DNA-only, and iC10 RNA-only). 

### PREAMBLE #####################################################################################
library(yaml)
require(survival)
require(survminer)
require(patchwork)
require(scales)
library(formattable)
library(stringr)
library(ggplot2)
library(ggpubr)

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
path.git <- paste0(basedir,"lisem/git_repo/brcarepred/")

### Comparison with iC10 DNA-only and RNA-only----
clinic = read.table(paste0(basedir,"lisem/bc_landscape/github/clinic_METABRIC.txt"), header = T, sep = "\t") # see: clinic_METABRIC.R
table(clinic$eniclust, useNA = "always")
table(is.na(clinic$eniclust))
clinic = clinic[which(!is.na(clinic$eniclust)),]

clinic$type <- NA
clinic$type[which(clinic$eniclust %in% c("ic3","ic4ER+","ic7","ic8","Other"))] <- "Typical-risk"
clinic$type[which(clinic$eniclust %in% c("ic1","ic2","ic6","ic9"))] <- "High-risk"
table(clinic$type, useNA = "always")

fit_drs_1 <- survfit(Surv(TDR_20, DR_20) ~ type, data = clinic)
dfs_1 <- ggsurvplot(fit_drs_1,conf.int = T,pval = F,risk.table = F,
                    palette = c("#ff6a00","#7aa6c2"), legend.title="", legend.labs=c("High-risk","Typical-risk")) + 
  ggtitle("ER+ ICs - ENiClust") + xlab("Time (years)") + ylab("DRF Survival")

clinic[,"ICs"] <- clinic$type
clinic$ICs <- factor(clinic$ICs, levels = c("Typical-risk","High-risk"))
fit_coxph_1 = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", "ICs+AGE+GRADE+SIZE+ER.Status_LM+LN")), data = clinic)
plot_forest_1 = ggforest(fit_coxph_1, main = 'Time to distant relapse / ENiClust (METABRIC)', data = clinic)

p_cox = as.character(summary(fit_coxph_1)$coefficients["ICsHigh-risk","Pr(>|z|)"])
hr_cox = summary(fit_coxph_1)$coefficients["ICsHigh-risk","exp(coef)"]
dfs_1 <- dfs_1$plot + 
  ggplot2::annotate(geom="text", x=3, y=0.2, label=paste0("P = " ,str_sub(p_cox,1,4),"x10-",unlist(strsplit(p_cox,"e-0"))[2]), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=3, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

# iC10 CNA-only
clinic$type <- NA
clinic$type[which(clinic$iC10_DNA %in% c("ic3","ic4ER+","ic7","ic8","Other"))] <- "Typical-risk"
clinic$type[which(clinic$iC10_DNA %in% c("ic1","ic2","ic6","ic9"))] <- "High-risk"
clinic$type[which(is.na(clinic$eniclust))] <- NA
table(clinic$type, useNA = "always")

fit_drs_2 <- survfit(Surv(TDR_20, DR_20) ~ type, data = clinic)
dfs_2 <- ggsurvplot(fit_drs_2,conf.int = T,pval = F,risk.table = F,
                    palette = c("#ff6a00","#7aa6c2"), legend.title="", legend.labs=c("High-risk","Typical-risk")) + 
  ggtitle("ER+ ICs - iC10 DNA-only") + xlab("Time (years)") + ylab("DRF Survival")

clinic[,"ICs"] <- clinic$type
clinic$ICs <- factor(clinic$ICs, levels = c("Typical-risk","High-risk"))
fit_coxph_2 = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", "ICs+AGE+GRADE+SIZE+ER.Status_LM+LN")), data = clinic)
plot_forest_2 = ggforest(fit_coxph_2, main = 'Time to distant relapse / iC10 DNA-only (METABRIC)', data = clinic)

p_cox = as.character(summary(fit_coxph_2)$coefficients["ICsHigh-risk","Pr(>|z|)"]*10^4)
hr_cox = summary(fit_coxph_2)$coefficients["ICsHigh-risk","exp(coef)"]
dfs_2 <- dfs_2$plot + 
  ggplot2::annotate(geom="text", x=3, y=0.2, label=paste0("P = " ,str_sub(p_cox,1,4),"x10-4"), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=3, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

# iC10 RNA-only
clinic$type <- NA
clinic$type[which(clinic$iC10_RNA %in% c("ic3","ic4ER+","ic7","ic8","Other"))] <- "Typical-risk"
clinic$type[which(clinic$iC10_RNA %in% c("ic1","ic2","ic6","ic9"))] <- "High-risk"
clinic$type[which(is.na(clinic$eniclust))] <- NA
table(clinic$type, useNA = "always")

fit_drs_3 <- survfit(Surv(TDR_20, DR_20) ~ type, data = clinic)
dfs_3 <- ggsurvplot(fit_drs_3,conf.int = T,pval = F,risk.table = F,
                    palette = c("#ff6a00","#7aa6c2"), legend.title="", legend.labs=c("High-risk","Typical-risk")) + 
  ggtitle("ER+ ICs - iC10 RNA-only") + xlab("Time (years)") + ylab("DRF Survival")

clinic[,"ICs"] <- clinic$type
clinic$ICs <- factor(clinic$ICs, levels = c("Typical-risk","High-risk"))
fit_coxph_3 = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", "ICs+AGE+GRADE+SIZE+ER.Status_LM+LN")), data = clinic)
plot_forest_3 = ggforest(fit_coxph_3, main = 'Time to distant relapse / iC10 RNA-only (METABRIC)', data = clinic)

p_cox = as.character(summary(fit_coxph_3)$coefficients["ICsHigh-risk","Pr(>|z|)"])
hr_cox = summary(fit_coxph_3)$coefficients["ICsHigh-risk","exp(coef)"]
dfs_3 <- dfs_3$plot + 
  ggplot2::annotate(geom="text", x=3, y=0.2, label=paste0("P = " ,str_sub(p_cox,1,4),"x10-",unlist(strsplit(p_cox,"e-0"))[2]), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=3, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

# iC10 DNA+RNA
clinic$type <- NA
clinic$type[which(paste0("ic",clinic$iC10) %in% c("ic3","ic4ER+","ic7","ic8","Other"))] <- "Typical-risk"
clinic$type[which(paste0("ic",clinic$iC10) %in% c("ic1","ic2","ic6","ic9"))] <- "High-risk"
clinic$type[which(is.na(clinic$eniclust))] <- NA
table(clinic$type, useNA = "always")

fit_drs_4 <- survfit(Surv(TDR_20, DR_20) ~ type, data = clinic)
dfs_4 <- ggsurvplot(fit_drs_4,conf.int = T,pval = F,risk.table = F,
                    palette = c("#ff6a00","#7aa6c2"), legend.title="", legend.labs=c("High-risk","Typical-risk")) + 
  ggtitle("ER+ ICs - iC10 DNA+RNA") + xlab("Time (years)") + ylab("DRF Survival")

clinic[,"ICs"] <- clinic$type
clinic$ICs <- factor(clinic$ICs, levels = c("Typical-risk","High-risk"))
fit_coxph_4 = coxph(formula(paste0("Surv(TDR_20, DR_20) ~", "ICs+AGE+GRADE+SIZE+ER.Status_LM+LN")), data = clinic)
plot_forest_4 = ggforest(fit_coxph_4, main = 'Time to distant relapse / iC10 DNA+RNA (METABRIC)', data = clinic)

p_cox = as.character(summary(fit_coxph_4)$coefficients["ICsHigh-risk","Pr(>|z|)"])
hr_cox = summary(fit_coxph_4)$coefficients["ICsHigh-risk","exp(coef)"]
dfs_4 <- dfs_4$plot + 
  ggplot2::annotate(geom="text", x=3, y=0.2, label=paste0("P = " ,str_sub(p_cox,1,4),"x10-",unlist(strsplit(p_cox,"e-0"))[2]), size=5, hjust=0) +
  ggplot2::annotate(geom="text", x=3, y=0.25, label=paste0("HR = " ,round(hr_cox,digits=2)), size=5, hjust=0) 

### SAVE ##########################################################################################
## Figure----
svg(filename = paste0(basedir,"lisem/bc_landscape/submission/figures/Low_as_High_dfs.svg"), w=7, h=6.44)
dfs_1 + dfs_4 + dfs_2 + dfs_3
dev.off()

## Save CoxModels----
fit.coxph <- list("ENiClust"=fit_coxph_1, "iC10 DNA+RNA"=fit_coxph_4, "iC10 DNA-only"=fit_coxph_2, "iC10 RNA-only"=fit_coxph_3)
save(fit.coxph, file=paste0(basedir,"lisem/bc_landscape/submission/data/SuppFigure1f_fit.coxph.RData"))

## Save KM fits----
fit.drs <- list("ENiClust"=fit_drs_1, "iC10 DNA+RNA"=fit_drs_4, "iC10 DNA-only"=fit_drs_2, "iC10 RNA-only"=fit_drs_3)
save(fit.drs, file=paste0(basedir,"lisem/bc_landscape/submission/data/Figure1d_fit.drs.RData"))

## SourceData----
sourcetable = clinic[,c("METABRIC.ID","TDR_20", "DR_20","AGE","GRADE","SIZE","ER.Status_LM","LN","eniclust", "iC10", "iC10_DNA", "iC10_RNA")]
colnames(sourcetable) = c("METABRIC.ID","TDR_20", "DR_20","AGE","GRADE","SIZE","ER.Status","LN","ENiClust", "iC10 DNA+RNA", "iC10 DNA-only", "iC10 RNA-only")
sourcetable$ENiClust[which(sourcetable$ENiClust == "Other")] = "ic3/ic7"
sourcetable[,"iC10 DNA+RNA"] = paste0("ic", sourcetable[,"iC10 DNA+RNA"])

write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/submission/data/Figure1c_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
