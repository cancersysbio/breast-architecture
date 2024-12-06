### CREATE SUPPLEMENTARY FIGURE 1d #############################################################################
# creates supplementary figure 1d
# supplementary figure 1d provides the performance metrics of the DNA-based predictive models iC10 DNA-only and ENiClust compared to iC10 DNA+RNA or ENiClust compared to iC10 DNA-only on the Testing dataset (n=187 TCGA-WES, left), the with-held Validation  datasets: ICGC (n=189 WGS, middle) and METABRIC (n=986 arrays, right).

### PREAMBLE #####################################################################################
library(yaml)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(dplyr)
library(tibble)
library(tidyr)
library(tidyverse)
library(caret)
library(yardstick, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.2.2-Seurat/")
library(iC10TrainingData, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
library(iC10, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")

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
path_eniclust <- "/oak/stanford/groups/ccurtis2/users/lisem/eniclust/wes_wgs/dev/test/rand/"
file <- "rand_59_cv_noLOW_v14"
outdir <- "/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/github/"

data(train.CN)

### Load megatables: can use the last ones without the dup in Hartwig since Hartwig not showed here
mega_1 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_primary_megatable.txt"), header = T, sep = "\t")
mega_2 = read.table(paste0(basedir,"khoulaha/BreastLandscape/data/2024-09-05_metastatic_megatable.txt"), header = T, sep = "\t")
mega = rbind(mega_1[,intersect(colnames(mega_1), colnames(mega_2))], mega_2[,intersect(colnames(mega_1), colnames(mega_2))])

## 1. Performance
load(paste0(basedir,"lisem/bc_landscape/github/references_table.RData"))

balanced_accuracy = data.frame()
performance_table = list()
for(set in c("validation","icgc","metabric")){
  
  performance_table_tmp_all = data.frame()
    
  col_metrics = c(paste0(rep(c("F1 ", "P ", "R "),each=9),c(paste0("IC", c(1,2,6,9)),"Typical",paste0("IC", c(5,10,"4ER-")),"TNBC")), paste0("PR AUC IC", c(1,2,6,9)))
  row_references = c("iC10 DNA+RNA", "iC10 DNA-only", "iC10 DNA-only (All)")[ifelse(set == "icgc",-4,-3)]
  performance_table_tmp = matrix(NA, nrow = length(row_references), ncol = length(col_metrics), dimnames = list(row_references, col_metrics) )
  
  for(ref in c("iC10 DNA+RNA", "iC10 DNA-only", "iC10 DNA-only (All)")[ifelse(set == "icgc",-4,-3)]){
    
    predictions = read.table(paste0(path_eniclust, file, "/Pipeline_2/",set,"_predictions.txt"),sep="\t",header=T)
    
    if(set == "icgc"){
      predictions$Sample = str_sub(predictions$Sample,1,nchar("CURTIS_H000434"))
      er_status = mega %>% dplyr::select(Sample,ER) %>% bind_rows(er_status[,c("Sample","ER")])
    }else{ #for validation, METABRIC uses references_table data
      er_status = read.table(paste0(path_eniclust,file,"/Pipeline_2/validation_set.txt"),sep="\t",header=T)
    }
    
    if(set == "metabric" & ref == "iC10 DNA-only"){
      predictions = predictions[which(!predictions$Sample %in% gsub("[.]","-",colnames(train.CN))),]
    }
    
    if(ref == "iC10 DNA-only (All)"){
      reference = references_table[["iC10 DNA-only"]][["icgc_all"]]
    }else{
      reference = references_table[[ref]][[set]]
    }
    
    if(set == "validation"){
      predictions = predictions[which(grepl("TCGA",predictions$Sample)),]
    }
    if(set %in% c("icgc","metabric")){
      predictions = predictions[which(predictions$Sample %in% reference$Sample),]
      er_status$Sample = sapply(er_status$Sample, function(x) if(set == "icgc"){str_sub(x,1,nchar("CURTIS_H000434"))}else{x})
    }
    
    if(set == "metabric"){
      predictions$voting[which(predictions$voting == "ic4" & predictions$Sample %in% references_table[["iC10 DNA+RNA"]][["metabric"]]$Sample[which(references_table[["iC10 DNA+RNA"]][["metabric"]]$ER.Status_LM == "-")])] = "ic4ER-"
      predictions$voting[which(predictions$voting == "ic4" & predictions$Sample %in% references_table[["iC10 DNA+RNA"]][["metabric"]]$Sample[which(references_table[["iC10 DNA+RNA"]][["metabric"]]$ER.Status_LM == "+")])] = "ic4ER+"
      
      if(ref == "iC10 DNA-only"){
        reference$iC10[which(reference$iC10 == "ic4" & reference$Sample %in% references_table[["iC10 DNA+RNA"]][["metabric"]]$Sample[which(references_table[["iC10 DNA+RNA"]][["metabric"]]$ER.Status_LM == "-")])] = "ic4ER-"
        reference$iC10[which(reference$iC10 == "ic4" & reference$Sample %in% references_table[["iC10 DNA+RNA"]][["metabric"]]$Sample[which(references_table[["iC10 DNA+RNA"]][["metabric"]]$ER.Status_LM == "+")])] = "ic4ER+"
      }
    }else{
      predictions$voting[which(predictions$voting == "ic4" & predictions$Sample %in% er_status$Sample[which(er_status$ER == 0)])] = "ic4ER-"
      predictions$voting[which(predictions$voting == "ic4" & predictions$Sample %in% er_status$Sample[which(er_status$ER != 0)])] = "ic4ER+"
      
      reference$iC10[which(reference$iC10 == "ic4" & reference$Sample %in% er_status$Sample[which(er_status$ER == 0)])] = "ic4ER-"
      reference$iC10[which(reference$iC10 == "ic4" & reference$Sample %in% er_status$Sample[which(er_status$ER != 0)])] = "ic4ER+"
    }
    reference$iC10[which(reference$iC10 %in% c("ic3","ic7"))] = "Other"
    
    if(any(c(table(predictions$voting), table(reference$iC10)) < 4)){
      print(paste(set,ref))
      print(unique(c(names(table(predictions$voting))[which(table(predictions$voting) < 4)],
                     names(table(reference$iC10))[which(table(reference$iC10) < 4)])))
    }
    
    stopifnot(all(predictions$Sample %in% reference$Sample))
    reference = reference[which(reference$Sample %in% predictions$Sample),]
    reference = reference[order(match(reference$Sample, predictions$Sample)),]
    stopifnot(all(reference$Sample == predictions$Sample))
    
    predictions$type = predictions$voting
    predictions$type[which(predictions$voting %in% c("Other","ic4ER+","ic8"))] = "Typical"
    
    reference$type = reference$iC10
    reference$type[which(reference$iC10 %in% c("Other","ic4ER+","ic8"))] = "Typical"
    
    reference = reference[which(!reference$Sample %in% predictions$Sample[which(predictions$type == "ic4")]),]
    predictions = predictions[which(predictions$type != "ic4"),]
    
    # 1.0. F1 score for IC1, IC2, IC6, IC9, Typical (IC3, IC4ER+, IC7, IC8), IC5, IC10, IC4ER-
    
    mat = confusionMatrix(data=as.factor(predictions$type), reference = as.factor(reference$type))
    if(ref == "iC10 DNA+RNA"){
      balanced_accuracy_tmp = data.frame(set = set, ref= ref, version = file, 
                                         class = gsub("Class: ","",rownames(mat$byClass)), balanced_accuracy = mat$byClass[,"Balanced Accuracy"])
    }
    
    predictions$type[which(predictions$voting %in% c("ic10","ic4ER-"))] = "TNBC"
    reference$type[which(reference$iC10 %in% c("ic10","ic4ER-"))] = "TNBC"
    
    mat2 = confusionMatrix(data=as.factor(predictions$type), reference = as.factor(reference$type))
    
    performance_table_tmp[ref,grepl("F1",colnames(performance_table_tmp))] = sapply(gsub("IC","ic",colnames(performance_table_tmp)[grepl("F1",colnames(performance_table_tmp))]), function(x) if(!grepl("TNBC",x)){mat$byClass[which(endsWith(rownames(mat$byClass), gsub("F1 ","",x))),"F1"]}else{mat2$byClass[which(endsWith(rownames(mat2$byClass), gsub("F1 ","",x))),"F1"]})
    performance_table_tmp[ref,grepl("P ",colnames(performance_table_tmp))] = sapply(gsub("IC","ic",colnames(performance_table_tmp)[grepl("P ",colnames(performance_table_tmp))]), function(x) if(!grepl("TNBC",x)){mat$byClass[which(endsWith(rownames(mat$byClass), gsub("P ","",x))),"Precision"]}else{mat2$byClass[which(endsWith(rownames(mat2$byClass), gsub("P ","",x))),"Precision"]})
    performance_table_tmp[ref,startsWith(colnames(performance_table_tmp),"R ")] = sapply(gsub("IC","ic",colnames(performance_table_tmp)[startsWith(colnames(performance_table_tmp),"R ")]), function(x) if(!grepl("TNBC",x)){mat$byClass[which(endsWith(rownames(mat$byClass), gsub("R ","",x))),"Recall"]}else{mat2$byClass[which(endsWith(rownames(mat2$byClass), gsub("R ","",x))),"Recall"]})
    
    # 1.1. PR AUC for IC1, IC2, IC6, IC9
    auc = c()
    for(High_IC in paste0("ic",c(1,2,6,9))){
      data = data.frame(Sample = predictions$Sample, 
                        predictions = ifelse(predictions$voting == High_IC,High_IC,paste0("not-",High_IC)),
                        reference = ifelse(reference$iC10 == High_IC,High_IC,paste0("not-",High_IC)), 
                        IC = predictions[,paste0("voting_proba_", High_IC)], stringsAsFactors = F)
      data$reference = factor(data$reference)
      data$predictions = factor(data$predictions)
      
      auc_tmp = pr_auc(data, reference, IC)$.estimate
      auc = c(auc,auc_tmp)
    }
    performance_table_tmp[ref,grepl("PR AUC",colnames(performance_table_tmp))] <- auc
  }
  performance_table_tmp = as.data.frame(performance_table_tmp) %>% 
    tibble::rownames_to_column("reference") %>% 
    tidyr::pivot_longer(cols=!reference, names_to='metrics', values_to='value') %>% mutate(version = file)
  performance_table_tmp = as.data.frame(performance_table_tmp)
  performance_table_tmp_all = rbind(performance_table_tmp_all, performance_table_tmp)

  balanced_accuracy = rbind(balanced_accuracy, balanced_accuracy_tmp)
  performance_table[[length(performance_table) + 1]] <- performance_table_tmp_all
  names(performance_table)[[length(performance_table)]] = set
}

performance_table_iC10 = list()
for(set in c("validation","icgc","metabric")){
  
  col_metrics = c(paste0(rep(c("F1 ", "P ", "R "),each=9),c(paste0("IC", c(1,2,6,9)),"Typical",paste0("IC", c(5,10,"4ER-")), "TNBC")), paste0("PR AUC IC", c(1,2,6,9)))
  performance_table_tmp = matrix(NA, nrow = 1, ncol = length(col_metrics), dimnames = list("iC10 DNA+RNA", col_metrics))
  
  ref = "iC10 DNA+RNA"
  predictions = references_table[["iC10 DNA-only"]][[set]]
  reference = references_table[[ref]][[set]]
  
  er_path = ifelse(set == "validation", paste0(path_eniclust, file, "/Pipeline_2/validation_set.txt"), paste0(path_eniclust, file, "/Pipeline_2/",set,"_transformed_matrix.txt"))
  er_status = read.table(er_path,sep="\t",header=T)
  
  if(set == "validation"){
    predictions = predictions[which(grepl("TCGA",predictions$Sample)),]
    predictions = predictions[which(predictions$Sample %in% er_status$Sample),]
    stopifnot(nrow(predictions) == 203)
  }
  if(set %in% c("icgc","metabric")){
    predictions = predictions[which(predictions$Sample %in% reference$Sample),]
    er_status$Sample = sapply(er_status$Sample, function(x) if(set == "icgc"){str_sub(x,1,nchar("CURTIS_H000434"))}else{x})
  }
  
  if(set != "metabric"){
    predictions$iC10[which(predictions$iC10 == "ic4" & predictions$Sample %in% er_status$Sample[which(er_status$ER == 0)])] = "ic4ER-"
    predictions$iC10[which(predictions$iC10 == "ic4" & predictions$Sample %in% er_status$Sample[which(er_status$ER != 0)])] = "ic4ER+"
    
    reference$iC10[which(reference$iC10 == "ic4" & reference$Sample %in% er_status$Sample[which(er_status$ER == 0)])] = "ic4ER-"
    reference$iC10[which(reference$iC10 == "ic4" & reference$Sample %in% er_status$Sample[which(er_status$ER != 0)])] = "ic4ER+"
  }else{
    predictions$iC10[which(predictions$iC10 == "ic4" & predictions$Sample %in% references_table[["iC10 DNA+RNA"]][["metabric"]]$Sample[which(references_table[["iC10 DNA+RNA"]][["metabric"]]$ER.Status_LM == "-")])] = "ic4ER-"
    predictions$iC10[which(predictions$iC10 == "ic4" & predictions$Sample %in% references_table[["iC10 DNA+RNA"]][["metabric"]]$Sample[which(references_table[["iC10 DNA+RNA"]][["metabric"]]$ER.Status_LM == "+")])] = "ic4ER+"
    predictions = predictions[which(predictions$iC10 != "ic4"),]
  }
  if(any(c(table(predictions$voting), table(reference$iC10)) < 4)){
    print(paste(set,ref))
    print(unique(c(names(table(predictions$voting))[which(table(predictions$voting) < 4)],
                   names(table(reference$iC10))[which(table(reference$iC10) < 4)])))
  }
  
  stopifnot(all(predictions$Sample %in% reference$Sample))
  reference = reference[which(reference$Sample %in% predictions$Sample),]
  reference = reference[order(match(reference$Sample, predictions$Sample)),]
  stopifnot(all(reference$Sample == predictions$Sample))
  
  predictions$type = predictions$iC10
  predictions$type[which(predictions$iC10 %in% c("ic3","ic7","ic4ER+","ic8"))] = "Typical"
  
  reference$type = reference$iC10
  reference$type[which(reference$iC10 %in% c("ic3","ic7","ic4ER+","ic8"))] = "Typical"
  
  # 1.0. F1 score for IC1, IC2, IC6, IC9, Typical (IC3, IC4ER+, IC7, IC8), IC5, IC10, IC4ER-
  mat = confusionMatrix(data=as.factor(predictions$type), reference = as.factor(reference$type))
  
  predictions$type[which(predictions$iC10 %in% c("ic10","ic4ER-"))] = "TNBC"
  reference$type[which(reference$iC10 %in% c("ic10","ic4ER-"))] = "TNBC"
  
  mat2 = confusionMatrix(data=as.factor(predictions$type), reference = as.factor(reference$type))
  
  performance_table_tmp[ref,grepl("F1",colnames(performance_table_tmp))] = sapply(gsub("IC","ic",colnames(performance_table_tmp)[grepl("F1",colnames(performance_table_tmp))]), function(x) if(!grepl("TNBC",x)){mat$byClass[which(endsWith(rownames(mat$byClass), gsub("F1 ","",x))),"F1"]}else{mat2$byClass[which(endsWith(rownames(mat2$byClass), gsub("F1 ","",x))),"F1"]})
  performance_table_tmp[ref,grepl("P ",colnames(performance_table_tmp))] = sapply(gsub("IC","ic",colnames(performance_table_tmp)[grepl("P ",colnames(performance_table_tmp))]), function(x) if(!grepl("TNBC",x)){mat$byClass[which(endsWith(rownames(mat$byClass), gsub("P ","",x))),"Precision"]}else{mat2$byClass[which(endsWith(rownames(mat2$byClass), gsub("P ","",x))),"Precision"]})
  performance_table_tmp[ref,startsWith(colnames(performance_table_tmp),"R ")] = sapply(gsub("IC","ic",colnames(performance_table_tmp)[startsWith(colnames(performance_table_tmp),"R ")]), function(x) if(!grepl("TNBC",x)){mat$byClass[which(endsWith(rownames(mat$byClass), gsub("R ","",x))),"Recall"]}else{mat2$byClass[which(endsWith(rownames(mat2$byClass), gsub("R ","",x))),"Recall"]})
  
  # 1.1. PR AUC for IC1, IC2, IC6, IC9
  auc = c()
  for(High_IC in paste0("ic",c(1,2,6,9))){
    data = data.frame(Sample = predictions$Sample, 
                      predictions = ifelse(predictions$iC10 == High_IC,High_IC,paste0("not-",High_IC)),
                      reference = ifelse(reference$iC10 == High_IC,High_IC,paste0("not-",High_IC)), 
                      IC = predictions[,paste0("proba_", High_IC)], stringsAsFactors = F)
    data$reference = factor(data$reference)
    data$predictions = factor(data$predictions)
    
    auc_tmp = pr_auc(data, reference, IC)$.estimate
    auc = c(auc,auc_tmp)
  }
  performance_table_tmp[ref,grepl("PR AUC",colnames(performance_table_tmp))] <- auc
  
  performance_table_tmp = as.data.frame(performance_table_tmp) %>% 
    tibble::rownames_to_column("reference") %>% 
    tidyr::pivot_longer(cols=!reference, names_to='metrics', values_to='value') %>% mutate(version = "iC10 DNA-only")
  performance_table_tmp = as.data.frame(performance_table_tmp)
  
  performance_table_iC10[[length(performance_table_iC10) + 1]] <- performance_table_tmp
  names(performance_table_iC10)[[length(performance_table_iC10)]] = set
}

n_set = c()
p = data.frame()
for(set in c("validation","icgc","metabric")){
  if(set == "validation"){
    n = read.table(paste0(path_eniclust, file, "/Pipeline_2/",set,"_predictions.txt"),sep="\t",header=T) %>% dplyr::filter(grepl("TCGA",Sample)) %>% nrow()
  }
  if(set == "icgc"){
    mysamples = read.table(paste0(path_eniclust, file, "/Pipeline_2/",set,"_predictions.txt"),sep="\t",header=T) %>% mutate(Sample = str_sub(Sample,1,nchar("CURTIS_H001746"))) %>% dplyr::select(Sample) %>% pull()
    n = paste0(as.character(references_table[["iC10 DNA+RNA"]]$icgc %>% dplyr::filter(Sample %in% mysamples) %>% nrow()) , "-" , as.character(references_table[["iC10 DNA-only"]]$icgc_all %>% dplyr::filter(Sample %in% mysamples) %>% nrow()))
  }
  if(set == "metabric"){
    mysamples = read.table(paste0(path_eniclust, file, "/Pipeline_2/",set,"_predictions.txt"),sep="\t",header=T) %>% dplyr::filter(Sample %in% references_table[["iC10 DNA+RNA"]]$metabric$Sample) %>% pull(Sample)
    n = paste0(as.character(length(mysamples)) , "-" , as.character(length(which(!mysamples %in% gsub("[.]","-",colnames(train.CN))))))
  }
  n_set = c(n_set, n)
  p_tmp <- performance_table[[set]] %>% 
    bind_rows(performance_table_iC10[[set]]) %>%
    rowwise() %>% mutate(type = ifelse(sapply(strsplit(metrics," "),"[[",1) == "PR","PR AUC",sapply(strsplit(metrics," "),"[[",1)), 
                         metrics = ifelse(grepl("Typical",metrics),"Typical",ifelse(grepl("TNBC",metrics),"TNBC",paste0("IC",sapply(strsplit(metrics,"IC"),"[[",2)))))
  p_tmp$set = paste0(toupper(set)," n=",n)
  p <- bind_rows(p, p_tmp)
  
}
p$flag = ">3" 
p$flag[which(p$version == file & grepl("ICGC",p$set) & p$reference == "iC10 DNA+RNA" & grepl("IC4ER-|IC5",p$metrics))] = "<4"
p$flag[which(p$version == file & grepl("ICGC",p$set) & p$reference == "iC10 DNA-only" & grepl("IC4ER-|IC5",p$metrics))] = "<4"
p$flag[which(p$version == "iC10 DNA-only" & grepl("ICGC",p$set) & p$reference == "iC10 DNA+RNA" & grepl("IC4ER-",p$metrics))] = "<4"

### SAVE ##########################################################################################
## Figures----
pdf(paste0(outdir,"SuppFigure1d.pdf"), w=11.92, h=5.93)
p %>%
  dplyr::mutate(set = factor(set, levels=c("VALIDATION n=203","ICGC n=189-641","METABRIC n=1894-937"))) %>%
  ggplot(aes(x=factor(ifelse(version == "iC10 DNA-only","iC10","ENiClust"), levels = c("iC10","ENiClust")), 
             y=factor(metrics,levels=rev(c(paste0("IC",c(1,2,6,9)),"Typical","TNBC",paste0("IC",c(5,10,"4ER-"))))), 
             fill = value)) +
  geom_tile(aes(color = flag), linewidth = 1) + 
  scale_fill_gradient(low="white",high="#23005c", name = "Value", limits = c(0.15,1)) +
  geom_text(aes(label = as.character(round(value,digits = 2)), col = value < 0.5), size=3) +
  scale_color_manual(values=c(">3"="transparent","<4"="red","FALSE"="white","TRUE"="black")) +
  facet_grid(factor(type,levels=c("PR AUC","F1","P","R")) ~ 
               set + 
               factor(reference,levels=c("iC10 DNA+RNA","iC10 DNA-only","iC10 DNA-only (All)")), scales = "free", space='free') +
  xlab('') + ylab("Score") +
  theme_LM + 
  guides(color="none")
dev.off()

## SourceData----
performance_table$validation$set = "validation"
performance_table$validation$version = "ENiClust"
performance_table$icgc$set = "icgc"
performance_table$validation$version = "ENiClust"
performance_table$metabric$set = "metabric"
performance_table$validation$version = "ENiClust"

performance_table_iC10$validation$set = "validation"
performance_table_iC10$icgc$set = "icgc"
performance_table_iC10$metabric$set = "metabric"
sourcetable_ENiClust = rbind(rbind(performance_table$validation,performance_table$icgc),performance_table$metabric)
sourcetable_iC10 = rbind(rbind(performance_table_iC10$validation,performance_table_iC10$icgc),performance_table_iC10$metabric)
sourcetable = rbind(sourcetable_ENiClust, sourcetable_iC10)
colnames(sourcetable)[4] = "model"
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/github/SupplementaryFigure1d_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

balanced_accuracy$version = "ENiClust"
colnames(balanced_accuracy)[3] = "model"
write.table(balanced_accuracy, paste0(basedir,"/lisem/bc_landscape/github/SupplementaryFigure1d_balanced_accuracy.txt"), row.names = F, col.names = T, quote = F, sep="\t")

## For Main Text----
#ENiClust yielded high balanced accuracy for the four ER+ High-risk ICs (IC1, IC2, IC6, and IC9) subtypes 
#(TCGA-WES Validation: 0.77-0.91; Nik-Zainal-WGS Validation: 0.72-0.99; METABRIC-ASCAT: 0.74-0.87), 
# as well as IC5 (HER2+, TCGA-WES Validation: 0.95; Nik-Zainal-WGS Validation: 0.63; METABRIC-ASCAT: 0.91) 
# and IC10 (TNBC, TCGA-WES Validation: 0.88; Nik-Zainal-WGS Validation: 0.89; METABRIC-ASCAT: 0.84)
# TCGA-Validation
range(balanced_accuracy[which(balanced_accuracy$set == "validation" & balanced_accuracy$class %in% c("ic1","ic2","ic6","ic9")),"balanced_accuracy"]) # 0.77-0.91
balanced_accuracy[which(balanced_accuracy$set == "validation" & balanced_accuracy$class == "ic5"),] # 0.95
balanced_accuracy[which(balanced_accuracy$set == "validation" & balanced_accuracy$class == "ic10"),] # 0.88

# ICGC-Replication
range(balanced_accuracy$balanced_accuracy[which(balanced_accuracy$set == "icgc" & balanced_accuracy$class %in% c("ic1","ic2","ic6","ic9"))]) #0.72-0.99
balanced_accuracy$balanced_accuracy[which(balanced_accuracy$set == "icgc" & balanced_accuracy$class == "ic5")] #0.63
balanced_accuracy$balanced_accuracy[which(balanced_accuracy$set == "icgc" & balanced_accuracy$class == "ic10")] #0.89

# METABRIC-Replication
range(balanced_accuracy$balanced_accuracy[which(balanced_accuracy$set == "metabric" & balanced_accuracy$class %in% c("ic1","ic2","ic6","ic9"))]) #0.74-0.87
balanced_accuracy$balanced_accuracy[which(balanced_accuracy$set == "metabric" & balanced_accuracy$class == "ic5")] #0.91
balanced_accuracy$balanced_accuracy[which(balanced_accuracy$set == "metabric" & balanced_accuracy$class == "ic10")] #0.84

