#####################################################
### Preparing clinical table for METABRIC samples ###
#####################################################

### Set paths----
basedir <- "/oak/stanford/groups/ccurtis2/users/"
path.git <- paste0(basedir,"lisem/git_repo/brcarepred/")

## see: https://github.com/cclab-brca/brcarepred/blob/master/Pipeline/4.CreateModelIntCluster.R (from Rueda et al. Nature 2019)
clinic = read.table(file=paste0(path.git,"Tables/TableS6.txt"), header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)

## We remove Samples with no follow-up time
ids <- which(clinic$T==0)
if (length(ids)>0) clinic <- clinic[-ids,]

## We remove samples with stage 4
clinic <- clinic[-which(clinic$Stage==4),]

## We remove benign, DCIS or PHYL
unique(clinic$Histological.Type)
bad.hist <- which(clinic$Histological.Type %in% c("BENIGN", "PHYL", "DCIS"))
if (length(bad.hist)>0) clinic <- clinic[-bad.hist,]

clinic <- clinic[which(!is.na(clinic$T)),]
clinic <- clinic[which(!is.na(clinic$Death)),]

## We move a bit relapses from diagnosis
clinic$TLR[which(clinic$TLR==0)] <- 0.1 
clinic$TDR[which(clinic$TDR==0)] <- 0.1

clinic$TR = base::pmin(clinic$TDR,clinic$TLR)
clinic$R = ifelse(clinic$LR | clinic$DR,1,0)
clinic$T[which(clinic$T==clinic$TDR & clinic$DR==1)] <- clinic$T[which(clinic$T==clinic$TDR & clinic$DR==1)] + 0.1
clinic$T[which(clinic$T==clinic$TLR & clinic$LR==1)] <- clinic$T[which(clinic$T==clinic$TLR & clinic$LR==1)] + 0.1
clinic$TDR[which(clinic$TLR==clinic$TDR & clinic$LR==1 & clinic$DR==1)] <- clinic$TDR[which(clinic$TLR==clinic$TDR & clinic$LR==1 & clinic$DR==1)] + 0.1
clinic$TLR[which(clinic$LR==0)] <- clinic$T[which(clinic$LR==0)]
clinic$TDR[which(clinic$DR==0)] <- clinic$T[which(clinic$DR==0)]

clinic$LR[which(clinic$TLR>clinic$TDR)] <- 0 # If local relapse occurred after distant, we don't consider it
clinic$TLR[which(clinic$TLR>clinic$TDR)] <- clinic$T[which(clinic$TLR>clinic$TDR)]

clinic$T <- clinic$T/365.25
clinic$TLR <- clinic$TLR/365.25
clinic$TDR <- clinic$TDR/365.25
clinic$TR <- clinic$TR/365.25

clinic$T_20 = ifelse(clinic$T>=20, 20, clinic$T) #time stops at 20 years
clinic$DeathBreast_20 = ifelse(clinic$T>=20, 0, clinic$DeathBreast) #death stops after 20 years

clinic$Death_20 = ifelse(clinic$T>=20, 0, clinic$Death)

clinic$TLR_20 = ifelse(clinic$TLR>=20, 20, clinic$TLR)
clinic$LR_20 = ifelse(clinic$TLR>=20, 0, clinic$LR)

clinic$TDR_20 = ifelse(clinic$TDR>=20, 20, clinic$TDR)
clinic$DR_20 = ifelse(clinic$TDR>=20, 0, clinic$DR)

clinic$TR_20 = ifelse(clinic$TR>=20, 20, clinic$TR)
clinic$R_20 = ifelse(clinic$TR>=20, 0, clinic$R)

clinic$LN <- clinic$Lymph.Nodes.Positive
clinic$LN[which(clinic$LN>=10)] <- 10
clinic$AGE <- clinic$Age.At.Diagnosis
clinic$Grade[clinic$Grade=="null"]=NA
clinic$GRADE <- as.factor(as.character(clinic$Grade))
clinic$SIZE <- as.numeric(as.character(clinic$Size))
clinic$Stage[clinic$Stage=="null"]=NA
clinic$STAGE = factor(clinic$Stage,labels = c("0","I","II","III"))

clinic$HT <- 1 * (clinic$HT!="NO/NA")
clinic$HT <- factor(clinic$HT, levels=c(0, 1), labels=c("NO", "YES")) 
clinic$CT <- 1 * (clinic$CT!="NO/NA")
clinic$CT <- factor(clinic$CT, levels=c(0, 1), labels=c("NO", "YES"))
clinic$RT[which(clinic$RT %in% c("NO/NA", "NONE RECORDED IN LANTIS", "Nnne"))] <- 0
clinic$RT[which(clinic$RT!=0)] <- 1
clinic$RT <- factor(clinic$RT, levels=c(0, 1), labels=c("NO", "YES"))
clinic$BS <- factor(clinic$Breast.Surgery, levels=c("BREAST CONSERVING", "MASTECTOMY"), labels=c("BC", "M"))

## Add iC10 DNA-only and RNA-only calls
clinic$ER.Status_LM = ifelse(is.na(clinic$ER.Status), clinic$ER.Expr, clinic$ER.Status)
clinic$ER.Status_LM[which(clinic$ER.Status_LM == "pos")] = "+"
clinic$ER.Status_LM[which(clinic$ER.Status_LM == "neg")] = "-"

metabric_iC10_DNA = read.table(paste0(basedir,"lisem/bc_landscape/ic_subtypes/data/metabric_iC10_CNA-only_labels.txt"), header = T, sep = "")
metabric_iC10_RNA = read.table(paste0(basedir,"lisem/bc_landscape/ic_subtypes/data/metabric_iC10_RNA-only_labels.txt"), header = T, sep = "")

clinic$iC10_DNA = sapply(clinic$METABRIC.ID, function(x) metabric_iC10_DNA$iC10[which(metabric_iC10_DNA$Sample == x)])
clinic$iC10_DNA[which(clinic$iC10_DNA == "ic4" & clinic$ER.Status_LM == "+")] = "ic4ER+"
clinic$iC10_DNA[which(clinic$iC10_DNA == "ic4" & clinic$ER.Status_LM == "-")] = "ic4ER-"

clinic$iC10_RNA = sapply(clinic$METABRIC.ID, function(x) metabric_iC10_RNA$iC10[which(metabric_iC10_RNA$Sample == x)])
clinic$iC10_RNA[which(clinic$iC10_RNA == "ic4" & clinic$ER.Status_LM == "+")] = "ic4ER+"
clinic$iC10_RNA[which(clinic$iC10_RNA == "ic4" & clinic$ER.Status_LM == "-")] = "ic4ER-"

# ENiClust
predictions = read.table(paste0(basedir, "lisem/bc_landscape/ic_subtypes/data/eniclust/metabric_predictions.txt"),sep="\t",header=T)
clinic$eniclust = sapply(clinic$METABRIC.ID, function(x) if(x %in% predictions$Sample){predictions$voting[which(predictions$Sample == x)]}else{NA})
clinic$eniclust[which(clinic$eniclust == "ic4")] = paste0("ic4ER",ifelse(is.na(clinic$ER.Status[which(clinic$eniclust == "ic4")]),clinic$ER.Expr[which(clinic$eniclust == "ic4")],ifelse(clinic$ER.Status[which(clinic$eniclust == "ic4")] == "pos","+","-")) )

## NB: Covariates to consider: AGE+GRADE+SIZE+LN+ER.Expr
write.table(clinic, paste0(basedir,"lisem/bc_landscape/github/clinic_METABRIC.txt"), col.names = T, row.names = F, sep = "\t")
