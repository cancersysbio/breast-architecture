### CREATE EXTENDED DATA FIGURE 3m #############################################################################
# creates extended data figure 3m
# extended data figure 3m provides differential pattern of relapse across ER+ IC subtypes and by histology (IDC: invasive ductal carcinoma and ILC: invasive lobular carcinoma), illustrated by the cumulative (black) and annual (red) risk of relapse.

### PREAMBLE #####################################################################################
library(yaml)
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
basedir <- "/oak/stanford/groups/ccurtis2/users/"
path.git <- paste0(basedir,"lisem/git_repo/brcarepred/")
outdir <- paste0(basedir,"lisem/bc_landscape/submission/figures/")

### Load Cox Model data----
load(file=paste0(path.git,"ILC/IntClust_ILC_AllProbs.RData")) # internally called IntClust_AllProbs.RData

r <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psdo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psld)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psl)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psc)

r0 <- do.call("cbind", r)
rownames(r0) <- psc[[1]][,1]

clinic = read.table(paste0(basedir,"lisem/bc_landscape/github/clinic_METABRIC.txt"), header = T, sep = "\t") # see: clinic_METABRIC.R
clinic$eniclust[which(clinic$eniclust != "Other")] = sapply(strsplit(clinic$eniclust[which(clinic$eniclust != "Other")],"ic"),"[[",2)
table(clinic[,c("eniclust", "Histological.Type")])

clinic$group = NA
clinic$group[which(clinic$Histological.Type == "ILC" & clinic$eniclust %in% c("3","7","8","Other","4ER+"))] = "ILC-L"
clinic$group[which(clinic$Histological.Type == "ILC" & clinic$eniclust %in% c("1","2","6","9"))] = "ILC-H"

clinic$group[which(grepl("IDC",clinic$Histological.Type) & clinic$Histological.Type != "IDC+ILC" & clinic$eniclust %in% c("3","7","8","Other","4ER+"))] = "IDC-L"
clinic$group[which(grepl("IDC",clinic$Histological.Type) & clinic$Histological.Type != "IDC+ILC" & clinic$eniclust %in% c("1","2","6","9"))] = "IDC-H"
clinic$group[which(grepl("IDC",clinic$Histological.Type) & clinic$Histological.Type != "IDC+ILC" & clinic$eniclust == "4ER-")] = "IDC-4-"
clinic$group[which(grepl("IDC",clinic$Histological.Type) & clinic$Histological.Type != "IDC+ILC" & clinic$eniclust == "5")] = "IDC-5"
clinic$group[which(grepl("IDC",clinic$Histological.Type) & clinic$Histological.Type != "IDC+ILC" & clinic$eniclust == "10")] = "IDC-10"
table(clinic$group)

clinic = clinic[which(!is.na(clinic$group)),]

unique(clinic$group)
clinic$group <- factor(clinic$group, levels=c("IDC-L","ILC-L","IDC-H","ILC-H","IDC-4-","IDC-5","IDC-10"))
clinic <- clinic[which(clinic$METABRIC.ID %in% colnames(r0)),]
clinic <- clinic[match(colnames(r0), clinic$METABRIC.ID),]
mean(clinic$METABRIC.ID==colnames(r0))

r0.INTCLUST.mean <- t(apply(r0, 1, function(y) tapply(y, clinic$group, mean)))
r0.INTCLUST.SE <- t(apply(r0, 1, function(y) tapply(y, clinic$group, function(z) sd(z)/sqrt(length(z)))))

ss <- paste0("n=", sapply(r, function(x) ncol(x)))
res <- data.frame(X=c(r0.INTCLUST.mean), Year=rep(rownames(r0.INTCLUST.mean), length(levels(clinic$group))),
                  Predictions=rep("After Surgery", length(levels(clinic$group))*nrow(r0.INTCLUST.mean)),
                  IntClust=rep(colnames(r0.INTCLUST.mean), rep(nrow(r0.INTCLUST.mean), length(levels(clinic$group)))),
                  li=c(r0.INTCLUST.mean) - 1.96*c(r0.INTCLUST.SE),
                  ui=c(r0.INTCLUST.mean) + 1.96*c(r0.INTCLUST.SE))

names(ss) <- levels(res$IntClust)
res$pch <- 17
res$colr <- "black"
ss <- ss[levels(res$IntClust)]

res$Year <- as.numeric(as.character(res$Year))
data = res
data$Xbis <-as.numeric(sapply(unique(data$IntClust), function(ic) sapply(1:5, function(i) if(i == 1){data$X[which(data$IntClust == ic)][i]}else{data$X[which(data$IntClust == ic)][i] - data$X[which(data$IntClust == ic)][i-1]})))
data$libis <-as.numeric(sapply(unique(data$IntClust), function(ic) sapply(1:5, function(i) if(i == 1){data$li[which(data$IntClust == ic)][i]}else{data$li[which(data$IntClust == ic)][i] - data$li[which(data$IntClust == ic)][i-1]})))
data$uibis <-as.numeric(sapply(unique(data$IntClust), function(ic) sapply(1:5, function(i) if(i == 1){data$ui[which(data$IntClust == ic)][i]}else{data$ui[which(data$IntClust == ic)][i] - data$ui[which(data$IntClust == ic)][i-1]})))

data <- data.frame(X = c(data$X, data$Xbis),
                   li = c(data$li, data$libis),
                   ui = c(data$ui, data$uibis), 
                   Year = c(data$Year, data$Year), 
                   Predictions = c(rep("Cumulative", nrow(data)), rep("Annual", nrow(data))),
                   IntClust = c(data$IntClust, data$IntClust), stringsAsFactors = F)
data$group <- NA
data$group[which(grepl("-H",data$IntClust))] <- "IC1/IC2/IC6/IC9"
data$group[which(grepl("-L",data$IntClust))] <- "IC3/IC4/IC7/IC8"

p <- data %>%
  dplyr::filter(IntClust %in% c("IDC-L","ILC-L","IDC-H","ILC-H")) %>%
  dplyr::rowwise() %>% dplyr::mutate(IntClust = ifelse(grepl("IDC",IntClust),ifelse(grepl("-H",IntClust),"ER+ High - IDC","ER+ Typical - IDC"),
                                                       ifelse(grepl("-H",IntClust),"ER+ High - ILC","ER+ Typical - ILC"))) %>%
  ggplot(aes(x=Year, y=X, fill=Predictions, col=Predictions)) +
  geom_point() +
  geom_line() +
  scale_color_manual(name = "Risk", values = c("Annual" = "red", "Cumulative" = "black")) +
  scale_fill_manual(name = "Risk", values = c("Annual" = "red", "Cumulative" = "black")) +
  facet_wrap(factor(IntClust, levels = c("ER+ High - IDC","ER+ High - ILC","ER+ Typical - IDC","ER+ Typical - ILC"))~group, ncol=2, nrow = 2) + ylab('Probability of relapse') +
  theme_LM + theme(axis.text = element_text(size=14),
                   axis.title = element_text(size=18),
                   strip.text = element_text(size=20), 
                   legend.text = element_text(size=14),
                   legend.title = element_text(size=18))

### SAVE ##########################################################################################
## Figure----
png(filename = file.path(outdir, 'Figure2e.png'), res = 300, width=5.21*2, height=1.91*5, units = 'in')
p
dev.off()

## SourceData----
sourcetable = data
unique(sourcetable$IntClust)
sourcetable$IntClust[which(grepl("-L",sourcetable$IntClust))] = paste0(sapply(strsplit(sourcetable$IntClust[which(grepl("-L",sourcetable$IntClust))],"-"),"[[",1)," - ER+ Typical")
sourcetable$IntClust[which(grepl("-H",sourcetable$IntClust))] = paste0(sapply(strsplit(sourcetable$IntClust[which(grepl("-H",sourcetable$IntClust))],"-"),"[[",1)," - ER+ High")
sourcetable$IntClust[which(sourcetable$IntClust == "IDC-4-" )] = "IDC - IC4ER-"
sourcetable$IntClust[which(sourcetable$IntClust == "IDC-10" )] = "IDC - IC10"
sourcetable$IntClust[which(sourcetable$IntClust == "IDC-5" )] = "IDC - IC5"

colnames(sourcetable) = c("Probability","lower","upper","Year","Predictions","Subgroup","ENiClust")
sourcetable$ENiClust[which(is.na(sourcetable$ENiClust))] = sapply(strsplit(sourcetable$Subgroup[which(is.na(sourcetable$ENiClust))],"- "),"[[",2)
write.table(sourcetable, "/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/submission/data/Figure2e_sourcetable.txt", row.names = F, col.names = T, quote = F, sep="\t")

## Main Text----
sourcetable[which(sourcetable$Year == 5 & sourcetable$ENiClust == "IC1/IC2/IC6/IC9" & sourcetable$Predictions == "Cumulative"),]
sourcetable[which(sourcetable$Year == 20 & sourcetable$ENiClust == "IC1/IC2/IC6/IC9" & sourcetable$Predictions == "Cumulative"),]

sourcetable[which(sourcetable$Year == 5 & sourcetable$ENiClust == "IC3/IC4/IC7/IC8" & sourcetable$Predictions == "Cumulative"),]
sourcetable[which(sourcetable$Year == 20 & sourcetable$ENiClust == "IC3/IC4/IC7/IC8" & sourcetable$Predictions == "Cumulative"),]
