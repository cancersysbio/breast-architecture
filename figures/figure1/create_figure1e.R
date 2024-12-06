### CREATE FIGURE 1e #############################################################################
# creates figure 1e
# figure 1e provides differential pattern of relapse across IC subtypes and groups, illustrated by the cumulative (black) and annual (red) risk of relapse over time.

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
load(file=paste0(path.git,"Figure1/IntClust_AllProbs.RData"))

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

clinic = clinic[which(!is.na(clinic$eniclust)),]

clinic$eniclust <- factor(clinic$eniclust, levels=c(1:2, "4ER+", "4ER-", 5:6, 9:10, "Other"))
clinic <- clinic[which(clinic$METABRIC.ID %in% colnames(r0)),]
clinic <- clinic[match(colnames(r0), clinic$METABRIC.ID),]

r0.INTCLUST.mean <- t(apply(r0, 1, function(y) tapply(y, clinic$eniclust, mean)))
r0.INTCLUST.SE <- t(apply(r0, 1, function(y) tapply(y, clinic$eniclust, function(z) sd(z)/sqrt(length(z)))))

ss <- paste0("n=", sapply(r, function(x) ncol(x)))
res <- data.frame(X=c(r0.INTCLUST.mean), Year=rep(rownames(r0.INTCLUST.mean), length(levels(clinic$eniclust))),
                  Predictions=rep("After Surgery", length(levels(clinic$eniclust))*nrow(r0.INTCLUST.mean)),
                  IntClust=rep(colnames(r0.INTCLUST.mean), rep(nrow(r0.INTCLUST.mean), length(levels(clinic$eniclust)))),
                  li=c(r0.INTCLUST.mean) - 1.96*c(r0.INTCLUST.SE),
                  ui=c(r0.INTCLUST.mean) + 1.96*c(r0.INTCLUST.SE))

res$IntClust <- paste0("IC", res$IntClust)
res$IntClust[which(res$IntClust == "ICOther")] <- "IC3/IC7/IC8"

res$IntClust <- factor(res$IntClust, levels=c("IC1","IC2","IC4ER+","IC4ER-","IC5","IC6","IC9","IC10","IC3/IC7/IC8"))
names(ss) <- levels(res$IntClust)
res$pch <- 17
res$colr <- "black"

ss <- ss[levels(res$IntClust)]
res$Year <- as.numeric(as.character(res$Year))

data = res[which(res$Predictions == "After Surgery"),]
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
data$group[which(data$IntClust %in% c("IC3/IC7/IC8","IC4ER+"))] <- "ER+ Typical"
data$group[which(data$IntClust %in% c("IC1","IC2","IC6","IC9"))] <- "ER+ High"
data$group[which(data$IntClust %in% c("IC10","IC4ER-"))] <- "TNBC"
data$group[which(data$IntClust == "IC5")] <- "HER2+"
data$group <- factor(data$group, levels = c("TNBC","ER+ High","ER+ Typical","HER2+"))

p <- ggplot(data, aes(x=Year, y=X, fill=Predictions, col=Predictions)) +
  geom_point() +
  geom_line() +
  scale_color_manual(name = "Risk", values = c("Annual" = "red", "Cumulative" = "black")) +
  scale_fill_manual(name = "Risk", values = c("Annual" = "red", "Cumulative" = "black")) +
  facet_wrap(group~IntClust, ncol=2, nrow = 5) + ylab('Probability of relapse') +
  theme_LM + theme(axis.text = element_text(size=14),
                   axis.title = element_text(size=18),
                   strip.text = element_text(size=20), 
                   legend.text = element_text(size=14),
                   legend.title = element_text(size=18))

### SAVE ##########################################################################################
## Figure----
png(filename = file.path(outdir, 'Figure1e.png'), res = 300, width=5.21*2, height=1.91*5, units = 'in')
p
dev.off()

## SourceData----
sourcetable = data
colnames(sourcetable) = c("Probability","lower","upper","Year","Predictions","ENiClust","Subgroup")
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/submission/data/Figure1e_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
