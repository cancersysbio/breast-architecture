### CREATE SUPPLEMENTARY FIGURE 1h #############################################################################
# creates supplementary figure 1h
# supplementary figure 1h provides the PAM50 proliferation score from mRNA in ER+ High-risk and ER+ Typical-risk cases at pre-treatment and early on-treatment time points.

### PREAMBLE #####################################################################################
library(yaml)
library(scales)
library(dplyr)
library(tidyr)
library(lemon, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
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

PAM50_PROLIF = c("BIRC5","CCNB1","CDC20","NUF2","CEP55","KNTC2","MKI67","PTTG1","RRM2","TYMS","UBE2C")

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'
outdir <- paste0(basedir,"lisem/bc_landscape/revision/figures/")

load(file = paste0(basedir,"lisem/Selli/data/selli_merged_dataset.RData"))
load(file = paste0(basedir,"lisem/Selli/data/selli_exp_norm_genes_all.RData"))

mydata$treatment[which(endsWith(mydata$Sample,"-3"))] = "early"
mydata$PAM50_PROLIF <- sapply(mydata$Sample, function(x) mean(unlist(y_genes[which(rownames(y_genes) %in% PAM50_PROLIF), which(grepl(x,colnames(y_genes)))])))
mydata$ic_group <- sapply(mydata$Patient, function(x) if(any(mydata$Patient == x & mydata$treatment == "pre")){ifelse(mydata$iC10[which(mydata$Patient == x & mydata$treatment == "pre")] %in% paste0("ic",c(1,2,6,9)), "High-risk",ifelse(mydata$iC10[which(mydata$Patient == x & mydata$treatment == "pre")] %in% paste0("ic",c(3,4,7,8)), "Typical-risk", NA))}else{NA})

p <- mydata %>%
  dplyr::filter(treatment %in% c("pre","early") & !is.na(ic_group)) %>%
  group_by(Patient, treatment, ic_group) %>% summarize(PAM50_PROLIF = mean(PAM50_PROLIF)) %>%
  ggplot(aes(x = factor(treatment, levels = c("pre","early")), y = PAM50_PROLIF, col = ic_group)) + 
  geom_boxplot(aes(col = ic_group), width = 0.5) +
  lemon::geom_pointpath(aes(colour=ic_group, group=Patient), linesize = 1, size=4, alpha=0.7) + # position
  scale_color_manual(name = "", values = c("High-risk" = "#ff8400", "Typical-risk" = "#529ed4")) +
  facet_grid(. ~ ic_group) + ylim(4.95, 7) +
  xlab("") + ylab("Proliferation score") + scale_x_discrete( #expand = c(0.6,1.5), 
    labels=c('pre' = "Pre-tx",'early' = "Early-tx")) +
  theme_LM + ggtitle("Proliferation level (mRNA) from Pre-tx to Early-tx") +
  theme(legend.position = "none", axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14),
        title = element_text(size=16), 
        strip.text = element_text(size=16))

### SAVE ##########################################################################################
## Figures----
png(filename = file.path(outdir, 'SupplementaryFigure_1h.png'), res = 300, width = 6.53, height = 6.02, units = 'in')
p
dev.off()

## SourceData----
sourcedata <- mydata %>% dplyr::filter(treatment %in% c("pre","early") & !is.na(ic_group))
write.table(sourcedata, paste0(basedir,"lisem/bc_landscape/github/SupplementaryFigure1h_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")

