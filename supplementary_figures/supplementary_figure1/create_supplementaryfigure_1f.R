### CREATE SUPPLEMENTARY FIGURE 1f #############################################################################
# creates supplementary figure 1f
# supplementary figure 1f provides the forest plot of the hazard ratios of the ER+ Typical vs. ER+ High-risk classification from the four IC subgroups classifiers (ENiClust, iC10 DNA+RNA, iC10 DNA-only, and iC10 RNA-only) for DRF survival.

### PREAMBLE #####################################################################################
library(yaml)
require(survival)
require(survminer)
require(patchwork)
require(scales)
library(formattable, lib.loc = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.2.2-Seurat/")

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

ggforestplot_LM <- function(data, range=NULL, v_line=NULL, stat_lable=NULL, colour_lines = colour_lines, pt_shape=19, summary_shape=23, pad=0, x_step_delta=1, x_lab_pos = NULL){
  
  dat <- data
  
  # Set bounds for each horizontal line segment
  dat$x1 <- dat$lower
  dat$x2 <- dat$upper
  
  dat$arrow_low<-FALSE
  dat$arrow_high<-FALSE
  
  if (!is.null(range)){
    dat$arrow_low[dat$lower<range[1]]<-TRUE
    dat$x1[dat$arrow_low==TRUE]<-range[1]
    dat$arrow_high[dat$upper>range[2]]<-TRUE
    dat$x2[dat$arrow_high==TRUE]<-range[2]
  }
  
  n_dat <- dim(dat)[1] + 2 #This will be the "height" of the graph
  
  #Initiate the plot
  g <-ggplot() +
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x	= element_line(size=1),
          axis.title.y = element_blank(),
          axis.text.x=element_text(size=12)
    )
  
  x_range <- range[2]-range[1]
  x_delta <- x_range/4      
  
  x_lab_start <-range[1]-1*x_delta-pad
  x_lab_end <-range[2]-1*x_delta-pad
  
  l <- dim(dat)[1]
  y_pos <- seq((n_dat-1),(n_dat-l),-1)
  dat$y_pos <- y_pos
  
  # Add test lables to the left of the data points
  temp_labs <- dat[,"Group"]
  pval<- dat[,"pval"]
  g <- g + ggplot2::annotate("text",x=x_lab_start,y=y_pos,label=temp_labs, hjust=0)
  g <- g + ggplot2::annotate("text",x=range[2]+pad,y=y_pos,label=ifelse(pval < 0.01,paste0("italic(p) == ",scientific(pval)), paste0("italic(p) == ",digits(pval,3))), hjust=0,parse = TRUE) + coord_cartesian(xlim = c(x_lab_start, range[2]+x_delta+pad), clip = 'off')
  
  # Apply points and lines
  for (j in 1:l) {
    dd <- dat[j,c("x1","x2","arrow_low","arrow_high","y_pos")]
    if(dd$arrow_low==TRUE & dd$arrow_high==TRUE){
      g <- g + geom_segment(data=dd, aes(x=x1,xend=x2,y=y_pos, yend=y_pos), arrow=arrow(length = unit(0.125, "inches"), ends="both"),colour=colour_lines)}
    if(dd$arrow_low==TRUE & dd$arrow_high==FALSE){
      g <- g + geom_segment(data=dd, aes(x=x1, xend=x2,y=y_pos, yend=y_pos), arrow=arrow(length = unit(0.125, "inches"), ends="first"),colour=colour_lines)}
    if(dd$arrow_low==FALSE & dd$arrow_high==TRUE){
      g <- g + geom_segment(data=dd, aes(x=x1,xend=x2,y=y_pos, yend=y_pos), arrow=arrow(length = unit(0.125, "inches"), ends="last"),colour=colour_lines)}
    if(dd$arrow_low==FALSE & dd$arrow_high==FALSE){
      g <- g+ geom_segment(data=dd, aes(x=x1,xend=x2,y=y_pos, yend=y_pos), colour=colour_lines)}
  }
  g <- g + geom_point(data=dat, aes(x=test,y=y_pos, size=n_dat),shape=pt_shape)
  g <- g + geom_text(data=dat, aes(x=test,y=y_pos-0.3,label=digits(test,2)), hjust=0.5, size = 3)
  g <- g + geom_text(data=dat, aes(x=test,y=y_pos-0.5,label=paste0("(",digits(lower,2),"-",digits(upper,2),")")), hjust=0.5, size = 3)
  
  # Horizontal line
  if(!is.null(v_line)){
    g <- g + geom_vline(xintercept=v_line, linetype="dashed")
  }
  g <- g + xlim(x_lab_start[1],range[2]) + xlab(stat_lable)
  g <- g + xlim(x_lab_start[1],range[2]) + scale_x_continuous(breaks=seq(range[1],range[2],x_step_delta)) + xlab(stat_lable)
  w <- (range[2]-x_lab_pos)/(range[2]-x_lab_start)
  g <- g + theme(axis.title.x=element_text(size=14, face="bold", hjust=1-w, margin = ggplot2::margin(t=7, r=0, b=0, l=0)))
  
  return(g)
}

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'
outdir <- "/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/submission/figures/"

### Combined forest plot
load(file = paste0(basedir, "lisem/bc_landscape/submission/data/SuppFigure1f_fit.coxph.RData")) # see: create_figure1c.R

Group <- names(fit.coxph)
data <- as.data.frame(Group)
for(i in 1:length(Group)){
  group = Group[i]
  fit <- fit.coxph[group]
  data$test[i] = as.numeric(summary(fit[[1]])$coefficients["ICsHigh-risk","exp(coef)"])
  data$pval[i] = as.numeric(summary(fit[[1]])$coefficients["ICsHigh-risk","Pr(>|z|)"])
  data$lower[i] = as.numeric(summary(fit[[1]])$conf.int["ICsHigh-risk","lower .95"])
  data$upper[i] = as.numeric(summary(fit[[1]])$conf.int["ICsHigh-risk","upper .95"])
}
data$qval = data$pval

p <- ggforestplot_LM(data, range=c(0.9,2.5), v_line=1 , #Group
                     stat_lable="Hazard ratio", 
                     colour_lines = "black", pt_shape=15, summary_shape=9 ,pad = 0.1, 
                     x_step_delta=0.3, x_lab_pos =1)

### SAVE ##########################################################################################
## Figures----
png(filename = paste0(outdir,"SuppFigure1f.png"), res = 300, w=7.44,h=2.93, units = 'in')
p
dev.off()

## SourceData----
sourcetable = data[,c("Group","test","pval","lower","upper")]
colnames(sourcetable)[1] = "Model"

write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/submission/data/SupplementaryFigure1f_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
