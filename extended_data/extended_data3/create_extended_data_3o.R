<<<<<<< HEAD:extended_data/extended_data3/create_extended_data_3o.R
### CREATE EXTENDED DATA FIGURE 3L #############################################################################
# Ternary of IDC/ILC (Metabric)

### PREAMBLE #####################################################################################
=======
### CREATE EXTENDED DATA FIGURE 3k #############################################################################
# creates extended data figure 3k
# extended data figure 3k provides the association between recurrence in the ER+ Typical-risk samples and the distance to each archetype, transcriptomic proliferative, HRD LOH scores and histological type in the METABRIC dataset.

### PREAMBLE #####################################################################################
library(yaml)
>>>>>>> c398e855a1f91e3baab879e682684528cdf0bfff:extended_data/extended_data3/create_extended_data_3l.R
library(tidyverse)
library(ggpubr)
<<<<<<< HEAD:extended_data/extended_data3/create_extended_data_3o.R
library(yaml)
library(Ternary)
=======
>>>>>>> c398e855a1f91e3baab879e682684528cdf0bfff:extended_data/extended_data3/create_extended_data_3l.R

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

<<<<<<< HEAD:extended_data/extended_data3/create_extended_data_3o.R
### MAIN ##########################################################################################
tern <- read.delim(file.path(main_repo_path, "data", "ExtendedData3l_sourcetable.txt"), sep = ",")
tern_arc <- tern %>% 
  drop_na(histo_type) %>% 
  filter(histo_type %in% c("ILC", "IDC"))

pdf(file.path(main_repo_path, "plots", 'ExtendedData3l.pdf'),h=3*2, w=3*2)
TernaryPlot(point = "up", alab = "Non amp-associated SVs",
            blab = "Diploidy - Amp-associated SV", clab = "WGD",
            lab.cex = 0.8, grid.minor.lines = 0,
            grid.lines = 5,
            lab.offset = 0.18,
            grid.lty = 'solid', col = rgb(0.9, 0.9, 0.9), grid.col = 'white', 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            padding = 0.12)

TernaryLines(list(c(0, 255, 0), rep(50,3)), col = 'grey')
TernaryLines(list(c(0, 0, 255), rep(50,3)), col = 'grey')
TernaryLines(list(c(255, 0, 0), rep(50,3)), col = 'grey')

AddToTernary(points, tern_arc[which(tern_arc$histo_type == "IDC"),c(2,1,3)], pch = 21, cex = 2, bg = alpha("#a3cfffff", 0.4), col = "white")
AddToTernary(points, tern_arc[which(tern_arc$histo_type == "ILC"),c(2,1,3)], pch = 21, cex = 2, bg = alpha("#335e8dff", 0.7), col = "white")

legend('bottomright', 
       legend = c('IDC', 'ILC'),
       cex = 0.8, bty = 'n', pch = 21, pt.cex = 1.8,
       pt.bg = c(alpha("#a3cfffff", 0.7), alpha("#335e8dff", 0.7)))


data_points <- list(
  "1" = c(0, 0, 255),
  "2" = c(255, 0, 0),
  "3" = c(0, 255, 0))
xy <- CoordinatesToXY(data_points)
points(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), pch = 19, cex = 2.8, col = alpha(c("purple", "blue", "brown"),0.7))
text(xy[1, ] + c(-sqrt(0.1**2/2),0,sqrt(0.1**2/2)), xy[2, ]+c(-sqrt(0.1**2/2),0.1,-sqrt(0.1**2/2)), names(data_points), cex = 0.8, font = 2)
dev.off()

## TME association ----
mb_df <- tern %>%
  filter(Cohort == "METABRIC" &
           group == "ER+ Typical" &
           histo_type %in% c("IDC", "ILC")) %>%
  drop_na(TME_subtype) %>%
  dplyr::rename("Typical"= "Arc2", "High"="Arc1", "TNBC"="Arc3")

model <- lm(Fibrotic ~ Typical  + histo_type, data = mb_df)
summary(model)

model2 <- lm(Immune ~ TNBC  + histo_type, data = mb_df)
summary(model2)

test <- lm(Fibrotic ~ histo_type, data=mb_df)
summary(test)
=======
config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

### MAIN ##########################################################################################
basedir <- "/oak/stanford/groups/ccurtis2/users/"
outdir <- paste0(basedir, "lisem/bc_landscape/submission/figures/")
path.git <- "/oak/stanford/groups/ccurtis2/users/lisem/git_repo/brcarepred/"

# Load data from https://github.com/cclab-brca/brcarepred
clinic = read.table(file=paste0(path.git,"Tables/TableS6.txt"), header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)

tern_arc = read.table(paste0(basedir,"lisem/bc_landscape/github/ExtendedData3d_sourcetable.txt"), header = T, sep="\t")
tern_arc = tern_arc[which(tern_arc$Subgroup == "ER+ Typical"),]

# Add recurrent
tern_arc$DR = sapply(tern_arc$Sample, function(x) if(x %in% clinic$METABRIC.ID){clinic$DR[which(clinic$METABRIC.ID == x)]}else{NA})
tern_arc$LR = sapply(tern_arc$Sample, function(x) if(x %in% clinic$METABRIC.ID){clinic$LR[which(clinic$METABRIC.ID == x)]}else{NA})
tern_arc$recurrence = sapply(1:nrow(tern_arc), function(i) if(tern_arc$Sample[i] %in% clinic$METABRIC.ID){ifelse(tern_arc$DR[i] == 1 | tern_arc$LR[i] == 1, "Recurrent","not-Recurrent")}else{NA} )
table(tern_arc$recurrence)/nrow(tern_arc)

# Add histopathological type
tern_arc$Histological.Type = sapply(tern_arc$Sample, function(x) if(x %in% clinic$METABRIC.ID){clinic$Histological.Type[which(clinic$METABRIC.ID == x)]}else{NA})

x = length(which(tern_arc$Histological.Type == "ILC" & tern_arc$recurrence == "Recurrent"))
y = length(which(tern_arc$Histological.Type == "ILC" & tern_arc$recurrence == "not-Recurrent"))
z = length(which(grepl("IDC",tern_arc$Histological.Type) & tern_arc$Histological.Type != "IDC+ILC" & tern_arc$recurrence == "Recurrent"))
w = length(which(grepl("IDC",tern_arc$Histological.Type) & tern_arc$Histological.Type != "IDC+ILC" & tern_arc$recurrence == "not-Recurrent"))
mat = matrix(c(x,y,z,w), nrow = 2)
fisher.test(mat)$estimate

# Add HRD
mb_input = read.delim(paste0(basedir,"lisem/bc_landscape/ic_subtypes/data/eniclust/metabric_transformed_matrix.txt"))
tern_arc$HRD = sapply(tern_arc$Sample, function(x) if(x %in% mb_input$Sample){mb_input$hrd_loh[which(mb_input$Sample == x)]/4}else{NA}) # HRD_LOH feature was transformed by factor 4 in ENiClust

# Add proliferation rate
mb_tme_sig = read.table(paste0(basedir, "/khoulaha/BreastLandscape/data/2023-11-15_tme_megatable_primary_met.txt"), header = T, sep = "\t")
mb_tme_sig = mb_tme_sig[which(mb_tme_sig$Cohort == "METABRIC"),]
tern_arc$prolif = sapply(tern_arc$Sample, function(x) if(x %in% mb_tme_sig$Sample){mb_tme_sig$Proliferation_rate[which(mb_tme_sig$Sample == x)]}else{NA})

# Merge coefficients together
var_arc <- apply(tern_arc[,2:4],1,var)

summary(lm(var_arc~recurrence,tern_arc))
summary(lm(Arc1~recurrence + var_arc,tern_arc)) #TNBC
summary(lm(Arc2~recurrence + var_arc,tern_arc)) #Typical yes and neg
summary(lm(Arc3~recurrence + var_arc,tern_arc)) #High/HER2+ no but pos

coefficient <- c(summary(lm(Arc1~recurrence,tern_arc))$coefficient["recurrenceRecurrent","Estimate"], 
                 summary(lm(Arc2~recurrence,tern_arc))$coefficient["recurrenceRecurrent","Estimate"], 
                 summary(lm(Arc3~recurrence,tern_arc))$coefficient["recurrenceRecurrent","Estimate"],
                 fisher.test(mat)$estimate,
                 summary(lm(HRD~recurrence,tern_arc))$coefficient["recurrenceRecurrent","Estimate"],
                 summary(lm(prolif~recurrence,tern_arc))$coefficient["recurrenceRecurrent","Estimate"])
p_values <- c(lmp(lm(Arc1~recurrence,tern_arc)), 
              lmp(lm(Arc2~recurrence,tern_arc)), 
              lmp(lm(Arc3~recurrence,tern_arc)),
              fisher.test(mat)$p.value,
              lmp(lm(HRD~recurrence,tern_arc)),
              lmp(lm(prolif~recurrence,tern_arc)))

significance <- sapply(p_values, function(x) if(x > 0.05){""}else{ifelse(between(x,0.01,0.05),"*",ifelse(between(x,0.001,0.01),"**",ifelse(x<0.001,"***",NA)))})

data <- data.frame(variable = c("TNBC -enriched","ER+ Typical -enriched","ER+ High/HER2+ -enriched",
                                "Histology ILC","HRD LOH score","Proliferation"),
                   coefficient = coefficient, 
                   pvalue = p_values,
                   significance = significance,
                   group = c(rep("Archetype",3),rep("Other features",3)), stringsAsFactors = F)

### SAVE ##########################################################################################
## Figure----
data$variable <- factor(data$variable, levels = c("ER+ High/HER2+ -enriched", "ER+ Typical -enriched", "TNBC -enriched",
                                                  "Histology ILC","HRD LOH score","Proliferation"))

labels = paste(data$variable, paste0("(",data$significance,")"))
names(labels) = c("TNBC -enriched", "ER+ Typical -enriched", "ER+ High/HER2+ -enriched", 
                  "Histology ILC","HRD LOH score","Proliferation")

png(filename = file.path(outdir, 'ExtendedData3k.png'), res = 300, width = 6.11, height = 4.80, units = 'in')
ggplot(data, aes(x=variable, y=coefficient)) +
  geom_segment(aes(x=variable, xend=variable, y=0, yend=coefficient, col=variable), linewidth = 1.5) +
  geom_point(aes(col=variable, fill=variable), size=5, alpha=0.7, shape=21, stroke=2) +
  scale_color_manual(name = "", values = c("TNBC -enriched"="#bb62f5ff", "ER+ Typical -enriched"="#4d4dffff", "ER+ High/HER2+ -enriched"="#bf6969ff","Histology ILC"="grey","HRD LOH score"="grey","Proliferation"="grey")) +
  scale_fill_manual(name = "", values = c("TNBC -enriched"="#bb62f5ff", "ER+ Typical -enriched"="#4d4dffff", "ER+ High/HER2+ -enriched"="#bf6969ff","Histology ILC"="grey","HRD LOH score"="grey","Proliferation"="grey")) +
  facet_wrap(group~., scales = "free", nrow = 2) +
  xlab('') + ylab('Coefficients') +
  scale_x_discrete(labels = labels) +
  coord_flip() + theme_LM + theme(legend.position = "none",
                                  axis.text.x = element_text(size=14),
                                  axis.text.y = element_text(size=14, angle=25),
                                  axis.title = element_text(size=20),
                                  strip.text = element_text(size=20))
dev.off()

## SourceData---- 
write.table(data, paste0(basedir,"lisem/bc_landscape/github/ExtendedData3k_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
>>>>>>> c398e855a1f91e3baab879e682684528cdf0bfff:extended_data/extended_data3/create_extended_data_3l.R
