### CREATE SUPPLEMENTARY FIGURE 4e #############################################################################
# creates supplementary figure 4e
# supplementary figure 4e provides the correlation between various mutational features and the distance to each archetype with archetypes defined by Pareto or NMF.

### PREAMBLE #####################################################################################
library(yaml)
library(BoutrosLab.plotting.general, lib = "/home/lisem/R/x86_64-pc-linux-gnu-library/4.1/")
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
date <- Sys.Date()
### ADD JABBA AND FACETS EVENTS ###
add_jabba_and_facets_events <- function(arch, jabbafile, facetsfile) {
  # read in jabba 
  jabba <- read.delim(
    jabbafile,
    as.is = TRUE
  )
  # merge 
  arch <- merge(arch[,1:4], jabba[,1:15], by = 'sample', all.x = TRUE)
  # read in facets loh 
  facets <- read.delim(
    facetsfile,
    as.is = TRUE
  )
  arch <- merge(arch, facets[,c('sample','hrd_loh')], by = 'sample', all.x = TRUE)
  return(arch)
}

### ADD SBS SIGNATURES ###
add_sbs_signatures <- function(arch, sbsfile, megatable, primary = TRUE) {
  # read in sbs 
  sbs <- readRDS(sbsfile)
  sbs <- sbs[which(sbs$Sample %in% megatable$sample),]
  colnames(sbs) <- gsub("Sample","sample", colnames(sbs))
  # only keep signatures that have at least 10% contribution in 10% of patients
  if (primary) {
    nonzero <- colSums(sbs[,-1] > 10)
    nonzero_sig <- names(nonzero)[nonzero > 0.1*nrow(sbs)]
    sbs <- sbs[,c('sample',nonzero_sig)]
  }
  arch <- merge(
    arch,
    sbs,
    by = 'sample'
  )
  return(arch)
}

### TEST FEATURES AGAINST NMF LFs ###
test_features_against_archtypes <- function(arch, features) {
  # test correlation with NMF LFs
  res <- do.call(rbind, sapply(
    features,
    function(x) {
      print(x)
      arch1 <- cor.test(arch$LF1, arch[,x], method = 'spearman')
      arch2 <- cor.test(arch$LF2, arch[,x], method = 'spearman')
      arch3 <- cor.test(arch$LF3, arch[,x], method = 'spearman')
      data.frame(
        rho = c(arch1$estimate, arch2$estimate, arch3$estimate),
        p = c(arch1$p.value, arch2$p.value, arch3$p.value),
        archetype = c('LF1','LF2','LF3'),
        feature = x
      )
    },
    simplify = FALSE
  ))
  res$fdr <- p.adjust(res$p, method = 'fdr')
  return(res)
}

### SET FEATURE ORDER ###
set_feature_order <- function(plot_data) {
  tmp1 <- plot_data[plot_data$archetype == 'LF1',]
  feat_order1 <- tmp1[order(tmp1$rho),'feature']
  tmp3 <- plot_data[plot_data$archetype == 'LF3',]
  feat_order3 <- tmp3[order(tmp3$rho),'feature']
  tmp2 <- plot_data[plot_data$archetype == 'LF2',]
  feat_order2 <- tmp2[order(tmp2$rho),'feature']
  # minor edits to feature order
  feat_order <- c(
    grep(
      'SBS',
      unique(c(as.character(feat_order1), as.character(feat_order3),as.character(feat_order2))),
      value = TRUE
    ),
    grep(
      'SBS',
      unique(c(as.character(feat_order1), as.character(feat_order3),as.character(feat_order2))),
      invert = TRUE,
      value = TRUE
    )
  )
  # swapping dm and inv to match categorization 
  feat_order <- feat_order[c(1:10,17,12,11,13:16,18:27)]
  return(feat_order)
}

### MAIN ##########################################################################################
basedir <- "/oak/stanford/groups/ccurtis2/users"
kh_basedir <- file.path(basedir, "khoulaha/BreastLandscape")
lm_basedir <- file.path(basedir, "lisem/bc_landscape")
sp_basedir <- file.path(basedir, "syparkmd/Projects/BCgenomic")

# read NMF LFs 
arch <- read.delim(
  file.path(lm_basedir, 'revision/data/nmf_rank_prim.txt'),
  as.is = TRUE, sep= " "
)
arch$Sample = rownames(arch)
arch = arch[,c(4,1:3)]
colnames(arch) <- gsub('Sample','sample', colnames(arch))

# read in jabba 
arch <- add_jabba_and_facets_events(
  arch = arch,
  jabbafile = file.path(kh_basedir, 'svs/jabba/primary/2023-05-08_jabba_events.txt'),
  facetsfile = file.path(basedir, '..', 'isabl/analyses/65/27/16527/bestfit_merged.tsv')
)

# read in megatable 
megatable <- read.delim(
  file.path(kh_basedir, 'data/2024-09-05_primary_megatable.txt'),
  as.is = TRUE
)
megatable_subset <- megatable[,c('Sample','genome_doubled','fraction_cna','SVB','AA_ecdna','AA_complex','group')]
colnames(megatable_subset) <- gsub('Sample','sample', colnames(megatable_subset))
# merge with archtypes
arch <- merge(arch, megatable_subset, by = 'sample', all.x = TRUE)
# reformat genome doubling
arch$genome_doubled <- (arch$genome_doubled == 'True')*1
# add dm and cmplxdm together
arch$dm <- arch$dm+arch$cpxdm
# read in sbs 
arch <- add_sbs_signatures(
  arch = arch,
  sbsfile = file.path(sp_basedir, 'Meta_data/sig_cont_intermediate.rds'),
  megatable = megatable_subset
)

# test correlation with NMF LFs
features <- colnames(arch)[c(5:24,26:34)]
features <- features[features != 'cpxdm']
res <- test_features_against_archtypes(arch = arch, features = features)

### METASTATIC ####################################################################################
# read NMF LFs 
metarch <- read.delim(
  file.path(lm_basedir, 'revision/data/nmf_rank_met.txt'),
  as.is = TRUE, sep= " "
)
metarch$Sample = rownames(metarch)
metarch = metarch[,c(4,1:3)]
colnames(metarch) <- gsub('Sample','sample',colnames(metarch))

# read in jabba and facet events
metarch <- add_jabba_and_facets_events(
  arch = metarch,
  jabbafile = file.path(kh_basedir, 'svs/jabba/metastatic/2023-06-09_jabba_events.txt'),
  facetsfile = file.path(basedir, '..', 'isabl/analyses/32/93/23293/bestfit_merged.tsv')
)

# read in megatable 
metmegatable <- read.delim(
  file.path(kh_basedir, 'data/2024-09-05_metastatic_megatable.txt'),
  as.is = TRUE
)
metmegatable_subset <- metmegatable[,c('Sample','genome_doubled','fraction_cna','SVB','AA_ecdna','AA_complex','group')]
colnames(metmegatable_subset) <- gsub('Sample','sample', colnames(metmegatable_subset))
metarch <- merge(metarch, metmegatable_subset, by = 'sample')
# read in sbs 
metarch <- add_sbs_signatures(
  arch = metarch,
  sbsfile = file.path(sp_basedir, 'Meta_data/sig_cont_intermediate.rds'),
  megatable = metmegatable_subset,
  primary = FALSE
)

metarch$genome_doubled <- (metarch$genome_doubled == 'True')*1
metarch$dm <- metarch$dm+metarch$cpxdm
# test features
metres <- test_features_against_archtypes(arch = metarch, features = features)

### SAVE ##########################################################################################
## Figure----
# find sig features and create plot data
sig_features <- unique(res[res$fdr < 0.05,'feature'])
primplot_data <- res[res$feature %in% sig_features,]
metplot_data <- metres[metres$feature %in% sig_features,]
metplot_data$stage <- 'metastatic'
primplot_data$stage <- 'primary'
plot_data <- rbind(primplot_data, metplot_data)

plot_data$code <- NA
plot_data[plot_data$archetype == 'LF2' & plot_data$stage == 'metastatic','code'] <- '2. ER+ Typical -enriched'
plot_data[plot_data$archetype == 'LF3' & plot_data$stage == 'metastatic','code'] <- '3. HER2+/ER+ High -enriched'
plot_data[plot_data$archetype == 'LF1' & plot_data$stage == 'metastatic','code'] <- '1. TNBC -enriched'
plot_data[plot_data$archetype == 'LF1' & plot_data$stage == 'primary','code'] <- '1. TNBC -enriched'
plot_data[plot_data$archetype == 'LF2' & plot_data$stage == 'primary','code'] <- '2. ER+ Typical -enriched'
plot_data[plot_data$archetype == 'LF3' & plot_data$stage == 'primary','code'] <- '3. HER2+/ER+ High -enriched'

# keep features in main plot
main <- c("SBS1","SBS18","SBS8","SBS39","AA_ecdna","tyfonas","AA_complex","bfb","SVB","tic","hrd_loh","del","dup","genome_doubled","fraction_cna")
feat_xlab <- c("SBS1","SBS18","SBS8","SBS39","ecdna","tyfonas","complex amp","bfb","SVB","tic","HRD LOH","del","dup","wgd","fga")

main_plot_data <- plot_data[plot_data$feature %in% main,]

main_plot_data$index <- NA
for (i in 1:length(main)) {
  feattmp <- main[i]
  main_plot_data[main_plot_data$feature == feattmp,'index'] <- i
}

cov.legend <- list(
  legend = list(
    colours = c('#9ebff1ff','#fe85cfff'),
    labels = c('Primary', 'Metastatic'),
    border = 'transparent'
  )
)
cov.legend.grob <- legend.grob(legends = cov.legend)

create.scatterplot( # w=15.86, h=5.07
  rho ~ index | code,
  data = main_plot_data,
  #horizontal = TRUE,
  groups = main_plot_data$stage,
  type = c('h','p'),
  abline.h = 0,
  abline.v = c(4.5,8.5),
  abline.lty = 2,
  xaxis.rot = 45,
  xaxis.cex = 1.25,
  ylimits = c(-1,1),
  xlimits = c(0,16),
  #xaxis.rot = 45,
  yat = seq(-1,1,0.5),
  yaxis.cex = 1,
  col = c('#fe85cfff','#9ebff1ff'),
  ylab.label = 'Spearman\'s Rho',
  xlab.label = 'Genomic Features',
  add.text = TRUE,
  text.labels = c('Simple SVs','Complex\nAmplifications','SNV\nSignatures'),
  text.x = c(12.5,6.5,2.5),
  text.y = c(0.91,0.85,0.85),
  layout = c(3,1),
  legend = list(inside = list(fun = cov.legend.grob, x = 0, y = 0.12)),
  text.cex = 1,
  xaxis.lab = feat_xlab,
  xat = 1:length(main),
  width = 17,
  height = 5,
  #filename = paste0('/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/publication/figures/', 'nmf_svs_features_scatterplot.pdf'),
  resolution = 300
)

## SourceData---- 
sourcedata <- plot_data
colnames(sourcedata)[3] <- "nmf_lf"
write.table(sourcedata, paste0(basedir,"/lisem/bc_landscape/github/SupplementaryFigure4e_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
