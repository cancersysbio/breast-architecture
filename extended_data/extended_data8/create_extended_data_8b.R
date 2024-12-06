### CREATE EXTENDED DATA FIGURE 8b #############################################################################
# creates extended data figure 8b
# extended data figure 8b provides the 

### PREAMBLE #####################################################################################
library(yaml)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
require(scales)

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
path_rna_tme <- "/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data"
outdir <- paste0(basedir,"lisem/bc_landscape/submission/figures/")

### Load IMC data
scRNA = read.table(paste0(basedir,'lisem/bc_spatial/imc/Ali_Nat_Cancer_2020/20191219-ftp/METABRIC_IMC/to_public_repository/single_cell_data.csv'), header = T, sep = ',')

### Creating dataframe
data_sample = data.frame(sample = rep(unique(scRNA$metabricId), each = length(unique(scRNA$description))),
                         celltype = rep(unique(scRNA$description), length(unique(scRNA$metabricId))),
                         stringsAsFactors = F)

data_sample = data_sample %>% dplyr::rowwise() %>% mutate(proportion = length(which(scRNA$metabricId == sample & scRNA$description == celltype)) / length(which(scRNA$metabricId == sample)))

## Cross cell type proportions with TME signatures
mb_tme = read.table(paste0(path_rna_tme, "/2023-11-15_tme_megatable_primary_met.txt"), header = T, sep = "\t")
mb_tme = mb_tme[which(mb_tme$Cohort == "METABRIC"),]

data_sample = data_sample %>% 
  dplyr::rowwise() %>% dplyr::mutate(tme = ifelse(sample %in% mb_tme$Sample, mb_tme$TME_subtype[which(mb_tme$Sample == sample)] , NA))

data_sample$tme = factor(data_sample$tme, levels = rev(c("D","F","IE/F","IE")))
dodge <- position_dodge(width = 0.5)

sort(unique(data_sample$celltype))
data_sample_sub = data_sample[which(!is.na(data_sample$tme) & data_sample$celltype %in% c("B cells","Fibroblasts","Fibroblasts CD68+","Myofibroblasts","T cells","Vascular SMA+")),]
data_sample_sub$celltype = factor(data_sample_sub$celltype, levels = c("B cells","Fibroblasts","Fibroblasts CD68+","T cells","Myofibroblasts","Vascular SMA+"))
data_sample_sub$tme = factor(data_sample_sub$tme, levels = c("D","F","IE/F","IE"))

my_comparisons <- list( c("D", "F"), c("F", "IE"), c("F", "IE/F"), c("D","IE"), c("D","IE/F"), c("IE/F","IE"))

p <- ggplot(data = data_sample_sub, aes(x = tme, y = proportion, col = tme)) +
  geom_violin(aes(fill = tme), alpha = 0.5, position = dodge) +
  geom_boxplot(width = 0.2, position = dodge) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  scale_color_manual(values = c("F"="#b80045ff","D"="#bb8800ff","IE"="#572eccff","IE/F"="#2ec0baff")) + 
  scale_fill_manual(values = c("F"="#b80045ff","D"="#bb8800ff","IE"="#572eccff","IE/F"="#2ec0baff")) + 
  xlab("") + ylab("Mean proportion") +
  facet_wrap(~ celltype, ncol = 3, nrow = 2) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test") + # Add pairwise comparisons p-value
  theme_LM + ggtitle("Cell type proportion in each TME subtype") + theme(legend.position = "None", strip.background = element_blank(), 
                                                                         strip.text.x = element_text(size = 16, face = "bold"),
                                                                         axis.text.x = element_text(size=16),
                                                                         axis.text.y = element_text(size=16),
                                                                         axis.title.y = element_text(size=18), 
                                                                         title = element_text(size=20)) 

### SAVE ##########################################################################################
## Figure----
png(filename = file.path(outdir, 'ExtendedData8b.png'), res = 300, width = 7.82, height = 7.28, units = 'in')
p
dev.off()

## SourceData---- 
colnames(data_sample_sub) = c("Sample","CellType","Proportion","TME")
sourcetable = data_sample_sub
length(unique(sourcetable$Sample))
write.table(sourcetable, "/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/github/ExtendedData8b_sourcetable.txt", row.names = F, col.names = T, quote = F, sep="\t")
