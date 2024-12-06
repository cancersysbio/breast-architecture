# Extract alterations affecting immune escape pathways 
# Data used in Figure 4C, 4D, Extended Data 9A, 9E

library(tidyverse)

### PREAMBLE #####################################################################################
library(tidyverse)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################

# Read data 
primary_megatable <- read.delim(file.path(main_repo_path, 'data','primary_megatable.txt'))
metastatic_megatable <- read.delim(file.path(main_repo_path, 'data','metastatic_megatable.txt'))

samples <- rbind(
  primary_megatable %>% select(Sample, group, genome_doubled) %>% mutate(Stage = "Primary"), 
  metastatic_megatable %>% select(Sample, group, genome_doubled) %>% mutate(Stage = "Metastatic")
)

### Set data paths #################################################################
# Copy number data is merged gene level data as output by FACETS
cn_prim_path <- "/path/to/primary/CN/data.tsv"
cn_met_path <- "/path/to/metastatic/CN/data.tsv"

# List of individuals homozygous for HLA alleles 
# Columns Sample, A_hom, B_hom, C_hom
homozygous_individuals <- "/path/to/hla/homozygotes.txt"

# Damaging SVs 
sv_prim_path <- "/path/to/primary/damaging/SV/data.txt"
sv_met_path <- "/path/to/metastatic/damaging/SV/data.txt"

# Complex alterations 
complex_amp_prim_path <- "/path/to/primary/complex/alt/data.txt"
complex_amp_met_path <- "/path/to/metastatic/complex/alt/data.txt"




## Genes and signatures ----
gene_signature_path <- file.path(main_repo_path, 'data','gene_signatures.gmt')
if(!file.exists(gene_signature_path)) {
  stop("Error: TME gene signature path not found, download here: https://github.com/BostonGene/MFP/blob/master/signatures/gene_signatures.gmt")
} 

tme_fges <- GSEABase::getGmt(gene_signature_path)

immune_escape_genes_path <- file.path(main_repo_path, 'data','immune_selected_genes.tsv')
if(!file.exists(immune_escape_genes_path)) {
  stop("Error: Immune escape gene list path not found, download here: https://github.com/UMCUGenetics/Genetic-Immune-Escape/blob/main/external_data/immune_selected_genes.tsv")
} 

hartwig_pathways <- read.delim(immune_escape_genes_path) %>% select(Hugo_Symbol = Gene, Pathway.general)

# Add HLA-II genes minus CIITA (already in AP pathway)
hartwig_pathways <- rbind(
  hartwig_pathways, 
  data.frame(Hugo_Symbol = unname(unlist(GSEABase::geneIds(tme_fges["MHCII"])))[1:8], 
             Pathway.general = "HLA-II")
) %>% rename(Pathway = Pathway.general)


## CN data ------
# Primary
cols <- rep("NULL", 25)
cols[c(1, 2, 16:17, 24:25)] <- "character"
cn_prim <- read.delim(cn_prim_path, colClasses = cols)

# Metastatic 
cn_met <- read.delim(cn_met_path, colClasses = cols)

cn <- rbind(cn_prim, cn_met)

# Filtering for immune genes 
cn <- cn %>% filter(filter == "PASS", 
                    cn_state != "DIPLOID", 
                    gene %in% hartwig_pathways$Hugo_Symbol)

length(unique(cn$gene))

amp_del <- cn %>% mutate(Alteration = ifelse(grepl("AMP", cn_state), "Amplification", 
                                             ifelse(cn_state == "HOMDEL", "Homozygous deletion", NA))) %>% 
  filter(!is.na(Alteration))

hla_loh <- cn %>% filter(grepl("HLA", gene)) %>% mutate(tcn = as.numeric(tcn), lcn = as.numeric(lcn)) %>% 
  filter((tcn == 1)|((lcn == 0)&(tcn > 0))) %>% mutate(Alteration = "HLA LOH")

# Remove if individual was homozygous for that allele
hla_hom <- read.delim(homozygous_individuals, sep = ",")
hla_hom <- hla_hom %>% select(Individual = Sample, `HLA-A` = A_hom, `HLA-B` = B_hom, `HLA-C` = C_hom) %>% pivot_longer(cols = 2:4, names_to = "gene")

hla_loh_hom <- hla_loh %>% mutate(Individual = substr(sample, 1, 14)) %>% inner_join(hla_hom) %>% filter(value) 

hla_loh <- hla_loh %>% anti_join(hla_loh_hom)

cn <- rbind(amp_del, hla_loh) %>% select(Sample = sample, Hugo_Symbol = gene, Alteration)


## Damaging SVs ----------
sv_prim <- read.delim(sv_prim_path)
sv_met <- read.delim(sv_met_path)

sv <- rbind(sv_prim, sv_met)
sv <- sv %>% filter(Gene %in% hartwig_pathways$Hugo_Symbol)

length(unique(sv$Gene))

sv <- sv %>% select(Sample = sample, Hugo_Symbol = Gene) %>% mutate(Alteration = "Structural variant")

## Complex alterations ------------
# Considering only amplifications in PD-L1 and SETDB1 as all the others will get filtered
complex_amp_prim <- read.delim(complex_amp_prim_path)
complex_amp_met <- read.delim(complex_amp_met_path)

complex_amp <- rbind(complex_amp_prim, complex_amp_met)

complex_amp <- rbind(
  complex_amp %>% filter(grepl("CD274", immune_tme)) %>% distinct(sample, ampliconID, .keep_all = T) %>% mutate(Hugo_Symbol = "CD274"), 
  complex_amp %>% filter(grepl("SETDB1", oncogenes)) %>% distinct(sample, ampliconID, .keep_all = T) %>% mutate(Hugo_Symbol = "SETDB1")
) %>% select(Sample = sample, amplicon_type, Hugo_Symbol)

complex_amp <- complex_amp %>% filter(!amplicon_type %in% c("Linear amplification", "No amp/Invalid"))
complex_amp <- complex_amp %>% mutate(Alteration = paste0(amplicon_type, " amplification")) %>% 
  select(-amplicon_type)

## Removing duplicates -------
# For samples with HLA LOH or complex amplifications, remove calls from facets for same gene/sample
duplicated_alterations <- complex_amp %>% select(-Alteration) %>% inner_join(cn) 
cn <- cn %>% anti_join(duplicated_alterations)

## SNV data -- already together and filtered 
snv <- read.delim(config$snv$protein_coding_immune)
snv <- snv %>% select(Sample = Tumor_Sample_Barcode, Hugo_Symbol) %>% mutate(Alteration = "Single nucleotide variant")

## Combine all alterations ---------
alterations <- rbind(cn, sv, snv, complex_amp)

# Filter alterations to match Hartwig paper 
alterations <- alterations %>% filter(case_when((Hugo_Symbol %in% c("SETDB1", "CD274")) & (grepl("mplification", Alteration)) ~ T, 
                                                !(Hugo_Symbol %in% c("SETDB1", "CD274")) & (!grepl("mplification", Alteration)) ~ T, .default = F)) %>% 
  left_join(hartwig_pathways) 

# Filter and add stage
alterations <- rbind(
  alterations %>% filter(Sample %in% primary_megatable$Sample) %>% mutate(Stage = "Primary"), 
  alterations %>% filter(Sample %in% metastatic_megatable$Sample) %>% mutate(Stage = "Metastatic")
)

write.table(alterations, file.path(main_repo_path, "data", "immune_pathway_alterations.txt"), 
            sep = "\t", quote = F, row.names = F)



