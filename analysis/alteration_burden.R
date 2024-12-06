# Calculate alteration burden

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
# CN is given as merged FACETS CNCF output files 
cn_prim_path <- "/path/to/primary/CN/data.tsv"
cn_met_path <- "/path/to/metastatic/CN/data.tsv"

# Damaging SVs 
sv_prim_path <- "/path/to/primary/damaging/SV/data.txt"
sv_met_path <- "/path/to/metastatic/damaging/SV/data.txt"

# All SVs 
sv_prim_all_path <- "/path/to/primary/all/SV/data.txt"
sv_met_all_path <- "/path/to/metastatic/all/SV/data.txt"

# Protein coding SNVs -- merged primary and metastatic
coding_snv_path <- "/path/to/coding/snv.txt"

# Alterations -----
## CN --------
cn_prim <- read.delim(cn_prim_path) %>% 
  select(Sample = ID, chrom, loc.start, loc.end, tcn, lcn) %>% mutate(Sample = substr(Sample, 1, 26))

cn_met <- read.delim(cn_met_path) %>% 
  select(Sample = ID, chrom, loc.start, loc.end, tcn, lcn) %>% mutate(Sample = substr(Sample, 1, 26))

cn <- rbind(cn_prim, cn_met)

# CN segments need to be annotated, FACETS calls AMP when > 4 without WGD, or when > 5 with WGD
cn_filt <- cn %>% inner_join(samples %>% select(Sample, genome_doubled)) %>% mutate(Alteration = ifelse(
  tcn == 0, "HOMDEL", ifelse((genome_doubled == "True")&(tcn > 5),"AMP", ifelse((genome_doubled == "False")&(tcn > 4), "AMP", NA))  
)) %>% filter(!is.na(Alteration))

## Structural variants --------
# Damaging
sv_prim <- read.delim(sv_prim_path)
sv_met <- read.delim(sv_met_path)

sv <- rbind(sv_prim, sv_met)
sv_filt <- sv %>% select(Sample = sample, SV_ID) %>% distinct()

damaging_sv_ids <- sv_filt %>% mutate(SV_ID = stringr::str_split(SV_ID, "-", simplify = T)[,1]) %>% mutate(SV_ID = substr(SV_ID, 28, 100)) %>% pull(SV_ID) %>% unique()

# All 
sv_all_prim <- read.delim(sv_prim_all_path)
sv_all_met <- read.delim(sv_met_all_path)

sv_all <- rbind(sv_all_prim, sv_all_met)

sv_non_damaging <- sv_all %>% filter(!(name %in% damaging_sv_ids)) %>% select(Sample = sample, name) %>% distinct()

## SNVs ------
snvs <- read.delim(coding_snv_path)
snvs_interest <- c("Missense_Mutation", "Splice_Site", "Frame_Shift_Del", "Splice_Region", 
                   "Nonsense_Mutation", "Frame_Shift_Ins", "In_Frame_Ins", "In_Frame_Del", 
                   "Translation_Start_Site", "Nonstop_Mutation")

snvs_filt <- snvs %>% filter(Variant_Classification %in% snvs_interest)

# Combine alterations and save ----
alterations <- rbind(
  cn_filt %>% select(Sample, Alteration) %>% mutate(Alteration = recode(Alteration, "AMP" = "Amplification", "HOMDEL" = "Homozygous deletion")), 
  sv_filt %>% select(Sample) %>% mutate(Alteration = "Damaging structural variant"), 
  sv_non_damaging %>% select(Sample) %>% mutate(Alteration = "Non damaging structural variant"), 
  snvs_filt %>% select(Sample = Tumor_Sample_Barcode) %>% mutate(Alteration = "Protein coding SNV")
) %>% inner_join(samples %>% select(Sample, group, Stage)) 

write.table(alterations, file.path(main_repo_path, "data", "alteration_burden.txt"), 
            quote = F, sep = "\t", row.names = F)
