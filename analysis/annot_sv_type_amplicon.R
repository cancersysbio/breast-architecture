# Annotate SV types overlapping amplicon locations ----
# ED4E (primary) and SF5c (metastatic)

# PREAMBLE ----
library(rstudioapi)
library(tidyr)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(IRanges)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# MAIN ----
her2_samples_w_amplicons <- "list of HER2+ samples with amplicons, e.g. from coamp_long"

# Read data
sv <- fread("/path/to/primary/all/SV/data.txt") %>%
  mutate(Sample=substr(sample, 1, 14)) %>%
  select(-sample) %>%
  filter(Sample %in% her2_samples_w_amplicons) # subset samples with ic amplicon

hits_strand1 = list()
for (i in seq_along(her2_samples_w_amplicons)){
  coamp_temp = coamp_long[coamp_long$Sample == her2_samples_w_amplicons[i], ]
  
  query <- GRanges(
    seqnames = coamp_temp$gencode_seqnames_hg19,
    ranges = IRanges(start = coamp_temp$gencode_start_hg19, end = coamp_temp$gencode_end_hg19),
    strand = '+',
    metadata = coamp_temp$Sample
  )
  mcols(query)$Cytoband_hg19 <- coamp_temp$Cytoband_hg19
  
  sv_temp <- sv[sv$Sample == her2_samples_w_amplicons[i], ]
  
  subject_strand1 <- GRanges(
    seqnames = paste0('chr', sv_temp$chrom1),
    ranges = IRanges(start = sv_temp$pos1, end = sv_temp$pos1_end),
    strand = sv_temp$strand1,
    metadata = sv_temp$Sample
  )
  mcols(subject_strand1)$svclass <- sv_temp$svclass
  mcols(subject_strand1)$name <- sv_temp$name
  
  
  hits_temp <- IRanges::findOverlaps(query, subject_strand1,
                                     ignore.strand=F,
                                     type="any",
                                     select="all")
  
  hits_strand1[[i]] <- data_frame(
    sample = query[hits_temp@from]$metadata,
    cytoband = query[hits_temp@from]$Cytoband_hg19,
    sv_type = subject_strand1[hits_temp@to]$svclass,
    sv_name = subject_strand1[hits_temp@to]$name
  )
}
hits_strand1_df <- do.call(rbind, hits_strand1) %>%
  distinct(sample, sv_name, .keep_all = TRUE) %>%
  mutate(strand='strand1')


hits_strand2 = list()
for (i in seq_along(her2_samples_w_amplicons)){
  coamp_temp = coamp_long[coamp_long$Sample == her2_samples_w_amplicons[i], ]
  
  query <- GRanges(
    seqnames = coamp_temp$gencode_seqnames_hg19,
    ranges = IRanges(start = coamp_temp$gencode_start_hg19, end = coamp_temp$gencode_end_hg19),
    strand = '+',
    metadata = coamp_temp$Sample
  )
  mcols(query)$Cytoband_hg19 <- coamp_temp$Cytoband_hg19
  
  sv_temp <- sv[sv$Sample == her2_samples_w_amplicons[i], ]
  subject_strand2 <- GRanges(
    seqnames = paste0('chr', sv_temp$chrom2),
    ranges = IRanges(start = sv_temp$pos2, end = sv_temp$pos2_end),
    strand = sv_temp$strand2,
    metadata = sv_temp$Sample
  )
  mcols(subject_strand2)$svclass <- sv_temp$svclass
  mcols(subject_strand2)$name <- sv_temp$name
  
  
  hits_temp <- IRanges::findOverlaps(query, subject_strand2,
                                     ignore.strand=F,
                                     type="any",
                                     select="all")
  
  hits_strand2[[i]] <- data_frame(
    sample = query[hits_temp@from]$metadata,
    cytoband = query[hits_temp@from]$Cytoband_hg19,
    sv_type = subject_strand2[hits_temp@to]$svclass,
    sv_name = subject_strand2[hits_temp@to]$name
  )
  
}

hits_strand2_df <- do.call(rbind, hits_strand2) %>%
  distinct(sample, sv_name, .keep_all = TRUE) %>%
  mutate(strand='strand2')

hits_strand_df <- rbind(hits_strand1_df, hits_strand2_df)
samples_sv_type_in_amplicon <- hits_strand_df %>%
  dplyr::distinct(sample, sv_type, cytoband, .keep_all = TRUE)

samples_sv_type_in_amplicon %>%
  fwrite(file.path(main_repo_path, "data", "ExtendedData4e_sourcetable2.txt"))
