### CREATE EXTENDED DATA 2C + 2F #############################################################################
# 2C -- Alteration burden in metastatic ER+ samples by treatment received
# 2F -- Alteration burden in metastatic samples by metastatic site
### PREAMBLE #####################################################################################
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

# Read data -- generated in alteration_burden.R
alterations_data_path <- file.path(main_repo_path, 'data','alteration_burden.txt')
if(!file.exists(alterations_data_path)) {
  stop("Error: Alteration burden data not found, run alteration_burden.R analysis script. ")
} 

alterations <- read.delim(alterations_data_path)

primary_megatable <- read.delim(file.path(main_repo_path, 'data','primary_megatable.txt'))
metastatic_megatable <- read.delim(file.path(main_repo_path, 'data','metastatic_megatable.txt'))

samples <- rbind(
  primary_megatable %>% select(Sample, group, genome_doubled) %>% mutate(Stage = "Primary"), 
  metastatic_megatable %>% select(Sample, group, genome_doubled) %>% mutate(Stage = "Metastatic")
)


# Separate by met sites, Extended Data 2F ----
hartwig_clinical <- read.delim(file.path(main_repo_path, "data", "Hartwig_clinical.txt")) %>% 
  mutate(sample = substr(Individual.System.ID, 1, 21)) %>% 
  inner_join(metastatic_megatable %>% mutate(sample = substr(Sample, 1, 21))) %>% select(-Sample) %>%
  mutate(biopsySite_LM = gsub("Subcateneous", "Subcutaneous", biopsySite_LM))

sites_interest <- table(hartwig_clinical$biopsySite_LM)[table(hartwig_clinical$biopsySite_LM) > 10] %>% names()
alterations_met_site <- alterations %>% inner_join(hartwig_clinical %>% rename(Sample = 3)) %>% filter(biopsySite_LM %in% sites_interest, biopsySite_LM != "Unknown")

samples_per_group_site <- hartwig_clinical %>% filter(!is.na(group), biopsySite_LM %in% sites_interest) %>% group_by(group, biopsySite_LM) %>% summarize(samples_group = n()) 
samples_per_site <- hartwig_clinical %>% filter(!is.na(group), biopsySite_LM %in% sites_interest) %>% group_by(biopsySite_LM) %>% summarize(samples_site = n()) 

source_data_2f <- alterations_met_site %>% group_by(biopsySite_LM, Alteration) %>% summarize(n = n()) %>% 
  left_join(samples_per_site) %>% mutate(alterations_per_sample = n/samples_site) %>% 
  group_by(biopsySite_LM) %>% mutate(total_alts = sum(alterations_per_sample)) %>%
  arrange(total_alts)

## Save plot ----
pdf(file.path(main_repo_path, "plots", "ExtendedData2F.pdf"), height = 6, width = 6)
ggbarplot(source_data_2f, "biopsySite_LM", "alterations_per_sample", 
          color = "Alteration", fill = "Alteration", palette = alteration_colors,
          ylab = "Alterations per sample", xlab = "Biopsy site") +
  geom_text(aes(label = samples_site, y = total_alts*1.07), check_overlap = T) +
  theme_LM +
  rotate_x_text(35) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(ncol = 1), 
         color = guide_legend(ncol = 1))
dev.off()



# ER+ High risk vs ER+ Typical risk, Extended Data 2C -----
er <- c("ER+ High", "ER+ Typical")

alts_treatment <- alterations_met_site %>% filter(group %in% er) %>% 
  select(Sample, group, hasSystemicPreTreatment, hasRadiotherapyPreTreatment, Alteration) %>%
  filter(hasSystemicPreTreatment != "NULL", hasRadiotherapyPreTreatment != "NULL") %>%
  mutate(pretreatmentType = ifelse(hasSystemicPreTreatment == "Yes", 
                                   ifelse(hasRadiotherapyPreTreatment == "Yes", "Both", "Systemic"), 
                                   ifelse(hasRadiotherapyPreTreatment == "Yes", "Radiotherapy", "Neither"))) %>% 
  select(-hasSystemicPreTreatment, -hasRadiotherapyPreTreatment)

samples_per_treatment <- hartwig_clinical %>% filter(group %in% er, hasSystemicPreTreatment != "NULL", hasRadiotherapyPreTreatment != "NULL") %>%
  mutate(pretreatmentType = ifelse(hasSystemicPreTreatment == "Yes", 
                                   ifelse(hasRadiotherapyPreTreatment == "Yes", "Both", "Systemic"), 
                                   ifelse(hasRadiotherapyPreTreatment == "Yes", "Radiotherapy", "Neither")))%>% 
  group_by(group, pretreatmentType) %>% summarize(samples_group = n()) 

samples_with_primary <- rbind(
  primary_megatable %>% filter(group %in% er) %>% select(group) %>% mutate(pretreatmentType = "Primary") %>% 
    group_by(group, pretreatmentType) %>% summarize(samples_group = n()), 
  samples_per_treatment
)

alts_primary <- alterations %>% filter(Stage == "Primary", group %in% er) %>% dplyr::rename(pretreatmentType = Stage)

source_data_2c <- rbind(alts_treatment, alts_primary) %>% 
  mutate(pretreatmentType = factor(pretreatmentType, levels = c("Primary", "Both", "Radiotherapy", "Systemic", "Neither"))) %>% 
  group_by(group, pretreatmentType, Alteration) %>% summarize(n = n()) %>% 
  left_join(samples_with_primary) %>% mutate(alterations_per_sample = n/samples_group) %>% 
  group_by(group, pretreatmentType) %>% mutate(total_alts = sum(alterations_per_sample))

## Save plot ----
pdf(file.path(main_repo_path, "plots", "ExtendedData2C.pdf"), height = 6, width = 6)
ggbarplot(source_data_2c, "pretreatmentType", "alterations_per_sample", 
          color = "Alteration", fill = "Alteration", palette = alteration_colors,
          ylab = "Alterations per sample", xlab = "Pretreatment type", facet.by = "group") +
  geom_text(aes(label = samples_group, y = total_alts*1.07), check_overlap = T) +
  theme_LM +
  rotate_x_text(35) +
  theme(legend.position = "top") +
  guides(fill = guide_legend(ncol = 1), 
         color = guide_legend(ncol = 1)) 
dev.off()


### SAVE ##########################################################################################
source_data_2f %>% 
  select(Biopsy_Site = biopsySite_LM, Alteration, 
         Alterations_per_Sample = alterations_per_sample, 
         Samples_per_Site = samples_site) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData2f_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")

source_data_2c %>% 
  select(Subtype = group, Pretreatment_Type = pretreatmentType, Alteration,
         Alterations_per_Sample = alterations_per_sample, Samples_per_Group = samples_group) %>% 
  write.table(file.path(main_repo_path, "data", "ExtendedData2c_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
