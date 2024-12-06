# CREATE sWGS SUMMARY FILE ----
# Create sWGS summary file for TCGA samples
# overlap with WGS

# PREAMBLE ----
library(rstudioapi)
library(data.table)
library(dplyr)
library(yaml)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

# MAIN ----
id_map <- fread(config$map$to_curtis_ids)
ecdna_wgs <- fread(config$aa$proj7) %>%
  mutate(sample_name = substr(sample_name, 1, 14)) %>%
  left_join(id_map, by=c('sample_name' = 'Individual.System.ID')) %>%
  filter(!is.na(TCGA)) %>%
  mutate(type = 'wgs') %>%
  distinct(sample_name, .keep_all = TRUE)

tcga_sample_ids <- ecdna_wgs %>%
  distinct(TCGA) %>%
  pull()

swgs_dir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/ecDNA/TCGA/sWGS/'

dfs <- list()
for (id in tcga_sample_ids) {
  filename <- paste0(swgs_dir, id, '/', id, '_amplicon_classification_profiles.tsv')
  
  if (file.exists(filename)) {
    dfs[[id]] <- fread(filename)
  } else {
    print(paste("File", filename, "does not exist. Skipping..."))
  }
}

ecdna_swgs <- do.call(rbind, dfs) %>%
  mutate(type='swgs') %>%
  distinct(sample_name, .keep_all = TRUE)

# SAVE ----
fwrite(ecdna_swgs, '../data/Summary_TCGA_sWGS.tsv')
