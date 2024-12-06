### CREATE ECDNA SUMMARY TABLE ####################################################################
# create ecdna summary table

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(argparse)


setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape')
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-p', '--project', type = 'character', help = 'project, options are primary or metastatic');

args <- parser$parse_args();
### FIND SAMPLES WITH NO SEEDS OR NO AMPLICONS ####################################################
find_samples_no_seeds_or_amplicons <- function(summary) {
	# check for seeds 
	seeds <- do.call(rbind, sapply(
		summary,
		function(x) {
			analysis <- strsplit(x, '\t')[[1]][1]
			print(analysis)
			if (nchar(analysis) == 4) {
				a1 <- substr(analysis, 1, 2)
				a2 <- substr(analysis, 3, 4)
			} else if (nchar(analysis) == 5) {
				a1 <- substr(analysis, 2, 3)
				a2 <- substr(analysis, 4, 5)
			}
			dir <- file.path(
				'/oak/stanford/groups/ccurtis2/isabl/analyses', a1, a2, analysis
				)
			bed <- list.files(path = dir, pattern = 'WG01.bed', full.name = TRUE)
			size <- file.info(bed)$size
			if (size == 0) {
				data.frame(
					analysis = analysis,
					sample = gsub('.bed', '', basename(bed)),
					size = size,
					amplicons = NA
					)
			} else {
				summaryfile <- list.files(path = dirname(bed), pattern = 'summary', recursive = TRUE, full.name = TRUE)
				if (file.info(summaryfile)$size > 0) {
					data.frame(
						analysis = analysis,
						sample = gsub('.bed', '', basename(bed)),
						size = size,
						amplicons = system(paste('head -n 1', summaryfile), intern = TRUE)
						)
				} else {
					data.frame(
						analysis = analysis,
						sample = gsub('.bed', '', basename(bed)),
						size = size,
						amplicons = NA
						)
					}
				}
			},
		simplify = FALSE
		))
	# find samples with no seeds 
	noseeds <- seeds[which(seeds$size == 0),'sample']
	# find samples with no amplicons 
	noamplicons <- seeds[which(seeds$amplicons == '#Amplicons = 0'),'sample']
	# create output list 
	output <- list()
	output[['noseeds']] <- as.character(noseeds)
	output[['noamplicons']] <- as.character(noamplicons)
	return(output)
}

### CREATE SUMMARY OF ECDNA #######################################################################
create_summary_of_ecdna_profile <- function(ecdna) {
	# create a file with number of amplicons, number of ecDNA, number of complex linear
	amp <- do.call(rbind, sapply(
		unique(ecdna$sample_name),
		function(x) {
			tmp <- ecdna[ecdna$sample_name == x,]
			tmp <- tmp[tmp$amplicon_decomposition_class != 'No amp/Invalid',]
			data.frame(
				Individual.System.ID = as.character(x),
				Experiment.System.ID = as.character(x),
				amp = nrow(tmp),
				ecdna = sum(tmp$ecDNA. == 'Positive'),
				complex = sum(tmp$amplicon_decomposition_class == 'Complex non-cyclic')
				)
			},
		simplify = FALSE
		))
	amp$Individual.System.ID <- substr(amp$Individual.System.ID, 1, 14)
	noamp <- data.frame(
		Individual.System.ID = as.character(substr(c(noseeds, noamplicons),1,14)),
		Experiment.System.ID = as.character(c(noseeds, noamplicons)),
		amp = 0,
		ecdna = 0,
		complex = 0
		)
	amp <- rbind(amp, noamp)
	return(amp)
}

### MAIN ##########################################################################################
# set project number
if (args$project == 'primary') {
	project_num <- 7
	ecdna_file <- '/oak/stanford/groups/ccurtis2/isabl/analyses/48/40/4840/classification_profiles_merged.tsv'
	gene_file <- '/oak/stanford/groups/ccurtis2/isabl/analyses/48/40/4840/gene_list_merged.tsv'
	cna_file <- '/oak/stanford/groups/ccurtis2/isabl/analyses/65/27/16527/gene_level_merged.tsv'
} else if (args$project == 'metastatic') {
	project_num <- 17
	ecdna_file <- '/oak/stanford/groups/ccurtis2/isabl/analyses/44/33/34433/classification_profiles_merged.tsv'
	gene_file <- '/oak/stanford/groups/ccurtis2/isabl/analyses/44/33/34433/gene_list_merged.tsv'
	cn_file <- '/oak/stanford/groups/ccurtis2/isabl/analyses/32/93/23293/gene_level_merged.tsv'
	}

# find all summary files
summary <- system(
		paste('isabl get-results --result-key summary -fi projects', project_num, '-fi application.pk 84 --verbose'),
		intern = TRUE
		)
# find samples with no seeds or amplicons
seeds <- find_samples_no_seeds_or_amplicons(summary)
noseeds <- seeds$noseeds
noamplicons <- seeds$noamplicons

# read in ecDNA summary files 
ecdna <- read.delim(
	ecdna_file,
	as.is = TRUE
	)
ecdna$Individual.ID <- substr(ecdna$sample_name, 1, 14)

# create summary of ecdna profile 
amp <- create_summary_of_ecdna_profile(ecdna)

# write to file 
write.table(
	amp,
	file = file.path('data', paste0(date, '_', args$project, '_ecdna_summary.txt')),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
