### CREATE QC SUMMARY TABLE ########################################################################
# create qc summary table

### PREAMBLE ######################################################################################
date <- Sys.Date()

setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape')
### REMOVE DUPLICATE DONORS #######################################################################
# read in all experiments 
experiments <- read.csv(
	'isabl_configs/ExportedExperiments_2022-10-19.csv',
	as.is = TRUE
	)

# read in all pcawg samples 
# file downloaded from https://dcc.icgc.org/releases/PCAWG/clinical_and_histology
pcawg <- read.csv(
	'../ICGC/pcawg_specimen_histology_August2016_v9.csv',
	as.is = TRUE
	)
pcawg <- pcawg[grep('BRCA', pcawg$project_code),c('submitted_sample_id','icgc_donor_id','donor_wgs_included_excluded')]
# find pcawg sample ids that overlap with tcga 
pcawgtcga <- pcawg[grep('TCGA', pcawg$submitted_sample_id),]
pcawgtcga$tcga_donor <- substr(pcawgtcga$submitted_sample_id, 1, 12)
all(pcawgtcga$tcga_donor %in% experiments$Individual.ID)

# find experiment ids that overlap 
expids <- experiments[experiments$Individual.ID %in% pcawgtcga$tcga_donor,'Experiment.System.ID']
write.table(
	data.frame(expids),
	file = 'isabl_configs/2022-10-19_overlapping_donors_to_remove.txt',
	row.names = FALSE,
	quote = FALSE,
	col.names = FALSE
	)

### QC ############################################################################################
# find all qc result dirs 
multiqc <- system(
		'isabl get-results --result-key multiqc_stats -fi projects 7 -fi application.pk 1',
		intern = TRUE
		)

# read in and concat 
qc <- do.call(rbind, sapply(
	multiqc,
	function(x) {
		tmp <- read.delim(x, as.is = TRUE)
		analysis <- substr(x, 52, 56)
		print(tmp[1,'Sample'])
		colnames <- c(
			'Sample',
			'QualiMap_mqc.generalstats.qualimap.avg_gc',
			'QualiMap_mqc.generalstats.qualimap.median_coverage',
			'QualiMap_mqc.generalstats.qualimap.mean_coverage',
			'QualiMap_mqc.generalstats.qualimap.percentage_aligned',
			'QualiMap_mqc.generalstats.qualimap.mapped_reads'
			)
		if (all(colnames %in% colnames(tmp))) {
			tmp2 <- tmp[,colnames]
			tmp2$analysis <- analysis
			return(tmp2)
		}
	},
	simplify = FALSE
	))
rownames(qc) <- 1:nrow(qc)

# remove duplicate samples 
qc <- qc[!qc$Sample %in% expids,]

# write to file 
write.table(
	qc,
	file = paste0('qc/', date, '_multiqc_stats.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# reformat 
colnames(qc) <- c('Experiment.System.ID','avg_qc','median_coverage','mean_coverage',
	'percentage_aligned','mapped_reads','analysis')
tumor <- qc[grep('T01', qc$Experiment.System.ID),]
normal <- qc[grep('N01', qc$Experiment.System.ID),]

# add purity
tumor <- merge(tumor, facets[,c('Experiment.System.ID','purity')], by = 'Experiment.System.ID', all = TRUE) 

tumor$Individual.System.ID <- substr(tumor$Experiment.System.ID, 1, 14)
normal$Individual.System.ID <- substr(normal$Experiment.System.ID, 1, 14)

# read in all samples 
map <- read.delim(
	'isabl_configs/2022-10-24_mapping_donor_ids.txt',
	as.is = TRUE
	)

# find sample missing 
missing <- map$Individual.System.ID[which(!map$Individual.System.ID %in% normal$Individual.System.ID)]
normal <- rbind(
	normal,
	data.frame(Experiment.System.ID = paste0(missing, '_N01_01_WG01'), avg_qc = NA, median_coverage = NA,
		mean_coverage = NA, percentage_aligned = NA, mapped_reads = NA, analysis = NA, Individual.System.ID = missing)
	)

qc_rf <- merge(tumor, normal, by = 'Individual.System.ID', all = TRUE)
colnames(qc_rf) <- gsub('.x', '_tumor', colnames(qc_rf), fixed = TRUE)
colnames(qc_rf) <- gsub('.y', '_normal', colnames(qc_rf), fixed = TRUE)
qc_rf <- qc_rf[,grep('analysis', colnames(qc_rf), invert = TRUE)]

### ADD QC FLAGS ##################################################################################
qc_rf$flag <- ''
qc_rf[which(qc_rf$median_coverage_tumor < 25),'flag'] <- 'Median tumor coverage < 25'
qc_rf[which(qc_rf$median_coverage_tumor > 25 & qc_rf$median_coverage_normal < 10),'flag'] <- 'Median normal coverage < 10'

exp_subset <- experiments[,c('Experiment.System.ID','Individual.ID','Sample.ID')]
colnames(exp_subset) <- c('Experiment.System.ID_tumor','Patient_ID','Tumor_Sample_Barcode')
qc_rf <- merge(qc_rf, exp_subset, by = 'Experiment.System.ID_tumor', all.x = TRUE)
qc_rf2 <- qc_rf[,c(1,16,17,9,3:8,10:15)]
colnames(qc_rf2) <- c('sample','Patient_ID','Tumor_Sample_Barcode','sample_normal', colnames(qc_rf2)[-c(1:4)])

# add metadata for one sample
qc_rf2[qc_rf2$sample_normal == 'CURTIS_H000544_N01_01_WG01','sample'] <- 'CURTIS_H000544_T01_01_WG01'
qc_rf2[qc_rf2$sample_normal == 'CURTIS_H000544_N01_01_WG01','Patient_ID'] <- 'DO44103'
qc_rf2[qc_rf2$sample_normal == 'CURTIS_H000544_N01_01_WG01','Tumor_Sample_Barcode'] <- 'SA491668'
qc_rf2[qc_rf2$sample_normal == 'CURTIS_H002081_N01_01_WG01','sample'] <- 'CURTIS_H002081_T01_01_WG01'
qc_rf2[qc_rf2$sample_normal == 'CURTIS_H002081_N01_01_WG01','Patient_ID'] <- 'PD7211'
qc_rf2[qc_rf2$sample_normal == 'CURTIS_H002081_N01_01_WG01','Tumor_Sample_Barcode'] <- 'PD7211a'

# add purity flags
qc_rf2[qc_rf2$sample == 'CURTIS_H001425_T01_01_WG01','flag'] <- 'Cannot estimate purity'
qc_rf2[qc_rf2$sample == 'CURTIS_H001749_T01_01_WG01','flag'] <- 'Cannot estimate purity'
# add alignment and import flags
qc_rf2[qc_rf2$sample == 'CURTIS_H002081_T01_01_WG01','flag'] <- 'Alignment errors'
qc_rf2[qc_rf2$sample == 'CURTIS_H000544_T01_01_WG01','flag'] <- 'Import errors'

write.table(
	qc_rf2,
	file = paste0('data/', date, '_primary_qc_metrics.txt'),
	sep = '\t', 
	row.names = FALSE,
	quote = FALSE 
	)
