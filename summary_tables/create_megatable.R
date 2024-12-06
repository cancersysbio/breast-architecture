### CREATE MEGA TABLE #############################################################################
# create mega table of all mutations

### FORMAT ENICLUST AND HORMONE ###################################################################
format_eniclust_and_hormone <- function(filter_primary = FALSE) {
	basedir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape'
	# read in all primary samples 
	primary <- read.delim(
		file.path(basedir, 'isabl_configs','2022-10-24_mapping_donor_ids.txt'),
		as.is = TRUE
		)
	# add hartwig primary samples 
	hartwig <- read.delim(
		file.path(basedir, 'clinical', 'Hartwig_clinical.txt'),
		as.is = TRUE
		)
	if (filter_primary) {
		hartwig <- hartwig[which(hartwig$biopsySite == 'Primary'),]
	}
	# read in IC10 labels
	p7 <- read.table(
		'/oak/stanford/groups/ccurtis2/users/lisem/ICGC/icgc_iC10_CNA-only_labels.txt',
		as.is = TRUE
		)
	p17 <- read.table(
		'/oak/stanford/groups/ccurtis2/users/lisem/Hartwig/iC10/hmf_iC10_CNA-only_labels_Isabl.txt',
		as.is = TRUE
		)
	ic10 <- rbind(p7[,c('Sample','iC10')], p17[,c('Sample','iC10')])

	# find hormone 
	icgc <- read.csv(
		file.path(basedir, '../GEB/ICGC/Supplementary Table 1 CLINICAL.PATHOLOGY.DATA.FREEZE.ANALYSIS.v4.032015.csv'),
		as.is = TRUE
		)
	icgc <- icgc[,c('sample_name','final.ER','final.HER2','final.PR')]
	colnames(icgc) <- c('NikZainal','ER','HER2','PR')
	icgc[icgc == 'negative'] <- 0
	icgc[icgc == 'positive'] <- 1
	icgc <- merge(primary[,c('Individual.System.ID','NikZainal')], icgc, by = 'NikZainal')

	# read in clinical status tcga 
	tcga_ihc <- read.csv(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/GEB/IBC/TCGA/germline/epitope_prediction/Thennavan_brca_annotations.csv',
		as.is = TRUE
		)
	tcga_ihc <- tcga_ihc[tcga_ihc$Tumor.or.Normal == 'Tumor',c('CLID','er_status_by_ihc')]
	tcga_ihc$TCGA <- substr(tcga_ihc$CLID, 1, 12)
	tcga_ihc <- tcga_ihc[grep('-06A', tcga_ihc$CLID, invert = TRUE),]

	load(file.path(basedir, '../GEB/IBC/TCGA/genefuAnnotations.RData'))
	tcga <- BRCA_annotation[,c('ID','ER','HER2','PR')]
	tcga$ID <- substr(tcga$ID, 1, 12)
	colnames(tcga) <- c('TCGA','ER_rna','HER2','PR')
	tcga <- merge(tcga, tcga_ihc, by = 'TCGA')
	tcga$ER <- NA
	tcga[which(tcga$er_status_by_ihc == 'Positive'),'ER'] <- 1
	tcga[which(tcga$er_status_by_ihc == 'Negative'),'ER'] <- 0
	tcga[which(tcga$er_status_by_ihc %in% c('[Not Evaluated]','Indeterminate') & tcga$ER_rna == 1),'ER'] <- 1
	tcga[which(tcga$er_status_by_ihc %in% c('[Not Evaluated]','Indeterminate') & tcga$ER_rna == 0),'ER'] <- 0
	tcga <- merge(primary[,c('Individual.System.ID','TCGA')], tcga[,c('TCGA','ER','HER2','PR')], by = 'TCGA')

	hormones <- merge(primary, rbind(tcga[,-1], icgc[,-1]), by = 'Individual.System.ID', all.x = TRUE)
	hormones$Individual.System.ID <- paste0(hormones$Individual.System.ID, '_T01_01_WG01')
	hormones$TNBC <- (hormones$ER == 0 & hormones$PR == 0 & hormones$HER2 == 0)*1
	hormones$Hartwig <- NA
	# add hartwig
	hartwig_hormones <- hartwig[,c('Individual.System.ID','sampleId','cancerSubtype')]
	colnames(hartwig_hormones) <- gsub('sampleId','Hartwig', colnames(hartwig_hormones))
	hartwig_hormones$ER <- (hartwig_hormones$cancerSubtype %in% c('ER-positive/HER2 unknown','ER-positive/HER2-negative','ER-positive/HER2-positive'))*1
	hartwig_hormones$HER2 <- (hartwig_hormones$cancerSubtype %in% c('ER-negative/HER2-positive','ER-positive/HER2-positive'))*1
	hartwig_hormones$TNBC <- (hartwig_hormones$cancerSubtype == 'Triple negative')*1
	hartwig_hormones[hartwig_hormones$cancerSubtype %in% c('', 'Subtype unknown'),'ER'] <- NA
	hartwig_hormones[hartwig_hormones$cancerSubtype %in% c('', 'Subtype unknown'),'HER2'] <- NA
	hartwig_hormones[hartwig_hormones$cancerSubtype %in% c('', 'Subtype unknown'),'TNBC'] <- NA
	hartwig_hormones$PR <- NA
	hartwig_hormones$TCGA <- NA
	hartwig_hormones$PCAWG <- NA
	hartwig_hormones$NikZainal <- NA
	hormones <- rbind(hormones, hartwig_hormones[,colnames(hormones)])
	colnames(hormones) <- gsub('Individual.System.ID','Sample', colnames(hormones))

	hormones <- merge(hormones, ic10, by = 'Sample', all.x = TRUE)

	# add eniclust
	icgc_eniclust <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v14/Pipeline_2/icgc_predictions.txt',
		#'/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/ic_subtypes/data/eniclust/icgc_predictions.txt',
		as.is = TRUE
		)
	hartwig_eniclust <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/lisem/eniclust/wes_wgs/dev/test/rand/rand_59_cv_noLOW_v14/Pipeline_2/hmf_predictions.txt',
		#'/oak/stanford/groups/ccurtis2/users/lisem/bc_landscape/ic_subtypes/data/eniclust/hmf_predictions.txt',
		as.is = TRUE
		)
	eniclust <- rbind(icgc_eniclust[,c('Sample','voting')], hartwig_eniclust[,c('Sample','voting')])
	colnames(eniclust) <- c('Sample','ENiClust')
	hormones <- merge(hormones, eniclust, by = 'Sample',all.x = TRUE)

	hormones <- hormones[,c("Sample","iC10","ENiClust","ER","HER2","PR","TNBC")]
	hormones <- unique(hormones)
	return(hormones)
}

### FORMAT PLOIDY, WGD AND PGA ####################################################################
format_ploidy_wgd_pga <- function() {
	# read in facets best fit 
	pfacets <- read.delim(
		'/oak/stanford/groups/ccurtis2/isabl/analyses/65/27/16527/bestfit_merged.tsv',
		as.is = TRUE
		)
	pfacets <- pfacets[,c('sample','purity','ploidy','genome_doubled','fraction_cna', 'fraction_loh')]
	mfacets <- read.delim(
		'/oak/stanford/groups/ccurtis2/isabl/analyses/32/93/23293/bestfit_merged.tsv',
		as.is = TRUE
		)
	mfacets <- mfacets[,c('sample','purity','ploidy','genome_doubled','fraction_cna', 'fraction_loh')]
	facets <- rbind(pfacets, mfacets)
	colnames(facets) <- gsub('sample','Sample', colnames(facets))
	return(facets)
}

### FORMAT CNA SIGNATURES #########################################################################
format_cna_signatures <- function() {
	pcnasig <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/cnas/signatures/primary/COSMIC/Assignment_Solution/Activities/Assignment_Solution_Activities.txt',
		as.is = TRUE
		)
	mcnasig <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/cnas/signatures/metastatic/COSMIC/Assignment_Solution/Activities/Assignment_Solution_Activities.txt',,
		as.is = TRUE
		)
	mcnasig$Samples <- gsub('_hisens','', mcnasig$Samples)
	cnasig <- rbind(pcnasig, mcnasig)
	colnames(cnasig) <- gsub("Samples", "Sample", colnames(cnasig))
	return(cnasig)
}

### FORMAT SV SIGNATURES ##########################################################################
format_primary_sv_signatures <- function() {
	svsig <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/svs/sv_signatures/primary/clustered_vs_nonclustered_primary_only_updated_nolowqual/SV32/All_Solutions/SV32_6_Signatures/Activities/SV32_S6_NMF_Activities.txt',
		as.is = TRUE
		)
	colnames(svsig) <- gsub("Samples", "Sample", colnames(svsig))
	return(svsig)
}
format_metastatic_sv_signatures <- function() {
	svsig <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/svs/sv_signatures/metastatic/clustered_vs_nonclustered_metastatic_only_updated/SV32/All_Solutions/SV32_6_Signatures/Activities/SV32_S6_NMF_Activities.txt',
		as.is = TRUE
		)
	colnames(svsig) <- gsub("Samples", "Sample", colnames(svsig))
	return(svsig)
}

### FORMAT SV SUMMARY #############################################################################
format_sv_summary <- function() {
	ptotal_sv <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/svb_project7.txt',
		as.is = TRUE
		)
	mtotal_sv <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/svb_project17.txt',
		as.is = TRUE
		)
	total_sv <- rbind(ptotal_sv, mtotal_sv)
	colnames(total_sv) <- gsub('sample',"Sample", colnames(total_sv))
	pdamaging_sv <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/svb_damaging_project7.txt',
		as.is = TRUE
		)
	mdamaging_sv <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/svb_damaging_project17.txt',
		as.is = TRUE
		)
	damaging_sv <- rbind(pdamaging_sv, mdamaging_sv)	
	colnames(damaging_sv) <- gsub('sample',"Sample", colnames(damaging_sv))
	psv_clustering <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/2023-01-17_primary_sv_clusters_summary.txt',
		as.is = TRUE
		)
	msv_clustering <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/2023-01-23_metastatic_sv_clusters_summary.txt',
		as.is = TRUE
		)
	sv_clustering <- rbind(psv_clustering, msv_clustering)
	colnames(sv_clustering) <- gsub('Experiment.System.ID','Sample', colnames(sv_clustering))
	# merge 
	sv <- merge(total_sv[,c('Sample','SVB')], damaging_sv[,c('Sample','SV_dama','SV_dama_prop')], 
		by = 'Sample', all = TRUE)
	sv <- merge(sv, sv_clustering[,-2], by = 'Sample', all = TRUE)
	# add 0s to samples with 0 SVs
	zero_svs <- c('CURTIS_H000503_T01_01_WG01','CURTIS_H001469_T01_01_WG01','CURTIS_H001474_T01_01_WG01',
		'CURTIS_H001748_T01_01_WG01','CURTIS_H002072_T01_01_WG01','CURTIS_H001685_T01_01_WG01')
	sv[sv$Sample %in% zero_svs,'SVB'] <- 0
	sv[sv$Sample %in% zero_svs,'SV_dama'] <- 0
	sv[sv$Sample %in% zero_svs,'SV_dama_prop'] <- 0
	sv[sv$Sample %in% zero_svs,'number_sv_clusters'] <- 0
	sv[sv$Sample %in% zero_svs,'prop_sv_in_clusters'] <- 0
	sv[sv$Sample %in% zero_svs,'median_num_sv_in_clusters'] <- 0
	return(sv)
}

### FORMAT CHROMOTHRIPSIS #########################################################################
format_chromothripsis <- function() {
	# read in chromothripsis 
	pcmt <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/chromothripsis_project7.txt',
		as.is = TRUE
		)
	mcmt <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/chromothripsis_project17.txt',
		as.is = TRUE
		)
	cmt <- rbind(pcmt, mcmt)
	colnames(cmt) <- gsub('sample', 'Sample', colnames(cmt))
	colnames(cmt) <- gsub('confidence', 'confidence_chromothripsis', colnames(cmt))
	return(cmt[,c('Sample','high_confidence_chromothripsis','low_confidence_chromothripsis')])
}

### FORMAT ECDNA SUMMARY ##########################################################################
format_ecdna_summary <- function() {
	pecdna <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/2023-01-17_primary_ecdna_summary.txt',
		as.is = TRUE
		)
	pecdna <- pecdna[,-1]
	mecdna <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/2023-01-17_metastatic_ecdna_summary.txt',
		as.is = TRUE
		)
	mecdna <- mecdna[,-1]
	ecdna <- rbind(pecdna, mecdna)
	colnames(ecdna) <- c('Sample','AA_amplicons','AA_ecdna','AA_complex')
	return(ecdna)
}

### FORMAT HLA ####################################################################################
format_hla <- function(qc) {
	phla <- read.csv(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/2023-01-18_primary_hla_summary.txt', 
		as.is = TRUE
		)
	phla <- phla[phla$Sample %in% qc$sample_normal,-c(7:8)]
	colnames(phla) <- gsub("Sample", "sample_normal", colnames(phla))
	phla <- merge(phla, qc[,c("sample_normal", "sample")], by = 'sample_normal')
	phla <- phla[,-1]
	colnames(phla) <- gsub("sample","Sample", colnames(phla))
	mhla <- read.csv(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data/2023-01-29_met_hla_summary.txt', 
		as.is = TRUE
		)
	mhla$Sample <- gsub('N01','T01',mhla$Sample)
	hla <- rbind(phla, mhla[,colnames(phla)])
	return(hla)
}

### ADD FOUR MAIN GROUPS ##########################################################################
add_four_main_groups <- function(megatable) {
	# set four main groups
	megatable$group <- NA
	megatable[which(megatable$ENiClust %in% c('ic1','ic2','ic6','ic9')),'group'] <- 'ER+ High'
	megatable[which(megatable$ENiClust == 'ic5'),'group'] <- 'HER2+'
	megatable[which(megatable$ENiClust == 'ic10'),'group'] <- 'IC10'
	megatable[which(megatable$ENiClust == 'ic4' & megatable$ER == 0),'group'] <- 'IC4ER-'
	megatable[which(megatable$ENiClust %in% c('Other','ic8') | (megatable$ENiClust == 'ic4' & megatable$ER == 1)),'group'] <- 'ER+ Typical'
	return(megatable)
}

### FILTER PATIENTS WITHOUT ENOUGH INFORMED CONSENT ###############################################
filter_for_consent <- function(megatable, hartwig) {
	# read in list of patients without enough informated consent
	c1 <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/GEB/hartwig/HMFIDs_insufficient_IC_06-12-22.txt',
		header = FALSE,
		as.is = TRUE
		)
	c2 <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/GEB/hartwig/HMFIDs_insufficient_IC_15-06-22.txt',
		as.is = TRUE
		)
	c3 <- read.delim(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/GEB/hartwig/samples_to_remove_insufficient_IC.txt',
		header = FALSE,
		as.is = TRUE
		)
	# find system ids 
	c2_sid <- hartwig[hartwig$hmfPatientId %in% c2$hmfPatientId,'Individual.System.ID']
	c1_sid <- hartwig[hartwig$hmfSampleId %in% c1$V1,'Individual.System.ID']
	c3_sid <- hartwig[hartwig$X.patientId %in% c3$V1,'Individual.System.ID']
	# remove samples 
	megatable <- megatable[!megatable$Sample %in% c(c2_sid, c1_sid, c3_sid),]
	return(megatable)
}

### ADD COHORT ####################################################################################
add_cohort <- function(megatable) {
	basedir <- '/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape'
	# read in all primary samples 
	primary <- read.delim(
		file.path(basedir, 'isabl_configs','2022-10-24_mapping_donor_ids.txt'),
		as.is = TRUE
		)
	pcawg <- paste0(primary[!is.na(primary$PCAWG),'Individual.System.ID'], '_T01_01_WG01')
	nz <- paste0(primary[is.na(primary$PCAWG) & !is.na(primary$NikZainal),'Individual.System.ID'], '_T01_01_WG01')
	tcga <- paste0(primary[is.na(primary$PCAWG) & !is.na(primary$TCGA),'Individual.System.ID'], '_T01_01_WG01')
	# add hartwig primary samples 
	hartwig <- read.delim(
		file.path(basedir, 'clinical', 'Hartwig_clinical.txt'),
		as.is = TRUE
		)
	# add cohort
	megatable$cohort <- NA
	megatable[megatable$Sample %in% hartwig$Individual.System.ID,'cohort'] <- "Hartwig"
	megatable[megatable$Sample %in% pcawg,'cohort'] <- "PCAWG"
	megatable[megatable$Sample %in% nz,'cohort'] <- "NikZainal"
	megatable[megatable$Sample %in% tcga,'cohort'] <- "TCGA"
	return(megatable)
}


### MAIN ##########################################################################################
setwd('/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape')
# read in qc table to start 
qc <- read.delim(
	'data/2022-11-15_primary_qc_metrics.txt',
	as.is = TRUE
	)
# add hartwig primary samples 
hartwig <- read.delim(
		file.path('clinical', 'Hartwig_clinical.txt'),
		as.is = TRUE
		)
phartwig <- hartwig[which(hartwig$biopsySite == 'Primary'),]
mhartwig <- hartwig[which(hartwig$biopsySite != 'Primary'),]

# read in eniclust and hormones 
hormones <- format_eniclust_and_hormone(filter_primary = TRUE)
mhormones <- format_eniclust_and_hormone()

# read in ploidy, wgd and pga
facets <- format_ploidy_wgd_pga()

# read in cna signatures
cnasig <- format_cna_signatures() 

# read in sv signatures 
svsig <- format_primary_sv_signatures()
msvsig <- format_metastatic_sv_signatures()

# read in sv summary
sv <- format_sv_summary()

# read in chromothripsis
chromothripsis <- format_chromothripsis() 

# read in ecdna summary
ecdna <- format_ecdna_summary()

# read in hla alleles
hla <- format_hla(qc)

# merge into megatable 
mtable <-  Reduce(
	function(dtf1, dtf2) merge(dtf1, dtf2, by = "Sample", all = TRUE),
	list(hormones, facets, cnasig, svsig, sv, chromothripsis, ecdna, hla)
	)
# only keep samples that pass qc 
mtable <- unique(mtable[mtable$Sample %in% c(qc[qc$flag == '','sample'], phartwig$Individual.System.ID),])
# CURTIS_H000525_T01_01_WG01 included twice due to two samples in TCGA, only difference is PR status so keeping the PR = 1 
mtable <- mtable[!duplicated(mtable$Sample),]
# remove two pcawg samples that are mets and one hartiwg sample (CURTIS_H004312_T01_01_WG01) without enough consent
mtable <- mtable[!mtable$Sample %in% c('CURTIS_H000508_T01_01_WG01','CURTIS_H000525_T01_01_WG01','CURTIS_H004312_T01_01_WG01'),]
# add subtype groups
mtable <- add_four_main_groups(mtable)

# MAKE SURE TO REMOVE NON CONSENT PATIENTS
mtable <- filter_for_consent(mtable, hartwig)

# remove samples with high number of inversions suggesting low quality - 7 samples
lowqual <- read.delim(
	'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/qc/low_qual_removed_sv_signature.txt',
	as.is = TRUE,
	header = FALSE
	)
mtable <- mtable[!mtable$Sample %in% lowqual$V1,]

# add cohort 
mtable <- add_cohort(mtable)

# set SV signatures to 0 if too few SVs
mtable[which(is.na(mtable$SV32A) & mtable$SVB < 10),c('SV32A','SV32B','SV32C','SV32D','SV32E','SV32F')] <- 0
# remove duplicated donor coming from hartwig CURTIS_H004115_T02_01_WG01 (chosing sample with higher purity)
mtable <- mtable[mtable$Sample != 'CURTIS_H004115_T02_01_WG01',]


# write to table
write.table(
	mtable,
	file = file.path(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data',
		paste0(date, '_primary_megatable.txt')
		),
	sep = '\t',
	quote = FALSE,
	row.names = FALSE
	)


# merge into megatable 
met_mtable <-  Reduce(
	function(dtf1, dtf2) merge(dtf1, dtf2, by = "Sample", all = TRUE),
	list(mhormones, facets, cnasig, msvsig, sv, chromothripsis, ecdna, hla)
	)
# only keep samples that pass qc 
met_mtable <- unique(met_mtable[which(met_mtable$Sample %in% c(mhartwig$Individual.System.ID,'CURTIS_H000508_T01_01_WG01',
	'CURTIS_H000525_T01_01_WG01')),])
met_mtable <- met_mtable[!duplicated(met_mtable$Sample),]
# add subtype groups
met_mtable <- add_four_main_groups(met_mtable)
# MAKE SURE TO REMOVE NON CONSENT PATIENTS
met_mtable <- filter_for_consent(met_mtable, hartwig)
# add cohort
met_mtable$cohort <- 'Hartwig'
met_mtable[met_mtable$Sample %in% c('CURTIS_H000508_T01_01_WG01','CURTIS_H000525_T01_01_WG01'),'cohort'] <- 'PCAWG'
# remove duplicated samples from the same donor
donor_ids <- substr(met_mtable$Sample, 1, 14)
unique_donors <- unique(donor_ids)
duplicate_donors <- unique(donor_ids[duplicated(donor_ids)])
# find sample with highest purity
toremove <- c()
for (i in duplicate_donors) {
	tmp <- met_mtable[grep(i, met_mtable$Sample),]
	tokeep <- tmp$Sample[which(tmp$purity == max(tmp$purity))]
	toremove <- c(toremove, tmp$Sample[!tmp$Sample %in% tokeep])
}
met_mtable2 <- met_mtable[!met_mtable$Sample %in% toremove,]

# write to table
write.table(
	met_mtable2,
	file = file.path(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data',
		paste0(date, '_metastatic_megatable.txt')
		),
	sep = '\t',
	quote = FALSE,
	row.names = FALSE
	)

