### CREATE SUPPLEMENTARY FIGURE 4KL  ##############################################################
# project PCs onto METABRIC 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(randomForest)

set.seed(1234)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### CALCULATE CN PROP #############################################################################
calculate_cn_prop <- function(input_dtf) {
	# calculate proportions
	cn_prop <- input_dtf[,grep('CN', colnames(input_dtf))]
	cn_prop <- cn_prop[,!apply(cn_prop, 2, function(x) {all(x == 0)})]
	cn_prop <- scale(cn_prop/rowSums(cn_prop))
	return(cn_prop)
}

### MAP PCAS ######################################################################################
map_pcas <- function(primarypcas, prim_cn_prop, meta_cn_prop, metabric) {
	# map using random forest
	cn1 <- randomForest(prim_cn_prop, primarypcas$PC1, mtry = 4)
	cn2 <- randomForest(prim_cn_prop, primarypcas$PC2, mtry = 4)
	# here we map the centroid using the trained svm
	cn1.proj <- predict(cn1,meta_cn_prop)
	cn2.proj <- predict(cn2,meta_cn_prop)
	# create data frame
	metapca <- data.frame(
		PC1=cn1.proj,
		PC2=cn2.proj,
		sample = gsub('.', '-', metabric$sample, fixed = TRUE)
		)
	return(metapca)
}

### ANNOTATE WITH ENICLUST ########################################################################
annotate_with_eniclut <- function(pca) {
	eniclust <- read.delim(
		file.path(main_repo_path, 'data', 'metabric_predictions.txt'),
		as.is = TRUE
		)
	eniclust <- eniclust[,c('Sample','voting')]
	colnames(eniclust) <- c('sample','voting')
	clinical <- read.delim(
		'clinical/Metabric_clinical.txt',
		skip = 4,
		as.is = TRUE
		)
	clinical <- clinical[,c('PATIENT_ID','INTCLUST','OS_MONTHS','OS_STATUS','AGE_AT_DIAGNOSIS','ER_STATUS')]
	colnames(clinical) <- c('sample','intclust','time_to_event','event','age','ER')
	pca <- merge(pca, clinical, by = 'sample')
	pca <- merge(pca, eniclust, by = 'sample')
	pca$eniclust <- NA
	pca[which((pca$voting %in% c('ic8','Other')) | (pca$voting == 'ic4' & pca$ER == 'Positive')),'eniclust'] <- 'ER+ Typical'
	pca[which(pca$voting %in% c('ic2','ic6','ic9','ic1')),'eniclust'] <- 'ER+ High'
	pca[which(pca$voting == 'ic5'),'eniclust'] <- 'HER2+'
	pca[which(pca$voting == 'ic10'),'eniclust'] <- 'IC10'
	pca[which(pca$voting == 'ic4' & pca$ER == 'Negative'),'eniclust'] <- 'IC4ER-'
	return(pca)
}

### MAIN ##########################################################################################
# read in primary pcas 
primarypcas <- read.delim(
	file.path(main_repo_path, 'data', 'SupplementaryFigure3a_sourcetable.txt'),
	as.is = TRUE
	)
# read in megatable
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
megatable <- megatable[megatable$Sample %in% primarypcas$Sample,]
prim_cn_prop <- calculate_cn_prop(megatable)
primarypcas <- primarypcas[primarypcas$Sample %in% megatable$Sample,]

# read in metabric 
metabric <- read.delim(
	file.path(main_repo_path, 'data', 'metabric_cn_signatures.txt'),
	as.is = TRUE
	)
meta_cn_prop <- calculate_cn_prop(metabric)

# find cn sig in botn
cnsig <- intersect(colnames(prim_cn_prop), colnames(meta_cn_prop))
prim_cn_prop <- prim_cn_prop[,cnsig]
meta_cn_prop <- meta_cn_prop[,cnsig]

# map primary PCAs to metabric
pca <- map_pcas(primarypcas, prim_cn_prop, meta_cn_prop, metabric)
pca <- annotate_with_eniclut(pca)

# set colour plalette
total_cols <- c('#fa954eff','#904dbdff','#bd854dff','#a1c1d4ff','#c2b7dbff')
names(total_cols) <- c('ER+ High','IC10','HER2+','ER+ Typical','IC4ER-')

pca$eniclust <- gsub("Low","Typical", pca$eniclust)

# create plot
col <- rep(NA, nrow(pca))
for (i in names(total_cols)) {
		col[which(pca$eniclust == i)] <- total_cols[i]
		}

# create scatterplot
create.scatterplot(
	PC2 ~ PC1,
	data = pca,
	col = col,
	xlimits = c(-4.5,5),
	filename = 'supplementary_figure4l.pdf',
	key = list(
		text = list(
			lab = names(total_cols),
			cex = 1,
			col = 'black'
			),
		points = list(
			pch = 19,
			col = total_cols,
			cex = 1
			),
		x = 0.99,
		y = 0.01,
		title = 'Subtype',
		corner = c(1,0),
		padding.text = 2
		),
	resolution = 300
	)

### SUPPLEMENTARY FIGURE 4K #######################################################################
# map primary PCAs to metabric
colnames(megatable) <- gsub('Sample','sample', colnames(megatable))
pca_tcga <- map_pcas(primarypcas, prim_cn_prop, prim_cn_prop, megatable)
colnames(pca_tcga) <- c("PC1","PC2", "sample")
pca_tcga <- merge(pca_tcga, megatable[,c("sample","group")], by = "sample")
colnames(pca_tcga) <- gsub('group','eniclust', colnames(pca_tcga))
pca_tcga$eniclust <- gsub('Low','Typical', pca_tcga$eniclust)

# create plot
col <- rep(NA, nrow(pca_tcga))
for (i in names(total_cols)) {
		col[which(pca_tcga$eniclust == i)] <- total_cols[i]
		}

# create scatterplot
create.scatterplot(
	PC2 ~ PC1,
	data = pca_tcga,
	col = col,
	xlimits = c(-5,5.5),
	ylmits = c(-4,4),
	filename = 'supplementary_figure4k.pdf',
	key = list(
		text = list(
			lab = names(total_cols),
			cex = 1,
			col = 'black'
			),
		points = list(
			pch = 19,
			col = total_cols,
			cex = 1
			),
		x = 0.99,
		y = 0.01,
		title = 'Subtype',
		corner = c(1,0),
		padding.text = 2
		),
	resolution = 300
	)
