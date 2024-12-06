### CREATE SUPPLEMENTARY FIGURE 4CD ################################################################
# create supplementary figure 4cd
# metastatic pca plot

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(tidyr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# read in megatable
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
	as.is = TRUE
	)

# run pca on cna and sv signatures
input_dtf <- megatable[,grep('^SV32|^CN|group|Sample', colnames(megatable))]
input_dtf <- input_dtf[!apply(input_dtf[,grep('^SV32|^CN', colnames(input_dtf))], 1, function(x) any(is.na(x))),]
input_dtf2 <- input_dtf[,-which(colnames(input_dtf) %in% c('group','Sample'))]
# remove columns that are all zero 
input_dtf2 <- input_dtf2[,!apply(input_dtf2, 2, function(x) {all(x == 0)})]

# calculate proportions
cn_prop <- input_dtf2[,grep('CN', colnames(input_dtf2))]
cn_prop <- cn_prop/rowSums(cn_prop)
sv_prop <- input_dtf2[,grep('SV32', colnames(input_dtf2))]
sv_prop <- sv_prop/rowSums(sv_prop)
input_dtf_prop <- cbind(cn_prop, sv_prop)

# run pca on proportions
pca.res <- prcomp(input_dtf_prop, scale = TRUE)
pca <- as.data.frame(pca.res$x)
pca$eniclust <- gsub('Low','Typical',input_dtf$group)
pca$Sample <- input_dtf$Sample

# set colour plalette
total_cols <- c('#fa954eff','#904dbdff','#bd854dff','#a1c1d4ff','#c2b7dbff')
names(total_cols) <- c('ER+ High','IC10','HER2+','ER+ Typical','IC4ER-')

col <- rep(NA, nrow(pca))
for (i in names(total_cols)) {
		col[which(pca$eniclust == i)] <- total_cols[i]
		}

# create scatterplot
pca$PC2 <- -1*pca$PC2
create.scatterplot(
	PC2 ~ PC1,
	data = pca,
	col = col,
	xlimits = c(-10,10),
	filename = 'supplementary_figure4c.pdf',
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
		y = 0.02,
		title = 'Subtype',
		corner = c(1,0),
		padding.text = 2
		),
	resolution = 300
	)

### SAVE TO FILE ##################################################################################
# write table 
write.table(
	pca,
	file.path(main_repo_path, 'data', 'SupplementaryFigure4c_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### SUPPLEMENTARY FIGURE 4D #######################################################################
activities <- read.delim(
	file.path(main_repo_path, 'data', 'sv_signature_composition_primary_metastatic.txt'),
	as.is = TRUE
	)
# read in megatable
pmtable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
pmtable$disease <- 'primary'
mmtable <- read.delim(
	file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
	as.is = TRUE
	)
mmtable$disease <- 'metastatic'

mtable <- rbind(
	pmtable[,c('Sample','ENiClust','group','disease', 'SVB', grep('^CN', colnames(pmtable), value = TRUE))],
	mmtable[,c('Sample','ENiClust','group','disease', 'SVB', grep('^CN', colnames(pmtable), value = TRUE))]
	)

colnames(activities) <- gsub('Samples','Sample', colnames(activities))
input_dtf <- merge(activities, mtable, by = 'Sample', all.y = TRUE)
input_dtf[which(is.na(input_dtf$SV32A) & input_dtf$SVB <= 6),c('SV32A','SV32B','SV32C','SV32D','SV32E','SV32F')] <- 0

# run pca on cna and sv signatures
input_dtf <- input_dtf[,grep('^SV32|^CN|group|Sample', colnames(input_dtf))]
input_dtf[which(input_dtf$Sample %in% c("CURTIS_H000503_T01_01_WG01","CURTIS_H001469_T01_01_WG01","CURTIS_H001474_T01_01_WG01","CURTIS_H001685_T01_01_WG01",
	"CURTIS_H001748_T01_01_WG01","CURTIS_H002072_T01_01_WG01")),grep('^SV32', colnames(input_dtf))] <- 0
input_dtf <- input_dtf[!apply(input_dtf[,grep('^SV32|^CN', colnames(input_dtf))], 1, function(x) any(is.na(x))),]
input_dtf2 <- input_dtf[,-which(colnames(input_dtf) %in% c('group','Sample'))]
# remove columns that are all zero 
input_dtf2 <- input_dtf2[,!apply(input_dtf2, 2, function(x) {all(x == 0)})]

# calculate proportions
cn_prop <- input_dtf2[,grep('CN', colnames(input_dtf2))]
cn_prop <- cn_prop/rowSums(cn_prop)
sv_prop <- input_dtf2[,grep('SV32', colnames(input_dtf2))]
#sv_prop <- sv_prop[,-which(colnames(sv_prop) == 'SV32F')]
sv_prop <- sv_prop/rowSums(sv_prop)
sv_prop[is.na(sv_prop)] <- 0
#input_dtf_prop <- cn_prop
input_dtf_prop <- cbind(cn_prop, sv_prop)

# run pca on proportions
pca.res <- prcomp(input_dtf_prop, scale = TRUE)
pca <- as.data.frame(pca.res$x)
pca$eniclust <- input_dtf$group
pca$Sample <- input_dtf$Sample
pca$metastatic <- (pca$Sample %in% mmtable$Sample)*1


# set colour plalette
col <- rep('grey', nrow(pca))
col[which(pca$metastatic == 1)] <- '#fe85cfff'
col[which(pca$metastatic == 0)] <- '#9ebff1ff'


#pca$PC2 <- pca$PC2*-1
pca$PC1 <- pca$PC1*-1
# create scatterplot
create.scatterplot(
	PC2 ~ PC1,
	data = pca,
	col = col,
	filename = 'supplementary_figure4d.pdf',
	key = list(
		text = list(
			lab = c('Primary','Metastatic'),
			cex = 1,
			col = 'black'
			),
		points = list(
			pch = 19,
			col = c('#9ebff1ff', '#fe85cfff'),
			cex = 1
			),
		x = 0.01,
		y = 0.01,
		corner = c(0,0),
		padding.text = 2
		),
	resolution = 300
	)

### SAVE TO FILE ##################################################################################
# write to file 
write.table(
	pca,
	file = file.path(main_repo_path, 'data', 'SupplementaryFigure4d_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
