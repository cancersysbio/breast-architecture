### SUPPLEMENTARY FIGURE 3ABCD  ##############################################################
# create supplementary figure 3abcd

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
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)

### SUPPLEMENTARY FIGURE 3A #######################################################################
# run pca on cna and sv signatures
input_dtf <- megatable[,grep('^SV32|^CN|group|Sample', colnames(megatable))]
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
sv_prop <- sv_prop/rowSums(sv_prop)
sv_prop[is.na(sv_prop)] <- 0
input_dtf_prop <- cbind(cn_prop, sv_prop)

# rename SV signatures
colnames(input_dtf_prop) <- gsub('SV32A','RS6', colnames(input_dtf_prop))
colnames(input_dtf_prop) <- gsub('SV32D','RS4', colnames(input_dtf_prop))
colnames(input_dtf_prop) <- gsub('SV32B','RS2', colnames(input_dtf_prop))
colnames(input_dtf_prop) <- gsub('SV32F','RS1', colnames(input_dtf_prop))
colnames(input_dtf_prop) <- gsub('SV32C','RS3', colnames(input_dtf_prop))
colnames(input_dtf_prop) <- gsub('SV32E','RS5', colnames(input_dtf_prop))


# run pca on proportions
pca.res <- prcomp(input_dtf_prop, scale = TRUE)
pca <- as.data.frame(pca.res$x)
pca$eniclust <- input_dtf$group
pca$Sample <- input_dtf$Sample


# set colour plalette
total_cols <- c('#fa954eff','#904dbdff','#bd854dff','#a1c1d4ff','#c2b7dbff')
names(total_cols) <- c('ER+ High','IC10','HER2+','ER+ Typical','IC4ER-')

col <- rep('grey', nrow(pca))
for (i in names(total_cols)) {
		col[which(pca$eniclust == i)] <- total_cols[i]
		}

pca$PC2 <- pca$PC2*-1
pca$PC1 <- pca$PC1*-1
# create scatterplot
create.scatterplot(
	PC2 ~ PC1,
	data = pca,
	col = col,
	filename = 'supplementary_figure3a.pdf',
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
		y = 0.99,
		title = 'Subtype',
		corner = c(1,1),
		padding.text = 2
		),
	resolution = 300
	)

### SAVE TO FILE ##################################################################################
# write to file 
write.table(
	input_dtf_prop,
	file = file.path(main_repo_path, 'data', 'SupplementaryFigure3a_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### SUPPLEMENTARY FIGURE 3C #######################################################################
library("factoextra")
# create factor loading map
pdf(paste0(date, '_primary_PC1_PC3_pca_factor_loadings_plot.pdf'), width = 7, height = 7)
fviz_pca_var(pca.res,axes = c(1,3),col.var="contrib")+
 scale_color_gradient2(low="white", mid="grey",
           high="black", midpoint=5) +
 theme_minimal()+ xlim(c(1,-1))
dev.off()

### SUPPLEMENTARY FIGURE 3D #######################################################################
# extract loadings 
loadings <- as.data.frame(pca.res$rotation)
pc1_cn <- loadings[grep('CN', rownames(loadings)),1]
pc1_cn <- pc1_cn[order(-pc1_cn)]
pc2_cn <- loadings[grep('CN', rownames(loadings)),2]
pc2_cn <- pc2_cn[order(-pc2_cn)]

pc1_cn <- pc1_cn[abs(pc1_cn) > 0.15]
pc2_cn <- pc2_cn[abs(pc2_cn) > 0.15]

contrib <- fviz_contrib(pca.res, choice = 'var', axes = c(1,2)) 

pdf(paste0(date, '_primary_pca_factor_barplot.pdf'), width = 7, height = 7)
fviz_contrib(pca.res, choice = 'var', axes = c(1,2), top = 10,  fill = "lightgray", color = "black")+
	theme_minimal()+xlab("Signature")
dev.off()

### SUPPLEMENTARY FIGURE 3B #######################################################################
# run pca on raw activities
pca.res <- prcomp(input_dtf2, scale = TRUE)
pca <- as.data.frame(pca.res$x)
pca$eniclust <- gsub('Low','Typical', input_dtf$group)
pca$Sample <- input_dtf$Sample

# set colour plalette
total_cols <- c('#fa954eff','#904dbdff','#bd854dff','#a1c1d4ff','#c2b7dbff')
names(total_cols) <- c('ER+ High','IC10','HER2+','ER+ Typical','IC4ER-')

col <- rep(NA, nrow(pca))
for (i in names(total_cols)) {
		col[which(pca$eniclust == i)] <- total_cols[i]
		}

pca$PC1 <- pca$PC1*-1
pca$PC2 <- pca$PC2*-1
# create scatterplot
create.scatterplot(
	PC1 ~ PC2,
	data = pca,
	col = col,
	xlimits = c(-7,7),
	ylimits = c(-8,8),
	yat = seq(-5,5,5),
	xat = seq(-10,5,5),
	filename = 'supplementary_figure3b.pdf',
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
		x = 0.98,
		y = 0.98,
		title = 'Subtype',
		corner = c(1,1),
		padding.text = 2
		),
	resolution = 300
	)

### SAVE TO FILE ##################################################################################
# write to file 
write.table(
	pca,
	file = file.path(main_repo_path, 'data', 'SupplementaryFigure3b_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
