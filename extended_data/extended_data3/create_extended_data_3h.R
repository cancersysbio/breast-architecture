### CREATE EXTENDED DATA 3H #######################################################################
# create barplot showing ratio of primary vs metastatic sv signatur activity
### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### ACTIVITIES COMPARISON #########################################################################
activities_file <- file.path(main_repo_path, 'data', 'sv_signature_primary_metastatic_activities.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put sv signature activities for both primary and metastatic tumours in data directory ...")
}
activities <- read.delim(
	activities_file,
	as.is = TRUE
	)
# read in megatable
pmtable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
pmtable <- pmtable[which(pmtable$group != ''),]
pmtable$disease <- 'primary'
mmtable <- read.delim(
	file.path(main_repo_path, 'data', 'metastatic_megatable.txt'),
	as.is = TRUE
	)
mmtable <- mmtable[which(mmtable$group != ''),]
mmtable$disease <- 'metastatic'

mtable <- rbind(
	pmtable[,c('Sample','ENiClust','group','disease')],
	mmtable[,c('Sample','ENiClust','group','disease')]
	)

colnames(activities) <- gsub('Samples','Sample', colnames(activities))
input_dtf <- merge(activities, mtable, by = 'Sample', all.y = TRUE)

colnames(input_dtf) <- gsub('SV32A','RS6', colnames(input_dtf))
colnames(input_dtf) <- gsub('SV32B','RS2', colnames(input_dtf))
colnames(input_dtf) <- gsub('SV32C','RS4', colnames(input_dtf))
colnames(input_dtf) <- gsub('SV32D','RS5', colnames(input_dtf))
colnames(input_dtf) <- gsub('SV32E','RS3', colnames(input_dtf))
colnames(input_dtf) <- gsub('SV32F','RS1', colnames(input_dtf))
# # map tumor and met sv signatures 

pplot_data <- gather(
	input_dtf[input_dtf$disease == 'primary',c('group','RS1','RS2','RS3','RS4','RS5','RS6')],
	value = 'SV',
	key = 'Sig',
	-group
	)
pplot_data <- pplot_data[!is.na(pplot_data$group),]

mplot_data <- gather(
	input_dtf[input_dtf$disease == 'metastatic',c('group','RS1','RS2','RS3','RS4','RS5','RS6')],
	value = 'SV',
	key = 'Sig',
	-group
	)
mplot_data <- mplot_data[!is.na(mplot_data$group),]

pmedian_data <- aggregate(pplot_data$SV, pplot_data[,c('group','Sig')], median, na.rm = TRUE)
mmedian_data <- aggregate(mplot_data$SV, mplot_data[,c('group','Sig')], median, na.rm = TRUE)
median_data <- merge(pmedian_data, mmedian_data, by = c('group','Sig'))
colnames(median_data) <- c('group','Sig','primary','metastatic')
median_data$ratio <- (median_data$primary+1)/median_data$metastatic
median_data$logratio <- log(median_data$ratio) 
median_data$index <- paste(median_data$group, median_data$Sig, sep = '_')

res <- list()
for (i in unique(pplot_data$group[!is.na(pplot_data$group)])) {
	for (j in unique(pplot_data$Sig)) {
		tmp <- wilcox.test(
			pplot_data[which(pplot_data$group == i & pplot_data$Sig == j),'SV'],
			mplot_data[which(mplot_data$group == i & mplot_data$Sig == j),'SV']
			)
		res[[paste(i, j)]] <- data.frame(
			group = i,
			sig = j,
			p = tmp$p.value
			)
	}
}
res <- do.call(rbind, res)
res$fdr <- p.adjust(res$p, method = 'fdr')

res <- res[order(as.character(res$group), as.character(res$sig)),]

col <- rep('black', nrow(median_data))
col[which(res$fdr < 0.05)] <- 'firebrick3'

create.barplot(
	logratio ~ index,
	data = median_data,
	ylimits = c(-5.5,5.5),
	col = col,
	filename = 'extended_data3h.pdf',
	xaxis.lab = rep(c('RS1','RS2','RS3','RS4','RS5','RS6'), 5),
	xaxis.cex = 1,
	xaxis.rot = 90,
	add.rectangle = TRUE,
	xleft.rectangle = c(0,12.5, 24.5),
	#xleft.rectangle = c(0,6.5,12.5),
	ybottom.rectangle = -5.5,
	xright.rectangle = c(6.5,18.5, 32),
	#xright.rectangle = c(3.5,9.5,16),
	ytop.rectangle = 5.5,
	col.rectangle = 'grey50',
	alpha.rectangle = 0.5,
	add.text = TRUE,
	ylab.label = 'log Ratio',
	xlab.label = 'SV Signature',
	text.labels = c('ER+ High','ER+ Typical','HER2+','IC10','IC4ER-'),
	text.x = seq(3.5,32,6),
	text.y = rep(5.2,5),
	height = 6,
	resolution = 400
	)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  median_data,
  file = file.path(main_repo_path, 'data', 'ExtendedData3h_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )
