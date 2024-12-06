### CREATE EXTENDED DATA 3H #######################################################################
# create barplot of co-amplified genes

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

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

# read in amplicons 
amplicons_file <- file.path(main_repo_path, 'data', 'amplicon_segments.txt')
if (!file.exists(amplicons_file)) {
	stop("Please put amplicon segments as outputted by AA in data directory ...")
}
amplicons <- read.delim(amplicons_files, as.is = TRUE)
pamplicons <- amplicons[which(amplicons$sample %in% megatable$sample),]
ecdna_all <- pamplicons[pamplicons$ecDNA == 'Positive',]

# find oncogenes of interest
ecdna_all$FGFR1 <- grepl('FGFR1', ecdna_all$oncogenes)*1
ecdna_all$CCND1 <- grepl('CCND1', ecdna_all$oncogenes)*1
ecdna_all$RPS6KB1 <- grepl('RPS6KB1', ecdna_all$oncogenes)*1
ecdna_all$MYC <- grepl('MYC', ecdna_all$oncogenes)*1
ecdna_all$ERBB2 <- grepl('ERBB2', ecdna_all$oncogenes)*1

ecdnas <- unique(ecdna_all[,c('sample','ampliconID','FGFR1','CCND1','RPS6KB1','MYC','ERBB2')])
num_ecdnas <- as.data.frame(table(ecdnas$sample))
num_ecdnas <- num_ecdnas[order(-num_ecdnas$Freq),]

ecdnas_coamp <- ecdnas[which(rowSums(ecdnas[,-c(1,2)]) > 1),]
ecdnas_coamp_df <- as.data.frame(table(ecdnas_coamp[,-c(1,2)]))
ecdnas_coamp_df <- ecdnas_coamp_df[ecdnas_coamp_df$Freq != 0,]
ecdnas_coamp_df$genes <- apply(
	ecdnas_coamp_df,
	1,
	function(x) {
		i <- which(x[1:5] != 0)
		paste(names(x)[i], collapse = '|')
		})
ecdnas_coamp_df <- ecdnas_coamp_df[order(-ecdnas_coamp_df$Freq),]
ecdnas_coamp_df$index <- 1:nrow(ecdnas_coamp_df)

create.barplot(
	Freq ~ index,
	data = ecdnas_coamp_df,
	ylimits = c(0,14),
	yat = seq(0,12,4),
	ylab.label = 'Number of ecDNA\nharboring \u2265 2 oncogenes',
	filename = 'extended_data4h.pdf',
	xlab.label = 'Genes Co-amplified',
	xaxis.lab = ecdnas_coamp_df$genes,
	xaxis.fontface = 4,
	xaxis.rot = 45,
	xaxis.cex = 0.8,
	resolution = 300
	)

### SAVE ##########################################################################################
# write source data to file 
write.table(
  ecdnas_coamp_df,
  file = file.path(main_repo_path, 'data', 'ExtendedData4h_sourcetable.txt'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE
  )