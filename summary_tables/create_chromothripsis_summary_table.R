### SUMMARIZE CHROMOTHRIPSIS ######################################################################
# summarize output from chromothripsis

### PREAMBLE ######################################################################################
# set date
date <- Sys.Date()

### CATEGORIZE CHROMOTHRIPSIS CONFIDENCE ##########################################################
# confidence categorizations from shatterseek tutorial 
categorize_chromothripsis_confidence <- function(dtf) {
	SVsix <- rowSums(dtf[,c('number_DEL','number_DUP','number_h2hINV','number_t2tINV')]) >= 6
	SVthree <- rowSums(dtf[,c('number_DEL','number_DUP','number_h2hINV','number_t2tINV')]) >= 3
	SVinter <- dtf$number_TRA >= 4
	CNseven <- dtf$max_number_oscillating_CN_segments_2_states >= 7
	CNfour <- dtf$max_number_oscillating_CN_segments_2_states >= 4
	fragp <- dtf$pval_fragment_joins < 0.05
	chrp <- dtf$pval_exp_chr < 0.05
	expp <- dtf$pval_exp_cluster < 0.05

	dtf$confidence <- NA
	dtf[which(SVsix & CNseven & fragp & (chrp | expp)),'confidence'] <- 'High'
	dtf[which(SVthree & SVinter & CNseven & fragp),'confidence'] <- 'High'
	dtf[which(SVsix & CNfour & fragp & (chrp | expp)),'confidence'] <- 'Low'
	return(dtf)
}

summarize_regions <- function(dtf) {
	if (nrow(dtf) > 0) {
		regions <- paste0(dtf$chrom, ':', dtf$start, '-', dtf$end)
		regions <- paste(regions, collapse = ';')
		return(regions)
	} else {
		return(NA)
		}
}

### READ AND SUMMARIZE SHATTERSEEK ################################################################
read_and_summarize_shatterseek <- function(summary) {
	dtf <- do.call(rbind, sapply(
		summary,
		function(x) {
			tmp <- read.delim(
				strsplit(x, '\t')[[1]][2],
				as.is = TRUE
				)
			tmp <- categorize_chromothripsis_confidence(tmp)
			# split into high and low 
			high <- tmp[which(tmp$confidence == 'High'),]
			low <- tmp[which(tmp$confidence == 'Low'),]
			# summarize regions
			high_regions <- summarize_regions(high)
			low_regions <- summarize_regions(low)
			# return dataframe of number of low and high confident chromothripsis and chromosomes 
			out <- data.frame(
				sample = gsub('.chromothripsis.summary.tsv', '', strsplit(x, '/')[[1]][11]),
				high_confidence = nrow(high),
				low_confidence = nrow(low),
				high_region = high_regions, 
				low_region = low_regions
				)
			return(out)
			},
		simplify = FALSE
		))
	return(dtf)
}
### MAIN ##########################################################################################
# find all summary files
summary <- system(
		'isabl get-results --result-key summary -fi projects 7 -fi application.pk 119 --verbose',
		intern = TRUE
		)
# summarize and categorize chromothripsis 
chromothripsis <- read_and_summarize_shatterseek(summary)

# write to file 
write.table(
	chromothripsis,
	file = file.path(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data',
		'chromothripsis_project7.txt'
		),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### METASTATIC ##############################################################################
# find all summary files
summary <- system(
		'isabl get-results --result-key summary -fi projects 17 -fi application.pk 119 --verbose',
		intern = TRUE
		)
# summarize and categorize chromothripsis 
chromothripsis <- read_and_summarize_shatterseek(summary[grep('No result', summary, invert = TRUE)])

# write to file 
write.table(
	chromothripsis,
	file = file.path(
		'/oak/stanford/groups/ccurtis2/users/khoulaha/BreastLandscape/data',
		'chromothripsis_project17.txt'
		),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

