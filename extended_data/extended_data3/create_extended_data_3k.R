### CREATE EXTENDED DATA FIGURE 3K ################################################################
# plot archetype prognosis
# create extended data figure 3k

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(survival)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

date <- Sys.Date()
### MAIN ##########################################################################################
# calculate read in clinical annotations
st6 <- read.delim(file.path(main_repo_path, 'data','clinic_METABRIC.txt'), as.is = TRUE)
st6$ER <- st6$ER.Status
st6 <- st6[,c('METABRIC.ID','STAGE','SIZE','GRADE','TDR_20','DR_20','AGE','ER','LN','BS')]
colnames(st6) <- c('sample','stage','size','grade','TDR','DR','age','ER','LN','BS')

# check if archetypes have been generated yet
## TODO: update archetypes file list
archetypes_file <- file.path(main_repo_path, 'data','archetypes_file.txt')
if (!file.exists(archetypes_file)) {
	stop("Please run script XXX ...")
}
# read in archetypes 
arch <- read.delim(
	archetypes_file,
	as.is = TRUE
	)
arch <- merge(arch, st6, by = 'sample')

res <- list()
for (i in c('Arc1','Arc2','Arc3')) {
		fit <- coxph(as.formula(paste('Surv(TDR, DR) ~', i, '+ ER + age + grade + size + LN')), data = arch)
		ci <- confint(fit)
		res[[i]] <- data.frame(
			arch = i,
			#subtype = j,
			hr = coef(fit)[[i]],
			p = summary(fit)$coefficients[i,5],
			l95 = ci[i,1],
			u95 = ci[i,2]
			)
	#}
}
res <- do.call(rbind, res)

# create plot data
plot_data <- res
plot_data$hr <- as.numeric(as.character(plot_data$hr))
plot_data$l95 <- as.numeric(as.character(plot_data$l95))
plot_data$u95 <- as.numeric(as.character(plot_data$u95))
plot_data$index <- 1:nrow(plot_data)

# create plot
create.scatterplot(
	hr ~ index,
	horizontal = TRUE,
	data = plot_data,
	xlimits = c(0.5,3.5),
	xat = 1:3,
	xaxis.lab = c('TNBC\n-enriched','ER+ Typical\n-enriched','ER+ High/HER2+\n-enriched'),
	filename = 'extended_data3k.pdf',
	y.error.up = plot_data$u95-plot_data$hr,
    y.error.down = plot_data$hr-plot_data$l95,
    yat = c(-0.69,0,0.92),
    yaxis.lab = c(0.5,1,2.5),
    cex = 1.2,
    error.bar.lwd = 1.5,
    xaxis.rot = 90,
    left.padding = 10,
    abline.h = 0,
    abline.lty = 2,
    ylab.label = 'Hazard Ratio',
    xlab.label = 'Archetype',
    ylimits = c(-1.5,1.5),
    height = 8,
    width = 5,
	resolution = 300
	)

### SAVE ##########################################################################################
# save source data
write.table(
	plot_data,
	file = file.path(main_repo_path, 'data', 'ExtendedData3k_sourcetable.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
