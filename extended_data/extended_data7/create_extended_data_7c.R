### CREATE EXTENDED DATA 7C #######################################################################
# create multi-analysis of ER association with eccDNA

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(metafor)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

date <- Sys.Date()
### ICGC ############################################################################################
# read in nikzainal mrna 
dtf_nk <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7c_sourcetable_icgc.txt'),
	as.is = TRUE
	)

# test association for each subgroup 
nk_res <- list()
for (j in c('AA_ecdna','AA_complex')) {
	fit <- glm(as.formula(paste(j, '~ ESR1 + purity + ENiClust')), data = dtf_nk[which(dtf_nk$group == 'ER+ High'),], family = 'binomial')
	ci <- confint(fit)
	nk_res[[j]] <- data.frame(
		subtype = 'ER+ High',
		cohort = 'NikZainal',
		type = j,
		coef = coef(fit)[['ESR1']],
		p = summary(fit)$coefficients['ESR1',4],
		se = summary(fit)$coefficients['ESR1',2],
		l95 = ci['ESR1',1],
		u95 = ci['ESR1',2]
		)
}
nk_res <- do.call(rbind, nk_res)


### TCGA ############################################################################################
# read in rna 
dtf_tcga <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7c_sourcetable_tcga.txt'),
	as.is = TRUE
	)

# test association for each subgroup 
tcga_res <- list()
for (i in c('ER+ High','HER2+')) {
	for (j in c('AA_ecdna','AA_complex')) {
		if (i == 'ER+ High') {
			fit <- glm(as.formula(paste(j, '~ ESR1 + purity + ENiClust')), data = dtf_tcga[dtf_tcga$group == i,], family = 'binomial')
		} else if (i == 'HER2+') {
			fit <- glm(as.formula(paste(j, '~ ESR1 + purity')), data = dtf_tcga[which(dtf_tcga$group == i & dtf_tcga$ER == 1),], family = 'binomial')
			}
		ci <- confint(fit)
		tcga_res[[paste0(i, j)]] <- data.frame(
			subtype = i,
			cohort = 'TCGA',
			type = j,
			coef = coef(fit)[['ESR1']],
			p = summary(fit)$coefficients['ESR1',4],
			se = summary(fit)$coefficients['ESR1',2],
			l95 = ci['ESR1',1],
			u95 = ci['ESR1',2]
			)
	}
}
tcga_res <- do.call(rbind, tcga_res)

### HARTWIG #######################################################################################
# read in hartwig rna data
dtf_hartwig <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData7c_sourcetable_hartwig.txt'),
	as.is = TRUE
	)

# test association for each subgroup 
hartwig_res <- list()
for (i in c('ER+ High','HER2+')) {
	for (j in c('AA_ecdna','AA_complex')) {
		if (i == 'HER2+') {
			# too few samples to subset down to only HER2+/ER+ so including as a term instead
			fit <- glm(as.formula(paste(j, '~ ESR1 + purity + ER')), data = dtf_hartwig[which(dtf_hartwig$group == i),], family = 'binomial')
		} else if (i == 'ER+ High') {
			fit <- glm(as.formula(paste(j, '~ ESR1 + purity')), data = dtf_hartwig[dtf_hartwig$group == i,], family = 'binomial')
		} 
		ci <- confint(fit)
		hartwig_res[[paste0(i, j)]] <- data.frame(
			subtype = i,
			cohort = 'Hartwig',
			type = j,
			coef = coef(fit)[['ESR1']],
			p = summary(fit)$coefficients['ESR1',4],
			se = summary(fit)$coefficients['ESR1',2],
			l95 = ci['ESR1',1],
			u95 = ci['ESR1',2]
			)
	}
}
hartwig_res <- do.call(rbind, hartwig_res)

### CALCULATE META-ANALYSIS #######################################################################
plot_data <- rbind(nk_res, tcga_res, hartwig_res)

meta_plot_data <- list()
for (subtype in c('ER+ High','HER2+')) {
	for (type in c('AA_ecdna','AA_complex')) {
	meta_data <- escalc(measure="OR", yi = coef, sei = se, 
                        data = plot_data[plot_data$subtype == subtype & plot_data$type == type,])
        res <- rma(yi, vi, data=meta_data)
        meta_plot_data[[paste0(subtype, type)]] <- data.frame(
                        subtype = subtype,
                        type = type,
                        coef = res$beta[[1]],
                        l95 = res$ci.lb,
                        u95 = res$ci.ub,
                        p = res$pval
                        )
	}
}
meta_plot_data <- do.call(rbind, meta_plot_data)
meta_plot_data$group <- paste(meta_plot_data$subtype, meta_plot_data$type)
meta_plot_data <- meta_plot_data[order(meta_plot_data$subtype, meta_plot_data$type),]
meta_plot_data$index <- 1:nrow(meta_plot_data)

meta_plot_data_cyclic <- meta_plot_data[meta_plot_data$type == 'AA_ecdna',]
meta_plot_data_cyclic$index <- 1:nrow(meta_plot_data_cyclic)
cyclicplot <- create.scatterplot(
        index ~ coef,
        data = meta_plot_data_cyclic,
        horizontal = TRUE,
        xlimits = c(-8,8),
        xaxis.lab = rep('', 10),
        #xat = c(log(c(0.1,0.5)),0, log(c(2,8))),
        #xaxis.lab = c('0.1','0.5', '1.0', '2.0', '8.0'),
        #xlab.label = 'Odds Ratio',
        #ylab.label = '',
        #yaxis.lab = gsub('Her2','HER2+', unique(plot_data_subtype$subtype)),
        yat = 1:nrow(meta_plot_data_cyclic),
        yaxis.lab = c('ER+ High','HER2+/ER+'),
        ylimits = c(0.5,2.5),
        #yaxis.lab = plot_data_cyclic$cohort,
        #ylimits = c(0.5, nrow(plot_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        xlab.label= '',
        ylab.label = '',
        add.text = TRUE,
        text.labels = 'Cyclic',
        text.y = 0.75,
        text.x = 5.25,
        #text.labels = c('ER+ High','HER2+/ER+','Cyclic'),
        # text.x = c(-5.25,-5.25,5.25),
        # text.y = c(1.25,4,1.25),
        #filename = paste0(date, '_ecdna_esr1_scatterplot.png'),
        # legend = list(
        #         right = list(fun = cov.grob),
        #         inside = list(fun = cov.legend.grob, corner = c(1,0), x = 0.99, y = 0.01)
        #         ),
        add.rectangle = TRUE,
        xleft.rectangle = -8,
        ybottom.rectangle = c(0),
        xright.rectangle = 8,
        ytop.rectangle = c(1.5),
        col.rectangle = 'mediumpurple2',
        alpha.rectangle = 0.25,
        # width = 8,
        # height = 5,
        x.error.right = meta_plot_data_cyclic$u95-meta_plot_data_cyclic$coef,
        x.error.left = meta_plot_data_cyclic$coef-meta_plot_data_cyclic$l95,
        top.padding = 2,
        resolution = 300
        )
meta_plot_data_complex <- meta_plot_data[meta_plot_data$type == 'AA_complex',]
meta_plot_data_complex$index <- 1:nrow(meta_plot_data_complex)
complexplot <- create.scatterplot(
        index ~ coef,
        data = meta_plot_data_complex,
        horizontal = TRUE,
        xlimits = c(-8,8),
        xat = log(c(0.1,1,8)),
        xaxis.lab = c('0.1','1.0', '8.0'),
        #xlab.label = 'Odds Ratio',
        #ylab.label = '',
        #yaxis.lab = gsub('Her2','HER2+', unique(plot_data_subtype$subtype)),
        yat = 1:nrow(meta_plot_data_complex),
        yaxis.lab = c('ER+ High','HER2+/ER+'),
        ylimits = c(0.5,2.5),
        #ylimits = c(0.5, nrow(plot_data_subtype)+0.75),
        main.cex = 2,
        abline.v = 0,
        key = NULL,
        add.text = TRUE,
        xlab.label = 'Odds Ratio',
        ylab.label = '',
        text.labels = 'Non-cyclic',
        text.x = 5.25,
        text.y = 0.75,
        #text.labels = c('ER+ High','HER2+/ER+','Non-Cyclic'),
        #text.x = c(-5.25,-5.25,5.25),
        #text.y = c(1.25,4,1.25),
        #filename = paste0(date, '_ecdna_esr1_scatterplot.png'),
        # legend = list(
        #         right = list(fun = cov.grob),
        #         inside = list(fun = cov.legend.grob, corner = c(1,0), x = 0.99, y = 0.01)
        #         ),
        add.rectangle = TRUE,
        xleft.rectangle = -8,
        ybottom.rectangle = c(0),
        xright.rectangle = 8,
        ytop.rectangle = c(1.5),
        col.rectangle = 'mediumpurple2',
        alpha.rectangle = 0.25,
        # width = 8,
        # height = 5,
        x.error.right = meta_plot_data_complex$u95-meta_plot_data_complex$coef,
        x.error.left = meta_plot_data_complex$coef-meta_plot_data_complex$l95,
        top.padding = 2,
        resolution = 300
        )

create.multipanelplot(
	list(cyclicplot,complexplot),
	filename = 'extended_data7c.pdf',
	plot.objects.heights = c(0.45, 0.55),
	height = 5.5,
	main = 'ESR1 mRNA',
	main.cex = 2,
	resolution = 300
	)