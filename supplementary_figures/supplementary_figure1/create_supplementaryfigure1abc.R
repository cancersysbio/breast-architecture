### CREATE SUPPLEMENTARY FIGURES 1ABC #############################################################
# compare mutational burden against PCAWG

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

# Set the main path for repo
main_repo_path <- "/Users/khoulaha/git/breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}

###################################################################################################
# read in megatable 
megatable <- read.delim(
	file.path(main_repo_path, 'data', 'primary_megatable.txt'),
	as.is = TRUE
	)
megatable$Individual.System.ID <- substr(megatable$Sample, 1, 14)

## CREATE SUPPLEMENTARY FIGURE 1A ################################################################# 
# read in number of snvs
snvs <- read.delim(
    file.path(main_repo_path, 'data', 'SupplementaryFigure1a_sourcetable.txt'),
    as.is = TRUE
    )

# create plot 
create.scatterplot(
    pcawg_snvs ~ our_snvs,
    data = snvs,
    filename = 'supplementary_figure1a.png',
    ylab.label = 'PCAWG SNV Burden',
    xlab.label = 'SNV Burden',
    yat = seq(0,80000,20000),
    xat = seq(0,80000,20000),
    xlimits = c(0,80000),
    ylimits = c(0,80000),
    legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = snvs$pcawg_snvs,
                         y = snvs$our_snvs,
                         label.items = c('spearman','spearman.p'),
                         alpha.background = 0,
                         key.cex = 1.2
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1)
                 )
             ),
    resolution = 300
    )

### CREATE SUPPLEMENTARY FIGURE 1B ################################################################
# read in number of sv
pcawg_sv <- read.delim(
    file.path(main_repo_path, 'data', 'SupplementaryFigure1b_sourcetable.txt'),
    as.is = TRUE
    )

create.scatterplot(
	pcawg_sv ~ our_sv,
	data = pcawg_sv,
	filename = 'supplementary_figure1b.png',
	ylab.label = 'PCAWG SV Burden',
	xlab.label = 'SV Burden',
	legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = pcawg_sv$our_sv,
                         y = pcawg_sv$pcawg_sv,
                         label.items = c('spearman','spearman.p'),
                         alpha.background = 0,
                         key.cex = 1.2
                         )
                     ),
                 x = 0.01,
                 y = 0.99,
                 corner = c(0,1)
                 )
             ),
	resolution = 300
	)


### CREATE SUPPLEMENTARY FIGURE 1C ################################################################
# read in number of snvs
pga_df <- read.delim(
    file.path(main_repo_path, 'data', 'SupplementaryFigure1c_sourcetable.txt'),
    as.is = TRUE
    )

# create plot of fraction genome altered
pga_df$pcawg_pga <- pga_df$pcawg_pga/100
create.scatterplot(
    pcawg_pga ~ our_fraction_cna,
    data = pga_df,
    filename = 'supplementary_figure1c.png',
    ylab.label = 'PCAWG FGA',
    xlab.label = 'FGA',
    legend = list(
             inside = list(
                 fun = draw.key,
                 args = list(
                     key = get.corr.key(
                         x = pga_df$pcawg_pga,
                         y = pga_df$our_fraction_cna,
                         label.items = c('spearman','spearman.p'),
                         alpha.background = 0,
                         key.cex = 1.2
                         )
                     ),
                 x = 0.99,
                 y = 0.01,
                 corner = c(1,0)
                 )
             ),
    resolution = 300
    )

