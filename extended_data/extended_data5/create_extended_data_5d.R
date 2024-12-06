### CREATE EXTENDED DATA 5D #######################################################################
# create oncogene ecDNA barplot for DCIS

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(GenomicRanges)
library(biovizBase)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/germline-epitopes' and try again.")
}


date <- Sys.Date()
### MAIN ##########################################################################################
# read in plot data
oncogene_dist <- read.delim(
	file.path(main_repo_path, 'data', 'ExtendedData5d_sourcetable.txt'),
	as.is = TRUE
	)

# create barplot
create.barplot(
	Prop ~ subtype,
	groups = oncocode,
	stack = TRUE,
	ylimits = c(0,1),
	yat = seq(0,1,0.2),
	xaxis.cex = 1.1,
	#xaxis.rot = 45,
	main = 'DCIS',
	xaxis.lab = c('ER+\nHigh','ER+\nTypical','HER2+','IC10','IC4ER-'),
	main.cex = 2,
	 legend = list(
             right = list(
                 fun = draw.key,
                 args = list(
                     key = list(
                         points = list(
                             col = 'black',
                             pch = 22,
                             cex = 3,
                             fill = c('white', 'darkgrey','black')
                             ),
                         text = list(
                             lab = c('No\nOncogene','Alternative\nOncogene','IC-specific\nOncogene')
                             ),
                         padding.text = 5,
                         cex = 1
                         )
                     )
                 )
             ),
	data = oncogene_dist,
	add.text = TRUE,
	text.labels = oncogene_dist$Freq[c(3,2,1,6,5,4,9,8,7,15,13)],
	text.x = c(1,1,1,2,2,2,3,3,3,5,5),
	text.y = c(0.82, 0.62, 0.35, 0.96, 0.7, 0.35, 0.85,0.96,
		0.35, 0.62, 0.15),
	#text.col = 'black',
	text.col = c(rep('black',8),'white',rep('black',2)),
	#col = c('black','darkgrey','white'),
	col = c(
		'#fa954eff','#fccaa7','white',
		'#a1c1d4ff','#d0e0e9','white',
		'#8b0100ff','#c58087','white',
		'#7c26ccff','#7c26cc6c','white',
		'#c2b7dbff','#e6e2ff','white'
		),
	ylab.label = 'Proportion',
	xlab.label = 'Subtype',
	filename = 'extended_data5d.pdf',
	resolution = 300
	)