### CREATE FIGURE 3C #############################################################################
# creates figure 2c
# plot example ecDNA copy number and supporting discordant reads

### PLOT GGRAPH ECDNA #############################################################################
library(plyranges)
library(fishHook)
library(gGnome)

date <- Sys.Date()
### MAIN ##########################################################################################
# read in jabba
jabbafile <- ''
if (jabbafile == '') {
	stop("Please specify JaBbA output file, genome graph with integer copy numbers ...")
}
jabba <- gG(jabba = jabbafile)
jabba_amp <- amp(jabba)

# read in coverage file 
covfile <- ''
if (covfile == '') {
	stop("Please specify the tumor/normal ratio coverage file. Same as the file that used to run JaBbA ...")
}
cov <- read_bigwig(
	covfile
	)
sg.cov.gt <- gTrack(cov, y.field = 'score', max.ranges = 1e4, lwd.border = 0.3, circles = TRUE, name = 'Read Depth')

# read in gencode annotation
gencodefile <- ''
if (gencodefile == '') {
	stop("Please specify gencode annotation file. We used Gencode v39lift37 ...")
}
gencode = track.gencode(
	gencode = gencodefile,
	stack.gap = 1e5, cex.label = 0.8, height = 20,
	grep = c('CCND1','FGFR1','ZNF703','NARS2','PAK1'))

regions_to_plot <- makeGRangesFromDataFrame(data.frame(
	chr = c(11,11,8),
	start = c(69010222, 73488024, 35000000),
	end = c(71517158, 80018321,43750000)
	))

pdf('figure3c.pdf')
plot(c(gencode, sg.cov.gt, jabba_amp$gtrack(y.field = 'cn')), regions_to_plot)
dev.off()
