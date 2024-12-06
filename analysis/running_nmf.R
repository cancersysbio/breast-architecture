### CREATE SUPPLEMENTARY FIGURE 4e #############################################################################
# creates supplementary figure 4e
# supplementary figure 4e provides the correlation between various mutational features and the distance to each archetype with archetypes defined by Pareto or NMF.

### PREAMBLE #####################################################################################
library(yaml)
library(NMF)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
basedir = "/oak/stanford/groups/ccurtis2/users/"
outdir = paste0(basedir,"lisem/bc_landscape/revision/figures/")

## PCA Primary-only
x = read.table(paste0(basedir,'khoulaha/BreastLandscape/data/primary_sv_cn_sig_proportions.txt'), sep = "\t", header = T)
rownames(x) = x$sample
x = x[,-which(colnames(x) == "sample")]

# Estimating the factorization rank
V.random <- randomize(x)
estim.r.random <- nmf(V.random, 2:6, nrun=200, seed=123456)
save(estim.r.random, file = paste0(basedir,"lisem/bc_landscape/revision/data/estim.r.random_prim.RData"))

estim.r <- nmf(x, 2:6, nrun = 200, seed=123456) 
save(estim.r, file = paste0(basedir,"lisem/bc_landscape/revision/data/estim.r_prim.RData"))

png(filename = file.path(outdir, 'nmf_estim_rank_prim.png'), res = 300, width = 6.62, height = 6.52, units = 'in')
plot(estim.r, estim.r.random) # to estimate the factorization rank (or number of clusters)
dev.off()

# Multiple runs with rank = 3
res.multirun <- nmf(x, 3, nrun = 200, seed=123456)
save(res.multirun, file = paste0(basedir,"lisem/bc_landscape/revision/data/nmf_res_prim.RData"))

## PCA Metastatic-only
x = read.table(paste0(basedir,'khoulaha/BreastLandscape/data/metastatic_sv_cn_sig_proportions.txt'), sep = "\t", header = T)
rownames(x) = x$sample
x = x[,-which(colnames(x) == "sample")]

# Estimating the factorization rank
V.random <- randomize(x)
estim.r.random <- nmf(V.random, 2:6, nrun=200, seed=123456)
save(estim.r.random, file = paste0(basedir,"lisem/bc_landscape/revision/data/estim.r.random_mets.RData"))

estim.r <- nmf(x, 2:6, nrun = 200, seed=123456) 
save(estim.r, file = paste0(basedir,"lisem/bc_landscape/revision/data/estim.r_mets.RData"))

png(filename = file.path(outdir, 'nmf_estim_rank_mets.png'), res = 300, width = 6.62, height = 6.52, units = 'in')
plot(estim.r, estim.r.random) # to estimate the factorization rank (or number of clusters)
dev.off()

# Multiple runs with rank = 3
res.multirun <- nmf(x, 3, nrun = 200, seed=123456)
save(res.multirun, file = paste0(basedir,"lisem/bc_landscape/revision/data/nmf_res_mets.RData"))

## PCA Primary-Metastatic
x = read.table(paste0(basedir,'khoulaha/BreastLandscape/data/primary_metastatic_sv_cn_sig_proportions.txt'), sep = "\t", header = T)
rownames(x) = x$sample
x = x[,-which(colnames(x) == "sample")]

# Estimating the factorization rank
V.random <- randomize(x)
estim.r.random <- nmf(V.random, 2:6, nrun=200, seed=123456)
save(estim.r.random, file = paste0(basedir,"lisem/bc_landscape/revision/data/estim.r.random_all.RData"))

estim.r <- nmf(x, 2:6, nrun = 200, seed=123456) 
save(estim.r, file = paste0(basedir,"lisem/bc_landscape/revision/data/estim.r_all.RData"))

png(filename = file.path(outdir, 'nmf_estim_rank_all.png'), res = 300, width = 6.62, height = 6.52, units = 'in')
plot(estim.r, estim.r.random) # to estimate the factorization rank (or number of clusters)
dev.off()

# Multiple runs with rank = 3
res.multirun <- nmf(x, 3, nrun = 200, seed=123456)
save(res.multirun, file = paste0(basedir,"lisem/bc_landscape/revision/data/nmf_res_all.RData"))

### SAVE ##########################################################################################
load(file = paste0(basedir,"lisem/bc_landscape/revision/data/nmf_res_prim.RData"))
nmf_rank_prim <- as.data.frame(res.multirun@fit@W)
colnames(nmf_rank_prim) <- paste0("LF",1:3)

load(file = paste0(basedir,"lisem/bc_landscape/revision/data/nmf_res_mets.RData"))
nmf_rank_met <- as.data.frame(res.multirun@fit@W)
colnames(nmf_rank_met) <- paste0("LF",1:3)

load(file = paste0(basedir,"lisem/bc_landscape/revision/data/nmf_res_all.RData"))
nmf_rank_all <- as.data.frame(res.multirun@fit@W)
colnames(nmf_rank_all) <- paste0("LF",1:3)

all(rownames(nmf_rank_prim) %in% tern_arc_prim$Sample)
all(rownames(nmf_rank_met) %in% tern_arc_met$Sample)
all(rownames(nmf_rank_all) %in% tern_arc_all$Sample)

all(rownames(nmf_rank_prim) == tern_arc_prim$Sample)
all(rownames(nmf_rank_met) == tern_arc_met$Sample)
all(rownames(nmf_rank_all) == tern_arc_all$Sample)

## From comparison Pareto vs. NMF----
colnames(nmf_rank_prim) = paste0("LF",c(3,1,2))
colnames(nmf_rank_met) = paste0("LF",c(3,1,2))
colnames(nmf_rank_all) = paste0("LF",c(1,3,2))

write.table(nmf_rank_prim, paste0(basedir,"lisem/bc_landscape/revision/data/nmf_rank_prim.txt"), row.names = TRUE,  quote=F)
write.table(nmf_rank_met, paste0(basedir,"lisem/bc_landscape/revision/data/nmf_rank_met.txt"), row.names = TRUE,  quote=F)
write.table(nmf_rank_all, paste0(basedir,"lisem/bc_landscape/revision/data/nmf_rank_all.txt"), row.names = TRUE,  quote=F)
