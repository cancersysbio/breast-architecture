# Plot configs and color palettes
library(RColorBrewer)

theme_LM <- theme_classic() + grids()

alteration_colors <- c("Homozygous deletion" = "blue3", "Amplification" = "firebrick2", 
                       "Damaging structural variant" = "darkorchid3", "Protein coding SNV" = "goldenrod2", 
                       "Complex amplification" = "darkorange2", "Non damaging structural variant" = "turquoise4")

immune_alt_colors <- c("Homozygous deletion" = "blue3", "Amplification" = "firebrick2", 
                       "Structural variant" = "darkorchid3", "Single nucleotide variant" = "goldenrod2", 
                       "HLA LOH" = "dodgerblue2", "Cyclic amplification" = "darkorange2", 
                       "Complex non-cyclic amplification" = "firebrick4", "Multiple" = "darkgreen")

ic_colors <- c("ic1"= "#FF5500", "ic2"="#00EE76", "ic3" = "#cd32caff", 
               "ic4"="#bbe0e2ff", "ic5"="#8B0000", "ic6"="#FFFF40", "ic7" = "blue", 
               "ic8"="#FFAA00", "ic9"="#EE82EE", "ic10"="#7D26CD", "IC3/IC7" = "#cd32caff", 
               "ic4ER-"="#c2b7dbff", "ic4ER+"= "#00c4cdff")
names(ic_colors) <- gsub("ic", "IC", names(ic_colors))

group_colors <- c("ER+ High" = "#fa954eff", "ER+ Typical"="#a1c1d4ff", "HER2+" = "#8b0100ff", "IC10" = "#904dbdff", "IC4ER-" = "#c2b7dbff")

tme_colors <- c("IE/F"="#5CBDB9", "IE"="#5036C4", "F"="#AA1946", "D"="#B48913")

immunophenotypes_colors <- rev(brewer.pal(n = 6, name = "RdBu"))
names(immunophenotypes_colors) <- 1:6

