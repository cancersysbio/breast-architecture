### CREATE FIGURE 1d #############################################################################
# creates figure 1d
# figure 1d provides the difference in distance relapse-free (DRF) survival probability between ER+ Typical and ER+ High-risk classes detected by the four different IC subtypes classifiers, at 1, 5, 10, 15, and 20 years
# and the delta in Cox proportional model Hazard Ratio (HR) between ER+ Typical and ER+ High-risk classes compared with the reference model (iC10 DNA+RNA), at 1, 5, 10, 15, and 20 years.

### PREAMBLE #####################################################################################
library(yaml)
library(ggplot2)
library(ggpubr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

config = yaml.load_file(file.path(main_repo_path, 'config.yml'))

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

### MAIN ##########################################################################################
basedir <- '/oak/stanford/groups/ccurtis2/users/'

## see create_figure1c.R
load(file=paste0(basedir,"lisem/bc_landscape/submission/data/Figure1d_fit.drs.RData"))
names(fit.drs)

seq_time <- c(1, 5, 10, 15, 20)
for(fit in fit.drs){
  print(fit)
  print(summary(fit, times = seq_time)$surv)
}

## Scatter plot difference of KM estimate at different time point
df2 <- data.frame(Time = rep(seq_time,4),
                  Model = rep(names(fit.drs), each = length(seq_time)), 
                  HR = as.numeric(sapply(fit.drs, function(fit) summary(fit, times = seq_time)$surv[(length(seq_time)+1):(length(seq_time)*2)] - summary(fit, times = seq_time)$surv[1:length(seq_time)])),
                  upper = as.numeric(sapply(fit.drs, function(fit) summary(fit, times = seq_time)$upper[(length(seq_time)+1):(length(seq_time)*2)] - summary(fit, times = seq_time)$upper[1:length(seq_time)])),
                  lower = as.numeric(sapply(fit.drs, function(fit) summary(fit, times = seq_time)$lower[(length(seq_time)+1):(length(seq_time)*2)] - summary(fit, times = seq_time)$lower[1:length(seq_time)])), stringsAsFactors = F) 

df2$Model <- factor(df2$Model, levels = names(fit.drs))
p <- ggplot(df2, aes(x=Time, y=HR, group=Model, color=Model)) + 
  geom_line(linetype = 'dashed') + xlab('Time (years)') + ylab('DRF Survival Difference') +
  geom_point(size = 4, alpha = 0.7) + scale_color_manual(name = "", values = c("ENiClust"="#825CA6","iC10 DNA+RNA"="black","iC10 DNA-only"="#317EC2","iC10 RNA-only"="red")) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                position=position_dodge(0.05)) + 
  ggtitle('ER+ Typical - ER+ High-risk') +
  theme_LM + theme(axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12),
                   axis.title.x = element_text(size = 14),
                   axis.title.y = element_text(size = 14),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 16))

df2_ref = df2 %>% 
  dplyr::filter(Model %in% c("iC10 DNA+RNA")) %>%
  dplyr::filter(Time > 1)

p2 <- df2 %>% 
  dplyr::filter(Model %in% names(fit.drs)) %>%
  dplyr::filter(Time > 1) %>%
  dplyr::group_by(Model) %>%
  summarize(Model = Model, Time = Time, HR = HR - df2_ref$HR, lower = lower - df2_ref$lower, upper = upper - df2_ref$upper) %>%
  ggplot(aes(x=Time, y=HR, group=Model, color=Model)) + 
  geom_line(linetype = 'dashed') + xlab('Time (years)') + ylab('Delta HR') +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_point(size = 4) + 
  scale_color_manual(name = "", values = c("ENiClust"="#825CA6","iC10 DNA+RNA"="black","iC10 DNA-only"="#317EC2","iC10 RNA-only"="red")) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2, position=position_dodge(0.05)) + 
  ggtitle('ER+ Typical - ER+ High-risk') +
  theme_LM + theme(axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12),
                   axis.title.x = element_text(size = 14),
                   axis.title.y = element_text(size = 14),
                   legend.text = element_text(size = 14),
                   title = element_text(size = 16))

### SAVE ##########################################################################################
## Figure----
svg(paste0(basedir,"lisem/bc_landscape/submission/figures/Figure1d.svg"), w = 8.98, h = 4.50)
(p + theme(legend.position = "none")) + p2
dev.off()

## SourceData----
sourcetable = df2
write.table(sourcetable, paste0(basedir,"lisem/bc_landscape/submission/data/Figure1d_sourcetable.txt"), row.names = F, col.names = T, quote = F, sep="\t")
