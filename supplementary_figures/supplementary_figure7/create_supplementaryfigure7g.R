# CREATE SUPPLEMENTARY FIGURE 6G ----
# Forest plot neoantigen and GIE

# PREAMBLE ----
library(rstudioapi)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(yaml)
library(Publish)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# FUNCTIONS ----
theme_LM = theme_classic() + grids()
theme_LM = theme_LM + theme(axis.text.x = element_text(size = 12),
                            axis.text.y = element_text(size = 12),
                            axis.title.x = element_text(size = 14),
                            axis.title.y = element_text(size = 14),
                            legend.text = element_text(size = 14),
                            title = element_text(size = 16),
                            strip.text.x = element_text(size = 14), 
                            strip.text.y = element_text(size = 14))

group_colors <- c("ER+ High" = "#fa954eff", 
                  "ER+ Typical"="#a1c1d4ff", 
                  "HER2+" = "#bd854dff", 
                  "IC10" = "#904dbdff", 
                  "IC4ER-" = "#c2b7dbff")

ic_colours <- c("IC1"="#ff5500ff", "IC2"="#00ee77ff", "IC6"="#fffe41ff",
                "IC9"="#ee82edff", "IC5"="#8b0100ff", "IC3"="#cd32caff", "IC3/IC7"="#cd32caff",
                "IC8"="#feaa00ff", "IC10"="#7c26ccff", "IC4ER+"="#00c4cdff", "IC4ER-"="#c2b7dbff")

plot_forest <- function(df) {
  df %>%
    ggplot(data=., 
           aes(x=group, y=slope, ymin=cl, ymax=cu, color=model)) +
    geom_pointrange(position=position_dodge(width = 0.5)) + 
    geom_hline(yintercept=0, lty=2) + 
    coord_flip() + 
    labs(x="", y="", color="Model") +
    theme_bw()  +
    scale_color_manual(values = c("grey", "black")) +
    facet_wrap(~Stage, scales = "free") +
    theme_LM +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 14),
          title = element_text(size = 16))
}

get_model_data <- function(data, snv_col, neoantigen_col, n=4) {
  data_model <- data %>% 
    filter(group !="IC4ER-", !!as.name(snv_col) != 0) %>%
    mutate(log_snv_col = log(!!as.name(snv_col)))
  
  # Filter for q1 - n*(q3 -q1), q3 + n(q3 - q1) by group
  data_model <- data_model %>% group_by(group) %>%
    mutate(uplim_na = quantile(!!as.name(neoantigen_col))[[4]] + n*IQR(!!as.name(neoantigen_col)),
           lowlim_na = quantile(!!as.name(neoantigen_col))[[2]] - n*IQR(!!as.name(neoantigen_col)),
           uplim_snv = mean(log_snv_col) + n*IQR(log_snv_col),
           lowlim_snv = mean(log_snv_col) - n*IQR(log_snv_col)) %>%
    filter(!!as.name(neoantigen_col) > lowlim_na, !!as.name(neoantigen_col) < uplim_na,
           log_snv_col > lowlim_snv, log_snv_col < uplim_snv)
  
  models <- lapply(groups, function(x) {
    return(lm(formula(paste0(neoantigen_col, "~ log_snv_col*Altered")),
              data = data_model %>% filter(group == x)))
  })
  
  
  model_summary <- data.frame(
    group = groups, 
    p = sapply(models, function(x) summary(x)$coefficients[4, 4]), 
    coef = sapply(models, function(x) summary(x)$coefficients[4, 1])
  ) %>% mutate(padj = p.adjust(p, method = "fdr"))
  
  return(data_model %>% left_join(model_summary))
}
get_model_data_robust <- function(data, snv_col, neoantigen_col) {
  data_model <- data %>% 
    # filter(group !="IC4ER-", !!as.name(snv_col) != 0) %>%
    filter(!!as.name(snv_col) != 0) %>%
    mutate(log_snv_col = log(!!as.name(snv_col)))
  
  models <- lapply(groups, function(x) {
    return(MASS::rlm(formula(paste0(neoantigen_col, "~ log_snv_col*Altered")),
                     data = data_model %>% filter(group == x)))
  })
  
  res_models <- lapply(models, function(x) {(add_residual(x))})
  
  model_summary <- data.frame(
    group = groups, 
    p = sapply(models, function(x) sfsmisc::f.robftest(x, var = "log_snv_col:Altered1")$p.value), 
    coef = sapply(models, function(x) summary(x)$coefficients[4, 1]), # logSNVs:Altered1
    intercept_noGIE = sapply(models, function(x) summary(x)$coefficients[1, 1]), # Intercept
    slope_noGIE = sapply(models, function(x) summary(x)$coefficients[2, 1]), # logSNVs
    intercept_yesGIE = sapply(models, function(x) {summary(x)$coefficients[1, 1] + summary(x)$coefficients[3, 1]}),  # Intercept  + Altered1
    slope_yesGIE = sapply(models, function(x) {summary(x)$coefficients[2, 1] + summary(x)$coefficients[4, 1]}), # log(SNVs) + logSVS:Altered
    ci_lower = sapply(res_models, function(x) {gsub('\\[', '', gsub('\\]', '', str_split(Publish::publish(x)[[1]]$CI.95[3], ";")[[1]][1]))}),
    ci_upper = sapply(res_models, function(x) {gsub('\\[', '', gsub('\\]', '', str_split(Publish::publish(x)[[1]]$CI.95[3], ";")[[1]][2]))})
  ) %>% mutate(padj = p.adjust(p, method = "fdr"))
  
  return(data_model %>% left_join(model_summary))
}
get_model_data_robust_nointeraction <- function(data, snv_col, neoantigen_col) {
  data_model <- data %>% 
    filter(group !="IC4ER-", !!as.name(snv_col) != 0) %>%
    # filter(!!as.name(snv_col) != 0) %>%
    mutate(log_snv_col = log(!!as.name(snv_col)))
  
  models <- lapply(groups, function(x) {
    return(MASS::rlm(formula(paste0(neoantigen_col, "~ log_snv_col + Altered")),
                     data = data_model %>% filter(group == x)))
  })
  
  res_models <- lapply(models, function(x) {(add_residual(x))})
  
  model_summary <- data.frame(
    group = groups, 
    p = sapply(models, function(x) sfsmisc::f.robftest(x, var = "Altered1")$p.value), 
    slope_yesGIE = sapply(models, function(x) summary(x)$coefficients[3, 1]), # Altered1
    ci_lower = sapply(res_models, function(x) {gsub('\\[', '', gsub('\\]', '', str_split(Publish::publish(x)[[1]]$CI.95[4], ";")[[1]][1]))}),
    ci_upper = sapply(res_models, function(x) {gsub('\\[', '', gsub('\\]', '', str_split(Publish::publish(x)[[1]]$CI.95[4], ";")[[1]][2]))})
  ) %>% mutate(padj = p.adjust(p, method = "fdr"))
  
  return(data_model %>% left_join(model_summary))
}
get_model_data_robust_nosnv <- function(data, snv_col, neoantigen_col) {
  # Same samples--filter similarly
  data_model <- data %>% 
    filter(group !="IC4ER-", !!as.name(snv_col) != 0) %>%
    # filter(!!as.name(snv_col) != 0) %>%
    mutate(log_snv_col = log(!!as.name(snv_col)))
  
  models <- lapply(groups, function(x) {
    return(MASS::rlm(formula(paste0(neoantigen_col, "~ Altered")),
                     data = data_model %>% filter(group == x)))
  })
  
  res_models <- lapply(models, function(x) {(add_residual(x))})
  
  model_summary <- data.frame(
    group = groups, 
    p = sapply(models, function(x) sfsmisc::f.robftest(x, var = "Altered1")$p.value), 
    slope_yesGIE = sapply(models, function(x) {summary(x)$coefficients[2, 1]}),
    ci_lower = sapply(res_models, function(x) {gsub('\\[', '', gsub('\\]', '', str_split(Publish::publish(x)[[1]]$CI.95[3], ";")[[1]][1]))}),
    ci_upper = sapply(res_models, function(x) {gsub('\\[', '', gsub('\\]', '', str_split(Publish::publish(x)[[1]]$CI.95[3], ";")[[1]][2]))})
  ) %>% mutate(padj = p.adjust(p, method = "fdr"))
  
  return(data_model %>% left_join(model_summary))
}

# MAIN----
## Primary ----
primary_data <- fread(file.path(main_repo_path, "data", "SupplementaryFigure6g_sourcetable.txt")) %>%
  filter(Stage == "Primary")
col <- brewer.pal(4, "Paired")[3:4]
names(col) <- c("No", "Yes")

add_residual <- function(object, ...){
  object$df.residual <- MASS:::summary.rlm(object)$df[2]
  return(object)
}

groups <- c("ER+ High", "ER+ Typical", "HER2+", "IC10")
primary_clonal <- get_model_data_robust(primary_data, "Clonal_SNVs", "Clonal_Neoantigens")
primary_subclonal <- get_model_data_robust(primary_data, "Subclonal_SNVs", "Subclonal_Neoantigens")
primary_all <- get_model_data_robust(primary_data, "SNVs", "Neoantigens")


## Metastatic ----
metastatic_data <- fread(file.path(main_repo_path, "data", "SupplementaryFigure6g_sourcetable.txt")) %>%
  filter(Stage == "Metastatic")

metastatic_clonal <- get_model_data_robust(metastatic_data, "Clonal_SNVs", "Clonal_Neoantigens")
metastatic_subclonal <- get_model_data_robust(metastatic_data, "Subclonal_SNVs", "Subclonal_Neoantigens")
metastatic_all <- get_model_data_robust(metastatic_data, "SNVs", "Neoantigens")

## All ----
all_data <- rbind(primary_all %>% mutate(Stage ="Primary"), 
                  metastatic_all %>% mutate(Stage = "Metastatic")) %>% 
  mutate(Stage = factor(Stage, levels = c("Primary", "Metastatic"))) %>% select(-padj)

all_data <- all_data %>% 
  left_join(all_data %>% select(ER_Status, p,coef, Stage) %>% 
              distinct() %>% 
              mutate(padj = p.adjust(p)))

# plot
plots <- list()
quantiles = c(0, 0.2, 0.4, 0.6, 0.8)

primary_all_gie <- get_model_data_robust_nosnv(primary_data, "SNVs", "Neoantigens")
metastatic_all_gie <- get_model_data_robust_nosnv(metastatic_data, "SNVs", "Neoantigens")
all_data_gie <- rbind(primary_all_gie %>% mutate(Stage ="Primary"), 
                      metastatic_all_gie %>% mutate(Stage = "Metastatic")) %>% 
  mutate(Stage = factor(Stage, levels = c("Primary", "Metastatic"))) %>% select(-padj)
all_data_gie <- all_data_gie %>% left_join(all_data_gie %>% select(group, p, Stage) %>% distinct() %>% mutate(padj = p.adjust(p)))

models_all_gie <- all_data_gie %>% 
  select(group, Stage, slope_yesGIE, ci_lower, ci_upper) %>%
  dplyr::rename("cl_yesGIE"="ci_lower", "cu_yesGIE"="ci_upper") %>%
  mutate(cl_yesGIE = as.double(cl_yesGIE),
         cu_yesGIE = as.double(cu_yesGIE)) %>%
  distinct() %>%
  pivot_longer(!c(group, Stage)) %>%
  separate(name, into = c('intercept', 'GIE'), sep="_") %>%
  mutate(GIE=case_when(grepl("no", GIE) ~ "No",
                       T ~ "Yes")) %>%
  pivot_wider(names_from=intercept, values_from = value)

primary_all_gie_nointer <- get_model_data_robust_nointeraction(primary_data, "SNVs", "Neoantigens")
metastatic_all_gie_nointer <- get_model_data_robust_nointeraction(metastatic_data, "SNVs", "Neoantigens")

all_data_gie_nointer <- rbind(primary_all_gie_nointer %>% mutate(Stage ="Primary"), 
                              metastatic_all_gie_nointer %>% mutate(Stage = "Metastatic")) %>% 
  mutate(Stage = factor(Stage, levels = c("Primary", "Metastatic"))) %>% select(-padj)
all_data_gie_nointer <- all_data_gie_nointer %>% left_join(all_data_gie_nointer %>% select(group, p, Stage) %>% distinct() %>% mutate(padj = p.adjust(p)))

models_all_nointer <- all_data_gie_nointer  %>%
  select(group, Stage, slope_yesGIE, ci_lower, ci_upper) %>%
  dplyr::rename("cl_yesGIE"="ci_lower", "cu_yesGIE"="ci_upper") %>%
  mutate(cl_yesGIE = as.double(cl_yesGIE),
         cu_yesGIE = as.double(cu_yesGIE)) %>%
  distinct() %>%
  pivot_longer(!c(group, Stage)) %>%
  separate(name, into = c('intercept', 'GIE'), sep="_") %>%
  mutate(GIE=case_when(grepl("no", GIE) ~ "No",
                       T ~ "Yes")) %>%
  pivot_wider(names_from=intercept, values_from = value)

forest_plot_df <- rbind(models_all_gie %>% mutate(model="GIE"), models_all_nointer %>% filter(GIE=="Yes") %>% mutate(model="GIE + log(SNV)"))

# SAVE ----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure6g.pdf"), height = 3, width = 10)
forest_plot_df %>% 
  mutate(group=factor(group, levels=c("IC10", "HER2+", "ER+ Typical", "ER+ High"))) %>%
  plot_forest(.)
dev.off()
