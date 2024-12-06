# CREATE SUPPLEMENTARY FIGURE 2g ----
# Create boxplot of RS4 and RS6 activity by alteration in SV type

# PREAMBLE ----
library(rstudioapi)
library(tidyr)
library(data.table)
library(dplyr)
library(yaml)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(broom)
library(scales)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set file directory as working directory
main_repo_path <- "../../../breast-architecture"
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

# FUNCTIONS ----
plot_sv_sig_sv_type <- function(plot_df, df_counts) {
  p <- plot_df %>%
    mutate(Signature = factor(Signature, levels = c("RS6", "RS4"))) %>%
    ggboxplot(data = .,
              x = 'Altered',
              y = 'Proportion',
              add = 'jitter',
              color = 'SVType',
              add.params = list(size = 1)) +
    labs(x="Altered", y = "Activity") +
    theme_LM +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          title = element_text(size = 16),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 14),
          legend.position = "top",
          legend.title = element_blank()) +
    scale_color_manual(values = colour_palette) +
    guides(color = guide_legend(nrow = 1)) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^floor(x)),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = expansion(mult = c(0.05, 0.1))
    )
  
  p <- p +
    geom_text(data = plot_df %>% distinct(Signature, SVType, .keep_all = TRUE),
              aes(x = 1.5, y = y, label = paste0("p=", p.val)),
              size = 3,
              hjust = 0.5)
  
  p <- facet(p, facet.by = c('Signature', 'SVType'),
             scales = 'free_y')
  
  p_type <- p +
    stat_summary(data = df_counts %>%
                   filter(Signature == 'RS6') %>%
                   top_n(1, Proportion),
                 aes(x = Altered,
                     label = paste0("n=", n),
                     y = 10^3+800),
                 geom = "text",
                 size = 3,
                 hjust = 0.5)
  
  return(p_type)
}

theme_LM = theme_classic() + grids()

colour_palette <- c(
  "BFB" = "#f9a41bff",
  "Chromoplexy" = "#428c40ff",
  "Chromothripsis" = "#682c8bff",
  "CPXDM" = "#fbd605ff",
  "DEL" = "#aa1a91ff",
  "DM" = "#4a86c5ff",
  "DUP" = "#63aaafff",
  "INV" = "#000000ff",
  "INVDUP" = "#cf2528ff",
  "Pyrgo" = "#9aca3bff",
  "Rigma" = "#f47f21ff",
  "TIC" = "#463c8dff",
  "TRA" = "#55bc7cff",
  "Tyfona" = "#ce2f79ff"
)


# MAIN ----
input_df <- fread(file.path(main_repo_path, "data", "SupplementaryFigure2fg_sourcetable.txt"))
jabba_types <- c(
  'BFB', 'Chromoplexy', 'Chromothripsis', 'CPXDM',
  'DEL', 'DM', 'DUP', 'INV', 'INVDUP', 'Pyrgo', 
  'Rigma', 'TIC', 'TRA', 'Tyfona'
)

med_res_df <- data.frame(jabba_type = character(),
                         b_0_modelA = numeric(),
                         p_0_modelA = numeric(), 
                         b_1_modelB = numeric(),
                         p_1_modelB = numeric(),
                         b_0_modelC = numeric(),
                         p_0_modelC = numeric(), 
                         b_1_modelC = numeric(),
                         p_1_modelC = numeric(),
                         stringsAsFactors = FALSE)


for (jabba_type in jabba_types) {
  formula_0 <- as.formula(paste(jabba_type, "~ RS6"))
  model.0 <- lm(formula_0, data = input_df)
  
  formula_2 <- as.formula(paste(jabba_type, "~ RS4"))
  model.2 <- lm(formula_2, data = input_df)
  
  formula_1 <- as.formula(paste(jabba_type, "~ RS6 + RS4"))
  model.1 <- lm(formula_1, data = input_df)
  
  # p-vals
  p_value_0 <- summary(model.0)$coefficients[, "Pr(>|t|)"][['RS6']]
  p_value_3 <- summary(model.1)$coefficients[, "Pr(>|t|)"][['RS4']]
  p_value_1 <- summary(model.1)$coefficients[, "Pr(>|t|)"][['RS6']]
  p_value_2 <- summary(model.1)$coefficients[, "Pr(>|t|)"][['RS4']]
  
  # b
  coef_0 <- coef(model.0)["RS6"]
  coef_3 <- coef(model.1)["RS4"]
  coef_1 <- coef(model.1)["RS6"]
  coef_2 <- coef(model.1)["RS4"]
  
  med_res_df <- rbind(med_res_df, 
                      data.frame(jabba_type = jabba_type, 
                                 b_0_modelA = coef_0,
                                 p_0_modelA = p_value_0, 
                                 b_1_modelB = coef_3,
                                 p_1_modelB = p_value_3,
                                 b_0_modelC = coef_1,
                                 p_0_modelC = p_value_1,
                                 b_1_modelC = coef_2,
                                 p_1_modelC = p_value_2))
  
  
}

# plot specific sv types 
svoi <- med_res_df %>%
  filter(p_0_modelA < 0.05 & p_0_modelC < 0.05 & p_1_modelC > 0.05 |
           p_1_modelB < 0.05 & p_1_modelC < 0.05 & p_0_modelC > 0.05) %>%
  pull(jabba_type)

plot_df <- input_df %>%
  select(all_of(svoi), RS4, RS6) %>%
  mutate_at(vars(BFB:Tyfona), ~ ifelse(is.na(.), NA, as.logical(.))) %>%
  pivot_longer(cols = starts_with('RS'),
               names_to = 'Signature',
               values_to = 'Proportion') %>%
  pivot_longer(cols = matches("BFB|Chromothripsis|CPXDM|DUP|INV|Tyfona"),
               names_to = "SVType",
               values_to = "Altered") %>%
  group_by(Signature, SVType) %>%
  filter(!is.na(Altered)) %>%
  mutate(Altered=ifelse(Altered == TRUE, 'Yes', 'No')) %>%
  mutate(Altered=factor(Altered, levels = c('Yes', 'No'))) %>%
  mutate(Signature=factor(Signature, levels = c('RS6', 'RS4'))) %>%
  ungroup()

df_counts <- plot_df %>%
  group_by(Signature, SVType, Altered) %>%
  count() %>%
  right_join(plot_df)

test_results <- med_res_df %>%
  select(jabba_type, p_0_modelC, p_1_modelC) %>%
  pivot_longer(cols = starts_with("p_"), 
               names_to = "model", 
               values_to = "p.value") %>%
  mutate(SVType = jabba_type,
         Signature = if_else(model == "p_0_modelC", "RS6", "RS4"),
         statistic = NA,  # Placeholder for statistic
         method = "Wilcoxon rank sum test with continuity correction",
         alternative = "two.sided",
         p.val = formatC(p.value, format = "e", digits = 2)) %>%
  mutate(Signature = factor(Signature, levels = c("RS6", "RS4")))

plot_df <- input_df %>%
  select(all_of(svoi), RS4, RS6) %>%
  mutate_at(vars(BFB:Tyfona), ~ ifelse(is.na(.), NA, as.logical(.))) %>%
  pivot_longer(cols = starts_with('RS'),
               names_to = 'Signature',
               values_to = 'Proportion') %>%
  pivot_longer(cols = matches("BFB|Chromothripsis|CPXDM|DUP|INV|Tyfona"),
               names_to = "SVType",
               values_to = "Altered") %>%
  mutate(Signature = factor(Signature, levels = c('RS6', 'RS4'))) %>%
  group_by(Signature, SVType) %>%
  filter(!is.na(Altered)) %>%
  mutate(Altered=ifelse(Altered == TRUE, 'Yes', 'No')) %>%
  mutate(Altered=factor(Altered, levels = c('Yes', 'No'))) %>%
  mutate(Signature=factor(Signature, levels = c('RS6', 'RS4'))) %>%
  ungroup()

plot_df <- plot_df %>%
  left_join(test_results, by = c("Signature", "SVType")) %>%
  mutate(y = ifelse(Signature == 'RS4', 10^4, 10^3))

# SAVE ----
med_res_df %>%
  filter(p_0_modelA < 0.05 & p_0_modelC < 0.05 & p_1_modelC > 0.05 |
           p_1_modelB < 0.05 & p_1_modelC < 0.05 & p_0_modelC > 0.05) %>%
  fwrite(file.path(main_repo_path, "plots", "SupplementaryFigure2g_stats.tsv"))

pdf(file.path(main_repo_path, "plots", "SupplementaryFigure2g.pdf"), width=6.5, height=5)
plot_sv_sig_sv_type(plot_df, df_counts)
dev.off()