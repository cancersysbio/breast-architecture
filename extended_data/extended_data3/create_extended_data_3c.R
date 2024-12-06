### CREATE EXTENDED DATA FIGURE 3C #############################################################################
# Barplot of SV types in HRD-like tumors

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(patchwork)
library(yaml)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### FUNCTIONS #####################################################################################
theme_LM = theme_classic() + grids()

plot_sv_prop <- function(plot_df, label_group, label_hrd, sv_type, col_pal) {
  plot_df %>%
    ggplot(aes(x=!!sym(label_hrd), y=Proportion, fill=!!sym(sv_type))) +
    geom_bar(stat='identity', position='stack') +
    facet_wrap(as.formula(paste0("~", label_group)), scales = 'free_x') + 
    labs(x = "Group", y = "Proportion of Events") +
    scale_fill_manual(values = col_pal) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_LM +
    theme(axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 12),
          title = element_text(size = 20),
          strip.text.x = element_text(size = 18), 
          strip.text.y = element_text(size = 18),
          legend.position = "right",
          legend.title = element_blank())
}

fisher_test_single_category <- function(contingency_table, hrd_labels, category1, category2) {
  fisher_matrix <- as.matrix(contingency_table[,-c(1,2)])
  rownames(fisher_matrix) <- hrd_labels
  
  fisher_test <- fisher.test(fisher_matrix)
  
  return(list(
    category1 = category1, 
    category2 = category2, 
    p_value = fisher_test$p.value, 
    odds_ratio = fisher_test$estimate
  ))
}

sv_prop_test <- function(hrd_prop_df, hrd_total_df, label_group, label_hrd, group2test, hrd_labels2test, sv_type_col, sv_type2test) {
  contingency_table <- hrd_prop_df %>%
    # left_join(hrd_total_df, by = c(label_group, label_hrd)) %>%
    left_join(hrd_total_df) %>%
    filter(!!sym(label_group) == group2test) %>%
    filter(!!sym(label_hrd) %in% hrd_labels2test) %>%
    filter(!!sym(sv_type_col) == sv_type2test) %>%
    mutate(countYes = prop * total_events,
           countNo = (1 - prop) * total_events) %>%
    select(!!sym(label_group), !!(label_hrd), countYes, countNo)
  if (nrow(contingency_table) == 0) {
    return(NULL)
  }
  else {
    fisher_test <- fisher_test_single_category(contingency_table, hrd_labels2test, 'countYes', 'countNo')
    
    return(list(
      group2test = group2test, 
      sv2test = sv_type, 
      p_value = fisher_test$p_value, 
      or = fisher_test$odds_ratio
    ))
    
  }
  
  
}


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
  "Tyfona" = "#ce2f79ff",
  "Other" = "#D3D3D3"
)

colour_palette_custom <- c(
  "BFB" = "#f9a41bff",
  "Chromoplexy" = "#428c40ff",
  "Chromothripsis" = "#682c8bff",
  "CPXDM" = "#fbd605ff",
  # "DEL" = "#aa1a91ff",
  "DEL (<10Kb)" = "#F5E4F1",
  "DEL (>10Kb)" = "#E7BDDF",
  "DM" = "#4a86c5ff",
  # "DUP" = "#63aaafff",
  "DUP (<10Kb)" = "#ECF4F5",
  "DUP (10-100Kb)" = "#B1D5D8",
  "DUP (>100Kb)" = "#95C7CC",
  "INV" = "#000000ff",
  "INVDUP" = "#cf2528ff",
  "Pyrgo" = "#9aca3bff",
  "Rigma" = "#f47f21ff",
  "TIC" = "#463c8dff",
  "TRA" = "#55bc7cff",
  "Tyfona" = "#ce2f79ff",
  "Other" = "#D3D3D3"
)

### MAIN ##########################################################################################
sv_sig_total <- read.delim(file.path(main_repo_path, "data", "ExtendedData3c_sourcetable.txt"), sep = ",")

non_clustered_del_cols <- grep("del", names(sv_sig_total), value = TRUE)
for (col in non_clustered_del_cols) {
  sv_sig_total[[paste0("proportion_", col)]] <- sv_sig_total[[col]] / sv_sig_total$total_non_clustered_del
}

non_clustered_tds_cols <- grep("tds", names(sv_sig_total), value = TRUE)
for (col in non_clustered_tds_cols) {
  sv_sig_total[[paste0("proportion_", col)]] <- sv_sig_total[[col]] / sv_sig_total$total_non_clustered_tds
}

sv_sig_total_prop <- sv_sig_total %>%
  select(label_chord_only, label_group, starts_with('proportion_')) %>%
  select(-c('proportion_total_non_clustered_del', 'proportion_total_non_clustered_tds')) %>%
  pivot_longer(cols=starts_with('proportion_'),
               values_to = 'prop',
               names_to='length') %>%
  mutate(sv_type = case_when(grepl('del', length) ~ 'DEL',
                             grepl('tds', length) ~ 'DUP')) %>%
  mutate(length = str_split(length, '_') %>% sapply(tail, 1))


hrd_label_group_df_long <- read.delim(file.path(main_repo_path, "data", "ExtendedData3c_hrd_sourcetable.txt"), sep = ",")

hrd_sv_prop_df <- hrd_label_group_df_long %>%
  group_by(label_chord_only, label_group, sv_type_tidy) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

hrd_sv_plot_df <- hrd_sv_prop_df %>% 
  full_join(sv_sig_total_prop, by=c('label_chord_only', 'label_group', 'sv_type_tidy' = 'sv_type')) %>%
  mutate(Proportion = case_when(is.na(`prop.y`) ~ `prop.x`,
                                !is.na(`prop.y`) ~ `prop.x` * `prop.y`)) %>%
  dplyr::mutate(sv_type_tidy = case_when(!(is.na(length)) ~ paste0(sv_type_tidy, ' (', length, ')'),
                                         T ~ sv_type_tidy)) %>%
  dplyr::mutate(sv_type_tidy = factor(sv_type_tidy, levels = c('DEL (<10Kb)', 'DEL (>10Kb)', 
                                                               'DUP (<10Kb)', 'DUP (10-100Kb)', 'DUP (>100Kb)', 
                                                               'DM', 'INVDUP', 'Pyrgo', 'Other')))


# Stats test 
groups_to_test <- c('ER+ High', 'ER+ Typical', 'TNBC')
chord_labels <- list('ER+ High' = c('BRCA2-like', 'non-HRD-like'),
                     'ER+ Typical' = c('BRCA2-like', 'non-HRD-like'))
sv_types_to_test <- c('DM', 'INVDUP')

hrd_total_df <- hrd_label_group_df_long %>%
  group_by(label_chord_only, label_group, sv_type_tidy) %>%
  summarise(n = n()) %>%
  group_by(label_chord_only, label_group) %>%
  summarise(total_events = sum(n)) %>%
  select(label_chord_only, label_group, total_events)

stat_results <- data.frame(group2test = character(),
                           sv2test = character(),
                           p_value = numeric(),
                           or = numeric(),
                           stringsAsFactors = FALSE)

for (group in groups_to_test) {
  for (sv_type in sv_types_to_test) {
    if (group == 'TNBC') {
      # have to test the BRCA1/2 individually
      for (hrd in c('BRCA1-like', 'BRCA2-like')) {
        hrd_labels2test = c(hrd, 'non-HRD-like')
        result <- sv_prop_test(hrd_prop_df = hrd_sv_prop_df,
                               hrd_total_df = hrd_total_df,
                               label_group = 'label_group',
                               label_hrd = 'label_chord_only',
                               group2test = group,
                               hrd_labels2test = hrd_labels2test,
                               sv_type_col = 'sv_type_tidy',
                               sv_type2test = sv_type)
        if (!is.null(result)) {
          stat_results <- rbind(stat_results,
                                data.frame(group2test = paste0('TNBC (', hrd, ')'),
                                           sv2test = sv_type,
                                           p_value = result$p_value,
                                           or = result$or,
                                           stringsAsFactors = FALSE))
        } else {
          break
        }
      }
    } else {
      result <- sv_prop_test(hrd_prop_df = hrd_sv_prop_df,
                             hrd_total_df = hrd_total_df,
                             label_group = 'label_group',
                             label_hrd = 'label_chord_only',
                             group2test = group,
                             hrd_labels2test = chord_labels[[group]],
                             sv_type_col = 'sv_type_tidy',
                             sv_type2test = sv_type)
      stat_results <- rbind(stat_results, result)
    }
  }
}


sv_sig_total_df <- sv_sig_total %>%
  select(label_chord_only, label_group, total_non_clustered_del, total_non_clustered_tds) %>%
  pivot_longer(cols = starts_with('total'),
               names_to = 'sv_type',
               values_to = 'total_events') %>%
  mutate(sv_type=substr(sv_type, nchar(sv_type) - 2, nchar(sv_type)))

sv_sig_prop_df <- sv_sig_total_prop %>%
  mutate(sv_type=case_when(sv_type == 'DEL' ~ 'del',
                           sv_type == 'DUP' ~ 'tds')) %>%
  mutate(sv_type_tidy = paste0(sv_type, '_', length))


groups_to_test <- c('ER+ High', 'ER+ Typical', 'TNBC')
lengths2test <- c('del_<10Kb', 'del_>10Kb', 'tds_<10Kb', 'tds_>10Kb')
for (group in groups_to_test) {
  for (lengths in lengths2test) {
    if (group == 'TNBC') {
      # have to test the BRCA1/2 individually
      for (hrd in c('BRCA1-like', 'BRCA2-like')) {
        hrd_labels2test = c(hrd, 'non-HRD-like')
        result <- sv_prop_test(hrd_prop_df = sv_sig_prop_df,
                               hrd_total_df = sv_sig_total_df,
                               label_group = 'label_group',
                               label_hrd = 'label_chord_only',
                               group2test = group,
                               hrd_labels2test = hrd_labels2test,
                               sv_type_col = 'sv_type_tidy',
                               sv_type2test = lengths)
        if (!is.null(result)) {
          stat_results <- rbind(stat_results,
                                data.frame(group2test = paste0('TNBC (', hrd, ')'),
                                           sv2test = lengths,
                                           p_value = result$p_value,
                                           or = result$or,
                                           stringsAsFactors = FALSE))
        } else {
          break
        }
      }
    } else {
      result <- sv_prop_test(hrd_prop_df = sv_sig_prop_df,
                             hrd_total_df = sv_sig_total_df,
                             label_group = 'label_group',
                             label_hrd = 'label_chord_only',
                             group2test = group,
                             hrd_labels2test = chord_labels[[group]],
                             sv_type_col = 'sv_type_tidy',
                             sv_type2test = lengths)
      stat_results <- rbind(stat_results,
                            data.frame(group2test = group,
                                       sv2test = lengths,
                                       p_value = result$p_value,
                                       or = result$or,
                                       stringsAsFactors = FALSE))
    }
  }
}


pdf(file.path(main_repo_path, 'plots', "ExtendedData3c.pdf"), height=6, width=7.5)
plot_sv_prop(plot_df=hrd_sv_plot_df, 
             label_group='label_group', 
             label_hrd='label_chord_only', 
             sv_type='sv_type_tidy', 
             col_pal=colour_palette_custom)
dev.off()
stat_results %>%
  write.table(file.path(main_repo_path, "plots", "ExtendedData3c_stats.txt"), 
              quote = F, row.names = F, sep = "\t")