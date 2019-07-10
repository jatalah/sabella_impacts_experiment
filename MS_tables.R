library(tidyverse)
library(vegan)

# 01 PERMANOVA tables ---------------
all_permanova_tables <-
  bind_rows(
    "Macrofauna" = macrofauna_permanova_table,
    "Eukaryotes" = euk_permanova_table,
    "Bacteria" = bact_permanova_table,
    .id = "Dataset"
  ) %>%
  filter(Terms != 'Total') %>% 
  mutate_at(vars(SS:R2), ~round(.,digits = 2)) %>% 
  mutate(P = signif(P,2)) %>% 
  write_csv(., 'tables/all_permanova_tables.csv', na = "")


# CAP tables-------------------------------
all_cap_tables <-
  bind_rows(
    "Macrofauna" = macrofauna_cap_table,
    "Eukaryotes" = euk_cap_table,
    "Bacteria" = bact_cap_table,
    .id = "Dataset"
  ) %>%
  mutate_at(vars(SS:F), ~round(.,digits = 2)) %>% 
  mutate(P = signif(P,2)) %>% 
  write_csv(., 'tables/all_cap_tables.csv', na = "")


# SIMPER tables -------------
SIMPER_table <-
  bind_rows(
    Macrofauna = simper_macrofauna,
    Eukaryote = simper_euk,
    Bacteria = simper_bact,
    .id = 'Dataset'
  ) %>%
  rename(P = p) %>% 
  mutate_at(vars(average:cusum), ~ round(., digits = 3)) %>%
  mutate(P = round(P, 3)) %>%
  write_csv(., 'tables/SIMPER_tables.csv', na = "")

## ANCOVA tables--------------
ANCOVA_table <- 
bind_rows(
  "Macrofauna" = ancova_table_macro,
  "Eukaryotes" = ancova_table_euk,
  "Bacteria" = ancova_table_bact,
  .id = 'Dataset'
) %>% 
  filter(term != "(Intercept)") %>% 
  rename(SS = sumsq,
         F = statistic,
         P = p.value) %>% 
  mutate_at(vars(SS:F), ~ round(., digits = 2)) %>%
  mutate(P = round(P, 3),
         term = fct_recode(term, `Final density` = "final_density")) %>%
  write_csv('tables/ANCOVA_tables.csv', na = '')
