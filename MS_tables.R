library(tidyverse)
library(vegan)

# 01 PERMANOVA tables ---------------
macrofauna_permanova_table <- 
  data.frame(permanova$aov.tab) %>% 
  rownames_to_column() %>% 
  rename(Terms = rowname,
         df = Df,
         SS = SumsOfSqs,
         MS = MeanSqs,
         Pseudo_F = F.Model,
         P = Pr..F.)

euk_permanova_table <- 
  data.frame(permanova_euk$aov.tab) %>% 
  rownames_to_column() %>% 
  rename(Terms = rowname,
         df = Df,
         SS = SumsOfSqs,
         MS = MeanSqs,
         Pseudo_F = F.Model,
         P = Pr..F.)


bact_permanova_table <- 
  data.frame(permanova_bact$aov.tab) %>% 
  rownames_to_column() %>% 
  rename(Terms = rowname,
         df = Df,
         SS = SumsOfSqs,
         MS = MeanSqs,
         Pseudo_F = F.Model,
         P = Pr..F.)


all_permanova_tables <-
  bind_rows(
    "Macrofauna" = macrofauna_permanova_table,
    "Eukaryotes" = euk_permanova_table,
    "Bacteria" = bact_permanova_table,
    .id = "Dataset"
  ) %>%
  write_csv(., 'tables/all_permanova_tables.csv', na = "")


# CAP tables-------------------------------
macrofauna_cap_table <- 
  data.frame(cap_anova_macrofauna) %>% 
  rownames_to_column() %>% 
  rename(Terms = rowname,
         df = Df,
         SS = SumOfSqs,
         P = Pr..F.)

euk_cap_table <- 
  data.frame(cap_anova_euk) %>% 
  rownames_to_column() %>% 
  rename(Terms = rowname,
         df = Df,
         SS = SumOfSqs,
         P = Pr..F.)

bact_cap_table <- 
  data.frame(cap_anova_bact) %>% 
  rownames_to_column() %>% 
  rename(Terms = rowname,
         df = Df,
         SS = SumOfSqs,
         P = Pr..F.)



all_cap_tables <-
  bind_rows(
    "Macrofauna" = macrofauna_cap_table,
    "Eukaryotes" = euk_cap_table,
    "Bacteria" = bact_cap_table,
    .id = "Dataset"
  ) %>%
  write_csv(., 'tables/all_cap_tables.csv', na = "")


# SIMPER tables -------------
SIMPER_table <-
  bind_rows(
    Macrcofauna = simper_macrofauna,
    Eukaryote = simper_euk,
    Bacteria = simper_bact,
    .id = 'Dataset'
  ) %>%
  write_csv(., 'tables/SIMPER_tables.csv', na = "")



## ANCOVA tables--------------
ANCOVA_table <- 
bind_rows(
  "Macrofauna" = ancova_table_macro,
  "Eukaryotes" = ancova_table_euk,
  "Bacteria" = ancova_table_bact,
  .id = 'Dataset'
) %>% 
  write_csv('tables/ANCOVA_tables.csv', na = '')
