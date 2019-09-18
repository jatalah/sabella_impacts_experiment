# 01 load libraries ----------------
suppressWarnings(suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(vegan)
  library(ggord)
  library(ggpubr)
  library(ggvegan)
  library(ggrepel)
  library(viridis)
  library(MASS)
}))
source('H:/javiera/Stats/R folder/theme_javier.R')
theme_set(theme_javier())

# read controls and ambient data-------------------------------
controls <-
  read_excel('data/cleaned_data/macrofauna_data_edited.xlsx', sheet = 'R') %>%
  dplyr::filter(Sampling == 2 &
                  Treatment == "Ambient" |
                  Treatment == 'Control') %>%
  dplyr::filter(sample != 25 &
                  sample != 30)
controls[is.na(controls)] <- 0

present_taxa_con <- 
  controls %>% 
  gather(taxa, value, 7:ncol(.)) %>% 
  group_by(taxa) %>% 
  summarise_at(vars(value), funs(sum), na.rm =T) %>%
  arrange(desc(value)) %>% 
  dplyr::filter(value>0) %>% 
  dplyr::select(taxa) %>% 
  as.list()

cont_dat <- 
  controls %>% 
  dplyr::select(sample,
                Treatment, 
                Density,
                present_taxa_con$taxa) %>% #dom_taxa$taxa #7:100
  arrange(sample) 

## permanvoa controls vs ambient---------------------
permanova_con <-
  adonis(
    sqrt(sqrt(cont_dat[, -c(1:6)] )) ~ Treatment ,
    data = controls,
    method =  "bray",
    permutations = 999
  )
permanova_con

mds_cont <- metaMDS(sqrt(sqrt(cont_dat[, -c(1:6)])), distance = 'bray')

# extract taxa with absoulte scores >0.6 in either MDS axes-------------- 
driving_taxa <- names(which(abs(mds_cont$species[,1])>.6|abs(mds_cont$species[,2])>.6))
ggord(
  mds_cont,
  grp_in = cont_dat$Treatment,
  poly = F,
  alpha = 1,
  arrow = F,
  txt = 3,
  var_sub = driving_taxa,
  repel = T
)

# Diversity indices macrofaunal controls--------------
indices_con <-
  cont_dat %>%
  mutate(N = rowSums(cont_dat[,-c(1:6)]),
         S = specnumber(cont_dat[,-c(1:6)]))

summary(lm(log(N) ~ Treatment, data = indices_con))
summary(lm(S ~ Treatment, data = indices_con))

# 02 euk_controls---------------
# read data
euk_cont <- 
  left_join(
    read_csv('data/cleaned_data/euks_rare.csv'),
    read_excel('data/cleaned_data/molecular_treatments.xlsx'),
    by = 'sample'
  ) %>%
  dplyr::select(sample, Plot:Density, everything()) %>% 
  filter(Treatment == "Ambient"| Treatment == "Control")
         

# get data present on these treatments only
euk_present_taxa <- 
  euk_cont %>% 
  gather(taxa, value, 4:ncol(.)) %>% 
  group_by(taxa) %>% 
  summarise_at(vars(value), funs(sum), na.rm =T) %>%
  dplyr::filter(value>0) %>% 
  dplyr::select(taxa) %>% 
  as.list()

# euk data selection 
euk_cont_dat <- 
  euk_cont %>% 
  dplyr::select(sample,
                Treatment, 
                Density,
                euk_present_taxa$taxa)

# permanova euk controls-----
permanova_euk_cont <-
  adonis(
    log10(euk_cont_dat[, 4:ncol(euk_cont_dat)] + 1) ~ Treatment ,
    data = euk_cont_dat,
    method =  "bray",
    permutations = 1e4
  )
permanova_euk_cont


mds_cont_euk <- metaMDS(log10(euk_cont_dat[, -c(1:4)]+1), distance = 'bray')

# extract taxa with absoulte scores >0.6 in either MDS axes
ggord(
  mds_cont_euk,
  grp_in = euk_cont_dat$Treatment,
  poly = F,
  alpha = 1,
  arrow = F,
  txt = NULL,
  vec_ext= 0
)

dis_euk_control <- vegdist(log10(euk_cont_dat[,-c(1:4)] + 1))
beta_disp_euk <-
  betadisper(
    dis_euk_control, euk_cont_dat$Treatment
  )
anova(beta_disp_euk)
plot(beta_disp_euk)

# Diversity indices euk controls-----------------
indices_con_euk <-
  euk_cont_dat %>%
  mutate(N = rowSums(euk_cont_dat[,-c(1:4)]),
         S = specnumber(euk_cont_dat[,-c(1:4)]))

anova(lm(log(N) ~ Treatment, data = indices_con_euk))
anova(lm(S ~ Treatment, data = indices_con_euk))

boxplot(S ~ Treatment, data = indices_con_euk)

# 03 bacterial data ---------------------
# data prep
bact_cont <- 
  left_join(bac,
            treat,
            by = 'sample') %>%
  dplyr::select(sample, Treatment:Density, everything()) %>% 
  filter(Treatment == "Ambient"| Treatment == "Control")

# get bact data present on these treatments only
bact_present_taxa <- 
  bact_cont %>% 
  gather(taxa, value, 5:ncol(.)) %>% 
  group_by(taxa) %>% 
  summarise_at(vars(value), funs(sum), na.rm =T) %>%
  dplyr::filter(value>0) %>% 
  dplyr::select(taxa) %>% 
  as.list()

# euk data selection -------
bact_cont_dat <- 
  bact_cont %>% 
  dplyr::select(sample,
                Treatment, 
                Density,
                bact_present_taxa$taxa)

# permanova euk controls----------------
permanova_bact_cont <-
  adonis(
    log10(bact_cont_dat[, 4:ncol(bact_cont_dat)] + 1) ~ Treatment ,
    data = bact_cont_dat,
    method =  "bray",
    permutations = 1e4
  )
permanova_bact_cont

mds_cont_bact <- metaMDS(log10(bact_cont_dat[, -c(1:4)]+1), distance = 'bray')

# extract taxa with absoulte scores >0.6 in either MDS axes-------------- 
ggord(
  mds_cont_bact,
  grp_in = bact_cont_dat$Treatment,
  poly = F,
  alpha = 1,
  arrow = F,
  txt = NULL,
  vec_ext= 0
)


dis_bact_control <- vegdist(log10(bact_cont_dat[,-c(1:4)] + 1))
beta_disp_bact <-
  betadisper(
    dis_bact_control, bact_cont_dat$Treatment
  )
anova(beta_disp_bact)
plot(beta_disp_bact)

# t-test on S and N for controls vs ambient-----------------
indices_con_bact<-
  bact_cont_dat %>%
  mutate(N = rowSums(bact_cont_dat[,-c(1:6)]),
         S = specnumber(bact_cont_dat[,-c(1:6)]))

anova(lm(N ~ Treatment, data = indices_con_bact))
anova(lm(S ~ Treatment, data = indices_con_bact))

boxplot(N ~ Treatment, data = indices_con_bact)
boxplot(S ~ Treatment, data = indices_con_bact)
