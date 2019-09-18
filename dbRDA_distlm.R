# rda
library(vegan)
library(ggord)
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(caret)


# Macrofauna distlm-----
env <-  
  read_csv('data/env_vars.csv') %>% 
  arrange(Plot)

fauna_dat <- 
  read_csv('data/fauna_dat.csv') %>% 
  rename(Plot = sample) %>% 
  arrange(Plot)

# 3. Transformation and scaling-------------
# t_fauna <- sqrt(sqrt(fauna_dat[,-c(1:4)]))

t_fauna <- log(fauna_dat[,-c(1:4)]+1)
# st_para <-
#   env %>%
#   select(Gravel:Phaeo) %>%
#   data.frame %>%
#   preProcess(method = "YeoJohnson")
# 
# st_preds <-
#   predict(st_para, data.frame(env)) %>%
#   select(Gravel:Phaeo)
# 
# env1 <- 
#   env %>% 
#   dplyr::select(Gravel:Phaeo) 
# names(env1)

env2 <- 
  env1 %>% 
  dplyr::select(Gravel, Medium, Fine, Organics:Phaeo) 

# VIF -----
diag(solve(cor(env2))) %>% enframe()

env2 %>% 
  gather() %>% 
  ggplot(aes(value))+
  geom_histogram()+
  facet_wrap(~key,scales = 'free')

# macrofauna distlm-----
macro_distlm <- 
  adonis2(t_fauna ~ .,
          data = env2,
          method = 'bray',
          by = 'margin')

macro_distlm

adonis2(t_fauna ~ .,
        data = env2,
        method = 'bray',
        by = NULL)

# macrofauna rda----
set.seed(123)
macro_rda <- dbrda(t_fauna ~ .,data = env2)
saveRDS(macro_rda, 'data/output_data/macro_rda.RDS')

# macro_ordisted <- ordistep(macro_rda)
# summary(macro_ordisted, axes = 0)

RsquareAdj(macro_rda)$adj.r.squared
RsquareAdj(macro_rda)$r.squared
anova(macro_rda, permutations = how(nperm = 999))
macro_ordisted

# Eukaryote RDA biplot----------
macro_rda_plot <- 
  ggord(
  macro_rda,
  grp_in = env$Treatment,
  poly = F,
  repel = T,
  txt = 3,
  addsize = NA,
  ellipse = F,
  arrow = 0,
  alpha = .7,
  veclsz = .2,
  vec_ext = 0.5,
  grp_title = 'Treatment'
) +
  theme_javier(base_size = 10) +
  labs(subtitle = 'a. Macrofauna') +
  scale_color_viridis_d() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"), aspect.ratio = 1)


# Bacteria- ------
bact_rare_log_mean  <-  
  readRDS("data/output_data/bact_rare_log.RDS") %>% 
  merge_samples(., group = "Plot", fun = mean)

bact_rare_log_mean_df  <-
  as(otu_table(bact_rare_log_mean), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column('Plot') %>% 
  mutate(Plot = as.numeric(sub("_.*", "", Plot))) %>% 
  arrange(Plot)

# macrofauna distlm-----
distlm_bact_margin <- 
  adonis2(
    bact_rare_log_mean_df ~ ., data = env2 ,
    method = 'bray',
    by = 'margin',
    permutations = 999
  )

distlm_bact_margin

adonis2(
  bact_rare_log_mean_df ~ ., data = env2 ,
  method = 'bray',
  by = NULL,
  permutations = 999
)

# bacteria dbRDA-------------
bact_rda <- dbrda(bact_rare_log_mean_df ~ . , data = env2)
saveRDS(bact_rda, 'data/output_data/bact_rda.RDS')

anova(bact_rda)

# bact_ordisted <- ordistep(bact_rda)

# bacteria dbRDA biplot -------------
bact_rda_plot <- 
  ggord(
  bact_rda,
  grp_in = env$Treatment,
  poly = F,
  repel = T,
  txt = 3,
  addsize = NA,
  ellipse = F,
  arrow = 0,
  alpha = .7,
  veclsz = .2,
  vec_ext = 0.5,
  grp_title = 'Treatment'
) +
  theme_javier(base_size = 10) +
  labs(subtitle = 'c. Bacteria') +
  scale_color_viridis_d() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"), aspect.ratio = 1)


# Eukaryotes -------
euk_rare_log <- readRDS("data/euk_rare_log.RDS")
sample_data(euk_rare_log)$Plot <- as.character(sample_data(euk_rare_log)$Plot) 
euk_rare_log_mean <- merge_samples(euk_rare_log, group = "Plot", fun = mean)

euk_rare_log_mean_df  <-
  as(otu_table(euk_rare_log_mean), "matrix") %>%
  as.data.frame() %>%
  rownames_to_column('Plot') %>% 
  mutate(Plot = as.numeric(sub("_.*", "", Plot))) %>% 
  arrange(Plot)

# eukaryote distlm-----
distlm_euk_margin <- 
  adonis2(
    euk_rare_log_mean_df ~ ., data = env2 ,
    method = 'bray',
    by = 'margin',
    permutations = 999
  )

adonis2(
  euk_rare_log_mean_df ~ ., data = env2 ,
  method = 'bray',
  by = NULL
)

# Eukaryote RDA-----
euk_rda <- dbrda(euk_rare_log_mean_df ~. , env2)
saveRDS(euk_rda, 'data/output_data/euk_rda.RDS')

# euk_ordisted <- ordistep(euk_rda)
# euk_ordisted
# anova(euk_ordisted)


## dbRDA biplot euk------------
euk_rda_plot <-
  ggord(
    euk_rda,
    grp_in = env$Treatment,
    poly = F,
    ellipse = F,
    repel = T,
    txt = 3,
    addsize = NA,
    arrow = 0,
    alpha = .7,
    veclsz = .2,
    vec_ext = 0.5,
    grp_title = 'Treatment'
  ) +
  theme_javier(base_size = 10) +
  labs (subtitle = "b. Eukaryotes") +
  scale_color_viridis_d() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "lines"), aspect.ratio = 1)
