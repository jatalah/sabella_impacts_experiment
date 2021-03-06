# 01 load libraries ----------------
suppressWarnings(suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(vegan)
  library(ggord)
  library(ggpubr)
  library(ggvegan)
  library(ggrepel)
  library(emmeans)
}))
source('theme_javier.R')
theme_set(theme_javier())
set.seed(1234)
# rm(list = ls())

# 02 read and prep data -----------------------------------
macrofauna_dat <-
  read_excel('data/macrofauna_data.xlsx', sheet = 'R') %>%
  dplyr::filter(Sampling == 2 & Treatment != "Ambient")# excluding ambient
 
macrofauna_dat[is.na(macrofauna_dat)] <- 0
# read large macrofaunal data 

# large_fauna <-
#   read_excel('data/sabella_large_macrofauna_data.xlsx', sheet = 'R') %>%
#   dplyr::select(-c(code, Treatment:Density))
# 
# large_fauna[is.na(large_fauna)] <- 0

# filter taxa with zero abundance 
present_taxa <- 
  macrofauna_dat %>% 
  gather(taxa, value, 7:100) %>% 
  group_by(taxa) %>% 
  summarise_at(vars(value), funs(sum), na.rm =T) %>%
  arrange(desc(value)) %>% 
  dplyr::filter(value>0) %>% 
  dplyr::select(taxa) %>% 
  as.list()

# treat_dat <-
#   macrofauna_dat %>%
#   dplyr::select(code:Density)

treat_dat <- 
  read_csv('data/treatments.csv') %>%
  dplyr::select(-weight)


# impute missing value for mimic final_density usinbg the mean ----
treat_dat %>% 
  filter(Treatment=='Mimic') %>% 
  group_by(Density) %>% 
  summarise(mean(final_density, na.rm = T))

# use 12.5 as estimate of missing value------------
treat_dat <- 
  treat_dat %>% 
  mutate(final_density = ifelse(is.na(final_density), 12.5, final_density))


fauna_dat <- 
  macrofauna_dat %>% 
  dplyr::select(sample, present_taxa$taxa) %>% 
  full_join(treat_dat, by = 'sample') %>% 
  dplyr::select(sample, Treatment, Density, final_density, everything()) %>% 
  arrange(sample) %>% 
  write_csv('data/fauna_dat.csv')


# transform data----------------
log_fauna <- log10(fauna_dat[,-c(1:4)]+1)
fauna_4th_root <- sqrt(sqrt(fauna_dat[,-c(1:4)]))

# calculate dissimmilarity matrices--------
d_4th_root_fauna<- vegdist(fauna_4th_root , method = 'bray')
d_log_fauna<- vegdist(log_fauna, method = 'bray')

# 03 mean + CI plots of dominant taxa-----------
long_fauna_dat <- 
  fauna_dat %>% 
  gather(taxa, abun,`Theora lubrica`:`Taeniogyrus dendyi`)


taxa_bar_plot <- 
  long_fauna_dat %>%
  group_by(taxa) %>%
  filter(sum(abun, na.rm = T) > 8) %>%
  ungroup() %>%
  ggplot(aes(x = Treatment, y = abun, fill = Treatment)) +
  # stat_summary(fun.data = mean_se, size = .5) +
  stat_summary(fun.y = mean,
               geom = "bar",
               color = 1) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.4) +
  facet_wrap(~ taxa, scales = 'free_y') +
  scale_fill_viridis_d() +
  theme_javier() +
  labs('x = Treatment', y = 'Abundance')

ggsave(
  taxa_bar_plot,
  filename = 'figures/bar_plot_macrofauna_taxa.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 6,
  width = 10
)



# 04 MDS ordination -------------
mds <- metaMDS(fauna_4th_root, distance = 'bray')

mds_dat <-
  fortify(mds) %>%
  rename(sample = Label) %>%
  dplyr::filter(Score == 'sites') %>% 
  bind_cols(treat_dat)


# get taxa with higher scores 
vec_taxa <-
  mds$species[, 1:2] %>% 
  data.frame %>% 
  mutate(Taxa = rownames(.)) %>% dplyr::filter(abs(MDS1) >.8 | abs(MDS2) >.8) %>% 
  dplyr::select(Taxa) %>% 
  as.list()

## MDS biplot 
ggord(
  mds,
  grp_in = fauna_dat$Treatment,
  poly = F,
  alpha = 1,
  ellipse = T,
  # arrow = .3,
  arrow = 0,
  repel = T,
  text = .01,
  vec_ext = 1,
  # size = fauna_dat$Density
  var_sub = vec_taxa$Taxa
) +
  theme_javier() +
  scale_shape_manual(values = 15:19) 

# 05 PCO ordination --------------------------------------------------
pco <-
  capscale(
    log10(fauna_dat[, -c(1:4)] + 1) ~ 1 ,
    data = fauna_dat,
    metaMDSdist = TRUE,
    sqrt.dist = T,
    add = T,
    distance = "bray",
    permutations = 1e6,
    na.action = na.omit
  )
pco_dat <- fortify(pco)

pco_ordination <-
  pco_dat %>%
  dplyr::filter(Score == 'sites') %>%
  ggplot(.,
         aes(
           x = MDS1,
           y =  MDS2,
           color = fauna_dat$Treatment,
           # shape = fauna_dat$Treatment,
           size = fauna_dat$Density
         )) +
  geom_point(alpha = .7) +
  # scale_color_discrete(name = "Treatment") +
  scale_shape_manual(values = 15:18, name = "Treatment") +
  labs(x = "PCO1", y = "PCO2")+
  scale_radius(limits = c(0, 50),
               range = c(3, 7),
               name = 'Density') +
  scale_color_viridis_d(name = NULL,
                      option = 'D')
pco_ordination

# 06 CAPscale-----------------------
cap <-
  capscale(
    log(fauna_dat[, -c(1:4)]+1) ~ Treatment + final_density,
    # sqrt(sqrt(fauna_dat[, -c(1:4)])) ~ Treatment * final_density,
    data = fauna_dat,
    sqrt.dist = T,
    add = F,
    distance =  "bray",
    permutations = 1e6,
    na.action = na.omit
  )

# anova table Cap
cap_anova_macrofauna <- anova(cap, permutations = 999, by = 'terms')
cap_anova_macrofauna


macrofauna_cap_table <- 
  data.frame(cap_anova_macrofauna) %>% 
  rownames_to_column() %>% 
  rename(Terms = rowname,
         df = Df,
         SS = SumOfSqs,
         P = Pr..F.)

# get data for plots 
cap_dat <- fortify(cap) 

cap_sites <- 
  cap_dat %>% 
  dplyr::filter(Score == 'sites') %>% 
  bind_cols(treat_dat)

cap_taxa <-
  cap_dat %>%
  dplyr::filter(Score == 'species') %>%
  dplyr::filter(abs(CAP1) > .15 | abs(CAP2) > .15) %>%
  mutate(plot_label = if_else(
    condition = str_detect(Label, " "),
    true = paste0('italic(', Label, ')'),
    false = as.character(Label)
  )) %>%
  mutate(plot_label = str_replace(plot_label," ","~"))

# CAP biplot-------
cap_plot <-
  ggplot(cap_sites,
         aes(
           x = CAP1,
           y =  CAP2
         )) +
  geom_point(alpha = .7, aes(color = Treatment, 
                             # shape = fauna_dat$Treatment,
                             size = Density)) + 
  # scale_color_discrete(name = "Treatment") +
  # scale_shape_manual(values = 15:17, name = "Treatment") +
  scale_radius(limits = c(0, 50),
               range = c(3, 7),
               name = 'Density') +
  scale_color_viridis_d(name = "Treatment",
                      option = 'D')
cap_plot

# biplot with taxa scores
cap_biplot <-
  cap_plot +
  geom_text_repel(
    data = cap_taxa,
    aes(label = plot_label, x = CAP1 * 4, y =  CAP2 * 4),
    size = 4,
    parse = TRUE
  )
cap_biplot

ggsave(
  cap_biplot,
  filename = 'figures/CAP_biplot.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 4,
  width = 5
)

# Both plot in one
ggarrange(pco_ordination, cap_plot, common.legend = T, labels = "auto")
ggsave(
  last_plot(),
  filename = 'figures/PCO_CAP.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 4,
  width = 8
)

# 07 PERMANOVA--------------------------------------
permanova <-
  adonis(
    log(fauna_dat[, -c(1:4)]+1) ~  Treatment + final_density,
    data = fauna_dat,
    method =  "bray",
    permutations = 9999
  )
permanova

macrofauna_permanova_table <- 
  data.frame(permanova$aov.tab) %>% 
  rownames_to_column() %>%
  rename(Terms = rowname,
         df = Df,
         SS = SumsOfSqs,
         MS = MeanSqs,
         "Pseudo-F" = F.Model,
         P = Pr..F.)



# 08 PERMANOVA pairwise comparisons-------------
source('pairwise.adonis.R')
set.seed(99) 
# .Random.seed
pairwise.adonis(
  log10(fauna_dat[, -c(1:4)]+1),
  # sqrt(sqrt(fauna_dat[, -c(1:4)] )),
                fauna_dat$Treatment,
                sim.method = 'bray')

# pairs  F.Model         R2 p.value p.adjusted
# 1 Control vs Sabella 1.484889 0.11011501   0.091      0.273
# 2   Control vs Mimic 1.234271 0.09326323   0.204      0.612
# 3   Sabella vs Mimic 1.486769 0.07629632   0.046      0.138

# 09 SIMPER ----
simp <- simper(fauna_dat[, -c(1:4)], fauna_dat$Treatment, permutations = 999)
# simp <- simper(fauna_dat[,-c(1:4)], fauna_dat$Treatment, permutations = 999)
summary(simp, digits = 2, ordered = T)

simper_macrofauna <-
  bind_rows("Control vs. Sabella" = simp$Control_Sabella %>% data.frame(),
            "Control vs.Mimic" = simp$Control_Mimic %>% data.frame(),
            "Sabella vs. Mimic" = simp$Sabella_Mimic %>% data.frame(),
            .id = "Groups") %>% 
  dplyr::select(-overall, -ord) %>% 
  dplyr::filter(cusum<.7) %>% 
  rename(Taxa = "species") %>% 
  print(digits = 2)


# 10 PERMDISP ----------------------
dis <- vegdist(sqrt(sqrt(fauna_dat[,-c(1:3)] )))
beta_disp <-
  betadisper(
    dis, fauna_dat$Treatment
  )
beta_disp
anova(beta_disp)
permutest(beta_disp, pairwise = TRUE, permutations = 999)
plot(beta_disp)
boxplot(beta_disp)

# 11 Univariate diversity index -------------------
indices <- 
  fauna_dat %>% 
  transmute(N = rowSums(fauna_dat[, -c(1:4)]),
            H = diversity(fauna_dat[, -c(1:4)]),
            S = specnumber(fauna_dat[, -c(1:4)]),
            J = H/log(S),
            F = fisher.alpha(fauna_dat[, -c(1:4)]),
            ES = rarefy(fauna_dat[, -c(1:4)], min(N))) %>% 
  bind_cols(treat_dat) %>% 
  write_csv('data/macro_indices.csv')

boxplots <- 
  indices %>% 
  gather(index, value, c("N", "S", "J")) %>% 
  mutate(index = fct_relevel(index, c("N", "S", "J"))) %>% 
  ggplot(., aes(x = Treatment, y = value, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  facet_wrap(~index,scales = 'free') +
  scale_fill_viridis_d(guide = F)
boxplots


# ancova indices------
indices_dat <-
  indices %>%
  mutate(N = log(N +1)) %>% 
  gather(index, value, c("N", "S", "J")) %>%
  group_by(index) %>%
  nest() %>%
  mutate(
    ancova = map(.x = data, ~ aov(value ~ Treatment + final_density, data = .x)),
    anova_tab = map(ancova, ~car::Anova(., type = 2)),
    ancova_table = map(anova_tab, broom::tidy),
    residuals = map(ancova, broom::augment),
    emm = map(ancova,~emmeans(., ~Treatment)),
    posthoc = map(emm,pairs))


ancova_table_macro <- 
  indices_dat %>% 
  select(index,ancova_table) %>% 
  unnest()

indices_dat %>% 
  select(index,posthoc) %>% 
  flatten()

## model validation----
res <- 
  indices_dat %>% 
  select(index, residuals) %>% 
  unnest(residuals, .drop = TRUE)

# Plot of fitted vs. residual values by index to check the assumptions heteroscedasticity or any other pattern,
#  e.g. non-linearity or missing covariates
ggplot(res) +
  geom_point(aes(x = .fitted, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2)

# boxplot by type
ggplot(res) +
  geom_boxplot(aes(x = Treatment, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2)

# qqplot of the normalised residuals to check the assumption of normality------------
ggplot(res) +
  stat_qq(aes(sample = .std.resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_abline(
    intercept =  0,
    slope = 1,
    lty = 2,
    col = 2
  )
