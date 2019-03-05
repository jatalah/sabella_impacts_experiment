# 01 load libraries ----------------
suppressWarnings(suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(vegan)
  library(ggord)
  library(ggpubr)
  library(ggvegan)
  library(ggrepel)
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
d_4th_root_fauna<- vegdist(fauna_4th_root , method = 'bray')

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
    sqrt(sqrt(fauna_dat[, -c(1:4)])) ~ Treatment * final_density,
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

# CAP biplot
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
    size = 3,
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
    sqrt(sqrt(fauna_dat[, -c(1:3)])) ~ Treatment * final_density,
    data = fauna_dat,
    method =  "bray",
    permutations = 9999
  )
permanova


data_mimic <- filter(fauna_dat, Treatment == 'Control' | Treatment == "Mimic")
data_sabella <- filter(fauna_dat, Treatment == 'Control' | Treatment == "Sabella")

adonis(
  sqrt(sqrt(data_mimic[, -c(1:4)])) ~ final_density,
  data = data_mimic,
  method =  "bray",
  permutations = 999
)

adonis(
  sqrt(sqrt(data_sabella[, -c(1:4)])) ~ final_density,
  data = data_sabella,
  method =  "bray",
  permutations = 999
)
 


# 08 PERMANOVA pairwise comparisons-------------
source('H:/javiera/Stats/R folder/pairwise.adonis.R')

pairwise.adonis(sqrt(sqrt(fauna_dat[, -c(1:3)] )),
                fauna_dat$Treatment,
                sim.method = 'bray')

# pairs  F.Model         R2 p.value p.adjusted
# 1   Control vs Mimic 1.311313 0.09851117   0.148      0.444
# 2 Control vs Sabella 1.571497 0.11579390   0.065      0.195
# 3   Mimic vs Sabella 1.537155 0.07867856   0.046      0.138

# 09 SIMPER ----
simp <-
  simper(fauna_4th_root, fauna_dat$Treatment, permutations = 999)
summary(simp, digits = 2, ordered = T)

sim_mc_cont_sabell <-
  simp$Control_Sabella %>%
  data.frame() %>%
  filter(p < 0.1) %>%
  arrange(p) %>%
  dplyr::select(species, ratio, ava, avb, cusum, p) %>%
  rename(Taxa = "species") %>% 
  print(digits = 2)

sim_mc_cont_mimic <-
  simp$Control_Mimic %>%
  data.frame() %>%
  filter(p < 0.2) %>%
  arrange(p) %>%
  dplyr::select(species, ratio, ava, avb, cusum, p) %>%
  rename(Taxa = "species") %>% 
  print(digits = 2)

simper_macrofauna <-
  bind_rows("Sabella vs. Control" = sim_mc_cont_sabell,
            "Mimic vs. Control" = sim_mc_cont_mimic,.id = "Groups")

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
  write_csv(., 'data/macro_indices.csv')

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
    ancova = map(.x = data, ~ lm(value ~ Treatment * final_density, data = .x)),
    anova_tab = map(ancova, anova),
    ancova_table = map(anova_tab, broom::tidy),
    residuals = map(ancova, broom::augment)
  ) 

ancova_macro <- 
  indices_dat %>% 
  select(index,ancova_table) %>% 
  unnest()



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

# ## other approaches----------
# library(caret)
# ancova_S <- 
#   train(
#   S ~ Treatment * final_density,
#   method = 'lm',
#   data = indices,
#   preProcess = "BoxCox",
#   trControl = trainControl(method = "none")
# )
# 
# anova(lm(log(N)~Treatment*final_density, indices))
# anova(lm(S~Treatment*final_density, indices))
# anova(lm(J~Treatment*final_density, indices))
# 
# # regression plots

# library(ggpmisc)
# lm_plots <-
#   indices %>%
#   filter(Treatment == 'Control' | Treatment == "Mimic") %>%
#   gather(index, value, c("N", "S", "J")) %>%
#   mutate(index = fct_relevel(index, c("N", "S", "J"))) %>%
#   ggplot(aes(x = final_density, y = value)) +
#   geom_point(alpha = .8) +
#   facet_wrap( ~ index, scales = 'free') +
#   scale_fill_viridis_d(guide = F) +
#   geom_smooth(method = 'lm') +
#   stat_poly_eq(formula = y~x, parse = T)
# 
# 
# indices %>%
#   filter(Treatment == 'Control' | Treatment == "Sabella") %>%
#   gather(index, value, c("N", "S", "J")) %>%
#   mutate(index = fct_relevel(index, c("N", "S", "J"))) %>%
#   ggplot(aes(x = final_density, y = value)) +
#   geom_point(alpha = .8) +
#   facet_wrap( ~ index, scales = 'free') +
#   scale_fill_viridis_d(guide = F) +
#   geom_smooth(method = 'lm') +
#   stat_poly_eq(formula = y~x, parse = T)
# 
# indices %>%
#   filter(Treatment == 'Control' | Treatment == "Mimic") %>%
#   lm(S~final_density, data = .) %>% 
#   summary(.)
# 
# # 12 Diversity including excavated samples --------------------------
# all_fauna_dat <- 
#   bind_rows(fauna_dat[,-c(2:3)], large_fauna) %>% 
#   group_by(sample) %>% 
#   summarise_all(sum, na.rm = T) %>% 
#   left_join(., treat_dat) %>% 
#   dplyr::select(sample,Treatment, Density, everything(), - (code:Treat))
# 
# indices_all <- 
#   all_fauna_dat %>% 
#   transmute(N = rowSums(all_fauna_dat[, -c(1:3)]),
#             H = diversity(all_fauna_dat[, -c(1:3)]),
#             S = specnumber(all_fauna_dat[, -c(1:3)]),
#             J = H/log(S),
#             F = fisher.alpha(all_fauna_dat[, -c(1:3)]),
#             ES = rarefy(all_fauna_dat[, -c(1:3)], min(N))) %>% 
#   bind_cols(treat_dat)
# 
# boxplots_all <- 
#   indices_all %>% 
#   gather(index, value, c("N", "S", "J", "ES")) %>% 
#   mutate(index = fct_relevel(index, c("N", "S", "J"))) %>% 
#   ggplot(., aes(x = Treatment, y = value, fill = Treatment)) +
#   geom_boxplot(alpha = .7) +
#   facet_wrap(~index,scales = 'free') +
#   scale_fill_viridis(guide = F, discrete = T)
# boxplots_all
