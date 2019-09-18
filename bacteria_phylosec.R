# Analyses of bacterial data from the field Sabella experiment 

# 01 load libraries--------------
library(phyloseq)
library(tidyverse)
library(ggvegan)
library(vegan)
library(ggrepel)
library(emmeans)
source('theme_javier.R')
theme_set(theme_javier())

# read OTU table with OTU > 5 reads--------------
otu_bact <-
  read.csv("data/bacteria_data.csv",
           row.names = 1) %>%
  select(-"SS.Ext.Blank_S94",-"SS.PCRBlank_S95" ,-"SS.water_S96") %>%
  as.matrix(.) %>% 
  otu_table(., taxa_are_rows = TRUE)

# taxonomy -------
tax_bact <-
  read.csv('data/taxonomy_bacteria.csv', row.names = 1) %>%
  as.matrix(.) %>%
  tax_table(.)

# read sample data-------
sample_bact <-
  read.csv('data/molecular_treatments.csv', row.names = 1) %>%
  # mutate(Plot = as.character(Plot)) %>% 
  sample_data(.)

# merge data ---------
bact_dat <- 
  phyloseq(otu_bact, tax_bact, sample_bact) %>% 
  subset_samples(Treatment != "Ambient") %>% 
  filter_taxa(function(x) mean(x) > 1e-5, TRUE) %>% # remove OTU with very small no. of reads
  subset_taxa(Phylum!="0" & Family !='Chloroplast') # remove Phylum 0 and Chloroplast

# pre-process-------------
table(tax_table(bact_dat)[, "Phylum"], exclude = NULL) # remove unassigned Phylum
ntaxa(bact_dat)

# An option is to combine replicates using 
# merge_samples(bact_dat, "Replicate", fun=mean)

# rarerefy data-----------
set.seed(1234)

# get the total read abundace for each sample in descending order
samp_sum <- 
  sample_sums(bact_dat) %>% 
  data_frame() %>% 
  rename(N = '.') %>% 
  arrange(N) %>% 
  print()

# rarefaction curves----
rarecurve(
  t(otu_table(bact_dat)),
  step = 50,
  cex = 0.5,
  label = F,
  main = 'Bacteria samples',
  ylab = 'AVS'
)
abline(v = 5000, lty = 2)


min(sample_sums(bact_dat))

bact_rare  <- rarefy_even_depth(bact_dat,
                                sample.size = 5000,
                                rngseed = 123)
# 2 samples removed because they contained fewer reads than `sample.size`.
# SS.Bac.13c_S39
# SS.Bac.22c_S66

bact_rare_log <- transform_sample_counts(bact_rare, function(x) log10(x+1))

saveRDS(bact_rare_log, "data/output_data/bact_rare_log.RDS")
saveRDS(bact_rare, "data/output_data/bact_rare.RDS")

## filter data---
# bact_rare_filt  <- filter_taxa(bact_rare, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
# bact_rare_filt_log <- transform_sample_counts(bact_rare_filt, function(x) log10(x+1))

# check data names and taxonomy descrition -------
ntaxa(bact_rare_log)
nsamples(bact_rare_log)
rank_names(bact_rare_log)

length(get_taxa_unique(bact_rare_log, "Phylum")) # 44
length(get_taxa_unique(bact_rare_log, "Class")) # 92
length(get_taxa_unique(bact_rare_log, "Sub_class")) # 206
length(get_taxa_unique(bact_rare_log, "Family"))# 281



# PERMANOVA------------------
library(vegan)
taxono_bact <- 
  read_csv('data/taxonomy_bacteria.csv')

metadata_bact <- as(sample_data(bact_rare_log), "data.frame")

permanova_bact <-
  adonis(distance(bact_rare_log, method = "bray") ~ Treatment + final_density,
         data = metadata_bact, permutations = 999)
permanova_bact

bact_permanova_table <- 
  data.frame(permanova_bact$aov.tab) %>% 
  rownames_to_column() %>% 
  dplyr::rename(Terms = rowname,
         df = Df,
         SS = SumsOfSqs,
         MS = MeanSqs,
         "Pseudo-F" = F.Model,
         P = Pr..F.)


# 08 PERMANOVA pairwise comparisons-------------
adonis(
  distance(subset_samples(
    bact_rare_log, Treatment %in% c("Sabella", "Control")
  ), method = "bray") ~ Treatment ,
  data = filter(metadata_bact, Treatment != "Mimic"),
)

adonis(
  distance(subset_samples(
    bact_rare_log, Treatment %in% c("Mimic", "Control")
  ), method = "bray") ~ Treatment ,
  data = filter(metadata_bact, Treatment != "Sabella")
)

adonis(
  distance(subset_samples(
    bact_rare_log, Treatment %in% c("Sabella", "Mimic")
  ), method = "bray") ~ Treatment ,
  data = filter(metadata_bact, Treatment != "Control")
)


# CAP analysis---------
bact_CAP <-
  ordinate(
    bact_rare_log,
    "CAP",
    distance =  'bray',
    formula = ~ Treatment + final_density
  )


cap_anova_bact <- anova(bact_CAP, by = 'terms', permutations = 999)
cap_anova_bact

bact_cap_table <- 
  data.frame(cap_anova_bact) %>% 
  rownames_to_column() %>% 
  dplyr::rename(Terms = rowname,
                df = Df,
                SS = SumOfSqs,
                P = Pr..F.)

metadata_bact <- as(sample_data(bact_rare_log), "data.frame")

cap_bact_dat <- fortify(bact_CAP)

cap_bact_sites <- 
  cap_bact_dat %>% 
  dplyr::filter(Score == 'sites') %>% 
  bind_cols(metadata_bact) %>% 
  select(-Score,-Label) %>% 
  write_csv('data/output_data/cap_bact_sites.csv')

cap_bact_taxa <-
  cap_bact_dat %>%
  dplyr::filter(Score == 'species') %>%
  dplyr::filter(abs(CAP1) > .075 | abs(CAP2) > .075) %>%
  dplyr::rename(OTU = "Label") %>% 
  left_join(taxono_bact, by = "OTU")

# CAP biplot------
cap_bact_plot <- 
  ggplot(cap_bact_sites,
       aes(
         x = CAP1,
         y =  CAP2
       )) +
  geom_point(alpha = .7, aes(color = Treatment, 
                             size = Density)) + 
  scale_radius(limits = c(0, 50),
               range = c(3, 7),
               name = 'Density') +
  scale_color_viridis_d(name = "Treatment",
                        option = 'D')

cap_bact_biplot <-
  cap_bact_plot +
  geom_text_repel(
    data = cap_bact_taxa,
    aes(label = Family, x = CAP1*10 , y =  CAP2*10 ),
    size = 4
  )

cap_bact_biplot

dat_bact_raw <- psmelt(bact_rare) 


# # agglomarate at family level------------
# bact_rare_fam  <- tax_glom(bact_rare, "Family")
# physeq_bact_fam  <- tax_glom(physeq_bact2, "Family")
physeq_bact_subc  <- tax_glom(physeq_bact2, "Sub_class")
physeq_bact_class  <- tax_glom(bact_rare, "Class")
physeq_bact_phyl  <- tax_glom(bact_rare, "Phylum")

# Indicator spp--------
bact_rare_filt <- filter_taxa(bact_rare, function(x) mean(x) > .9, TRUE)
# bact_rare_filt  <- filter_taxa(bact_rare, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
bact_rare_df  <-
  as(otu_table(bact_rare_filt), "matrix") %>%
  t(.) %>% 
  as.data.frame()

ind_sp_bact <- multipatt(bact_rare_df, metadata_bact$Treatment, duleg=TRUE)

summary(ind_sp_bact, indvalcomp = TRUE)


# indicator taxa table--------------
taxo_table_bact <- 
  read.csv('data/taxonomy_bacteria.csv') %>% 
  rename(OTU =  "ï..OTU")

ind_bact_table <- 
  ind_sp_bact$sign %>% 
  data.frame() %>% 
  rownames_to_column('OTU') %>% 
  filter(p.value<0.05) %>% 
  mutate(Treatment = fct_recode(
    factor(index),
    Sabella = "3",
    Mimic = "2",
    Control = "1"
  )) %>%
  select(Treatment, OTU, stat, p.value) %>%
  arrange(Treatment, p.value) %>% 
  left_join(taxo_table_bact) %>% 
  write_csv('tables/inc_bact_table.csv')


ind_ASV_bact <- 
  ind_bact_table %>%
  select(OTU) %>% 
  as.list()

# relative abundance plot of indicator taxa -----
bact_rare_rel_df <- 
  bact_rare %>% 
  filter_taxa(function(x) mean(x) > 1, TRUE) %>% 
  transform_sample_counts(.,function(x) (x/sum(x))*100) %>% 
  psmelt(.) 

bact_rare_melt <- 
  bact_rare %>% 
  filter_taxa(function(x) mean(x) > 1, TRUE) %>% 
  psmelt(.) 

ind_bact_plot_data <- 
  bact_rare_melt %>%
  filter(OTU %in% ind_ASV_bact$OTU) %>%
  mutate(Taxa = Phylum) %>% 
  write_csv('data/output_data/ind_bact_plot_data.csv')

ind_bact_plot <- 
  ggplot(ind_bact_plot_data, aes(x = Phylum,
             y = Abundance,
             color = Treatment)) +
  stat_summary(fun.data = mean_se,
               size = .5,
               position = position_dodge(width = .5), alpha = 0.5) +
  labs(y = '', x = '') +
  coord_flip() +
  scale_color_viridis_d() +
  theme(legend.position = c(.8, .2))
  # scale_y_sqrt(breaks = c(1,5, 10, 20, 30))
ind_bact_plot

# SIMPER--------------------
# ## Simper analysis, extract abundance matrix from a phyloseq object
# bact_rare_phyl_df  <-
#   as(otu_table(physeq_bact_phyl), "matrix") %>%
#   t(.) %>% 
#   as.data.frame()
# 
# # running the simper analysis on the dataframe and the variable of interest "time"
# simper_bact_fam <- simper(bact_rare_phyl_df, metadata_bact$Treatment, permutations = 100)
# 
# simper_bact <-
#   bind_rows("Sabella vs. Control " = simper_bact_fam$Sabella_Control %>% data.frame(),
#             "Mimic vs. Control" = simper_bact_fam$Mimic_Control %>% data.frame(),
#             "Mimic vcs. Sabella" = simper_bact_fam$Mimic_Sabella %>% data.frame(),
#             .id = "Groups") %>% 
#   dplyr::filter(cusum<.7) %>%
#   dplyr::rename(OTU = "species") %>%
#   left_join(taxono_bact, by = "OTU") %>%
#   dplyr::rename(Taxa = "Phylum") %>%
#   dplyr::select(Groups, Taxa, average, sd, ratio, ava, avb, cusum, p) %>%
#   print(digits = 2)


# # relative abundance plots -----
# bact_rare_rel <- 
#   bact_rare %>% 
#   filter_taxa(function(x) mean(x) > 1, TRUE) %>% 
#   transform_sample_counts(. ,function(x) (x/sum(x))*100) 
# 
# dat_bact <- 
#   psmelt(bact_rare_rel) %>% 
#   write_csv('data/bact_rel_abund.csv')
# 
# # Phylum level------------------
# ggplot(dat_bact %>% filter(Phylum != 0),
#        aes(
#          x = Phylum,
#          y = Abundance,
#          color = Treatment
#        )) +
#   stat_summary(fun.data = mean_se, size = .5,
#                position = position_dodge(width = .5)) +
#   scale_color_discrete() +
#   labs(y = 'Abundance') +
#   coord_flip() +
#   theme(legend.position = c(.8,.8))


# Diversity indices bacteria----------
dat_bact_raw <- 
  bact_rare %>% 
  psmelt(.)

bact_indices <-
  dat_bact_raw %>%
  group_by(Sample) %>%
  mutate(
    # N = sum(Abundance),
    H = diversity(Abundance),
    S = specnumber(Abundance),
    J = H / log(S)
  )

treat_bact_df <- sampledata %>% mutate(Sample = rownames(.))

bact_indices_raw <- 
  bact_indices %>%
  group_by(Sample) %>%
  summarise(S = mean(S),
            J = mean(J)) %>%
  left_join(treat_bact, by = 'Sample') %>% 
  write_csv('data/output_data/bact_indices_raw.csv')

# N_bact <- 
#   ggplot(bact_indices_raw, aes(x = Treatment, y = N, fill = Treatment)) +
#   geom_boxplot(alpha = .8) +
#   scale_fill_viridis_d(guide = F, discrete = T) +
#   labs(x='', y = 'Number of individuals') + scale_y_log10()

S_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = S, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis_d(guide = F) +
  labs(x='', y = 'Number of taxa')

J_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = J, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis_d(guide = F) +
  labs(x='', y = 'Community evenness')


# anova indices---
# ancova indices------
indices_bact_dat <-
  bact_indices_raw %>%
  gather(index, value, S) %>%
  group_by(index) %>%
  nest() %>%
  mutate(
    ancova = map(.x = data, ~ lm(value ~ Treatment  + final_density, data = .x)),
    ancova_sum = map(ancova, summary),
    anova_tab = map(ancova, ~car::Anova(., type = "II")),
    ancova_table = map(anova_tab, broom::tidy),
    residuals = map(ancova, broom::augment),
    emm = map(ancova,~emmeans(., ~Treatment)),
    posthoc = map(emm,pairs))

ancova_table_bact <- 
  indices_bact_dat %>% 
  select(index,ancova_table) %>% 
  unnest()

indices_bact_dat %>%
  select(index, posthoc) %>%
  flatten()

## model validation----
res_bact <- 
  indices_bact_dat %>% 
  select(index, residuals) %>% 
  unnest(residuals, .drop = TRUE)

# Plot of fitted vs. residual values by index to check the assumptions heteroscedasticity or any other pattern,
#  e.g. non-linearity or missing covariates
ggplot(res_bact) +
  geom_point(aes(x = .fitted, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2)

# boxplot by type
ggplot(res_bact) +
  geom_boxplot(aes(x = Treatment, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2)

# qqplot of the normalised residuals to check the assumption of normality------------
ggplot(res_bact) +
  stat_qq(aes(sample = .std.resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_abline(
    intercept =  0,
    slope = 1,
    lty = 2,
    col = 2
  )
