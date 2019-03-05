# Analyses of eukaryote data from the field Sabella experiment 

# 01 load libraries--------------
library(phyloseq)
library(tidyverse)
library(ggvegan)
library(vegan)
library(ggrepel)
source('theme_javier.R')
theme_set(theme_javier())

# 02 read feature table and taxonomic table--------------
otumat <- read.csv('data/eukaryote_data.csv', row.names = 1) %>% as.matrix(.)
taxmat<- read.csv('data/taxonomy_eukaryote.csv',row.names = 1 ) %>% as.matrix(.)

OTU  <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)

# 03 join and create phylosec object
physeq  <- phyloseq(OTU, TAX)

# read sample data-------
sampledata <- 
  read.csv('data/molecular_treatments.csv', row.names = 1) %>% 
  sample_data(.)

# merge data ---------------
physeq1  <- merge_phyloseq(physeq, sampledata)

# 04 pre-process-------------
# remove ambient data---
physeq3 <- subset_samples(physeq1, Treatment != "Ambient")
physeq4 <- filter_taxa(physeq3, function(x) mean(x) > 1e-5, TRUE)

# check data structur
ntaxa(physeq4)
nsamples(physeq4)
rank_names(physeq4)

# agglomarate at family level------------
physeq_fam  <- tax_glom(physeq4, "Family")
physeq_order  <- tax_glom(physeq4, "Order")
physeq_class  <- tax_glom(physeq4, "Class")

# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
physeq_filt <- filter_taxa(physeq4, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_fam_filt  <- filter_taxa(physeq_fam, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_order_filt  <- filter_taxa(physeq_order, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_class_filt  <- filter_taxa(physeq_class, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# Transformation ----
# filtered data
physeq_filt_log <- transform_sample_counts(physeq_filt, function(x) log(x+1))
physeq_fam_filt_log <- transform_sample_counts(physeq_fam_filt, function(x) log(x+1))
physeq_fam_filt_4throot <- transform_sample_counts(physeq_fam_filt, function(x) sqrt(sqrt(x)))
physeq_order_filt_log <- transform_sample_counts(physeq_order_filt, function(x) log(x+1))
physeq_class_filt_log <- transform_sample_counts(physeq_class_filt, function(x) log(x+1))

# unfiltered data
physeq_log <- transform_sample_counts(physeq4, function(x) log(x+1))
physeq_fam_log <- transform_sample_counts(physeq_fam, function(x) log(x+1))
physeq_order_log <- transform_sample_counts(physeq_order, function(x) log(x+1))
physeq_class_log <- transform_sample_counts(physeq_class, function(x) log(x+1))

# Standardize abundances to the median sequencing depth
total  <- median(sample_sums(physeq_fam_filt))
standf <- function(x, t=total) round(t * (x / sum(x)))
physeq_fam_filt_std <- transform_sample_counts(physeq_fam_filt, standf)

# Filter the taxa using a cutoff of 1 for the Coefficient of Variation
physeq_fam_filt_std2 <- filter_taxa(physeq_fam_log, function(x) sd(x)/mean(x) > 1.0, TRUE)


# bar plot of family by treatments------
bar_plot_phylum_by_treat <- 
  plot_bar(
  physeq_class_filt,
  x = "Phylum",
  fill = 'Phylum',
  y = 'Abundance',
  facet_grid = . ~ Treatment
) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack", color = 1) +
  scale_fill_viridis_d() +
  coord_flip()

ggsave(
  bar_plot_phylum_by_treat,
  filename = 'figures/bar_plot_phylum_by_treat.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 5,
  width = 6
)


# calculate ordination of log transformed sample counts
log_ord <- ordinate(physeq_fam_filt_std2, method = "RDA", distance = 'bray')

plot_ordination(
  physeq_fam_filt_std2,
  log_ord,
  type = "samples",
  color = "Treatment"
)

plot_ordination(
  physeq_fam_filt_std2,
  log_ord,
  type = "split",
  shape = "Treatment",
  color = 'Density',
  label = "Family",
  title = "PCO"
)

# CAP -------------
log_ord_CAP <-
  ordinate(
    physeq_log,
    "CAP",
    distance =  'bray',
    formula = ~ Treatment * Density
  )

cap_plot <-
  plot_ordination(physeq_fam,
                  log_ord_CAP,
                  color = 'Treatment') +
  geom_point(aes(size = Density), alpha = 0.5) +
  scale_radius(limits = c(0, 50),
               range = c(3, 7),
               name = 'Density') +
  scale_color_viridis_d(name = "Treatment", option = 'D')

cap_plot


plot_ordination(
  physeq_log,
  log_ord_CAP,
  type = "taxa",
  color = "Phylum",
  title = "taxa"
)


# PERMANOVA on phylosec objects------------------
metadata <- as(sample_data(physeq_fam_filt_log), "data.frame")

permanova_euk <-
  adonis(distance(physeq_fam_filt_log, method = "bray") ~ Treatment * final_density,
         data = metadata)
permanova_euk


euk_tax <- 
  read.csv('data/taxonomy_eukaryote.csv') %>% 
  rename(OTU = "X")

# Simper analysis--------------------
# extract abundance matrix from a phyloseq object
physeq_fam_OTUs  <-
  as(otu_table(physeq_fam_filt_log), "matrix") %>%
  t(.) %>% 
  as.data.frame()

# running the simper analysis on the dataframe and the variable of interest "time"
simper_euk_fam <-
  simper(physeq_fam_OTUs, metadata$Treatment, permutations = 100)

simper_euk_sabe_cont <-
  simper_euk_fam$Sabella_Control %>%
  data.frame() %>%
  rename(OTU = "species") %>%
  left_join(., euk_tax, by = "OTU") %>%
  filter(p < 0.05) %>%
  arrange(p) %>%
  dplyr::select(Family, ratio, ava, avb, cusum, p) %>%
  rename(Taxa = "Family") %>% 
  print(digits = 2)

simper_euk_mim_cont <-
  simper_euk_fam$Mimic_Control %>%
  data.frame() %>%
  rename(OTU = "species") %>%
  left_join(., euk_tax, by = "OTU") %>%
  filter(p < 0.05) %>%
  arrange(p)  %>%
  dplyr::select(Family, ratio, ava, avb, cusum, p) %>%
  rename(Taxa = "Family") %>% 
  print(digits = 2)

simper_euk <-
  bind_rows(Sabella_vs_Control = simper_euk_sabe_cont,
            Mimic_vs_Control = simper_euk_mim_cont,.id = "Groups")

# CAP biplot with phyloseq dat-----------------
physeq_fam_filt_log_CAP <-
  ordinate(
    physeq_fam_filt_log,
    "CAP",
    distance =  'bray',
    formula = ~ Treatment * final_density
  )

# CAP table
cap_anova_euk <- anova(physeq_fam_filt_log_CAP, by = 'terms')
cap_anova_euk


metadata <- as(sample_data(physeq_fam_filt_log), "data.frame")

cap_euk_dat <- fortify(physeq_fam_filt_log_CAP)

# eukaryote capscale biplot
cap_euk_sites <- 
  cap_euk_dat %>% 
  dplyr::filter(Score == 'sites') %>% 
  bind_cols(metadata)

cap_euk_taxa <-
  cap_euk_dat %>%
  dplyr::filter(Score == 'species') %>%
  dplyr::filter(abs(CAP1) > .15 | abs(CAP2) > .15) %>%
  rename(OTU = "Label") %>% 
  left_join(euk_tax, by = "OTU")


# CAP biplot
cap_euk_biplot <-
  ggplot(cap_euk_sites,
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
                        option = 'D') +
  geom_text_repel(
    data = cap_euk_taxa,
    aes(label = Family, x = CAP1 * 5, y =  CAP2 * 5),
    size = 3
  )
cap_euk_biplot


# without taxa---------
cap_euk_plot <- 
  ggplot(cap_euk_sites,
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


# Plot richness----
plot_richness(physeq_fam_filt_std2, x = "Treatment")

# heatmaps-----
x <- subset_taxa(physeq4, Class == "Polychaeta")
y  <- tax_glom(x, taxrank = "Family")
plot_heatmap(
  y,
  method = "NMDS",
  distance = "bray",
  taxa.label = "Family",
  taxa.order = "Family",
  sample.order = "Treatment",
  sample.label = "Treatment",
  low = "#000033",
  high = "#FF3300",
  trans = scales::log10_trans()
)


# relative abundance plots -----
physeq_rel <- transform_sample_counts(physeq4,function(x) x/sum(x))
dat_euk <- psmelt(physeq_rel)
dat_euk_raw <- psmelt(physeq4)

# Phylum
ggplot(dat_euk_raw %>% drop_na(Phylum),
       aes(
         x = Phylum,
         y = Abundance,
         fill = Phylum
       )) +
  stat_summary(fun.y = mean, geom = "bar", color = 'gray20') +  
  stat_summary(fun.data = mean_se, geom = "errorbar",width = 0.2, color = 'gray20')+
  scale_fill_viridis_d(guide = F) +
  coord_flip() +
  facet_wrap(~Treatment) +
  labs(y = 'Relative abundance (%)')

# or-
# or at Phylum level
phylum_relabund_euk <- 
  ggplot(dat_euk %>% drop_na(Phylum),
       aes(
         x = Phylum,
         y = Abundance*100,
         color = Treatment
       )) +
  stat_summary(fun.data = mean_se, size = .5,
               position = position_dodge(width = .5)) +
  scale_color_viridis_d() +
  labs(y = 'Relative abundance (%)', x='') +
  coord_flip() +
  theme(legend.position = c(.8,.8))

## Diversity indices---------------
euk_indices_raw <-
  dat_euk_raw %>%
  group_by(Sample) %>%
  mutate(
    N = sum(Abundance),
    H = diversity(Abundance),
    S = specnumber(Abundance),
    J = H / log(S)
  ) 

treat_euk <-
  sampledata %>% mutate(Sample = rownames(.))

euk_indices_raw <-
  euk_indices_raw %>%
  group_by(Sample) %>%
  summarise(S = mean(S),
            J = mean(J),
            N = mean(N)) %>%
  left_join(treat_euk, by = 'Sample') %>%
  write_csv(., 'data/cleaned_data/euk_indices_raw')

N_euk <- 
  ggplot(euk_indices_raw, aes(x = Treatment, y = N, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis_d(guide = F) +
  labs(x='', y = 'Number of individuals') + scale_y_log10()

S_euk <- 
  ggplot(euk_indices_raw, aes(x = Treatment, y = S, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis_d(guide = F) +
  labs(x='', y = 'Number of taxa')

J_euk <- 
  ggplot(euk_indices_raw, aes(x = Treatment, y = J, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis_d(guide = F) +
  labs(x='', y = 'Community evenness')

# ancova indices------
indices_euk_dat <-
  euk_indices_raw %>%
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

ancova_table_euk <- 
  indices_euk_dat %>% 
  select(index,ancova_table) %>% 
  unnest()

## model validation----
res_euk <- 
  indices_euk_dat %>% 
  select(index, residuals) %>% 
  unnest(residuals, .drop = TRUE)

# Plot of fitted vs. residual values by index to check the assumptions heteroscedasticity or any other pattern,
#  e.g. non-linearity or missing covariates
ggplot(res_euk) +
  geom_point(aes(x = .fitted, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2)

# boxplot by type
ggplot(res_euk) +
  geom_boxplot(aes(x = Treatment, y = .resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_hline(yintercept = 0,
             lty = 2,
             col = 2)

# qqplot of the normalised residuals to check the assumption of normality------------
ggplot(res_euk) +
  stat_qq(aes(sample = .std.resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_abline(
    intercept =  0,
    slope = 1,
    lty = 2,
    col = 2
  )
