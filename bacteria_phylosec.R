library(tidyverse)
library(phyloseq)

# read OTU table with OTU > 5 reads--------------
otu_bact <-
  read.csv("data/bacteria_data.csv",
           row.names = 1,
           check.names = FALSE) %>% 
  select(-"SS.Ext.Blank_S94", - "SS.PCRBlank_S95" , -"SS.water_S96" ) %>% 
  as.matrix(.)

tax_bact<- read.csv('data/taxonomy_bacteria.csv',row.names = 1 )%>% as.matrix(.)

sampledata <- 
  sample_data(read.csv('data/molecular_treatments.csv', row.names = 1))


OTU_bact  <- otu_table(otu_bact, taxa_are_rows = TRUE)
TAX_bact <- tax_table(tax_bact)

# merge data ---------
physeq_bact <- phyloseq(OTU_bact, TAX_bact, sampledata)

# pre-process-------------
physeq_bact1 <- subset_samples(physeq_bact, Treatment != "Ambient")
physeq_bact2 <- filter_taxa(physeq_bact1, function(x) mean(x) > 1e-5, TRUE)

# rarerefy data
physeq_bact_rare  <- rarefy_even_depth(physeq_bact2)
ntaxa(physeq_bact_rare)
physeq_bact_rare_fam  <- tax_glom(physeq_bact_rare, "Family")
physeq_bact_rare_fam_filt  <- filter_taxa(physeq_bact_rare_fam, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_bact_rare_fam_filt_log <- transform_sample_counts(physeq_bact_rare_fam_filt, function(x) log(x+1))

# check data names and taxa
ntaxa(physeq_bact2)
nsamples(physeq_bact2)
rank_names(physeq_bact2)

# agglomarate at family level------------
physeq_bact_fam  <- tax_glom(physeq_bact2, "Family")
physeq_bact_subc  <- tax_glom(physeq_bact2, "Sub_class")
physeq_bact_class  <- tax_glom(physeq_bact2, "Class")
physeq_bact_phyl  <- tax_glom(physeq_bact2, "Phylum")

# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
physeq_bact_fam_filt  <- filter_taxa(physeq_bact_fam, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_bact_subc_filt  <- filter_taxa(physeq_bact_subc, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_bact_class_filt  <- filter_taxa(physeq_bact_class, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_bact_phyl_filt  <- filter_taxa(physeq_bact_phyl, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
physeq_bact_filt  <- filter_taxa(physeq_bact2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# log transform ----
# filtered data
physeq_bact_phyl_filt_log <- transform_sample_counts(physeq_bact_phyl_filt, function(x) log(x+1))
physeq_bact_class_filt_log <- transform_sample_counts(physeq_bact_class_filt, function(x) log(x+1))
physeq_bact_fam_filt_log <- transform_sample_counts(physeq_bact_fam_filt, function(x) log(x+1))
physeq_bact_subc_filt_log <- transform_sample_counts(physeq_bact_subc_filt, function(x) log(x+1))
physeq_bact_filt_log <- transform_sample_counts(physeq_bact_class_filt, function(x) log(x+1))

# unfiltered data-------------
physeq_bact_log <- transform_sample_counts(physeq_bact2, function(x) log(x+1))
physeq_bact_fam_log <- transform_sample_counts(physeq_bact_fam, function(x) log(x+1))
physeq_bact_subc_log <- transform_sample_counts(physeq_bact_subc, function(x) log(x+1))
physeq_bact_class_log <- transform_sample_counts(physeq_bact_class, function(x) log(x+1))



# bar plot of family by treatments------
bar_plot_all_bact<- 
  plot_bar(
    physeq_bact_fam_filt_log,
    x = "Phylum",
    fill = 'Phylum',
    y = 'Abundance',
    facet_grid = . ~ Treatment
  ) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack") +
  scale_fill_viridis_d(guide = F) +
  coord_flip()

ggsave(
  bar_plot_all_bact,
  filename = 'figures/bar_plot_all_bact.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 5,
  width = 10
)

# bar plot proteo------
physeq_proteo  <-  subset_taxa(physeq_bact_fam_filt_log, Phylum=="Proteobacteria")
bar_plot_proteo <- 
  plot_bar(
    physeq_proteo,
    x = "Sub_class",
    fill = 'Sub_class',
    y = 'Abundance',
    facet_grid = . ~ Treatment
  ) +
  geom_bar(aes(fill = Sub_class), stat = "identity", position = "stack") +
  scale_fill_viridis_d(guide = F) +
  coord_flip()

ggsave(
  bar_plot_proteo,
  filename = 'figures/bar_plot_proteo.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 5,
  width = 10
)

physeq_bacteroidetes  <-  subset_taxa(physeq_bact_fam_filt_log, Phylum=="Bacteroidetes")
bar_plot_bacteroidetes <- 
  plot_bar(
    physeq_bacteroidetes,
    x = "Sub_class",
    fill = 'Sub_class',
    y = 'Abundance',
    facet_grid = . ~ Treatment
  ) +
  geom_bar(aes(fill = Sub_class), stat = "identity", position = "stack") +
  scale_fill_viridis_d(guide = F) +
  coord_flip()

ggsave(
  bar_plot_bacteroidetes,
  filename = 'figures/bar_plot_bacteroidetes.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 5,
  width = 10
)


# relative abundance plots -----
physeq_bact_rel <- transform_sample_counts(physeq_bact_filt,function(x) x/sum(x))
dat_bact_rel <- psmelt(physeq_bact_rel)
dat_bact <- psmelt(physeq_bact_filt)

# Phylum level------------------
ggplot(dat_bact %>% filter(Phylum != 0),
         aes(
           x = Phylum,
           y = Abundance,
           color = Treatment
         )) +
  stat_summary(fun.data = mean_se, size = .5,
               position = position_dodge(width = .5)) +
  scale_color_viridis_d() +
  labs(y = 'Abundance') +
  coord_flip() +
  theme(legend.position = c(.8,.8))



# PERMANOVA on phylosec objects------------------
library(vegan)
taxono_bact <- 
  read_csv('data/taxonomy_bacteria.csv') %>%
  rename(OTU = "#OTU ID")

metadata_bact <- as(sample_data(physeq_bact_class_filt_log), "data.frame")

permanova_bact <-
  adonis(distance(physeq_bact_class_filt_log, method = "bray") ~ Treatment * final_density,
         data = metadata_bact, permutations = 999)
permanova_bact


# SIMPER--------------------
## Simper analysis, extract abundance matrix from a phyloseq object
physeq_bact_fam_OTUs  <-
  as(otu_table(physeq_bact_class_filt_log), "matrix") %>%
  t(.) %>% 
  as.data.frame()

# running the simper analysis on the dataframe and the variable of interest "time"
simper_bact_fam <- simper(physeq_bact_fam_OTUs, metadata_bact$Treatment, permutations = 100)

simper_bact_cont_sab <- 
  simper_bact_fam$Sabella_Control %>%
  data.frame() %>%
  rename(OTU = "species") %>%
  left_join(., taxono_bact, by = "OTU") %>%
  filter(p < 0.05) %>%
  # arrange(p) %>%
  dplyr::select(Class, ratio, ava, avb, cusum, p) %>%
  rename(Taxa = "Class") %>% 
  print(digits = 2)

simper_bact_mim_cont <- 
simper_bact_fam$Mimic_Control %>%
  data.frame() %>% 
  rename(OTU = "species") %>% 
  left_join(., taxono_bact, by = "OTU") %>% 
  filter(p<0.05) %>% 
  arrange(p)  %>% 
  dplyr::select(Class, ratio, ava, avb,cusum,p) %>% 
  rename(Taxa = "Class") %>% 
  print(digits = 2)

simper_bact <-
  bind_rows("Sabella vs. Control" = simper_bact_cont_sab,
            "Mimic vs. Control" = simper_bact_mim_cont,.id = "Groups")

# CAP analysis-----------------
physeq_bact_fam_filt_log_CAP <-
  ordinate(
    physeq_bact_fam_filt_log,
    "CAP",
    distance =  'bray',
    formula = ~ Treatment * final_density
  )


cap_anova_bact <- anova(physeq_bact_fam_filt_log_CAP, by = 'terms', permutations = 9999)
cap_anova_bact

metadata_bact <- as(sample_data(physeq_bact_class_filt_log), "data.frame")

cap_bact_dat <- fortify(physeq_bact_class_filt_log_CAP)

# bacteria capscale biplot------
cap_bact_sites <- 
  cap_bact_dat %>% 
  dplyr::filter(Score == 'sites') %>% 
  bind_cols(metadata_bact)

cap_bact_taxa <-
  cap_bact_dat %>%
  dplyr::filter(Score == 'species') %>%
  dplyr::filter(abs(CAP1) > .1 | abs(CAP2) > .1) %>%
  rename(OTU = "Label") %>% 
  left_join(taxono_bact, by = "OTU")

# CAP biplot
cap_bact_biplot <-
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
                        option = 'D') +
  geom_text_repel(
    data = cap_bact_taxa,
    aes(label = Class, x = CAP1*10 , y =  CAP2*10 ),
    size = 3
  )
cap_bact_biplot

# Diversity indices bacteria----------
dat_bact_raw <- psmelt(physeq_bact2) 
# dat_bact_raw <- psmelt(physeq_bact_fam)

bact_indices_raw <-
  dat_bact_raw %>%
  group_by(Sample) %>%
  mutate(
    N = sum(Abundance),
    H = diversity(Abundance),
    S = specnumber(Abundance),
    J = H / log(S)
  )

treat_bact <- sampledata %>% mutate(Sample = rownames(.))
bact_indices_raw <- 
  bact_indices_raw %>%
  group_by(Sample) %>%
  summarise(S = mean(S),
            J = mean(J),
            N = mean(N)) %>%
  left_join(treat_bact, by = 'Sample') %>% 
  write_csv(., 'data/bact_indices_raw')

N_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = N, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'Number of individuals') + scale_y_log10()

S_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = S, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'Number of taxa')

J_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = J, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'Community evenness')


# anova indices---
# ancova indices------
indices_bact_dat <-
  bact_indices_raw %>%
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

ancova_table_bact <- 
  indices_bact_dat %>% 
  select(index,ancova_table) %>% 
  unnest()

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
