# Analyses of eukaryote data from the field Sabella experiment 

# 01 load libraries--------------
library(phyloseq)
library(tidyverse)
library(ggvegan)
library(vegan)
library(ggrepel)
library(emmeans)
library(indicspecies)
source('theme_javier.R')
theme_set(theme_javier())

# 02 read feature table and taxonomic table--------------
otumat <- read.csv('data/eukaryote_data.csv', row.names = 1) %>% as.matrix(.)

# taxonomy data for phylosec
taxmat<- read.csv('data/taxonomy_eukaryote.csv',row.names = 1 ) %>% as.matrix(.)

# taxonomy as dataframe
# euk_tax <- 
#   read.csv('data/taxonomy_eukaryote.csv') %>% 
#   dplyr::rename(OTU = "X")

OTU  <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)

# 03 join and create phylosec object
physeq  <- phyloseq(OTU, TAX)

# read sample data-------
sampledata <- read.csv('data/molecular_treatments.csv', row.names = 1) %>% 
  sample_data(.)

# merge data ---------------
physeq1  <- merge_phyloseq(physeq, sampledata)

# 04 pre-process-------------
# remove ambient data, data with mean< 1e-4, unassigned  and  non-eukaryote taxa------
euk_data <- 
  subset_samples(physeq1, Treatment != "Ambient") %>% 
  subset_taxa(Domain=="Eukaryota") %>% 
  filter_taxa(function(x) mean(x) > 1e-5, TRUE) %>% 
  subset_taxa(Phylum!="")

table(tax_table(euk_data)[, "Phylum"], exclude = NULL)

# check data structure
ntaxa(euk_data)
nsamples(euk_data)
rank_names(euk_data)

# Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
# euk_data_filt <- filter_taxa(euk_data, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# # Transformation ----
# # filtered data
# physeq_filt_log <- transform_sample_counts(physeq_filt, function(x) log(x+1))
# physeq_fam_filt_log <- transform_sample_counts(physeq_fam_filt, function(x) log(x+1))

# # Filter the taxa using a cutoff of 1 for the Coefficient of Variation
# physeq_fam_filt_std2 <- filter_taxa(physeq_fam_log, function(x) sd(x)/mean(x) > 1.0, TRUE)

# rarerefy data------------
# get the total read abundace for each sample in descending order
samp_sum_euk <- 
  sample_sums(euk_data) %>% 
  data_frame() %>% 
  rename(N = '.') %>% 
  arrange(N) %>% 
  print(n  = 15)

rarecurve(
  t(otu_table(euk_data)),
  step = 50,
  cex = 0.5,
  label = F,
  main = 'Eukaryote samples',
  ylab = 'AVS'
)
abline(v = 5000, lty = 2)

euk_rare  <-
  rarefy_even_depth(euk_data,
                    sample.size = 5000,
                  rngseed = 123)
# 1 samples removed because they contained fewer reads than `sample.size`.
# SS.Bac.22c_S66	

ntaxa(euk_rare)

# taxonomy description-------
rank_names(euk_rare)
ntaxa(euk_rare) #1210
length(get_taxa_unique(euk_rare, "Phylum")) # 29
length(get_taxa_unique(euk_rare, "Class")) # 64
length(get_taxa_unique(euk_rare, "Order")) # 137
length(get_taxa_unique(euk_rare, "Family"))# 232
length(get_taxa_unique(euk_rare, "Genus")) # 304
length(get_taxa_unique(euk_rare, "Species")) #339


# physeq_euk_rare_filt  <- filter_taxa(physeq_euk_rare, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
euk_rare_log <- transform_sample_counts(euk_rare, function(x) log(x+1))
saveRDS(euk_rare_log, "data/output_data/euk_rare_log.RDS")

metadata_euk <- as(sample_data(euk_rare_log), "data.frame")

# PERMANOVA------------------
treat_euk <- as(sample_data(euk_rare_log), "data.frame")

permanova_euk <-
  adonis(distance(euk_rare_log, method = "bray") ~ Treatment +  final_density,
         data = treat_euk)
permanova_euk

euk_permanova_table <- 
  data.frame(permanova_euk$aov.tab) %>% 
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
    euk_rare_log, Treatment %in% c("Sabella", "Control")
  ), method = "bray") ~ Treatment ,
  data = filter(treat_euk, Treatment != "Mimic")
)

adonis(
  distance(subset_samples(
    euk_rare_log, Treatment %in% c("Mimic", "Control")
  ), method = "bray") ~ Treatment ,
  data = filter(treat_euk, Treatment != "Sabella")
)

adonis(
  distance(subset_samples(
    euk_rare_log, Treatment %in% c("Sabella", "Mimic")
  ), method = "bray") ~ Treatment ,
  data = filter(treat_euk, Treatment != "Control")
)

# SIMPER analysis--------------------
# extract abundance matrix from a phyloseq object
# agglomarate data----------------
euk_fam  <- tax_glom(euk_rare, "Family")
euk_phylum  <- tax_glom(euk_rare, "Phylum")

euk_rare_df  <-
  as(otu_table(euk_fam), "matrix") %>%
  t(.) %>% 
  as.data.frame()

# indicator spp--------
euk_rare_df  <-
  as(otu_table(euk_rare), "matrix") %>%
  t(.) %>% 
  as.data.frame()

ind_sp_euk <- multipatt(euk_rare_df, treat_euk$Treatment, duleg=TRUE) 

summary(ind_sp_euk, indvalcomp = TRUE)


# indicator taxa table--------------
taxo_table_euk <- 
  read.csv('data/taxonomy_eukaryote.csv') %>%
  rename(OTU = X)

ind_euk_table <- 
  ind_sp_euk$sign %>% 
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
  left_join(taxo_table_euk) %>% 
  write_csv('tables/ind_euk_table.csv')


ind_ASV_euk <- 
  ind_euk_table %>%
  select(OTU) %>% 
  as.list()
  
# relative abundance plot of indicator taxa -----
euk_rare_melt <- 
  euk_rare %>% 
  filter_taxa(function(x) mean(x) > 1, TRUE) %>% 
  # transform_sample_counts(.,function(x) (x/sum(x))*100) %>% 
  psmelt(.) 
  # write_csv('data/euk_rare_rel_df.csv')

ind_euk_plot_dat <- 
  euk_rare_melt %>%
  filter(OTU %in% ind_ASV_euk$OTU) %>% 
  mutate(Taxa = Family) %>% 
  write_csv('data/output_data/ind_euk_plot_dat.csv')

ind_euk_plot <-
  ggplot(ind_euk_plot_dat, aes(x = Family,
             y = Abundance+1,
             color = Treatment)) +
  stat_summary(fun.data = mean_se,
               size = .5,
               position = position_dodge(width = .5),
               alpha = 0.5) +
  labs(y = 'Abundance', x = '') +
  coord_flip() +
  scale_color_viridis_d() +
  theme(legend.position = c(.8, .2)) +
  scale_y_log10()

ind_euk_plot
# # SIMPER
# simper_euk_fam <- simper(euk_rare_df, treat_euk$Treatment, permutations = 100)
# 
# simper_euk <-
#   bind_rows("Sabella vs. Control " = simper_euk_fam$Sabella_Control %>% data.frame(),
#             "Mimic vs. Control" = simper_euk_fam$Mimic_Control %>% data.frame(),
#             "Mimic vcs. Sabella" = simper_euk_fam$Mimic_Sabella %>% data.frame(),
#             .id = "Groups") %>% 
#   dplyr::filter(cusum<.5) %>% 
#   dplyr::rename(OTU = "species") %>%
#   left_join(euk_tax, by = "OTU") %>%
#   dplyr::rename(Taxa = "Family") %>%
#   dplyr::select(Groups, Taxa, average, sd,ratio, ava, avb, cusum, p) %>%
#   print(digits = 2)


# CAP analysis with phyloseq dat-----------------
euk_CAP <-
  ordinate(
    euk_rare_log,
    "CAP",
    distance =  'bray',
    formula = ~ Treatment + final_density
  )

# CAP table
cap_anova_euk <- anova(euk_CAP, by = 'terms')
cap_anova_euk


euk_cap_table <- 
  data.frame(cap_anova_euk) %>% 
  rownames_to_column() %>% 
  dplyr::rename(Terms = rowname,
         df = Df,
         SS = SumOfSqs,
         P = Pr..F.)


cap_euk_dat <- fortify(euk_CAP)

# eukaryote capscale biplot----------
cap_euk_sites <- 
  cap_euk_dat %>% 
  dplyr::filter(Score == 'sites') %>% 
  bind_cols(metadata_euk) %>% 
  select(-Score,-Label) %>% 
  write_csv('data/output_data/cap_euk_sites.csv')


cap_euk_taxa <-
  cap_euk_dat %>%
  dplyr::filter(Score == 'species') %>%
  dplyr::filter(abs(CAP1) > .2 | abs(CAP2) > .2) %>%
  dplyr::rename(OTU = "Label") %>% 
  left_join(euk_tax, by = "OTU")


# CAP biplot
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


cap_euk_biplot <-
  cap_euk_plot   +
  geom_text_repel(
    data = cap_euk_taxa,
    aes(label = Family, x = CAP1 * 5, y =  CAP2 * 5),
    size = 4
  )
cap_euk_biplot



# relative abundance plots -----
euk_rare_rel_df <- 
  euk_rare %>% 
  filter_taxa(function(x) mean(x) > 1, TRUE) %>% 
  transform_sample_counts(.,function(x) (x/sum(x))*100) %>% 
  psmelt(.) %>% 
  write_csv('data/output_data/euk_rel_abund.csv')

# # or at Phylum level
# phylum_relabund_euk <-
#   ggplot(euk_rare_rel_df %>% drop_na(Phylum),
#          aes(x = Phylum,
#              y = Abundance,
#              color = Treatment)) +
#   stat_summary(fun.data = mean_se,
#                size = .5,
#                position = position_dodge(width = .5)) +
#   scale_color_viridis_d() +
#   labs(y = 'Relative abundance (%)', x = '') +
#   coord_flip() +
#   theme(legend.position = c(.8, .8))

## Diversity indices---------------
dat_euk_rare_raw <- 
  euk_rare %>% 
  psmelt(.)


euk_indices <-
  dat_euk_rare_raw %>%
  group_by(Sample) %>%
  mutate(
    # N = sum(Abundance),
    H = diversity(Abundance),
    S = specnumber(Abundance),
    J = H / log(S)
  ) 

treat_euk_df <- treat_euk %>% mutate(Sample = rownames(.))

euk_indices_raw <-
  euk_indices %>%
  group_by(Sample) %>%
  summarise(S = mean(S),
            J = mean(J)) %>%
  left_join(treat_euk_df, by = 'Sample') %>%
  write_csv('data/output_data/euk_indices_raw.csv')

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
  gather(index, value, c("S")) %>%
  group_by(index) %>%
  nest() %>%
  mutate(
    ancova = map(.x = data, ~ aov(value ~ Treatment + final_density, data = .x)),
    anova_tab = map(ancova, ~car::Anova(., type = 2)),
    ancova_table = map(anova_tab, broom::tidy),
    residuals = map(ancova, broom::augment),
    emm = map(ancova,~emmeans(., ~Treatment)),
    posthoc = map(emm,pairs))

ancova_table_euk <-
  indices_euk_dat %>%
  select(index, ancova_table) %>%
  unnest()

indices_euk_dat %>%
  select(index, posthoc) %>%
  flatten()

## model validation----
res_euk <- 
  indices_euk_dat %>% 
  select(index, residuals) %>% 
  unnest(residuals, .drop = TRUE)

# Plot of fitted vs. residual values by index to check the assumptions
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

# qqplot of the normalised residuals to check normality
ggplot(res_euk) +
  stat_qq(aes(sample = .std.resid), alpha = .3) +
  facet_wrap( ~ index, scale = 'free') +
  geom_abline(
    intercept =  0,
    slope = 1,
    lty = 2,
    col = 2
  )
