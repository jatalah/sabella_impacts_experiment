# load libraries------------
library(phyloseq)
library(tidyverse)
rm(list = ls())
source('other scripts/theme_javier.R')
theme_set(theme_javier())

# data preparation------------
phylosec_mac <- 
  read_csv('data/fauna_dat.csv') %>% ### from macrofauna.R script
  dplyr::select(-Treatment, -Density) %>% 
  mutate(sample = paste0("sa",sample)) %>% 
  gather(Taxa, abund, -sample) %>%
  spread(sample, abund) %>% 
  column_to_rownames("Taxa") %>% 
  as.matrix(.)


taxmat_mac <- 
  read.csv('data/taxonomy_mac.csv',row.names = 1 ) %>% 
  dplyr::select(-query, -Genus, -Species) %>% 
  as.matrix(.)

OTU_mac <- otu_table(phylosec_mac, taxa_are_rows = TRUE)
TAX_mac <- tax_table(taxmat_mac)

sampledat_mac <-
  sample_data(
    read.csv('data/fauna_dat.csv', row.names = 1) %>%
      dplyr::select(Treatment, Density)
  ) %>% 
  arrange(row.names(.)) %>% 
  as.data.frame()

sampledata_mac <- sample_data(sampledat_mac)
physeq_mac  <- phyloseq(OTU_mac, TAX_mac)

# merge data------
physeq_mac1  <- merge_phyloseq(physeq_mac, sampledata_mac)

# average samples by treatment---------
physeq_mac1_treat_mean <- merge_samples(physeq_mac1, "Treatment", fun=mean)
physeq_mac1_rel <- transform_sample_counts(physeq_mac1,function(x) x/sum(x))
physeq_mac1_log <- transform_sample_counts(physeq_mac1,function(x) log10(x +1))

dat <- psmelt(physeq_mac1) # create data frame
dat_rel <- psmelt(physeq_mac1_rel)

# Figure xx macrofauna barplots by group ------
# class-----
ggplot(dat %>% drop_na(Class),
       aes(
         x = Class,
         y = Abundance,
         fill = Class
       )) +
  stat_summary(fun.y = mean, geom = "bar", color = 'gray20') +  
  stat_summary(fun.data = mean_se, geom = "errorbar",width = 0.2, color = 'gray20')+
  scale_fill_viridis_d(guide = F) +
  coord_flip() +
  facet_wrap(~Treatment) +
  labs(y = 'Relative abundance')

# or ----
# ggplot(dat %>% drop_na(Class),
#        aes(
#          x = Class,
#          y = Abundance,
#          fill = Treatment,
#          group =  Treatment
#        )) +
#   stat_summary(
#     fun.y = mean,
#     geom = "bar",
#     color = 'gray20',
#     position = position_dodge(width = .9)
#   ) +
#   stat_summary(
#     fun.data = mean_se,
#     geom = "errorbar",
#     width = 0.2,
#     color = 'gray20',
#     position = position_dodge(width = .9)
#   ) +
#   scale_fill_viridis_d() +
#     labs(y = 'Relative abundance')

# or at Phylum level------------------
ggplot(dat %>% drop_na(Phylum),
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


# family ----
ggplot(dat,
       aes(
         x = Family,
         y = Abundance,
         color = Treatment
       )) +
  stat_summary(fun.data = mean_se, size = .5,
               position = position_dodge(width = .5)) +
  scale_color_viridis_d() +
  labs(y = 'Abundance') +
  coord_flip() +
  theme(legend.position = c(.8,.8))


plot_bar(
  physeq_mac1_rel,
  x = "Phylum",
  fill = 'Phylum',
  y = 'Abundance',
  facet_grid = .~Treatment
) +
  geom_bar(aes(fill = Phylum), stat = "identity") +
  scale_fill_viridis_d(guide =F) +
  coord_flip()


plot_bar(
  physeq_mac1,
  x = "Class",
  fill = 'Class',
  y = 'Abundance',
  facet_grid = .~Treatment
) +
  geom_bar(aes(fill = Class), stat = "identity", position = "stack") +
  scale_fill_viridis_d(guide =F) +
  coord_flip()


physeq_mac1_polyc <- subset_taxa(physeq_mac1, Class=="Polychaeta")

bar_plot_ploychaeta <- 
  plot_bar(
  physeq_mac1_polyc,
  x = "Family",
  fill = 'Family',
  y = 'Abundance',
  facet_grid = .~Treatment
) +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack") +
  scale_fill_viridis_d(guide =F) +
  coord_flip()

ggsave(
  bar_plot_ploychaeta,
  filename = 'figures/bar_plot_ploychaeta.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 5,
  width = 10
)

physeq_mac1_bivalvia <- subset_taxa(physeq_mac1, Class=="Bivalvia")

bar_plot_bivalvia <- 
  plot_bar(
    physeq_mac1_bivalvia,
    x = "Family",
    fill = 'Family',
    y = 'Abundance',
    facet_grid = .~Treatment
  ) +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack") +
  scale_fill_viridis_d(guide =F) +
  coord_flip()
bar_plot_bivalvia


plot_bar(
  physeq_mac1_polyc,
  x = "Family",
  fill = 'Family',
  y = 'Abundance',
  facet_grid = .~Treatment
) +geom_bar(aes(fill = Family), stat = "identity", position = "stack") +
  scale_fill_viridis_d(guide =F) +
  coord_flip() +
  geom_boxplot()

library(scales)
plot_heatmap(
  physeq_mac1,
  "NMDS",
  "bray",
  taxa.order = "Class",
  sample.order = "Treatment",
  sample.label = "Treatment",
  # na.value="black",
  # low="#FFFFCC", high="#FF3300",
  low="#000033", high="#FF3300"
)

## rel abund calcs 
dat %>%
  group_by(Sample) %>%
  summarise(n = sum(Abundance)) %>%
  right_join(., dat, by = 'Sample') %>%
  group_by(Phylum, Treatment, Sample) %>%
  mutate(rel = Abundance / n * 100) %>% 
  group_by(Phylum, Treatment) %>% 
  summarise(mean_rel = mean(rel)) %>% 
  print(n = 30)


