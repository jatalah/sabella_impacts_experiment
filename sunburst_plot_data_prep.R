library(plotly)
library(phyloseq)
library(tidyverse)
library(ggpubr)

# bacteria ---
bact_rare_rel_df <- read_csv('data/bact_rel_abund.csv')

data_bact_sunburst <-
  bact_rare_rel_df %>% 
  group_by(Phylum, Class,  Sub_class, Family) %>%
  summarise(N = sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(prop = N/sum(N)*100) %>% 
  select(N, everything()) %>% 
  arrange(desc(prop)) %>% 
  print()


data_bact_sunburst <- 
  data_bact_sunburst %>% 
  group_by(Phylum) %>%
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  mutate(prop = N/sum(N)*100) %>% 
  filter(prop>.5) %>% 
  arrange(desc(prop)) %>% 
  write_csv('data/output_data/data_bact_sunburst.csv')

data_bact_sunburst <- read_csv('data/output_data/data_bact_sunburst.csv')
bact_donut <- 
  ggplot(data_bact_sunburst, aes(x = 2, fill = Phylum, y = prop)) +
  geom_bar(stat = 'identity', color = "black") +
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  # scale_fill_brewer(palette = "Spectral") +
  scale_fill_viridis_d(option = 'B') +
  xlim(0.5, 2.5) +
  labs(subtitle = 'Bacteria')
bact_donut

# eukaryote---
data_euk_sunburst <-
  euk_rare_rel_df %>%
  group_by(Phylum) %>%
  summarise(N = sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(prop = N/sum(N)*100) %>% 
  select(N, everything())

data_euk_sunburst <-
  data_euk_sunburst %>%
  arrange(desc(Phylum)) %>%
  arrange(desc(prop)) %>%
  write_csv('data/output_data/data_euk_sunburst.csv')

# macrofauna -------------
data_mac_sunburst <-
  psmelt(physeq_mac) %>%
  group_by(Phylum) %>%
  summarise(N = sum(Abundance)) %>% 
  ungroup() %>% 
  mutate(prop = N/sum(N)*100) %>% 
  select(N, everything()) %>% 
  arrange(desc(prop)) %>% 
  write_csv('data/output_data/data_mac_sunburst.csv')

all_sunburst_data <- 
  bind_rows(Macrofauna = data_mac_sunburst, Eukaryote = data_euk_sunburst,.id = 'dataset') %>% 
  group_by(dataset) %>% 
  arrange(desc(Phylum)) %>%
  mutate(pos = cumsum(prop)- prop/2) %>% 
  filter(prop>.2) %>% 
  write_csv('data/output_data/all_sunburst_data.csv')

all_sunburst_data <-  read_csv('data/output_data/all_sunburst_data.csv')
euk_mac_donut <- 
  ggplot(all_sunburst_data, aes(x = 2, fill = Phylum, y = prop)) +
  geom_bar(stat = 'identity', color = "black") +
  coord_polar(theta = "y", start = 0) +
  facet_wrap(~dataset, ncol = 1) +
  theme_void() +
  scale_fill_viridis_d() +
  xlim(0.5, 2.5) 

# save all plot ----
all_donuts <- 
  ggarrange(euk_mac_donut, bact_donut)

ggsave(
  all_donuts,
  filename = 'figures/all_dounuts.tiff',
  device = 'tiff',
  compression = 'lzw',width =6, height = 6 
)

ggsave(
  all_donuts,
  filename = 'figures/all_dounuts.svg',
  device = 'svg',
  width =6, height = 6 
)