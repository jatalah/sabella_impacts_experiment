library(plotly)
library(phyloseq)
library(tidyverse)
library(ggpubr)

# bacteria ---
bact_rare_rel_df <- read_csv('data/bact_rel_abund.csv')

data_bact_sunburst <-
  bact_rare_rel_df %>%
  group_by(Phylum) %>%
  summarise(N = sum(Abundance)) %>%
  ungroup() %>%
  mutate(prop = N / sum(N) * 100) %>%
  arrange(desc(prop))

# eukaryote---
data_euk_sunburst <-
  euk_rare_rel_df %>%
  group_by(Phylum) %>%
  summarise(N = sum(Abundance)) %>%
  ungroup() %>%
  mutate(prop = N / sum(N) * 100) %>%
  select(N, everything()) %>%
  arrange(desc(Phylum)) %>%
  arrange(desc(prop))

# macrofauna -------------
data_mac_sunburst <-
  psmelt(physeq_mac) %>%
  group_by(Phylum) %>%
  summarise(N = sum(Abundance)) %>%
  ungroup() %>%
  mutate(prop = N / sum(N) * 100) %>%
  arrange(desc(prop))

all_sunburst_data <-
  bind_rows(
    `a. Macrofauna` = data_mac_sunburst,
    `b. Eukaryote` = data_euk_sunburst,
    `c. Bacteria` = data_bact_sunburst,
    .id = 'dataset'
  ) %>%
  group_by(dataset) %>%
  arrange(desc(Phylum)) %>%
  mutate(pos = cumsum(prop) - prop / 2) %>%
  write_csv('data/output_data/all_sunburst_data.csv')

euk_mac_donut <-
  all_sunburst_data %>%
  filter(dataset != "c. Bacteria") %>%
  ggplot(aes(x = 2, fill = Phylum, y = prop)) +
  geom_bar(stat = 'identity', color = "gray40") +
  coord_polar(theta = "y", start = 0) +
  theme_void(base_size = 9) +
  scale_fill_viridis_d(option = 'D', alpha = .8) +
  xlim(0.5, 2.5) +
  facet_wrap( ~ dataset) +
  theme(
    legend.position = c(.5,-.1),
    legend.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "lines"),
    legend.spacing.y = unit(0, 'line'),
    legend.spacing.x = unit(.1, 'line'),
    legend.key.size = unit(.3, "cm")
  ) +
  guides(fill = guide_legend(nrow = 5,
                             byrow = TRUE))

bact_donut <-
  all_sunburst_data %>%
  filter(dataset == "c. Bacteria" & prop > .2) %>%
  ggplot(aes(x = 2, fill = Phylum, y = prop)) +
  geom_bar(stat = 'identity', color = "gray40") +
  coord_polar(theta = "y", start = 0) +
  theme_void(base_size = 9) +
  scale_fill_viridis_d(option = 'B', alpha = .8) +
  xlim(0.5, 2.5) +
  facet_wrap( ~ dataset) +
  theme(
    legend.position = c(.4,-.1),
    legend.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "lines"),
    legend.spacing.y = unit(0, 'line'),
    legend.spacing.x = unit(.1, 'line'),
    legend.key.size = unit(.3, "cm")
  ) +
  guides(fill = guide_legend(nrow = 5,
                             byrow = TRUE))


# save all plot ----
all_donuts <- ggarrange(euk_mac_donut, bact_donut, widths = c(2, 1))


