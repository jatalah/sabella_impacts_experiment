# load libraries -------------
library(tidyverse)
library(ggpubr)
library(ggrepel)

source('H:/javiera/Stats/R folder/theme_javier.R')
theme_set(theme_javier())

# Figure 1 Sunburst diagrams --------------------------
all_sunburst_data <-  read_csv('data/output_data/all_sunburst_data.csv')

euk_mac_donut <-
  all_sunburst_data %>%
  filter(dataset != "c. Bacteria") %>%
  ggplot(aes(x = 2, fill = Phylum, y = prop)) +
  geom_bar(stat = 'identity', color = "gray40") +
  coord_polar(theta = "y", start = 0) +
  theme_void(base_size = 10) +
  scale_fill_viridis_d(option = 'D', alpha = .8) +
  xlim(0.5, 2.5) +
  facet_wrap( ~ dataset) +
  theme(
    legend.position = c(.35,-.05),
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
  theme_void(base_size = 10) +
  scale_fill_viridis_d(option = 'B', alpha = .8) +
  xlim(0.5, 2.5) +
  facet_wrap( ~ dataset) +
  theme(
    legend.position = c(.35,-.05),
    legend.title = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "lines"),
    legend.spacing.y = unit(0, 'line'),
    legend.spacing.x = unit(.1, 'line'),
    legend.key.size = unit(.3, "cm")
  ) +
  guides(fill = guide_legend(nrow = 5,
                             byrow = TRUE))


all_donuts <- ggarrange(euk_mac_donut, bact_donut, widths = c(2, 1))

ggsave(
  all_donuts,
  filename = 'figures/Figure 1.tiff',
  device = 'tiff',
  compression = 'lzw',
  width = 18,
  height = 9,
  units = 'cm',
  dpi = 300
)


# Figure 2 CAP ordination bubble plots -------------------
cap_bact_sites <- read_csv('data/output_data/cap_bact_sites.csv')
cap_euk_sites <- read_csv('data/output_data/cap_euk_sites.csv')
cap_mac_sites <- read_csv('data/output_data/cap_mac_sites.csv')

all_cap_data <-
  bind_rows(
    `a. Macrofauna` = cap_mac_sites,
   `b. Eukaryote` = cap_euk_sites,
    `c. Bacteria` = cap_bact_sites,
    .id = 'dataset'
  )

CAP_biplots <- 
  ggplot(all_cap_data,
         aes(
           x = CAP1,
           y =  CAP2
         )) +
  geom_point(alpha = .5, aes(color = Treatment, 
                             size = Density)) + 
  scale_radius(limits = c(0, 50),
               range = c(3, 7),
               name = 'Density') +
  scale_color_viridis_d(name = "Treatment",
                        option = 'D') +
  facet_wrap(~dataset, scales = 'free') +
  theme(legend.position = 'bottom')
  
ggsave(
  CAP_biplots,
  filename = 'figures/Figure 2.tiff',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 80,
  width = 169,
  units = 'mm'
)

# Figure 3 Abundace plots of indicator taxa--------------------
ind_bact_plot_dat <-read_csv('data/output_data/ind_bact_plot_data.csv')
ind_euk_plot_dat <- read_csv('data/output_data/ind_euk_plot_dat.csv')
ind_mac_plot_dat <- 
  read_csv('data/output_data/ind_mac_plot_dat.csv') %>% 
  mutate(Taxa = if_else(
    condition = str_detect(Taxa, " "),
    true = paste0('italic(', Taxa, ')'),
    false = as.character(Taxa)
  )) %>%
  mutate(Taxa = str_replace(Taxa," ","~"))

ind_taxa_all_data <-
  bind_rows(
    `a. Macrofauna` = ind_mac_plot_dat,
    `b. Eukaryote` = ind_euk_plot_dat,
    `c. Bacteria` = ind_bact_plot_dat,
    .id = 'dataset'
  )

phylum_relabund_plots <- 
  ggplot(ind_taxa_all_data,
         aes(x = Taxa,
             y = Abundance + 1,
             color = Treatment)) +
  stat_summary(
    fun.data = mean_se,
    size = .5,
    position = position_dodge(width = .5),
    alpha = 0.5
  ) +
  labs(y = '', x = '') +
  coord_flip() +
  scale_color_viridis_d() +
  theme_javier(base_size = 10) +
  facet_wrap( ~ dataset, scales = 'free') +
  theme(legend.position = 'bottom') +
  scale_y_log10() + 
  scale_x_discrete(labels = ggplot2:::parse_safe)

ggsave(
  phylum_relabund_plots,
  filename = 'figures/Figure 3.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 120,
  width = 189,
  units = 'mm'
)

# Figure 4 Boxplots of taxa richness ----------------------
indices_mac <- read_csv('data/output_data/macro_indices.csv')
euk_indices <- read_csv('data/output_data/euk_indices_raw.csv')
bact_indices <- read_csv('data/output_data/bact_indices_raw.csv')

indices_all_data <- 
  bind_rows(
    `a. Macrofauna` = indices_mac,
    `b. Eukaryote` = euk_indices,
    `c. Bacteria` = bact_indices,
    .id = 'dataset'
  )

S_plots <- 
  ggplot(indices_all_data, aes(x = Treatment, y = S, fill = Treatment)) +
  geom_boxplot(alpha = .7) +
  scale_fill_viridis_d(guide = F) +
  labs(x='', y = 'No. of taxa')  +
  theme_javier(base_size = 10) +
  facet_wrap(~dataset, scales = 'free') +
  theme(legend.position = 'bottom')

S_plots

ggsave(
  S_plots,
  filename = 'figures/Figure 4.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 65,
  width = 169, 
  units = 'mm',
  dpi = 600
)


# Figure 5 boxplots of env variables------------
env_long <-
  read_csv('data/output_data/env_long.csv') %>%
  mutate(
    key = fct_relevel(
      key,
      "Clay",
      "Silt",
      "Very fine sand",
      'Fine sand',
      "Medium sand",
      "Coarse sand",
      "Gravel"
    )
  )

env_boxplots <- 
  ggplot(env_long,
       aes(x = Treatment, y = value, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_viridis_d(guide = F, alpha = .7) +
  theme_javier(base_size = 10) +
  theme(axis.title = element_blank()) +
  facet_wrap( ~ key, scales = 'free', nrow = 2)

ggsave(
  env_boxplots,
  filename = 'figures/env_boxplots.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 65,
  width = 169, 
  units = 'mm',
  dpi = 600
)

# Figure 6 dbRDA biplots ----------------
# this biplots come from "scripts/dbRDA_distlm.R"



rda_all_plots <-
  ggarrange(
    macro_rda_plot,
    euk_rda_plot,
    bact_rda_plot,
    common.legend = T,
    nrow = 1,
    legend = 'bottom',
    align = 'h'
  )

ggsave(
  rda_all_plots,
  filename = 'figures/Figure 5.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 65,
  width = 169, 
  units = 'mm',
  dpi = 600
)
