# load libraries -------------
library(tidyverse)
library(ggpubr)
library(ggrepel)

source('H:/javiera/Stats/R folder/theme_javier.R')
theme_set(theme_javier())

# 01 Figure 1 Sunburst diagrams --------------------------
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
all_donuts <- ggarrange(euk_mac_donut, bact_donut)

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

# Figure 3 Abundace plots of indicator taxa----------------------
ind_bact_plot_dat <-read_csv('data/output_data/ind_bact_plot_data.csv')
ind_euk_plot_dat <- read_csv('data/output_data/ind_euk_plot_dat.csv')
ind_mac_plot_dat <- read_csv('data/output_data/ind_mac_plot_dat.csv')

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
  scale_y_log10()


ggsave(
  phylum_relabund_plots,
  filename = 'figures/Figure 3.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 100,
  width = 169,
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

# Figure 5 dbRDA biplots ----------------
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
