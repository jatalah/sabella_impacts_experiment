# load libraries -------------
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(viridis)
source('H:/javiera/Stats/R folder/theme_javier.R')
theme_set(theme_javier())

# 01 Figure ordinations for PCO and CAP--------------------------
# pcos_caps <- 
# ggarrange(
#   pco_ordination,
#   cap_plot,
#   pco_euk_ordination,
#   cap_euk_plot,
#   pco_plot_bact, 
#   cap_plot_bact,
#   common.legend = T,
#   labels = "auto",
#   ncol = 2,
#   nrow = 3,
#   legend = 'bottom' 
# )
# 
# ggsave(
#   pcos_caps,
#   filename = 'figures/all_PCO_CAP.tif',
#   device = 'tiff',
#   compression = 'lzw',
#   height = 8,
#   width = 5
# )


## CAP plots for macrofauna.R , aukaryotes.R and bacteria.R------
CAP_biplots <- 
  ggarrange(
    cap_biplot + ggtitle('Macrofauna'),
    cap_euk_biplot + ggtitle('Eukaryotes'),
    cap_bact_biplot + ggtitle('Bacteria'),
    common.legend = T,
    # labels = "auto",
    ncol = 3,
    nrow = 1,
    legend = 'bottom'
  )
CAP_biplots

ggsave(
  CAP_biplots,
  filename = 'figures/all_CAP.tif',
  device = 'tiff',
  compression = 'lzw',
  dpi = 600,
  height = 6*.9,
  width = 16*.9
)

## Boxplots diversity
indices_mac <- read_csv('data/macro_indices.csv')
 
S_mac <- 
  ggplot(indices_mac, aes(x = Treatment, y = S, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'S') 
N_mac <- 
  ggplot(indices_mac, aes(x = Treatment, y = N, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'N') + scale_y_log10()
  
J_mac <- 
  ggplot(indices_mac, aes(x = Treatment, y = J, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'J')

# eukaryote indices-------
euk_indices <- read_csv('data/euk_indices_raw')

N_euk <- 
  ggplot(euk_indices_raw, aes(x = Treatment, y = N, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'N') + scale_y_log10()

S_euk <- 
  ggplot(euk_indices_raw, aes(x = Treatment, y = S, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'S')

J_euk <- 
  ggplot(euk_indices_raw, aes(x = Treatment, y = J, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'J')

# bacteria
bact_indices_raw <- read_csv('data/cleaned_data/bact_indices_raw')

N_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = N, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'N') + scale_y_log10()

S_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = S, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'S')

J_bact <- 
  ggplot(bact_indices_raw, aes(x = Treatment, y = J, fill = Treatment)) +
  geom_boxplot(alpha = .8) +
  scale_fill_viridis(guide = F, discrete = T) +
  labs(x='', y = 'J')


indices_plots <- 
  ggarrange(
    S_mac + ggtitle('Macrofauna') + ylab('Richness'),
    N_mac + ggtitle('') + ylab('Abundance'),
    J_mac + ggtitle('') + ylab('Evenness'),
    S_euk + ggtitle('Eukaryotes')  + ylab('Richness'),
    N_euk + ggtitle('') + ylab('Abundance'),
    J_euk + ggtitle('') + ylab('Evenness'),
    S_bact + ggtitle('Bacteria') + ylab('Richness'),
    N_bact + ggtitle('') + ylab('Total abundance'),
    J_bact + ggtitle('') + ylab('Evenness'),
    common.legend = T,
    # labels = "auto",
    ncol = 3,
    nrow = 3,
    legend = 'bottom'
  )
indices_plots

ggsave(
  indices_plots,
  filename = 'figures/_all_indices_plots.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 5,
  dpi = 900,
  width = 7
)


# relative abundace by phylum----------------------
# macrofauna
phylum_relabund_mac <- 
  ggplot(dat %>% drop_na(Phylum),
         aes(
           x = Phylum,
           y = Abundance,
           color = Treatment
         )) +
  stat_summary(fun.data = mean_se, size = .5,
               position = position_dodge(width = .5), alpha = .8) +
  scale_color_viridis_d() +
  labs(y = '', x = '') +
  coord_flip() +
  theme(legend.position = c(.8,.8))

# euk
phylum_relabund_euk <- 
  ggplot(dat_euk %>% drop_na(Phylum),
         aes(
           x = Phylum,
           y = Abundance,
           color = Treatment
         )) +
  stat_summary(fun.data = mean_se, size = .5,
               position = position_dodge(width = .5), alpha = .8) +
  scale_color_viridis_d() +
  labs(y = 'Relative abundance (%)', x='') +
  coord_flip() +
  theme(legend.position = c(.8,.8))


# bacteria------
phylum_relabund_bact <- 
  ggplot(dat_bact %>% filter(Phylum !=0),
         aes(
           x = Phylum,
           y = Abundance,
           color = Treatment
         )) +
  stat_summary(fun.data = mean_se, size = .5,
               position = position_dodge(width = .5), alpha = .8) +
  scale_color_viridis_d() +
  labs(y = '', x = '') +
  coord_flip() +
  theme(legend.position = c(.8,.8))



phylum_relabund_plots <- 
  ggarrange(
    phylum_relabund_mac + ggtitle('Macrofauna'),
    phylum_relabund_euk + ggtitle('Eukaryotes'),
    phylum_relabund_bact + ggtitle('Bacteria'),
    common.legend = T,
    # labels = "auto",
    ncol = 3,
    nrow = 1,
    legend = 'bottom'
  )


phylum_relabund_plots

ggsave(
  phylum_relabund_plots,
  filename = 'figures/phylum_relabund_plots.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 6 ,
  dpi = 900,
  width = 9
)
