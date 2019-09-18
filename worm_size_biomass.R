library(tidyverse)
library(readxl)
## read data----
worms <- read_excel('data/sabella_worms_data.xlsx', 2)
names(worms)


long_worm <- 
  worms %>% 
  group_by(sample) %>% 
  gather(id, size, length:length__34) %>% 
  select(-id)

size_hist <- 
  ggplot(drop_na(long_worm)) +
  geom_histogram(aes(x = size * 10),
                 color = 1,
                 fill = 'gray80',
                 bins = 15) +
  theme_minimal() +
  labs(x = "Size (mm)" , y = 'Number of individuals') +
  scale_x_continuous(breaks = seq(100,600,100)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(
  size_hist,
  filename = 'figures/size_histogram.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 2,
  width = 3
)


summary(long_worm$size)
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
se(long_worm$size)

density <- 
  long_worm %>%
  drop_na(size) %>% 
  group_by(sample) %>%
  summarise_all(funs(n()))


worm_survival_dat <- 
  read_csv('data/treatments.csv') %>% 
  mutate(final_density = ifelse(is.na(final_density), 12.5, final_density)) %>% 
  mutate(Survival = (final_density/Density)*100,
         loss = 100 - Survival)
  
worm_survival_dat %>% 
  group_by(Treatment) %>% 
  summarise(mean_loss = mean(loss),
            se_loss = se(loss))

worm_survival_dat %>% 
  group_by(Treatment, Density) %>% 
  summarise(mean_loss = mean(loss),
            se_loss = se(loss))

worm_loss_plot <- 
  worm_survival_dat %>% 
  filter(Density>0) %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(Density), loss, fill = Treatment)) +
  theme_bw() +
  labs(y = "Percentage loss (%)", x = 'Density', parse = T)

ggsave(
  worm_loss_plot,
  filename = 'figures/worm_loss_boxplot.tif',
  device = 'tiff',
  compression = 'lzw',
  height = 3,
  width = 4
)
