library(tidyverse)
source('theme_javier.R')
theme_set(theme_javier())


# env data exploration-----
env <-  
  read_csv('data/env_vars.csv') %>% 
  arrange(Plot) 
  


env %>% 
  gather(key, value, Gravel:Phaeo) %>% 
  mutate(Treatment = fct_recode(Treatment, C = "Control", M = "Mimic", S = "Sabella")) %>% 
  ggplot(aes(x = Treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~key, scales = 'free',nrow = 2)

summary(env)


env %>% 
  mutate(silt_clay = Silt+Clay) %>% 
  gather(key, value, Gravel:silt_clay) %>% 
  group_by(key) %>% 
  summarise(mean = mean(value),
            se = sd(value) /sqrt(length(value)))
