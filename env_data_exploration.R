library(tidyverse)
library(ggord)
library(emmeans)
source('theme_javier.R')
theme_set(theme_javier())


# env data exploration-----
env <-
  read_csv('data/env_vars.csv') %>%
  arrange(Plot) %>%
  rename(
    Phaeopigments = Phaeo,
    `Very fine sand` = Very_fine,
    `Chlorophyll-a` = Chla,
    `Fine sand` = Fine,
    `Medium sand` = Medium,
    `Coarse sand` = Coarse
  )

env_ord <-
  prcomp(
    env %>% select(Gravel:Phaeopigments) %>% log10(),
    scale. = T,
    center = T
  )
ggord(
  env_ord,
  grp_in = env$Treatment,
  poly = F,
  ellipse = F,
  repel = T,
  txt = 3,
  addsize = NA,
  arrow = 0,
  alpha = .7,
  veclsz = .2,
  vec_ext = 5,
  grp_title = NULL
) +
  theme_javier(base_size = 10) +
  scale_color_viridis_d() +
  theme(aspect.ratio = 1)

env_long <-
  env %>%
  gather(key, value, Gravel:Phaeopigments) %>%
  mutate(
    Treatment = fct_recode(Treatment,
                           C = "Control",
                           M = "Mimic",
                           S = "Sabella"),
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
  ) %>%
  write_csv('data/output_data/env_long.csv')

ggplot(env_long,
       aes(x = Treatment, y = value, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_viridis_d(guide = F, alpha = .7) +
  theme_javier(base_size = 10) +
  theme(axis.title = element_blank()) +
  facet_wrap(~ key, scales = 'free', nrow = 2)

summary(env)


env %>%
  mutate(silt_clay = Silt + Clay) %>%
  gather(key, value, Gravel:silt_clay) %>%
  group_by(key) %>%
  summarise(mean = mean(value),
            se = sd(value) / sqrt(length(value)))


# ancova indices------
env_ancova <-
  env_long %>%
  group_by(key) %>%
  nest() %>%
  mutate(
    ancova = map(.x = data, ~ lm(value ~ Treatment  + Density, data = .x)),
    ancova_sum = map(ancova, summary),
    anova_tab = map(ancova, ~ car::Anova(., type = "III")),
    ancova_table = map(anova_tab, broom::tidy),
    residuals = map(ancova, broom::augment),
    emm = map(ancova,  ~ emmeans(., ~ Treatment)),
    posthoc = map(emm, pairs)
  )

env_ancova %>%
  select(key, ancova_table) %>%
  unnest(ancova_table) %>%
  print(n = 50)


env_ancova %>%
  select(key, posthoc) %>%
  flatten()


ggplot(env_long, aes(x = Density, y = value)) +
  geom_point() +
  geom_smooth() +
  facet_wrap( ~ key, scales = 'free')
