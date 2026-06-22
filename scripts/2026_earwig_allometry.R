#rm(list=ls())
# if git is ahead by X commits do this: git reset --soft HEAD~1 (8=# of commits)
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))


# We show that weapon allometry is structured across nested biological levels, with sex determining 
# scaling plasticity, morphs diverging primarily in intercept, environmental conditions modifying scaling parameters, 
# and families exhibiting integrated intercept–slope covariance.

## ---- analysis ----
library(scales)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(broman)
library(cowplot)
library(flextable)
library(officer)
library(knitr)
library(officedown)
library(broom.mixed)
library(stringr)
library(performance)
library(brms)
library(tidybayes)
library(mixtools)
library(mixsmsn)
library(cmdstanr)
library(ggside)
library(posterior)
library(emmeans)

# packageVersion("mixsmsn")

earwig_data_raw <- read.csv(file="data/raw/earwig_allometry.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()

earwig_data <- earwig_data_raw %>%
  mutate(sex = recode(sex, femelle = "female",
                      male = "male")) %>%
  rename(id=Id_adulte) %>%
  mutate(density = as.factor(density), diet = as.factor(diet), sex = as.factor(sex), id_mere = as.factor(id_mere))

# select complete replicates: density-1 with 1, density-4 with 3-4, density-8 with 6, 7,or 8
earwig_two <- earwig_data %>%
  group_by(boite_petri, diet, density) %>%
  summarise(num_individ = n()) %>%
  filter(density=="8" & num_individ >= 6 | density=="4" & num_individ >= 3 | density=="1" & num_individ == 1) %>%
  mutate(valid=1) %>%
  ungroup() %>%
  dplyr::select(boite_petri, valid) 

dat.new <- left_join(earwig_data,earwig_two, by="boite_petri") %>%
  filter(valid==1) %>%
  mutate(treat_date = dmy(date_treatment), adult_date=dmy(date_adulte)) %>%
  mutate(dev_time= as.numeric(difftime(adult_date,treat_date,units = "days"))) %>%
  filter(sex!="") %>%
  droplevels()

# ================================
# Morph analysis
# ================================

## Males ##
# males <- dat.new %>%
#   filter(sex=="male")
# 
# bimodal.analysis <- smsn.mix(males$forceps_L, nu=3, g=2, get.init = TRUE, criteria = TRUE,
#                              group = TRUE, family = "Skew.normal", iter.max = 1000, calc.im=TRUE, obs.prob = TRUE)
# # bimodal.analysis$aic
# mixsmsn::mix.print(bimodal.analysis, digits=3)
# saveRDS(bimodal.analysis, file ="data/processed/bimodal.analysis.rda")
# str(bimodal.analysis)
# 
# mix.hist(males$forceps_L,bimodal.analysis)
# 
# # stick the group number from the mixed analysis onto the "males" data tibble
# # and convert group to morph
# ## bimodal
# males$group <- bimodal.analysis$group
# males <- males %>% mutate(
#   morph = case_when(
#     group == 2 ~ "minor",
#     group == 1 ~ "major")
# )
# 
# males$group <- as.factor(males$group)

# save(males,file="data/processed/males.Rda")

# males %>%
# arrange(forceps_L) %>%
#   as.tibble() %>%
#   print(n=1000) # to determine brachy vs macro cut-off

## FEMALES ##
# females <- dat.new %>%
#   filter(sex=="female") %>%
#   mutate(group=3, morph="female")
# str(females)
# 
# # save(females,file="data/processed/females.Rda")
# 
# females$group <- as.factor(females$group)
# 
# ## COMBINE DATA FRAMES
# dat.morphs <- rbind(males, females)

# save(dat.morphs, file = "data/processed/dat.morphs.Rda")
load(file ="data/processed/dat.morphs.Rda")

# # smsn.mix not categorizing individual correctly: he's clearly a minor
# dat.morphs <- dat.morphs %>%
#   mutate(group=replace(group, id=="DJODI", "2"))
# # 
# dat.morphs %>%
#   filter(morph=="CCNJO")

# # earwig_data_complete %>%
# #   filter(id_mere=="MNLDO")
# 
# # save(dat.morphs,file="data/processed/dat.morphs.Rda")

# ============================
# Allometry: three-group model
# ============================

dat.morphs <- dat.morphs %>%
  mutate(
    logF = log(forceps_L),
    logP = log(pronotum),
    
    group3 = case_when(
      sex == "female" ~ "female",
      sex == "male" & morph == "minor" ~ "brachylabic",
      sex == "male" & morph == "major" ~ "macrolabic"
    ),
    
    group3  = factor(group3, levels = c("female", "brachylabic", "macrolabic")),
    diet    = factor(diet, levels = c("GOOD", "POOR")),
    density = factor(density, levels = c("1", "4", "8")),
    id_mere = factor(id_mere)
  ) %>%
  filter(!is.na(group3), !is.na(logF), !is.na(logP))

str(dat.morphs)

summary_group3 <- dat.morphs %>%
  group_by(group3, diet, density) %>%
  summarise(
    n = n(),
    mean_forceps = mean(forceps_L, na.rm = TRUE),
    se_forceps = sd(forceps_L, na.rm = TRUE) / sqrt(n),
    mean_pronotum = mean(pronotum, na.rm = TRUE),
    se_pronotum = sd(pronotum, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

summary.table.2 <- dat.morphs %>%
  group_by(group3) %>%
  summarise(
    mean_forceps  = mean(forceps_L, na.rm = TRUE),
    sd_forceps    = sd(forceps_L, na.rm = TRUE),
    n             = n(),
    se_forceps    = sd_forceps / sqrt(n),
    
    mean_pronotum = mean(pronotum, na.rm = TRUE),
    sd_pronotum   = sd(pronotum, na.rm = TRUE),
    se_pronotum   = sd_pronotum / sqrt(n),
    
    .groups = "drop"
  ) %>%
  mutate(
    group3 = factor(
      group3,
      levels = c("female", "macrolabic", "brachylabic")
    )
  ) %>%
  arrange(group3)

formula_all <- bf(
  logF ~ logP * group3 * diet * density +
    (1 + logP | id_mere)
)

priors_all <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(student_t(3, 0, 1), class = "sigma"),
  prior(student_t(3, 0, 1), class = "sd"),
  prior(lkj(2), class = "cor")
)

# mod_all <- brm(
#   formula = formula_all,
#   data    = dat.morphs,
#   family  = gaussian(),
#   prior   = priors_all,
#   chains  = 4,
#   cores   = 4,
#   iter    = 4000,
#   backend = "cmdstanr",
#   file    = "data/processed/mod_all.Rds",
#   control = list(adapt_delta = 0.97)
# )

mod_all <- readRDS(file = "data/processed/mod_all.Rds")

summary(mod_all)
print(tidy(mod_all), n = Inf)

conditional_effects(mod_all)

# ==============================
# Test for non-linearity
# ==============================

formula_all_quad <- bf(
  logF ~ logP * group3 * diet * density +
    I(logP^2) * group3 * diet * density +
    (1 + logP | id_mere)
)

priors_all_quad <- priors_all

# mod_all_quad <- brm(
#   formula = formula_all_quad,
#   data    = dat.morphs,
#   family  = gaussian(),
#   prior   = priors_all_quad,
#   backend = "cmdstanr",
#   chains  = 4,
#   cores   = 4,
#   iter    = 4000,
#   file    = "data/processed/mod_all_quad.Rds",
#   control = list(adapt_delta = 0.97)
# )

mod_all_quad <- readRDS(file = "data/processed/mod_all_quad.Rds")

# Add LOO criteria
mod_all <- add_criterion(mod_all, "loo")
mod_all_quad <- add_criterion(mod_all_quad, "loo")

# Compare linear vs quadratic model
loo_compare(mod_all, mod_all_quad)

loo_comp <- loo_compare(mod_all, mod_all_quad)

loo_table <- as.data.frame(loo_comp) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(
    Model = recode(Model,
                   mod_all = "Linear model",
                   mod_all_quad = "Quadratic model")
  )

loo_table

# ============================
# Condition-specific slopes
# ============================

slopes_all <- emtrends(
  mod_all,
  ~ group3 * diet * density,
  var = "logP"
)

slope_table <- as.data.frame(slopes_all)

# ============================
# Group comparisons within each environment
# ============================

group_slope_contrasts <- pairs(
  emtrends(mod_all, ~ group3 | diet * density, var = "logP")
)

# ============================
# Environmental slope plasticity within each group
# ============================

diet_slope_contrasts <- pairs(
  emtrends(mod_all, ~ diet | group3 * density, var = "logP")
)

density_slope_contrasts <- pairs(
  emtrends(mod_all, ~ density | group3 * diet, var = "logP")
)

summary(diet_slope_contrasts)
summary(density_slope_contrasts)

# ============================
# Posterior probabilities for slope contrasts
# ============================

x_seq <- seq(
  min(dat.morphs$logP, na.rm = TRUE),
  max(dat.morphs$logP, na.rm = TRUE),
  length.out = 100
)

new_all <- expand.grid(
  logP    = x_seq,
  group3  = levels(dat.morphs$group3),
  diet    = levels(dat.morphs$diet),
  density = levels(dat.morphs$density),
  id_mere = NA
)

epred_all <- posterior_epred(
  mod_all,
  newdata = new_all,
  re_formula = NA
)

new_all <- new_all %>%
  mutate(
    fit = apply(epred_all, 2, median),
    lwr = apply(epred_all, 2, quantile, 0.025),
    upr = apply(epred_all, 2, quantile, 0.975),
    panel = recode(group3,
                   female = "Female",
                   brachylabic = "Brachylabic male",
                   macrolabic = "Macrolabic male")
  )

raw_data <- dat.morphs %>%
  mutate(
    panel = recode(group3,
                   female = "Female",
                   brachylabic = "Brachylabic male",
                   macrolabic = "Macrolabic male")
  )

dat.morphs %>%
  count(group3, diet, density)

## ---- end

# ==============================
# PLOTS
# ==============================

library(patchwork)

#### figure 2 ####

dat.morphs.2 <- dat.morphs %>% mutate(
  morph = case_when(
    group == "2" ~ "brachylabic",
    group == "1" ~ "macrolabic",
    group == "3" ~ "female")
)

dat.morphs.2$morph <- factor(dat.morphs.2$morph, levels=c('female', 'macrolabic', 'brachylabic'))
dat.morphs.2$diet <- factor(dat.morphs.2$diet, levels=c('POOR', 'GOOD'))
library(ggrepel)

density.labs <- c("Low", "Medium","High")
names(density.labs) <- c("1", "4", "8")

sex.labs <- c("Females", "Males")
names(sex.labs) <- c("female", "male")

forceps.body.plot.both <- ggplot(dat.morphs.2, aes(x=pronotum, y=forceps_L,label=id_mere,shape=interaction(morph,diet), colour=interaction(morph,diet))) +
  geom_point(size=2,alpha=0.7) +
  geom_hline(data = dat.morphs.2 %>% filter(sex == "male"),
             aes(yintercept = 4.725), colour="grey", linetype="dashed") +
  facet_grid(sex~density,labeller = labeller(density = density.labs, sex = sex.labs)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_colour_manual("",values=c("black", "black", "black", "red", "red","black"), 
                      labels = c(
                        "female.POOR" = "Female - Poor diet",
                        "female.GOOD" = "Female - Good diet",
                        "brachylabic.POOR" = "Brachylabic male - Poor diet",
                        "brachylabic.GOOD" = "Brachylabic male - Good diet",
                        "macrolabic.POOR" = "Macrolabic male - Poor diet",
                        "macrolabic.GOOD"= "Macrolabic male - Good diet"),
                      limits = c("female.POOR","female.GOOD","brachylabic.POOR", "brachylabic.GOOD", "macrolabic.POOR", "macrolabic.GOOD")) +
  scale_shape_manual("",values=c(16,1,16,16,1,1), 
                     labels = c(
                       "female.POOR" = "Female - Poor diet",
                       "female.GOOD" = "Female - Good diet",
                       "brachylabic.POOR" = "Brachylabic male - Poor diet",
                       "brachylabic.GOOD" = "Brachylabic male - Good diet",
                       "macrolabic.POOR" = "Macrolabic male - Poor diet",
                       "macrolabic.GOOD"= "Macrolabic male - Good diet"),
                     limits = c("female.POOR","female.GOOD","brachylabic.POOR", "brachylabic.GOOD", "macrolabic.POOR", "macrolabic.GOOD")) +
  xlab("Pronotum length (mm)") +
  ylab("Forceps length (mm)")

figure_2 <- forceps.body.plot.both + geom_ysidedensity(aes(x=after_stat(density),group=sex, colour="black")) +
  ggside(collapse="y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(
    legend.position = c(0.02, 0.98),   # x, y in [0,1]
    legend.justification = c(0, 1)     # anchor legend's top-left corner
  ) +
  theme(legend.text = element_text(size = 12)) +
  theme(strip.text.y = element_text(size = 14, face="bold")) +
  theme(strip.text.x = element_text(size = 14, face="bold")) +
  theme(axis.title.x = element_text(size = 16, face="bold")) +
  theme(axis.title.y = element_text(size = 16, face="bold")) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(strip.background = element_rect(fill="white")) +
  theme(ggside.axis.text.x = element_blank()) +
  theme(ggside.axis.ticks.x = element_blank())

ggsave(figure_2,filename="Figure_2.jpg", width=14.83, height=8.83, dpi=300,antialias="default")


#### Figure 3: posterior slope estimates ####

slopes_all <- emtrends(
  mod_all,
  ~ group3 * diet * density,
  var = "logP"
)

slope_df <- as.data.frame(slopes_all) %>%
  mutate(
    group3 = recode(group3,
                    "female" = "Female",
                    "brachylabic" = "Brachylabic male",
                    "macrolabic" = "Macrolabic male"),
    
    diet = recode(diet,
                  "GOOD" = "Good",
                  "POOR" = "Poor"),
    
    density = recode(density,
                     "1" = "Low",
                     "4" = "Medium",
                     "8" = "High"),
    
    density = factor(density, levels = c("Low", "Medium", "High")),
    
    group3 = factor(group3,
                    levels = c("Female", "Brachylabic male", "Macrolabic male"))
  )

figure_3 <- ggplot(
  slope_df,
  aes(x = logP.trend,
      y = density,
      colour = diet)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dotted",
    linewidth = 0.6,
    alpha = 0.4
  ) +
  geom_errorbarh(
    aes(xmin = lower.HPD, xmax = upper.HPD),
    width = 0.2,
    linewidth = 0.6,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    size = 2.5,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ group3, nrow = 1) +
  scale_color_manual(
    values = c(
      "Good" = "#D55E00",
      "Poor" = "#0072B2"
    )
  ) +
  #coord_cartesian(xlim = c(-0.2, 1.8)) +
  labs(
    x = expression(bold(paste("Allometric slope (", beta, ")"))),
    y = "Rearing density",
    colour = "Diet"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = "right"
  )

figure_3

ggsave(
  filename = "figure_3.jpg",
  plot = figure_3,
  width = 12.83,
  height = 8.83,
  dpi = 300,
  antialias = "default"
)

# ==============================
# Figure S1: predicted allometries
# Three-group model
# ==============================

library(dplyr)
library(ggplot2)
library(brms)

# Ensure variables match the fitted model
dat.morphs <- dat.morphs %>%
  mutate(
    group3 = factor(group3,
                    levels = c("female", "brachylabic", "macrolabic")),
    diet = factor(diet, levels = c("GOOD", "POOR")),
    density = factor(density, levels = c("1", "4", "8")),
    panel = factor(
      recode(group3,
             female = "Female",
             brachylabic = "Brachylabic",
             macrolabic = "Macrolabic"),
      levels = c("Female", "Brachylabic", "Macrolabic")
    )
  )

# Prediction grid: only across observed logP range per cell
x_ranges <- dat.morphs %>%
  group_by(group3, diet, density) %>%
  summarise(
    min_logP = min(logP, na.rm = TRUE),
    max_logP = max(logP, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

new_all <- x_ranges %>%
  rowwise() %>%
  do(data.frame(
    group3 = .$group3,
    diet = .$diet,
    density = .$density,
    logP = seq(.$min_logP, .$max_logP, length.out = 50),
    n = .$n
  )) %>%
  ungroup() %>%
  mutate(
    group3 = factor(group3, levels = levels(dat.morphs$group3)),
    diet = factor(diet, levels = levels(dat.morphs$diet)),
    density = factor(density, levels = levels(dat.morphs$density)),
    id_mere = NA,
    panel = factor(
      recode(group3,
             female = "Female",
             brachylabic = "Brachylabic",
             macrolabic = "Macrolabic"),
      levels = c("Female", "Brachylabic", "Macrolabic")
    )
  )

# Posterior fitted values, excluding family-level effects
epred_all <- posterior_epred(
  mod_all,
  newdata = new_all,
  re_formula = NA
)

new_all <- new_all %>%
  mutate(
    fit = apply(epred_all, 2, median),
    lwr = apply(epred_all, 2, quantile, 0.025),
    upr = apply(epred_all, 2, quantile, 0.975)
  )

# Plot
figure_S1 <- ggplot(
  new_all,
  aes(
    x = logP,
    y = fit,
    colour = diet,
    linetype = density,
    group = interaction(diet, density)
  )
) +
  geom_point(
    data = dat.morphs,
    aes(x = logP, y = logF),
    inherit.aes = FALSE,
    colour = "grey50",
    alpha = 0.12,
    size = 0.7
  ) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ panel, nrow = 1) +
  labs(
    x = "log(Pronotum length, mm)",
    y = "log(Forceps length, mm)",
    colour = "Diet",
    linetype = "Density"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    strip.background = element_rect(fill = "white"),
    legend.position = "right"
  )

figure_S1

ggsave(
  filename = "figure_S1.jpg",
  plot = figure_S1,
  width = 10,
  height = 6,
  dpi = 300
)



