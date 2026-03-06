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
# # earwig_data_complete %>%
# #   filter(id=="DJODI")
# # earwig_data_complete %>%
# #   filter(id_mere=="MNLDO")
# 
# # save(dat.morphs,file="data/processed/dat.morphs.Rda")

# ============================
# Allometry
# ============================

# Transform variables
dat.morphs <- dat.morphs %>%
  mutate(
    logF  = log(forceps_L),
    logP  = log(pronotum),
    logPc = scale(logP, scale = FALSE)[,1],  # mean center
    sex   = factor(sex),
    group = factor(morph),
    diet  = factor(diet),
    density = factor(density)
  )

# Reference levels
dat.morphs$sex     <- relevel(dat.morphs$sex, ref = "female")
dat.morphs$diet    <- relevel(dat.morphs$diet, ref = "GOOD")
dat.morphs$density <- relevel(dat.morphs$density, ref = "1")

summary.table <- dat.morphs %>%
  group_by(morph,diet,density) %>%
  summarise(n=n())

summary.table.2 <- dat.morphs %>%
  group_by(morph) %>%
  summarize(mean_for = mean(forceps_L, na.rm = TRUE),
            sd.forceps = sd(forceps_L, na.rm = TRUE),
            n = n(), # counts the number of observations in the current group
            se.forceps = sd.forceps / sqrt(n), # standard error formula
            mean.pr = mean(pronotum,na.rm=TRUE),
            sd.pro = sd(pronotum, na.rm=TRUE),
            se.pro = sd.pro/sqrt(n),
            .groups = 'drop' # drops the grouping structure after summarizing
    )

# sex difference
formula_sex <- bf(
  logF ~ logPc * sex * diet +
    logPc * sex * density +
    (1 + logPc | id_mere)
)

priors_sex <- c(
  
  # Fixed effects (slopes & intercept shifts)
  prior(normal(0, 1), class = "b"),
  
  # Intercept
  prior(normal(0, 2), class = "Intercept"),
  
  # Residual SD
  prior(student_t(3, 0, 1), class = "sigma"),
  
  # Random effect SDs
  prior(student_t(3, 0, 1), class = "sd"),
  
  # Random effect correlation (intercept-slope)
  prior(lkj(2), class = "cor")
)

# mod_sex <- brm(
#   formula = formula_sex,
#   data    = dat.morphs,
#   family  = gaussian(),
#   prior   = priors_sex,
#   chains  = 4,
#   cores   = 4,
#   iter    = 4000,
#   backend = "cmdstanr",
#      file = "data/processed/mod_sex.Rds",
#   control = list(adapt_delta = 0.97)
# )

# saveRDS(mod_sex, file = "mod_sex.Rds")
mod_sex <- readRDS(file = "data/processed/mod_sex.Rds")

#### Compute Posterior Slopes per Sex × Environment

# Slopes for all combinations
slopes_sex <- emtrends(
  mod_sex,
  ~ sex * diet * density,
  var = "logPc"
)

summary(slopes_sex)
slopes_draws <- as.data.frame(slopes_sex)

#### Posterior Probability Male Slope > Female Slope
male_vs_female <- pairs(
  emtrends(mod_sex, ~ sex | diet * density, var = "logPc")
)

summary(male_vs_female)

# Extract posterior draws
draws <- as_draws_df(mod_sex)

# Female slope baseline
beta_female_base <- draws$b_logPc

# Male slope baseline
beta_male_base <- draws$b_logPc + draws$`b_logPc:sexmale`

# Posterior probability male > female (GOOD, density1)
mean(beta_male_base > beta_female_base)

#### Probability Male Slope Changes Under Poor Diet
beta_male_good <- draws$b_logPc + draws$`b_logPc:sexmale`

beta_male_poor <- beta_male_good +
  draws$`b_logPc:sexmale:dietPOOR`

mean(beta_male_poor < beta_male_good)

delta_male_diet <- beta_male_poor - beta_male_good
mean(abs(delta_male_diet) > 0.1)

#### Probability Female Slope Changes Under Poor Diet
beta_female_good <- draws$b_logPc
beta_female_poor <- draws$b_logPc + draws$`b_logPc:dietPOOR`

delta_female_diet <- beta_female_poor - beta_female_good
mean(delta_female_diet > 0)
mean(delta_female_diet < 0)

quantile(delta_female_diet, c(.025, .5, .975))
mean(abs(delta_female_diet) > 0.1)

#### Full Table of Posterior Slopes
# Female slopes
female_good_d1 <- draws$b_logPc
female_poor_d1 <- draws$b_logPc + draws$`b_logPc:dietPOOR`
female_good_d4 <- draws$b_logPc + draws$`b_logPc:density4`
female_good_d8 <- draws$b_logPc + draws$`b_logPc:density8`

quantile(female_good_d1, probs = c(.025, .5, .975))
quantile(female_poor_d1, probs = c(.025, .5, .975))
quantile(female_good_d4, probs = c(.025, .5, .975))
quantile(female_good_d8, probs = c(.025, .5, .975))

# Male slopes
male_good_d1 <- draws$b_logPc + draws$`b_logPc:sexmale`
male_poor_d1 <- male_good_d1 + draws$`b_logPc:sexmale:dietPOOR`
male_good_d4 <- male_good_d1 + draws$`b_logPc:sexmale:density4`
male_good_d8 <- male_good_d1 + draws$`b_logPc:sexmale:density8`

quantile(male_good_d1, probs = c(.025, .5, .975))
quantile(male_poor_d1, probs = c(.025, .5, .975))
quantile(male_good_d4, probs = c(.025, .5, .975))
quantile(male_good_d8, probs = c(.025, .5, .975))

mean(abs(delta_male_diet) > abs(delta_female_diet))

## density
delta_male_density4  <- male_good_d4 - male_good_d1
delta_male_density8  <- male_good_d8 - male_good_d1

mean(abs(delta_male_density4) > 0.1)
mean(abs(delta_male_density8) > 0.1)

delta_female_density4  <- female_good_d4 - female_good_d1
delta_female_density8  <- female_good_d8 - female_good_d1

mean(abs(delta_female_density4) > 0.1)
mean(abs(delta_female_density8) > 0.1)

### MALES ONLY
males <- subset(dat.morphs, sex == "male")
males$morph <- factor(males$morph)
males$morph <- relevel(males$morph, ref = "minor")  # brachylabic as reference

formula_morph <- bf(
  logF ~ logPc * morph * diet +
    logPc * morph * density +
    (1 + logPc | id_mere)
)

priors_morph <- c(
  
  # Population-level fixed effects (main effects + interactions)
  prior(normal(0, 0.5), class = "b"),
  
  # Intercept (log-scale forceps size)
  prior(normal(1.2, 0.5), class = "Intercept"),
  
  # Residual SD
  prior(student_t(3, 0, 0.5), class = "sigma"),
  
  # Random effect SDs (family variation)
  prior(student_t(3, 0, 0.3), class = "sd"),
  
  # Correlation between random intercept and slope
  prior(lkj(2), class = "cor")
)

# mod_morph <- brm(
#   formula = formula_morph,
#   data    = males,
#   family  = gaussian(),
#   prior   = priors_morph,
#   backend = "cmdstanr",
#   chains  = 4,
#   cores   = 4,
#   iter    = 4000,
#   file = "data/processed/mod_morph.Rds",
#   control = list(adapt_delta = 0.97)
# )
# 
# saveRDS(mod_morph, file = "mod_morph.Rds")

mod_morph <- readRDS(file = "data/processed/mod_morph.Rds")

draws <- as_draws_df(mod_morph)

# Baseline Slopes (GOOD diet, density 1):
# Do macrolabic and brachylabic males differ in their allometric slope under the reference environmental conditions?
beta_brachy_good_d1 <- draws$b_logPc

beta_macro_good_d1 <- draws$b_logPc + draws$`b_logPc:morphmajor`
delta_baseline <- beta_macro_good_d1 - beta_brachy_good_d1

mean(delta_baseline > 0)
quantile(delta_baseline, c(.025, .5, .975)) # baseline slopes are not different

# Diet Slope Plasticity Within Each Morph
beta_brachy_poor_d1 <- beta_brachy_good_d1 + draws$`b_logPc:dietPOOR`
beta_macro_poor_d1 <- beta_macro_good_d1 +
  draws$`b_logPc:dietPOOR` +
  draws$`b_logPc:morphmajor:dietPOOR`
delta_brachy_diet <- beta_brachy_poor_d1 - beta_brachy_good_d1
delta_macro_diet  <- beta_macro_poor_d1  - beta_macro_good_d1

quantile(delta_brachy_diet, c(.025, .5, .975))
quantile(delta_macro_diet,  c(.025, .5, .975))

mean(abs(delta_brachy_diet) > 0.1)
mean(abs(delta_macro_diet)  > 0.1)

mean(abs(delta_macro_diet) > abs(delta_brachy_diet))

# Density slope plasticity
# density 4
beta_brachy_d4 <- beta_brachy_good_d1 + draws$`b_logPc:density4`
beta_macro_d4  <- beta_macro_good_d1  +
  draws$`b_logPc:density4` +
  draws$`b_logPc:morphmajor:density4`

delta_brachy_d4 <- beta_brachy_d4 - beta_brachy_good_d1
delta_macro_d4  <- beta_macro_d4  - beta_macro_good_d1

mean(abs(delta_macro_d4) > abs(delta_brachy_d4))

# Density 8
beta_brachy_d8 <- beta_brachy_good_d1 + draws$`b_logPc:density8`
beta_macro_d8  <- beta_macro_good_d1  +
  draws$`b_logPc:density8` +
  draws$`b_logPc:morphmajor:density8`

delta_brachy_d8 <- beta_brachy_d8 - beta_brachy_good_d1
delta_macro_d8  <- beta_macro_d8  - beta_macro_good_d1

mean(abs(delta_macro_d8) > abs(delta_brachy_d8))

summ <- function(x) {
  c(
    median = median(x),
    l95 = quantile(x, .025),
    u95 = quantile(x, .975),
    p_gt0 = mean(x > 0),
    p_abs_gt_01 = mean(abs(x) > 0.1)
  )
}

summ(delta_brachy_diet)
summ(delta_macro_diet)

## ---- end

# ==============================
# PLOTS
# ==============================

library(patchwork)

#### figure 1 ####

dat.morphs.2 <- dat.morphs %>% mutate(
  morph = case_when(
    group == 2 ~ "minor",
    group == 1 ~ "major",
    group == 3 ~ "female")
)

dat.morphs.2$morph <- factor(dat.morphs.2$morph, levels=c('female', 'major', 'minor'))
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
                        "female.POOR" = "Female - Poor",
                        "female.GOOD" = "Female - Good",
                        "minor.POOR" = "Brachylabic - Poor",
                        "minor.GOOD" = "Brachylabic - Good",
                        "major.POOR" = "Macrolabic - Poor",
                        "major.GOOD"= "Macrolabic - Good",
                        "old.low"="Old - Low",
                        "old.high"="Old - High"),
                      limits = c("female.POOR","female.GOOD","minor.POOR", "minor.GOOD", "major.POOR", "major.GOOD")) +
  scale_shape_manual("",values=c(16,1,16,16,1,1), 
                     labels = c(
                       "female.POOR" = "Female - Poor",
                       "female.GOOD" = "Female - Good",
                       "minor.POOR" = "Brachylabic - Poor",
                       "minor.GOOD" = "Brachylabic - Good",
                       "major.POOR" = "Macrolabic - Poor",
                       "major.GOOD"= "Macrolabic - Good",
                       "old.low"="Old - Low",
                       "old.high"="Old - High"),
                     limits = c("female.POOR","female.GOOD","minor.POOR", "minor.GOOD", "major.POOR", "major.GOOD")) +
  xlab("Pronotum length (mm)") +
  ylab("Forceps length (mm)")

figure_1 <- forceps.body.plot.both + geom_ysidedensity(aes(x=after_stat(density),group=sex, colour="black")) +
  ggside(collapse="y") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.text.y = element_text(size = 14, face="bold")) +
  theme(strip.text.x = element_text(size = 14, face="bold")) +
  theme(axis.title.x = element_text(size = 16, face="bold")) +
  theme(axis.title.y = element_text(size = 16, face="bold")) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(strip.background = element_rect(fill="white")) +
  theme(ggside.axis.text.x = element_blank()) +
  theme(ggside.axis.ticks.x = element_blank())

ggsave(figure_1,filename="Figure_1.jpg", width=14.83, height=8.83, dpi=800,antialias="default")

#### figure 2 ####

# Range of logPc
x_seq <- seq(
  min(dat.morphs$logPc, na.rm = TRUE),
  max(dat.morphs$logPc, na.rm = TRUE),
  length.out = 100
)

# Grid for females (from sex model)
new_female <- expand.grid(
  logPc   = x_seq,
  sex     = "female",
  diet    = c("GOOD", "POOR"),
  density = c("1", "4", "8"),
  id_mere = NA
)

# Grid for males (from morph model)
new_male <- expand.grid(
  logPc   = x_seq,
  group   = c("1","2"),   # 1=macro, 2=brachy
  diet    = c("GOOD", "POOR"),
  density = c("1", "4", "8"),
  id_mere = NA
)

# Females from sex model
epred_female <- posterior_epred(
  mod_sex,
  newdata = new_female,
  re_formula = NA   # exclude random effects
)

# Males from morph model
epred_male <- posterior_epred(
  mod_morph,
  newdata = new_male,
  re_formula = NA
)

# Female summaries
female_summary <- apply(epred_female, 2, median)

new_female$fit <- female_summary
new_female$panel <- "Female"

# Male summaries
male_summary <- apply(epred_male, 2, median)

new_male$fit <- male_summary
new_male$panel <- ifelse(new_male$group == "1",
                         "Macrolabic male",
                         "Brachylabic male")

plot_data <- bind_rows(
  new_female %>% select(logPc, diet, density, fit, panel),
  new_male %>% select(logPc, diet, density, fit, panel)
)

plot_data$panel <- factor(
  plot_data$panel,
  levels = c(
    "Female",
    "Brachylabic male",
    "Macrolabic male"
  )
)

raw_data <- dat.morphs %>%
  mutate(
    panel = case_when(
      sex == "female" ~ "Female",
      sex == "male" & group == "2" ~ "Brachylabic male",
      sex == "male" & group == "1" ~ "Macrolabic male"
    )
  )

raw_data$panel <- factor(
  raw_data$panel,
  levels = c("Female", "Brachylabic male", "Macrolabic male")
) 

figure_2 <- ggplot(plot_data,
       aes(x = logPc,
           y = fit,
           color = diet,
           linetype = density,
           group = interaction(diet, density))) +
  geom_point(
    data = raw_data,
    aes(x = logPc, y = logF),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.08,
    size = 0.7
  ) + 
  geom_line(size = 1.2) +
  facet_wrap(~ panel, nrow = 1, scales = "fixed") +
  scale_color_manual(values = c("GOOD" = "#D55E00",
                                "POOR" = "#0072B2")) +
  scale_linetype_manual(values = c("1" = "solid",
                                   "4" = "dashed",
                                   "8" = "dotted")) +
  labs(
    x = "log(Pronotum)",
    y = "log(Forceps length)",
    color = "Diet",
    linetype = "Density"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(size = 14, face="bold"),
    strip.text.x = element_text(size = 14, face="bold"),
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    ggside.axis.text.x = element_blank(),
    ggside.axis.ticks.x = element_blank(),
    legend.position = c(0.1, 0.75)
  )

ggsave(figure_2, filename="figure_2.jpg", width=12.83, height=8.83, dpi=300,antialias="default")


#### figure 3 ####

draws_sex   <- as_draws_df(mod_sex)
draws_morph <- as_draws_df(mod_morph)

# female slopes
female_good_d1 <- draws_sex$b_logPc
female_poor_d1 <- draws_sex$b_logPc + draws_sex$`b_logPc:dietPOOR`
female_good_d4 <- draws_sex$b_logPc + draws_sex$`b_logPc:density4`
female_good_d8 <- draws_sex$b_logPc + draws_sex$`b_logPc:density8`

female_poor_d4 <- female_poor_d1 + draws_sex$`b_logPc:density4`
female_poor_d8 <- female_poor_d1 + draws_sex$`b_logPc:density8`

# brachy slopes
brachy_good_d1 <- draws_morph$b_logPc
brachy_poor_d1 <- brachy_good_d1 + draws_morph$`b_logPc:dietPOOR`
brachy_good_d4 <- brachy_good_d1 + draws_morph$`b_logPc:density4`
brachy_good_d8 <- brachy_good_d1 + draws_morph$`b_logPc:density8`

brachy_poor_d4 <- brachy_poor_d1 + draws_morph$`b_logPc:density4`
brachy_poor_d8 <- brachy_poor_d1 + draws_morph$`b_logPc:density8`

# macro slopes
macro_good_d1 <- brachy_good_d1 + draws_morph$`b_logPc:group1`

macro_poor_d1 <- macro_good_d1 +
  draws_morph$`b_logPc:dietPOOR` +
  draws_morph$`b_logPc:group1:dietPOOR`

macro_good_d4 <- macro_good_d1 +
  draws_morph$`b_logPc:density4` +
  draws_morph$`b_logPc:group1:density4`

macro_good_d8 <- macro_good_d1 +
  draws_morph$`b_logPc:density8` +
  draws_morph$`b_logPc:group1:density8`

macro_poor_d4 <- macro_poor_d1 +
  draws_morph$`b_logPc:density4` +
  draws_morph$`b_logPc:group1:density4`

macro_poor_d8 <- macro_poor_d1 +
  draws_morph$`b_logPc:density8` +
  draws_morph$`b_logPc:group1:density8`

summ <- function(x, label, group) {
  data.frame(
    group = group,
    treatment = label,
    median = median(x),
    l95 = quantile(x, .025),
    u95 = quantile(x, .975)
  )
}

slope_df <- bind_rows(
  
  # Female
  summ(female_good_d1, "GOOD d1", "Female"),
  summ(female_good_d4, "GOOD d4", "Female"),
  summ(female_good_d8, "GOOD d8", "Female"),
  summ(female_poor_d1, "POOR d1", "Female"),
  summ(female_poor_d4, "POOR d4", "Female"),
  summ(female_poor_d8, "POOR d8", "Female"),
  
  # Brachy
  summ(brachy_good_d1, "GOOD d1", "Brachylabic male"),
  summ(brachy_good_d4, "GOOD d4", "Brachylabic male"),
  summ(brachy_good_d8, "GOOD d8", "Brachylabic male"),
  summ(brachy_poor_d1, "POOR d1", "Brachylabic male"),
  summ(brachy_poor_d4, "POOR d4", "Brachylabic male"),
  summ(brachy_poor_d8, "POOR d8", "Brachylabic male"),
  
  # Macro
  summ(macro_good_d1, "GOOD d1", "Macrolabic male"),
  summ(macro_good_d4, "GOOD d4", "Macrolabic male"),
  summ(macro_good_d8, "GOOD d8", "Macrolabic male"),
  summ(macro_poor_d1, "POOR d1", "Macrolabic male"),
  summ(macro_poor_d4, "POOR d4", "Macrolabic male"),
  summ(macro_poor_d8, "POOR d8", "Macrolabic male")
)

slope_df$treatment <- factor(
  slope_df$treatment,
  levels = c(
    "GOOD d1", "GOOD d4", "GOOD d8",
    "POOR d1", "POOR d4", "POOR d8"
  )
)

slope_df$treatment <- rev(levels(slope_df$treatment))

slope_df$group <- factor(
  slope_df$group,
  levels = c("Female", "Brachylabic male", "Macrolabic male")
) 


figure_3 <- ggplot(slope_df,
       aes(x = median,
           y = treatment)) +
  
  # vertical reference line at 0
  geom_vline(xintercept = 0,
             linetype = "dotted",
             size = 0.6) +
  
  # credible intervals
  geom_errorbarh(aes(xmin = l95, xmax = u95),
                 height = 0.2,
                 size = 0.6,
                 color = "black") +
  
  # posterior medians
  geom_point(size = 2.5,
             color = "black") +
  
  facet_wrap(~ group, scales = "fixed") +
  
  scale_y_discrete(limits = rev) +
  
  labs(
    x = expression(paste("Allometric slope (", beta, ")")),
    y = NULL
  ) +
  
  geom_hline(yintercept = 3.5, color = "grey80") +
  
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(size = 14, face="bold"),
    strip.text.x = element_text(size = 14, face="bold"),
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    ggside.axis.text.x = element_blank(),
    ggside.axis.ticks.x = element_blank(),
    legend.position = c(0.1, 0.75)
  )

ggsave(figure_3, filename="figure_3.jpg", width=12.83, height=8.83, dpi=300,antialias="default")

#### figure 4 ####
draws_sex <- as_draws_df(mod_sex)

# Female
delta_f_diet <- draws_sex$`b_logPc:dietPOOR`
delta_f_d4   <- draws_sex$`b_logPc:density4`
delta_f_d8   <- draws_sex$`b_logPc:density8`

# Male
delta_m_diet <- draws_sex$`b_logPc:dietPOOR` +
  draws_sex$`b_logPc:sexmale:dietPOOR`

delta_m_d4 <- draws_sex$`b_logPc:density4` +
  draws_sex$`b_logPc:sexmale:density4`

delta_m_d8 <- draws_sex$`b_logPc:density8` +
  draws_sex$`b_logPc:sexmale:density8`

summ <- function(x, group, env){
  data.frame(
    group = group,
    environment = env,
    median = median(x),
    l95 = quantile(x,.025),
    u95 = quantile(x,.975)
  )
}

delta_sex_df <- bind_rows(
  summ(delta_f_diet,"Female","POOR diet"),
  summ(delta_f_d4,"Female","Density 4"),
  summ(delta_f_d8,"Female","Density 8"),
  summ(delta_m_diet,"Male","POOR diet"),
  summ(delta_m_d4,"Male","Density 4"),
  summ(delta_m_d8,"Male","Density 8")
)

figure_4 <- ggplot(delta_sex_df,
       aes(x = median, y = group)) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept=c(-0.1,0.1),
             linetype="dotted", colour="grey50") +
  geom_errorbarh(aes(xmin=l95, xmax=u95),
                 height=.15, size=.9) +
  geom_point(size=3) +
  facet_wrap(~environment) +
  coord_cartesian(xlim = c(-0.85, 0.7)) +
  labs(x=expression(Delta*beta), y=NULL) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(size = 14, face="bold"),
    strip.text.x = element_text(size = 14, face="bold"),
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    ggside.axis.text.x = element_blank(),
    ggside.axis.ticks.x = element_blank(),
    legend.position = c(0.1, 0.75)
  )

ggsave(figure_4, filename="figure_4.jpg", width=12.83, height=8.83, dpi=300,antialias="default")

#### figure 5 ####
draws_morph <- as_draws_df(mod_morph)

# Brachylabic (reference morph)
delta_b_diet <- draws_morph$`b_logPc:dietPOOR`
delta_b_d4   <- draws_morph$`b_logPc:density4`
delta_b_d8   <- draws_morph$`b_logPc:density8`

# Macrolabic
delta_M_diet <- draws_morph$`b_logPc:dietPOOR` +
  draws_morph$`b_logPc:morphmajor:dietPOOR`

delta_M_d4 <- draws_morph$`b_logPc:density4` +
  draws_morph$`b_logPc:morphmajor:density4`

delta_M_d8 <- draws_morph$`b_logPc:density8` +
  draws_morph$`b_logPc:morphmajor:density8`

summ <- function(x, group, env){
  data.frame(
    group = group,
    environment = env,
    median = median(x),
    l95 = quantile(x,.025),
    u95 = quantile(x,.975)
  )
}

delta_morph_df <- bind_rows(
  summ(delta_b_diet,"Brachylabic","POOR diet"),
  summ(delta_b_d4,"Brachylabic","Density 4"),
  summ(delta_b_d8,"Brachylabic","Density 8"),
  summ(delta_M_diet,"Macrolabic","POOR diet"),
  summ(delta_M_d4,"Macrolabic","Density 4"),
  summ(delta_M_d8,"Macrolabic","Density 8")
)

factor(group, levels = c("Macrolabic","Brachylabic"))

figure_5 <- ggplot(delta_morph_df,
       aes(x = median, y = group)) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept=c(-0.1,0.1),
             linetype="dotted", colour="grey50") +
  geom_errorbarh(aes(xmin=l95, xmax=u95),
                 height=.15, size=.9) +
  geom_point(size=3) +
  facet_wrap(~environment) +
  coord_cartesian(xlim = c(-0.7, 0.7)) +
  labs(x=expression(Delta*beta), y=NULL) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(size = 14, face="bold"),
    strip.text.x = element_text(size = 14, face="bold"),
    axis.title.x = element_text(size = 16, face="bold"),
    axis.title.y = element_text(size = 16, face="bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    ggside.axis.text.x = element_blank(),
    ggside.axis.ticks.x = element_blank(),
    legend.position = c(0.1, 0.75)
  )

ggsave(figure_5, filename="figure_5.jpg", width=12.83, height=8.83, dpi=300,antialias="default")




