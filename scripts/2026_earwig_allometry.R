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
dat.morphs %>%
  filter(morph=="CCNJO")

# # earwig_data_complete %>%
# #   filter(id_mere=="MNLDO")
# 
# # save(dat.morphs,file="data/processed/dat.morphs.Rda")

# ============================
# Allometry
# ============================

# ============================
# sex difference
# ============================

# Transform variables
dat.morphs <- dat.morphs %>%
  mutate(
    logF  = log(forceps_L),
    logP  = log(pronotum),
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
  logF ~ logP * sex * diet * density +
    (1 + logP | id_mere)
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

#saveRDS(mod_sex, file = "mod_sex.Rds")
mod_sex <- readRDS(file = "data/processed/mod_sex.Rds")

print(tidy(mod_sex), n = "all")

slopes_sex <- emtrends(
  mod_sex,
  ~ sex * diet * density,
  var = "logP"
)

summary(slopes_sex)

male_vs_female <- pairs(
  emtrends(mod_sex, ~ sex | diet * density, var = "logP")
)

summary(male_vs_female)


# =============================
# test for non-linearity
# =============================

formula_sex_quad <- bf(
  logF ~ logP * sex * diet * density +
    I(logP^2) * sex * diet * density +
    (1 + logP | id_mere)
)

# mod_sex_quad <- brm(
#   formula = formula_sex_quad,
#   data    = dat.morphs,
#   family  = gaussian(),
#   prior   = priors_sex,
#   chains  = 4,
#   cores   = 4,
#   iter    = 4000,
#   backend = "cmdstanr",
#   file    = "data/processed/mod_sex_quad.Rds",
#   control = list(adapt_delta = 0.97)
# )

mod_sex_quad <- readRDS(file = "data/processed/mod_sex_quad.Rds")

summary(mod_sex_quad)

loo.comp <- loo(mod_sex, mod_sex_quad)


# ============================
# males only
# ============================
males <- subset(dat.morphs, sex == "male")
males$morph <- factor(males$morph)
males$morph <- relevel(males$morph, ref = "minor")  # brachylabic as reference

formula_morph <- bf(
  logF ~ logP * morph * diet * density +
    (1 + logP | id_mere)
)

priors_morph <- c(
  
  # Population-level fixed effects
  prior(normal(0, 0.5), class = "b"),
  
  # Intercept (log-scale forceps size)
  prior(normal(1.2, 0.5), class = "Intercept"),
  
  # Residual SD
  prior(student_t(3, 0, 0.5), class = "sigma"),
  
  # Random effect SDs
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
#   file    = "data/processed/mod_morph.Rds",
#   control = list(adapt_delta = 0.97)
# )

# saveRDS(mod_morph, file = "mod_morph.Rds")

mod_morph <- readRDS(file = "data/processed/mod_morph.Rds")

print(tidy(mod_morph), n = "all")

draws <- as_draws_df(mod_morph)

slopes_morph <- emtrends(
  mod_morph,
  ~ morph * diet * density,
  var = "logP"
)

summary(slopes_morph)

pairs(
  emtrends(mod_morph, ~ morph | diet * density, var = "logP")
)

pairs(
  emtrends(mod_morph, ~ diet | morph * density, var = "logP")
)


pairs(
  emtrends(mod_morph, ~ density | morph * diet, var = "logP")
)


## ---- end

# ==============================
# PLOTS
# ==============================

library(patchwork)

#### figure 2 ####

dat.morphs.2 <- dat.morphs %>% mutate(
  morph = case_when(
    group == "minor" ~ "brachylabic",
    group == "major" ~ "macrolabic",
    group == "female" ~ "female")
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

#### figure 3 ####


dat.morphs <- dat.morphs %>%
  mutate(
    group = case_when(
      sex == "female" ~ "female",
      sex == "male" & morph == "minor" ~ "minor",
      sex == "male" & morph == "major" ~ "major"
    )
  )

dat.morphs$group <- factor(dat.morphs$group,
                           levels = c("female", "minor", "major"))

formula_all <- bf(
  logF ~ logP * group * diet * density +
    (1 + logP | id_mere)
)

mod_all <- brm(
  formula = formula_all,
  data    = dat.morphs,
  family  = gaussian(),
  prior   = priors_sex,   # you can reuse or tweak
  backend = "cmdstanr",
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  file    = "data/processed/mod_all.Rds",
  control = list(adapt_delta = 0.97)
)

# -----------------------------
# 1. SEQUENCE FOR PREDICTION
# -----------------------------
x_seq <- seq(
  min(dat.morphs$logP, na.rm = TRUE),
  max(dat.morphs$logP, na.rm = TRUE),
  length.out = 100
)

# -----------------------------
# 2. NEW DATA (ALL GROUPS)
# -----------------------------
new_all <- expand.grid(
  logP    = x_seq,
  group   = c("female", "minor", "major"),
  diet    = c("GOOD", "POOR"),
  density = c("1", "4", "8"),
  id_mere = NA
)

# Ensure factor levels match model
new_all$group   <- factor(new_all$group, levels = levels(dat.morphs$group))
new_all$diet    <- factor(new_all$diet, levels = levels(dat.morphs$diet))
new_all$density <- factor(new_all$density, levels = levels(dat.morphs$density))

# -----------------------------
# 3. MODEL PREDICTIONS
# -----------------------------
epred_all <- posterior_epred(
  mod_all,
  newdata = new_all,
  re_formula = NA
)

# -----------------------------
# 4. SUMMARISE POSTERIOR
# -----------------------------
new_all$fit <- apply(epred_all, 2, median)
new_all$lwr <- apply(epred_all, 2, quantile, 0.025)
new_all$upr <- apply(epred_all, 2, quantile, 0.975)

# -----------------------------
# 5. PANEL LABELS
# -----------------------------
new_all$panel <- case_when(
  new_all$group == "female" ~ "Female",
  new_all$group == "minor"  ~ "Brachylabic male",
  new_all$group == "major"  ~ "Macrolabic male"
)

new_all$panel <- factor(
  new_all$panel,
  levels = c("Female", "Brachylabic male", "Macrolabic male")
)

# -----------------------------
# 6. RAW DATA
# -----------------------------
raw_data <- dat.morphs %>%
  mutate(
    panel = case_when(
      sex == "female" ~ "Female",
      sex == "male" & morph == "minor" ~ "Brachylabic male",
      sex == "male" & morph == "major" ~ "Macrolabic male"
    )
  )

raw_data$panel <- factor(
  raw_data$panel,
  levels = c("Female", "Brachylabic male", "Macrolabic male")
)

# -----------------------------
# 7. FINAL PLOT
# -----------------------------
figure_3 <- ggplot(new_all,
                   aes(x = logP,
                       y = fit,
                       color = diet,
                       fill = diet,
                       linetype = density,
                       group = interaction(diet, density))) +
  
  # prediction lines
  geom_line(linewidth = 1.2) +
  
  # raw data
  geom_point(
    data = raw_data,
    aes(x = logP, y = logF),
    inherit.aes = FALSE,
    color = "black",
    alpha = 0.08,
    size = 0.7
  ) +
  
  facet_wrap(~ panel, nrow = 1, scales = "fixed") +
  
  scale_color_manual(values = c(
    "GOOD" = "#D55E00",
    "POOR" = "#0072B2"
  )) +
  
  scale_fill_manual(values = c(
    "GOOD" = "#D55E00",
    "POOR" = "#0072B2"
  )) +
  
  scale_linetype_manual(values = c(
    "1" = "solid",
    "4" = "dashed",
    "8" = "dotted"
  )) +
  
  labs(
    x = "log(Pronotum length, mm)",
    y = "log(Forceps length, mm)",
    color = "Diet",
    fill = "Diet",
    linetype = "Density"
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    strip.background = element_rect(fill = "white"),
    legend.position = c(0.1, 0.75)
  )


figure_3
ggsave(figure_3, filename="figure_3.jpg", width=12.83, height=8.83, dpi=300,antialias="default")


#### figure 4 ####

slopes_all <- emtrends(
  mod_all,
  ~ group * diet * density,
  var = "logP"
)

slope_df <- as.data.frame(slopes_all)

slope_df <- slope_df %>%
  mutate(
    group = recode(group,
                   "female" = "Female",
                   "minor"  = "Brachylabic male",
                   "major"  = "Macrolabic male"),
    
    diet = recode(diet,
                  "GOOD" = "Good",
                  "POOR" = "Poor"),
    
    density = recode(density,
                     "1" = "Low",
                     "4" = "Medium",
                     "8" = "High"),
    
    density = factor(density, levels = c("Low", "Medium", "High")),
    
    group = factor(group,
                   levels = c("Female", "Brachylabic male", "Macrolabic male"))
  )

figure_4 <- ggplot(slope_df,
                   aes(x = logP.trend,
                       y = density,
                       color = diet)) +
  
  geom_vline(xintercept = 0,
             linetype = "dotted",
             linewidth = 0.6,
             alpha = 0.4) +
  
  geom_errorbarh(aes(xmin = lower.HPD, xmax = upper.HPD),
                 width = 0.2,
                 linewidth = 0.6,
                 position = position_dodge(width = 0.5)) +
  
  geom_point(size = 2.5,
             position = position_dodge(width = 0.5)) +
  
  facet_wrap(~ group, nrow = 1) +
  
  scale_color_manual(values = c(
    "Good" = "#D55E00",
    "Poor" = "#0072B2"
  )) +
  
  labs(
    x = expression(paste("Allometric slope (", beta, ")")),
    y = "Rearing density",
    color = "Diet"
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

ggsave(figure_4, filename="figure_4.jpg", width=12.83, height=8.83, dpi=300,antialias="default")

