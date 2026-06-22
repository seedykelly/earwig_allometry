# ==============================================================================
# RESEARCH PROJECT: 2026 Earwig Allometry
#
# DESCRIPTION:
# This script analyzes weapon allometry across nested biological levels. 
# It examines how sex, morph (major/minor), diet, and density influence 
# scaling plasticity and intercept/slope parameters.
#
# KEY STEPS:
# 1. Data cleaning and preparation
# 2. Bayesian linear and quadratic allometric modeling
# 3. Model comparison (LOO)
# 4. Post-hoc slope analysis and environmental plasticity
# 5. High-quality visualization for publication
# ==============================================================================

## ---- 1. Load Libraries ----

library(tidyverse)  # Core data manipulation and visualization
library(scales)     # For scale formatting (e.g., label_number)
library(broman)     # For emmeans/emtrends integration
library(cowplot)    # For plot arrangement
library(flextable)  # For table generation
library(officer)    # For document interaction
library(knitr)       # For table formatting
library(officedown) # For advanced RMarkdown/Quarto
library(broom.mixed) # For tidying model outputs
library(performance) # For model diagnostics
library(brms)        # Bayesian regression modeling
library(tidybayes)   # For posterior distribution visualization
library(mixtools)    # For mixture modeling (used in scratchpad)
library(mixsmsn)    # Used in scratchpad
library(cmdstanr)   # High-performance Stan backend
library(ggside)     # For marginal distributions in plots
library(posterior) # For posterior sample manipulation
library(emmeans)     # For marginal means and trends
library(patchwork)  # For combining plots
library(ggrepel)    # For non-overlapping labels
library(here)

## ---- 2. Data Loading and Cleaning ----

# Load raw observations
earwig_data_raw <- read.csv(file="data/raw/earwig_allometry.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()

# Clean basic variables and convert to factors
earwig_data <- earwig_data_raw %>%
  mutate(
    sex = recode(sex, femelle = "female", male = "male"),
    id = Id_adulte,
    density = as.factor(density), 
    diet = as.factor(diet), 
    sex = as.factor(sex), 
    id_mere = as.factor(id_mere)
  )

# Filter for complete replicates to ensure statistical consistency:
# density-1 must have 1, density-4 must have 3-4, density-8 must have 6, 7, or 8
earwig_two <- earwig_data %>%
  group_by(boite_petri, diet, density) %>%
  summarise(num_individ = n(), .groups = "drop") %>%
  filter(
    (density == "8" & num_individ >= 6) | 
    (density == "4" & num_individ >= 3) | 
    (density == "1" & num_individ == 1)
  ) %>%
  mutate(valid = 1) %>%
  dplyr::select(boite_petri, valid) 

# Merge valid replicates back and process dates/development time
dat.new <- left_join(earwig_data, earwig_two, by = "boite_petri") %>%
  filter(valid == 1) %>%
  mutate(
    treat_date = dmy(date_treatment), 
    adult_date = dmy(date_adulte),
    dev_time = as.numeric(difftime(adult_date, treat_date, units = "days"))
  ) %>%
  filter(sex != "") %>%
  droplevels()

## ---- 3. Morph Integration ----

# The morph categorization (identifying major/minor males) was performed 
# in the scratchpad: scripts/2026_earwig_morph_analysis_scratchpad.R
load(file = "data/processed/dat.morphs.Rda")

# Prepare allometric variables and standardized factors for modeling and plotting
# We consolidate factor levels and labels here to avoid repetitive 'recode' in plots.
dat.morphs <- dat.morphs %>%
  mutate(
    logF = log(forceps_L),
    logP = log(pronotum),
    
    # Grouping for modeling: Combining sex and morph into one factor
    group3 = case_when(
      sex == "female" ~ "female",
      sex == "male" & morph == "minor" ~ "brachylabic",
      sex == "male" & morph == "major" ~ "macrolabic"
    ),
    
    # Define factor levels and labels once for consistency
    group3  = factor(group3, levels = c("female", "brachylabic", "macrolabic")),
    diet    = factor(diet, levels = c("GOOD", "POOR")),
    density = factor(density, levels = c("1", "4", "8")),
    id_mere = factor(id_mere)
  ) %>%
  filter(!is.na(group3), !is.na(logF), !is.na(logP))

## ---- 4. Bayesian Allometric Modeling ----

# Load pre-fitted models and comparison results
message("Loading pre-fitted models from models/...")
mod_all <- readRDS(here::here("models", "mod_all.rds"))
mod_all_quad <- readRDS(here::here("models", "mod_all_quad.rds"))
loo_table <- readRDS(here::here("models", "loo_table.rds"))

print(loo_table)

## ---- 5. Post-hoc Statistical Analysis ----

# Generate descriptive summary table for Results (used in .qmd)
summary.table.2 <- dat.morphs %>%
  group_by(group3) %>%
  summarise(
    mean_forceps = mean(forceps_L, na.rm = TRUE),
    sd_forceps   = sd(forceps_L, na.rm = TRUE),
    n            = n(),
    se_forceps   = sd_forceps / sqrt(n),
    mean_pronotum = mean(pronotum, na.rm = TRUE),
    sd_pronotum   = sd(pronotum, na.rm = TRUE),
    se_pronotum   = sd_pronotum / sqrt(n),
    .groups = "drop"
  )

# Re-calculate conditional slopes and contrasts from loaded models
# to ensure all subsequent visualization objects are available.

# 5.1 Calculate conditional slopes (emtrends)
slopes_all <- emtrends(mod_all, ~ group3 * diet * density, var = "logP")

# Convert emtrends object to a data frame for the Quarto document
slope_table <- as.data.frame(slopes_all)
# Rename columns to match .qmd expectations
if("estimate" %in% colnames(slope_table)) {
  slope_table <- rename(slope_table, logP.trend = estimate)
}
if("lower.CL" %in% colnames(slope_table)) {
  slope_table <- rename(slope_table, lower.HPD = lower.CL)
}
if("upper.CL" %in% colnames(slope_table)) {
  slope_table <- rename(slope_table, upper.HPD = upper.CL)
}

# 5.2 Group comparisons within environments
group_slope_contrasts <- pairs(emtrends(mod_all, ~ group3 | diet * density, var = "logP"))

# 5.3 Environmental plasticity (slope differences)
diet_slope_contrasts <- pairs(emtrends(mod_all, ~ diet | group3 * density, var = "logP"))
density_slope_contrasts <- pairs(emtrends(mod_all, ~ density | group3 * diet, var = "logP"))

## ---- 6. Visualizations ----

# Pre-format data for all plots to ensure uniform labeling and aesthetics
plot_data_master <- dat.morphs %>%
  mutate(
    # Readable labels for facets and legends
    Phenotype = factor(recode(group3, 
                               female = "Female", 
                               brachylabic = "Brachylabic male", 
                               macrolabic = "Macrolabic male"),
                       levels = c("Female", "Brachylabic male", "Macrolabic male")),
    Diet_lab = factor(recode(diet, GOOD = "Good", POOR = "Poor"), levels = c("Poor", "Good")),
    Density_lab = factor(recode(density, "1" = "Low", "4" = "Medium", "8" = "High"), 
                         levels = c("Low", "Medium", "High"))
  )

# #### Figure 2: Observed Data and Marginal Distributions ####

figure_2 <- ggplot(plot_data_master, aes(x = pronotum, y = forceps_L, colour = Diet_lab)) +
  geom_point(size = 3, alpha = 0.3, stroke = NA) +
  facet_grid(group3 ~ density, 
             labeller = labeller(density = c("1"="Low", "4"="Medium", "8"="High"), 
                                 group3 = c("female"="Females", "brachylabic"="Brachylabic males", "macrolabic"="Macrolabic males"))) +
  geom_ysidedensity(aes(x = after_stat(density), group = group3), colour = "black") +
  ggside(collapse = "y") +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_colour_manual("Diet", values = c("black", "red"), labels = c("Poor", "Good")) +
  labs(x = "Pronotum length (mm)", y = "Forceps length (mm)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    strip.background = element_rect(fill = "white"),
    ggside.axis.text.x = element_blank(),
    ggside.axis.ticks.x = element_blank()
  )

ggsave(figure_2, filename = "Figure_2.jpg", width = 14.83, height = 8.83, dpi = 300)

# #### Figure 3: Posterior Slope Estimates ####

slope_df <- as.data.frame(slopes_all) %>%
  mutate(
    Phenotype = factor(recode(group3, female = "Female", brachylabic = "Brachylabic male", macrolabic = "Macrolabic male"),
                       levels = c("Female", "Brachylabic male", "Macrolabic male")),
    Diet_lab = factor(recode(diet, GOOD = "Good", POOR = "Poor"), levels = c("Good", "Poor")),
    Density_lab = factor(recode(density, "1" = "Low", "4" = "Medium", "8" = "High"), 
                         levels = c("Low", "Medium", "High"))
  )

figure_3 <- ggplot(slope_df, aes(x = logP.trend, y = Density_lab, colour = Diet_lab)) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.6, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lower.HPD, xmax = upper.HPD), width = 0.2, linewidth = 0.6, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Phenotype, nrow = 1) +
  scale_color_manual(values = c("Good" = "#D55E00", "Poor" = "#0072B2")) +
  labs(x = expression(bold("Allometric slope (" * beta * ")")), y = "Rearing density", colour = "Diet") +
  theme_classic() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

ggsave(filename = "figure_3.jpg", plot = figure_3, width = 12.83, height = 8.83, dpi = 300)

# #### Figure S1: Predicted Allometries (Population-Level) ####

# Generate prediction grid for population-level averages
prediction_grid_s1 <- dat.morphs %>%
  group_by(group3, diet, density) %>%
  summarise(min_logP = min(logP), max_logP = max(logP), .groups = "drop") %>%
  mutate(logP = map2(min_logP, max_logP, ~seq(.x, .y, length.out = 50))) %>%
  unnest(logP) %>%
  mutate(id_mere = NA, boite_petri = NA)

# Get posterior predictions (excluding random effects)
epred_s1 <- posterior_epred(mod_all, newdata = prediction_grid_s1, re_formula = NA)

new_all_s1 <- prediction_grid_s1 %>%
  mutate(
    fit = apply(epred_s1, 2, median),
    lwr = apply(epred_s1, 2, quantile, 0.025),
    upr = apply(epred_s1, 2, quantile, 0.975),
    Phenotype = factor(recode(group3, female="Female", brachylabic="Brachylabic", macrolabic="Macrolabic"),
                       levels = c("Female", "Brachylabic", "Macrolabic"))
  )

figure_S1 <- ggplot(new_all_s1, aes(x = logP, y = fit, colour = diet, linetype = density, group = interaction(diet, density))) +
  geom_point(data = dat.morphs, aes(x = logP, y = logF), inherit.aes = FALSE, colour = "grey50", alpha = 0.12, size = 0.7) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ Phenotype, nrow = 1) +
  labs(x = "log(Pronotum length, mm)", y = "log(Forceps length, mm)", colour = "Diet", linetype = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 12))

ggsave("figure_S1.jpg", plot = figure_S1, width = 10, height = 6, dpi = 300)

# #### Figure S2: Linear vs. Quadratic Comparison ####

# Helper to extract population-level predictions
get_pop_preds <- function(fit, newdata, label) {
  epred <- posterior_epred(fit, newdata = newdata, re_formula = NA)
  newdata %>%
    mutate(
      logF_med = apply(epred, 2, median),
      logF_low = apply(epred, 2, quantile, probs = 0.025),
      logF_high = apply(epred, 2, quantile, probs = 0.975),
      Model = label
    )
}

# Create common prediction grid
prediction_grid_s2 <- dat.morphs %>%
  group_by(group3, diet, density) %>%
  summarise(min_logP = min(logP), max_logP = max(logP), .groups = "drop") %>%
  mutate(logP = map2(min_logP, max_logP, ~seq(.x, .y, length.out = 100))) %>%
  unnest(logP) %>%
  mutate(id_mere = dat.morphs$id_mere[1], boite_petri = dat.morphs$boite_petri[1])

set.seed(123)
lin_preds <- get_pop_preds(mod_all, prediction_grid_s2, "Linear")
set.seed(123)
quad_preds <- get_pop_preds(mod_all_quad, prediction_grid_s2, "Quadratic")

# Combine and back-transform for plotting
plot_data_s2 <- bind_rows(lin_preds, quad_preds) %>%
  mutate(
    Phenotype = factor(recode(as.character(group3), female="Female", brachylabic="Brachylabic male", macrolabic="Macrolabic male"),
                       levels = c("Female", "Brachylabic male", "Macrolabic male")),
    Diet_lab = factor(recode(as.character(diet), GOOD="Good", POOR="Poor"), levels = c("Good", "Poor")),
    Density_lab = factor(recode(as.character(density), "1"="Low", "4"="Medium", "8"="High"), 
                       levels = c("Low", "Medium", "High")),
    pronotum_pred = exp(logP),
    forceps_med = exp(logF_med),
    forceps_low = exp(logF_low),
    forceps_high = exp(logF_high),
    Model = factor(Model, levels = c("Linear", "Quadratic"))
  )

obs_data_s2 <- dat.morphs %>%
  mutate(
    Phenotype = factor(recode(group3, female="Female", brachylabic="Brachylabic male", macrolabic="Macrolabic male"),
                       levels = c("Female", "Brachylabic male", "Macrolabic male")),
    Diet_lab = factor(recode(as.character(diet), GOOD="Good", POOR="Poor"), levels = c("Good", "Poor")),
    Density_lab = factor(recode(as.character(density), "1"="Low", "4"="Medium", "8"="High"), 
                       levels = c("Low", "Medium", "High")),
    pronotum_obs = exp(logP),
    forceps_obs = exp(logF)
  )

figure_s_quadratic <- ggplot() +
  geom_point(data = obs_data_s2, aes(x = pronotum_obs, y = forceps_obs, colour = Diet_lab), alpha = 0.18, size = 0.7) +
  geom_ribbon(data = filter(plot_data_s2, Model == "Quadratic"), 
              aes(x = pronotum_pred, ymin = forceps_low, ymax = forceps_high, fill = Diet_lab), alpha = 0.18, colour = NA) +
  geom_line(data = filter(plot_data_s2, Model == "Linear"), 
            aes(x = pronotum_pred, y = forceps_med, colour = Diet_lab, linetype = Model), colour = "grey30", linewidth = 0.65) +
  geom_line(data = filter(plot_data_s2, Model == "Quadratic"), 
            aes(x = pronotum_pred, y = forceps_med, colour = Diet_lab, linetype = Model), linewidth = 1) +
  facet_grid(rows = vars(Phenotype), cols = vars(Density_lab), scales = "free_y") +
  scale_x_log10() + scale_y_log10() +
  scale_linetype_manual(name = "Model", values = c("Linear" = "dotdash", "Quadratic" = "solid"),
                       labels = c("Linear model", "Quadratic model")) +
  labs(x = "Pronotum length (mm)", y = "Forceps length (mm)", colour = "Diet") +
  theme_classic() +
  theme(legend.position = "top", strip.text = element_text(face = "bold"), axis.title = element_text(face = "bold"))

ggsave("figure_S2.jpg", plot = figure_s_quadratic, width = 10, height = 6, dpi = 300)

# ---- end ----
