# ==============================================================================
# RESEARCH PROJECT: 2026 Earwig Allometry
#
# DESCRIPTION:
# This script handles the heavy lifting of Bayesian model fitting.
# It loads the necessary data, performs required transformations,
# defines the models, and saves the results to the 'models/' directory.
# ==============================================================================

## ---- 1. Load Libraries ----

library(tidyverse)
library(brms)
library(cmdstanr)
library(here)

## ---- 2. Load and Prepare Data ----

# Load the processed morph data
data_path <- here::here("data", "processed", "dat.morphs.Rda")

# Handle potential naming discrepancies in .Rda files
temp_env <- new.env()
load(data_path, envir = temp_env)
actual_object_name <- ls(temp_env)[1]
dat.morphs <- temp_env[[actual_object_name]]

message(paste0("Successfully loaded data object '", actual_object_name, "' as 'dat.morphs'"))

# Perform necessary transformations required for the models
# Replicating logic from the main analysis script
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
    
    # Ensure factors are set correctly for brms
    group3  = factor(group3, levels = c("female", "brachylabic", "macrolabic")),
    diet    = factor(diet, levels = c("GOOD", "POOR")),
    density = factor(density, levels = c("1", "4", "8")),
    id_mere = factor(id_mere)
  ) %>%
  filter(!is.na(group3), !is.na(logF), !is.na(logP))

message("Data transformations complete.")

## ---- 3. Model Definitions ----

# Bayesian Linear Allometry Model
formula_all <- bf(
  logF ~ logP * group3 * diet * density +
    (1 + logP | id_mere) + (1 | boite_petri)
)

# Informative priors
priors_all <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(student_t(3, 0, 1), class = "sigma"),
  prior(student_t(3, 0, 1), class = "sd"),
  prior(lkj(2), class = "cor")
)

# Quadratic version of the model
formula_all_quad <- bf(
  logF ~ logP * group3 * diet * density +
    I(logP^2) * group3 * diet * density +
    (1 + logP | id_mere) + (1 | boite_petri)
)

## ---- 4. Fit Models ----

# Set seed for reproducibility
set.seed(4821)

# 4.1 Fit the linear model
message("Fitting linear model (mod_all)...")
mod_all <- brm(
  formula = formula_all,
  data    = dat.morphs,
  family  = gaussian(),
  prior   = priors_all,
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97)
)

# 4.2 Fit the quadratic model
message("Fitting quadratic model (mod_all_quad)...")
mod_all_quad <- brm(
  formula = formula_all_quad,
  data    = dat.morphs,
  family  = gaussian(),
  prior   = priors_all,
  backend = "cmdstanr",
  chains  = 4,
  cores   = 4,
  iter    = 4000,
  control = list(adapt_delta = 0.97)
)

## ---- 5. Model Comparison ----

message("Performing LOO comparison...")
mod_all <- add_criterion(mod_all, "loo")
mod_all_quad <- add_criterion(mod_all_quad, "loo")

loo_table <- as.data.frame(loo_compare(mod_all, mod_all_quad)) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(Model = recode(Model, mod_all = "Linear model", mod_all_quad = "Quadratic model"))

print(loo_table)

## ---- 6. Save Outputs ----

message("Saving models and comparison table...")

saveRDS(mod_all, here::here("models", "mod_all.rds"))
saveRDS(mod_all_quad, here::here("models", "mod_all_quad.rds"))
saveRDS(loo_table, here::here("models", "loo_table.rds"))

message("All models fitted and saved successfully.")
