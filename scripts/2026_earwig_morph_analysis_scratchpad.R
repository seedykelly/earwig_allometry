# ==============================================================================
# SCRATCHPAD: Morph Analysis (One-time execution)
# Purpose: Perform bimodal analysis on male forceps length to identify morphs 
#          (minor vs. major) and combine with female data.
# ==============================================================================

# This script contains code that is only required once to generate the 
# processed morph data. Once 'data/processed/dat.morphs.Rda' is created, 
# this section is no longer needed for the main analysis pipeline.

library(dplyr)
library(mixsmsn)

# 1. Bimodal analysis for Males
# We use the smsn.mix function to identify two distinct modes in the forceps length distribution.
males <- dat.new %>%
  filter(sex == "male")

bimodal.analysis <- smsn.mix(males$forceps_L, nu=3, g=2, get.init = TRUE, criteria = TRUE,
                             group = TRUE, family = "Skew.normal", iter.max = 1000, calc.im=TRUE, obs.prob = TRUE)

# mixsmsn::mix.print(bimodal.analysis, digits=3)
# saveRDS(bimodal.analysis, file ="data/processed/bimodal.analysis.rda")

# mix.hist(males$forceps_L, bimodal.analysis)

# 2. Assign morphs to Males
males$group <- bimodal.analysis$group
males <- males %>% mutate(
  morph = case_when(
    group == 2 ~ "minor",
    group == 1 ~ "major"
  )
)
males$group <- as.factor(males$group)

# save(males, file="data/processed/males.Rda")

# 3. Prepare Females
females <- dat.new %>%
  filter(sex == "female") %>%
  mutate(group = 3, morph = "female")

females$group <- as.factor(females$group)

# 4. Combine Datasets
dat.morphs <- rbind(males, females)

# 5. Manual Correction
# Note: One individual was manually corrected based on visual inspection of the distribution.
dat.morphs <- dat.morphs %>%
  mutate(group = replace(group, id == "DJODI", "2"))

# Save the finalized morph-categorized dataset
save(dat.morphs, file = "data/processed/dat.morphs.Rda")
