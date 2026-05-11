
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)

OUT <- "results/exp2/"

# Load data
df <- read.csv("../data/exp2_fold_source_error_rates.csv")

# Convert variables
df$Source <- factor(df$Source)
df$fold <- factor(df$fold)

# Set Human as reference
df$Source <- relevel(df$Source, ref = "Human")

# --------------------------
# Model 1: False Negative Rate
# --------------------------

model_fnr <- lmer(FNR ~ Source + (1|fold), data=df)

summary(model_fnr)
anova(model_fnr)

# Pairwise comparisons vs Human
fnr_contrast <- emmeans(model_fnr, pairwise ~ Source, adjust="fdr")
fnr_contrast$contrasts

# Export FNR results
write.csv(tidy(model_fnr, effects = "fixed"),
          paste0(OUT, "Exp_2_fnr_fixed_effects.csv"), row.names = FALSE)
write.csv(as.data.frame(anova(model_fnr)),
          paste0(OUT, "Exp_2_fnr_anova.csv"), row.names = TRUE)
write.csv(as.data.frame(fnr_contrast$contrasts),
          paste0(OUT, "Exp_2_fnr_contrasts_fdr.csv"), row.names = FALSE)
write.csv(as.data.frame(emmeans(model_fnr, ~ Source)),
          paste0(OUT, "Exp_2_fnr_emm.csv"), row.names = FALSE)


# --------------------------
# Model 2: False Positive Rate
# --------------------------

model_fpr <- lmer(FPR ~ Source + (1|fold), data=df)

summary(model_fpr)
anova(model_fpr)

# Pairwise comparisons vs Human
fpr_contrast <- emmeans(model_fpr, pairwise ~ Source, adjust="fdr")
fpr_contrast$contrasts

# Export FPR results
write.csv(tidy(model_fpr, effects = "fixed"),
          paste0(OUT, "Exp_2_fpr_fixed_effects.csv"), row.names = FALSE)
write.csv(as.data.frame(anova(model_fpr)),
          paste0(OUT, "Exp_2_fpr_anova.csv"), row.names = TRUE)
write.csv(as.data.frame(fpr_contrast$contrasts),
          paste0(OUT, "Exp_2_fpr_contrasts_fdr.csv"), row.names = FALSE)
write.csv(as.data.frame(emmeans(model_fpr, ~ Source)),
          paste0(OUT, "Exp_2_fpr_emm.csv"), row.names = FALSE)

message("Exp 2 results saved to: ", OUT)