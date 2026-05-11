

df <- read.csv("../data/exp3_filewise_predictions_long.csv")

# ---------------------------
# Experiment 3: Biomarker distortion stats
# ---------------------------

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)

# 1) Load your file-wise biomarker table
#df <- read.csv("exp3_filewise_biomarkers.csv", stringsAsFactors = FALSE)

# 2) Clean and set factor levels
df <- df %>%
  mutate(
    Filename = factor(Filename),
    Source   = factor(Source),
    Group    = factor(Group)
  )

# Ensure Human is the reference transcription source
df <- df %>% mutate(Source = relevel(Source, ref = "Human"))

# Ensure Control/Healthy is the reference group (OLD = 0 baseline)
# Your labels: DEM and OLD
df <- df %>% mutate(Group = relevel(Group, ref = "OLD"))

# Optional: check for missing values
stopifnot(sum(is.na(df$ttr)) == 0,
          sum(is.na(df$pronoun_ratio)) == 0,
          sum(is.na(df$content_ratio)) == 0,
          sum(is.na(df$mean_word_len)) == 0)

# Helper function to run model + exports for one biomarker
run_biomarker_model <- function(data, outcome, out_prefix) {
  
  form <- as.formula(paste0(outcome, " ~ Group * Source + (1|Filename)"))
  
  m <- lmer(form, data = data, REML = TRUE)
  
  cat("\n============================\n")
  cat("Outcome:", outcome, "\n")
  cat("============================\n")
  print(summary(m))
  cat("\n--- Type III ANOVA ---\n")
  print(anova(m))  # includes Group, Source, Group:Source
  
  # Estimated dementia effect within each Source:
  # (DEM - OLD) computed per source, then compare to Human
  emm <- emmeans(m, ~ Group | Source)
  
  dem_effect_by_source <- contrast(emm, method = "revpairwise") %>%
    as.data.frame()
  # 'revpairwise' gives DEM - OLD if Group levels are OLD, DEM (check output).
  # If sign is flipped, swap to "pairwise" or relevel.
  
  # Compare dementia effect (DEM-OLD) in each ASR vs Human (interaction interpretation)
  # This is the core "distortion" test.
  # Build a table of DEM-OLD per Source first, then compare those effects vs Human.
  dem_effect_emm <- emmeans(m, ~ Group * Source)
  dem_effects <- contrast(dem_effect_emm, interaction = "pairwise")  # all interactions
  
  # A simpler and clearer approach:
  # Get dementia effect per source: (DEM-OLD) at each Source
  dem_eff <- contrast(emmeans(m, ~ Group | Source), method = "pairwise")  # DEM-OLD per Source
  dem_eff_df <- as.data.frame(dem_eff)
  
  # Compare those dementia effects to Human: ASR vs Human (FDR)
  # We do this by fitting an emmeans on the "Group difference" and contrasting.
  # Use emtrends-like approach via 'contrast' twice:
  dem_eff_emm <- emmeans(m, ~ Source, by = NULL,
                         weights = "equal",
                         # This 'Source' EMM isn't the group effect; we use the computed diffs instead.
  )
  # Instead: directly contrast DEM-OLD effects using 'dem_eff' object
  # dem_eff is already per Source; we can contrast its estimates.
  #dem_eff_vs_human <- contrast(dem_eff, method = "trt.vs.ctrl", ref = "Human", adjust = "fdr")
  dem_eff_vs_human <- contrast(dem_eff, method = "trt.vs.ctrl", ref = 1, adjust = "fdr")
  dem_eff_vs_human_df <- as.data.frame(dem_eff_vs_human)
  
  # Save outputs
  write.csv(tidy(m, effects = "fixed"), paste0(out_prefix, "_fixed_effects.csv"), row.names = FALSE)
  write.csv(as.data.frame(anova(m)), paste0(out_prefix, "_anova.csv"))
  write.csv(dem_eff_df, paste0(out_prefix, "_dem_effect_by_source.csv"), row.names = FALSE)
  write.csv(dem_eff_vs_human_df, paste0(out_prefix, "_dem_effect_vs_human_fdr.csv"), row.names = FALSE)
  
  invisible(list(model = m,
                 dem_effect_by_source = dem_eff_df,
                 dem_effect_vs_human = dem_eff_vs_human_df))
}

# 3) Run for each biomarker
res_ttr   <- run_biomarker_model(df, "ttr",           "exp3_ttr")
res_pro   <- run_biomarker_model(df, "pronoun_ratio", "exp3_pronoun_ratio")
res_cont  <- run_biomarker_model(df, "content_ratio", "exp3_content_ratio")
res_wlen  <- run_biomarker_model(df, "mean_word_len", "exp3_mean_word_len")

# 4) Optional: sanity check distribution of word counts by Source
# (Useful as a covariate if needed)
summary(df$n_words)