# ---------------------------
# Experiment 1: Stats in R
# ---------------------------

# Packages
library(tidyverse)
library(lme4)
library(lmerTest)   # p-values for lmer
library(emmeans)    # marginal means + contrasts
library(broom.mixed)
library(Matrix)

# ---- 1) Load data ----
df <- read.csv("../data/exp1_filewise_predictions_long.csv", stringsAsFactors = FALSE)
# ---- 2) OPTIONAL: rename columns if needed ----
# Uncomment/edit if your columns differ.
# df <- df %>%
#   rename(
#     y_true = label,        # or Group -> make label separately
#     y_prob = prob,         # predicted probability for DEM class
#     y_pred = pred          # predicted class (0/1)
#   )

# ---- 3) Basic cleaning ----
df <- df %>%
  mutate(
    Filename = as.factor(Filename),
    Source   = as.factor(Source),
    Task     = as.factor(Task),
    Group    = as.factor(Group)
  )

# Ensure Group is coded as 0/1 for modeling (DEM = 1, OLD = 0)
# Adjust if your labels differ.
df <- df %>%
  mutate(
    y_true = case_when(
      Group %in% c("DEM", "Dementia", "dementia") ~ 1L,
      Group %in% c("OLD", "Healthy", "healthy", "control") ~ 0L,
      TRUE ~ NA_integer_
    )
  )

stopifnot(sum(is.na(df$y_true)) == 0)

# Set reference levels (important for interpretation)
df <- df %>%
  mutate(
    Source = relevel(Source, ref = "Human")
  )

# ---- 4) Create outcome variable(s) ----
# Preferred: use predicted probability if available
has_prob <- "y_prob" %in% names(df)

# If probability is not present but y_pred exists, use correctness model
if (!has_prob) {
  if (!("y_pred" %in% names(df))) stop("Need either y_prob or y_pred in the CSV.")
  df <- df %>%
    mutate(correct = as.integer(y_pred == y_true))
} else {
  # Ensure y_prob is numeric and within (0,1)
  df$y_prob <- as.numeric(df$y_prob)
  if (any(df$y_prob < 0 | df$y_prob > 1, na.rm = TRUE)) stop("y_prob must be between 0 and 1.")
  # logit transform to make linear mixed model reasonable
  eps <- 1e-6
  df <- df %>%
    mutate(
      y_prob_clip = pmin(pmax(y_prob, eps), 1 - eps),
      logit_prob  = qlogis(y_prob_clip)
    )
}

# ---------------------------
# 5) Main model(s)
# ---------------------------

# (A) If you have probabilities -> linear mixed model on logit(prob)
if (has_prob) {
  m1 <- lmer(logit_prob ~ Source + Task + (1 | Filename), data = df)
  m2 <- lmer(logit_prob ~ Source * Task + (1 | Filename), data = df)  # interaction check
  
  cat("\n--- Mixed model: logit(prob) ~ Source + Task + (1|Filename) ---\n")
  print(summary(m1))
  cat("\n--- ANOVA (Type III) for main effects ---\n")
  print(anova(m1))
  
  cat("\n--- Interaction check: Source*Task ---\n")
  print(anova(m1, m2))  # likelihood ratio test
  
  # Estimated marginal means (back-transformed to probability)
  emm_src <- emmeans(m1, ~ Source, type = "response") # type=response returns prob scale
  print(emm_src)
  
  # Pairwise comparisons: each ASR vs Human, with FDR adjustment
  comp_vs_human <- contrast(emm_src, method = "trt.vs.ctrl", ref = "Human", adjust = "fdr")
  comp_vs_human_df <- as.data.frame(comp_vs_human)
  print(comp_vs_human_df)
  
  # Optional: all pairwise Source comparisons
  # all_pairs <- pairs(emm_src, adjust="fdr")
  
  # Task effect summary (marginal means per Task)
  emm_task <- emmeans(m1, ~ Task, type = "response")
  print(emm_task)
  
  # Export tables
  write.csv(as.data.frame(emm_src), "exp1_emm_source.csv", row.names = FALSE)
  write.csv(comp_vs_human_df, "exp1_source_vs_human_fdr.csv", row.names = FALSE)
  write.csv(as.data.frame(emm_task), "exp1_emm_task.csv", row.names = FALSE)
  
} else {
  
  # (B) If you only have predicted class -> mixed-effects logistic regression on correctness
  m1 <- glmer(correct ~ Source + Task + (1 | Filename), data = df, family = binomial)
  m2 <- glmer(correct ~ Source * Task + (1 | Filename), data = df, family = binomial)
  
  cat("\n--- Mixed model: correct ~ Source + Task + (1|Filename) ---\n")
  print(summary(m1))
  cat("\n--- Likelihood ratio test for interaction: Source*Task ---\n")
  print(anova(m1, m2))
  
  emm_src <- emmeans(m1, ~ Source, type = "response")
  comp_vs_human <- contrast(emm_src, method = "trt.vs.ctrl", ref = "Human", adjust = "fdr")
  comp_vs_human_df <- as.data.frame(comp_vs_human)
  
  emm_task <- emmeans(m1, ~ Task, type = "response")
  
  write.csv(as.data.frame(emm_src), "exp1_emm_source.csv", row.names = FALSE)
  write.csv(comp_vs_human_df, "exp1_source_vs_human_fdr.csv", row.names = FALSE)
  write.csv(as.data.frame(emm_task), "exp1_emm_task.csv", row.names = FALSE)
}

# ---------------------------
# 6) Optional: quick plots
# ---------------------------
# Visualise marginal mean performance by Source (probability scale)

# If you ran the probability model:
# plot(emm_src) + ggtitle("Experiment 1: Estimated mean P(DEM) by Source")