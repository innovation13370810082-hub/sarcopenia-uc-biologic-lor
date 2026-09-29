# ============================================================
# Model Increment Analysis
#
# This script compares a base Cox model with an extended Cox model
# containing one additional exposure variable. Variable names have been
# anonymized for public release.
#
# Required anonymized columns:
#   follow_up_time, event, exposure, and covariate_01 to covariate_12
# ============================================================

library(survival)
library(boot)

input_file <- file.path("data", "analysis_data.csv")
figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. Read the analysis data
# ============================================================

data <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 2. Define categorical variables
# ============================================================

factor_vars <- c(
  "covariate_XX",
  "exposure"
)

data[factor_vars] <- lapply(
  data[factor_vars],
  as.factor
)


# ============================================================
# 3. Base clinical model
#
# Contains only the base covariates
# Does not contain the additional exposure variable
# ============================================================

predictors_base <- c(
  "covariate_XX"
)


# ============================================================
# 4. Extended model
#
# Adds only the anonymized exposure variable
# ============================================================

predictors_extended <- c(
  predictors_base,
  "exposure"
)


# ============================================================
# 5. Complete-case dataset
#
# Both models must use the same observations.
# ============================================================

model_vars <- unique(
  c(
    "follow_up_time",
    "event",
    predictors_extended
  )
)

model_data <- data[
  complete.cases(data[, model_vars]),
  ,
  drop = FALSE
]


cat(
  "Analysis sample size:",
  nrow(model_data),
  "\n"
)

cat(
  "Number of events:",
  sum(model_data$event == 1),
  "\n"
)

cat(
  "Number censored:",
  sum(model_data$event == 0),
  "\n"
)


# ============================================================
# 6. Construct the two Cox model formulas
# ============================================================

formula_base <- as.formula(
  paste(
    "Surv(follow_up_time, event) ~",
    paste(
      predictors_base,
      collapse = " + "
    )
  )
)


formula_extended <- as.formula(
  paste(
    "Surv(follow_up_time, event) ~",
    paste(
      predictors_extended,
      collapse = " + "
    )
  )
)


cat(
  "\nBase model formula:\n"
)

print(formula_base)


cat(
  "\nBase + exposure model formula:\n"
)

print(formula_extended)


# ============================================================
# 7. Fit both models using the complete-case dataset
# ============================================================

fit_base <- coxph(
  formula_base,
  data = model_data,
  x = TRUE,
  y = TRUE,
  model = TRUE
)


fit_extended <- coxph(
  formula_extended,
  data = model_data,
  x = TRUE,
  y = TRUE,
  model = TRUE
)


# ============================================================
# 8. Display the model results
# ============================================================

cat(
  "\n====================================\n",
  "Base model\n",
  "====================================\n"
)

print(
  summary(fit_base)
)


cat(
  "\n====================================\n",
  "Base + exposure model\n",
  "====================================\n"
)

print(
  summary(fit_extended)
)


# ============================================================
# 9. HR, 95% CI, and P value for the exposure variable
# ============================================================

coef_extended <- summary(
  fit_extended
)$coefficients

ci_extended <- summary(
  fit_extended
)$conf.int


# Identify the coefficient corresponding to the exposure variable
extended_row <- grep(
  "^exposure",
  rownames(coef_extended)
)


if (length(extended_row) == 1) {
  
  HR_extended <- ci_extended[
    extended_row,
    "exp(coef)"
  ]
  
  Lower_extended <- ci_extended[
    extended_row,
    "lower .95"
  ]
  
  Upper_extended <- ci_extended[
    extended_row,
    "upper .95"
  ]
  
  P_extended <- coef_extended[
    extended_row,
    "Pr(>|z|)"
  ]
  
  
  cat(
    "\nexposure HR =",
    round(HR_extended, 3),
    "\n95% CI =",
    round(Lower_extended, 3),
    "to",
    round(Upper_extended, 3),
    "\nP =",
    round(P_extended, 4),
    "\n"
  )
}


# ============================================================
# 10. Nested model likelihood-ratio test
#
# Test whether adding the exposure variable improves overall model fit
# relative to the same base model.
#
# This P value describes the incremental improvement in model fit;
# it is not a P value for the change in C-index.
# ============================================================

lrt_result <- anova(
  fit_base,
  fit_extended,
  test = "LRT"
)

cat(
  "\n====================================\n",
  "Likelihood ratio test\n",
  "Base vs Base + exposure\n",
  "====================================\n"
)

print(
  lrt_result
)


# ============================================================
# 11. Apparent C-index
# ============================================================

Y <- Surv(
  model_data$follow_up_time,
  model_data$event
)


# Base risk score
lp_base_app <- predict(
  fit_base,
  newdata = model_data,
  type = "lp"
)


# Base + exposure risk score
lp_extended_app <- predict(
  fit_extended,
  newdata = model_data,
  type = "lp"
)


# C-index
C_base_app <- concordance(
  Y ~ lp_base_app,
  reverse = TRUE
)$concordance


C_extended_app <- concordance(
  Y ~ lp_extended_app,
  reverse = TRUE
)$concordance


Delta_C_app <-
  C_extended_app -
  C_base_app


cat(
  "\n====================================\n",
  "Apparent performance\n",
  "====================================\n",
  
  "Base apparent C-index =",
  round(C_base_app, 4),
  "\n",
  
  "Base + exposure apparent C-index =",
  round(C_extended_app, 4),
  "\n",
  
  "Apparent Delta C-index =",
  round(Delta_C_app, 4),
  "\n"
)


# ============================================================
# 12. Bootstrap function for optimism correction
#
# In each bootstrap iteration:
#
# 1) Sample observations with replacement.
# 2) Refit the base Cox model in the bootstrap sample.
# 3) Refit the extended Cox model in the bootstrap sample.
# 4) Calculate the C-index in the bootstrap sample.
# 5) Apply the same fitted model to the original dataset.
# 6) Calculate optimism as the difference between the two estimates.
# ============================================================

bootstrap_optimism_cox <- function(
    data,
    indices
) {
  
  # ----------------------------------------------------------
  # A. Bootstrap sample
  # ----------------------------------------------------------
  
  boot_data <- data[
    indices,
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # B. Check event status
  # ----------------------------------------------------------
  
  event_b <- boot_data$event
  
  if (
    length(unique(event_b)) < 2 ||
    sum(event_b == 1) < 5 ||
    sum(event_b == 0) < 5
  ) {
    
    return(
      rep(
        NA,
        9
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # C. Refit the base Cox model in the bootstrap sample
  # ----------------------------------------------------------
  
  fit_base_b <- try(
    coxph(
      formula_base,
      data = boot_data,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    silent = TRUE
  )
  
  
  # ----------------------------------------------------------
  # D. Refit the extended Cox model in the bootstrap sample
  # ----------------------------------------------------------
  
  fit_extended_b <- try(
    coxph(
      formula_extended,
      data = boot_data,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    silent = TRUE
  )
  
  
  if (
    inherits(
      fit_base_b,
      "try-error"
    ) ||
    inherits(
      fit_extended_b,
      "try-error"
    )
  ) {
    
    return(
      rep(
        NA,
        9
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # E. Risk scores in the bootstrap sample
  # ----------------------------------------------------------
  
  lp_base_boot <- try(
    predict(
      fit_base_b,
      newdata = boot_data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  
  lp_extended_boot <- try(
    predict(
      fit_extended_b,
      newdata = boot_data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  
  # ----------------------------------------------------------
  # F. Risk scores from the same models in the original dataset
  # ----------------------------------------------------------
  
  lp_base_orig <- try(
    predict(
      fit_base_b,
      newdata = data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  
  lp_extended_orig <- try(
    predict(
      fit_extended_b,
      newdata = data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  
  if (
    inherits(
      lp_base_boot,
      "try-error"
    ) ||
    inherits(
      lp_extended_boot,
      "try-error"
    ) ||
    inherits(
      lp_base_orig,
      "try-error"
    ) ||
    inherits(
      lp_extended_orig,
      "try-error"
    )
  ) {
    
    return(
      rep(
        NA,
        9
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # G. Bootstrap sample C-index
  # ----------------------------------------------------------
  
  Y_boot <- Surv(
    boot_data$follow_up_time,
    boot_data$event
  )
  
  
  C_base_boot <- concordance(
    Y_boot ~ lp_base_boot,
    reverse = TRUE
  )$concordance
  
  
  C_extended_boot <- concordance(
    Y_boot ~ lp_extended_boot,
    reverse = TRUE
  )$concordance
  
  
  # ----------------------------------------------------------
  # H. C-index in the original dataset
  # ----------------------------------------------------------
  
  Y_orig <- Surv(
    data$follow_up_time,
    data$event
  )
  
  
  C_base_orig <- concordance(
    Y_orig ~ lp_base_orig,
    reverse = TRUE
  )$concordance
  
  
  C_extended_orig <- concordance(
    Y_orig ~ lp_extended_orig,
    reverse = TRUE
  )$concordance
  
  
  # ----------------------------------------------------------
  # I. Optimism
  # ----------------------------------------------------------
  
  optimism_base <-
    C_base_boot -
    C_base_orig
  
  
  optimism_extended <-
    C_extended_boot -
    C_extended_orig
  
  
  # ----------------------------------------------------------
  # J. Delta C
  # ----------------------------------------------------------
  
  delta_boot <-
    C_extended_boot -
    C_base_boot
  
  
  delta_orig <-
    C_extended_orig -
    C_base_orig
  
  
  optimism_delta <-
    delta_boot -
    delta_orig
  
  
  return(
    c(
      C_base_boot = C_base_boot,
      C_base_orig = C_base_orig,
      optimism_base = optimism_base,
      
      C_extended_boot = C_extended_boot,
      C_extended_orig = C_extended_orig,
      optimism_extended = optimism_extended,
      
      delta_boot = delta_boot,
      delta_orig = delta_orig,
      optimism_delta = optimism_delta
    )
  )
}


# ============================================================
# 13. Run the bootstrap procedure
# ============================================================

set.seed(2026)

boot_cox <- boot(
  data = model_data,
  statistic = bootstrap_optimism_cox,
  R = 2000
)


# ============================================================
# 14. Organize the bootstrap results
# ============================================================

boot_values <- as.data.frame(
  boot_cox$t
)


names(
  boot_values
) <- c(
  "C_base_boot",
  "C_base_orig",
  "optimism_base",
  
  "C_extended_boot",
  "C_extended_orig",
  "optimism_extended",
  
  "delta_boot",
  "delta_orig",
  "optimism_delta"
)


# Remove failed bootstrap iterations
boot_values <- boot_values[
  complete.cases(
    boot_values
  ),
  ,
  drop = FALSE
]


cat(
  "\nNumber of valid bootstrap iterations:",
  nrow(boot_values),
  "\n"
)


# ============================================================
# 15. Calculate mean optimism
# ============================================================

mean_optimism_base <- mean(
  boot_values$optimism_base
)


mean_optimism_extended <- mean(
  boot_values$optimism_extended
)


mean_optimism_delta <- mean(
  boot_values$optimism_delta
)


cat(
  "\nMean optimism - Base =",
  round(
    mean_optimism_base,
    4
  ),
  "\n"
)


cat(
  "Mean optimism - Base + exposure =",
  round(
    mean_optimism_extended,
    4
  ),
  "\n"
)


cat(
  "Mean optimism - Delta C =",
  round(
    mean_optimism_delta,
    4
  ),
  "\n"
)


# ============================================================
# 16. Optimism-corrected C-index
# ============================================================

C_base_corrected <-
  C_base_app -
  mean_optimism_base


C_extended_corrected <-
  C_extended_app -
  mean_optimism_extended


Delta_C_corrected <-
  C_extended_corrected -
  C_base_corrected


# Optional verification:
Delta_C_corrected_check <-
  Delta_C_app -
  mean_optimism_delta


cat(
  "\n====================================\n",
  "Optimism-corrected performance\n",
  "====================================\n",
  
  "Base apparent C-index =",
  round(
    C_base_app,
    4
  ),
  "\n",
  
  "Base optimism-corrected C-index =",
  round(
    C_base_corrected,
    4
  ),
  "\n\n",
  
  "Base + exposure apparent C-index =",
  round(
    C_extended_app,
    4
  ),
  "\n",
  
  "Base + exposure optimism-corrected C-index =",
  round(
    C_extended_corrected,
    4
  ),
  "\n\n",
  
  "Apparent Delta C-index =",
  round(
    Delta_C_app,
    4
  ),
  "\n",
  
  "Optimism-corrected Delta C-index =",
  round(
    Delta_C_corrected,
    4
  ),
  "\n",
  
  "Check Delta C-index =",
  round(
    Delta_C_corrected_check,
    4
  ),
  "\n"
)


# ============================================================
# 17. Create the final results table
# ============================================================

result_cox <- data.frame(
  
  Model = c(
    "Base clinical model",
    "Base clinical model + exposure",
    "Incremental value"
  ),
  
  Apparent_C_index = c(
    C_base_app,
    C_extended_app,
    Delta_C_app
  ),
  
  Mean_optimism = c(
    mean_optimism_base,
    mean_optimism_extended,
    mean_optimism_delta
  ),
  
  Optimism_corrected_C_index = c(
    C_base_corrected,
    C_extended_corrected,
    Delta_C_corrected
  )
)


print(
  result_cox,
  digits = 4
)
# ============================================================
# 18. Calculate the 95% CI for the optimism-corrected change in C-index
#     Location-shifted bootstrap CI
# ============================================================

# ------------------------------------------------------------
# 1. Point estimate
# ------------------------------------------------------------

cat(
  "Apparent Delta C-index =",
  round(Delta_C_app, 5),
  "\n"
)

cat(
  "Optimism-corrected Delta C-index =",
  round(Delta_C_corrected, 5),
  "\n"
)


# ------------------------------------------------------------
# 2. Mean optimism in the change in C-index
# ------------------------------------------------------------

optimism_delta <- mean(
  boot_values$optimism_delta,
  na.rm = TRUE
)

cat(
  "Mean optimism of Delta C =",
  round(optimism_delta, 5),
  "\n"
)


# ------------------------------------------------------------
# 3. Apparent change in C-index from each bootstrap iteration
# ------------------------------------------------------------

delta_boot <- boot_values$delta_boot

delta_boot <- delta_boot[
  is.finite(delta_boot)
]

cat(
  "Valid bootstrap replicates =",
  length(delta_boot),
  "\n"
)


# ------------------------------------------------------------
# 4. Location shift
#
# corrected Delta C
# = apparent Delta C - optimism
#
# Shift the bootstrap sampling distribution by the mean optimism.
# ------------------------------------------------------------

delta_boot_shifted <-
  delta_boot -
  optimism_delta


# ------------------------------------------------------------
# 5. Calculate the 95% CI
# ------------------------------------------------------------

CI_delta_corrected <- quantile(
  delta_boot_shifted,
  probs = c(0.025, 0.975),
  na.rm = TRUE,
  names = FALSE
)


# ------------------------------------------------------------
# 6. Display the results
# ------------------------------------------------------------

cat(
  "\n====================================\n",
  "Optimism-corrected Delta C-index\n",
  "====================================\n",
  
  "Delta C-index =",
  sprintf("%.4f", Delta_C_corrected),
  "\n",
  
  "95% CI =",
  sprintf("%.4f", CI_delta_corrected[1]),
  "to",
  sprintf("%.4f", CI_delta_corrected[2]),
  "\n"
)
# ============================================================
# Bootstrap-based P value for corrected Delta C-index
# H0: Delta C = 0
# ============================================================

x <- delta_boot_shifted[
  is.finite(delta_boot_shifted)
]

B <- length(x)

# Proportion of corrected bootstrap estimates less than or equal to zero
p_lower <- (
  sum(x <= 0) + 1
) / (B + 1)

# Proportion of corrected bootstrap estimates greater than or equal to zero
p_upper <- (
  sum(x >= 0) + 1
) / (B + 1)

# Two-sided P value
P_delta <- min(
  1,
  2 * min(
    p_lower,
    p_upper
  )
)

cat(
  "Corrected Delta C-index =",
  sprintf("%.4f", Delta_C_corrected),
  "\n95% CI =",
  sprintf("%.4f", CI_delta_corrected[1]),
  "to",
  sprintf("%.4f", CI_delta_corrected[2]),
  "\nBootstrap P =",
  ifelse(
    P_delta < 0.001,
    "<0.001",
    sprintf("%.3f", P_delta)
  ),
  "\n"
)

# ============================================================
# 18. Export the final analysis results to Excel
# ============================================================

library(openxlsx)


# ============================================================
# A. Model performance
# ============================================================

model_performance <- data.frame(
  Model = c(
    "Base clinical model",
    "Base clinical model + exposure",
    "Incremental value"
  ),
  
  Apparent_C_index = c(
    C_base_app,
    C_extended_app,
    Delta_C_app
  ),
  
  Mean_optimism = c(
    mean_optimism_base,
    mean_optimism_extended,
    mean_optimism_delta
  ),
  
  Optimism_corrected_C_index = c(
    C_base_corrected,
    C_extended_corrected,
    Delta_C_corrected
  )
)


# ============================================================
# B. Incremental prognostic value
# ============================================================

# Extract the likelihood-ratio test results
LRT_chisq <- lrt_result$Chisq[2]
LRT_df <- lrt_result$Df[2]
LRT_p <- lrt_result$`Pr(>|Chi|)`[2]


incremental_value <- data.frame(
  Metric = c(
    "Base optimism-corrected C-index",
    "Base + exposure optimism-corrected C-index",
    "Optimism-corrected Delta C-index",
    "Delta C-index 95% CI lower",
    "Delta C-index 95% CI upper",
    "Likelihood-ratio Chi-square",
    "Likelihood-ratio df",
    "Likelihood-ratio P value"
  ),
  
  Value = c(
    C_base_corrected,
    C_extended_corrected,
    Delta_C_corrected,
    CI_delta_corrected[1],
    CI_delta_corrected[2],
    LRT_chisq,
    LRT_df,
    LRT_p
  )
)


# ============================================================
# C. Effect of the exposure variable in the extended Cox model
# ============================================================

extendedopenia_effect <- data.frame(
  Variable = "exposure",
  HR = HR_extended,
  Lower_95CI = Lower_extended,
  Upper_95CI = Upper_extended,
  P_value = P_extended
)


# ============================================================
# D. Summary of bootstrap internal validation
# ============================================================

bootstrap_summary <- data.frame(
  Metric = c(
    "Number of requested bootstrap resamples",
    "Number of valid bootstrap resamples",
    "Mean optimism - Base",
    "Mean optimism - Base + exposure",
    "Mean optimism - Delta C"
  ),
  
  Value = c(
    2000,
    nrow(boot_values),
    mean_optimism_base,
    mean_optimism_extended,
    mean_optimism_delta
  )
)


# ============================================================
# E. Full model coefficients
# ============================================================

base_coef <- data.frame(
  Variable = rownames(summary(fit_base)$coefficients),
  HR = exp(summary(fit_base)$coefficients[, "coef"]),
  Lower_95CI = summary(fit_base)$conf.int[, "lower .95"],
  Upper_95CI = summary(fit_base)$conf.int[, "upper .95"],
  P_value = summary(fit_base)$coefficients[, "Pr(>|z|)"],
  row.names = NULL
)


extended_coef <- data.frame(
  Variable = rownames(summary(fit_extended)$coefficients),
  HR = exp(summary(fit_extended)$coefficients[, "coef"]),
  Lower_95CI = summary(fit_extended)$conf.int[, "lower .95"],
  Upper_95CI = summary(fit_extended)$conf.int[, "upper .95"],
  P_value = summary(fit_extended)$coefficients[, "Pr(>|z|)"],
  row.names = NULL
)


# ============================================================
# F. Create the Excel workbook
# ============================================================

wb <- createWorkbook()


# ------------------------------------------------------------
# Sheet 1
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Model_performance"
)

writeData(
  wb,
  "Model_performance",
  model_performance
)


# ------------------------------------------------------------
# Sheet 2
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Incremental_value"
)

writeData(
  wb,
  "Incremental_value",
  incremental_value
)


# ------------------------------------------------------------
# Sheet 3
# ------------------------------------------------------------

addWorksheet(
  wb,
  "exposure_effect"
)

writeData(
  wb,
  "exposure_effect",
  extendedopenia_effect
)


# ------------------------------------------------------------
# Sheet 4
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Bootstrap_summary"
)

writeData(
  wb,
  "Bootstrap_summary",
  bootstrap_summary
)


# ------------------------------------------------------------
# Sheet 5
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Base_model_Cox"
)

writeData(
  wb,
  "Base_model_Cox",
  base_coef
)


# ------------------------------------------------------------
# Sheet 6
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Extended_model_Cox"
)

writeData(
  wb,
  "Extended_model_Cox",
  extended_coef
)

# ============================================================
# Organize the likelihood-ratio test results
# ============================================================

LRT_chisq <- lrt_result$Chisq[2]
LRT_df    <- lrt_result$Df[2]
LRT_p     <- lrt_result$`Pr(>|Chi|)`[2]


lrt_table <- data.frame(
  Comparison = "Base clinical model vs Base + exposure",
  Chi_square = LRT_chisq,
  df = LRT_df,
  P_value = LRT_p
)


# ============================================================
# Write results to Excel
# ============================================================

addWorksheet(
  wb,
  "Model_comparison_LRT"
)

writeData(
  wb,
  "Model_comparison_LRT",
  lrt_table
)


# ============================================================
# G. Apply consistent header formatting
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (sheet_name in names(wb)) {
  
  addStyle(
    wb,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = 1:20,
    gridExpand = TRUE,
    stack = TRUE
  )
  
  setColWidths(
    wb,
    sheet = sheet_name,
    cols = 1:20,
    widths = "auto"
  )
}


# ============================================================
# H. Round numerical values to four decimal places
# ============================================================

num_style <- createStyle(
  numFmt = "0.0000"
)

for (sheet_name in names(wb)) {
  
  sheet_data <- switch(
    sheet_name,
    "Model_performance" = model_performance,
    "Incremental_value" = incremental_value,
    "exposure_effect" = extendedopenia_effect,
    "Bootstrap_summary" = bootstrap_summary,
    "Base_model_Cox" = base_coef,
    "Extended_model_Cox" = extended_coef
  )
  
  if (!is.null(sheet_data)) {
    
    numeric_cols <- which(
      sapply(
        sheet_data,
        is.numeric
      )
    )
    
    if (length(numeric_cols) > 0) {
      
      addStyle(
        wb,
        sheet = sheet_name,
        style = num_style,
        rows = 2:(nrow(sheet_data) + 1),
        cols = numeric_cols,
        gridExpand = TRUE,
        stack = TRUE
      )
    }
  }
}


# ============================================================
# I. Save the workbook
# ============================================================

saveWorkbook(
  wb,
  file.path(table_dir, "model_increment_results.xlsx"),
  overwrite = TRUE
)

cat(
  "\nResults saved to results/tables/model_increment_results.xlsx\n"
)
#================================================================================
# ============================================================
# Bootstrap C-index distribution
# Violin + boxplot + optimism-corrected point estimate
# ============================================================
# Plot model performance
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

optimism_plot_data <- boot_values %>%
  select(
    optimism_base,
    optimism_extended
  ) %>%
  rename(
    "Base clinical model" = optimism_base,
    "Base + exposure" = optimism_extended
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Model",
    values_to = "Optimism"
  )

optimism_plot_data$Model <- factor(
  optimism_plot_data$Model,
  levels = c(
    "Base clinical model",
    "Base + exposure"
  )
)

mean_points <- optimism_plot_data %>%
  group_by(Model) %>%
  summarise(
    Mean_optimism = mean(Optimism, na.rm = TRUE),
    .groups = "drop"
  )

p_optimism <- ggplot(
  optimism_plot_data,
  aes(
    x = Model,
    y = Optimism,
    fill = Model
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.35,
    width = 0.75,
    linewidth = 0.6
  ) +
  geom_boxplot(
    width = 0.20,
    outlier.shape = NA,
    alpha = 0.65,
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    data = mean_points,
    aes(
      x = Model,
      y = Mean_optimism
    ),
    inherit.aes = FALSE,
    shape = 18,
    size = 4
  ) +
  geom_text(
    data = mean_points,
    aes(
      x = Model,
      y = Mean_optimism,
      label = sprintf("%.3f", Mean_optimism)
    ),
    inherit.aes = FALSE,
    vjust = -1.2,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "Base clinical model" = "#7F8C8D",
      "Base + exposure" = "#4C78A8"
    )
  ) +
  labs(
    x = NULL,
    y = "Optimism in C-index"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 11.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 13
    )
  )

p_optimism
#========================================================================ROC
# Install the package before first use if necessary.
# install.packages("timeROC")

library(timeROC)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. Risk scores from the two Cox models
# ============================================================

lp_base <- predict(
  fit_base,
  newdata = model_data,
  type = "lp"
)

lp_extended <- predict(
  fit_extended,
  newdata = model_data,
  type = "lp"
)

# ============================================================
# 2. Define evaluation times
# ============================================================

auc_times <- seq(
  6,
  36,
  by = 1
)


# ============================================================
# 3. Base model time-dependent ROC
# ============================================================

roc_base <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = lp_base,
  cause = 1,
  weighting = "marginal",
  times = auc_times,
  iid = TRUE
)


# ============================================================
# 4. Base + exposure
# ============================================================

roc_extended <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = lp_extended,
  cause = 1,
  weighting = "marginal",
  times = auc_times,
  iid = TRUE
)


# ============================================================
# 5. Organize the AUC results
# ============================================================

auc_data <- data.frame(
  Time = auc_times,
  Base = roc_base$AUC,
  exposure = roc_extended$AUC
)


auc_long <- auc_data %>%
  tidyr::pivot_longer(
    cols = c(
      Base,
      exposure
    ),
    names_to = "Model",
    values_to = "AUC"
  )


auc_long$Model <- factor(
  auc_long$Model,
  levels = c(
    "Base",
    "exposure"
  ),
  labels = c(
    "Base clinical model",
    "Base + exposure"
  )
)
# ============================================================
# 6. Plot the time-dependent AUC curves
# ============================================================

p_auc <- ggplot(
  auc_long,
  aes(
    x = Time,
    y = AUC,
    color = Model
  )
) +
  
  geom_line(
    linewidth = 1.2
  ) +
  
  scale_x_continuous(
    breaks = seq(
      6,
      36,
      by = 6
    )
  ) +
  
  scale_y_continuous(
    limits = c(
      0.50,
      0.90
    ),
    breaks = seq(
      0.50,
      0.90,
      by = 0.10
    )
  ) +
  
  labs(
    x = "Time (months)",
    y = "Time-dependent AUC",
    color = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = c(
      0.72,
      0.15
    ),
    
    legend.background = element_blank(),
    
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      color = "black"
    )
  )

p_auc
#=================================================================
# ============================================================
# 18. Optimism-corrected time-dependent AUC
# ============================================================

library(timeROC)
library(boot)
library(ggplot2)
library(dplyr)
library(tidyr)


# ============================================================
# 1. Define evaluation times
# ============================================================

auc_times <- seq(
  6,
  36,
  by = 1
)


# ============================================================
# 2. Risk scores in the original complete-case dataset
# ============================================================

lp_base_full <- predict(
  fit_base,
  newdata = model_data,
  type = "lp"
)

lp_extended_full <- predict(
  fit_extended,
  newdata = model_data,
  type = "lp"
)


# ============================================================
# 3. Apparent time-dependent AUC
# ============================================================

tdroc_base_app <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = lp_base_full,
  cause = 1,
  weighting = "marginal",
  times = auc_times,
  iid = FALSE
)

tdroc_extended_app <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = lp_extended_full,
  cause = 1,
  weighting = "marginal",
  times = auc_times,
  iid = FALSE
)


AUC_base_app <- tdroc_base_app$AUC
AUC_extended_app <- tdroc_extended_app$AUC

Delta_AUC_app <-
  AUC_extended_app -
  AUC_base_app

# ============================================================
# 4. Bootstrap optimism function for AUC(t)
# ============================================================

bootstrap_auc_optimism <- function(
    data,
    indices
) {
  
  # ----------------------------------------------------------
  # A. Bootstrap sample
  # ----------------------------------------------------------
  
  boot_data <- data[
    indices,
    ,
    drop = FALSE
  ]
  
  
  # Ensure that both event states are represented
  if (
    length(unique(boot_data$event)) < 2
  ) {
    
    return(
      rep(
        NA,
        length(auc_times) * 9
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # B. Refit both Cox models in each bootstrap sample
  # ----------------------------------------------------------
  
  fit_base_b <- try(
    coxph(
      formula_base,
      data = boot_data,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    silent = TRUE
  )
  
  
  fit_extended_b <- try(
    coxph(
      formula_extended,
      data = boot_data,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    silent = TRUE
  )
  
  
  if (
    inherits(fit_base_b, "try-error") ||
    inherits(fit_extended_b, "try-error")
  ) {
    
    return(
      rep(
        NA,
        length(auc_times) * 9
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # C. Bootstrap sample risk scores
  # ----------------------------------------------------------
  
  lp_base_boot <- try(
    predict(
      fit_base_b,
      newdata = boot_data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  lp_extended_boot <- try(
    predict(
      fit_extended_b,
      newdata = boot_data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  
  # ----------------------------------------------------------
  # D. Apply each bootstrap-fitted model to the original dataset
  # ----------------------------------------------------------
  
  lp_base_orig <- try(
    predict(
      fit_base_b,
      newdata = data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  lp_extended_orig <- try(
    predict(
      fit_extended_b,
      newdata = data,
      type = "lp"
    ),
    silent = TRUE
  )
  
  
  if (
    inherits(lp_base_boot, "try-error") ||
    inherits(lp_extended_boot, "try-error") ||
    inherits(lp_base_orig, "try-error") ||
    inherits(lp_extended_orig, "try-error")
  ) {
    
    return(
      rep(
        NA,
        length(auc_times) * 9
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # E. AUC(t) in the bootstrap sample
  # ----------------------------------------------------------
  
  roc_base_boot <- try(
    timeROC(
      T = boot_data$follow_up_time,
      delta = boot_data$event,
      marker = lp_base_boot,
      cause = 1,
      weighting = "marginal",
      times = auc_times,
      iid = FALSE
    ),
    silent = TRUE
  )
  
  
  roc_extended_boot <- try(
    timeROC(
      T = boot_data$follow_up_time,
      delta = boot_data$event,
      marker = lp_extended_boot,
      cause = 1,
      weighting = "marginal",
      times = auc_times,
      iid = FALSE
    ),
    silent = TRUE
  )
  
  
  # ----------------------------------------------------------
  # F. AUC(t) in the original dataset
  # ----------------------------------------------------------
  
  roc_base_orig <- try(
    timeROC(
      T = data$follow_up_time,
      delta = data$event,
      marker = lp_base_orig,
      cause = 1,
      weighting = "marginal",
      times = auc_times,
      iid = FALSE
    ),
    silent = TRUE
  )
  
  
  roc_extended_orig <- try(
    timeROC(
      T = data$follow_up_time,
      delta = data$event,
      marker = lp_extended_orig,
      cause = 1,
      weighting = "marginal",
      times = auc_times,
      iid = FALSE
    ),
    silent = TRUE
  )
  
  
  if (
    inherits(roc_base_boot, "try-error") ||
    inherits(roc_extended_boot, "try-error") ||
    inherits(roc_base_orig, "try-error") ||
    inherits(roc_extended_orig, "try-error")
  ) {
    
    return(
      rep(
        NA,
        length(auc_times) * 9
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # G. Extract AUC estimates
  # ----------------------------------------------------------
  
  base_boot <- roc_base_boot$AUC
  base_orig <- roc_base_orig$AUC
  
  extended_boot <- roc_extended_boot$AUC
  extended_orig <- roc_extended_orig$AUC
  
  
  # ----------------------------------------------------------
  # H. Optimism
  # ----------------------------------------------------------
  
  optimism_base <-
    base_boot -
    base_orig
  
  optimism_extended <-
    extended_boot -
    extended_orig
  
  
  # ----------------------------------------------------------
  # I. Delta AUC
  # ----------------------------------------------------------
  
  delta_boot <-
    extended_boot -
    base_boot
  
  delta_orig <-
    extended_orig -
    base_orig
  
  optimism_delta <-
    delta_boot -
    delta_orig
  
  
  # ----------------------------------------------------------
  # J. Return results
  # ----------------------------------------------------------
  
  return(
    c(
      base_boot,
      base_orig,
      optimism_base,
      
      extended_boot,
      extended_orig,
      optimism_extended,
      
      delta_boot,
      delta_orig,
      optimism_delta
    )
  )
}

# ============================================================
# 5. Run bootstrap
# ============================================================

set.seed(2026)

boot_auc <- boot(
  data = model_data,
  statistic = bootstrap_auc_optimism,
  R = 2000
)

cat(
  "Bootstrap completed.\n"
)

# ============================================================
# 6. Organize bootstrap AUC results
# ============================================================

n_time <- length(
  auc_times
)

boot_auc_mat <- boot_auc$t


# ------------------------------------------------------------
# Base
# ------------------------------------------------------------

idx_base_boot <-
  1:n_time

idx_base_orig <-
  (n_time + 1):(2 * n_time)

idx_opt_base <-
  (2 * n_time + 1):(3 * n_time)


# ------------------------------------------------------------
# Base + exposure
# ------------------------------------------------------------

idx_extended_boot <-
  (3 * n_time + 1):(4 * n_time)

idx_extended_orig <-
  (4 * n_time + 1):(5 * n_time)

idx_opt_extended <-
  (5 * n_time + 1):(6 * n_time)


# ------------------------------------------------------------
# Delta
# ------------------------------------------------------------

idx_delta_boot <-
  (6 * n_time + 1):(7 * n_time)

idx_delta_orig <-
  (7 * n_time + 1):(8 * n_time)

idx_opt_delta <-
  (8 * n_time + 1):(9 * n_time)

# ============================================================
# 7. Mean optimism at each time point
# ============================================================

mean_opt_base_auc <- apply(
  boot_auc_mat[, idx_opt_base, drop = FALSE],
  2,
  mean,
  na.rm = TRUE
)


mean_opt_extended_auc <- apply(
  boot_auc_mat[, idx_opt_extended, drop = FALSE],
  2,
  mean,
  na.rm = TRUE
)


mean_opt_delta_auc <- apply(
  boot_auc_mat[, idx_opt_delta, drop = FALSE],
  2,
  mean,
  na.rm = TRUE
)

# ============================================================
# 8. Optimism-corrected AUC(t)
# ============================================================

AUC_base_corrected <-
  AUC_base_app -
  mean_opt_base_auc


AUC_extended_corrected <-
  AUC_extended_app -
  mean_opt_extended_auc


Delta_AUC_corrected <-
  AUC_extended_corrected -
  AUC_base_corrected

auc_corrected_results <- data.frame(
  Time = auc_times,
  
  Base_apparent =
    AUC_base_app,
  
  Base_corrected =
    AUC_base_corrected,
  
  Sarc_apparent =
    AUC_extended_app,
  
  Sarc_corrected =
    AUC_extended_corrected,
  
  Delta_apparent =
    Delta_AUC_app,
  
  Delta_corrected =
    Delta_AUC_corrected
)


print(
  auc_corrected_results,
  digits = 4
)

# ============================================================
# 9. Corrected time-dependent AUC curve
# ============================================================

auc_plot_data <- data.frame(
  Time = auc_times,
  
  `Base clinical model` =
    AUC_base_corrected,
  
  `Base + exposure` =
    AUC_extended_corrected,
  
  check.names = FALSE
)


auc_plot_long <- auc_plot_data %>%
  pivot_longer(
    cols = -Time,
    names_to = "Model",
    values_to = "AUC"
  )


auc_plot_long$Model <- factor(
  auc_plot_long$Model,
  levels = c(
    "Base clinical model",
    "Base + exposure"
  )
)


# ============================================================
# Publication-style optimism-corrected time-dependent AUC
# ============================================================

# ============================================================
# Publication-style optimism-corrected time-dependent AUC
# ============================================================

library(ggplot2)

p_auc_corrected <- ggplot(
  auc_plot_long,
  aes(
    x = Time,
    y = AUC,
    color = Model
  )
) +
  
  # AUC curves
  geom_line(
    linewidth = 1.35,
    lineend = "round"
  ) +
  
  # Colors
  scale_color_manual(
    values = c(
      "Base clinical model" = "#D95F59",
      "Base + exposure"   = "#168AAD"
    )
  ) +
  
  # X axis
  scale_x_continuous(
    breaks = c(6, 12, 18, 24, 30, 36),
    limits = c(5.5, 36.5),
    expand = c(0.01, 0.01)
  ) +
  
  # Y axis
  scale_y_continuous(
    breaks = seq(0.50, 0.80, by = 0.05),
    limits = c(0.50, 0.80),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Time (months)",
    y = "Optimism-corrected AUC",
    color = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    # Legend
    legend.position = "top",
    legend.justification = "center",
    
    legend.text = element_text(
      size = 11.5,
      color = "black"
    ),
    
    legend.key.width = grid::unit(
      1.5,
      "cm"
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 14,
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    
    # Axis text
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    
    # Axis lines
    axis.line = element_line(
      linewidth = 0.7,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    axis.ticks.length = grid::unit(
      0.15,
      "cm"
    ),
    
    # Margins
    plot.margin = margin(
      t = 5,
      r = 12,
      b = 5,
      l = 5
    )
  )

p_auc_corrected


# ============================================================
# Apparent time-dependent ROC curves
# 12 / 24 / 36 months
# ============================================================

library(timeROC)
library(ggplot2)
library(dplyr)
library(patchwork)


# ============================================================
# 1. Risk scores from the original fitted models
# ============================================================

lp_base_app <- predict(
  fit_base,
  newdata = model_data,
  type = "lp"
)

lp_extended_app <- predict(
  fit_extended,
  newdata = model_data,
  type = "lp"
)


# ============================================================
# 2. Calculate time-dependent ROC
# ============================================================

roc_times <- c(12, 24, 36)

roc_base <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = lp_base_app,
  cause = 1,
  weighting = "marginal",
  times = roc_times,
  iid = TRUE
)

roc_extended <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = lp_extended_app,
  cause = 1,
  weighting = "marginal",
  times = roc_times,
  iid = TRUE
)


# ============================================================
# 3. Function for drawing each ROC panel
# ============================================================

make_roc_plot <- function(
    time_index,
    time_month
) {
  
  # ----------------------------------------------------------
  # ROC coordinates
  # ----------------------------------------------------------
  
  df_base <- data.frame(
    FPR = roc_base$FP[, time_index],
    TPR = roc_base$TP[, time_index],
    Model = "Base clinical model"
  )
  
  df_extended <- data.frame(
    FPR = roc_extended$FP[, time_index],
    TPR = roc_extended$TP[, time_index],
    Model = "Base + exposure"
  )
  
  roc_df <- bind_rows(
    df_base,
    df_extended
  )
  
  roc_df$Model <- factor(
    roc_df$Model,
    levels = c(
      "Base clinical model",
      "Base + exposure"
    )
  )
  
  
  # ----------------------------------------------------------
  # AUC values
  # ----------------------------------------------------------
  
  auc_base <- roc_base$AUC[time_index]
  auc_extended <- roc_extended$AUC[time_index]
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    roc_df,
    aes(
      x = FPR,
      y = TPR,
      color = Model
    )
  ) +
    
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.7,
      color = "grey55"
    ) +
    
    geom_line(
      linewidth = 1.25,
      lineend = "round"
    ) +
    
    scale_color_manual(
      values = c(
        "Base clinical model" = "#D95F59",
        "Base + exposure"   = "#168AAD"
      )
    ) +
    
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      expand = c(0, 0)
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(
        0,
        1,
        by = 0.2
      ),
      expand = c(0, 0)
    ) +
    
    coord_equal() +
    
    labs(
      title = paste0(
        time_month,
        "-month ROC"
      ),
      x = "1 - Specificity",
      y = "Sensitivity",
      color = NULL
    ) +
    
    # --------------------------------------------------------
  # AUC annotation
  # --------------------------------------------------------
  
  annotate(
    "text",
    x = 0.97,
    y = 0.22,
    hjust = 1,
    label = paste0(
      "Base: AUC = ",
      sprintf("%.3f", auc_base)
    ),
    size = 3.7,
    color = "#D95F59"
  ) +
    
    annotate(
      "text",
      x = 0.97,
      y = 0.13,
      hjust = 1,
      label = paste0(
        "Base + exposure: AUC = ",
        sprintf("%.3f", auc_extended)
      ),
      size = 3.7,
      color = "#168AAD"
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    theme(
      plot.title = element_text(
        size = 13.5,
        hjust = 0.5,
        face = "bold"
      ),
      
      axis.title = element_text(
        size = 12.5
      ),
      
      axis.text = element_text(
        size = 11,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.65,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.55,
        color = "black"
      ),
      
      legend.position = "none",
      
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
  
  return(p)
}


# ============================================================
# 4. Generate 12 / 24 / 36 month ROC curves
# ============================================================

p_roc12 <- make_roc_plot(
  time_index = 1,
  time_month = 12
)

p_roc24 <- make_roc_plot(
  time_index = 2,
  time_month = 24
)

p_roc36 <- make_roc_plot(
  time_index = 3,
  time_month = 36
)


# ============================================================
# 5. Display separately
# ============================================================

p_roc12
p_roc24
p_roc36

# ============================================================
# 6. Combine the three ROC panels
# ============================================================

p_roc_combined <-
  p_roc12 +
  p_roc24 +
  p_roc36 +
  plot_annotation(
    tag_levels = "A"
  )

p_roc_combined

ggsave(
  file.path(figure_dir, "time_dependent_roc.pdf"),
  plot = p_roc_combined,
  width = 12,
  height = 4.3,
  units = "in",
  device = cairo_pdf
)

roc_auc_results <- data.frame(
  Time_months = c(12, 24, 36),
  
  Base_apparent_AUC =
    roc_base$AUC,
  
  Base_plus_exposure_apparent_AUC =
    roc_extended$AUC,
  
  Apparent_Delta_AUC =
    roc_extended$AUC -
    roc_base$AUC
)

print(
  roc_auc_results,
  digits = 4
)
# ============================================================
# Calculate 95% confidence intervals
# ============================================================
# 12 / 24 / 36-month AUC + 95% CI
# Apparent + optimism-corrected
# ============================================================

library(timeROC)
library(dplyr)
library(openxlsx)

key_times <- c(12, 24, 36)

# Previously defined evaluation times:
# auc_times <- seq(6, 36, by = 1)

key_idx <- match(
  key_times,
  auc_times
)

if (any(is.na(key_idx))) {
  stop("12, 24 or 36 months not found in auc_times.")
}
# ============================================================
# Apparent AUC at 12 / 24 / 36 months
# ============================================================

roc_base_key <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = predict(fit_base, model_data, type = "lp"),
  cause = 1,
  weighting = "marginal",
  times = key_times,
  iid = TRUE
)

roc_extended_key <- timeROC(
  T = model_data$follow_up_time,
  delta = model_data$event,
  marker = predict(fit_extended, model_data, type = "lp"),
  cause = 1,
  weighting = "marginal",
  times = key_times,
  iid = TRUE
)


# 95% CI
ci_base_app <- confint(
  roc_base_key,
  level = 0.95
)

ci_extended_app <- confint(
  roc_extended_key,
  level = 0.95
)
# ============================================================
# Paired comparison of apparent AUC
# ============================================================

compare_app <- compare(
  roc_base_key,
  roc_extended_key,
  adjusted = FALSE
)

print(compare_app)

# ============================================================
# Corrected AUC bootstrap distributions
# ============================================================

base_boot_mat <- boot_auc_mat[
  ,
  idx_base_boot,
  drop = FALSE
]

extended_boot_mat <- boot_auc_mat[
  ,
  idx_extended_boot,
  drop = FALSE
]


# Delta bootstrap distribution
delta_boot_mat <- (
  extended_boot_mat -
    base_boot_mat
)

# ============================================================
# Containers
# ============================================================

Base_corr_lower <- numeric(length(key_times))
Base_corr_upper <- numeric(length(key_times))

Sarc_corr_lower <- numeric(length(key_times))
Sarc_corr_upper <- numeric(length(key_times))

Delta_corr_lower <- numeric(length(key_times))
Delta_corr_upper <- numeric(length(key_times))


# ============================================================
# Bootstrap percentile CI after location correction
# ============================================================

for (k in seq_along(key_times)) {
  
  j <- key_idx[k]
  
  
  # ----------------------------------------------------------
  # Base model
  # ----------------------------------------------------------
  
  x_base <- base_boot_mat[, j]
  
  x_base <- x_base[
    is.finite(x_base)
  ]
  
  # Shift bootstrap distribution so that its center corresponds
  # to the optimism-corrected estimate
  shift_base <-
    AUC_base_corrected[j] -
    AUC_base_app[j]
  
  x_base_corrected <-
    x_base +
    shift_base
  
  ci_base <- quantile(
    x_base_corrected,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  Base_corr_lower[k] <- ci_base[1]
  Base_corr_upper[k] <- ci_base[2]
  
  
  # ----------------------------------------------------------
  # Base + exposure
  # ----------------------------------------------------------
  
  x_extended <- extended_boot_mat[, j]
  
  x_extended <- x_extended[
    is.finite(x_extended)
  ]
  
  shift_extended <-
    AUC_extended_corrected[j] -
    AUC_extended_app[j]
  
  x_extended_corrected <-
    x_extended +
    shift_extended
  
  ci_extended <- quantile(
    x_extended_corrected,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  Sarc_corr_lower[k] <- ci_extended[1]
  Sarc_corr_upper[k] <- ci_extended[2]
  
  
  # ----------------------------------------------------------
  # Corrected Delta AUC
  # ----------------------------------------------------------
  
  x_delta <- delta_boot_mat[, j]
  
  x_delta <- x_delta[
    is.finite(x_delta)
  ]
  
  delta_app_j <-
    AUC_extended_app[j] -
    AUC_base_app[j]
  
  shift_delta <-
    Delta_AUC_corrected[j] -
    delta_app_j
  
  x_delta_corrected <-
    x_delta +
    shift_delta
  
  ci_delta <- quantile(
    x_delta_corrected,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  Delta_corr_lower[k] <- ci_delta[1]
  Delta_corr_upper[k] <- ci_delta[2]
}

# ============================================================
# Corrected results
# ============================================================

corrected_auc_table <- data.frame(
  
  Time_months = key_times,
  
  Base_corrected_AUC =
    AUC_base_corrected[key_idx],
  
  Base_corrected_Lower95 =
    Base_corr_lower,
  
  Base_corrected_Upper95 =
    Base_corr_upper,
  
  
  exposure_corrected_AUC =
    AUC_extended_corrected[key_idx],
  
  exposure_corrected_Lower95 =
    Sarc_corr_lower,
  
  exposure_corrected_Upper95 =
    Sarc_corr_upper,
  
  
  Corrected_Delta_AUC =
    Delta_AUC_corrected[key_idx],
  
  Delta_Lower95 =
    Delta_corr_lower,
  
  Delta_Upper95 =
    Delta_corr_upper
)


print(
  corrected_auc_table,
  digits = 4
)

apparent_auc_table <- data.frame(
  
  Time_months = key_times,
  
  Base_apparent_AUC =
    roc_base_key$AUC,
  
  Base_Lower95 =
    ci_base_app$CI_AUC[, 1],
  
  Base_Upper95 =
    ci_base_app$CI_AUC[, 2],
  
  
  exposure_apparent_AUC =
    roc_extended_key$AUC,
  
  exposure_Lower95 =
    ci_extended_app$CI_AUC[, 1],
  
  exposure_Upper95 =
    ci_extended_app$CI_AUC[, 2],
  
  
  Apparent_Delta_AUC =
    roc_extended_key$AUC -
    roc_base_key$AUC
)


print(
  apparent_auc_table,
  digits = 4
)


# ============================================================
# Publication-friendly corrected table
# ============================================================

corrected_auc_formatted <- data.frame(
  
  Time = paste0(
    key_times,
    " months"
  ),
  
  `Base clinical model` = sprintf(
    "%.3f (%.3f–%.3f)",
    AUC_base_corrected[key_idx],
    Base_corr_lower,
    Base_corr_upper
  ),
  
  `Base + exposure` = sprintf(
    "%.3f (%.3f–%.3f)",
    AUC_extended_corrected[key_idx],
    Sarc_corr_lower,
    Sarc_corr_upper
  ),
  
  `Delta AUC` = sprintf(
    "%.3f (%.3f–%.3f)",
    Delta_AUC_corrected[key_idx],
    Delta_corr_lower,
    Delta_corr_upper
  ),
  
  check.names = FALSE
)


print(
  corrected_auc_formatted
)

# ============================================================
# Export everything to Excel
# ============================================================

wb_auc <- createWorkbook()


# ------------------------------------------------------------
# Sheet 1: Corrected AUC
# ------------------------------------------------------------

addWorksheet(
  wb_auc,
  "Corrected_AUC"
)

writeData(
  wb_auc,
  "Corrected_AUC",
  corrected_auc_table
)


# ------------------------------------------------------------
# Sheet 2: Publication table
# ------------------------------------------------------------

addWorksheet(
  wb_auc,
  "Corrected_AUC_formatted"
)

writeData(
  wb_auc,
  "Corrected_AUC_formatted",
  corrected_auc_formatted
)


# ------------------------------------------------------------
# Sheet 3: Apparent AUC
# ------------------------------------------------------------

addWorksheet(
  wb_auc,
  "Apparent_AUC"
)

writeData(
  wb_auc,
  "Apparent_AUC",
  apparent_auc_table
)


# ------------------------------------------------------------
# Sheet 4: Full corrected AUC trajectory
# ------------------------------------------------------------

addWorksheet(
  wb_auc,
  "Full_AUC_trajectory"
)

writeData(
  wb_auc,
  "Full_AUC_trajectory",
  auc_corrected_results
)


# ============================================================
# Formatting
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (
  sheet_name in names(wb_auc)
) {
  
  addStyle(
    wb_auc,
    sheet = sheet_name,
    style = header_style,
    rows = 1,
    cols = 1:ncol(
      readWorkbook(
        wb_auc,
        sheet = sheet_name
      )
    ),
    gridExpand = TRUE
  )
  
  setColWidths(
    wb_auc,
    sheet = sheet_name,
    cols = 1:ncol(
      readWorkbook(
        wb_auc,
        sheet = sheet_name
      )
    ),
    widths = "auto"
  )
}


# ============================================================
# Save
# ============================================================

saveWorkbook(
  wb_auc,
  file = file.path(table_dir, "time_dependent_auc_results.xlsx"),
  overwrite = TRUE
)

cat(
  "\nSaved: results/tables/time_dependent_auc_results.xlsx\n"
)
