# ============================================================
# Overlap Weighting Analysis
#
# Description:
#   Perform propensity score overlap weighting (OW) to compare
#   exposed and unexposed groups. The script evaluates
#   covariate balance, produces diagnostic plots, fits an OW-weighted
#   Cox model for loss of response (event), and exports a weighted
#   baseline table.
#
# Required anonymized columns:
#   exposure, event, follow_up_time, and covariate_01 to covariate_16
# ============================================================

# Install the packages before first use if necessary:
# install.packages(c(
#   "WeightIt", "cobalt", "survey", "survival", "tableone",
#   "openxlsx", "ggplot2", "dplyr", "tibble", "tidyr", "patchwork"
# ))
library(WeightIt)
library(cobalt)
library(survey)
library(survival)
library(tableone)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)

# ============================================================
# 1. File paths
# ============================================================

# Run this script from the root directory of the GitHub repository.
input_file <- file.path("data", "analysis_data.csv")
figure_dir <- file.path("results", "figures")
table_dir <- file.path("results", "tables")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ============================================================
# 2. Define categorical variables
# ============================================================

factor_vars <- c(
  "covariate_02",
  "covariate_03",
  "covariate_04",
  "covariate_16",
  "covariate_08",
  "covariate_12"
)

data[factor_vars] <- lapply(
  data[factor_vars],
  as.factor
)

data$exposure <- as.factor(data$exposure)


# ============================================================
# 3. Specify baseline covariates for overlap weighting
#    Continuous variables are retained in their original continuous form.
# ============================================================

covariates <- c(
  "covariate_01",
  "covariate_02",
  "covariate_03",
  "covariate_04",
  "covariate_05",
  "covariate_13",
  "covariate_14",
  "covariate_15",
  "covariate_16",
  "covariate_06",
  "covariate_07",
  "covariate_08",
  "covariate_09",
  "covariate_10",
  "covariate_11",
  "covariate_12"
)


# ============================================================
# 4. Construct the propensity score formula
# ============================================================

formula_ps <- as.formula(
  paste(
    "exposure ~",
    paste(covariates, collapse = " + ")
  )
)

print(formula_ps)


# ============================================================
# 4. Overlap Weighting
#    estimand = "ATO"
# ============================================================

set.seed(2000)

w_ow <- weightit(
  formula_ps,
  data = data,
  method = "glm",
  estimand = "ATO"
)

summary(w_ow)
summary(w_ow$weights)

hist(
  w_ow$weights,
  breaks = 30,
  main = "Distribution of Overlap Weights",
  xlab = "Overlap weight"
)

max(w_ow$weights)
min(w_ow$weights)

bal_ow <- bal.tab(
  w_ow,
  un = TRUE,
  thresholds = c(m = 0.1)
)

print(bal_ow)
# ============================================================
# 5. Assess covariate balance and create a Love plot
# ============================================================

# ============================================================
# Extract SMDs before and after overlap weighting
# ============================================================

bal_ow <- bal.tab(
  w_ow,
  un = TRUE
)

love_df <- as.data.frame(bal_ow$Balance) %>%
  rownames_to_column("Variable") %>%
  select(
    Variable,
    Diff.Un,
    Diff.Adj
  ) %>%
  mutate(
    Before = abs(Diff.Un),
    After  = abs(Diff.Adj)
  )


# ============================================================
# Use anonymized variable names as plot labels
# ============================================================

love_df$Variable_label <- love_df$Variable

love_df$Variable_label[
  love_df$Variable == "prop.score"
] <- "Propensity score"

# ============================================================
# Exclude the propensity score itself from the Love plot
# ============================================================

love_df <- love_df %>%
  filter(Variable != "prop.score")


# ============================================================
# Sort variables by their absolute SMD before weighting
# ============================================================

love_df <- love_df %>%
  arrange(Before)

love_df$Variable_label <- factor(
  love_df$Variable_label,
  levels = love_df$Variable_label
)


# ============================================================
# Convert the balance results to long format for plotting
# ============================================================

love_long <- love_df %>%
  select(
    Variable_label,
    Before,
    After
  ) %>%
  tidyr::pivot_longer(
    cols = c(Before, After),
    names_to = "Stage",
    values_to = "SMD"
  ) %>%
  mutate(
    Stage = factor(
      Stage,
      levels = c("Before", "After"),
      labels = c(
        "Before weighting",
        "After overlap weighting"
      )
    )
  )


# ============================================================
# Create a publication-style Love plot
# ============================================================

p_love <- ggplot() +
  
  # Connect the estimates before and after weighting
  geom_segment(
    data = love_df,
    aes(
      x = Before,
      xend = After,
      y = Variable_label,
      yend = Variable_label
    ),
    linewidth = 0.55,
    color = "grey70"
  ) +
  
  # Plot the SMD estimates
  geom_point(
    data = love_long,
    aes(
      x = SMD,
      y = Variable_label,
      color = Stage,
      shape = Stage
    ),
    size = 3.2
  ) +
  
  # Add the conventional SMD threshold of 0.1
  geom_vline(
    xintercept = 0.1,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey35"
  ) +
  
  # Configure the horizontal axis
  scale_x_continuous(
    breaks = seq(0, 1.2, 0.2),
    limits = c(0, 1.2),
    expand = expansion(
      mult = c(0.01, 0.03)
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Before weighting" = "#D55E00",
      "After overlap weighting" = "#0072B2"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Before weighting" = 16,
      "After overlap weighting" = 17
    )
  ) +
  
  labs(
    x = "Absolute standardized mean difference",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    legend.position = "top",
    
    legend.justification = "left",
    
    legend.text = element_text(
      size = 11.5
    ),
    
    axis.text.x = element_text(
      size = 11.5,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 11.5,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 10)
    ),
    
    axis.line = element_line(
      linewidth = 0.6,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5
    ),
    
    plot.margin = margin(
      10, 20, 10, 10
    )
  )

p_love

ggsave(
  file.path(figure_dir, "overlap_weighting_love_plot.pdf"),
  plot = p_love,
  width = 8.5,
  height = 7,
  units = "in",
  device = cairo_pdf
)
# ============================================================
# 6. Plot propensity score distributions before and after OW
# ============================================================

# Add propensity scores and overlap weights to the analysis data
data$PS <- w_ow$ps
data$OW <- w_ow$weights

# Create descriptive group labels for plotting
data$exposure_plot <- factor(
  data$exposure,
  levels = c(0, 1),
  labels = c("Unexposed", "Exposed")
)

# ============================================================
# A. Propensity score distribution before overlap weighting
# ============================================================

p_before <- ggplot(
  data,
  aes(
    x = PS,
    fill = exposure_plot,
    color = exposure_plot
  )
) +
  geom_density(
    alpha = 0.30,
    linewidth = 0.9,
    adjust = 1
  ) +
  scale_fill_manual(
    values = c(
      "Unexposed" = "#4C78A8",
      "Exposed" = "#E58B65"
    )
  ) +
  scale_color_manual(
    values = c(
      "Unexposed" = "#4C78A8",
      "Exposed" = "#E58B65"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(add = c(0.02, 0.02))
  ) +
  labs(
    title = "Before overlap weighting",
    x = "Propensity score",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
      size = 13,
      face = "bold"
    ),
    axis.text = element_text(color = "black")
  )


# ============================================================
# B. Propensity score distribution after overlap weighting
# ============================================================

p_after <- ggplot(
  data,
  aes(
    x = PS,
    weight = OW,
    fill = exposure_plot,
    color = exposure_plot
  )
) +
  geom_density(
    alpha = 0.30,
    linewidth = 0.9,
    adjust = 1
  ) +
  scale_fill_manual(
    values = c(
      "Unexposed" = "#4C78A8",
      "Exposed" = "#E58B65"
    )
  ) +
  scale_color_manual(
    values = c(
      "Unexposed" = "#4C78A8",
      "Exposed" = "#E58B65"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(add = c(0.02, 0.02))
  ) +
  labs(
    title = "After overlap weighting",
    x = "Propensity score",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 13,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_blank(),
    axis.text = element_text(color = "black")
  )


# ============================================================
# Combine the two panels and use a shared legend
# ============================================================

p_combined <- p_before + p_after +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

p_combined

ggsave(
  file.path(figure_dir, "propensity_score_distribution_ow.pdf"),
  plot = p_combined,
  width = 11,
  height = 5,
  units = "in",
  device = cairo_pdf
)

# ============================================================
# 7. Calculate the Jensen-Shannon divergence (JSD)
#    Kernel densities are estimated on the same evaluation grid.
# ============================================================
# ============================================================

calc_jsd <- function(x1, x2, w1 = NULL, w2 = NULL,
                     from = 0, to = 1, n = 512, adjust = 1) {
  
  # Assign equal weights when no weights are supplied
  if (is.null(w1)) {
    w1 <- rep(1, length(x1))
  }
  
  if (is.null(w2)) {
    w2 <- rep(1, length(x2))
  }
  
  # Normalize the weights within each group
  w1 <- w1 / sum(w1)
  w2 <- w2 / sum(w2)
  
  # Use a common bandwidth for both groups
  pooled_x <- c(x1, x2)
  
  bw_common <- density(
    pooled_x,
    from = from,
    to = to,
    n = n,
    adjust = adjust
  )$bw
  
  # Estimate the density for the first group
  d1 <- density(
    x1,
    weights = w1,
    bw = bw_common,
    from = from,
    to = to,
    n = n
  )
  
  # Estimate the density for the second group
  d2 <- density(
    x2,
    weights = w2,
    bw = bw_common,
    from = from,
    to = to,
    n = n
  )
  
  p <- d1$y
  q <- d2$y
  
  # Convert the densities into discrete probability distributions
  p <- p / sum(p)
  q <- q / sum(q)
  
  # Add a small constant to avoid log(0)
  eps <- 1e-12
  
  p <- p + eps
  q <- q + eps
  
  p <- p / sum(p)
  q <- q / sum(q)
  
  # Calculate the midpoint distribution
  m <- (p + q) / 2
  
  # KL divergence
  kl_pm <- sum(p * log2(p / m))
  kl_qm <- sum(q * log2(q / m))
  
  # Jensen-Shannon divergence
  jsd <- 0.5 * kl_pm + 0.5 * kl_qm
  
  return(jsd)
}

# ============================================================
# Define the two exposure groups
# ============================================================

idx_non <- data$exposure_plot == "Unexposed"
idx_exposed <- data$exposure_plot == "Exposed"

# ============================================================
# Calculate the JSD after overlap weighting
# ============================================================

JSD_after <- calc_jsd(
  x1 = data$PS[idx_non],
  x2 = data$PS[idx_exposed],
  
  w1 = data$OW[idx_non],
  w2 = data$OW[idx_exposed]
)

JSD_after
# ============================================================
# 8. Fit the overlap-weighted Cox proportional hazards model
# ============================================================

# Add overlap weights to the analysis data
data$OW <- w_ow$weights

# Define the overlap-weighted survey design
design_ow <- svydesign(
  ids = ~1,
  weights = ~OW,
  data = data
)

# Fit the overlap-weighted Cox model
fit_ow <- svycoxph(
  Surv(follow_up_time, event) ~ exposure,
  design = design_ow
)

summary(fit_ow)
HR <- exp(coef(fit_ow))

# 95% CI
CI <- exp(confint(fit_ow))

# Extract the P value
P <- summary(fit_ow)$coefficients[, "Pr(>|z|)"]

result_ow <- data.frame(
  HR = HR,
  Lower95CI = CI[, 1],
  Upper95CI = CI[, 2],
  P_value = P
)

print(result_ow)
# ============================================================
# 9. Create the overlap-weighted baseline characteristics table
# ============================================================

# ============================================================
# Specify variables included in the baseline table
# ============================================================

vars_table <- c(
  "covariate_01",
  "covariate_02",
  "covariate_03",
  "covariate_04",
  "covariate_05",
  "covariate_13",
  "covariate_14",
  "covariate_15",
  "covariate_16",
  "covariate_06",
  "covariate_07",
  "covariate_08",
  "covariate_09",
  "covariate_10",
  "covariate_11",
  "covariate_12"
)

# Identify categorical variables
factor_vars_table <- c(
  "covariate_02",
  "covariate_03",
  "covariate_04",
  "covariate_16",
  "covariate_08",
  "covariate_12"
)

# ============================================================
# Convert categorical variables to factors
# ============================================================

data[factor_vars_table] <- lapply(
  data[factor_vars_table],
  as.factor
)

data$exposure <- as.factor(data$exposure)


# ============================================================
# Add overlap weights to the analysis data
# ============================================================

data$OW <- w_ow$weights


# ============================================================
# Define the overlap-weighted survey design
# ============================================================

design_ow <- svydesign(
  ids = ~1,
  weights = ~OW,
  data = data
)


# ============================================================
# Create the baseline characteristics table after overlap weighting
# test = TRUE includes between-group P values.
# ============================================================

table_after_ow <- svyCreateTableOne(
  vars = vars_table,
  strata = "exposure",
  data = design_ow,
  factorVars = factor_vars_table,
  test = TRUE
)


# ============================================================
# Specify non-normally distributed continuous variables
#
# covariate_05 and covariate_10 are summarized using mean (SD).
# The remaining continuous variables are summarized using median (IQR).
# ============================================================

nonnormal_vars <- c(
  "covariate_01",
  "covariate_13",
  "covariate_14",
  "covariate_15",
  "covariate_06",
  "covariate_07",
  "covariate_09",
  "covariate_11"
)

# ============================================================
# Calculate the weighted median and IQR for a variable in each group.
# Median (IQR)
# IQR = Q3 - Q1
# ============================================================

get_weighted_median_iqr <- function(var, group) {
  
  # Subset the survey design to the selected group
  design_sub <- subset(
    design_ow,
    exposure == group
  )
  
  # Calculate the weighted first quartile, median, and third quartile
  q <- svyquantile(
    as.formula(paste0("~", var)),
    design = design_sub,
    quantiles = c(0.25, 0.50, 0.75),
    ci = FALSE,
    na.rm = TRUE
  )
  
  q <- as.numeric(q)
  
  Q1 <- q[1]
  Median <- q[2]
  Q3 <- q[3]
  
  IQR <- Q3 - Q1
  
  # Return the result in Median (IQR) format
  sprintf(
    "%.2f (%.2f)",
    Median,
    IQR
  )
}

# ============================================================
# Generate the final formatted baseline table
#
# smd = FALSE excludes SMDs from this exported table.
# test = TRUE retains the P values.
# ============================================================

df_after_ow <- print(
  table_after_ow,
  nonnormal = nonnormal_vars,
  smd = FALSE,
  test = TRUE,
  printToggle = FALSE,
  noSpaces = TRUE,
  quote = FALSE
)

df_after_ow


# ============================================================
# Export the weighted baseline table to Excel
# ============================================================

wb <- createWorkbook()

addWorksheet(
  wb,
  "After_OW"
)

writeData(
  wb,
  sheet = "After_OW",
  x = df_after_ow,
  rowNames = TRUE
)

# Format the table header
headerStyle <- createStyle(
  textDecoration = "bold",
  border = "Bottom"
)

addStyle(
  wb,
  sheet = "After_OW",
  style = headerStyle,
  rows = 1,
  cols = 1:(ncol(df_after_ow) + 1),
  gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "After_OW",
  cols = 1:(ncol(df_after_ow) + 1),
  widths = "auto"
)

saveWorkbook(
  wb,
  file.path(table_dir, "overlap_weighted_baseline_table.xlsx"),
  overwrite = TRUE
)

cat(
  "Overlap-weighted baseline table saved to:",
  file.path(table_dir, "overlap_weighted_baseline_table.xlsx"),
  "\n"
)
