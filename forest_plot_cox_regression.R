# ============================================================
# Cox Regression Forest Plot
#
# Description:
#   Cox regression analyses were performed using IBM SPSS Statistics.
#   This script imports the Cox regression results exported from SPSS
#   and creates a publication-style forest plot in R.
#
# Required input columns:
#   Characteristic: variable name
#   HR: hazard ratio
#   lower: lower limit of the 95% confidence interval
#   upper: upper limit of the 95% confidence interval
#   P: P value
# ============================================================
# Load required packages
# Install forestploter first if needed: install.packages("forestploter")
library(forestploter)
library(grid)


# ============================================================
# 1. File paths
# ============================================================

# Relative paths are used so that the script can run on different computers.
# Run this script from the root directory of the GitHub repository.
input_file <- file.path("data", "cox_results.csv")
output_dir <- file.path("results", "figures")
output_file <- file.path(output_dir, "cox_forest_plot.pdf")

# Create the output directory if it does not already exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read Cox regression results exported from SPSS
# ============================================================

# Cox regression was performed in IBM SPSS Statistics.
# The resulting HRs, 95% CIs, and P values were organized in a CSV file.
df <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 3. Format values displayed in the forest plot
# ============================================================

# Combine the HR and 95% confidence interval into one text column
df$`HR (95% CI)` <- sprintf(
  "%.2f (%.2f–%.2f)",
  df$HR,
  df$lower,
  df$upper
)

# Format P values to three decimal places; values below 0.001 are shown as <0.001
# The code also preserves values that are already stored in the form "<0.001".
df$P <- as.character(df$P)

df$`P value` <- ifelse(
  grepl("<", df$P),
  df$P,
  ifelse(
    as.numeric(df$P) < 0.001,
    "<0.001",
    sprintf("%.3f", as.numeric(df$P))
  )
)

# forestploter uses this blank column as the plotting area for HRs and 95% CIs.
# Increasing the number of spaces makes the plotting area wider.
df$` ` <- "                    "

# Select and arrange the columns displayed in the final figure
plot_data <- df[, c(
  "Characteristic",
  " ",
  "HR (95% CI)",
  "P value"
)]


# ============================================================
# 4. Define the forest plot theme
# ============================================================

tm <- forest_theme(
  base_size = 11,

  # Point estimates and 95% confidence intervals
  ci_pch = 15,
  ci_col = "black",
  ci_fill = "black",
  ci_alpha = 1,
  ci_lty = 1,
  ci_lwd = 1.2,
  ci_Theight = 0.15,

  # Reference line at HR = 1
  refline_lwd = 1,
  refline_lty = "dashed",
  refline_col = "black",

  # Axis and table-header formatting
  xaxis_cex = 0.9,
  header_gp = gpar(
    fontface = "bold",
    fontsize = 11
  ),

  # Table-body formatting
  core = list(
    bg_params = list(
      fill = "white",
      col = NA
    ),
    fg_params = list(
      fontsize = 10.5
    ),
    padding = unit(c(5, 5), "mm")
  )
)


# ============================================================
# 5. Create the forest plot
# ============================================================

p <- forest(
  plot_data,

  # Cox regression estimates
  est = df$HR,
  lower = df$lower,
  upper = df$upper,

  # Draw the forest plot in the second (blank) table column
  ci_column = 2,
  ref_line = 1,

  # Display HRs on a logarithmic scale
  x_trans = "log",
  xlim = c(0.5, 5),
  ticks_at = c(0.5, 1, 2, 3, 5),
  xlab = "Hazard ratio",
  theme = tm
)

# Display the figure in the active graphics device
plot(p)


# ============================================================
# 6. Export the figure as a PDF
# ============================================================

pdf(
  output_file,
  width = 9,
  height = 4.2,
  family = "Helvetica"
)

plot(p)
dev.off()
