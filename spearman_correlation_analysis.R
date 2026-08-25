# ============================================================
# Spearman Correlation Analysis
#
# Description:
#   Calculate pairwise Spearman correlation coefficients and P values
#   among anonymized variables. The script creates a correlation heatmap
#   and selected scatterplots involving covariate_17.
#
# Required input columns:
#   covariate_XX
# ============================================================

# Install the packages before first use if necessary:
# install.packages(c("ggplot2", "corrplot"))
library(ggplot2)
library(corrplot)


# ============================================================
# 1. File paths
# ============================================================

# Run this script from the root directory of the GitHub repository.
input_file <- file.path("data", "analysis_data.csv")
figure_dir <- file.path("results", "figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read and prepare the analysis data
# ============================================================

data_cor <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove accidental leading or trailing spaces from column names
names(data_cor) <- trimws(names(data_cor))

# Variables included in the correlation matrix
correlation_vars <- c(
  "covariate_XX
)

cor_data <- data_cor[, correlation_vars]


# ============================================================
# 3. Calculate Spearman correlation and P-value matrices
# ============================================================

# Pairwise complete observations are used for correlations involving
# variables with missing values.
cor_mat <- cor(
  cor_data,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Calculate a P value for every pair of variables
cor_test_matrix <- function(data) {
  n_vars <- ncol(data)

  p_mat <- matrix(
    NA_real_,
    nrow = n_vars,
    ncol = n_vars,
    dimnames = list(colnames(data), colnames(data))
  )

  for (i in seq_len(n_vars)) {
    for (j in seq_len(n_vars)) {
      if (i == j) {
        p_mat[i, j] <- 0
        next
      }

      pair_data <- data[, c(i, j)]
      pair_data <- pair_data[complete.cases(pair_data), , drop = FALSE]

      test_result <- cor.test(
        pair_data[[1]],
        pair_data[[2]],
        method = "spearman",
        exact = FALSE
      )

      p_mat[i, j] <- test_result$p.value
    }
  }

  p_mat
}

p_mat <- cor_test_matrix(cor_data)

# ============================================================
# 4. Create the correlation heatmap
# ============================================================

# Convert P values into conventional significance symbols
star_mat <- matrix(
  "",
  nrow = nrow(p_mat),
  ncol = ncol(p_mat),
  dimnames = dimnames(p_mat)
)

star_mat[p_mat < 0.05] <- "*"
star_mat[p_mat < 0.01] <- "**"
star_mat[p_mat < 0.001] <- "***"
diag(star_mat) <- ""

# Combine each Spearman coefficient with its significance symbol
label_mat <- matrix(
  "",
  nrow = nrow(cor_mat),
  ncol = ncol(cor_mat),
  dimnames = dimnames(cor_mat)
)

for (i in seq_len(nrow(cor_mat))) {
  for (j in seq_len(ncol(cor_mat))) {
    label_mat[i, j] <- paste0(
      sprintf("%.2f", cor_mat[i, j]),
      star_mat[i, j]
    )
  }
}

corrplot(
  cor_mat,
  method = "color",
  type = "full",
  diag = TRUE,
  order = "original",
  col = colorRampPalette(c("#3B6FB6", "white", "#D95F59"))(200),
  tl.col = "black",
  tl.cex = 0.9,
  tl.srt = 45,
  cl.lim = c(-1, 1),
  mar = c(0, 0, 1, 0)
)

# Add coefficients and significance symbols at the center of each cell
n_vars <- ncol(cor_mat)

for (i in seq_len(n_vars)) {
  for (j in seq_len(n_vars)) {
    text(
      x = j,
      y = n_vars - i + 1,
      labels = label_mat[i, j],
      cex = 0.65,
      col = "black"
    )
  }
}


# ============================================================
# 5. Create selected correlation scatterplots
# ============================================================

# This helper function calculates Spearman's rho and creates a scatterplot
# with a LOESS smooth and 95% confidence band.
create_correlation_plot <- function(data, outcome, y_label, output_name) {
  test_result <- cor.test(
    data$covariate_17,
    data[[outcome]],
    method = "spearman",
    exact = FALSE
  )

  rho_value <- unname(test_result$estimate)
  p_value <- test_result$p.value

  p_text <- if (p_value < 0.001) {
    paste0("Spearman ρ = ", sprintf("%.2f", rho_value), ", P < 0.001")
  } else {
    paste0(
      "Spearman ρ = ", sprintf("%.2f", rho_value),
      ", P = ", sprintf("%.3f", p_value)
    )
  }

  scatter_plot <- ggplot(
    data,
    aes(x = covariate_17, y = .data[[outcome]])
  ) +
    geom_point(
      size = 2.3,
      alpha = 0.65
    ) +
    geom_smooth(
      method = "loess",
      se = TRUE,
      linewidth = 0.9
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = p_text,
      hjust = 1.08,
      vjust = 1.5,
      size = 4.2
    ) +
    labs(
      x = "covariate_XX",
      y = y_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(size = 13),
      plot.margin = margin(12, 18, 12, 12)
    )

  print(scatter_plot)

  ggsave(
    filename = file.path(figure_dir, output_name),
    plot = scatter_plot,
    width = 6,
    height = 5,
    units = "in",
    device = cairo_pdf
  )

  invisible(scatter_plot)
}

# Generate the prespecified covariate_17 association plots
p_covariate_10 <- create_correlation_plot(
  data_cor, "covariate_XX", "covariate_XX", "covariate_XX_vs_covariate_XX.pdf"
)

