# ============================================================
# Optimism-Corrected C-index Dumbbell Plot
#
# Description:
#   Compare the optimism-corrected C-index of a base model and an
#   extended model containing an anonymized exposure variable.
#
# Required input columns:
#   Metric, Value
# ============================================================

library(ggplot2)


# ============================================================
# 1. File paths
# ============================================================

input_file <- file.path("results", "tables", "model_increment_results.csv")
output_dir <- file.path("results", "figures")
output_file <- file.path(output_dir, "corrected_cindex_dumbbell.pdf")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read the model performance results
# ============================================================

cindex_result <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

print(cindex_result)


# ============================================================
# 3. Extract C-index results
# ============================================================

C_base <- cindex_result$Value[
  cindex_result$Metric == "Base optimism-corrected C-index"
]

C_extended <- cindex_result$Value[
  cindex_result$Metric ==
    "Base + exposure optimism-corrected C-index"
]

Delta_C <- cindex_result$Value[
  cindex_result$Metric ==
    "Optimism-corrected Delta C-index"
]

Delta_lower <- cindex_result$Value[
  cindex_result$Metric ==
    "Delta C-index 95% CI lower"
]

Delta_upper <- cindex_result$Value[
  cindex_result$Metric ==
    "Delta C-index 95% CI upper"
]

# Display the extracted values
cat(
  "Base =", C_base,
  "\nBase + exposure =", C_extended,
  "\nDelta C =", Delta_C,
  "\n95% CI =", Delta_lower, "to", Delta_upper,
  "\n"
)


# ============================================================
# 4. Create the dumbbell plot
# ============================================================

p_dumbbell <- ggplot() +

  # Line connecting the two models
  geom_segment(
    aes(
      x = 1,
      xend = 2,
      y = C_base,
      yend = C_extended
    ),
    linewidth = 1.3,
    color = "grey55"
  ) +

  # Base model
  geom_point(
    aes(
      x = 1,
      y = C_base
    ),
    size = 5.5,
    color = "#D95F59"
  ) +

  # Extended model
  geom_point(
    aes(
      x = 2,
      y = C_extended
    ),
    size = 5.5,
    color = "#168AAD"
  ) +

  # Base-model estimate
  annotate(
    "text",
    x = 1,
    y = C_base - 0.012,
    label = sprintf("%.3f", C_base),
    size = 4.5,
    fontface = "bold"
  ) +

  # Extended-model estimate
  annotate(
    "text",
    x = 2,
    y = C_extended + 0.012,
    label = sprintf("%.3f", C_extended),
    size = 4.5,
    fontface = "bold"
  ) +

  # Change in C-index
  annotate(
    "text",
    x = 1.5,
    y = max(C_base, C_extended) + 0.040,
    label = paste0(
      "ΔC-index = ",
      sprintf("%.3f", Delta_C)
    ),
    size = 4.4,
    fontface = "bold"
  ) +

  # Confidence interval for the change in C-index
  annotate(
    "text",
    x = 1.5,
    y = max(C_base, C_extended) + 0.025,
    label = paste0(
      "95% CI ",
      sprintf("%.3f", Delta_lower),
      "–",
      sprintf("%.3f", Delta_upper)
    ),
    size = 4.1
  ) +

  # Horizontal axis
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c(
      "Base model",
      "Base + exposure"
    ),
    limits = c(0.75, 2.25)
  ) +

  # Vertical axis
  scale_y_continuous(
    limits = c(0.55, 0.70),
    breaks = seq(
      0.55,
      0.70,
      by = 0.05
    )
  ) +

  labs(
    x = NULL,
    y = "Optimism-corrected C-index"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 11.5,
      color = "black"
    ),

    axis.title.y = element_text(
      size = 13.5,
      margin = margin(r = 10)
    ),

    axis.line = element_line(
      linewidth = 0.7,
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.6,
      color = "black"
    ),

    plot.margin = margin(
      15, 15, 15, 15
    )
  )

p_dumbbell

ggsave(
  output_file,
  plot = p_dumbbell,
  width = 6.5,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)
