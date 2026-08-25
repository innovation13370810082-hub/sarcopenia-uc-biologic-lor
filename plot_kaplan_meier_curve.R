# ============================================================
# Kaplan-Meier Curve for Loss of Response
#
# Description:
#   Compare LOR-free survival between patients with and without
#   sarcopenia using Kaplan-Meier analysis and the log-rank test.
#
# Required input columns:
#   Sarcopenia: sarcopenia status (0 = no, 1 = yes)
#   LOR: loss-of-response event indicator (0 = censored, 1 = event)
#   Follow-up time: follow-up duration in months
# ============================================================

# Install these packages before first use if necessary:
# install.packages(c("survival", "survminer"))
library(survival)
library(survminer)


# ============================================================
# 1. File paths
# ============================================================

# Run this script from the root directory of the GitHub repository.
input_file <- file.path("data", "analysis_data.csv")
output_dir <- file.path("results", "figures")
output_file <- file.path(output_dir, "kaplan_meier_sarcopenia_lor.pdf")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read and prepare the analysis data
# ============================================================

df <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Assign descriptive labels to the sarcopenia groups
df$Sarcopenia <- factor(
  df$Sarcopenia,
  levels = c(0, 1),
  labels = c("Non-sarcopenia", "Sarcopenia")
)

# Optional checks of group sizes and event counts
table(df$LOR)
table(df$Sarcopenia)
table(df$Sarcopenia, df$LOR)


# ============================================================
# 3. Fit Kaplan-Meier curves
# ============================================================

km_fit <- survfit(
  Surv(`Follow-up time`, LOR) ~ Sarcopenia,
  data = df
)

# Calculate the two-sided log-rank P value
logrank_test <- survdiff(
  Surv(`Follow-up time`, LOR) ~ Sarcopenia,
  data = df
)

logrank_p <- 1 - pchisq(
  logrank_test$chisq,
  df = length(logrank_test$n) - 1
)

# Create a plotmath-compatible label with an italic P
p_label <- if (logrank_p < 0.0001) {
  "italic(P) < 0.0001"
} else {
  paste0("italic(P) == ", format.pval(logrank_p, digits = 3, eps = 0.0001))
}


# ============================================================
# 4. Create the Kaplan-Meier plot and risk table
# ============================================================

km_plot <- ggsurvplot(
  km_fit,
  data = df,

  # Statistical display
  pval = FALSE,
  pval.method = FALSE,
  conf.int = TRUE,

  # Number-at-risk table
  risk.table = TRUE,
  risk.table.height = 0.24,
  risk.table.title = "Number at risk",
  risk.table.y.text.col = FALSE,
  risk.table.y.text = TRUE,

  # Axes
  xlim = c(0, 40),
  break.x.by = 12,
  ylim = c(0, 1),
  break.y.by = 0.2,
  xlab = "Time (months)",
  ylab = "LOR-free survival probability",

  # Legend
  legend.title = NULL,
  legend.labs = c("Non-sarcopenia", "Sarcopenia"),
  legend = c(0.72, 0.86),

  # Survival curves and censoring marks
  size = 1.1,
  censor = TRUE,
  censor.shape = 124,
  censor.size = 3,

  # Figure style
  ggtheme = theme_classic(base_size = 13),
  font.x = c(13, "plain"),
  font.y = c(13, "plain"),
  font.tickslab = c(11),
  font.legend = c(11),
  surv.median.line = "none"
)


# ============================================================
# 5. Refine the main plot and risk table
# ============================================================

# Add the automatically calculated log-rank P value
km_plot$plot <- km_plot$plot +
  annotate(
    "text",
    x = 1.5,
    y = 0.15,
    label = p_label,
    parse = TRUE,
    hjust = 0,
    size = 4.5
  ) +
  theme(
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.6),
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      size = 13,
      margin = margin(r = 8)
    ),
    axis.text = element_text(size = 11),
    legend.background = element_blank(),
    legend.key = element_blank(),
    plot.margin = margin(8, 10, 4, 8)
  )

# Refine the appearance of the number-at-risk table
km_plot$table <- km_plot$table +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_text(
      size = 12,
      margin = margin(t = 6)
    ),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(
      size = 12,
      face = "plain",
      hjust = 0
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(0, 10, 5, 8)
  )

# Display the combined Kaplan-Meier plot and risk table
km_plot


# ============================================================
# 6. Export the figure as a PDF
# ============================================================

pdf(
  file = output_file,
  width = 7,
  height = 6.5,
  useDingbats = FALSE
)

print(km_plot)
dev.off()
