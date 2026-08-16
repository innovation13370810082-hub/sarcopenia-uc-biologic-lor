# ============================================================
# Univariable Cox Regression Forest Plot
# SCI-style forest plot
# ============================================================
install.packages("forestploter")
# 1. 加载R包
library(ggplot2)
library(dplyr)
library(patchwork)
library(forestploter)
library(grid)

# 2. 设置工作目录
setwd("D:/身体成分分析IBD/分析")

# 3. 读取CSV
df <- read.csv(
  "单因素Cox回归.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

View(df)

# ============================================================
# 3. 数据整理
# ============================================================

# HR和95%CI
df$`HR (95% CI)` <- sprintf(
  "%.2f (%.2f–%.2f)",
  df$HR,
  df$lower,
  df$upper
)

# P值
# 推荐CSV中的P尽量存成真实数值
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

# forestploter需要一个空白列来放森林图
# 空格越多，森林图区越宽
df$` ` <- "                    "


# ============================================================
# 4. 创建用于显示的表格
# ============================================================

plot_data <- df[, c(
  "Characteristic",
  " ",
  "HR (95% CI)",
  "P value"
)]


# ============================================================
# 5. 设置森林图主题
# ============================================================

tm <- forest_theme(
  
  base_size = 11,
  
  # HR点和95%CI
  ci_pch = 15,
  ci_col = "black",
  ci_fill = "black",
  ci_alpha = 1,
  ci_lty = 1,
  ci_lwd = 1.2,
  ci_Theight = 0.15,
  
  # HR = 1参考线
  refline_lwd = 1,
  refline_lty = "dashed",
  refline_col = "black",
  
  # X轴
  xaxis_cex = 0.9,
  
  # 表头
  header_gp = gpar(
    fontface = "bold",
    fontsize = 11
  ),
  
  # 正文：全部白色背景
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
# 6. 绘制森林图
# ============================================================

p <- forest(
  
  plot_data,
  
  # HR和CI
  est = df$HR,
  lower = df$lower,
  upper = df$upper,
  
  # 第2列为空白列，因此森林图画在第2列
  ci_column = 2,
  
  # HR = 1
  ref_line = 1,
  
  # 对数坐标
  x_trans = "log",
  
  # X轴范围
  xlim = c(0.5, 5),
  
  # 刻度
  ticks_at = c(0.5, 1, 2, 3, 5),
  
  # X轴标题
  xlab = "Hazard ratio",
  
  # 主题
  theme = tm
)


# ============================================================
# 7. 显示
# ============================================================

plot(p)


# ============================================================
# 8. 导出PDF
# ============================================================

pdf(
  "Univariable_Cox_Forest_Plot2.pdf",
  width = 9,
  height = 4.2,
  family = "Helvetica"
)

plot(p)

dev.off()