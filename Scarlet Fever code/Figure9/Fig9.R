# Load required packages
library(ggplot2)
library(reshape2)
library(corrplot)
library(Hmisc)
library(imputeTS)          # Kalman filter imputation
library(dplyr)             # Data manipulation

# Read data
data <- read.csv("E:/scarletfever_2013-2020_2.csv", 
                 fileEncoding = "GB18030")

# Define required variable columns
data_columns <- c("PM2.5", "PM10", "SO2", "NO2", "O3_8h", "CO", 
                  "meantemp", "humi", "rain", "sunlight", "beta")

# Keep necessary columns (remove ghost columns)
keep_cols <- c("name", "year", "month", data_columns)
data <- data[, keep_cols]

# Sort by time
data <- data[with(data, order(name, year, month)), ]

# Set working directory
work_dir <- "E:/Scarlet Fever/fig9"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# Check column names
cat("Column names in data:\n")
print(colnames(data))

# Ensure all analysis columns are numeric (ignore warnings, NAs will be imputed later)
for (col in data_columns) {
  if (!is.numeric(data[[col]])) {
    suppressWarnings(data[[col]] <- as.numeric(as.character(data[[col]])))
    cat("Converted column", col, "to numeric type\n")
  }
}

# ========== Kalman filter imputation ==========
# Group by name and apply Kalman filter imputation to each variable's time series
data_interp <- data %>%
  group_by(name) %>%
  arrange(year, month, .by_group = TRUE) %>%
  mutate(across(all_of(data_columns), 
                ~ na_kalman(., model = "StructTS", smooth = TRUE))) %>%
  ungroup()

# Check for remaining missing values after imputation
na_after_interp <- sapply(data_interp[data_columns], function(x) sum(is.na(x)))
cat("\nMissing values after Kalman filter imputation:\n")
print(na_after_interp)

# If any missing values remain, remove them (assume all are filled here)
data_subset <- data_interp[, data_columns]
data_subset <- na.omit(data_subset)   # Safety removal
cat("Rows removed due to remaining missing values:", nrow(data_interp) - nrow(data_subset), "\n")

# Check if subset is empty
if (nrow(data_subset) == 0) {
  stop("Data subset is empty, please check data quality")
}

# ---------- Compute correlation matrix and significance ----------
rcorr_result <- rcorr(as.matrix(data_subset), type = "pearson")
cor_matrix <- rcorr_result$r          # Correlation coefficient matrix
p_matrix  <- rcorr_result$P           # p-value matrix

# Save correlation matrix
write.csv(cor_matrix, 
          file.path(work_dir, "correlation_matrix.csv"), 
          row.names = TRUE, 
          fileEncoding = "GB18030")

# Save p-value matrix
write.csv(p_matrix,
          file.path(work_dir, "p_value_matrix.csv"),
          row.names = TRUE,
          fileEncoding = "GB18030")
cat("\np-value matrix saved to:", file.path(work_dir, "p_value_matrix.csv"), "\n")

# ---------- Prepare data for plotting ----------
melted_cor <- melt(cor_matrix, varnames = c("Var1", "Var2"), value.name = "cor")
melted_p   <- melt(p_matrix,  varnames = c("Var1", "Var2"), value.name = "p")
plot_data <- merge(melted_cor, melted_p, by = c("Var1", "Var2"))

# Create labels: add asterisk for p < 0.01
plot_data$label <- ifelse(
  plot_data$p < 0.01,
  paste0(round(plot_data$cor, 2), "*"),
  as.character(round(plot_data$cor, 2))
)

# No asterisk on diagonal
plot_data$label <- ifelse(
  as.character(plot_data$Var1) == as.character(plot_data$Var2),
  as.character(round(plot_data$cor, 2)),
  plot_data$label
)

# ---------- Custom label function for axes ----------
create_labels <- function(labels) {
  label_list <- list(
    "PM2.5" = expression(PM[2.5]),
    "PM10"  = expression(PM[10]),
    "SO2"   = expression(SO[2]),
    "NO2"   = expression(NO[2]),
    "O3_8h" = expression(O[3]),
    "CO"    = "CO",
    "meantemp" = "Mean temperature",
    "humi"     = "Relative humidity",
    "rain"     = "Precipitation",
    "sunlight" = "Sunlight",
    "beta"     = "Transmission rate"
  )
  sapply(labels, function(x) ifelse(x %in% names(label_list), label_list[[x]], x))
}

# ---------- Draw heatmap ----------
p_heatmap <- ggplot(plot_data, aes(x = Var1, y = Var2, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 5, color = "black") +
  scale_fill_gradient2(
    low = "#2b8cbe",
    mid = "white",
    high = "#e41a1c",
    midpoint = 0, 
    limits = c(-1, 1),
    name = ""
  ) +
  scale_x_discrete(labels = create_labels) +
  scale_y_discrete(labels = create_labels) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(x = "", y = "")

# ---------- Save plots ----------
ggsave(
  file.path(work_dir, "correlation_heatmap.eps"),
  plot = p_heatmap,
  device = cairo_ps,
  width = 12, 
  height = 10, 
  dpi = 600
)

ggsave(
  file.path(work_dir, "correlation_heatmap.png"),
  plot = p_heatmap,
  width = 12, 
  height = 10, 
  dpi = 300
)

# Output completion messages
cat("\nCorrelation matrix saved to:", file.path(work_dir, "correlation_matrix.csv"))
cat("\np-value matrix saved to:", file.path(work_dir, "p_value_matrix.csv"))
cat("\nCorrelation heatmap (EPS format) saved to:", file.path(work_dir, "correlation_heatmap.eps"))
cat("\nCorrelation heatmap (PNG format) saved to:", file.path(work_dir, "correlation_heatmap.png"))