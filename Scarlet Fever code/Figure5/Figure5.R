library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster)
library(factoextra)
library(patchwork)
library(stringr)
library(htmlwidgets)

output_dir <- "E:/Scarlet Fever/Figure5"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
  message("Created output directory: ", output_dir)
} else {
  message("Output directory already exists: ", output_dir)
}

data <- read.csv("E:/scarletfever_2013-2020_1.csv", 
                 fileEncoding = "GB18030")

data <- data %>%
  rename(
    Province_EN = province,
    Province_CN = name,
    Month = month,
    Year = year,
    Beta = beta
  )

province_to_region <- function(province) {
  case_when(
    province %in% c("Beijing", "Tianjin", "Hebei", "Shanxi", "Inner Mongolia") ~ "North China",
    province %in% c("Liaoning", "Jilin", "Heilongjiang") ~ "Northeast",
    province %in% c("Shanghai", "Jiangsu", "Zhejiang", "Anhui", "Fujian", "Jiangxi", "Shandong") ~ "East China",
    province %in% c("Henan", "Hubei", "Hunan") ~ "Central China",
    province %in% c("Guangdong", "Guangxi", "Hainan") ~ "South China",
    province %in% c("Chongqing", "Sichuan", "Guizhou", "Yunnan", "Tibet") ~ "Southwest",
    province %in% c("Shaanxi", "Gansu", "Qinghai", "Ningxia", "Xinjiang") ~ "Northwest"
  )
}

skewness <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  sqrt(n) * sum(x^3) / (sum(x^2)^(3/2))
}

kurtosis <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  n * sum(x^4) / (sum(x^2)^2) - 3
}

processed_data <- data %>%
  mutate(Year = as.numeric(as.character(Year))) %>%
  filter(Year >= 2013 & Year <= 2020) %>%
  mutate(Month = as.integer(Month)) %>%
  mutate(Month_str = sprintf("%02d", Month)) %>%
  mutate(Date = as.Date(paste(Year, Month_str, "01", sep = "-"))) %>%
  mutate(Beta = ifelse(is.na(Beta), median(Beta, na.rm = TRUE), Beta)) %>%
  mutate(Region = province_to_region(Province_EN))

province_order <- processed_data %>%
  group_by(Province_EN) %>%
  summarise(Median_Beta = median(Beta)) %>%
  arrange(Median_Beta) %>%
  pull(Province_EN)

processed_data$Province_EN <- factor(processed_data$Province_EN, levels = province_order)

province_stats <- processed_data %>%
  group_by(Province_EN, Province_CN, Region) %>%
  summarise(
    Observations = n(),
    Mean = mean(Beta),
    Median = median(Beta),
    SD = sd(Beta),
    IQR = IQR(Beta),
    Min = min(Beta),
    Max = max(Beta),
    Skewness = skewness(Beta),
    Kurtosis = kurtosis(Beta),
    .groups = 'drop'
  ) %>%
  arrange(desc(Median))

write.csv(province_stats, paste0(output_dir, "beta_province_statistics.csv"), row.names = FALSE)

beta_violin <- ggplot(processed_data, aes(x = Province_EN, y = Beta, fill = Region)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, size = 0.8, shape = 16) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "",
       x = "Province",
       y = "Transmission Rate (Beta)",
       fill = "Region") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = "bottom"
  ) +
  coord_flip()

ggsave(paste0(output_dir, "beta_violin.png"), 
       plot = beta_violin,
       device = "png",
       width = 14,
       height = 10,
       units = "in",
       dpi = 600)

ggsave(paste0(output_dir, "beta_violin.eps"), 
       plot = beta_violin,
       device = "eps",
       width = 14,
       height = 10,
       units = "in",
       dpi = 600)

message("Analysis completed! All results saved to: ", output_dir)
message("Generated files:")
message("- beta_violin.png (high-resolution PNG)")
message("- beta_violin.eps (vector EPS)")
message("- beta_province_statistics.csv (descriptive statistics)")