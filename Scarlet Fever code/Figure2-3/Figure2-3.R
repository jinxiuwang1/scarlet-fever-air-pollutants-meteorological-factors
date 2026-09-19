# Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(imputeTS)

output_dir <- "E:/Scarlet Fever/Figure2-3"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = TRUE)
  message("Created output directory: ", output_dir)
} else {
  message("Output directory already exists: ", output_dir)
}

data <- read.csv("E:/scarletfever_2013-2020_1.csv", 
                 fileEncoding = "GB18030")

preprocess_data <- function(df) {
  processed <- df %>%
    rename(
      Province = province,
      Year = year,
      Month = month
    ) %>%
    select(Province, Year, Month, 
           meantemp, rain, humi, sunlight,
           SO2, NO2, PM10, CO, O3_8h, PM2.5)
  
  processed$Year <- as.numeric(as.character(processed$Year))
  processed <- processed %>% filter(Year >= 2013 & Year <= 2020)
  
  vars_to_check <- c("meantemp", "rain", "humi", "sunlight",
                     "SO2", "NO2", "PM10", "CO", "O3_8h", "PM2.5")
  
  for (var in vars_to_check) {
    if (!is.numeric(processed[[var]])) {
      char_vals <- as.character(processed[[var]])
      char_vals <- gsub("[^0-9\\.eE-]", "", char_vals)
      char_vals[char_vals == ""] <- NA
      processed[[var]] <- as.numeric(char_vals)
    }
  }
  
  processed <- processed %>%
    group_by(Province) %>%
    mutate(
      across(all_of(vars_to_check), 
             ~ if (all(is.na(.))) {
               rep(median(., na.rm = TRUE), n())
             } else {
               na_kalman(.)
             }
      )
    ) %>%
    ungroup()
  
  return(processed)
}

processed_data <- preprocess_data(data)

create_monthly_boxplot <- function(data, var_name, title_label, unit_label, is_pollutant = FALSE) {
  p <- ggplot(data, aes(x = factor(Month), y = .data[[var_name]])) +
    geom_boxplot(fill = "#4E84C4", color = "#1F3552", alpha = 0.7, 
                 outlier.color = "#D55E00", outlier.shape = 16, outlier.size = 1.5) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
    labs(
      title = title_label,
      x = "Month",
      y = if (is_pollutant) {
        paste0("Concentration ", unit_label)
      } else {
        paste0("Value ", unit_label)
      }
    ) +
    scale_x_discrete(labels = month.abb) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

meteo_vars <- list(
  list(var = "meantemp", title = "Mean Temperature", unit = "(°C)"),
  list(var = "rain", title = "Precipitation", unit = "(mm)"),
  list(var = "humi", title = "Relative Humidity", unit = "(%)"),
  list(var = "sunlight", title = "Sunlight", unit = "(hours)")
)

meteo_plots <- lapply(meteo_vars, function(v) {
  create_monthly_boxplot(processed_data, v$var, v$title, v$unit, is_pollutant = FALSE)
})

combined_meteo <- wrap_plots(meteo_plots, ncol = 2) +
  plot_annotation(
    title = "Monthly Distribution of Meteorological Variables (2013-2020)",
    subtitle = "Data from 31 Provinces in China",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
    )
  )

ggsave(paste0(output_dir, "meteo_monthly_boxplots.eps"),
       plot = combined_meteo,
       device = "eps",
       width = 14,
       height = 10,
       units = "in",
       dpi = 600)

pollutant_vars <- list(
  list(var = "SO2", title = expression(SO[2]), unit = "(μg/m³)"),
  list(var = "NO2", title = expression(NO[2]), unit = "(μg/m³)"),
  list(var = "PM10", title = expression(PM[10]), unit = "(μg/m³)"),
  list(var = "CO", title = expression(CO), unit = "(mg/m³)"),
  list(var = "O3_8h", title = expression(O[3]), unit = "(μg/m³)"),
  list(var = "PM2.5", title = expression(PM[2.5]), unit = "(μg/m³)")
)

pollutant_plots <- lapply(pollutant_vars, function(v) {
  create_monthly_boxplot(processed_data, v$var, v$title, v$unit, is_pollutant = TRUE)
})

combined_pollutants <- wrap_plots(pollutant_plots, ncol = 3) +
  plot_annotation(
    title = "Monthly Distribution of Air Pollutants in China (2013-2020)",
    subtitle = "Data from 31 Provinces",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30")
    )
  ) +
  theme(
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(paste0(output_dir, "pollutant_monthly_boxplots.eps"),
       plot = combined_pollutants,
       device = "eps",
       width = 16,
       height = 12,
       units = "in",
       dpi = 600)

print(combined_meteo)
print(combined_pollutants)

message("Analysis completed! All results saved to: ", output_dir)
message("Generated EPS files:")
message("- meteo_monthly_boxplots.eps")
message("- pollutant_monthly_boxplots.eps")