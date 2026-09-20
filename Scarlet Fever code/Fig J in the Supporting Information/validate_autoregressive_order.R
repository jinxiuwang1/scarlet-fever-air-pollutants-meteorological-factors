library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

output_dir <- "E:/Scarlet Fever/validate_autoregressive_order"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data_file <- "E:/scarletfever_2013-2020_2.csv"

data <- read.csv(
  data_file,
  fileEncoding = "GB18030"
)

is_empty_column <- function(x) {
  x_chr <- trimws(as.character(x))
  all(is.na(x_chr) | x_chr == "")
}

data <- data[, !sapply(data, is_empty_column)]

pollution_vars <- c("SO2", "CO", "NO2", "O3_8h")
meteo_vars <- c("rain", "sunlight", "humi", "meantemp")

required_vars <- c(
  "name", "year", "month",
  "beta", "population",
  pollution_vars,
  meteo_vars
)

missing_vars <- setdiff(required_vars, names(data))

if (length(missing_vars) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_vars, collapse = ", ")
  )
}

data <- data[, required_vars]

missing_tokens <- c("NA", "N/A", "NaN", "missing", "", " ")

count_missing_like <- function(x) {
  if (is.character(x) || is.factor(x)) {
    x_chr <- trimws(as.character(x))
    return(sum(is.na(x_chr) | x_chr %in% missing_tokens))
  }
  
  return(sum(is.na(x)))
}

numeric_vars <- c(
  "year", "month",
  "beta", "population",
  pollution_vars,
  meteo_vars
)

na_conversion_log <- data.frame(
  variable = numeric_vars,
  NA_before_conversion = NA_integer_,
  NA_after_conversion = NA_integer_,
  newly_created_NA = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(numeric_vars)) {
  
  var <- numeric_vars[i]
  na_before <- count_missing_like(data[[var]])
  
  if (!is.numeric(data[[var]])) {
    x_chr <- trimws(as.character(data[[var]]))
    x_chr[x_chr %in% missing_tokens] <- NA
    x_chr <- gsub(",", "", x_chr)
    data[[var]] <- suppressWarnings(as.numeric(x_chr))
  }
  
  na_after <- sum(is.na(data[[var]]))
  
  na_conversion_log$NA_before_conversion[i] <- na_before
  na_conversion_log$NA_after_conversion[i] <- na_after
  na_conversion_log$newly_created_NA[i] <- na_after - na_before
}

cat("\nNA counts before and after numeric conversion:\n")
print(na_conversion_log)

write.csv(
  na_conversion_log,
  file.path(output_dir, "NA_conversion_log.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

if (any(is.na(data$beta))) {
  stop("beta has missing values. beta should not be imputed.")
}

if (any(data$beta < 0, na.rm = TRUE)) {
  stop("beta has negative values. The current quasipoisson log-link model is not suitable.")
}

if (any(is.na(data$population))) {
  stop("population has missing values and cannot be used in offset(log(population)).")
}

if (any(data$population <= 0, na.rm = TRUE)) {
  stop("population has zero or negative values and cannot be used in log(population).")
}

data <- data[with(data, order(name, year, month)), ]
rownames(data) <- NULL

data <- data %>%
  group_by(name) %>%
  arrange(year, month, .by_group = TRUE) %>%
  mutate(seq = row_number()) %>%
  ungroup()

safe_na_kalman <- function(x, global_median) {
  
  x <- as.numeric(x)
  
  if (!any(is.na(x))) {
    return(x)
  }
  
  if (all(is.na(x))) {
    return(rep(global_median, length(x)))
  }
  
  if (sum(!is.na(x)) < 3) {
    local_median <- median(x, na.rm = TRUE)
    
    if (!is.finite(local_median)) {
      local_median <- global_median
    }
    
    x[is.na(x)] <- local_median
    return(x)
  }
  
  result <- tryCatch(
    {
      na_kalman(
        x,
        model = "StructTS",
        smooth = TRUE
      )
    },
    error = function(e) {
      tryCatch(
        {
          na_interpolation(x, option = "linear")
        },
        error = function(e2) {
          local_median <- median(x, na.rm = TRUE)
          
          if (!is.finite(local_median)) {
            local_median <- global_median
          }
          
          x[is.na(x)] <- local_median
          return(x)
        }
      )
    }
  )
  
  if (any(is.na(result))) {
    result[is.na(result)] <- global_median
  }
  
  return(result)
}

cat("\nMissing values before meteorological imputation:\n")
print(sapply(data[meteo_vars], function(x) sum(is.na(x))))

for (var in meteo_vars) {
  
  new_var <- paste0(var, "_new")
  global_median <- median(data[[var]], na.rm = TRUE)
  
  if (!is.finite(global_median)) {
    stop(var, " is completely missing and cannot be imputed.")
  }
  
  data[[new_var]] <- ave(
    data[[var]],
    data$name,
    FUN = function(x) {
      safe_na_kalman(x, global_median)
    }
  )
}

meteo_new_vars <- paste0(meteo_vars, "_new")

cat("\nMissing values after meteorological imputation:\n")
print(sapply(data[meteo_new_vars], function(x) sum(is.na(x))))

write.csv(
  data,
  file.path(output_dir, "data_after_meteorological_imputation.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

get_valid_knots <- function(x) {
  
  x_valid <- x[is.finite(x)]
  
  if (length(unique(x_valid)) < 3) {
    return(numeric(0))
  }
  
  x_range <- range(x_valid, na.rm = TRUE)
  
  q <- quantile(
    x_valid,
    probs = c(0.1, 0.5, 0.9),
    na.rm = TRUE,
    names = FALSE
  )
  
  q <- unique(as.numeric(q))
  q <- q[q > x_range[1] & q < x_range[2]]
  
  return(q)
}

create_crossbasis <- function(var_name, df, lag_value = 3) {
  
  x <- df[[var_name]]
  x_valid <- x[is.finite(x)]
  
  if (length(unique(x_valid)) < 2) {
    stop(var_name, " has insufficient variation.")
  }
  
  x_range <- range(x_valid, na.rm = TRUE)
  var_knots <- get_valid_knots(x)
  
  if (length(var_knots) >= 1) {
    argvar_list <- list(
      fun = "ns",
      knots = var_knots,
      Boundary.knots = x_range
    )
  } else {
    argvar_list <- list(
      fun = "lin"
    )
  }
  
  cb <- crossbasis(
    x,
    lag = lag_value,
    argvar = argvar_list,
    arglag = list(
      fun = "ns",
      knots = logknots(lag_value, nk = 2)
    ),
    group = df$name
  )
  
  return(cb)
}

max_dlnm_lag <- 3
max_acf_lag <- 12
minimum_region_n <- 24

cb1.so2 <- create_crossbasis("SO2", data, lag_value = max_dlnm_lag)
cb1.co <- create_crossbasis("CO", data, lag_value = max_dlnm_lag)
cb1.no2 <- create_crossbasis("NO2", data, lag_value = max_dlnm_lag)
cb1.o3 <- create_crossbasis("O3_8h", data, lag_value = max_dlnm_lag)

cb1.sun <- create_crossbasis("sunlight_new", data, lag_value = max_dlnm_lag)
cb1.humi <- create_crossbasis("humi_new", data, lag_value = max_dlnm_lag)
cb1.rain <- create_crossbasis("rain_new", data, lag_value = max_dlnm_lag)
cb1.temp <- create_crossbasis("meantemp_new", data, lag_value = max_dlnm_lag)

model_no_ar <- glm(
  beta ~
    cb1.so2 + cb1.co + cb1.no2 + cb1.o3 +
    cb1.sun + cb1.humi + cb1.rain + cb1.temp +
    ns(seq, 8) +
    offset(log(population)) +
    factor(name),
  family = quasipoisson(),
  data = data,
  na.action = na.exclude
)

capture.output(
  summary(model_no_ar),
  file = file.path(output_dir, "model_without_AR_summary.txt")
)

resid_no_ar <- as.numeric(
  residuals(model_no_ar, type = "deviance")
)

if (length(resid_no_ar) == nrow(data)) {
  data$resid_no_ar <- resid_no_ar
} else {
  model_rows <- as.integer(rownames(model.frame(model_no_ar)))
  data$resid_no_ar <- NA_real_
  data$resid_no_ar[model_rows] <- resid_no_ar
}

compute_region_acf_pacf <- function(df, resid_var, max_lag, min_n) {
  
  region_list <- unique(df$name)
  out <- data.frame()
  
  for (reg in region_list) {
    
    sub <- df[df$name == reg, ]
    sub <- sub[with(sub, order(year, month)), ]
    
    r <- sub[[resid_var]]
    r <- r[!is.na(r)]
    
    if (length(r) < min_n || length(r) <= max_lag + 1 || sd(r) == 0) {
      next
    }
    
    acf_obj <- acf(
      r,
      lag.max = max_lag,
      plot = FALSE
    )
    
    pacf_obj <- pacf(
      r,
      lag.max = max_lag,
      plot = FALSE
    )
    
    acf_values <- as.numeric(acf_obj$acf)[-1]
    pacf_values <- as.numeric(pacf_obj$acf)
    ci_value <- qnorm(0.975) / sqrt(length(r))
    
    one_region <- data.frame(
      name = reg,
      lag = seq_len(max_lag),
      acf = acf_values[seq_len(max_lag)],
      pacf = pacf_values[seq_len(max_lag)],
      ci = ci_value,
      n = length(r),
      stringsAsFactors = FALSE
    )
    
    out <- rbind(out, one_region)
  }
  
  return(out)
}

region_acf_pacf <- compute_region_acf_pacf(
  df = data,
  resid_var = "resid_no_ar",
  max_lag = max_acf_lag,
  min_n = minimum_region_n
)

if (nrow(region_acf_pacf) == 0) {
  stop("No region has enough residual observations for ACF/PACF analysis.")
}

acf_pacf_summary <- region_acf_pacf %>%
  group_by(lag) %>%
  summarise(
    n_regions = n(),
    mean_acf = mean(acf, na.rm = TRUE),
    median_acf = median(acf, na.rm = TRUE),
    mean_pacf = mean(pacf, na.rm = TRUE),
    median_pacf = median(pacf, na.rm = TRUE),
    mean_ci = mean(ci, na.rm = TRUE),
    proportion_significant_acf = mean(abs(acf) > ci, na.rm = TRUE),
    proportion_significant_pacf = mean(abs(pacf) > ci, na.rm = TRUE),
    .groups = "drop"
  )

first_significant_lag <- region_acf_pacf %>%
  group_by(name) %>%
  summarise(
    first_significant_acf_lag = {
      sig_lags <- lag[abs(acf) > ci]
      
      if (length(sig_lags) == 0) {
        NA_integer_
      } else {
        min(sig_lags)
      }
    },
    first_significant_pacf_lag = {
      sig_lags <- lag[abs(pacf) > ci]
      
      if (length(sig_lags) == 0) {
        NA_integer_
      } else {
        min(sig_lags)
      }
    },
    .groups = "drop"
  )

write.csv(
  region_acf_pacf,
  file.path(output_dir, "region_ACF_PACF_values.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

write.csv(
  acf_pacf_summary,
  file.path(output_dir, "ACF_PACF_summary_by_lag.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

write.csv(
  first_significant_lag,
  file.path(output_dir, "first_significant_lag_by_region.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

plot_acf_pacf_summary <- function(summary_df) {
  
  ci_plot <- median(summary_df$mean_ci, na.rm = TRUE)
  
  y_limit <- range(
    c(
      summary_df$median_acf,
      summary_df$median_pacf,
      -ci_plot,
      ci_plot
    ),
    na.rm = TRUE
  )
  
  y_abs <- max(abs(y_limit), na.rm = TRUE)
  y_limit <- c(-y_abs, y_abs)
  
  par(
    mfrow = c(1, 2),
    mar = c(5, 5, 4, 2)
  )
  
  plot(
    summary_df$lag,
    summary_df$median_acf,
    type = "h",
    lwd = 4,
    xlab = "Lag (months)",
    ylab = "Median ACF",
    main = "Residual ACF by region",
    ylim = y_limit
  )
  
  points(
    summary_df$lag,
    summary_df$median_acf,
    pch = 16
  )
  
  abline(h = 0, lwd = 1)
  abline(h = c(-ci_plot, ci_plot), lty = 2)
  
  plot(
    summary_df$lag,
    summary_df$median_pacf,
    type = "h",
    lwd = 4,
    xlab = "Lag (months)",
    ylab = "Median PACF",
    main = "Residual PACF by region",
    ylim = y_limit
  )
  
  points(
    summary_df$lag,
    summary_df$median_pacf,
    pch = 16
  )
  
  abline(h = 0, lwd = 1)
  abline(h = c(-ci_plot, ci_plot), lty = 2)
}

tiff(
  file = file.path(output_dir, "ACF_PACF_by_region.tiff"),
  width = 3000,
  height = 1500,
  res = 300
)

plot_acf_pacf_summary(acf_pacf_summary)

dev.off()

setEPS()

postscript(
  file = file.path(output_dir, "ACF_PACF_by_region.eps"),
  width = 10,
  height = 5,
  horizontal = FALSE,
  onefile = FALSE,
  paper = "special"
)

plot_acf_pacf_summary(acf_pacf_summary)

dev.off()

acf_first_lag_table <- table(
  first_significant_lag$first_significant_acf_lag,
  useNA = "ifany"
)

pacf_first_lag_table <- table(
  first_significant_lag$first_significant_pacf_lag,
  useNA = "ifany"
)

sink(
  file.path(output_dir, "ACF_PACF_results.txt")
)

cat("ACF/PACF analysis for autoregressive order selection\n")
cat("====================================================\n\n")

cat("Model used for residual extraction:\n")
cat("Multivariable DLNM without autoregressive term\n")
cat("beta ~ cb1.so2 + cb1.co + cb1.no2 + cb1.o3 + cb1.sun + cb1.humi + cb1.rain + cb1.temp + ns(seq, 8) + offset(log(population)) + factor(name)\n\n")

cat("Number of regions included in region-specific ACF/PACF analysis:\n")
cat(length(unique(region_acf_pacf$name)), "\n\n")

cat("Distribution of the first significant ACF lag by region:\n")
print(acf_first_lag_table)

cat("\nDistribution of the first significant PACF lag by region:\n")
print(pacf_first_lag_table)

cat("\nSummary by lag:\n")
print(acf_pacf_summary)

cat("\nInterpretation:\n")
cat("Residual ACF/PACF was calculated within each region based on residuals from the multivariable DLNM without the autoregressive term.\n")
cat("If lag 1 shows the most consistent autocorrelation signal across regions and higher lags are less consistent, including a first-order autoregressive term is reasonable.\n")
cat("The final conclusion should be based on the ACF/PACF summary table and the distribution of the first significant lag across regions.\n")

sink()

cat("\nAnalysis completed.\n")
cat("TIFF figure saved to: ",
    file.path(output_dir, "ACF_PACF_by_region.tiff"),
    "\n",
    sep = "")
cat("EPS figure saved to: ",
    file.path(output_dir, "ACF_PACF_by_region.eps"),
    "\n",
    sep = "")
cat("Results saved to: ",
    file.path(output_dir, "ACF_PACF_results.txt"),
    "\n",
    sep = "")

print(acf_first_lag_table)