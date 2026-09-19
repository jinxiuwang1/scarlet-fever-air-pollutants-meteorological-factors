library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

output_dir <- "E:/Scarlet Fever/multivariable_model_diagnostics"

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
  stop("beta has missing values. beta should not be imputed for model diagnostics.")
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
  mutate(
    lag.value1 = dplyr::lag(beta, n = 1, default = NA),
    seq = row_number()
  ) %>%
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
  
  new_var <- paste0(var, "_imp")
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

meteo_imp_vars <- paste0(meteo_vars, "_imp")

cat("\nMissing values after meteorological imputation:\n")
print(sapply(data[meteo_imp_vars], function(x) sum(is.na(x))))

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

calc_pseudo_r2 <- function(model) {
  
  if (is.null(model$deviance) || is.null(model$null.deviance)) {
    return(NA_real_)
  }
  
  if (!is.finite(model$null.deviance) || model$null.deviance == 0) {
    return(NA_real_)
  }
  
  return(1 - model$deviance / model$null.deviance)
}

calc_lag1_correlation_by_region <- function(df, resid_var) {
  
  out <- df %>%
    group_by(name) %>%
    arrange(year, month, .by_group = TRUE) %>%
    summarise(
      n_residuals = sum(!is.na(.data[[resid_var]])),
      lag1_correlation = {
        r <- .data[[resid_var]]
        r <- r[!is.na(r)]
        
        if (length(r) < 3 || sd(r) == 0) {
          NA_real_
        } else {
          cor(r[-length(r)], r[-1])
        }
      },
      .groups = "drop"
    )
  
  return(out)
}

summarise_region_acf <- function(region_acf_df, col_name) {
  
  x <- region_acf_df[[col_name]]
  x <- x[is.finite(x)]
  
  data.frame(
    n_regions = length(x),
    mean_lag1_correlation = mean(x, na.rm = TRUE),
    median_lag1_correlation = median(x, na.rm = TRUE),
    mean_abs_lag1_correlation = mean(abs(x), na.rm = TRUE),
    median_abs_lag1_correlation = median(abs(x), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

analysis_vars <- c(
  "SO2",
  "CO",
  "NO2",
  "O3_8h",
  "sunlight_imp",
  "humi_imp",
  "rain_imp",
  "meantemp_imp"
)

cb_names <- setNames(
  paste0("cb_", make.names(analysis_vars)),
  analysis_vars
)

model_env <- new.env(parent = globalenv())

for (var in analysis_vars) {
  
  cb_obj <- create_crossbasis(
    var_name = var,
    df = data,
    lag_value = 3
  )
  
  assign(
    cb_names[[var]],
    cb_obj,
    envir = model_env
  )
}

formula_with_lag <- as.formula(
  paste(
    "beta ~",
    paste(cb_names, collapse = " + "),
    "+ ns(seq, 8) + offset(log(population)) + factor(name) + lag.value1"
  ),
  env = model_env
)

formula_without_lag <- as.formula(
  paste(
    "beta ~",
    paste(cb_names, collapse = " + "),
    "+ ns(seq, 8) + offset(log(population)) + factor(name)"
  ),
  env = model_env
)

full_model <- glm(
  formula_with_lag,
  family = quasipoisson(),
  data = data,
  na.action = na.exclude
)

model_without_lag <- glm(
  formula_without_lag,
  family = quasipoisson(),
  data = data,
  na.action = na.exclude
)

pseudo_r2 <- calc_pseudo_r2(full_model)
dispersion <- summary(full_model)$dispersion

cat("\n========== Model diagnostic results ==========\n")
cat("Deviance-based pseudo R2: ", round(pseudo_r2, 4), "\n", sep = "")
cat("Dispersion parameter: ", round(dispersion, 4), "\n", sep = "")

if (dispersion > 1.5) {
  cat("Overdispersion is present. The quasipoisson family accounts for this.\n")
} else {
  cat("No strong overdispersion was detected based on the dispersion parameter.\n")
}

capture.output(
  summary(full_model),
  file = file.path(output_dir, "full_model_summary.txt")
)

capture.output(
  summary(model_without_lag),
  file = file.path(output_dir, "model_without_lag_summary.txt")
)

data$resid_with_lag <- as.numeric(residuals(full_model, type = "deviance"))
data$fitted_with_lag <- as.numeric(fitted(full_model))

data$resid_without_lag <- as.numeric(residuals(model_without_lag, type = "deviance"))
data$fitted_without_lag <- as.numeric(fitted(model_without_lag))

diag_data <- data[!is.na(data$resid_with_lag), , drop = FALSE]

if (nrow(diag_data) < 5) {
  stop("Too few residuals are available for model diagnostics.")
}

resid_with_lag <- diag_data$resid_with_lag

acf_with_lag_obj <- acf(
  resid_with_lag,
  lag.max = 20,
  plot = FALSE,
  na.action = na.omit
)

pacf_with_lag_obj <- pacf(
  resid_with_lag,
  lag.max = 20,
  plot = FALSE,
  na.action = na.omit
)

acf_lag1 <- as.numeric(acf_with_lag_obj$acf[2])
pacf_lag1 <- as.numeric(pacf_with_lag_obj$acf[1])
acf_ci <- qnorm(0.975) / sqrt(length(resid_with_lag))

cat("\nResidual autocorrelation based on pooled ordered residuals:\n")
cat("Lag-1 ACF: ", round(acf_lag1, 4), "\n", sep = "")
cat("Lag-1 PACF: ", round(pacf_lag1, 4), "\n", sep = "")
cat("Approximate 95% CI: ±", round(acf_ci, 4), "\n", sep = "")

if (abs(acf_lag1) > acf_ci) {
  cat("Lag-1 ACF is outside the approximate 95% CI.\n")
} else {
  cat("Lag-1 ACF is within the approximate 95% CI.\n")
}

region_acf_with_lag <- calc_lag1_correlation_by_region(
  df = data,
  resid_var = "resid_with_lag"
)

region_acf_without_lag <- calc_lag1_correlation_by_region(
  df = data,
  resid_var = "resid_without_lag"
)

region_acf_comparison <- region_acf_with_lag %>%
  rename(
    lag1_correlation_with_lag = lag1_correlation,
    n_residuals_with_lag = n_residuals
  ) %>%
  left_join(
    region_acf_without_lag %>%
      rename(
        lag1_correlation_without_lag = lag1_correlation,
        n_residuals_without_lag = n_residuals
      ),
    by = "name"
  ) %>%
  mutate(
    abs_lag1_with_lag = abs(lag1_correlation_with_lag),
    abs_lag1_without_lag = abs(lag1_correlation_without_lag),
    reduced_abs_lag1 = abs_lag1_with_lag < abs_lag1_without_lag
  )

region_summary_with_lag <- summarise_region_acf(
  region_acf_comparison,
  "lag1_correlation_with_lag"
)

region_summary_without_lag <- summarise_region_acf(
  region_acf_comparison,
  "lag1_correlation_without_lag"
)

acf_reduction_rate <- mean(
  region_acf_comparison$reduced_abs_lag1,
  na.rm = TRUE
)

cat("\nRegion-specific lag-1 residual autocorrelation:\n")
cat(
  "Median absolute lag-1 correlation with autoregressive term: ",
  round(region_summary_with_lag$median_abs_lag1_correlation, 4),
  "\n",
  sep = ""
)
cat(
  "Median absolute lag-1 correlation without autoregressive term: ",
  round(region_summary_without_lag$median_abs_lag1_correlation, 4),
  "\n",
  sep = ""
)
cat(
  "Proportion of regions with reduced absolute lag-1 correlation: ",
  round(acf_reduction_rate, 4),
  "\n",
  sep = ""
)

if (
  is.finite(region_summary_with_lag$median_abs_lag1_correlation) &&
    is.finite(region_summary_without_lag$median_abs_lag1_correlation) &&
    region_summary_with_lag$median_abs_lag1_correlation <
    region_summary_without_lag$median_abs_lag1_correlation
) {
  cat("The autoregressive term reduced the median absolute lag-1 residual correlation across regions.\n")
} else {
  cat("The autoregressive term did not clearly reduce the median absolute lag-1 residual correlation across regions.\n")
}

diagnostic_summary <- data.frame(
  Metric = c(
    "N_used",
    "Deviance_based_pseudo_R2",
    "Dispersion_parameter",
    "Pooled_residual_lag1_ACF",
    "Pooled_residual_lag1_PACF",
    "Approximate_ACF_95CI",
    "Median_abs_region_lag1_ACF_with_AR",
    "Median_abs_region_lag1_ACF_without_AR",
    "Proportion_regions_with_reduced_abs_lag1_ACF"
  ),
  Value = c(
    nobs(full_model),
    pseudo_r2,
    dispersion,
    acf_lag1,
    pacf_lag1,
    acf_ci,
    region_summary_with_lag$median_abs_lag1_correlation,
    region_summary_without_lag$median_abs_lag1_correlation,
    acf_reduction_rate
  ),
  stringsAsFactors = FALSE
)

write.csv(
  diagnostic_summary,
  file.path(output_dir, "model_diagnostic_summary.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

write.csv(
  region_acf_comparison,
  file.path(output_dir, "region_lag1_residual_autocorrelation.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

plot_residual_diagnostics <- function(diag_data) {
  
  par(
    mfrow = c(2, 2),
    mar = c(5, 5, 4, 2)
  )
  
  plot(
    diag_data$fitted_with_lag,
    diag_data$resid_with_lag,
    xlab = "Fitted values",
    ylab = "Deviance residuals",
    main = "Residuals vs fitted values"
  )
  abline(h = 0, col = "red", lwd = 2)
  
  qqnorm(
    diag_data$resid_with_lag,
    main = "Q-Q plot of deviance residuals"
  )
  qqline(
    diag_data$resid_with_lag,
    col = "red",
    lwd = 2
  )
  
  plot(
    diag_data$resid_with_lag,
    type = "l",
    xlab = "Ordered observation index",
    ylab = "Deviance residuals",
    main = "Ordered residual series"
  )
  abline(h = 0, col = "red", lwd = 2)
  
  acf(
    diag_data$resid_with_lag,
    lag.max = 20,
    main = "ACF of deviance residuals",
    na.action = na.omit
  )
}

pdf(
  file = file.path(output_dir, "Residual_Diagnostics.pdf"),
  width = 10,
  height = 8
)

plot_residual_diagnostics(diag_data)

dev.off()

png(
  filename = file.path(output_dir, "Residual_Diagnostics.png"),
  width = 3000,
  height = 2400,
  res = 300
)

plot_residual_diagnostics(diag_data)

dev.off()

setEPS()

postscript(
  file = file.path(output_dir, "Residual_Diagnostics.eps"),
  width = 10,
  height = 8,
  horizontal = FALSE,
  onefile = FALSE,
  paper = "special"
)

plot_residual_diagnostics(diag_data)

dev.off()

cat("\nDiagnostic files saved to: ", output_dir, "\n", sep = "")
cat("PDF figure saved: ", file.path(output_dir, "Residual_Diagnostics.pdf"), "\n", sep = "")
cat("PNG figure saved: ", file.path(output_dir, "Residual_Diagnostics.png"), "\n", sep = "")
cat("EPS figure saved: ", file.path(output_dir, "Residual_Diagnostics.eps"), "\n", sep = "")
cat("Model diagnostics completed.\n")