# ============================================================
# DLNM sensitivity analysis for different beta smoothing windows
# Outcome variants:
#   beta_cases      = no smoothing
#   beta_smooth_2   = 2-month moving average
#   beta_smooth_3   = 3-month moving average (MAIN MODEL)
#   beta_smooth_4   = 4-month moving average
#   beta_smooth_5   = 5-month moving average
#
# All other DLNM specifications are FIXED at the main-model values:
#   maximum lag = 3 months
#   time spline df = 8
#   population adjustment = offset(log(population))
#   province fixed effects = factor(name)
#   first-order AR term = lag.value1, recalculated for each beta series
# ============================================================

library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

# ------------------------
# 1. File paths
# ------------------------
# Change these two paths to your own Windows paths.
data_file <- "E:/scarletfever_2013-2020_3.csv"
output_dir <- "E:/Scarlet Fever/smoothing_window_sensitivity_analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------
# 2. Read and clean data
# ------------------------
data <- read.csv(
  data_file,
  fileEncoding = "GB18030",
  check.names = FALSE
)

is_empty_column <- function(x) {
  x_chr <- trimws(as.character(x))
  all(is.na(x_chr) | x_chr == "")
}

data <- data[, !sapply(data, is_empty_column), drop = FALSE]

pollution_vars <- c("SO2", "CO", "NO2", "O3_8h")
meteo_vars <- c("rain", "sunlight", "humi", "meantemp")

beta_vars <- c(
  "beta_cases",
  "beta_smooth_2",
  "beta_smooth_3",
  "beta_smooth_4",
  "beta_smooth_5"
)

beta_labels <- c(
  beta_cases    = "No smoothing",
  beta_smooth_2 = "2-month moving average",
  beta_smooth_3 = "3-month moving average (main model)",
  beta_smooth_4 = "4-month moving average",
  beta_smooth_5 = "5-month moving average"
)

required_vars <- c(
  "name", "province", "year", "month", "population",
  pollution_vars, meteo_vars, beta_vars
)

missing_vars <- setdiff(required_vars, names(data))
if (length(missing_vars) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_vars, collapse = ", ")
  )
}

data <- data[, required_vars, drop = FALSE]

missing_tokens <- c("NA", "N/A", "NaN", "missing", "", " ", "-")

numeric_vars <- c(
  "year", "month", "population",
  pollution_vars, meteo_vars, beta_vars
)

for (var in numeric_vars) {
  if (!is.numeric(data[[var]])) {
    x_chr <- trimws(as.character(data[[var]]))
    x_chr[x_chr %in% missing_tokens] <- NA
    x_chr <- gsub(",", "", x_chr)
    data[[var]] <- suppressWarnings(as.numeric(x_chr))
  }
}

# Basic checks
if (any(is.na(data$population))) {
  stop("population contains missing values.")
}
if (any(data$population <= 0, na.rm = TRUE)) {
  stop("population contains zero or negative values.")
}

for (b in beta_vars) {
  if (any(is.na(data[[b]]))) {
    stop(b, " contains missing values. beta should not be imputed.")
  }
  if (any(data[[b]] < 0, na.rm = TRUE)) {
    stop(b, " contains negative values; current quasipoisson log-link model is unsuitable.")
  }
}

# Check panel structure before modeling
panel_count <- data %>%
  count(name, province, name = "n_months")

if (nrow(panel_count) != 31) {
  warning("Expected 31 provinces, but found ", nrow(panel_count), ".")
}
if (any(panel_count$n_months != 96)) {
  warning("At least one province does not contain exactly 96 monthly observations.")
}

if (anyDuplicated(data[, c("name", "year", "month")]) > 0) {
  stop("Duplicate name-year-month records were found.")
}

# Sort by province and time
# The same order is used for every beta smoothing specification.
data <- data[with(data, order(name, year, month)), , drop = FALSE]
rownames(data) <- NULL

# ------------------------
# 3. Meteorological imputation
#    Same procedure as the original sensitivity-analysis code
# ------------------------
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
          x
        }
      )
    }
  )

  if (any(is.na(result))) {
    result[is.na(result)] <- global_median
  }

  result
}

# Save pre-imputation missing counts for checking
missing_before_imputation <- data.frame(
  Variable = meteo_vars,
  Missing_N = sapply(meteo_vars, function(v) sum(is.na(data[[v]]))),
  stringsAsFactors = FALSE
)

write.csv(
  missing_before_imputation,
  file.path(output_dir, "missing_values_before_imputation.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

for (var in meteo_vars) {
  new_var <- paste0(var, "_new")
  global_median <- median(data[[var]], na.rm = TRUE)

  if (!is.finite(global_median)) {
    stop(var, " is completely missing and cannot be imputed.")
  }

  data[[new_var]] <- ave(
    data[[var]],
    data$name,
    FUN = function(x) safe_na_kalman(x, global_median)
  )
}

# ------------------------
# 4. Exposure variables and increments
# ------------------------
analysis_vars <- c(
  "SO2",
  "CO",
  "NO2",
  "O3_8h",
  "sunlight_new",
  "humi_new",
  "rain_new",
  "meantemp_new"
)

display_names <- c(
  SO2 = "SO2",
  CO = "CO",
  NO2 = "NO2",
  O3_8h = "O3",
  sunlight_new = "Sunlight",
  humi_new = "Humidity",
  rain_new = "Precipitation",
  meantemp_new = "Temperature"
)

increment_map <- c(
  SO2 = 10,
  CO = 1,
  NO2 = 10,
  O3_8h = 10,
  sunlight_new = 5,
  humi_new = 10,
  rain_new = 20,
  meantemp_new = 5
)

increment_label_map <- c(
  SO2 = "+10 ug/m3",
  CO = "+1 mg/m3",
  NO2 = "+10 ug/m3",
  O3_8h = "+10 ug/m3",
  sunlight_new = "+5 h",
  humi_new = "+10%",
  rain_new = "+20 mm",
  meantemp_new = "+5 C"
)

# ------------------------
# 5. DLNM helper functions
# ------------------------
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
  q
}

# Reviewer-7 smoothing sensitivity keeps max lag FIXED at 3 months.
MAX_LAG <- 3
TIME_DF <- 8

make_arglag <- function(lag_value) {
  if (lag_value >= 3) {
    return(
      list(
        fun = "ns",
        knots = logknots(lag_value, nk = 2)
      )
    )
  }

  list(fun = "lin")
}

create_crossbasis <- function(var_name, df, lag_value = MAX_LAG) {
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
    argvar_list <- list(fun = "lin")
  }

  crossbasis(
    x,
    lag = lag_value,
    argvar = argvar_list,
    arglag = make_arglag(lag_value),
    group = df$name
  )
}

calc_pseudo_r2 <- function(model) {
  if (is.null(model$deviance) || is.null(model$null.deviance)) {
    return(NA_real_)
  }

  if (!is.finite(model$null.deviance) || model$null.deviance == 0) {
    return(NA_real_)
  }

  1 - model$deviance / model$null.deviance
}

get_overall_p <- function(model, term_name) {
  out <- tryCatch(
    drop1(model, test = "F"),
    error = function(e) NULL
  )

  if (is.null(out) || !(term_name %in% rownames(out))) {
    return(NA_real_)
  }

  as.numeric(out[term_name, "Pr(>F)"])
}

format_rr <- function(rr, low, high) {
  if (any(is.na(c(rr, low, high)))) {
    return(NA_character_)
  }
  sprintf("%.3f (%.3f, %.3f)", rr, low, high)
}

format_p <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }
  if (p < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p)
}

# ------------------------
# 6. Fit one model for one beta smoothing specification
# ------------------------
fit_one_smoothing_model <- function(beta_var, smoothing_label) {

  data_used <- data

  # IMPORTANT:
  # For each smoothing window, beta and lag.value1 must be rebuilt together.
  data_used$beta <- data_used[[beta_var]]

  data_used <- data_used %>%
    group_by(name) %>%
    arrange(year, month, .by_group = TRUE) %>%
    mutate(
      lag.value1 = dplyr::lag(beta, n = 1, default = NA),
      seq = row_number()
    ) %>%
    ungroup()

  # Exposure medians are identical across smoothing settings because only beta changes.
  center_map_used <- sapply(
    analysis_vars,
    function(v) median(data_used[[v]], na.rm = TRUE)
  )

  model_env <- new.env(parent = globalenv())

  cb_names <- setNames(
    paste0("cb_", make.names(analysis_vars)),
    analysis_vars
  )

  for (var in analysis_vars) {
    cb_obj <- create_crossbasis(
      var_name = var,
      df = data_used,
      lag_value = MAX_LAG
    )

    assign(
      cb_names[[var]],
      cb_obj,
      envir = model_env
    )
  }

  model_terms <- c(
    cb_names,
    paste0("ns(seq, ", TIME_DF, ")"),
    "offset(log(population))",
    "factor(name)",
    "lag.value1"
  )

  formula_text <- paste(
    "beta ~",
    paste(model_terms, collapse = " + ")
  )

  fit <- glm(
    as.formula(formula_text, env = model_env),
    family = quasipoisson(),
    data = data_used,
    na.action = na.exclude
  )

  assign("fit", fit, envir = model_env)

  # Actual rows used by the fitted model
  mf <- model.frame(fit)
  model_rows <- suppressWarnings(as.integer(rownames(mf)))

  if (length(model_rows) != nrow(mf) || any(is.na(model_rows))) {
    # Fallback; usually not needed with the current data structure
    model_data <- data_used[complete.cases(mf), , drop = FALSE]
  } else {
    model_data <- data_used[model_rows, , drop = FALSE]
  }

  pseudo_r2 <- calc_pseudo_r2(fit)
  dispersion <- summary(fit)$dispersion

  run_crosspred <- function(basis_name, cen_value, at_values) {
    expr <- substitute(
      crosspred(
        BASIS,
        fit,
        cen = CEN,
        at = AT
      ),
      list(
        BASIS = as.name(basis_name),
        CEN = cen_value,
        AT = at_values
      )
    )

    eval(expr, envir = model_env)
  }

  model_result <- data.frame()

  for (var in analysis_vars) {
    cb_name <- cb_names[[var]]
    x <- model_data[[var]]
    x_range <- range(x, na.rm = TRUE)

    cen_value <- as.numeric(center_map_used[[var]])
    increment <- as.numeric(increment_map[[var]])
    target_value <- cen_value + increment

    overall_p <- get_overall_p(
      model = fit,
      term_name = cb_name
    )

    rr_values <- c(
      RR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_
    )

    note_text <- ""

    if (any(!is.finite(x_range)) || diff(x_range) == 0) {
      note_text <- "Invalid exposure range."
    } else if (target_value > x_range[2]) {
      note_text <- "Target outside observed range."
    } else {
      rr_values <- tryCatch(
        {
          pred_obj <- run_crosspred(
            basis_name = cb_name,
            cen_value = cen_value,
            at_values = sort(unique(c(cen_value, target_value)))
          )

          idx <- which.min(abs(pred_obj$predvar - target_value))

          c(
            RR = as.numeric(pred_obj$allRRfit[idx]),
            CI_low = as.numeric(pred_obj$allRRlow[idx]),
            CI_high = as.numeric(pred_obj$allRRhigh[idx])
          )
        },
        error = function(e) {
          note_text <<- conditionMessage(e)
          c(RR = NA_real_, CI_low = NA_real_, CI_high = NA_real_)
        }
      )
    }

    model_result <- rbind(
      model_result,
      data.frame(
        Smoothing = smoothing_label,
        Beta_variable = beta_var,
        Is_main_model = ifelse(beta_var == "beta_smooth_3", "Yes", "No"),
        Model_formula = formula_text,
        Max_lag = MAX_LAG,
        Time_df = TIME_DF,
        Population_handling = "offset(log(population))",
        AR_term = "With lag.value1 (recalculated for this beta)",
        N_used = nobs(fit),
        Pseudo_R2 = round(pseudo_r2, 4),
        Dispersion = round(dispersion, 4),
        Variable = unname(display_names[[var]]),
        Exposure = var,
        Increment = unname(increment_label_map[[var]]),
        Overall_P = as.numeric(overall_p),
        P_value = format_p(overall_p),
        RR = as.numeric(rr_values["RR"]),
        CI_low = as.numeric(rr_values["CI_low"]),
        CI_high = as.numeric(rr_values["CI_high"]),
        RR_95CI = format_rr(
          as.numeric(rr_values["RR"]),
          as.numeric(rr_values["CI_low"]),
          as.numeric(rr_values["CI_high"])
        ),
        Note = note_text,
        stringsAsFactors = FALSE
      )
    )
  }

  model_result
}

# ------------------------
# 7. Run the five smoothing specifications
# ------------------------
all_results <- data.frame()
error_log <- data.frame(
  Smoothing = character(),
  Beta_variable = character(),
  Error_message = character(),
  stringsAsFactors = FALSE
)

for (b in beta_vars) {
  label <- unname(beta_labels[[b]])

  cat("Fitting: ", label, " [", b, "]\n", sep = "")

  result_i <- tryCatch(
    fit_one_smoothing_model(
      beta_var = b,
      smoothing_label = label
    ),
    error = function(e) {
      error_log <<- rbind(
        error_log,
        data.frame(
          Smoothing = label,
          Beta_variable = b,
          Error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      )
      NULL
    }
  )

  if (!is.null(result_i)) {
    all_results <- rbind(all_results, result_i)
  }
}

if (nrow(all_results) == 0) {
  stop("All smoothing-window models failed. Check the error log/data.")
}

# Preserve desired row order
all_results$Smoothing <- factor(
  all_results$Smoothing,
  levels = unname(beta_labels)
)
all_results <- all_results[order(all_results$Smoothing, all_results$Variable), ]
all_results$Smoothing <- as.character(all_results$Smoothing)
rownames(all_results) <- NULL

# ------------------------
# 8. Save full long-format results
# ------------------------
write.csv(
  all_results,
  file.path(output_dir, "smoothing_sensitivity_results_long.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# ------------------------
# 9. Create reviewer-friendly RR tables
# ------------------------
rr_wide <- reshape(
  all_results[, c("Smoothing", "Variable", "RR_95CI")],
  idvar = "Smoothing",
  timevar = "Variable",
  direction = "wide"
)

names(rr_wide) <- gsub("^RR_95CI\\.", "", names(rr_wide))

# Reorder settings after reshape
rr_wide$Smoothing <- factor(rr_wide$Smoothing, levels = unname(beta_labels))
rr_wide <- rr_wide[order(rr_wide$Smoothing), ]
rr_wide$Smoothing <- as.character(rr_wide$Smoothing)

write.csv(
  rr_wide,
  file.path(output_dir, "smoothing_sensitivity_RR_95CI_wide.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# Pollutants only (similar to Supplementary Table H)
pollutant_names <- c("SO2", "CO", "NO2", "O3")
pollutant_table <- rr_wide[, c("Smoothing", pollutant_names), drop = FALSE]

write.csv(
  pollutant_table,
  file.path(output_dir, "smoothing_sensitivity_pollutants_Table.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# Meteorological variables only (similar to Supplementary Table I)
weather_names <- c("Sunlight", "Humidity", "Precipitation", "Temperature")
weather_table <- rr_wide[, c("Smoothing", weather_names), drop = FALSE]

write.csv(
  weather_table,
  file.path(output_dir, "smoothing_sensitivity_weather_Table.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# Overall P-value table
p_wide <- reshape(
  all_results[, c("Smoothing", "Variable", "P_value")],
  idvar = "Smoothing",
  timevar = "Variable",
  direction = "wide"
)

names(p_wide) <- gsub("^P_value\\.", "", names(p_wide))
p_wide$Smoothing <- factor(p_wide$Smoothing, levels = unname(beta_labels))
p_wide <- p_wide[order(p_wide$Smoothing), ]
p_wide$Smoothing <- as.character(p_wide$Smoothing)

write.csv(
  p_wide,
  file.path(output_dir, "smoothing_sensitivity_overall_P_wide.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# Model-fit summary (one row per smoothing window)
fit_summary <- unique(
  all_results[, c(
    "Smoothing",
    "Beta_variable",
    "Is_main_model",
    "Max_lag",
    "Time_df",
    "Population_handling",
    "AR_term",
    "N_used",
    "Pseudo_R2",
    "Dispersion"
  )]
)

fit_summary$Smoothing <- factor(fit_summary$Smoothing, levels = unname(beta_labels))
fit_summary <- fit_summary[order(fit_summary$Smoothing), ]
fit_summary$Smoothing <- as.character(fit_summary$Smoothing)

write.csv(
  fit_summary,
  file.path(output_dir, "smoothing_sensitivity_model_fit.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

if (nrow(error_log) > 0) {
  write.csv(
    error_log,
    file.path(output_dir, "smoothing_sensitivity_error_log.csv"),
    row.names = FALSE,
    fileEncoding = "GB18030"
  )
}

# ------------------------
# 10. Console output
# ------------------------
cat("\n============================================================\n")
cat("Smoothing-window sensitivity analysis completed.\n")
cat("IMPORTANT: 3-month moving average is the main model.\n")
cat("All other DLNM specifications were held fixed.\n")
cat("lag.value1 was recalculated separately for every beta series.\n")
cat("============================================================\n\n")

print(fit_summary)
cat("\nRR (95% CI) table:\n")
print(rr_wide)

cat("\nOutput directory: ", output_dir, "\n", sep = "")

if (nrow(error_log) > 0) {
  cat("WARNING: Some models failed. Check smoothing_sensitivity_error_log.csv\n")
}
