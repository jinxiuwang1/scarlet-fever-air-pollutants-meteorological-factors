library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

data_file <- "E:/scarletfever_2013-2020_2.csv"
work_dir <- "E:/Scarlet Fever/Figure11"

if (!dir.exists(work_dir)) {
  dir.create(work_dir, recursive = TRUE)
}

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

setwd(work_dir)

write.csv(
  na_conversion_log,
  file.path(work_dir, "NA_conversion_log.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

write.csv(
  data,
  file.path(work_dir, "data_after_meteorological_imputation.csv"),
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

cb1.so2 <- create_crossbasis("SO2", data, lag_value = 3)
cb1.co <- create_crossbasis("CO", data, lag_value = 3)
cb1.no2 <- create_crossbasis("NO2", data, lag_value = 3)
cb1.o3 <- create_crossbasis("O3_8h", data, lag_value = 3)

cb1.sun <- create_crossbasis("sunlight_new", data, lag_value = 3)
cb1.humi <- create_crossbasis("humi_new", data, lag_value = 3)
cb1.rain <- create_crossbasis("rain_new", data, lag_value = 3)
cb1.temp <- create_crossbasis("meantemp_new", data, lag_value = 3)

full_model <- glm(
  beta ~
    cb1.so2 + cb1.co + cb1.no2 + cb1.o3 +
    cb1.sun + cb1.humi + cb1.rain + cb1.temp +
    ns(seq, 8 * 1) +
    offset(log(population)) +
    factor(name) +
    lag.value1,
  family = quasipoisson(),
  data = data,
  na.action = na.exclude
)

pseudo_r2_value <- calc_pseudo_r2(full_model)

cat("\nMultivariable model pseudo R2: ", round(pseudo_r2_value, 4), "\n", sep = "")

capture.output(
  summary(full_model),
  file = file.path(work_dir, "multivariable_model_summary.txt")
)

model_rows <- as.integer(rownames(model.frame(full_model)))

if (any(is.na(model_rows))) {
  model_rows <- seq_len(nrow(data))
}

model_data <- data[model_rows, , drop = FALSE]

cat("\nOriginal rows: ", nrow(data), "\n", sep = "")
cat("Rows used in the multivariable DLNM model: ", nrow(model_data), "\n", sep = "")

get_overall_p <- function(model, term_name) {
  
  out <- tryCatch(
    {
      drop1(model, test = "F")
    },
    error = function(e) NULL
  )
  
  if (is.null(out)) {
    return(NA_real_)
  }
  
  if (!(term_name %in% rownames(out))) {
    return(NA_real_)
  }
  
  return(as.numeric(out[term_name, "Pr(>F)"]))
}

overall_p_values <- data.frame(
  Variable = character(),
  Model_term = character(),
  Overall_P = numeric(),
  stringsAsFactors = FALSE
)

term_map <- list(
  SO2 = "cb1.so2",
  CO = "cb1.co",
  NO2 = "cb1.no2",
  O3 = "cb1.o3",
  Sunlight = "cb1.sun",
  Humidity = "cb1.humi",
  Precipitation = "cb1.rain",
  Temperature = "cb1.temp"
)

for (nm in names(term_map)) {
  
  overall_p_values <- rbind(
    overall_p_values,
    data.frame(
      Variable = nm,
      Model_term = term_map[[nm]],
      Overall_P = get_overall_p(full_model, term_map[[nm]]),
      stringsAsFactors = FALSE
    )
  )
}

write.csv(
  overall_p_values,
  file.path(work_dir, "overall_significance_results.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

plot_contour <- function(pred_obj, file_name, xlab_value) {
  
  tiff(
    file = file.path(work_dir, file_name),
    width = 3000,
    height = 2500,
    res = 300
  )
  
  on.exit(dev.off())
  
  plot(
    pred_obj,
    "contour",
    xlab = xlab_value,
    ylab = "Lag (months)",
    key.title = title("RR", cex.main = 1.2),
    plot.axes = {
      axis(1, cex.axis = 1.8)
      axis(2, cex.axis = 1.8)
    },
    key.axes = axis(4, cex.axis = 1.8),
    main = "Contour plot",
    cex.main = 2.0,
    cex.lab = 1.8
  )
}

run_crosspred <- function(basis_name, cen_value, at_values, bylag_value = 0.2) {
  
  expr <- substitute(
    crosspred(
      BASIS,
      full_model,
      cen = CEN,
      at = AT,
      bylag = BYLAG
    ),
    list(
      BASIS = as.name(basis_name),
      CEN = cen_value,
      AT = at_values,
      BYLAG = bylag_value
    )
  )
  
  eval(expr, envir = .GlobalEnv)
}

make_prediction_and_rr <- function(
    variable_name,
    basis_name,
    display_name,
    file_name,
    xlab_value,
    increment,
    increment_label,
    p_value
) {
  
  x <- model_data[[variable_name]]
  x_range <- range(x, na.rm = TRUE)
  cen_value <- median(x, na.rm = TRUE)
  target_value <- cen_value + increment
  
  at_values <- seq(
    x_range[1],
    x_range[2],
    length.out = 50
  )
  
  pred_obj <- run_crosspred(
    basis_name = basis_name,
    cen_value = cen_value,
    at_values = at_values,
    bylag_value = 0.2
  )
  
  plot_contour(
    pred_obj = pred_obj,
    file_name = file_name,
    xlab_value = xlab_value
  )
  
  # ----- Peak (global maximum from the contour surface) -----
  rr_mat <- pred_obj$matRRfit
  max_pos <- which(rr_mat == max(rr_mat, na.rm = TRUE), arr.ind = TRUE)[1, , drop = FALSE]
  peak_exposure <- pred_obj$predvar[max_pos[1, 1]]
  peak_lag <- pred_obj$lag[max_pos[1, 2]]
  peak_RR <- rr_mat[max_pos[1, 1], max_pos[1, 2]]
  peak_CI_low <- pred_obj$matRRlow[max_pos[1, 1], max_pos[1, 2]]
  peak_CI_high <- pred_obj$matRRhigh[max_pos[1, 1], max_pos[1, 2]]
  
  peak_info <- data.frame(
    Variable = display_name,
    Increment = increment_label,
    Peak_exposure = peak_exposure,
    Peak_lag_months = peak_lag,
    Peak_RR = peak_RR,
    CI_low = peak_CI_low,
    CI_high = peak_CI_high,
    Peak_RR_95CI = sprintf("%.3f (%.3f, %.3f)", peak_RR, peak_CI_low, peak_CI_high),
    stringsAsFactors = FALSE
  )
  
  # ----- Cumulative RR at fixed increment -----
  if (target_value > x_range[2]) {
    cum_row <- data.frame(
      Variable = display_name,
      Increment = increment_label,
      Overall_P = as.numeric(p_value),
      Center = cen_value,
      Target = target_value,
      RR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      RR_95CI = NA_character_,
      Note = "Target value is outside the observed range.",
      stringsAsFactors = FALSE
    )
    lag_rows <- data.frame()
  } else {
    pred_rr <- tryCatch(
      {
        run_crosspred(
          basis_name = basis_name,
          cen_value = cen_value,
          at_values = sort(unique(c(cen_value, target_value))),
          bylag_value = 0.2
        )
      },
      error = function(e) NULL
    )
    
    if (is.null(pred_rr)) {
      cum_row <- data.frame(
        Variable = display_name,
        Increment = increment_label,
        Overall_P = as.numeric(p_value),
        Center = cen_value,
        Target = target_value,
        RR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        RR_95CI = NA_character_,
        Note = "RR calculation failed.",
        stringsAsFactors = FALSE
      )
      lag_rows <- data.frame()
    } else {
      idx <- which.min(abs(pred_rr$predvar - target_value))
      rr_value <- as.numeric(pred_rr$allRRfit[idx])
      rr_low <- as.numeric(pred_rr$allRRlow[idx])
      rr_high <- as.numeric(pred_rr$allRRhigh[idx])
      rr_text <- sprintf("%.3f (%.3f, %.3f)", rr_value, rr_low, rr_high)
      
      cum_row <- data.frame(
        Variable = display_name,
        Increment = increment_label,
        Overall_P = as.numeric(p_value),
        Center = cen_value,
        Target = target_value,
        RR = rr_value,
        CI_low = rr_low,
        CI_high = rr_high,
        RR_95CI = rr_text,
        Note = "",
        stringsAsFactors = FALSE
      )
      
      # ----- Lag-specific RR (lag0-lag3) -----
      pred_lag <- tryCatch(
        {
          run_crosspred(
            basis_name = basis_name,
            cen_value = cen_value,
            at_values = target_value,
            bylag_value = 1
          )
        },
        error = function(e) NULL
      )
      
      if (!is.null(pred_lag)) {
        lag_rr <- pred_lag$matRRfit[1, ]
        lag_low <- pred_lag$matRRlow[1, ]
        lag_high <- pred_lag$matRRhigh[1, ]
        lag_text <- sprintf("%.3f (%.3f, %.3f)", lag_rr, lag_low, lag_high)
        
        lag_rows <- data.frame(
          Variable = display_name,
          Increment = increment_label,
          Lag = 0:3,
          RR = lag_rr,
          CI_low = lag_low,
          CI_high = lag_high,
          RR_95CI = lag_text,
          stringsAsFactors = FALSE
        )
      } else {
        lag_rows <- data.frame()
      }
    }
  }
  
  return(list(cumulative = cum_row, lag_specific = lag_rows, peak = peak_info))
}

analysis_list <- list(
  SO2 = list(
    variable = "SO2",
    basis_name = "cb1.so2",
    display = "SO2",
    file = "FIGURE_5A_SO2.tiff",
    xlab = expression(SO[2] ~ "(" * mu * g / m^3 * ")"),
    increment = 10,
    increment_label = "SO2 (+10 ug/m³)"
  ),
  CO = list(
    variable = "CO",
    basis_name = "cb1.co",
    display = "CO",
    file = "FIGURE_5B_CO.tiff",
    xlab = expression(CO ~ "(mg/m"^3*")"),
    increment = 1,
    increment_label = "CO (+1 mg/m³)"
  ),
  NO2 = list(
    variable = "NO2",
    basis_name = "cb1.no2",
    display = "NO2",
    file = "FIGURE_5C_NO2.tiff",
    xlab = expression(NO[2] ~ "(" * mu * g / m^3 * ")"),
    increment = 10,
    increment_label = "NO2 (+10 ug/m³)"
  ),
  O3 = list(
    variable = "O3_8h",
    basis_name = "cb1.o3",
    display = "O3",
    file = "FIGURE_5D_O3.tiff",
    xlab = expression(O[3] ~ "(" * mu * g / m^3 * ")"),
    increment = 10,
    increment_label = "O3 (+10 ug/m³)"
  ),
  Sunlight = list(
    variable = "sunlight_new",
    basis_name = "cb1.sun",
    display = "Sunlight",
    file = "FIGURE_5E_Sunlight.tiff",
    xlab = "Sunlight (hours)",
    increment = 5,
    increment_label = "Sunlight (+5 h)"
  ),
  Humidity = list(
    variable = "humi_new",
    basis_name = "cb1.humi",
    display = "Humidity",
    file = "FIGURE_5F_Humidity.tiff",
    xlab = "Relative humidity (%)",
    increment = 10,
    increment_label = "Relative humidity (+10%)"
  ),
  Precipitation = list(
    variable = "rain_new",
    basis_name = "cb1.rain",
    display = "Precipitation",
    file = "FIGURE_5G_Precipitation.tiff",
    xlab = "Precipitation (mm)",
    increment = 20,
    increment_label = "Precipitation (+20 mm)"
  ),
  Temperature = list(
    variable = "meantemp_new",
    basis_name = "cb1.temp",
    display = "Mean temperature",
    file = "FIGURE_5H_Temperature.tiff",
    xlab = "Mean temperature (°C)",
    increment = 5,
    increment_label = "Mean temperature (+5 °C)"
  )
)

cumulative_results <- data.frame()
lag_specific_results <- data.frame()
peak_results <- data.frame()

for (nm in names(analysis_list)) {
  
  info <- analysis_list[[nm]]
  
  p_value <- overall_p_values$Overall_P[
    overall_p_values$Variable == nm
  ]
  
  result <- make_prediction_and_rr(
    variable_name = info$variable,
    basis_name = info$basis_name,
    display_name = info$display,
    file_name = info$file,
    xlab_value = info$xlab,
    increment = info$increment,
    increment_label = info$increment_label,
    p_value = p_value
  )
  
  cumulative_results <- rbind(cumulative_results, result$cumulative)
  if (nrow(result$lag_specific) > 0) {
    lag_specific_results <- rbind(lag_specific_results, result$lag_specific)
  }
  peak_results <- rbind(peak_results, result$peak)
}

cat("\nAdjusted cumulative RR estimates with overall significance:\n")
print(cumulative_results)

write.csv(
  cumulative_results,
  file.path(work_dir, "adjusted_cumulative_RR_with_P_values.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

if (nrow(lag_specific_results) > 0) {
  write.csv(
    lag_specific_results,
    file.path(work_dir, "lag_specific_RR_results.csv"),
    row.names = FALSE,
    fileEncoding = "GB18030"
  )
  cat("\nLag-specific RR results (lag0-lag3) saved.\n")
}

if (nrow(peak_results) > 0) {
  write.csv(
    peak_results,
    file.path(work_dir, "peak_lag_summary.csv"),
    row.names = FALSE,
    fileEncoding = "GB18030"
  )
  cat("\nPeak lag summary saved.\n")
}

saveRDS(
  full_model,
  file.path(work_dir, "multivariable_dlnm_model.rds")
)

cat("\nAnalysis completed.\n")
cat("Contour plots saved to: ", work_dir, "\n", sep = "")
cat(
  "Cumulative RR table saved to: ",
  file.path(work_dir, "adjusted_cumulative_RR_with_P_values.csv"),
  "\n",
  sep = ""
)
cat(
  "Lag-specific RR table saved to: ",
  file.path(work_dir, "lag_specific_RR_results.csv"),
  "\n",
  sep = ""
)
cat(
  "Peak lag table saved to: ",
  file.path(work_dir, "peak_lag_summary.csv"),
  "\n",
  sep = ""
)