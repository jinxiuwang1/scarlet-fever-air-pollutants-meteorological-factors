library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

data_file <- "E:/scarletfever_2013-2020_2.csv"
output_dir <- "E:/Scarlet Fever/multivariable_sensitivity_analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
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

numeric_vars <- c(
  "year", "month",
  "beta", "population",
  pollution_vars,
  meteo_vars
)

for (var in numeric_vars) {
  
  if (!is.numeric(data[[var]])) {
    x_chr <- trimws(as.character(data[[var]]))
    x_chr[x_chr %in% missing_tokens] <- NA
    x_chr <- gsub(",", "", x_chr)
    data[[var]] <- suppressWarnings(as.numeric(x_chr))
  }
}

if (any(is.na(data$beta))) {
  stop("beta has missing values. beta should not be imputed.")
}

if (any(data$beta < 0, na.rm = TRUE)) {
  stop("beta has negative values. The current quasipoisson log-link model is not suitable.")
}

if (any(is.na(data$population))) {
  stop("population has missing values and cannot be used in population sensitivity analysis.")
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

data$log_population <- log(data$population)

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

make_arglag <- function(lag_value) {
  
  if (lag_value >= 3) {
    return(
      list(
        fun = "ns",
        knots = logknots(lag_value, nk = 2)
      )
    )
  }
  
  return(
    list(
      fun = "lin"
    )
  )
}

create_crossbasis <- function(var_name, df, lag_value) {
  
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

fit_sensitivity_model <- function(
    setting,
    max_lag,
    time_df,
    population_handling,
    include_ar,
    year_filter = "All"
) {
  
  data_used <- data
  
  if (!is.null(year_filter) && year_filter != "All") {
    if (year_filter == "<=2019") {
      data_used <- data[data$year <= 2019, , drop = FALSE]
    } else {
      stop("Unknown year_filter: ", year_filter)
    }
    rownames(data_used) <- NULL
  }
  
  center_map_used <- sapply(
    analysis_vars,
    function(v) {
      median(data_used[[v]], na.rm = TRUE)
    }
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
      lag_value = max_lag
    )
    
    assign(
      cb_names[[var]],
      cb_obj,
      envir = model_env
    )
  }
  
  model_terms <- c(
    cb_names,
    paste0("ns(seq, ", time_df, ")")
  )
  
  if (population_handling == "offset") {
    model_terms <- c(model_terms, "offset(log(population))")
  }
  
  if (population_handling == "covariate") {
    model_terms <- c(model_terms, "log_population")
  }
  
  model_terms <- c(
    model_terms,
    "factor(name)"
  )
  
  if (include_ar) {
    model_terms <- c(model_terms, "lag.value1")
  }
  
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
  
  model_rows <- as.integer(rownames(model.frame(fit)))
  
  if (any(is.na(model_rows))) {
    model_rows <- seq_len(nrow(data_used))
  }
  
  model_data <- data_used[model_rows, , drop = FALSE]
  
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
          c(
            RR = NA_real_,
            CI_low = NA_real_,
            CI_high = NA_real_
          )
        }
      )
    }
    
    model_result <- rbind(
      model_result,
      data.frame(
        Setting = setting,
        Year_filter = year_filter,
        Model_formula = formula_text,
        Max_lag = max_lag,
        Time_df = time_df,
        Population_handling = population_handling,
        AR_term = ifelse(include_ar, "With lag.value1", "Without lag.value1"),
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
          rr = as.numeric(rr_values["RR"]),
          low = as.numeric(rr_values["CI_low"]),
          high = as.numeric(rr_values["CI_high"])
        ),
        Note = note_text,
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(model_result)
}

config_table <- data.frame(
  Setting = c(
    "Main model",
    "Lag 2 months",
    "Lag 4 months",
    "Lag 5 months",
    "Lag 6 months",
    "Time df 6",
    "Time df 10",
    "Time df 12",
    "No population adjustment",
    "log(population) covariate",
    "Without lag.value1",
    "Exclude 2020"
  ),
  Max_lag = c(
    3,
    2, 4, 5, 6,
    3, 3, 3,
    3, 3,
    3,
    3
  ),
  Time_df = c(
    8,
    8, 8, 8, 8,
    6, 10, 12,
    8, 8,
    8,
    8
  ),
  Population_handling = c(
    "offset",
    "offset", "offset", "offset", "offset",
    "offset", "offset", "offset",
    "none", "covariate",
    "offset",
    "offset"
  ),
  Include_AR = c(
    TRUE,
    TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE,
    TRUE, TRUE,
    FALSE,
    TRUE
  ),
  Year_filter = c(
    "All",
    "All", "All", "All", "All",
    "All", "All", "All",
    "All", "All",
    "All",
    "<=2019"
  ),
  stringsAsFactors = FALSE
)

all_results <- data.frame()

error_log <- data.frame(
  Setting = character(),
  Error_message = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(config_table))) {
  
  cfg <- config_table[i, ]
  
  cat("Fitting: ", cfg$Setting, "\n", sep = "")
  
  result_i <- tryCatch(
    {
      fit_sensitivity_model(
        setting = cfg$Setting,
        max_lag = cfg$Max_lag,
        time_df = cfg$Time_df,
        population_handling = cfg$Population_handling,
        include_ar = cfg$Include_AR,
        year_filter = cfg$Year_filter
      )
    },
    error = function(e) {
      
      error_log <<- rbind(
        error_log,
        data.frame(
          Setting = cfg$Setting,
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

write.csv(
  all_results,
  file.path(output_dir, "sensitivity_results_long.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

rr_wide <- reshape(
  all_results[, c("Setting", "Year_filter", "Variable", "RR_95CI")],
  idvar = c("Setting", "Year_filter"),
  timevar = "Variable",
  direction = "wide"
)

names(rr_wide) <- gsub("^RR_95CI\\.", "", names(rr_wide))

write.csv(
  rr_wide,
  file.path(output_dir, "sensitivity_RR_95CI_wide.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

formula_table <- unique(
  all_results[, c(
    "Setting",
    "Year_filter",
    "Model_formula",
    "Max_lag",
    "Time_df",
    "Population_handling",
    "AR_term",
    "N_used",
    "Pseudo_R2",
    "Dispersion"
  )]
)

write.csv(
  formula_table,
  file.path(output_dir, "sensitivity_model_formulas.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

if (nrow(error_log) > 0) {
  write.csv(
    error_log,
    file.path(output_dir, "sensitivity_error_log.csv"),
    row.names = FALSE,
    fileEncoding = "GB18030"
  )
}

cat("\nSensitivity analysis completed.\n")
cat("Long results saved to: ",
    file.path(output_dir, "sensitivity_results_long.csv"),
    "\n",
    sep = "")
cat("Wide RR table saved to: ",
    file.path(output_dir, "sensitivity_RR_95CI_wide.csv"),
    "\n",
    sep = "")
cat("Model formulas saved to: ",
    file.path(output_dir, "sensitivity_model_formulas.csv"),
    "\n",
    sep = "")

if (nrow(error_log) > 0) {
  cat("Error log saved to: ",
      file.path(output_dir, "sensitivity_error_log.csv"),
      "\n",
      sep = "")
}

cat("All output files saved in: ", output_dir, "\n", sep = "").