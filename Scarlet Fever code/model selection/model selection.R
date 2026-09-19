library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

data_file <- "E:/scarletfever_2013-2020_2.csv"
output_dir <- "E:/Scarlet Fever code/model_selection"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data <- read.csv(
  data_file,
  fileEncoding = "GB18030"
)

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

data <- data %>%
  group_by(name) %>%
  arrange(year, month, .by_group = TRUE) %>%
  mutate(
    lag.beta = dplyr::lag(beta, n = 1, default = NA),
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
    FUN = function(x) safe_na_kalman(x, global_median)
  )
}

cat("\nMissing values after meteorological imputation:\n")
print(
  sapply(
    data[paste0(meteo_vars, "_imp")],
    function(x) sum(is.na(x))
  )
)

write.csv(
  na_conversion_log,
  file.path(output_dir, "NA_conversion_log.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

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

create_crossbasis <- function(var_name, data, lag_value, lag_knots) {
  
  x <- data[[var_name]]
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
      knots = lag_knots
    ),
    group = data$name
  )
  
  return(cb)
}

lag_value <- 3
lag_knots <- logknots(lag_value, nk = 2)

cb.so2 <- create_crossbasis("SO2", data, lag_value, lag_knots)
cb.co <- create_crossbasis("CO", data, lag_value, lag_knots)
cb.no2 <- create_crossbasis("NO2", data, lag_value, lag_knots)
cb.o3 <- create_crossbasis("O3_8h", data, lag_value, lag_knots)
cb.sun <- create_crossbasis("sunlight_imp", data, lag_value, lag_knots)
cb.humi <- create_crossbasis("humi_imp", data, lag_value, lag_knots)
cb.rain <- create_crossbasis("rain_imp", data, lag_value, lag_knots)
cb.temp <- create_crossbasis("meantemp_imp", data, lag_value, lag_knots)

base_vars <- c(
  "beta", "lag.beta", "population",
  "seq", "name"
)

common_rows <- complete.cases(data[, base_vars]) &
  complete.cases(cb.so2) &
  complete.cases(cb.co) &
  complete.cases(cb.no2) &
  complete.cases(cb.o3) &
  complete.cases(cb.sun) &
  complete.cases(cb.humi) &
  complete.cases(cb.rain) &
  complete.cases(cb.temp)

analysis_data <- data[common_rows, ]

cb.so2 <- cb.so2[common_rows, , drop = FALSE]
cb.co <- cb.co[common_rows, , drop = FALSE]
cb.no2 <- cb.no2[common_rows, , drop = FALSE]
cb.o3 <- cb.o3[common_rows, , drop = FALSE]
cb.sun <- cb.sun[common_rows, , drop = FALSE]
cb.humi <- cb.humi[common_rows, , drop = FALSE]
cb.rain <- cb.rain[common_rows, , drop = FALSE]
cb.temp <- cb.temp[common_rows, , drop = FALSE]

cat("\nOriginal rows: ", nrow(data), "\n", sep = "")
cat("Analysis rows after applying a common complete-case sample: ", nrow(analysis_data), "\n", sep = "")

if (nrow(analysis_data) == 0) {
  stop("No rows remain after complete-case filtering.")
}

base_formula <- beta ~ ns(seq, 8 * 1) +
  offset(log(population)) +
  factor(name) +
  lag.beta

model_list <- list(
  "Base" = base_formula,
  "+ SO2" = update(base_formula, ~ . + cb.so2),
  "+ CO" = update(base_formula, ~ . + cb.so2 + cb.co),
  "+ NO2" = update(base_formula, ~ . + cb.so2 + cb.co + cb.no2),
  "+ O3" = update(base_formula, ~ . + cb.so2 + cb.co + cb.no2 + cb.o3),
  "+ Sunlight" = update(base_formula, ~ . + cb.so2 + cb.co + cb.no2 + cb.o3 + cb.sun),
  "+ Humidity" = update(base_formula, ~ . + cb.so2 + cb.co + cb.no2 + cb.o3 + cb.sun + cb.humi),
  "+ Precip" = update(base_formula, ~ . + cb.so2 + cb.co + cb.no2 + cb.o3 + cb.sun + cb.humi + cb.rain),
  "+ Temperature" = update(base_formula, ~ . + cb.so2 + cb.co + cb.no2 + cb.o3 + cb.sun + cb.humi + cb.rain + cb.temp)
)

pseudo_r2 <- function(model) {
  if (is.null(model$deviance) || is.null(model$null.deviance)) {
    return(NA_real_)
  }
  
  if (!is.finite(model$null.deviance) || model$null.deviance == 0) {
    return(NA_real_)
  }
  
  return(1 - model$deviance / model$null.deviance)
}

results <- data.frame(
  Model = character(),
  N = integer(),
  PseudoR2 = numeric(),
  Deviance = numeric(),
  DF_resid = numeric(),
  Dev_change = numeric(),
  DF_change = numeric(),
  F_value = numeric(),
  p_value_numeric = numeric(),
  p_value = character(),
  stringsAsFactors = FALSE
)

fit_list <- list()
prev_fit <- NULL

for (i in seq_along(model_list)) {
  
  model_name <- names(model_list)[i]
  form <- model_list[[i]]
  
  fit <- glm(
    form,
    family = quasipoisson(),
    data = analysis_data,
    na.action = na.fail
  )
  
  fit_list[[model_name]] <- fit
  
  r2 <- pseudo_r2(fit)
  dev <- fit$deviance
  df_resid <- fit$df.residual
  n_model <- nobs(fit)
  
  if (is.null(prev_fit)) {
    
    dev_change <- NA_real_
    df_change <- NA_real_
    f_val <- NA_real_
    p_val_numeric <- NA_real_
    p_val_text <- NA_character_
    
  } else {
    
    dev_change <- prev_fit$deviance - fit$deviance
    df_change <- prev_fit$df.residual - fit$df.residual
    
    if (is.finite(dev_change) && is.finite(df_change) && df_change > 0) {
      dispersion <- dev / df_resid
      f_val <- (dev_change / df_change) / dispersion
      p_val_numeric <- pf(f_val, df_change, df_resid, lower.tail = FALSE)
      p_val_text <- formatC(p_val_numeric, format = "e", digits = 2)
    } else {
      f_val <- NA_real_
      p_val_numeric <- NA_real_
      p_val_text <- NA_character_
    }
  }
  
  results <- rbind(
    results,
    data.frame(
      Model = model_name,
      N = n_model,
      PseudoR2 = round(r2, 4),
      Deviance = round(dev, 1),
      DF_resid = df_resid,
      Dev_change = ifelse(is.na(dev_change), NA, round(dev_change, 1)),
      DF_change = df_change,
      F_value = ifelse(is.na(f_val), NA, round(f_val, 3)),
      p_value_numeric = p_val_numeric,
      p_value = p_val_text,
      stringsAsFactors = FALSE
    )
  )
  
  cat("\nModel: ", model_name, "\n", sep = "")
  cat("  N = ", n_model, "\n", sep = "")
  cat("  Pseudo R2 = ", round(r2, 4), "\n", sep = "")
  cat("  Deviance = ", round(dev, 1), " (df = ", df_resid, ")\n", sep = "")
  
  if (!is.null(prev_fit)) {
    cat(
      "  F test vs previous model: F = ",
      round(f_val, 3),
      " (df1 = ",
      df_change,
      ", df2 = ",
      df_resid,
      "), p = ",
      p_val_text,
      "\n",
      sep = ""
    )
  }
  
  prev_fit <- fit
}

cat("\nModel selection results:\n")
print(results)

sig_rows <- which(!is.na(results$p_value_numeric) & results$p_value_numeric < 0.05)

if (length(sig_rows) > 0) {
  best_model_index <- sig_rows[length(sig_rows)]
} else {
  best_model_index <- 1
}

best_model_name <- results$Model[best_model_index]
best_model <- fit_list[[best_model_name]]

cat("\nSelected model:\n")
cat(best_model_name, "\n")

write.csv(
  results,
  file.path(output_dir, "model_selection_results.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

saveRDS(
  best_model,
  file.path(output_dir, "best_model.rds")
)

cat("\nResults saved to:\n")
cat(file.path(output_dir, "model_selection_results.csv"), "\n")
cat(file.path(output_dir, "best_model.rds"), "\n")