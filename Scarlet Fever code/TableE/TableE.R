library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

data_file <- "E:/scarletfever_2013-2020_2.csv"
output_dir <- "E:/Scarlet Fever/Lag_Selection"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data <- read.csv(data_file, fileEncoding = "GB18030")

is_empty_column <- function(x) {
  x_chr <- trimws(as.character(x))
  all(is.na(x_chr) | x_chr == "")
}

data <- data[, !sapply(data, is_empty_column)]

pollution_vars <- c("SO2", "CO", "NO2", "O3_8h")
meteo_vars <- c("rain", "sunlight", "humi", "meantemp")

required_vars <- c(
  "name", "year", "month", "beta", "population",
  pollution_vars, meteo_vars
)

data <- data[, required_vars]

missing_tokens <- c("NA", "N/A", "NaN", "missing", "", " ")

numeric_vars <- c(
  "year", "month", "beta", "population",
  pollution_vars, meteo_vars
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
  stop("beta has missing values.")
}

if (any(data$beta < 0, na.rm = TRUE)) {
  stop("beta has negative values.")
}

if (any(is.na(data$population))) {
  stop("population has missing values.")
}

if (any(data$population <= 0, na.rm = TRUE)) {
  stop("population has zero or negative values.")
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
          na_interpolation(
            x,
            option = "linear"
          )
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
  new_var <- paste0(var, "_imp")
  global_median <- median(data[[var]], na.rm = TRUE)

  data[[new_var]] <- ave(
    data[[var]],
    data$name,
    FUN = function(x) {
      safe_na_kalman(x, global_median)
    }
  )
}

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
  } else {
    return(
      list(
        fun = "lin"
      )
    )
  }
}

create_crossbasis <- function(var_name, df, lag_value) {
  x <- df[[var_name]]
  x_valid <- x[is.finite(x)]

  if (length(unique(x_valid)) < 2) {
    stop(
      var_name,
      " has insufficient variation."
    )
  }

  x_range <- range(
    x_valid,
    na.rm = TRUE
  )

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
    arglag = make_arglag(lag_value),
    group = df$name
  )

  return(cb)
}

pseudo_r2 <- function(model) {
  if (
    is.null(model$deviance) ||
    is.null(model$null.deviance)
  ) {
    return(NA_real_)
  }

  if (
    !is.finite(model$null.deviance) ||
    model$null.deviance == 0
  ) {
    return(NA_real_)
  }

  return(
    1 - model$deviance / model$null.deviance
  )
}

calc_qaic <- function(model) {
  D <- deviance(model)
  phi <- summary(model)$dispersion
  k <- length(coef(model))

  return(
    D / phi + 2 * phi * k
  )
}

calc_qbic <- function(model) {
  D <- deviance(model)
  phi <- summary(model)$dispersion
  k <- length(coef(model))
  n <- nobs(model)

  return(
    D / phi + log(n) * phi * k
  )
}

candidate_lags <- 1:6

dlnm_vars <- c(
  "SO2",
  "CO",
  "NO2",
  "O3_8h",
  "sunlight_imp",
  "humi_imp",
  "rain_imp",
  "meantemp_imp"
)

cb_by_lag <- list()

for (lag_val in candidate_lags) {
  lag_key <- as.character(lag_val)

  cb_by_lag[[lag_key]] <- list()

  for (var in dlnm_vars) {
    cb_by_lag[[lag_key]][[var]] <- create_crossbasis(
      var,
      data,
      lag_val
    )
  }
}

base_complete <- complete.cases(
  data[, c(
    "beta",
    "lag.beta",
    "population",
    "seq",
    "name"
  )]
)

cb_complete_all <- rep(
  TRUE,
  nrow(data)
)

for (lag_val in candidate_lags) {
  lag_key <- as.character(lag_val)

  for (var in dlnm_vars) {
    cb_complete_all <- cb_complete_all &
      complete.cases(
        cb_by_lag[[lag_key]][[var]]
      )
  }
}

common_rows <- base_complete & cb_complete_all

analysis_data <- data[
  common_rows,
  ,
  drop = FALSE
]

cat(
  "Number of complete observations used for model comparison:",
  nrow(analysis_data),
  "\n"
)

fit_full_model_by_lag <- function(lag_val) {
  lag_key <- as.character(lag_val)

  model_env <- new.env(
    parent = globalenv()
  )

  assign(
    "cb_so2",
    cb_by_lag[[lag_key]][["SO2"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  assign(
    "cb_co",
    cb_by_lag[[lag_key]][["CO"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  assign(
    "cb_no2",
    cb_by_lag[[lag_key]][["NO2"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  assign(
    "cb_o3",
    cb_by_lag[[lag_key]][["O3_8h"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  assign(
    "cb_sun",
    cb_by_lag[[lag_key]][["sunlight_imp"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  assign(
    "cb_humi",
    cb_by_lag[[lag_key]][["humi_imp"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  assign(
    "cb_rain",
    cb_by_lag[[lag_key]][["rain_imp"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  assign(
    "cb_temp",
    cb_by_lag[[lag_key]][["meantemp_imp"]][
      common_rows,
      ,
      drop = FALSE
    ],
    envir = model_env
  )

  full_formula <- as.formula(
    paste(
      "beta ~ cb_so2 + cb_co + cb_no2 + cb_o3 +",
      "cb_sun + cb_humi + cb_rain + cb_temp +",
      "ns(seq, 8) + offset(log(population)) + factor(name) + lag.beta"
    ),
    env = model_env
  )

  glm(
    full_formula,
    family = quasipoisson(),
    data = analysis_data,
    na.action = na.fail
  )
}

results_table <- data.frame()
fit_by_lag <- list()

for (lag_val in candidate_lags) {
  cat(
    "\nFitting model with maximum lag =",
    lag_val,
    "months...\n"
  )

  fit <- fit_full_model_by_lag(
    lag_val
  )

  fit_by_lag[[as.character(lag_val)]] <- fit

  r2 <- pseudo_r2(fit)
  qaic <- calc_qaic(fit)
  qbic <- calc_qbic(fit)

  results_table <- rbind(
    results_table,
    data.frame(
      Lag = lag_val,
      N = nobs(fit),
      PseudoR2 = round(r2, 4),
      QAIC = round(qaic, 2),
      QBIC = round(qbic, 2),
      Deviance = round(
        deviance(fit),
        1
      ),
      DF_resid = fit$df.residual,
      stringsAsFactors = FALSE
    )
  )

  cat(
    sprintf(
      "  Pseudo-R2 = %.4f, QAIC = %.2f, QBIC = %.2f\n",
      r2,
      qaic,
      qbic
    )
  )
}

cat(
  "\n========== Model Comparison Results ==========\n"
)

print(results_table)

write.csv(
  results_table,
  file.path(
    output_dir,
    "model_selection_QAIC_QBIC.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

saveRDS(
  fit_by_lag,
  file.path(
    output_dir,
    "fitted_models_by_lag.rds"
  )
)

cat(
  "\nAll results have been saved to:",
  output_dir,
  "\n"
)