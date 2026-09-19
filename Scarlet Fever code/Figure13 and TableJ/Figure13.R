library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

output_dir <- "E:/Scarlet Fever/Figure13"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

setwd(output_dir)

data <- read.csv(
  "E:/scarletfever_2013-2020_north.csv",
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

cat("\nMultivariable model pseudo R2: ", round(pseudo_r2_value, 4), "\n\n", sep = "")

capture.output(
  summary(full_model),
  file = file.path(output_dir, "multivariable_model_summary.txt")
)

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

run_crossreduce <- function(basis_name, cen_value, at_value = NULL) {
  
  if (is.null(at_value)) {
    expr <- substitute(
      crossreduce(
        BASIS,
        full_model,
        cen = CEN
      ),
      list(
        BASIS = as.name(basis_name),
        CEN = cen_value
      )
    )
  } else {
    expr <- substitute(
      crossreduce(
        BASIS,
        full_model,
        cen = CEN,
        at = AT
      ),
      list(
        BASIS = as.name(basis_name),
        CEN = cen_value,
        AT = at_value
      )
    )
  }
  
  eval(expr, envir = .GlobalEnv)
}

run_crosspred <- function(basis_name, cen_value, at_values) {
  
  expr <- substitute(
    crosspred(
      BASIS,
      full_model,
      cen = CEN,
      at = AT
    ),
    list(
      BASIS = as.name(basis_name),
      CEN = cen_value,
      AT = at_values
    )
  )
  
  eval(expr, envir = .GlobalEnv)
}

plot_cumulative_curve <- function(basis_name, cen_value, xlab_value, panel_title) {
  
  cr_obj <- run_crossreduce(
    basis_name = basis_name,
    cen_value = cen_value
  )
  
  plot(
    cr_obj,
    xlab = xlab_value,
    ylab = "RR",
    col = 2,
    lwd = 2,
    cex.lab = 1.8,
    cex.axis = 1.8,
    main = ""
  )
  
  mtext(
    panel_title,
    cex = 1.6,
    line = 1.4
  )
  
  return(cr_obj)
}

get_increment_rr <- function(basis_name, cen_value, target_value) {
  
  pred_obj <- tryCatch(
    {
      run_crosspred(
        basis_name = basis_name,
        cen_value = cen_value,
        at_values = sort(unique(c(cen_value, target_value)))
      )
    },
    error = function(e) NULL
  )
  
  if (is.null(pred_obj)) {
    return(
      c(
        RR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
    )
  }
  
  idx <- which.min(abs(pred_obj$predvar - target_value))
  
  return(
    c(
      RR = as.numeric(pred_obj$allRRfit[idx]),
      CI_low = as.numeric(pred_obj$allRRlow[idx]),
      CI_high = as.numeric(pred_obj$allRRhigh[idx])
    )
  )
}

var_info <- list(
  SO2 = list(
    var = "SO2",
    basis_name = "cb1.so2",
    term = "cb1.so2",
    xlab = expression(SO[2] ~ "(" * mu * g / m^3 * ")"),
    increment = 10,
    increment_label = "+10 ug/m3"
  ),
  CO = list(
    var = "CO",
    basis_name = "cb1.co",
    term = "cb1.co",
    xlab = expression(CO ~ "(mg/m"^3*")"),
    increment = 1,
    increment_label = "+1 mg/m³"
  ),
  NO2 = list(
    var = "NO2",
    basis_name = "cb1.no2",
    term = "cb1.no2",
    xlab = expression(NO[2] ~ "(" * mu * g / m^3 * ")"),
    increment = 10,
    increment_label = "+10 ug/m3"
  ),
  O3 = list(
    var = "O3_8h",
    basis_name = "cb1.o3",
    term = "cb1.o3",
    xlab = expression(O[3] ~ "(" * mu * g / m^3 * ")"),
    increment = 10,
    increment_label = "+10 ug/m3"
  ),
  Sunlight = list(
    var = "sunlight_new",
    basis_name = "cb1.sun",
    term = "cb1.sun",
    xlab = "Sunlight (hours)",
    increment = 5,
    increment_label = "+5 h"
  ),
  Humidity = list(
    var = "humi_new",
    basis_name = "cb1.humi",
    term = "cb1.humi",
    xlab = "Relative humidity (%)",
    increment = 10,
    increment_label = "+10%"
  ),
  Precipitation = list(
    var = "rain_new",
    basis_name = "cb1.rain",
    term = "cb1.rain",
    xlab = "Precipitation (mm)",
    increment = 20,
    increment_label = "+20 mm"
  ),
  `Mean temperature` = list(
    var = "meantemp_new",
    basis_name = "cb1.temp",
    term = "cb1.temp",
    xlab = "Mean temperature (°C)",
    increment = 5,
    increment_label = "+5 °C"
  )
)

draw_all_cumulative_curves <- function() {
  
  par(
    mfrow = c(2, 4),
    mar = c(6.6, 6.4, 6.2, 2.6),
    mgp = c(3.2, 1.1, 0),
    oma = c(0.6, 1.2, 1.0, 1.2)
  )
  
  curve_objects <- list()
  
  for (nm in names(var_info)) {
    
    info <- var_info[[nm]]
    cen_value <- median(data[[info$var]], na.rm = TRUE)
    
    curve_objects[[nm]] <- plot_cumulative_curve(
      basis_name = info$basis_name,
      cen_value = cen_value,
      xlab_value = info$xlab,
      panel_title = "Overall cumulative association"
    )
  }
  
  invisible(curve_objects)
}

setEPS()

postscript(
  file = file.path(output_dir, "Multivariate cumulative curve graph north.eps"),
  width = 16,
  height = 8,
  horizontal = FALSE
)

draw_all_cumulative_curves()

dev.off()

png(
  filename = file.path(output_dir, "Multivariate cumulative curve graph north.png"),
  width = 4800,
  height = 2400,
  res = 300
)

draw_all_cumulative_curves()

dev.off()

cat("EPS figure saved: ",
    file.path(output_dir, "Multivariate cumulative curve graph north.eps"),
    "\n",
    sep = "")

cat("PNG figure saved: ",
    file.path(output_dir, "Multivariate cumulative curve graph north.png"),
    "\n",
    sep = "")

results <- data.frame(
  Variable = character(),
  Increment = character(),
  Overall_P = numeric(),
  Center = numeric(),
  Target = numeric(),
  RR = numeric(),
  CI_low = numeric(),
  CI_high = numeric(),
  RR_95CI = character(),
  Note = character(),
  stringsAsFactors = FALSE
)

for (nm in names(var_info)) {
  
  info <- var_info[[nm]]
  
  x <- data[[info$var]]
  x_range <- range(x, na.rm = TRUE)
  cen_value <- median(x, na.rm = TRUE)
  target_value <- cen_value + info$increment
  
  p_value <- get_overall_p(
    model = full_model,
    term_name = info$term
  )
  
  if (target_value > x_range[2]) {
    
    results <- rbind(
      results,
      data.frame(
        Variable = nm,
        Increment = info$increment_label,
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
    )
    
  } else {
    
    rr_values <- get_increment_rr(
      basis_name = info$basis_name,
      cen_value = cen_value,
      target_value = target_value
    )
    
    rr_text <- ifelse(
      any(is.na(rr_values)),
      NA_character_,
      sprintf(
        "%.3f (%.3f, %.3f)",
        rr_values["RR"],
        rr_values["CI_low"],
        rr_values["CI_high"]
      )
    )
    
    results <- rbind(
      results,
      data.frame(
        Variable = nm,
        Increment = info$increment_label,
        Overall_P = as.numeric(p_value),
        Center = cen_value,
        Target = target_value,
        RR = as.numeric(rr_values["RR"]),
        CI_low = as.numeric(rr_values["CI_low"]),
        CI_high = as.numeric(rr_values["CI_high"]),
        RR_95CI = rr_text,
        Note = "",
        stringsAsFactors = FALSE
      )
    )
  }
}

cat("\nFinal RR, P-value, and 95% CI table:\n")
print(results)

write.csv(
  results,
  file.path(output_dir, "RR_CI_P_table.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

saveRDS(
  full_model,
  file.path(output_dir, "multivariable_cumulative_model.rds")
)

cat("\nAnalysis completed.\n")
cat("Table saved to: ",
    file.path(output_dir, "RR_CI_P_table.csv"),
    "\n",
    sep = "")
cat("All results saved in: ", output_dir, "\n", sep = "")