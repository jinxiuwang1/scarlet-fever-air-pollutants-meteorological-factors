library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)
library(sandwich)

data_file <- "E:/scarletfever_2013-2020_2.csv"
output_dir <- "E:/Scarlet Fever/Province-clustered"

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
  "name", "province", "year", "month",
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

# -------------------------------------------------------------------------
# Province identifiers
# -------------------------------------------------------------------------
data$name <- trimws(as.character(data$name))
data$province <- trimws(as.character(data$province))

if (any(is.na(data$name) | data$name == "")) {
  stop("name contains missing/blank values.")
}

if (any(is.na(data$province) | data$province == "")) {
  stop("province contains missing/blank values.")
}

# Check one-to-one mapping between Chinese and English province names
name_to_province <- aggregate(
  province ~ name,
  data = data,
  FUN = function(x) length(unique(x))
)

if (any(name_to_province$province != 1)) {
  bad_names <- name_to_province$name[name_to_province$province != 1]
  stop(
    "Some Chinese province names map to multiple English province names: ",
    paste(bad_names, collapse = ", ")
  )
}

province_to_name <- aggregate(
  name ~ province,
  data = data,
  FUN = function(x) length(unique(x))
)

if (any(province_to_name$name != 1)) {
  bad_provinces <- province_to_name$province[province_to_name$name != 1]
  stop(
    "Some English province names map to multiple Chinese province names: ",
    paste(bad_provinces, collapse = ", ")
  )
}

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

# -------------------------------------------------------------------------
# Province-clustered Overall cross-basis P
# Joint Wald-F test of all coefficients in one cross-basis.
# The effective covariance rank is used as numerator df, which is safer when
# the clustered covariance matrix is close to singular.
# -------------------------------------------------------------------------
get_cluster_overall_p <- function(
    model,
    term_name,
    vcov_cluster,
    n_clusters
) {
  
  coef_names <- names(coef(model))
  idx <- which(startsWith(coef_names, term_name))
  
  if (length(idx) == 0 || n_clusters <= 1) {
    return(NA_real_)
  }
  
  b <- coef(model)[idx]
  V <- vcov_cluster[idx, idx, drop = FALSE]
  
  if (any(!is.finite(b)) || any(!is.finite(V))) {
    return(NA_real_)
  }
  
  q_eff <- qr(V, tol = 1e-10)$rank
  
  if (!is.finite(q_eff) || q_eff < 1) {
    return(NA_real_)
  }
  
  V_inv <- tryCatch(
    MASS::ginv(V),
    error = function(e) NULL
  )
  
  if (is.null(V_inv)) {
    return(NA_real_)
  }
  
  W <- as.numeric(t(b) %*% V_inv %*% b)
  
  if (!is.finite(W)) {
    return(NA_real_)
  }
  
  F_stat <- W / q_eff
  
  pf(
    F_stat,
    df1 = q_eff,
    df2 = n_clusters - 1,
    lower.tail = FALSE
  )
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

# -------------------------------------------------------------------------
# P value for the SPECIFIC cumulative RR contrast reported by crosspred().
# H0: log(RR) = 0, i.e. RR = 1.
# Uses crosspred()$allfit and $allse and therefore corresponds to the reported
# crosspred 95% CI, apart from display rounding.
# -------------------------------------------------------------------------
calc_rr_contrast_p <- function(log_rr, se_log_rr) {
  
  if (any(is.na(c(log_rr, se_log_rr)))) {
    return(NA_real_)
  }
  
  if (!is.finite(log_rr) || !is.finite(se_log_rr) || se_log_rr <= 0) {
    return(NA_real_)
  }
  
  z_value <- log_rr / se_log_rr
  2 * pnorm(-abs(z_value))
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
  
  # -----------------------------------------------------------------------
  # Province-clustered robust covariance matrix
  # name     = Chinese province name used for fixed effects/DLNM grouping
  # province = English province name used as cluster ID
  # -----------------------------------------------------------------------
  cluster_id <- model_data$province
  n_clusters <- length(unique(cluster_id))
  
  if (n_clusters < 2) {
    stop("Fewer than two province clusters are available.")
  }
  
  vcov_cluster <- sandwich::vcovCL(
    fit,
    cluster = cluster_id,
    type = "HC1",
    cadjust = TRUE,
    fix = TRUE
  )
  
  # Original/model-based crosspred
  run_crosspred_model <- function(basis_name, cen_value, at_values) {
    
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
  
  # Province-clustered crosspred
  run_crosspred_cluster <- function(basis_name, cen_value, at_values) {
    
    basis_obj <- get(
      basis_name,
      envir = model_env
    )
    
    coef_names <- names(coef(fit))
    coef_idx <- which(startsWith(coef_names, basis_name))
    
    if (length(coef_idx) == 0) {
      stop("No fitted coefficients found for cross-basis: ", basis_name)
    }
    
    beta_cb <- coef(fit)[coef_idx]
    vcov_cb <- vcov_cluster[
      coef_idx,
      coef_idx,
      drop = FALSE
    ]
    
    if (length(beta_cb) != ncol(basis_obj)) {
      stop(
        "Coefficient count mismatch for ", basis_name,
        ": fitted coefficients = ", length(beta_cb),
        ", cross-basis columns = ", ncol(basis_obj)
      )
    }
    
    if (any(!is.finite(beta_cb)) || any(!is.finite(vcov_cb))) {
      stop("Non-finite clustered coefficient/covariance values for ", basis_name)
    }
    
    crosspred(
      basis_obj,
      coef = beta_cb,
      vcov = vcov_cb,
      model.link = "log",
      cen = cen_value,
      at = at_values
    )
  }
  
  model_result <- data.frame()
  
  for (var in analysis_vars) {
    
    cb_name <- cb_names[[var]]
    x <- model_data[[var]]
    x_range <- range(x, na.rm = TRUE)
    
    cen_value <- as.numeric(center_map_used[[var]])
    increment <- as.numeric(increment_map[[var]])
    target_value <- cen_value + increment
    
    # ---------------------------------------------------------------
    # P value 1: Overall cross-basis P
    # ---------------------------------------------------------------
    model_overall_p <- get_overall_p(
      model = fit,
      term_name = cb_name
    )
    
    cluster_overall_p <- get_cluster_overall_p(
      model = fit,
      term_name = cb_name,
      vcov_cluster = vcov_cluster,
      n_clusters = n_clusters
    )
    
    rr_model <- c(
      RR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      LogRR = NA_real_,
      SE_LogRR = NA_real_,
      RR_contrast_P = NA_real_
    )
    
    rr_cluster <- c(
      RR = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      LogRR = NA_real_,
      SE_LogRR = NA_real_,
      RR_contrast_P = NA_real_
    )
    
    note_text <- ""
    
    if (any(!is.finite(x_range)) || diff(x_range) == 0) {
      note_text <- "Invalid exposure range."
    } else if (target_value > x_range[2]) {
      note_text <- "Target outside observed range."
    } else {
      
      # -------------------------------------------------------------
      # Original/model-based RR, CI and RR contrast P
      # -------------------------------------------------------------
      rr_model <- tryCatch(
        {
          pred_model <- run_crosspred_model(
            basis_name = cb_name,
            cen_value = cen_value,
            at_values = sort(unique(c(cen_value, target_value)))
          )
          
          idx_model <- which.min(
            abs(pred_model$predvar - target_value)
          )
          
          model_log_rr <- as.numeric(pred_model$allfit[idx_model])
          model_se_log_rr <- as.numeric(pred_model$allse[idx_model])
          model_contrast_p <- calc_rr_contrast_p(
            log_rr = model_log_rr,
            se_log_rr = model_se_log_rr
          )
          
          c(
            RR = as.numeric(pred_model$allRRfit[idx_model]),
            CI_low = as.numeric(pred_model$allRRlow[idx_model]),
            CI_high = as.numeric(pred_model$allRRhigh[idx_model]),
            LogRR = model_log_rr,
            SE_LogRR = model_se_log_rr,
            RR_contrast_P = model_contrast_p
          )
        },
        error = function(e) {
          note_text <<- paste0(
            note_text,
            ifelse(nchar(note_text) > 0, " | ", ""),
            "Model-based prediction: ",
            conditionMessage(e)
          )
          
          c(
            RR = NA_real_,
            CI_low = NA_real_,
            CI_high = NA_real_,
            LogRR = NA_real_,
            SE_LogRR = NA_real_,
            RR_contrast_P = NA_real_
          )
        }
      )
      
      # -------------------------------------------------------------
      # Province-clustered RR, CI and RR contrast P
      # -------------------------------------------------------------
      rr_cluster <- tryCatch(
        {
          pred_cluster <- run_crosspred_cluster(
            basis_name = cb_name,
            cen_value = cen_value,
            at_values = sort(unique(c(cen_value, target_value)))
          )
          
          idx_cluster <- which.min(
            abs(pred_cluster$predvar - target_value)
          )
          
          cluster_log_rr <- as.numeric(pred_cluster$allfit[idx_cluster])
          cluster_se_log_rr <- as.numeric(pred_cluster$allse[idx_cluster])
          cluster_contrast_p <- calc_rr_contrast_p(
            log_rr = cluster_log_rr,
            se_log_rr = cluster_se_log_rr
          )
          
          c(
            RR = as.numeric(pred_cluster$allRRfit[idx_cluster]),
            CI_low = as.numeric(pred_cluster$allRRlow[idx_cluster]),
            CI_high = as.numeric(pred_cluster$allRRhigh[idx_cluster]),
            LogRR = cluster_log_rr,
            SE_LogRR = cluster_se_log_rr,
            RR_contrast_P = cluster_contrast_p
          )
        },
        error = function(e) {
          note_text <<- paste0(
            note_text,
            ifelse(nchar(note_text) > 0, " | ", ""),
            "Clustered prediction: ",
            conditionMessage(e)
          )
          
          c(
            RR = NA_real_,
            CI_low = NA_real_,
            CI_high = NA_real_,
            LogRR = NA_real_,
            SE_LogRR = NA_real_,
            RR_contrast_P = NA_real_
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
        AR_term = ifelse(
          include_ar,
          "With lag.value1",
          "Without lag.value1"
        ),
        N_used = nobs(fit),
        Pseudo_R2 = round(pseudo_r2, 4),
        Dispersion = round(dispersion, 4),
        Variable = unname(display_names[[var]]),
        Exposure = var,
        Increment = unname(increment_label_map[[var]]),
        Cluster_variable = "province",
        N_clusters = n_clusters,
        
        # -----------------------------------------------------------
        # Original/model-based inference, retained for comparison
        # -----------------------------------------------------------
        Model_Overall_crossbasis_P = as.numeric(model_overall_p),
        Model_Overall_crossbasis_P_fmt = format_p(model_overall_p),
        Model_RR = as.numeric(rr_model["RR"]),
        Model_CI_low = as.numeric(rr_model["CI_low"]),
        Model_CI_high = as.numeric(rr_model["CI_high"]),
        Model_RR_95CI = format_rr(
          rr = as.numeric(rr_model["RR"]),
          low = as.numeric(rr_model["CI_low"]),
          high = as.numeric(rr_model["CI_high"])
        ),
        Model_RR_contrast_P = as.numeric(rr_model["RR_contrast_P"]),
        Model_RR_contrast_P_fmt = format_p(
          as.numeric(rr_model["RR_contrast_P"])
        ),
        
        # -----------------------------------------------------------
        # Province-clustered inference
        # -----------------------------------------------------------
        Cluster_Overall_crossbasis_P = as.numeric(cluster_overall_p),
        Cluster_Overall_crossbasis_P_fmt = format_p(cluster_overall_p),
        Cluster_RR = as.numeric(rr_cluster["RR"]),
        Cluster_CI_low = as.numeric(rr_cluster["CI_low"]),
        Cluster_CI_high = as.numeric(rr_cluster["CI_high"]),
        Cluster_RR_95CI = format_rr(
          rr = as.numeric(rr_cluster["RR"]),
          low = as.numeric(rr_cluster["CI_low"]),
          high = as.numeric(rr_cluster["CI_high"])
        ),
        Cluster_RR_contrast_P = as.numeric(rr_cluster["RR_contrast_P"]),
        Cluster_RR_contrast_P_fmt = format_p(
          as.numeric(rr_cluster["RR_contrast_P"])
        ),
        
        # -----------------------------------------------------------
        # Reviewer-facing primary columns = clustered inference
        # -----------------------------------------------------------
        RR = as.numeric(rr_cluster["RR"]),
        CI_low = as.numeric(rr_cluster["CI_low"]),
        CI_high = as.numeric(rr_cluster["CI_high"]),
        RR_95CI = format_rr(
          rr = as.numeric(rr_cluster["RR"]),
          low = as.numeric(rr_cluster["CI_low"]),
          high = as.numeric(rr_cluster["CI_high"])
        ),
        Overall_crossbasis_P = as.numeric(cluster_overall_p),
        Overall_crossbasis_P_fmt = format_p(cluster_overall_p),
        RR_contrast_P = as.numeric(rr_cluster["RR_contrast_P"]),
        RR_contrast_P_fmt = format_p(
          as.numeric(rr_cluster["RR_contrast_P"])
        ),
        
        Inference = paste0(
          "Province-clustered robust SE (HC1 + cluster adjustment); ",
          "Overall P = joint cross-basis Wald-F; ",
          "Contrast P = reported cumulative RR contrast"
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

# -------------------------------------------------------------------------
# OUTPUT 1: complete long table
# -------------------------------------------------------------------------
write.csv(
  all_results,
  file.path(
    output_dir,
    "sensitivity_results_long_two_Pvalues.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# -------------------------------------------------------------------------
# OUTPUT 2: province-clustered RR (95% CI), wide format
# -------------------------------------------------------------------------
rr_cluster_wide <- reshape(
  all_results[, c(
    "Setting",
    "Year_filter",
    "Variable",
    "Cluster_RR_95CI"
  )],
  idvar = c("Setting", "Year_filter"),
  timevar = "Variable",
  direction = "wide"
)

names(rr_cluster_wide) <- gsub(
  "^Cluster_RR_95CI\\.",
  "",
  names(rr_cluster_wide)
)

write.csv(
  rr_cluster_wide,
  file.path(
    output_dir,
    "sensitivity_RR_95CI_province_clustered_wide.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# -------------------------------------------------------------------------
# OUTPUT 3: Overall cross-basis P, wide format
# Joint significance of the entire exposure-lag cross-basis.
# -------------------------------------------------------------------------
overall_crossbasis_p_wide <- reshape(
  all_results[, c(
    "Setting",
    "Year_filter",
    "Variable",
    "Cluster_Overall_crossbasis_P_fmt"
  )],
  idvar = c("Setting", "Year_filter"),
  timevar = "Variable",
  direction = "wide"
)

names(overall_crossbasis_p_wide) <- gsub(
  "^Cluster_Overall_crossbasis_P_fmt\\.",
  "",
  names(overall_crossbasis_p_wide)
)

write.csv(
  overall_crossbasis_p_wide,
  file.path(
    output_dir,
    "sensitivity_Overall_crossbasis_P_wide.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# -------------------------------------------------------------------------
# OUTPUT 4: RR contrast P, wide format
# Corresponds to the specific cumulative RR (95% CI) shown in the table.
# -------------------------------------------------------------------------
rr_contrast_p_wide <- reshape(
  all_results[, c(
    "Setting",
    "Year_filter",
    "Variable",
    "Cluster_RR_contrast_P_fmt"
  )],
  idvar = c("Setting", "Year_filter"),
  timevar = "Variable",
  direction = "wide"
)

names(rr_contrast_p_wide) <- gsub(
  "^Cluster_RR_contrast_P_fmt\\.",
  "",
  names(rr_contrast_p_wide)
)

write.csv(
  rr_contrast_p_wide,
  file.path(
    output_dir,
    "sensitivity_RR_contrast_P_wide.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# -------------------------------------------------------------------------
# OUTPUT 5: direct model-based vs clustered comparison
# -------------------------------------------------------------------------
comparison_table <- all_results[, c(
  "Setting",
  "Year_filter",
  "Variable",
  "Increment",
  "N_used",
  "N_clusters",
  "Dispersion",
  "Model_RR_95CI",
  "Model_Overall_crossbasis_P_fmt",
  "Model_RR_contrast_P_fmt",
  "Cluster_RR_95CI",
  "Cluster_Overall_crossbasis_P_fmt",
  "Cluster_RR_contrast_P_fmt"
)]

write.csv(
  comparison_table,
  file.path(
    output_dir,
    "sensitivity_model_vs_clustered_two_Pvalues.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# -------------------------------------------------------------------------
# OUTPUT 6: detailed main-model table
# -------------------------------------------------------------------------
main_model_clustered <- all_results[
  all_results$Setting == "Main model",
  c(
    "Variable",
    "Increment",
    "N_used",
    "N_clusters",
    "Dispersion",
    "Cluster_RR",
    "Cluster_CI_low",
    "Cluster_CI_high",
    "Cluster_RR_95CI",
    "Cluster_Overall_crossbasis_P",
    "Cluster_Overall_crossbasis_P_fmt",
    "Cluster_RR_contrast_P",
    "Cluster_RR_contrast_P_fmt",
    "Model_RR_95CI",
    "Model_Overall_crossbasis_P_fmt",
    "Model_RR_contrast_P_fmt"
  ),
  drop = FALSE
]

write.csv(
  main_model_clustered,
  file.path(
    output_dir,
    "MAIN_MODEL_clustered_RR_and_two_Pvalues.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# -------------------------------------------------------------------------
# OUTPUT 7: compact table for manuscript/reviewer response
# -------------------------------------------------------------------------
main_model_for_paper <- all_results[
  all_results$Setting == "Main model",
  c(
    "Variable",
    "Increment",
    "RR_95CI",
    "Overall_crossbasis_P_fmt",
    "RR_contrast_P_fmt"
  ),
  drop = FALSE
]

names(main_model_for_paper) <- c(
  "Variable",
  "Increment",
  "RR_95CI",
  "Overall_crossbasis_P",
  "RR_contrast_P"
)

write.csv(
  main_model_for_paper,
  file.path(
    output_dir,
    "MAIN_MODEL_for_manuscript.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

# -------------------------------------------------------------------------
# OUTPUT 8: model specification table
# -------------------------------------------------------------------------
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
    "N_clusters",
    "Pseudo_R2",
    "Dispersion",
    "Inference"
  )]
)

write.csv(
  formula_table,
  file.path(
    output_dir,
    "sensitivity_model_formulas_with_cluster.csv"
  ),
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

cat("\n============================================================\n")
cat("Sensitivity analysis completed.\n")
cat("Fixed effects / DLNM grouping: name (Chinese province name)\n")
cat("Cluster ID: province (English province name)\n")
cat("Cluster covariance: HC1 + cluster adjustment\n")
cat("\nTwo different P values:\n")
cat("1) Overall_crossbasis_P = joint test of entire DLNM cross-basis\n")
cat("2) RR_contrast_P       = P value corresponding to RR (95% CI)\n")
cat("============================================================\n\n")

cat("MAIN MODEL - manuscript/reviewer table:\n\n")
print(main_model_for_paper, row.names = FALSE)

cat("\nDetailed MAIN MODEL results:\n\n")
print(main_model_clustered, row.names = FALSE)

cat("\nOutput files:\n")
cat("1. sensitivity_results_long_two_Pvalues.csv\n")
cat("2. sensitivity_RR_95CI_province_clustered_wide.csv\n")
cat("3. sensitivity_Overall_crossbasis_P_wide.csv\n")
cat("4. sensitivity_RR_contrast_P_wide.csv\n")
cat("5. sensitivity_model_vs_clustered_two_Pvalues.csv\n")
cat("6. MAIN_MODEL_clustered_RR_and_two_Pvalues.csv\n")
cat("7. MAIN_MODEL_for_manuscript.csv\n")
cat("8. sensitivity_model_formulas_with_cluster.csv\n")

if (nrow(error_log) > 0) {
  cat("9. sensitivity_error_log.csv\n")
}

cat("\nAll output files saved in: ", output_dir, "\n", sep = "")
