library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

has_sandwich <- requireNamespace("sandwich", quietly = TRUE)
has_MASS <- requireNamespace("MASS", quietly = TRUE)


# ============================================================
# 1. File paths
# ============================================================

# IMPORTANT:
# This must be the SAME nationwide 31-province file used for the main model.
national_file <- "E:/scarletfever_2013-2020_2.csv"

# These two files are used ONLY to identify the original North/South grouping.
north_file <- "E:/scarletfever_2013-2020_north.csv"
south_file <- "E:/scarletfever_2013-2020_south.csv"

output_dir <- "E:/Scarlet Fever/North_South_interaction"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ============================================================
# 2. Variables
# ============================================================

pollution_vars <- c("SO2", "CO", "NO2", "O3_8h")
meteo_vars <- c("rain", "sunlight", "humi", "meantemp")

required_vars <- c(
  "name", "year", "month",
  "beta", "population",
  pollution_vars,
  meteo_vars
)

numeric_vars <- c(
  "year", "month",
  "beta", "population",
  pollution_vars,
  meteo_vars
)

missing_tokens <- c(
  "NA", "N/A", "NaN", "nan",
  "missing", "Missing",
  "", " ",
  "-", "--", "—"
)


# ============================================================
# 3. Helper: read a CSV as character first
# ============================================================

read_csv_character <- function(file_path) {

  df <- read.csv(
    file_path,
    fileEncoding = "GB18030",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )

  # Remove completely empty columns
  is_empty_column <- function(x) {
    x_chr <- trimws(as.character(x))
    all(is.na(x_chr) | x_chr == "")
  }

  df <- df[, !vapply(df, is_empty_column, logical(1)), drop = FALSE]

  return(df)
}


# ============================================================
# 4. Read nationwide data
# ============================================================

data_all <- read_csv_character(national_file)

missing_vars <- setdiff(required_vars, names(data_all))

if (length(missing_vars) > 0) {
  stop(
    "Nationwide file is missing required variables: ",
    paste(missing_vars, collapse = ", ")
  )
}

data_all <- data_all[, required_vars, drop = FALSE]


# ============================================================
# 5. Recover the EXACT old North/South grouping
# ============================================================

north_map <- read_csv_character(north_file)
south_map <- read_csv_character(south_file)

if (!("name" %in% names(north_map))) {
  stop("The old North file does not contain a 'name' column.")
}

if (!("name" %in% names(south_map))) {
  stop("The old South file does not contain a 'name' column.")
}

north_provinces <- unique(trimws(as.character(north_map$name)))
south_provinces <- unique(trimws(as.character(south_map$name)))

north_provinces <- north_provinces[
  !is.na(north_provinces) & north_provinces != ""
]

south_provinces <- south_provinces[
  !is.na(south_provinces) & south_provinces != ""
]

# A province must not occur in both groups
overlap_provinces <- intersect(north_provinces, south_provinces)

if (length(overlap_provinces) > 0) {
  stop(
    "These provinces occur in BOTH the old North and South files: ",
    paste(overlap_provinces, collapse = ", ")
  )
}

cat("\n====================================================\n")
cat("ORIGINAL NORTH/SOUTH GROUPING\n")
cat("====================================================\n")

cat("\nNorth provinces (", length(north_provinces), "):\n", sep = "")
print(north_provinces)

cat("\nSouth provinces (", length(south_provinces), "):\n", sep = "")
print(south_provinces)


# ============================================================
# 6. Add regional stratum to nationwide data
# ============================================================

data_all$name <- trimws(as.character(data_all$name))

data_all$region <- ifelse(
  data_all$name %in% north_provinces,
  "North",
  ifelse(
    data_all$name %in% south_provinces,
    "South",
    NA_character_
  )
)

# Check whether every province in the nationwide dataset was mapped
unmapped_provinces <- unique(
  data_all$name[is.na(data_all$region)]
)

unmapped_provinces <- unmapped_provinces[
  !is.na(unmapped_provinces) & unmapped_provinces != ""
]

if (length(unmapped_provinces) > 0) {
  stop(
    "The following nationwide provinces were NOT found in either ",
    "the old North or South grouping file: ",
    paste(unmapped_provinces, collapse = ", ")
  )
}

# South = reference; North = 1
data_all$region <- factor(
  data_all$region,
  levels = c("South", "North")
)

data_all$North <- ifelse(
  data_all$region == "North",
  1,
  0
)


# ============================================================
# 7. Convert numeric variables
# ============================================================

count_missing_like <- function(x) {

  x_chr <- trimws(as.character(x))

  sum(
    is.na(x_chr) |
      x_chr %in% missing_tokens
  )
}


na_conversion_log <- data.frame(
  variable = numeric_vars,
  NA_before_conversion = NA_integer_,
  NA_after_conversion = NA_integer_,
  newly_created_NA = NA_integer_,
  stringsAsFactors = FALSE
)


for (i in seq_along(numeric_vars)) {

  var_name <- numeric_vars[i]

  before_n <- count_missing_like(
    data_all[[var_name]]
  )

  x_chr <- trimws(
    as.character(data_all[[var_name]])
  )

  x_chr[x_chr %in% missing_tokens] <- NA_character_

  # Remove thousands separators if any
  x_chr <- gsub(
    ",",
    "",
    x_chr,
    fixed = TRUE
  )

  x_num <- suppressWarnings(
    as.numeric(x_chr)
  )

  after_n <- sum(
    is.na(x_num)
  )

  data_all[[var_name]] <- x_num

  na_conversion_log$NA_before_conversion[i] <- before_n
  na_conversion_log$NA_after_conversion[i] <- after_n
  na_conversion_log$newly_created_NA[i] <- after_n - before_n
}


cat("\n====================================================\n")
cat("NUMERIC CONVERSION CHECK\n")
cat("====================================================\n")

print(na_conversion_log)


if (any(na_conversion_log$newly_created_NA > 0)) {
  warning(
    "Some previously unrecognized non-numeric tokens became NA. ",
    "Inspect NA_conversion_log.csv."
  )
}


# ============================================================
# 8. Validate nationwide data and regional mapping
# ============================================================

if (any(is.na(data_all$name) | data_all$name == "")) {
  stop("Some observations have missing province names.")
}

if (any(is.na(data_all$year))) {
  stop("year contains missing/non-numeric values.")
}

if (any(is.na(data_all$month))) {
  stop("month contains missing/non-numeric values.")
}

if (any(!data_all$month %in% 1:12)) {
  stop("month contains values outside 1-12.")
}

if (any(is.na(data_all$beta))) {
  stop("beta contains missing values. beta should not be imputed.")
}

if (any(data_all$beta < 0, na.rm = TRUE)) {
  stop(
    "beta contains negative values. ",
    "The current quasi-Poisson log-link model is not suitable."
  )
}

if (any(is.na(data_all$population))) {
  stop("population contains missing values.")
}

if (any(data_all$population <= 0, na.rm = TRUE)) {
  stop("population contains zero or negative values.")
}


# Duplicate province-month check
dup_key <- duplicated(
  data_all[, c("name", "year", "month")]
)

if (any(dup_key)) {

  print(
    data_all[
      dup_key,
      c("name", "year", "month", "region"),
      drop = FALSE
    ]
  )

  stop(
    "Duplicate province-year-month observations were found."
  )
}


# Each province must belong to exactly one region
strata_per_province <- tapply(
  as.character(data_all$region),
  data_all$name,
  function(x) length(unique(x))
)

if (any(strata_per_province != 1)) {
  stop(
    "At least one province belongs to more than one regional stratum."
  )
}


# Sort by province and month
data_all <- data_all[
  order(
    data_all$name,
    data_all$year,
    data_all$month
  ),
  ,
  drop = FALSE
]

rownames(data_all) <- NULL


# Province count check
province_region_table <- unique(
  data_all[, c("name", "region")]
)

cat("\n====================================================\n")
cat("REGIONAL MAPPING CHECK\n")
cat("====================================================\n")

cat("\nProvince counts by region:\n")
print(
  table(province_region_table$region)
)

cat("\nObservation counts by region:\n")
print(
  table(data_all$region)
)

cat(
  "\nTotal unique provinces: ",
  length(unique(data_all$name)),
  "\n",
  sep = ""
)

cat(
  "Total observations: ",
  nrow(data_all),
  "\n",
  sep = ""
)


if (length(unique(data_all$name)) != 31) {
  warning(
    "Nationwide data do not contain exactly 31 unique provinces."
  )
}


# Check expected 96 months/province
months_per_province <- table(
  data_all$name
)

if (any(months_per_province != 96)) {

  warning(
    "At least one province does not have exactly 96 monthly records."
  )

  print(
    months_per_province[months_per_province != 96]
  )
}


# ============================================================
# 9. Recalculate lag.value1 and within-province time sequence
# ============================================================

data_all <- data_all %>%
  group_by(name) %>%
  arrange(year, month, .by_group = TRUE) %>%
  mutate(
    lag.value1 = dplyr::lag(
      beta,
      n = 1,
      default = NA_real_
    ),
    seq = dplyr::row_number()
  ) %>%
  ungroup()


# ============================================================
# 10. Meteorological imputation
#     Same province-specific approach as the original code
# ============================================================

safe_na_kalman <- function(x, global_median) {

  x <- as.numeric(x)

  if (!any(is.na(x))) {
    return(x)
  }

  if (all(is.na(x))) {
    return(
      rep(global_median, length(x))
    )
  }

  if (sum(!is.na(x)) < 3) {

    local_median <- median(
      x,
      na.rm = TRUE
    )

    if (!is.finite(local_median)) {
      local_median <- global_median
    }

    x[is.na(x)] <- local_median

    return(x)
  }

  result <- tryCatch(
    {
      imputeTS::na_kalman(
        x,
        model = "StructTS",
        smooth = TRUE
      )
    },
    error = function(e1) {

      tryCatch(
        {
          imputeTS::na_interpolation(
            x,
            option = "linear"
          )
        },
        error = function(e2) {

          local_median <- median(
            x,
            na.rm = TRUE
          )

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

  return(
    as.numeric(result)
  )
}


cat("\n====================================================\n")
cat("METEOROLOGICAL IMPUTATION\n")
cat("====================================================\n")

cat("\nMissing values before imputation:\n")
print(
  sapply(
    data_all[meteo_vars],
    function(x) sum(is.na(x))
  )
)


for (var_name in meteo_vars) {

  new_var <- paste0(
    var_name,
    "_new"
  )

  global_median <- median(
    data_all[[var_name]],
    na.rm = TRUE
  )

  if (!is.finite(global_median)) {
    stop(
      var_name,
      " is completely missing and cannot be imputed."
    )
  }

  data_all[[new_var]] <- ave(
    data_all[[var_name]],
    data_all$name,
    FUN = function(x) {
      safe_na_kalman(
        x,
        global_median
      )
    }
  )
}


meteo_new_vars <- paste0(
  meteo_vars,
  "_new"
)


cat("\nMissing values after imputation:\n")
print(
  sapply(
    data_all[meteo_new_vars],
    function(x) sum(is.na(x))
  )
)


if (anyNA(data_all[meteo_new_vars])) {
  stop(
    "Meteorological imputation did not remove all missing values."
  )
}


# Convert to ordinary data.frame and create stable row IDs
data_all <- as.data.frame(data_all)

data_all$row_id <- seq_len(
  nrow(data_all)
)

rownames(data_all) <- as.character(
  data_all$row_id
)


# Save mapping and cleaned data for checking
write.csv(
  province_region_table,
  file.path(
    output_dir,
    "province_region_mapping.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

write.csv(
  na_conversion_log,
  file.path(
    output_dir,
    "NA_conversion_log.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

write.csv(
  data_all,
  file.path(
    output_dir,
    "nationwide_data_with_region_after_imputation.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)


# ============================================================
# 11. Build COMMON cross-bases from the nationwide dataset
#
# This is necessary for a direct North-South interaction test.
# ============================================================

get_valid_knots <- function(x) {

  x_valid <- x[
    is.finite(x)
  ]

  if (length(unique(x_valid)) < 3) {
    return(
      numeric(0)
    )
  }

  x_range <- range(
    x_valid,
    na.rm = TRUE
  )

  q <- quantile(
    x_valid,
    probs = c(0.1, 0.5, 0.9),
    na.rm = TRUE,
    names = FALSE
  )

  q <- unique(
    as.numeric(q)
  )

  q <- q[
    q > x_range[1] &
      q < x_range[2]
  ]

  return(q)
}


create_crossbasis <- function(
  var_name,
  df,
  lag_value = 3
) {

  x <- df[[var_name]]

  x_valid <- x[
    is.finite(x)
  ]

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

  var_knots <- get_valid_knots(
    x
  )

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

  cb <- dlnm::crossbasis(
    x,
    lag = lag_value,
    argvar = argvar_list,
    arglag = list(
      fun = "ns",
      knots = dlnm::logknots(
        lag_value,
        nk = 2
      )
    ),
    group = df$name
  )

  return(cb)
}


cb_so2 <- create_crossbasis(
  "SO2",
  data_all,
  lag_value = 3
)

cb_co <- create_crossbasis(
  "CO",
  data_all,
  lag_value = 3
)

cb_no2 <- create_crossbasis(
  "NO2",
  data_all,
  lag_value = 3
)

cb_o3 <- create_crossbasis(
  "O3_8h",
  data_all,
  lag_value = 3
)

cb_sun <- create_crossbasis(
  "sunlight_new",
  data_all,
  lag_value = 3
)

cb_humi <- create_crossbasis(
  "humi_new",
  data_all,
  lag_value = 3
)

cb_rain <- create_crossbasis(
  "rain_new",
  data_all,
  lag_value = 3
)

cb_temp <- create_crossbasis(
  "meantemp_new",
  data_all,
  lag_value = 3
)


basis_list <- list(
  SO2 = cb_so2,
  CO = cb_co,
  NO2 = cb_no2,
  O3 = cb_o3,
  Sunlight = cb_sun,
  Humidity = cb_humi,
  Precipitation = cb_rain,
  Temperature = cb_temp
)


# ============================================================
# 12. Base nationwide multivariable DLNM
#     Same structure as the main model
# ============================================================

base_model <- glm(
  beta ~
    cb_so2 +
    cb_co +
    cb_no2 +
    cb_o3 +
    cb_sun +
    cb_humi +
    cb_rain +
    cb_temp +
    splines::ns(seq, 8) +
    offset(log(population)) +
    factor(name) +
    lag.value1,
  family = quasipoisson(
    link = "log"
  ),
  data = data_all,
  na.action = na.exclude
)


base_n <- stats::nobs(
  base_model
)


cat("\n====================================================\n")
cat("BASE NATIONWIDE MODEL CHECK\n")
cat("====================================================\n")

cat(
  "N used = ",
  base_n,
  "\n",
  sep = ""
)

cat(
  "Residual df = ",
  base_model$df.residual,
  "\n",
  sep = ""
)

cat(
  "Dispersion = ",
  round(
    summary(base_model)$dispersion,
    4
  ),
  "\n",
  sep = ""
)


capture.output(
  summary(base_model),
  file = file.path(
    output_dir,
    "base_nationwide_multivariable_model_summary.txt"
  )
)


# ============================================================
# 13. Primary interaction test:
#     nested quasi-Poisson F test
# ============================================================

nested_interaction_F_test <- function(
  base_model,
  interaction_model
) {

  if (
    stats::nobs(base_model) !=
      stats::nobs(interaction_model)
  ) {
    stop(
      "Base and interaction models use different sample sizes."
    )
  }

  a <- anova(
    base_model,
    interaction_model,
    test = "F"
  )

  if (nrow(a) < 2) {
    stop(
      "Nested model comparison did not return two rows."
    )
  }

  f_col <- grep(
    "^F$|F value",
    names(a),
    value = TRUE
  )

  p_col <- grep(
    "Pr\\(>F\\)",
    names(a),
    value = TRUE
  )

  if (length(f_col) == 0) {
    f_value <- NA_real_
  } else {
    f_value <- as.numeric(
      a[2, f_col[1]]
    )
  }

  if (length(p_col) == 0) {
    p_value <- NA_real_
  } else {
    p_value <- as.numeric(
      a[2, p_col[1]]
    )
  }

  if ("Df" %in% names(a)) {
    df_added <- abs(
      as.numeric(
        a[2, "Df"]
      )
    )
  } else {
    df_added <- NA_real_
  }

  return(
    c(
      df = df_added,
      F = f_value,
      P = p_value
    )
  )
}


# ============================================================
# 14. Optional province-clustered robust Wald test
# ============================================================

get_used_clusters <- function(
  model,
  df
) {

  used_rows <- rownames(
    model.frame(model)
  )

  idx <- match(
    used_rows,
    rownames(df)
  )

  if (any(is.na(idx))) {
    stop(
      "Could not align fitted observations with province clusters."
    )
  }

  return(
    factor(
      df$name[idx]
    )
  )
}


clustered_interaction_Wald_test <- function(
  model,
  interaction_prefix,
  df
) {

  if (!has_sandwich) {

    return(
      c(
        df = NA_real_,
        F = NA_real_,
        P = NA_real_,
        clusters = NA_real_
      )
    )
  }

  b_all <- coef(
    model
  )

  coef_names <- names(
    b_all
  )

  idx_coef <- grep(
    paste0(
      "^",
      interaction_prefix
    ),
    coef_names
  )

  if (length(idx_coef) == 0) {

    warning(
      "No interaction coefficients found with prefix: ",
      interaction_prefix
    )

    return(
      c(
        df = NA_real_,
        F = NA_real_,
        P = NA_real_,
        clusters = NA_real_
      )
    )
  }

  cluster_id <- get_used_clusters(
    model,
    df
  )

  G <- nlevels(
    cluster_id
  )

  V_all <- tryCatch(
    sandwich::vcovCL(
      model,
      cluster = cluster_id,
      type = "HC1"
    ),
    error = function(e) {
      warning(
        "vcovCL failed: ",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(V_all)) {

    return(
      c(
        df = NA_real_,
        F = NA_real_,
        P = NA_real_,
        clusters = G
      )
    )
  }

  b <- b_all[
    idx_coef
  ]

  V <- V_all[
    idx_coef,
    idx_coef,
    drop = FALSE
  ]

  keep <- is.finite(
    b
  )

  if (nrow(V) > 0) {

    keep <- keep &
      apply(
        V,
        1,
        function(z) all(is.finite(z))
      )
  }

  b <- b[
    keep
  ]

  V <- V[
    keep,
    keep,
    drop = FALSE
  ]

  if (
    length(b) == 0 ||
      nrow(V) == 0
  ) {

    return(
      c(
        df = NA_real_,
        F = NA_real_,
        P = NA_real_,
        clusters = G
      )
    )
  }

  V <- (
    V + t(V)
  ) / 2

  q_eff <- qr(
    V,
    tol = 1e-10
  )$rank

  if (
    !is.finite(q_eff) ||
      q_eff < 1
  ) {

    return(
      c(
        df = NA_real_,
        F = NA_real_,
        P = NA_real_,
        clusters = G
      )
    )
  }

  V_inv <- tryCatch(
    solve(V),
    error = function(e) {

      if (has_MASS) {
        MASS::ginv(V)
      } else {
        qr.solve(
          V,
          diag(nrow(V)),
          tol = 1e-10
        )
      }
    }
  )

  W <- as.numeric(
    t(b) %*%
      V_inv %*%
      b
  )

  if (!is.finite(W) || W < 0) {

    return(
      c(
        df = q_eff,
        F = NA_real_,
        P = NA_real_,
        clusters = G
      )
    )
  }

  F_stat <- W / q_eff

  p_value <- stats::pf(
    F_stat,
    df1 = q_eff,
    df2 = G - 1,
    lower.tail = FALSE
  )

  return(
    c(
      df = q_eff,
      F = F_stat,
      P = p_value,
      clusters = G
    )
  )
}


format_p <- function(p) {

  if (
    length(p) == 0 ||
      is.na(p) ||
      !is.finite(p)
  ) {
    return(NA_character_)
  }

  if (p < 0.001) {
    return("<0.001")
  }

  return(
    sprintf(
      "%.3f",
      p
    )
  )
}


# ============================================================
# 15. Fit one interaction model per environmental variable
# ============================================================

interaction_results <- data.frame(
  Variable = character(),
  N_used = integer(),
  N_provinces = integer(),

  Interaction_df = numeric(),
  F_interaction = numeric(),
  P_interaction = numeric(),
  P_interaction_fmt = character(),

  Cluster_df = numeric(),
  Cluster_F = numeric(),
  P_interaction_cluster = numeric(),
  P_interaction_cluster_fmt = character(),

  Interpretation = character(),

  stringsAsFactors = FALSE
)


interaction_models <- list()


for (var_label in names(basis_list)) {

  cat("\n====================================================\n")
  cat(
    "Testing North-South interaction for: ",
    var_label,
    "\n",
    sep = ""
  )
  cat("====================================================\n")

  target_basis <- basis_list[[var_label]]

  # South = 0 -> interaction contribution = 0
  # North = 1 -> interaction contribution = target cross-basis
  int_target <- sweep(
    as.matrix(target_basis),
    MARGIN = 1,
    STATS = data_all$North,
    FUN = "*"
  )

  colnames(int_target) <- paste0(
    "b",
    seq_len(
      ncol(int_target)
    )
  )


  interaction_model <- glm(
    beta ~
      cb_so2 +
      cb_co +
      cb_no2 +
      cb_o3 +
      cb_sun +
      cb_humi +
      cb_rain +
      cb_temp +
      int_target +
      splines::ns(seq, 8) +
      offset(log(population)) +
      factor(name) +
      lag.value1,
    family = quasipoisson(
      link = "log"
    ),
    data = data_all,
    na.action = na.exclude
  )


  if (
    stats::nobs(interaction_model) !=
      base_n
  ) {
    stop(
      "Interaction model for ",
      var_label,
      " uses a different sample size from the base model."
    )
  }


  interaction_models[[var_label]] <- interaction_model


  # Primary formal interaction test
  test_standard <- nested_interaction_F_test(
    base_model = base_model,
    interaction_model = interaction_model
  )

  p_standard <- as.numeric(
    test_standard["P"]
  )


  # Optional province-clustered robust interaction test
  test_cluster <- clustered_interaction_Wald_test(
    model = interaction_model,
    interaction_prefix = "int_target",
    df = data_all
  )

  p_cluster <- as.numeric(
    test_cluster["P"]
  )


  interpretation <- if (
    is.finite(p_standard) &&
      p_standard < 0.05
  ) {

    "Evidence of North-South heterogeneity"

  } else if (
    is.finite(p_standard)
  ) {

    "No statistical evidence of North-South heterogeneity"

  } else {

    "Interaction test unavailable"
  }


  cluster_used <- get_used_clusters(
    interaction_model,
    data_all
  )


  interaction_results <- rbind(
    interaction_results,
    data.frame(
      Variable = var_label,

      N_used = stats::nobs(
        interaction_model
      ),

      N_provinces = nlevels(
        cluster_used
      ),

      Interaction_df = as.numeric(
        test_standard["df"]
      ),

      F_interaction = as.numeric(
        test_standard["F"]
      ),

      P_interaction = p_standard,

      P_interaction_fmt = format_p(
        p_standard
      ),

      Cluster_df = as.numeric(
        test_cluster["df"]
      ),

      Cluster_F = as.numeric(
        test_cluster["F"]
      ),

      P_interaction_cluster = p_cluster,

      P_interaction_cluster_fmt = format_p(
        p_cluster
      ),

      Interpretation = interpretation,

      stringsAsFactors = FALSE
    )
  )


  capture.output(
    summary(
      interaction_model
    ),
    file = file.path(
      output_dir,
      paste0(
        "interaction_model_",
        var_label,
        "_summary.txt"
      )
    )
  )
}


# ============================================================
# 16. Safety check: all 8 interaction models must finish
# ============================================================

if (
  nrow(interaction_results) !=
    length(basis_list)
) {
  stop(
    "Not all 8 interaction models were completed successfully."
  )
}


# ============================================================
# 17. Optional BH correction for the 8 interaction tests
# ============================================================

interaction_results$P_interaction_BH <- p.adjust(
  interaction_results$P_interaction,
  method = "BH"
)

interaction_results$P_interaction_BH_fmt <- vapply(
  interaction_results$P_interaction_BH,
  format_p,
  character(1)
)


if (
  any(
    is.finite(
      interaction_results$P_interaction_cluster
    )
  )
) {

  interaction_results$P_interaction_cluster_BH <- p.adjust(
    interaction_results$P_interaction_cluster,
    method = "BH"
  )

} else {

  interaction_results$P_interaction_cluster_BH <- rep(
    NA_real_,
    nrow(interaction_results)
  )
}


interaction_results$P_interaction_cluster_BH_fmt <- vapply(
  interaction_results$P_interaction_cluster_BH,
  format_p,
  character(1)
)


# ============================================================
# 18. Print and save final results
# ============================================================

cat("\n\n====================================================\n")
cat("FINAL NORTH-SOUTH CROSS-BASIS INTERACTION RESULTS\n")
cat("====================================================\n\n")

print(
  interaction_results
)


write.csv(
  interaction_results,
  file.path(
    output_dir,
    "North_South_crossbasis_interaction_tests_FULL.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)


# Compact reviewer-facing table
reviewer_table <- interaction_results[
  ,
  c(
    "Variable",
    "Interaction_df",
    "F_interaction",
    "P_interaction",
    "P_interaction_fmt",
    "P_interaction_cluster",
    "P_interaction_cluster_fmt",
    "Interpretation"
  ),
  drop = FALSE
]


write.csv(
  reviewer_table,
  file.path(
    output_dir,
    "North_South_crossbasis_interaction_reviewer_table.csv"
  ),
  row.names = FALSE,
  fileEncoding = "GB18030"
)


saveRDS(
  base_model,
  file.path(
    output_dir,
    "base_nationwide_multivariable_model.rds"
  )
)


saveRDS(
  interaction_models,
  file.path(
    output_dir,
    "North_South_interaction_models.rds"
  )
)


cat(
  "\nPrimary reviewer-facing result file:\n",
  file.path(
    output_dir,
    "North_South_crossbasis_interaction_reviewer_table.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nFull result file:\n",
  file.path(
    output_dir,
    "North_South_crossbasis_interaction_tests_FULL.csv"
  ),
  "\n",
  sep = ""
)

cat("\nAnalysis completed successfully.\n")
