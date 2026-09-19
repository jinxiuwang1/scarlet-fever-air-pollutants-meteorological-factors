library(dplyr)
library(imputeTS)

data_file <- "E:/scarletfever_2013-2020_2.csv"
output_dir <- "E:/Scarlet Fever/multicollinearity diagnosis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data <- read.csv(
  data_file,
  fileEncoding = "GB18030"
)

# Remove empty columns.
data[data == ""] <- NA
data <- data[, colSums(!is.na(data)) > 0]

pollution_vars <- c("PM2.5", "PM10", "SO2", "CO", "NO2", "O3_8h")
meteo_vars <- c("rain", "sunlight", "humi", "meantemp")

required_vars <- c(
  "name", "year", "month",
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

numeric_vars <- c("year", "month", pollution_vars, meteo_vars)

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

data <- data[with(data, order(name, year, month)), ]

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

meteo_imp_vars <- paste0(meteo_vars, "_imp")

cat("\nMissing values after meteorological imputation:\n")
print(sapply(data[meteo_imp_vars], function(x) sum(is.na(x))))

vif_vars <- c(
  "PM2.5", "PM10", "SO2", "CO", "NO2", "O3_8h",
  "meantemp_imp", "humi_imp", "rain_imp", "sunlight_imp"
)

vif_data <- data[, vif_vars]
vif_data <- vif_data[complete.cases(vif_data), ]

cat("\nRows used for VIF analysis: ", nrow(vif_data), "\n", sep = "")

if (nrow(vif_data) == 0) {
  stop("No complete rows remain for VIF analysis.")
}

calculate_vif <- function(df, vars) {
  
  out <- data.frame(
    Variable = vars,
    N = NA_integer_,
    R_squared = NA_real_,
    VIF = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(vars)) {
    
    y_var <- vars[i]
    x_vars <- setdiff(vars, y_var)
    
    tmp <- df[, c(y_var, x_vars), drop = FALSE]
    tmp <- tmp[complete.cases(tmp), , drop = FALSE]
    
    out$N[i] <- nrow(tmp)
    
    if (nrow(tmp) <= length(x_vars) + 1) {
      next
    }
    
    if (length(unique(tmp[[y_var]])) < 2) {
      next
    }
    
    formula_i <- as.formula(
      paste(y_var, "~", paste(x_vars, collapse = " + "))
    )
    
    fit_i <- tryCatch(
      {
        lm(formula_i, data = tmp)
      },
      error = function(e) NULL
    )
    
    if (is.null(fit_i)) {
      next
    }
    
    r2 <- summary(fit_i)$r.squared
    
    out$R_squared[i] <- r2
    
    if (!is.finite(r2)) {
      next
    }
    
    if (r2 >= 1) {
      out$VIF[i] <- Inf
    } else {
      out$VIF[i] <- 1 / (1 - r2)
    }
  }
  
  return(out)
}

vif_results <- calculate_vif(vif_data, vif_vars)

vif_results$Collinearity_level <- cut(
  vif_results$VIF,
  breaks = c(-Inf, 5, 10, Inf),
  labels = c("Acceptable", "Moderate", "High")
)

cat("\nVIF results:\n")
print(vif_results)

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

write.csv(
  vif_results,
  file.path(output_dir, "VIF_results.csv"),
  row.names = FALSE,
  fileEncoding = "GB18030"
)

cat("\nMulticollinearity diagnosis completed.\n")
cat("Results saved to: ", output_dir, "\n", sep = "")