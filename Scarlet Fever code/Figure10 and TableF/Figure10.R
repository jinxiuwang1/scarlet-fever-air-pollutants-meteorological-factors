library(dlnm)
library(splines)
library(dplyr)
library(imputeTS)

data_file <- "E:/scarletfever_2013-2020_2.csv"
work_dir <- "E:/Scarlet Fever/Figure10"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)

data <- read.csv(data_file, fileEncoding = "GB18030")
is_empty_column <- function(x) {
  x_chr <- trimws(as.character(x))
  all(is.na(x_chr) | x_chr == "")
}
data <- data[, !sapply(data, is_empty_column)]

pollution_vars <- c("NO2", "PM10", "O3_8h", "PM2.5", "SO2", "CO")
meteo_vars <- c("rain", "sunlight", "humi", "meantemp")
required_vars <- c("name", "year", "month", "beta", "population", pollution_vars, meteo_vars)
missing_vars <- setdiff(required_vars, names(data))
if (length(missing_vars) > 0) stop("Missing variables: ", paste(missing_vars, collapse=", "))
data <- data[, required_vars]

missing_tokens <- c("NA", "N/A", "NaN", "missing", "", " ")
count_missing_like <- function(x) {
  if (is.character(x) || is.factor(x)) {
    x_chr <- trimws(as.character(x))
    return(sum(is.na(x_chr) | x_chr %in% missing_tokens))
  } else return(sum(is.na(x)))
}
vars_to_numeric <- c("year", "month", "beta", "population", pollution_vars, meteo_vars)
na_conversion_log <- data.frame(variable=vars_to_numeric, NA_before=NA, NA_after=NA, new_NA=NA)
for (i in seq_along(vars_to_numeric)) {
  var <- vars_to_numeric[i]
  na_before <- count_missing_like(data[[var]])
  if (!is.numeric(data[[var]])) {
    x_chr <- trimws(as.character(data[[var]]))
    x_chr[x_chr %in% missing_tokens] <- NA
    x_chr <- gsub(",", "", x_chr)
    data[[var]] <- suppressWarnings(as.numeric(x_chr))
  }
  na_after <- sum(is.na(data[[var]]))
  na_conversion_log[i, 2:4] <- c(na_before, na_after, na_after - na_before)
}
if (any(is.na(data$beta))) stop("beta has NA")
if (any(data$beta<0, na.rm=TRUE)) stop("beta negative")
if (any(is.na(data$population))) stop("population NA")
if (any(data$population<=0, na.rm=TRUE)) stop("population <=0")

data <- data[order(data$name, data$year, data$month), ]
data <- data %>% group_by(name) %>% mutate(lag.value1 = dplyr::lag(beta, n=1, default=NA)) %>% ungroup()
data$seq <- 1:nrow(data)

safe_na_kalman <- function(x, global_median) {
  x <- as.numeric(x)
  if (!any(is.na(x))) return(x)
  if (all(is.na(x))) return(rep(global_median, length(x)))
  if (sum(!is.na(x)) < 3) {
    x[is.na(x)] <- median(x, na.rm=TRUE)
    if (any(is.na(x))) x[is.na(x)] <- global_median
    return(x)
  }
  tryCatch(na_kalman(x, model="StructTS", smooth=TRUE), 
           error=function(e) na_interpolation(x, option="linear"))
}
for (var in meteo_vars) {
  new_var <- paste0(var, "_new")
  global_median <- median(data[[var]], na.rm=TRUE)
  data[[new_var]] <- ave(data[[var]], data$name, FUN=function(x) safe_na_kalman(x, global_median))
}
setwd(work_dir)
write.csv(na_conversion_log, "NA_conversion_log.csv", row.names=FALSE, fileEncoding="GB18030")
write.csv(data, "data_after_meteorological_imputation.csv", row.names=FALSE, fileEncoding="GB18030")

calc_pseudo_r2 <- function(model) {
  if(is.null(model$deviance) || is.null(model$null.deviance)) return(NA)
  if(model$null.deviance==0) return(NA)
  1 - model$deviance/model$null.deviance
}
get_valid_knots <- function(x) {
  x_valid <- x[is.finite(x)]
  if(length(unique(x_valid))<3) return(numeric(0))
  q <- quantile(x_valid, probs=c(0.1,0.5,0.9), na.rm=TRUE)
  q <- unique(q)
  rng <- range(x_valid)
  q[q>rng[1] & q<rng[2]]
}
get_overall_p <- function(model, term="cb") {
  tryCatch(drop1(model, test="F")[term, "Pr(>F)"], error=function(e) NA_real_)
}

plot_dlnm_contour <- function(pred, file_name, xlab_value) {
  tryCatch({
    tiff(file.path(work_dir, file_name), width=3000, height=2500, res=300)
    on.exit(dev.off())
    plot(pred, "contour", xlab=xlab_value, ylab="Lag (months)",
         key.title=title("RR", cex.main=1.2),
         plot.axes={axis(1,cex.axis=1.8); axis(2,cex.axis=1.8)},
         key.axes=axis(4,cex.axis=1.8), main="Contour plot", cex.main=2.0, cex.lab=1.8)
    return(TRUE)
  }, error=function(e){warning("Contour plot failed: ", e); FALSE})
}

fit_dlnm_one <- function(var_name, display_name, file_name, xlab_value,
                         at_fun, increment, increment_label) {
  cat("\nProcessing ", display_name, "...\n", sep="")
  exposure <- data[[var_name]]
  x_range <- range(exposure, na.rm=TRUE)
  
  empty_cum <- data.frame(Variable=display_name, Exposure=var_name, Increment=increment_label,
                          Overall_P=NA, Center=NA, Target=NA, RR=NA, CI_low=NA, CI_high=NA,
                          RR_95CI=NA, Pseudo_R2=NA, Sig_inc_threshold=NA, Sig_dec_threshold=NA,
                          Figure=file_name, Note="", stringsAsFactors=FALSE)
  empty_peak <- data.frame(Variable=display_name, Exposure=var_name, Increment=increment_label,
                           Peak_exposure=NA, Peak_lag_months=NA, Peak_RR=NA,
                           CI_low=NA, CI_high=NA, Peak_RR_95CI=NA, stringsAsFactors=FALSE)
  if(any(!is.finite(x_range)) || diff(x_range)==0) {
    warning(display_name, " invalid range")
    return(list(cumulative=empty_cum, lag_specific=data.frame(), peak=empty_peak, model=NULL))
  }
  
  knots_var <- get_valid_knots(exposure)
  argvar <- if(length(knots_var)>=1) list(fun="ns", knots=knots_var, Boundary.knots=x_range) else list(fun="lin")
  lk <- logknots(3, nk=2)
  cb <- crossbasis(exposure, lag=3, argvar=argvar, arglag=list(fun="ns", knots=lk), group=data$name)
  
  model <- glm(beta ~ cb + ns(seq, 8*1) + offset(log(population)) + factor(name) + lag.value1,
               family=quasipoisson(), data=data, na.action=na.exclude)
  
  model_rows <- as.integer(rownames(model.frame(model)))
  if(any(is.na(model_rows))) model_rows <- seq_len(nrow(data))
  exposure_model <- data[[var_name]][model_rows]
  x_model_range <- range(exposure_model, na.rm=TRUE)
  cen_value <- median(exposure_model, na.rm=TRUE)
  target_value <- cen_value + increment
  
  at_plot <- tryCatch(at_fun(x_model_range), error=function(e) NULL)
  pred_plot <- NULL
  if(!is.null(at_plot) && length(at_plot)>=2) {
    pred_plot <- tryCatch(crosspred(cb, model, cen=cen_value, at=at_plot, bylag=0.2), error=function(e) NULL)
    if(!is.null(pred_plot)) {
      plot_dlnm_contour(pred_plot, file_name, xlab_value)
    } else {
      warning("Crosspred for contour failed for ", display_name)
    }
  } else {
    warning("Invalid at_values_plot for ", display_name)
  }
  
  pred_lag <- crosspred(cb, model, cen=cen_value, at=target_value, bylag=1)
  cum_RR <- pred_lag$allRRfit[1]; cum_low <- pred_lag$allRRlow[1]; cum_high <- pred_lag$allRRhigh[1]
  cum_text <- sprintf("%.3f (%.3f, %.3f)", cum_RR, cum_low, cum_high)
  lag_rr <- pred_lag$matRRfit[1,]; lag_low <- pred_lag$matRRlow[1,]; lag_high <- pred_lag$matRRhigh[1,]
  lag_text <- sprintf("%.3f (%.3f, %.3f)", lag_rr, lag_low, lag_high)
  
  if(!is.null(pred_plot)) {
    rr_mat <- pred_plot$matRRfit
    max_pos <- which(rr_mat == max(rr_mat, na.rm=TRUE), arr.ind=TRUE)[1, , drop=FALSE]
    peak_exposure <- pred_plot$predvar[max_pos[1,1]]
    peak_lag <- pred_plot$lag[max_pos[1,2]]
    peak_RR <- rr_mat[max_pos[1,1], max_pos[1,2]]
    peak_CI_low <- pred_plot$matRRlow[max_pos[1,1], max_pos[1,2]]
    peak_CI_high <- pred_plot$matRRhigh[max_pos[1,1], max_pos[1,2]]
  } else {
    peak_idx_int <- which.max(lag_rr)
    peak_exposure <- target_value
    peak_lag <- peak_idx_int - 1
    peak_RR <- lag_rr[peak_idx_int]
    peak_CI_low <- lag_low[peak_idx_int]
    peak_CI_high <- lag_high[peak_idx_int]
    cat("  -> WARNING: using integer lag peak\n")
  }
  
  peak_info <- data.frame(Variable=display_name, Exposure=var_name, Increment=increment_label,
                          Peak_exposure=as.numeric(peak_exposure), Peak_lag_months=as.numeric(peak_lag),
                          Peak_RR=as.numeric(peak_RR), CI_low=as.numeric(peak_CI_low),
                          CI_high=as.numeric(peak_CI_high),
                          Peak_RR_95CI=sprintf("%.3f (%.3f, %.3f)", peak_RR, peak_CI_low, peak_CI_high),
                          stringsAsFactors=FALSE)
  
  sig_inc <- sig_dec <- NA_real_
  if(!is.null(pred_plot)) {
    xv <- pred_plot$predvar
    low <- pred_plot$allRRlow
    high <- pred_plot$allRRhigh
    idx <- which(!is.na(low) & !is.na(high))
    if(length(idx)>0) {
      xv <- xv[idx]; low <- low[idx]; high <- high[idx]
      inc <- low > 1
      if(any(inc)) {
        rle_inc <- rle(inc)
        cum_len <- 0
        for(i in seq_along(rle_inc$values)) {
          if(rle_inc$values[i]) {
            sig_inc <- xv[cum_len+1]
            break
          }
          cum_len <- cum_len + rle_inc$lengths[i]
        }
      }
      dec <- high < 1
      if(any(dec)) {
        rle_dec <- rle(dec)
        cum_len <- 0
        last <- NA
        for(i in seq_along(rle_dec$values)) {
          if(rle_dec$values[i]) last <- cum_len + rle_dec$lengths[i]
          cum_len <- cum_len + rle_dec$lengths[i]
        }
        if(!is.na(last)) sig_dec <- xv[last]
      }
    }
  }
  
  pseudo_r2 <- calc_pseudo_r2(model)
  overall_p <- get_overall_p(model, "cb")
  cum_row <- data.frame(Variable=display_name, Exposure=var_name, Increment=increment_label,
                        Overall_P=as.numeric(overall_p), Center=cen_value, Target=target_value,
                        RR=cum_RR, CI_low=cum_low, CI_high=cum_high, RR_95CI=cum_text,
                        Pseudo_R2=pseudo_r2, Sig_inc_threshold=sig_inc, Sig_dec_threshold=sig_dec,
                        Figure=file_name, Note="", stringsAsFactors=FALSE)
  
  lag_rows <- data.frame(Variable=display_name, Exposure=var_name, Increment=increment_label,
                         Lag=0:3, RR=lag_rr, CI_low=lag_low, CI_high=lag_high, RR_95CI=lag_text,
                         stringsAsFactors=FALSE)
  
  return(list(cumulative=cum_row, lag_specific=lag_rows, peak=peak_info, model=model))
}

at_50 <- function(xr) seq(xr[1], xr[2], length.out=50)
at_rain <- function(xr) seq(0, xr[2], by=10)
at_temp <- function(xr) seq(floor(xr[1]), ceiling(xr[2]), by=1)

exposure_list <- list(
  list(var="NO2", display="NO2", file="FIGURE_5A_NO2.tiff", xlab=expression(NO[2]~"("*mu*g/m^3*")"), at_fun=at_50, increment=10, increment_label="NO2 (+10 ug/m³)"),
  list(var="PM10", display="PM10", file="FIGURE_5B_PM10.tiff", xlab=expression(PM[10]~"("*mu*g/m^3*")"), at_fun=at_50, increment=10, increment_label="PM10 (+10 ug/m³)"),
  list(var="O3_8h", display="O3", file="FIGURE_5C_O3.tiff", xlab=expression(O[3]~"("*mu*g/m^3*")"), at_fun=at_50, increment=10, increment_label="O3 (+10 ug/m³)"),
  list(var="PM2.5", display="PM2.5", file="FIGURE_5D_PM25.tiff", xlab=expression(PM[2.5]~"("*mu*g/m^3*")"), at_fun=at_50, increment=10, increment_label="PM2.5 (+10 ug/m³)"),
  list(var="SO2", display="SO2", file="FIGURE_5E_SO2.tiff", xlab=expression(SO[2]~"("*mu*g/m^3*")"), at_fun=at_50, increment=10, increment_label="SO2 (+10 ug/m³)"),
  list(var="CO", display="CO", file="FIGURE_5F_CO.tiff", xlab=expression(CO~"(mg/m"^3*")"), at_fun=at_50, increment=1, increment_label="CO (+1 mg/m³)"),
  list(var="sunlight_new", display="Sunlight", file="FIGURE_5G_Sunlight.tiff", xlab="Sunlight (hours)", at_fun=at_50, increment=5, increment_label="Sunlight (+5 h)"),
  list(var="humi_new", display="Humidity", file="FIGURE_5H_Humidity.tiff", xlab="Relative humidity (%)", at_fun=at_50, increment=10, increment_label="Relative humidity (+10%)"),
  list(var="rain_new", display="Rainfall", file="FIGURE_5I_Rainfall.tiff", xlab="Precipitation (mm)", at_fun=at_rain, increment=20, increment_label="Precipitation (+20 mm)"),
  list(var="meantemp_new", display="Mean temperature", file="FIGURE_5J_Temperature.tiff", xlab="Mean temperature (°C)", at_fun=at_temp, increment=5, increment_label="Mean temperature (+5 °C)")
)

cumulative_results <- data.frame()
lag_results_all <- data.frame()
peak_results_all <- data.frame()
fit_objects <- list()

for (item in exposure_list) {
  fit_result <- tryCatch({
    fit_dlnm_one(item$var, item$display, item$file, item$xlab, item$at_fun, item$increment, item$increment_label)
  }, error=function(e) {
    warning("Failed for ", item$display, ": ", e)
    empty_cum <- data.frame(Variable=item$display, Exposure=item$var, Increment=item$increment_label,
                            Overall_P=NA, Center=NA, Target=NA, RR=NA, CI_low=NA, CI_high=NA,
                            RR_95CI=NA, Pseudo_R2=NA, Sig_inc_threshold=NA, Sig_dec_threshold=NA,
                            Figure=item$file, Note=paste("Error:", e), stringsAsFactors=FALSE)
    empty_peak <- data.frame(Variable=item$display, Exposure=item$var, Increment=item$increment_label,
                             Peak_exposure=NA, Peak_lag_months=NA, Peak_RR=NA, CI_low=NA, CI_high=NA,
                             Peak_RR_95CI=NA, stringsAsFactors=FALSE)
    return(list(cumulative=empty_cum, lag_specific=data.frame(), peak=empty_peak, model=NULL))
  })
  cumulative_results <- rbind(cumulative_results, fit_result$cumulative)
  if(nrow(fit_result$lag_specific)>0) lag_results_all <- rbind(lag_results_all, fit_result$lag_specific)
  peak_results_all <- rbind(peak_results_all, fit_result$peak)
  fit_objects[[item$display]] <- fit_result
}

write.csv(cumulative_results, file.path(work_dir, "DLNM_cumulative_RR_P_95CI_summary.csv"), row.names=FALSE, fileEncoding="GB18030")
if(nrow(lag_results_all)==0) lag_results_all <- data.frame(Variable=character(), Exposure=character(), Increment=character(), Lag=integer(), RR=numeric(), CI_low=numeric(), CI_high=numeric(), RR_95CI=character())
write.csv(lag_results_all, file.path(work_dir, "DLNM_lag_specific_RR_95CI.csv"), row.names=FALSE, fileEncoding="GB18030")
if(nrow(peak_results_all)==0) peak_results_all <- data.frame(Variable=character(), Exposure=character(), Increment=character(), Peak_exposure=numeric(), Peak_lag_months=numeric(), Peak_RR=numeric(), CI_low=numeric(), CI_high=numeric(), Peak_RR_95CI=character())
write.csv(peak_results_all, file.path(work_dir, "DLNM_peak_lag_summary.csv"), row.names=FALSE, fileEncoding="GB18030")

saveRDS(fit_objects, file.path(work_dir, "DLNM_fitted_objects.rds"))

cat("\nAnalysis completed.\n")