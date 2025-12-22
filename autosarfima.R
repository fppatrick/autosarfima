#' @title SARFIMA Model Selection via Machine Learning and Robust Fractional Parameter Estimation
#'
#' @description
#' This function performs automated model selection and estimation for
#' SARFIMA (Seasonal Autoregressive Fractionally Integrated Moving Average)
#' models using a machine-learning–based validation scheme.
#'
#' It evaluates multiple SARFIMA specifications—defined by candidate
#' non-seasonal and seasonal orders and selects the best model according to
#' a user-defined accuracy metric (RMSE or MAE). Fractional memory parameters
#' (both seasonal and non-seasonal) are estimated using the robust
#' M-periodogram, where robustness is reached via Huber regression.  
#'
#' When no seasonal structure is detected, only the non-seasonal
#' fractional parameter is estimated. Cross-validated multi-step forecasts
#' are used to compute the error metric and guide model selection.
#'
#' @author Patrick Ferreira Patrocinio
#'
#' @param timeseries A univariate time series object.
#' @param bandw.exp Coefficient bandwich that will truncate the spectrum length.
#'   Default: `0.78`
#' @param candidates Integer vector of possible seasonal periods.  
#'   Default: `c(0, 7, 12)`.
#' @param include.mean Logical indicating whether the mean should be included
#'   in the SARFIMA specification. Default: `FALSE`.
#' @param nonseas_max Integer vector specifying the maximum non-seasonal
#'   AR and MA orders to evaluate, e.g., `c(p_max, q_max)`.  
#'   Default: `c(3, 3)`.
#' @param seas_max Integer vector specifying the maximum seasonal SAR and SMA
#'   orders to evaluate, e.g., `c(P_max, Q_max)`.  
#'   Default: `c(3, 3)`.
#' @param metric Character string indicating the accuracy metric to optimise.
#'   Options: `"RMSE"` or `"MAE"`.  
#'   Default: `"RMSE"`.
#' @param rm_out Logical indicating whether outliers should be removed from the
#'   input time series prior to model fitting.  
#'   Default: `FALSE`.
#' @param S Integer defining the known seasonal period. If `NULL` or `FALSE`,
#'   the function estimates periodicity using the robust M-periodogram.  
#'   Default: `7`.
#' @param h Integer specifying the forecast horizon used for cross-validated
#'   error computation.  
#'   Default: `1`.
#' @param k Integer specifying the number of folds used in the rolling
#'   cross-validation scheme.  
#'   Default: `45`.
#' @param method Character string defining the estimation method for
#'   SARFIMA parameters.  
#'   Default: `"CSS"`.
#'
#' @return A list containing the selected SARFIMA model, fractional parameter
#'   estimates, validation metrics, and forecasting outputs.
#'
#' @seealso \code{fracdiff}, \code{rlm}, \code{MperioReg}
#-------------------------------
# Necessary packages
#-------------------------------
library(MASS)
library(tidyverse)
library(lubridate)
library(forecast)
library(dplyr)
library(future.apply)
#-------------------------------
# Function to replace outliers by the median (take care about it!).
#-------------------------------
m_data <- function(timeseries){
  out <- boxplot.stats(timeseries, coef = 3)$out
  timeseries[timeseries %in% out] <- median(timeseries)
  return(timeseries)
}
#-------------------------------
# Function to pass the filters non-seasonal and seasonal.
#-------------------------------
frac_filter <- function(timeseries, d) {
  
  stopifnot((TT <- length(timeseries)) >= 2)
  #x <- x - mean(x)
  np2 <- nextn(TT+TT - 1L)
  
  pad <- rep.int(0, np2-TT)
  
  k <- seq_len(TT - 1L)
  
  b <- c(1, cumprod((k - (d+1))/ k), pad)
  ## ~= convolve(x, b, type = "filter") :
  
  dx <- fft(fft(b) * fft(c(timeseries, pad)), inverse =TRUE)[seq_len(TT)] / np2
  
  Re(dx)
}
#-------------------------------
# Function to differentiate the series using fractional parameters
#-------------------------------
frac_diff <- function(timeseries,d,D){
  
  y <- frac_filter(timeseries = timeseries, d = d)
  
  u <- frac_filter(timeseries = y, d = D)
  
  return(u)
}
#-------------------------------
# This function estimate the fractional long-memory parameters (non-seasonal and, 
# if applicable, seasonal).
#-------------------------------
RobustdSperio <- function(timeseries, bandw.exp = 0.78, tau = 0.5) {
  
  if (!is.numeric(timeseries)) timeseries <- as.numeric(timeseries)
  timeseries <- na.fail(as.ts(timeseries))
  
  if (!is.null(dim(timeseries))) stop("Only univariate time series are allowed")
  
  if (any(is.na(timeseries))) stop("There are NA values")
  
  n <- length(timeseries)
  g <- trunc(n^bandw.exp)
  j <- 1:g
  w <- 2 * pi * j / n
  
  timeseries <- timeseries - mean(timeseries, na.rm = TRUE)
  
  per <- mqper(timeseries = timeseries, tau = tau)$perior
  
  if (length(per) < g) stop("Series to short for the choosen bandw.exp")
  
  periodogram <- per[1:g]
  
  s <- Period(timeseries = timeseries)
  
  if (is.null(s) || s < 0) s <- 0
  
  y.reg <- log(periodogram)
  
  d.reg <- log((2 * sin(w / 2))^2)
  
  if (s != 0) {
    D.reg <- log((2 * sin(s * w / 2))^2)
    X <- cbind(1, -D.reg, -d.reg)
    colnames(X) <- c("Intercept", "D.reg", "d.reg")
  } else {
    X <- cbind(1, -d.reg)
    colnames(X) <- c("Intercept", "d.reg")
  }
  
  fit <- rlm(X, y.reg, psi = psi.huber)
  
  coefs <- coef(fit)
  
  d <- coefs["d.reg"]
  
  D <- if ("D.reg" %in% names(coefs)) coefs["D.reg"] else 0
  
  return(list(d = as.numeric(d), D = as.numeric(D)))
}
#-------------------------------
# Function to calculate the periodicity based on the Robust periodogram. 
# If you are sure about the time series features set the value in the
# fit_sarfima_cv function (parameter S).
#-------------------------------
Period <- function(timeseries, candidates = c(0,7,12), tau= 0.5) {
  
  modts = timeseries - mean(timeseries)
  
  robper <- mqper(timeseries = timeseries, tau = tau)
  
  freqs <- robper$freq
  
  spec_den <- robper$perior
  
  max_index <- order(spec_den, decreasing = TRUE)
  dom_freq <- freqs[max_index]
  s <- round(1 / dom_freq)
  s <- s[which(s!=0)]
  
  ntests <- length(candidates)
  
  tests <- sapply(1:ntests, function(i){which(s%%candidates[i]==0)})
  
  S <- NULL
  
  for (i in 1:ntests) {
    S[i] <- length(tests[[i]])
  }
  
  sea_per <- ifelse(max(S)==0, 0, candidates[which.max(S)])
  
  return(sea_per)
}
#-------------------------------
# Rolling cross-validation of the SARFIMA model. If you do not know 
# what you are doing would be better you use auto.arima function from 
# forecast package. If you do not have at least a basic backgorund of time series... 
# Please, stay away for now and call someone who has.
#-------------------------------
fit_sarfima_cv <- function(timeseries,
                           include.mean = FALSE,
                           nonseas_max = c(4,4),
                           seas_max = c(4,4),
                           metric = "RMSE",
                           rm_out = FALSE,
                           S = 7,
                           h = 1,
                           k = 35,
                           method = "CSS") {
  
  start_time <- Sys.time()
  
  if (isTRUE(rm_out)) {
    clean_series <- m_data(timeseries = timeseries)
  } else {
    clean_series <- timeseries
  }
  
  if (!(metric %in% c("MAE","RMSE")))
    stop("Metrics allowed: MAE or RMSE")
  
  TT <- length(clean_series)
  
  if (is.null(h) || h==FALSE) h <- round(0.1 * TT)
  
  if (is.null(S) || S==FALSE) {
    S <- Period(timeseries = clean_series)
  }
  
  if (!(S %in% c(0, 7, 12)))
    warning("Please confirm the seasonal period of the series.")
  
  diffs <- ndiffs(clean_series)
  
  sym <- lawstat::symmetry.test(clean_series)
  
  if (sym$p.value <= 0.1) {
    diff_series <- diff(clean_series)
  } else {
    diff_series <- clean_series
  }
  
  # if(!(S==0)) {
  #   sdiff_series <- diff(diff_series, S)
  # } else {
  #   sdiff_series <- diff_series
  # }
  
  
  d_frac <- RobustdSperio(timeseries = diff_series)
  
  d_vals <- unique(c(0:diffs))
  
  D_vals <- unique(c(0:1))
  
  frac_series <- frac_diff(timeseries = timeseries, d = d_frac$d, D = d_frac$D)
  
  par_grid <- expand_grid(p=0:nonseas_max[1], d=d_vals, q=0:nonseas_max[2],
                          P=0:seas_max[1], D=D_vals, Q=0:seas_max[2])
  
  plan(multisession, workers = max(1, parallel::detectCores() - 1))
  
  results <- future_lapply(1:nrow(par_grid), function(i) {
    nonseas <- as.numeric(par_grid[i, 1:3])
    seas    <- as.numeric(par_grid[i, 4:6])
    
    fold_size <- floor((TT - h) / k)
    folds <- seq(from = fold_size, to = TT - h, by = fold_size)
    
    if (length(folds) > k) {
      folds <- folds[1:k]
    }
    
    cv_errors <- rep(NA, length(folds))
    
    for (j in seq_along(folds)) {
      train_end <- folds[j]
      train_data <- frac_series[1:train_end]
      test_data <- frac_series[(train_end + 1):min(train_end + h, TT)]
      
      if (length(test_data) < 1) {
        cv_errors[j] <- NA
        next
      }
      
      fit <- tryCatch({
        Arima(train_data, 
              order = nonseas, 
              seasonal = list(order = seas, period = S),
              include.mean = include.mean, 
              method = method)
      }, error = function(e) {
        NULL
      })
      
      if (is.null(fit)) {
        cv_errors[j] <- NA
        next
      }
      
      fc <- tryCatch({
        forecast::forecast(fit, h = length(test_data))
      }, error = function(e) {
        NULL
      })
      
      if (is.null(fc)) {
        cv_errors[j] <- NA
        next
      }
      
      pred <- fc$mean
      
      min_len <- min(length(pred), length(test_data))
      if (min_len < 1) {
        cv_errors[j] <- NA
        next
      }
      
      pred <- pred[1:min_len]
      actual <- test_data[1:min_len]
      
      fold_error <- switch(metric,
                           RMSE = sqrt(mean((actual - pred)^2, na.rm = TRUE)),
                           MAE  = mean(abs(actual - pred), na.rm = TRUE))
      
      cv_errors[j] <- fold_error
    }
    
    valid_errors <- cv_errors[!is.na(cv_errors)]
    if (length(valid_errors) == 0) {
      return(NULL)
    }
    
    err <- mean(valid_errors, na.rm = TRUE)
    
    list(params = c(nonseas, seas), error = err, fold_errors = cv_errors)
    
  })
  
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0)
    stop("All models failed to converge.")
  
  errors <- sapply(results, `[[`, "error")
  
  if (all(is.infinite(errors))) {
    stop("All models produced infinite errors. Check parameter ranges.")
  }
  
  best <- which.min(errors)
  
  best_params <- results[[best]]$params
  
  best_model <- tryCatch({
    Arima(frac_series,
          order = best_params[1:3],
          seasonal = list(order = best_params[4:6], period = S),
          include.mean = include.mean,
          method = method)
  }, error = function(e) NULL)
  
  if (is.null(best_model)) {
    stop("Best model failed to fit on full dataset.")
  }
  
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  cat("Time elapsed: ", round(time_elapsed, 2), " ", attr(time_elapsed, "units"), "\n")
  
  list(
    model = best_model,
    params = best_params,
    d = d_frac$d,
    D = d_frac$D
    #cv_error = errors[best],
    #fc = fc,
    #all_errors = errors
  )
}

fc_sarfima <- function(obj, h, level){
  
  frac_fc <- forecast::forecast(obj = obj[["model"]], h = h, level = level)
  
  Mean <- c(obj[["model"]]$x, frac_fc$mean)
  
  Lower <- c(obj[["model"]]$x, frac_fc$lower)
  
  Upper <- c(obj[["model"]]$x, frac_fc$upper)
  
  df <- data.frame(Mean, Lower, Upper)
  
  TT <- nrow(df)
  train <- TT - h
  
  df_orig <- sapply(df, frac_diff, d = -obj[["d"]], D = -obj[["D"]])
  
  df_orig <- as.data.frame(df_orig[(train+1):TT,])
  
  fc <- data.frame(mean = df_orig$Mean, 
                   lower = df_orig$Lower, 
                   upper = df_orig$Upper)
  
  return(fc)
  
}