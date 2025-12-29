# SARFIMA Model Selection via Machine Learning and Robust Fractional Parameters Estimation

An R function for automated selection of the best SARFIMA model especification.

---

## Description
This function performs automated model selection and estimation for
SARFIMA (Seasonal Autoregressive Fractionally Integrated Moving Average)
models using a machine-learning–based validation scheme.
It evaluates multiple SARFIMA specifications—defined by candidate
non-seasonal and seasonal orders and selects the best model according to
a user-defined accuracy metric (RMSE or MAE). Fractional memory parameters
(both seasonal and non-seasonal) are estimated using the robust
M-periodogram, where robustness is reached via Huber regression.  
When no seasonal structure is detected, only the non-seasonal
fractional parameter is estimated. Cross-validated multi-step forecasts
are used to compute the error metric and guide model selection.

---

## Installation and example of use
Ensure you have the 'mquantile.R' and 'mqper.R' files
```R
rm(list = ls())
graphics.off()
set.seed(123)
# Define the SARFIMA model parameters
sarfima_model <- list(
  phi = 0.3,        # Non-seasonal AR parameter (phi_1)
  dfrac = 0.4,      # Non-seasonal fractional differencing parameter (d)
  seasonal = list(
    phi = 0.6,      # Seasonal AR parameter (Phi_1)
    dfrac = 0.25,    # Seasonal fractional differencing parameter (D)
    period = 12     # Seasonal period (s)
  )
)
# Simulate the time series
ts_data <- arfima::arfima.sim(
  n = 1000,          # Length of the series
  model = sarfima_model,
  sigma2 = 1        # Innovation variance (optional, default is 1)
)
# Plot the simulated data
plot(ts_data, main = "Simulated SARFIMA Process")
# Check the Autocorrelation Function (ACF)
acf(ts_data, lag.max = 100, main = "Simulated SARFIMA Process")
# Load the necessary functions
source('mquantile.R')
source('mqper.R')
source('autosarfima.R')
# Define the train series
TT <- length(ts_data)
train <- TT - 90
ts_train <- ts_data[1:train]
# Estimate the standard model
fit <- arfima::arfima(ts_train, order = c(1,0,1),
                      seasonal = list(
                        order = c(1,0,1),
                        period = 12
                      ))
# Predictions
fc_cla <- predict(fit, n.ahead = 90)
# Plot the results
par(mfrow = c(1,2))
plot.ts(c(ts_train, as.numeric(fc_cla[[1]]$Forecast)), col = "blue", 
        ylab = 'series', main = 'Non cross validation')
lines(as.numeric(ts_data))  
legend("topleft", 
       legend = c('Original series','Predicted series'),
       col = c('black', 'blue'), 
       bty = 'n', lty = c(1,1))
# Estimate the cross-validation model
fit_cv <- fit_sarfima_cv(timeseries = ts_train, nonseas_max = c(2,2), 
                         seas_max = c(2,2), S = 12, k = 30)
# Predictions
fc_cv <- fc_sarfima(obj = fit_cv, h = 90, level = 0, S = 12)
# Plot the results
plot.ts(c(ts_train, fc_cv$mean), main = 'Cros validation', col = 'blue',
        ylab = 'series')
lines(as.numeric(ts_data))
legend("topleft", 
       legend = c('Original series','Predicted series'),
       col = c('black', 'blue'), 
       bty = 'n', lty = c(1,1))