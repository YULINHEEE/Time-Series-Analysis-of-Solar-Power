# upload the following libraries 
library(TSA)
library(fUnitRoots)
library(forecast)
library(tseries)
library(lmtest)
library(zoo)
library(rugarch)
# ================================
# Descriptive Statistic analysis
# ================================

df1<-read.csv("C:/Users/Yulin He/Documents/Datascience/9.Time Series Analysis/Final Project/daily_total_energy.csv")

#summary statistics
summary_stats <- data.frame(
   Min = min(ts_daily),
   Max = max(ts_daily),
   Mean = mean(ts_daily),
   Median = median(ts_daily),
   SD = sd(ts_daily)
)

print(summary_stats)
#check days
length(df1$total_kwh)
#check class
class(df1)
#transform date and keep the right order
df1$date <- as.Date(df1$date)
ts_zoo <- zoo(df1$total_kwh, order.by = df1$date)
plot(ts_zoo, type = "l",
     xlab = "Date", ylab = "Energy (kWh)",
     main = "Daily Solar Energy Production",
     col = "blue")

# Convert data to time series object
ts_daily <- ts(coredata(ts_zoo), start = c(2019, as.numeric(format(start(ts_zoo), "%j"))), frequency = 365)
plot(ts_daily, type = "l",
     xlab = "Date", ylab = "Daily Solar Power Output(kWh)",
     col = "blue")

class(ts_daily)
head(ts_daily)

# Autocorrelation Analysis
# First-order autocorrelation
y<-ts_daily
x<- zlag(ts_daily)
index <- 2:length(x)
cor(y[index],x[index]) #0.7630352
# Lag plot for first-order autocorrelation
plot(y[index], x[index],
     ylab = "Solar Power(kWh)",
     xlab = "Previous Energy",
     main = "Scatter plot of the daily solar energy in consecutive days")

#Normally distribution
qqnorm(ts_daily, main="Q-Q Normal Plot of Daily Solar Energy")
qqline(ts_daily)
#Shapiro-Wilk normality test
shapiro.test(ts_daily)

#===Stationarity Testing===
# ADF test to check whether solar.ts is stationary or not
adf.test(ts_daily) 
pp.test(ts_daily)
#p-value = 0.01, stationarity

# ACF & PACF plot
par(mfrow=c(1,2))
acf(ts_daily,main = "ACF plot",lag.max = 100) #ACF有明显的趋势，所以q=0,意味着白噪点有信息
pacf(ts_daily, main = 'PACF plot',lag.max = 100)
par(mfrow=c(1,1))

##changing variance
r.ts_daily=diff(log(ts_daily))*100
plot(r.ts_daily, ylab='Return ts_daily', xlab='Time', type='o',
     main = "Time series plot of the differenced and transformed daily energy.")

adf.test(r.ts_daily)
pp.test(r.ts_daily)
McLeod.Li.test(y=r.ts_daily,main="McLeod-Li test statistics for daily energy series")

qqnorm(r.ts_daily, main="Q-Q Normal Plot of Daily Solar Energy")
qqline(r.ts_daily)
shapiro.test(r.ts_daily)


#We will first specify the orders of ARMA part. 
#Then we will use the residuals of ARMA part to specify the orders of GARCH part.
#We will use returns series to specify the orders of ARMA part following the model specification steps.
# ================================
# Model Identification - ARMA Part
# ================================
# ACF & PACF plot
par(mfrow=c(1,2))
acf(as.numeric(r.ts_daily),main = "log-differenced transformation ACF plot",lag.max = 20)
pacf(as.numeric(r.ts_daily), main = 'log-differenced transformation PACF plot',lag.max = 20)
par(mfrow=c(1,1))
#{ARMA(1,1), ARMA(2,1), ARMA(3,1), ARMA(4,1)}  can be identified

#EACF
eacf(r.ts_daily, ar.max = 5, ma.max = 5)
#{ARMA(0,2) ARIMA(1,2) ARIMA(2,2)}

# BIC table
par(mfrow=c(1,1))
res = armasubsets(y=r.ts_daily,nar=5,nma=5,y.name='p',ar.method='ols')
plot(res)

#{ARMA(0,1) ARIMA(0,2) ARIMA(2,1) ARIMA(2,2)  ARIMA(2,5)}

#totally: {ARMA(0,1), ARMA(0,2), ARMA(1,1), ARMA(1,2), ARMA(2,1), ARMA(2,2), ARMA(2,5), ARMA(3,1), ARMA(4,1)}

#====residual diagnostics====
residual.analysis <- function(model, std = TRUE, start = 2, 
                              class = c("ARIMA","GARCH","ARMA-GARCH","garch","fGARCH","rugarch"),
                              max.lag = 30, title = NULL) {
   library(TSA)
   class <- match.arg(class)
   
   if (class == "ARIMA") {
      res.model <- if (std) rstandard(model) else residuals(model)
   } else if (class %in% c("GARCH", "garch")) {
      res.model <- model$residuals[start:model$n.used]
   } else if (class == "ARMA-GARCH") {
      res.model <- model@fit$residuals
   } else if (class == "fGARCH") {
      res.model <- model@residuals
   } else if (class == "rugarch") {
      res.model <- residuals(model, standardize = std)
      res.model <- as.numeric(res.model)
   }else {
      stop("The argument 'class' must be one of the predefined types.")
   }
   
   # Ljung-Box test p-values
   lbq.pvals <- sapply(1:max.lag, function(lag) Box.test(res.model, lag = lag, type = "Ljung-Box")$p.value)
   
   par(mfrow = c(2, 3))
   
   # Diagnostic plots
   plot(res.model, type = 'o', ylab = 'Standardised residuals', main = "Time series plot of residuals")
   abline(h = 0)
   hist(res.model, main = "Histogram of residuals", col = "lightgray")
   qqnorm(res.model, main = "QQ plot of residuals")
   qqline(res.model, col = 2)
   acf(as.numeric(res.model), main = "ACF of residuals")
   pacf(as.numeric(res.model), main = "PACF of residuals")
   
   # Ljung-Box p-value plot only
   plot(1:max.lag, lbq.pvals, type = "b", pch = 16, col = "darkgreen",
        main = "Ljung-Box test p-values", xlab = "Lag", ylab = "p-value",
        ylim = c(0, 1))
   abline(h = 0.05, col = "red", lty = 2)
   
   if (!is.null(title)) {
      mtext(title, outer = TRUE, line = -1.5, cex = 1.5, font = 2)
   }
   
   par(mfrow = c(1, 1))
   
   # Shapiro-Wilk test
   cat("\nShapiro-Wilk test:\n")
   print(shapiro.test(res.model))
}


#====build a model fitting function====
fit_arima_models <- function(ts_data, orders, std = TRUE, start = 2) {
   
   #using ML and CSS methods to fit ARIMA model
   #print the coeftest results
   results <- list()
   
   for (order in orders) {
      order_str <- paste(order, collapse = "")
      
      cat("\n===== ARIMA(", paste(order, collapse = ","), ") with ML =====\n")
      model_ml <- arima(ts_data, order = order, method = "ML", include.mean = TRUE)
      print(coeftest(model_ml))
      results[[paste0("model", order_str, "_ML")]] <- model_ml
      
      cat("\n===== ARIMA(", paste(order, collapse = ","), ") with CSS =====\n")
      model_css <- arima(ts_data, order = order, method = "CSS", include.mean = TRUE)
      print(coeftest(model_css))
      results[[paste0("model", order_str, "_CSS")]] <- model_css
   
      cat("Residual diagnostics for ARIMA(", paste(order, collapse = ","), ") CSS:\n")
      residual.analysis(model = model_css, std = std, start = start, class = "ARIMA",
                        title = paste0("Residual Diagnostics for ARMA(", order[1], ",", order[3], ")"))
   }
   
   return(results)
}


orders_to_test <- list(c(0,0,1),c(0,0,2),c(1,0,1),c(1,0,2),c(2,0,1),c(2,0,2),c(2,0,5),c(3,0,1),c(4,0,1))

results <- fit_arima_models(r.ts_daily, orders_to_test)

# ================================
# Comparised AIC and BIC
# ================================
#AIC/BIC SCORE
extract_and_sort_ic <- function(model_list, score = c("aic", "bic")) {
   score <- match.arg(score)
   
   ic_values <- sapply(names(model_list), function(name) {
      model <- model_list[[name]]
      val <- tryCatch(
         if (score == "aic") AIC(model) else BIC(model),
         error = function(e) NA 
      )
      return(val)
   })
   # delete NA models
   ic_values <- ic_values[!is.na(ic_values)]
   # order the results
   df <- data.frame(Model = names(ic_values), IC = ic_values)
   df <- df[order(df$IC), ]
   rownames(df) <- NULL
   return(df)
}

extract_and_sort_ic(results[grep("_ML$", names(results))], score = "aic")
extract_and_sort_ic(results[grep("_ML$", names(results))], score = "bic")

#model_102 is the best
# ================================
# Accuracy Evaluation
# ================================
model_001_css = Arima(r.ts_daily,order=c(0,0,1), method='CSS')
model_002_css = Arima(r.ts_daily,order=c(0,0,2), method='CSS')
model_101_css = Arima(r.ts_daily,order=c(1,0,1), method='CSS')
model_102_css = Arima(r.ts_daily,order=c(1,0,2), method='CSS')
model_201_css = Arima(r.ts_daily,order=c(2,0,1), method='CSS')
model_202_css = Arima(r.ts_daily,order=c(2,0,2), method='CSS')
model_205_css = Arima(r.ts_daily,order=c(2,0,5), method='CSS')
model_301_css = Arima(r.ts_daily,order=c(3,0,1), method='CSS')
model_401_css = Arima(r.ts_daily,order=c(4,0,1), method='CSS')
# get accuracy results
Smodel_001_css <- accuracy(model_001_css)[1:7]
Smodel_002_css <- accuracy(model_002_css)[1:7]
Smodel_101_css <- accuracy(model_101_css)[1:7]
Smodel_102_css <- accuracy(model_102_css)[1:7]
Smodel_201_css <- accuracy(model_201_css)[1:7]
Smodel_202_css <- accuracy(model_202_css)[1:7]
Smodel_205_css <- accuracy(model_205_css)[1:7]
Smodel_301_css <- accuracy(model_301_css)[1:7]
Smodel_401_css <- accuracy(model_401_css)[1:7]
# cpllect the results in one data frame
df.Smodels <- data.frame(
   rbind(Smodel_001_css, Smodel_002_css, Smodel_101_css, Smodel_102_css,
         Smodel_201_css, Smodel_202_css, Smodel_205_css, Smodel_301_css, Smodel_401_css))

colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,0,1)", "ARIMA(0,0,2)", "ARIMA(1,0,1)","ARIMA(1,0,2)",
                          "ARIMA(2,0,1)", "ARIMA(2,0,2)", "ARIMA(2,0,5)", "ARIMA(3,0,1)", "ARIMA(4,0,1)")
round(df.Smodels,  digits = 3)

#ARMA(1,2) fitted to the energy return series is the best model

#Overfitting models are ARMA(2,2) and ARMA(1,3)
model_103 = arima(r.ts_daily,order=c(1,0,3), method='ML')
coeftest(model_103)

model_103_css = arima(r.ts_daily,order=c(1,0,3), method='CSS')
coeftest(model_103_css)
residual.analysis(model = model_103_css)

orders_to_test <- list(c(0,0,1),c(0,0,2),c(1,0,1),c(1,0,2),c(1,0,3),c(2,0,1),c(2,0,2),c(2,0,5),c(3,0,1),c(4,0,1))
results <- fit_arima_models(r.ts_daily, orders_to_test)
extract_and_sort_ic(results[grep("_ML$", names(results))], score = "aic")
extract_and_sort_ic(results[grep("_ML$", names(results))], score = "bic")


model_103_css = Arima(r.ts_daily,order=c(1,0,3), method='CSS')
Smodel_103_css <- accuracy(model_103_css)[1:7]
# cllect the results in one data frame
df.Smodels <- data.frame(
   rbind(Smodel_001_css, Smodel_002_css, Smodel_101_css, Smodel_102_css, Smodel_103_css,
         Smodel_201_css, Smodel_202_css, Smodel_205_css, Smodel_301_css, Smodel_401_css))

colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1")
rownames(df.Smodels) <- c("ARIMA(0,0,1)", "ARIMA(0,0,2)","ARIMA(1,0,1)","ARIMA(1,0,2)","ARIMA(1,0,3)","ARIMA(2,0,1)","ARIMA(2,0,2)","ARIMA(2,0,5)","ARIMA(3,0,1)","ARIMA(4,0,1)")
round(df.Smodels,  digits = 3)

#==================================
#Model Identification - GARCH Part
#====================================
#==== use the residuals of ARMA(1,2) model to identify the orders of GARCH.====
model_102 <- arima(r.ts_daily, order = c(1,0,2), method='ML')
m12residuals = model_102$residuals
abs.res = abs(m12residuals)
sq.res = m12residuals^2
par(mfrow=c(1,2))
acf(as.numeric(abs.res), ci.type = "ma", main="ACF plot for absolute return series")
pacf(as.numeric(abs.res), main="PACF plot for absolute return series ")
par(mfrow=c(1,1))
eacf(abs.res)
#GARCH(0,1) GARCH(1,1) GARCH(2,1) GARCH(0,2) , GARCH(0,2) , GARCH(2,2)
par(mfrow=c(1,2))
acf(as.numeric(sq.res), ci.type = "ma",main = "ACF plot for square return series")
pacf(as.numeric(sq.res), main = "PACF plot for square return series")
par(mfrow=c(1,1))
eacf(sq.res)

#Overall:{GARCH(0,1), GARCH(0,2), GARCH(1,1), GARCH(1,2), GARCH(2,0), GARCH(2,1),GARCH(2,2)}

#====fit these models using ugarchfit() function from rugarch package.====
#====fit ARIMA(1,0,2)+GARCH(p,q)====
fit_arma12_garch_models_minimal <- function(data,
                                            garch_orders = list(c(0,1), c(0,2), c(1,1), c(1,2), c(2,0), c(2,1), c(2,2)),
                                            plot_path = "garch_diagnostics") {
   dir.create(plot_path, showWarnings = FALSE)
   
   model_names <- c()
   aic_vec <- c()
   bic_vec <- c()
   
   for (order in garch_orders) {
      p <- order[1]
      q <- order[2]
      model_name <- paste0("ARMA(1,2)+GARCH(", p, ",", q, ")")
      cat("\n===========================\nFitting:", model_name, "\n")
      
      spec <- ugarchspec(
         variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
         mean.model = list(armaOrder = c(1, 2), include.mean = TRUE),
         distribution.model = "norm"
      )
      
      fit <- tryCatch({
         ugarchfit(spec, data = data)
      }, error = function(e) {
         warning(paste("Failed to fit", model_name, ":", e$message))
         return(NULL)
      })
      
      if (!is.null(fit)) {
         show(fit)
         

         cat("Residual diagnostics for", model_name, ":\n")
         par(mfrow = c(3, 2), oma = c(0, 0, 2, 0))
         residual.analysis(fit, class = "rugarch")
         mtext(model_name, outer = TRUE, cex = 1.3, font = 2)
         par(mfrow = c(1, 1))  # reset
         
         #AIC / BIC
         ic <- infocriteria(fit)
         aic_vec <- c(aic_vec, ic[1])
         bic_vec <- c(bic_vec, ic[2])
         model_names <- c(model_names, model_name)
      }
   }
   
   result_df <- data.frame(Model = model_names, AIC = aic_vec, BIC = bic_vec)
   result_df <- result_df[order(result_df$AIC), ]
   rownames(result_df) <- NULL
   return(result_df)
}


result_table <- fit_arma12_garch_models_minimal(r.ts_daily)
print(result_table)

# ================================
#ARMA(1,2)+GARCH(2,0) forecasting.
#==================================
library(xts)
r.xts <- xts(r.ts_daily, order.by = as.Date(df1$date[-1]))
spec <- ugarchspec(
   variance.model = list(model = "sGARCH",garchOrder = c(2, 0)),
   mean.model = list(armaOrder = c(1, 2), include.mean = TRUE),
   distribution.model = "norm")

model_12_20_fit <- ugarchfit(spec = spec, data = r.xts,
                            solver = "hybrid",
                            solver.control = list(trace=0))

par(mfrow=c(1,1))
plot(model_12_20_fit, which = 1)

# 10 day ahead forecast
frc <- ugarchforecast(model_12_20_fit,n.ahead=10,data=r.ts_daily)
frc
class(frc)
plot(frc, which = 1) ## Forecasted return series + conditional SD bands
plot(frc, which = 3) ## Forecasted conditional sigma (volatility)

#====better visualization====
plot_garch_forecast_with_ci <- function(
      return_past,
      forecast_object,
      date_vector,
      save_path = NULL
) {
   # date
   dates_past <- as.Date(date_vector)[-1]
   n_forecast <- length(forecast_object@forecast$seriesFor)
   dates_future <- seq(from = tail(dates_past, 1) + 1, by = "day", length.out = n_forecast)
   dates_all <- c(dates_past, dates_future)
   
   # data
   returns_full <- c(return_past, rep(NA, n_forecast))
   forecast_returns <- c(rep(NA, length(return_past)), forecast_object@forecast$seriesFor)
   sigma_forecast <- forecast_object@forecast$sigmaFor
   upper <- forecast_object@forecast$seriesFor + 2 * sigma_forecast
   lower <- forecast_object@forecast$seriesFor - 2 * sigma_forecast
   
   upper_full <- c(rep(NA, length(return_past)), upper)
   lower_full <- c(rep(NA, length(return_past)), lower)
   
   par(bg = "white", mar = c(4, 4, 3, 2))
   plot(dates_all, returns_full, type = "l", col = "blue", lwd = 1,
        ylab = "Return (%)", xlab = "Date", 
        main = "Actual + Forecasted Return with ±2σ Band",
        ylim = range(c(returns_full, upper, lower), na.rm = TRUE))
   
   polygon(
      x = c(dates_future, rev(dates_future)),
      y = c(lower, rev(upper)),
      col = rgb(1, 0, 0, 0.2), border = NA)
   lines(dates_all, forecast_returns, col = "red", lwd = 2)
   legend("topleft", legend = c("Actual", "Forecast", "±2σ Band"),
          col = c("blue", "red", rgb(1, 0, 0, 0.2)), lty = c(1, 1, NA), 
          lwd = c(2, 2, NA), pch = c(NA, NA, 15), pt.cex = 2, bty = "n")
   
   grid()
   
   if (!is.null(save_path)) dev.off()
}


plot_garch_forecast_with_ci(
   return_past = r.ts_daily,
   forecast_object = frc,
   date_vector = df1$date
)


# All the transformations:To get the actual forecasts, we need to take the tansformation and differencing back.
#r.ts_daily = diff(log(ts_daily))*100
last_value <- tail(ts_daily, 1)  #last observation
r_forecast <- frc@forecast$seriesFor  # return
r_forecast

n <- length(r_forecast)
x_pred <- numeric(n)
x_pred[1] <- last_value * exp(r_forecast[1] / 100)

for (i in 2:n) {
   x_pred[i] <- x_pred[i - 1] * exp(r_forecast[i] / 100)
}
print(x_pred)


# x-index
dates <- as.Date(df1$date)
last_date <- tail(dates, 1)
future_dates <- seq(from = last_date + 1, by = "day", length.out = n)
time_index <- c(dates, future_dates)

# original data + predicted data
ts_full <- c(ts_daily, rep(NA, n))
pred_full <- c(rep(NA, length(ts_daily)), x_pred)

par(bg = "white", mar = c(4, 4, 3, 2), cex.axis = 1.2, cex.lab = 1.3, cex.main = 1.4)
plot(time_index, ts_full, type = "l", col = "darkgrey", lwd = 2,
     ylab = "Electricity (MWh)", xlab = "Date",
     main = "Solar Power Output Forecast Using ARMA(1,2)-GARCH(2,0)",
     ylim = range(c(ts_full, pred_full), na.rm = TRUE))
lines(time_index, pred_full, col = "blue", lwd = 2)
legend("topleft", legend = c("Actual", "Forecast"),
       col = c("darkgrey", "blue"), lty = 1, lwd = 2, bty = "n")
grid(col = "lightgray", lty = "dotted")
box()





