# Solar Power Time Series Forecasting with ARIMA-GARCH

This project focuses on short-term forecasting and volatility modeling of daily solar power generation using ARIMA and GARCH models. It combines classical time series forecasting with financial-style risk modeling to support better energy planning and uncertainty management.

---

## Project Overview

- **Objective**:  
  To forecast the next 15 days of solar energy output and evaluate the dynamics of return volatility, helping energy providers better plan capacity and manage unpredictability.

- **Data**:  
  Raw hourly solar generation data was aggregated into daily totals (in kWh). A log return transformation was applied to stabilize variance and enable GARCH modeling.

---

## Methodology

1. **Exploratory Data Analysis (EDA)**  
   - Aggregated raw data to daily kWh.
   - Visualized trends, seasonality, and variance structure.

2. **Stationarity Testing**  
   - ADF and PP tests confirmed the stationarity of the return series after log transformation.

3. **ARIMA Model Fitting**  
   - Identified optimal (p,d,q) using AIC and BIC.
   - Performed residual diagnostics (Ljung-Box test, ACF/PACF).

4. **GARCH Modeling**  
   - Tested multiple GARCH(p,q) specifications on residuals.
   - Selected best-fit model based on AIC and volatility patterns.

5. **Forecasting**  
   - Forecasted log returns with the GARCH model and back-transformed to kWh.
   - Plotted prediction intervals to visualize uncertainty.

## Key Findings

- **Best ARIMA Model**: ARIMA(1,0,2) based on goodness-of-fit and residual checks.
- **Best GARCH Model**: GARCH(1,1), capturing volatility clustering in solar returns.
- **Uncertainty**: Forecast intervals widened over time, highlighting increasing uncertainty in future solar output.
- **Application**: Energy providers can use this model to anticipate risky or unstable days for power planning.

## Tools Used

- **Language**: R  
- **Packages**: `forecast`, `rugarch`, `tseries`, `TSA`
