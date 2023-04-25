library(fBasics)  # load the packages
library(tseries)
library(astsa)
library(TSA)
library(fUnitRoots)  # adfTest()
library(fGarch)
require(timeSeries)
library(rugarch)
library(FinTS)
library(forecast)
library(rmgarch)
library(depmixS4)


data1 <- Y_hep                     #lodaing in data
hist(data1[,"Goiania"], density = TRUE, main = "Density of Goiania Hepatitis Results")
hist(data1[,"Brasilia"], density = TRUE, main = "Density of Brasilia Hepatitis Results")
## Both locations exhibit right skew


data <- ts(data1, frequency = 12, start = 2001, end = 2018)
data_b <- cbind(data1[,"Goiania"],data1[,"Brasilia"])
fit_arima_x <- auto.arima(data[, "Goiania"], ic = "aic")
fit_arima_y <- auto.arima(data[, "Brasilia"], ic = "aic")


plot(data, main = "Time Series Plot of 2 Hepatitis Locations")
TSA::acf(data[,"Goiania"], lag.max = 12)
data <- diff(data, lag = 1)
data_c <- cbind(data[,"Goiania"],data[,"Brasilia"])

pacf(data[,"Goiania"])
TSA::acf(data[,"Brasilia"])
pacf(data[,"Brasilia"])


##I let the ugarchspec estimate the skew and shape params
##Trying a skewed student t with GARCH(1,1):
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                   mean.model = list(armaOrder = c(fit_arima_x$arma[1],
                  fit_arima_y$arma[1])), distribution.model = "sstd")
fit <- ugarchfit(spec, data = data_c, solver.control = list(trace=0))
fit


##GARCH(1,0)
spec2 <- ugarchspec(variance.model = list(model = "sGARCH",
        garchOrder = c(1,0)),mean.model = list(armaOrder = 
        c(fit_arima_x$arma[1], fit_arima_y$arma[1])), distribution.model = "sstd")
fit2 <- ugarchfit(spec2, data = data_b, solver.control = list(trace=0))
fit2


#GARCH(0,1)
spec3 <- ugarchspec(variance.model = list(model = "sGARCH",
      garchOrder = c(0,1)),mean.model = list(armaOrder = 
      c(fit_arima_x$arma[1], fit_arima_y$arma[1])), distribution.model = "sstd")
fit3 <- ugarchfit(spec3, data = data_b, solver.control = list(trace=0))

fit3
##Trying a skewed-ged with GARCH(1,1):
spec4 <- ugarchspec(variance.model = list(model = "sGARCH",
       garchOrder = c(1,1)),mean.model = list(armaOrder = 
       c(fit_arima_x$arma[1], fit_arima_y$arma[1])), distribution.model = "sged")
fit4 <- ugarchfit(spec4, data = data_b, solver.control = list(trace=0))
fit4

#GARCH(1,0)
spec5 <- ugarchspec(variance.model = list(model = "sGARCH",
      garchOrder = c(1,0)),mean.model = list(armaOrder = 
      c(fit_arima_x$arma[1], fit_arima_y$arma[1])), distribution.model = "sged")
fit5 <- ugarchfit(spec5, data = data_b, solver.control = list(trace=0))
fit5

#GARCH(0,1)
spec6 <- ugarchspec(variance.model = list(model = "sGARCH",
      garchOrder = c(0,1)),mean.model = list(armaOrder = 
      c(fit_arima_x$arma[1], fit_arima_y$arma[1])), distribution.model = "sged")
fit6 <- ugarchfit(spec6, data = data_b, solver.control = list(trace=0))
fit6

### Trying higher order GARCH models under skewed t-dist
spec7 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,2)),
                   mean.model = list(armaOrder = c(fit_arima_x$arma[1],
                    fit_arima_y$arma[1])), distribution.model = "sstd")
fit7 <- ugarchfit(spec, data = data_b, solver.control = list(trace=0))
fit7

### Returns the same AIC value as GARCH(1,1), we shall disregard models
### 2-7 as they did not perform as well as model 1 or parsimony dictates
### we use simpler models with similar results.


#Getting the ACF and PACF 
resid1 <- fit@fit$residuals
se1 <- fit@fit$sigma
std_re <- resid1/se1
adf.test(resid1)
adf.test(std_re)
plot(std_re, main = "Standardized Residual Plot")
TSA::acf(std_re)
pacf(std_re)

### There appears to be no unit root present in the data! Also the 
### standardized residuals show no obvious heteroskedasicity! 


### Keeping the below code for now because it may serve a purpose. 
### Using fit one since it had the lowest AIC and the difference series 
### looks pretty good relative to non-differenced!

resid2 <- fit2@fit$residuals
se2 <- fit2@fit$sigma
std_re2 <- resid2/se2
resid3 <- fit3@fit$residuals
se3 <- fit3@fit$sigma
std_re3 <- resid3/se3
resid4 <- fit4@fit$residuals
se4 <- fit4@fit$sigma
std_re4 <- resid4/se4
resid5 <- fit5@fit$residuals
se5 <- fit5@fit$sigma
std_re5 <- resid5/se5
resid6 <- fit6@fit$residuals
se6 <- fit6@fit$sigma
std_re6 <- resid6/se6

acf(resid1, lag.max = 12, main = "ACF for Fit 1")
acf(std_re, lag.max = 12, main = "ACF for stdzd Fit 1")
acf(resid2, lag.max = 12, main = "ACF for Fit 2")
acf(resid3, lag.max = 12, main = "ACF for Fit 3")
acf(resid4, lag.max = 12, main = "ACF for Fit 4")
acf(resid5, lag.max = 12, main = "ACF for Fit 5")
acf(resid6, lag.max = 12, main = "ACF for Fit 6")

pacf(resid1, lag.max = 12, main = "PACF for Fit 1")
pacf(resid2, lag.max = 12, main = "PACF for Fit 2")
pacf(resid3, lag.max = 12, main = "PACF for Fit 3")
pacf(resid4, lag.max = 12, main = "PACF for Fit 4")
pacf(resid5, lag.max = 12, main = "PACF for Fit 5")
pacf(resid6, lag.max = 12, main = "PACF for Fit 6")



adf.test(resid1)
adf.test(resid2)
adf.test(resid3)
adf.test(resid4)
adf.test(resid5)
adf.test(resid6)

