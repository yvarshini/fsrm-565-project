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

Y_hep
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


##I let the dccspec due to the correlated and bivariate nature of the data 
##estimate the skew and shape parameters


##Trying a skewed student t with GARCH(1,1), ARIMA model set at (0,1,1) due to unit root presence:
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                   mean.model = list(armaOrder=c(0,1,1)), distribution.model = "sstd")
dccspec1 <- dccspec(uspec = multispec(replicate(2, spec)), dccOrder = c(1,1), distribution = "mvt")

dccfit1 <- dccfit(dccspec, data = data1)
dccfit1

aic <- 16.601
##GARCH(1,0)
spec2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,0)),
                   mean.model = list(armaOrder=c(0,1,1)), distribution.model = "sstd")
dccspec2 <- dccspec(uspec = multispec(replicate(2, spec2)), dccOrder = c(1,1), distribution = "mvt")

dccfit2 <- dccfit(dccspec2, data = data1)
dccfit2
aic <- c(aic, 16.601)
#GARCH(0,1)
spec3 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0,1)),
                    mean.model = list(armaOrder=c(0,1,1)), distribution.model = "sstd")
dccspec3 <- dccspec(uspec = multispec(replicate(2, spec3)), dccOrder = c(1,1), distribution = "mvlaplace")

dccfit3 <- dccfit(dccspec3, data = data1)
dccfit3
aic <- c(aic, 17.022)
##Trying a skewed-ged with GARCH(1,1):
spec4 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                    mean.model = list(armaOrder=c(0,1,1)), distribution.model = "sged")
dccspec4 <- dccspec(uspec = multispec(replicate(2, spec4)), dccOrder = c(1,1), distribution = "mvt")

dccfit4 <- dccfit(dccspec4, data = data1)
dccfit4
aic <- c(aic, 16.503)
#GARCH(2,0)
spec5 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,0)),
                    mean.model = list(armaOrder=c(0,1,1)), distribution.model = "sged")
dccspec5 <- dccspec(uspec = multispec(replicate(2, spec5)), dccOrder = c(1,1), distribution = "mvt")
dccfit5 <- dccfit(dccspec5, data = data1)
dccfit5
aic <- c(aic, 16.494)
#GARCH(0,1)
spec6 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0,1)),
                    mean.model = list(armaOrder=c(0,1,1)), distribution.model = "sged")
dccspec6 <- dccspec(uspec = multispec(replicate(2, spec6)), dccOrder = c(1,1), distribution = "mvt")

dccfit6 <- dccfit(dccspec6, data = data1)
dccfit6
aic <- c(aic, 23.521)

### Returns the same AIC value as GARCH(1,1), we shall disregard models
### 2-7 as they did not perform as well as model 1 or parsimony dictates
### we use simpler models with similar results.


#Getting the ACF and PACF 
resid1 <- dccfit1@mfit$stdresid

adf.test(resid1)
adf.test(std_re)
plot(resid1, main = "Standardized Residual Plot")
TSA::acf(resid1)
pacf(resid1)





resid2 <- dccfit2@mfit$stdresid

resid3 <- dccfit3@mfit$stdresid

resid4 <- dccfit4@mfit$stdresid

resid5 <- dccfit5@mfit$stdresid

resid6 <- dccfit6@mfit$stdresid


acf(resid1, lag.max = 12, main = "ACF for Fit 1")
acf(resid2, lag.max = 12, main = "ACF for Fit 2")
acf(resid3, lag.max = 12, main = "ACF for Fit 3")
acf(resid4, lag.max = 12, main = "ACF for Fit 4")
acf(resid5, lag.max = 12, main = "ACF for Fit 5")
acf(resid6, lag.max = 12, main = "ACF for Fit 6")

pacf(resid1, lag.max = 12)
pacf(resid2, lag.max = 12)
pacf(resid3, lag.max = 12)
pacf(resid4, lag.max = 12)
pacf(resid5, lag.max = 12)
pacf(resid6, lag.max = 12)




