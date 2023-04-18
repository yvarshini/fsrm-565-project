library(fBasics)  # load the packages
library(tseries)
install.packages("astsa")
library(astsa)
library(TSA)
library(fUnitRoots)  # adfTest()
library(fGarch)
require(timeSeries)
library(rugarch)
library(FinTS)
Y_hep
install.packages("rmgarch")
library(rmgarch)
install.packages("depmixS4")
library(depmixS4)
data1 <- Y_hep                     #lodaing in data
hist(data1[,"Goiania"], density = TRUE, main = "Density of Goiania Hepatitis Results")
hist(data1[,"Brasilia"], density = TRUE, main = "Density of Brasilia Hepatitis Results")
## Both locations exhibit right skew


data1 <- ts(data1, frequency = 12, start = 2001, end = 2018)
plot(data1, main = "Time Series Plot of 2 Hepatitis Locations")

##I let the ugarchspec estimate the skew and shape params
##Trying a skewed student t with GARCH(1,1):
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder = c(0,0)), distribution = "sstd")
fit <- ugarchfit(spec, data = data1, solver.control = list(trace=0))
fit


##GARCH(1,0)
spec2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,0)), mean.model = list(armaOrder = c(0,0)), distribution = "sstd")
fit2 <- ugarchfit(spec, data = data1, solver.control = list(trace=0))
fit2


#GARCH(0,1)
spec3 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0,1)), mean.model = list(armaOrder = c(0,0)), distribution = "sstd")
fit3 <- ugarchfit(spec, data = data1, solver.control = list(trace=0))
fit3


##Trying a skewed-ged with GARCH(1,1):
spec4 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder = c(0,0)), distribution = "sged")
fit4 <- ugarchfit(spec2, data = data1, solver.control = list(trace=0))
fit4


#GARCH(1,0)
spec5 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,0)), mean.model = list(armaOrder = c(0,0)), distribution = "sged")
fit5 <- ugarchfit(spec2, data = data1, solver.control = list(trace=0))
fit5

#GARCH(0,1)
spec6 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(0,1)), mean.model = list(armaOrder = c(0,0)), distribution = "sged")
fit6 <- ugarchfit(spec2, data = data1, solver.control = list(trace=0))
fit6


#Getting the ACF and PACF   (MAYBE STANDARDIZE THE RESIDUALS I FORGOT TO DO THAT)
resid1 <- fit@fit$residuals
resid2 <- fit2@fit$residuals
resid3 <- fit3@fit$residuals
resid4 <- fit4@fit$residuals
resid5 <- fit5@fit$residuals
resid6 <- fit6@fit$residuals


acf(resid1, lag.max = 24, main = "ACF for Fit 1")
acf(resid2, lag.max = 24, main = "ACF for Fit 2")
acf(resid3, lag.max = 24, main = "ACF for Fit 3")
acf(resid4, lag.max = 24, main = "ACF for Fit 4")
acf(resid5, lag.max = 24, main = "ACF for Fit 5")
acf(resid6, lag.max = 24, main = "ACF for Fit 6")

pacf(resid1, lag.max = 24, main = "PACF for Fit 1")
pacf(resid2, lag.max = 24, main = "PACF for Fit 2")
pacf(resid3, lag.max = 24, main = "PACF for Fit 3")
pacf(resid4, lag.max = 24, main = "PACF for Fit 4")
pacf(resid5, lag.max = 24, main = "PACF for Fit 5")
pacf(resid6, lag.max = 24, main = "PACF for Fit 6")


##Why are all of these the same :/

## Perform the unit root test(all exact same results)

adf.test(resid1)
adf.test(resid2)
adf.test(resid3)
adf.test(resid4)
adf.test(resid5)
adf.test(resid6)

