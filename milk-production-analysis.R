library(astsa)
library(forecast)
library(ggplot2)
#library(MASS)
library(tseries)
library(xts)
#library(zoo)
library(LSTS)

# Importing the data

production.df <- read.csv('./milk_production_data.csv')
milk.ts <- ts(production.df$Milk.Prod, start=1995, freq=12)
View(milk.ts)
head(milk.ts, n=12)
summary(milk.ts)
autoplot(milk.ts, main='Milk Production Time Series (pounds per cow)', ylab='Milk Production (pounds per cow)')

autoplot(milk.ts, main='Milk Production Time Series (pounds per cow)') +
  geom_smooth() +
  labs("Milk production (in pounds per cow)",
       y = "Milk production (in pounds per cow)",
       x = NULL)

ggseasonplot(milk.ts)
ggseasonplot(milk.ts, polar=TRUE)
ggmonthplot(milk.ts, main="dwef")

acf2(milk.ts, main=paste("ACF & PACF for Series: Milk Production (pounds per cow)"))
gglagplot(milk.ts)
Box.test(milk.ts, lag = 4, fitdf = 0, type = "Lj")

# 1. Removing heteroscedasticity

milkLog <- log(milk.ts)
autoplot(milk.ts, main='Milk Production Time Series (pounds per cow)', ylab='Milk Production')
autoplot(milkLog, main='Log Milk Production Time Series (pounds per cow)', ylab='Log of Milk Production')
par(mfrow=c(2,1))
tsplot(milk.ts, main='Milk Production Time Series (pounds per cow)', ylab='Milk Production')
tsplot(milkLog, main='Log Milk Production Time Series (pounds per cow)', ylab='Log of Milk Production')
par(mfrow=c(1,1))
acf2(milkLog, main='ACF & PACF for Logged Series: Milk Production')

print(lambda <- BoxCox.lambda(milk.ts))  # find the optimal value for parameter lambda
plot.ts(BoxCox(milk.ts, lambda = lambda))

autoplot(milkLog, main='Milk Production Time Series transformed using log', ylab='Log of Milk Production')
autoplot(BoxCox(milk.ts, lambda = lambda), main='Milk Production Time Series transformed using Box-Cox with lambda = -0.2210074', ylab='Box-Cox with lambda = -0.2210074')

acf2(milkLog, main='ACF & PACF for Logged Milk Production Series')
acf2(BoxCox(milk.ts, lambda = lambda), main='ACF & PACF for Milk Production Series transformed using Box-Cox')

# 2. Estimate the trend

## 2.1 Using a linear model

summary(fit <- lm(milkLog~time(milkLog)))
tsplot(milkLog, main='Linear Regression over time applied to Milk Production Logged Series', ylab='Log of Milk Production', col=4, lwd=2)
abline(fit)
tsplot(resid(fit), main="Detrended Milk Production Time Series")
acf2(resid(fit), main="ACF & PACF of the Detrended Milk Production Time Series")

summary(fitPoly <- lm(milkLog ~ poly(time(milkLog), 2, raw=TRUE)))
tsplot(milkLog, main='Polinomial Regression over time applied to Milk Production Logged Series', ylab='Log of Milk Production', col=4, lwd=2)
tsplot(resid(fitPoly), main="Detrended Milk Production Time Series")
acf2(resid(fitPoly), main="ACF & PACF of the Detrended Milk Production Time Series")

# 3. Remove season effects

## 3.1 Using a linear model:

milkLog.stl <- stl(milkLog, "periodic")
plot(milkLog.stl)  # original series

milkLog.sa <- seasadj(milkLog.stl)
tsplot(milkLog, main="Milk Production Logged series", ylab="Milk production (pounds per cow)")  # original series
autoplot(milkLog.sa, main="Milk Production seasonal adjusted logged series",  ylab="Milk production (pounds per cow)") # seasonal adjusted
acf(milkLog.sa)
ggseasonplot(milkLog.sa, main="Seasonal plot: Milk production")

summary(fit <- lm(milkLog.sa ~ time(milkLog.sa)))
tsplot(resid(fit), main="Residuals")
acf(resid(fit))

summary(fit <- lm(milkLog.sa ~ poly(time(milkLog.sa), 2, raw=TRUE)))
tsplot(resid(fit), main="Residuals of detrended, seasonal adjusted series", ylab="Residuals")
acf2(resid(fit))

remainder.ts <- milkLog.stl$time.series[,'remainder']
acf2(remainder.ts, main="ACF & PACF remainder")
checkresiduals(remainder.ts)

milk.arima = arima(resid(fit), order = c(1,0,0), include.mean=F)
milk.arima

acf2(milk.arima$resid[-1])


# Analyse the residuals

residuals <- milk.arima$resid[-1]
plot(residuals)
plot.ts(residuals)
acf2(residuals, main="ACF & PACF AR(1) residuals")
hist(residuals)
checkresiduals(residuals)

# Estimação dos modelos:

str(milk.ts)
milk.ts=as.xts(milk.ts)
str(milk.ts)
index(milk.ts)

milkLog <- as.xts(milkLog)
test.set <- last(milkLog,'12 month')
train.set <- window(milkLog,start=index(first(milkLog)), end=index(last(test.set)) - 1)
plot(train.set, main="Time plot of the train set")
plot(test.set, main="Time plot of the test set")
acf2(train.set, main="ACF & PACF of the train set")

adf.test(train.set, alternative = "stationary")

var(as.ts(train.set))
ndiffs(as.ts(train.set), test="kpss")
ndiffs(as.ts(train.set), test="adf")
ndiffs(as.ts(train.set), test="pp")

tsplot(as.ts(train.set), main="Milk Production (pounds per cow)")
tsplot(diff(as.ts(train.set)), main="Milk Production - First Difference", ylab="First Difference Series")

var(diff(as.ts(train.set)))

par(mfrow=c(2,1))
acf(train.set, main="", ylab="ACF before 1st diff")
acf(diff(as.ts(train.set)), main="", ylab="ACF 1st diff")
par(mfrow=c(1,1))
acf(resid(fit))

nsdiffs(as.ts(train.set))

par(mfrow=c(2,1))
plot(diff(as.ts(train.set),12), main="Series after 12th diff", ylab="")
acf(diff(as.ts(train.set),12), main="", ylab="ACF 12th diff")
par(mfrow=c(1,1))
var(diff(as.ts(train.set),12))

par(mfrow=c(2,1))
plot(diff(diff(as.ts(train.set),12)), main="Series after 1st and 12th diffs", ylab="")
acf(diff(diff(as.ts(train.set),12)), main="", ylab="ACF 1st and 12th diffs")
par(mfrow=c(1,1))
var(diff(diff(as.ts(train.set),12)))

adf.test(diff(diff(as.ts(train.set),12)), alternative = "stationary")
# p-values is small => reject h0 => stationary

train.set.diff <- diff(diff(as.ts(train.set),12))
acf2(train.set.diff, main="ACF & PACF after 1st and 12th diffs", ylab="ACF 1st and 12th diffs")

test.set.diff <- diff(diff(as.ts(test.set),12))

# Looking at seasonal part of the ACF & PACF, we can see that there are high correlation at year one on the ACF, and then the other years (2, 3 and 4) are not that significant. The seasonality of the PACF decreases slowly. This indicates an MA(1) model for the seasonal part.

model.sma1 <- sarima(train.set, 0,1,0,0,1,1,12) # the residuals are not normal, the ACF shows some correlation and the p-values are very small indicating that the null hypothesis of the residuals being uncorrelated should be rejected.
hist(train.set)
model.sma1
acf2(model.sma1$fit$residuals)
checkresiduals(model.sma1$fit)

# This seams to explain the seasonal part of the series, now we will try to model the other part.

model.ma1sma1 <- sarima(train.set, 0,1,1,0,1,1,12)
model.ma1sma1
acf2(model.ma1sma1$fit$residuals)
checkresiduals(model.ma1sma1$fit)

model.ma2sma1 <- sarima(train.set, 0,1,2,0,1,1,12)
model.ma2sma1
acf2(model.ma2sma1$fit$residuals)
checkresiduals(model.ma2sma1$fit)

model.ma5sma1 <- sarima(train.set, 0,1,5,0,1,1,12)
model.ma5sma1
acf2(model.ma5sma1$fit$residuals)
checkresiduals(model.ma5sma1$fit)

model.ar1sma1 <- sarima(train.set, 1,1,0,0,1,1,12)
model.ar1sma1
acf2(model.ar1sma1$fit$residuals)
checkresiduals(model.ar1sma1$fit)

model.ar5sma1 <- sarima(train.set, 5,1,0,0,1,1,12)
model.ar5sma1
acf2(model.ar5sma1$fit$residuals)
checkresiduals(model.ar5sma1$fit)

model.ar1ma1sma1 <- sarima(train.set, 1,1,1,0,1,1,12)
model.ar1ma1sma1
acf2(model.ar1ma1sma1$fit$residuals)
checkresiduals(model.ar1ma1sma1$fit)


# Auto arima:
auto.arima(train.set)
model.auto <- sarima(as.ts(train.set), 0,1,1,0,1,1,12)
model.auto

# Forecasting with SARIMA
pred.auto <- sarima.for(train.set, 12, 0,1,0,0,1,1,12)
accuracy(pred.auto$pred, test.set)


pred.auto <- sarima.for(train.set, 12, 0,1,1,0,1,1,12)

## Results

accuracy(pred.auto$pred, test.set)

high <- pred.auto$pred + 1.96*pred.auto$se
low <- pred.auto$pred - 1.96*pred.auto$se

inter <- data.frame(low, high)

inter[,(ncol(inter)+1)] = (coredata(test.set) > inter$low & coredata(test.set) < inter$high)[,1]
inter[,(ncol(inter)+1)] = coredata(test.set)
inter[,(ncol(inter)+1)] = pred.auto$pred
inter[,(ncol(inter)+1)] = abs(coredata(test.set) - pred.auto$pred)/coredata(test.set)

colnames(inter)[3] = "Contains"
colnames(inter)[4] = "Observed value"
colnames(inter)[5] = "Forecast"
colnames(inter)[6] = "Percentual Error"
inter

pred.original <- exp(pred.auto$pred)  # Exp of predictions


# Exponential Smoothing:

# To compare the forecast accuracy of the three methods, we use cross-validation:
e1 <- tsCV(milkLog, ses, h=12)
e2 <- tsCV(milkLog, holt, h=12)
e3 <- tsCV(milkLog, holt, damped=TRUE, h=12)
e4 <- tsCV(milkLog, hw, seasonal="additive", h=12)
e5 <- tsCV(milkLog, hw, seasonal="multiplicative", h=12)
e6 <- tsCV(milkLog, hw, damped=TRUE, seasonal="additive", h=12)
e7 <- tsCV(milkLog, hw, damped=TRUE, seasonal="multiplicative", h=12)

# Compare MSE:
mean(e1^2, na.rm=TRUE)
mean(e2^2, na.rm=TRUE)
mean(e3^2, na.rm=TRUE)
mean(e4^2, na.rm=TRUE)
mean(e5^2, na.rm=TRUE)
mean(e6^2, na.rm=TRUE)
mean(e7^2, na.rm=TRUE)
# Compare MAE:
mean(abs(e1), na.rm=TRUE)
mean(abs(e2), na.rm=TRUE)
mean(abs(e3), na.rm=TRUE)
mean(abs(e4), na.rm=TRUE)
mean(abs(e5), na.rm=TRUE)
mean(abs(e6), na.rm=TRUE)
mean(abs(e7), na.rm=TRUE)

# Additive
hwAdd <- hw(subset(as.ts(milkLog),end=length(as.ts(milkLog))-12), seasonal="additive", h=12)
hwAdd$model
accuracy(subset(hwAdd$fitted,start=length(hwAdd$fitted)-11), test.set)

autoplot(as.ts(milkLog)) +
  autolayer(hwAdd, series="HW additive", PI=FALSE)+
  guides(colour=guide_legend(title="Milk production forecasts"))+
  ylab("Milk Production (pounds per cow)")
summary(hwAdd)
checkresiduals(hwAdd)
Box.Ljung.Test(hwAdd$residuals, lag = NULL, main = NULL)

# Multiplicative
hwMulti <- hw(subset(as.ts(milkLog),end=length(as.ts(milkLog))-12), seasonal="multiplicative", h=12)
hwMulti$model
accuracy(subset(hwMulti$fitted,start=length(hwMulti$fitted)-11), test.set)

autoplot(as.ts(milkLog)) +
  autolayer(hwMulti, series="HW multiplicative", PI=FALSE)+
  guides(colour=guide_legend(title="Milk production forecasts"))+
  ylab("Milk Production (pounds per cow)")
summary(hwMulti)
checkresiduals(hwMulti)
Box.Ljung.Test(hwMulti$residuals, lag = NULL, main = NULL)

Box.Ljung.Test(model.ma1sma1$fit$residuals, lag = NULL, main = NULL)

hwAddD <- hw(subset(as.ts(milkLog),end=length(as.ts(milkLog))-12), seasonal="additive", damped=TRUE, h=12)
hwAddD$model
hwMultiD <- hw(subset(as.ts(milkLog),end=length(as.ts(milkLog))-12), seasonal="multiplicative", damped=TRUE, h=12)
hwMultiD$model

## State space models

# For the logged series, the ets function picks a model with multiplicative errors and seasonal component and additive trend.
fit.ets <- ets(train.set)
summary(fit.ets)
checkresiduals(fit.ets)

fit.ets %>% forecast(h=12) %>%
  autoplot() +
  ylab("Milk Production Time Series (pounds per cow)")

# When using the original series (before the log transformation) ets picks a model with multiplicative errors and seasonal component and additive trend.
test.set.ts <- last(milk.ts,'12 month')
train.set.ts <- window(milk.ts,start=index(first(milk.ts)), end=index(last(test.set)) - 1)
fit.ets.original <- ets(train.set.ts)
summary(fit.ets.original)
checkresiduals(fit.ets.original)

fit.ets.original %>% forecast(h=12) %>%
  autoplot() +
  ylab("Milk Production Time Series (pounds per cow)")



# Comparing ARIMA and Exponential Smoothing models

tsplot(as.ts(milkLog), main="Real values")

## ARIMA

arima.model <- Arima(train.set, order=c(0,1,1), seasonal=c(0,1,1))
arima.model %>% forecast(h=12) %>%
  autoplot() +
  ylab("Milk Production Time Series (pounds per cow)")

## Exponential Smoothing

autoplot(as.ts(milkLog)) +
  autolayer(hwMulti, series="HW multiplicative", PI=FALSE)+
  guides(colour=guide_legend(title="Milk production forecasts"))+
  ylab("Milk Production (pounds per cow)")

## State Space Models

fit.etc %>% forecast(h=12) %>%
  autoplot() +
  ylab("Milk Production Time Series (pounds per cow)")
