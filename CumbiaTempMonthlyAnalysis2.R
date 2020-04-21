#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## INITIAL PLOT ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

y=scan('CumbriaTempMonthly.txt',skip=2,nlines=240)
y
date1 = seq(2000, 2019+11/12,1/12)
yts = ts(y[1:216]) #Holding back 24 data values for forecasting.
date = seq(2000, 2017+11/12, 1/12) #2017 due to data held back.
plot(date, yts, main = "Monthly Temperature at Newton Rigg (10% Data Omitted)", type = 'p', col = 2, xlab = "Year",ylab = "Temperature")
lines(date, yts)
n = length(yts)
time = seq(1:n)
ymod1 = lm(yts~time)
yfit1 = ymod1$fitted
lines(date, yfit1)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## RESIDUALS OF LINEAR REGRESSION ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

yres1 = ymod1$residuals
ser = ts(yres1)
plot(ser, main = "Residuals of monthly temperature at Newton Rigg", type = 'l', col = 2) #Notice x scale change

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## PERIODOGRAM OF RESIDUALS ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

yspec = spec.pgram(ser, pad = 4, main = "Periodogram of the residuals")
plot(yspec$freq, yspec$spec, type = 'l', col = 4, main = "Periodogram of the residuals") #Large peak at 1/12= 0.9=08333

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## SEASONAL REGRESSION ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

time = seq(1:n)
month = c(rep(seq(1:12), 18))
fmonth = as.factor(month)
ymod1 = lm(yts~time+fmonth)
summary(ymod1) 
yfit1 = ymod1$fitted
plot(date, yts, type = 'p')
lines(date, yfit1, col = 4) #Good fit
yres1 = ymod1$residuals
plot(date, yres1, main = "Residuals of the fitted model", xlab = "Residuals", ylab = "Date")
lines(date, yres1) #Stationary.
acf(yres1, lag.max = 40) #MA(4)
pacf(yres1, lag.max = 40)
yspec = spec.pgram(yres1, pad = 4, main = "Raw periodogram of residuals of seasonal regression")
plot(yspec$freq, yspec$spec, type = 'l', main = "Periodogram of residuals of seasonal regression")
p = rep(0, 10)
for (i in 1:10){
  p[i] = Box.test(yres3, i, type = "Ljung-Box")$p.value
}
p.min = 0.05 #level of significance
plot(p, ylim = c(0, 1), main = "P-values for Ljung-box statistic", xlab = "lag", ylab = "p-value")
abline(h = p.min, col = 4, lty =2)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## HARMONIC REGRESSION ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

ymod2 = lm(yts~time + cos(2*pi*time/12) + sin(2*pi*time/12))
summary(ymod2)
yfit2 = ymod2$fitted
plot(date, yts, main = "Newton Rigg temperature with fitted model")
lines(date, yfit2, col =4)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## RESIDUALS OF HARMONIC REGRESSION ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

yres2 = ymod2$residuals
ts.plot(ymod2$res, main = "Residuals of harmonic regression") #Small x scale

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## PERIODOGRAM OF RESIDUALS ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

yspec = spec.pgram(yres2, pad = 4, main = "Raw periodogram of residuals of harmonic regression")
plot(yspec$freq, yspec$spec, type = 'l', main = "Periodogram of residuals of harmonic regression") #Small x scale again, shows good fit
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## ACF AND PACF OF RESIDUALS ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

acf(yres2, lag.max = 20, main = "ACF of Residuals") #MA(1)?
pacf(yres2, lag.max = 20, main = "PACF of Residuals", ylim = c(-0.5, 1)) #AR(1)? or AR(4) or AR(5)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## LJUNG-BOX STATISTICS ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

p = rep(0, 10)
for (i in 1:10){
  p[i] = Box.test(yres2, i, type = "Ljung-Box")$p.value
}
p.min = 0.05 #level of significance
plot(p, ylim = c(0, 1), main = "P-values for Ljung-box statistic", xlab = "lag", ylab = "p-value")
abline(h = p.min, col = 4, lty =2) #All below line

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## MA(4) FIT ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

ma4fit = arima(yres2, order = c(0, 0, 4), include.mean = FALSE)
ma4fit #Standard errors quite large, ar1 is best for this
tsdiag(ma4fit) #Higher p values compared to ma3, ar 1 may be better

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## RESIDUALS OF FIT ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

x = residuals(ma4fit)
plot(date, x, main = "Residuals after MA fit") #Looks more like white noise
lines(date, x)
acf(x, main = "ACF after MA fit", lag.max = 50)
pacf(x, main = "PACF after MA fit", lag.max = 50) #Both look like white noise

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## LJUNG-BOX OF RESIDUALS ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

p = rep(0, 10)
for (i in 1:10){
  p[i] = Box.test(x, i, type = "Ljung-Box")$p.value
}
p.min = 0.05 #level of significance
plot(p, ylim = c(0, 1), main = "P-values for Ljung-box statistic", xlab = "lag", ylab = "p-value")
abline(h = p.min, col = 4, lty =2) # All high #note lag 4 on last few plots


#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## FORECASTING ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

arpred = predict(ma4fit, n.ahead = 30)
fcast = arpred$pred
low = fcast - 1.96*arpred$se
upp = fcast + 1.96*arpred$se
## Plot ##
ts.plot(yres2, col = 1, xlim = c(0, n + 30), main = "Forecasting the MA(4) component")
lines(fcast, col = 2)
lines(low, lty = 3, col = 3)
lines(upp, lty = 3, col = 3)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## FORECASTING THE ORIGINAL DATA ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

timef = seq((n+1),(n+30))
coefficients(ymod2)
fts = coefficients(ymod2)[1] + coefficients(ymod2)[2]*timef +
  coefficients(ymod2)[3]*cos(2*pi*timef/12) +
  coefficients(ymod2)[4]*sin(2*pi*timef/12) +
  fcast
date = seq(2000, 2019 + 11/12, 1/12)
plot(yts, main="Monthly Temperature Forecast at Newton Rigg", col=1, xlim = c(0, n + 30), type = 'p')
lines(timef,fts,col=2)
points(seq((n+1),(n+24)),y[(n+1):(n+24)],col=1)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## FORECAST UNCERTAINTY ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

ftslow = coefficients(ymod2)[1] + coefficients(ymod2)[2]*timef +
  coefficients(ymod2)[3]*cos(2*pi*timef/12) +
  coefficients(ymod2)[4]*sin(2*pi*timef/12) +
  low
ftshigh = coefficients(ymod2)[1] + coefficients(ymod2)[2]*timef +
  coefficients(ymod2)[3]*cos(2*pi*timef/12) +
  coefficients(ymod2)[4]*sin(2*pi*timef/12) +
  upp

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## ADDING THESE TO THE PLOT ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

lines(timef,ftslow,col=3)
lines(timef,ftshigh,col=3)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## ZOOMING IN ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

timew = seq((n-10), (n + 30)) #(n-10 to n+24 in lecture examples, requires some fiddling)
plot(timew, y[timew], main = "Monthly Temperature Forecast at Newton Rigg", ylim = c(0, 25), ylab = "Temperature", xlab = "Time")
lines(yts)
lines(timef, fts, col = 2)
points(seq((n + 1), (n + 10)), y[(n + 1): (n + 10)], col = 1)
lines(timef, ftslow, col = 3)
lines(timef, ftshigh, col = 3)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
## ADDING IN UNCERTAINTY FROM LINEAR TREND ##
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

ftslow2 = coefficients(ymod2)[1] + 
  coefficients(ymod2)[2]*timef +
  -1.96*summary(ymod1)$coef[2,2]*seq(1,30) +
  coefficients(ymod2)[3]*cos(2*pi*timef/12) +
  coefficients(ymod2)[4]*sin(2*pi*timef/12) +
  low 
ftshigh2 = coefficients(ymod2)[1] + 
  coefficients(ymod2)[2]*timef +
  +1.96*summary(ymod1)$coef[2,2]*seq(1,30) +
  coefficients(ymod2)[3]*cos(2*pi*timef/12) +
  coefficients(ymod2)[4]*sin(2*pi*timef/12) +
  upp
lines(timef,ftslow2,lty = 3, col=4)
lines(timef,ftshigh2,lty = 3, col=4)
lines(yts, col = 1)
