
library("xlsx")

A <- read.table("D:/Paveloom/Курсовая/output1")
B <- read.table("D:/Paveloom/Курсовая/output2")

# options(digits=16)

# Размер выборки
N = nrow(A)

# Вектора абсциссы
x = seq(1:N-1)

# Коррелограмма (код на Фортране) (весь период, 5110 дней)

par(mfrow=c(1,2))
acf(D[,2], lag.max = 5110, xlab = "Time Lag, Days, k", ylab = "Autocorrelation Coefficient, r(k)", main = "Correlogram of Total Solar Irradiance (TSI) at 1-AU, 5110 days, ACF Function")
plot(x, A[,2], xlim = c(0, 5110), type='l', xlab = "Time Lag, Days, k", ylab = "Autocorrelation Coefficient, r(k)", main = "Correlogram of Total Solar Irradiance (TSI) at 1-AU, 5110 days")
dev.off()

# Два месяц

acf(D[,2], lag.max = 62, xlab = "Time Lag, Days, k", ylab = "Autocorrelation Coefficient, r(k)", main = "Correlogram of Total Solar Irradiance (TSI) at 1-AU, 62 days")
lines(x, A[,2])

# Один месяца

acf(D[,2], lag.max = 31, xlab = "Time Lag, Days, k", ylab = "Autocorrelation Coefficient, r(k)", main = "Correlogram of Total Solar Irradiance (TSI) at 1-AU, 31 days")
lines(x, A[,2])

# Периодограмма (код на Фортране) (до 0.1 по частоте)

par(mfrow=c(1,2))
# TSA::periodogram(D[,2], xlim = c(0.17064846416382252730E-03,0.1), ylim = c(0,1e+6), lwd = 0.01, ylab = "Periodogram, I(v)", xlab = "Frequency, 1/Days, v", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 10 days, Periodogram Function")
plot(B[,1], B[,3], type='l', xlim = c(5860,10), ylim = c(0,1e+5), log = "x", ylab = "Periodogram, I(T)", xlab = "Period, Days, T, Log Scale", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 0.1- Frequencies")
plot(B[,2], B[,3], type='l', xlim = c(0.17064846416382252730E-03,0.1),  ylim = c(0,1e+5), ylab = "Periodogram, I(v)", xlab = "Frequency, 1/Days, v", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 0.1- Frequencies")
dev.off()

# 10 дней

par(mfrow=c(2,1))
plot(B[,1], B[,3], type='l', xlim = c(1,10), ylim = c(0,1.35e+3), ylab = "Periodogram, I(T)", xlab = "Period, Days, T", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 10 days")
plot(B[,1], B[,3], type='l', xlim = c(1,10), log = "x", ylim = c(0,1.35e+3), ylab = "Periodogram, I(T)", xlab = "Period, Days, T, Log Scale", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 10 days")
dev.off()

# 31 день

# par(mfrow=c(1,2))
# TSA::periodogram(D[,2], xlim = c(0,0.00529010238907849829351535836177), ylim = c(0, 3.8e+7), lwd = 0.01, ylab = "Periodogram, I(v)", xlab = "Frequency, 1/Days, v", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 31 days, Periodogram Function")

par(mfrow=c(2,1))
plot(B[,1], B[,3], type='l', xlim = c(1,31), ylim = c(0, 6e+3), ylab = "Periodogram, I(T)", xlab = "Period, Days, T", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 31 days")
plot(B[,1], B[,3], type='l', xlim = c(1,31), log = "x", ylim = c(0, 6e+3), ylab = "Periodogram, I(T)", xlab = "Period, Days, T, Log Scale", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 31 days")
dev.off()

# dev.off()

# 62 дня

# par(mfrow=c(1,2))
# TSA::periodogram(D[,2], xlim = c(0,0.01058020477815699658703071672355), ylim = c(0, 3.8e+7), lwd = 0.01, ylab = "Periodogram, I(v)", xlab = "Frequency, 1/Days, v", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 62 days, Periodogram Function")

par(mfrow=c(2,1))
plot(B[,1], B[,3], type='l', xlim = c(1,62), ylim = c(0, 2.2e+4), ylab = "Periodogram, I(T)", xlab = "Period, Days, T", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 62 days")
plot(B[,1], B[,3], type='l', xlim = c(1,62), log = "x", ylim = c(0, 2.2e+4), ylab = "Periodogram, I(T)", xlab = "Period, Days, T, Log Scale", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 62 days")
dev.off()

# Долгопериодичная (20 лет)

par(mfrow=c(3,1))

plot(B[,1], B[,3], type='l', xlim = c(1,58600), log = "x", ylim = c(0, 70), ylab = "Periodogram, I(T)", xlab = "Period, Days, T, Log Scale", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 58600 days")
abline(v = 4502, col = "blue")
text(3500, 65, "4502", cex = 1, col = "blue")

plot(B[,1], B[,3], type='l', xlim = c(3000,7000), log = "x", ylim = c(40, 60), ylab = "Periodogram, I(T)", xlab = "Period, Days, T, Log Scale", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 58600 days, Zoomed")
abline(v = 4502, col = "blue")

plot(B[,1], B[,3], type='l', xlim = c(1,58600), log = "xy" , ylim = c(1e-8, 7e+2), ylab = "Periodogram, I(T), Log Scale", xlab = "Period, Days, T, Log Scale", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 58600 days")
abline(v = 4502, col = "blue")

dev.off()

# Полный диапазон частот

par(mfrow=c(2,1))
plot(B[,2], B[,3], type='l', xlim = c(0,1), ylim = c(0,5e+3), ylab = "Periodogram, I(v)", xlab = "Frequency, 1/Days, v", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 1- Frequencies")
plot(B[,2], B[,3], type='l', xlim = c(0,1),  log = "y", ylim = c(0.0001,1e+7), ylab = "Periodogram, I(v), Log Scale", xlab = "Frequency, 1/Days, v", main = "Periodogram of Total Solar Irradiance (TSI) at 1-AU, 1- Frequencies")
dev.off()

# Общий вид исходных данных (TSI от Юлианских дней)

D <- read.table("D:/Paveloom/Курсовая/SORCE Total Solar Irradiance - Daily Average, Time Series - 2.txt")

plot(D[,1], D[,2], ylim = c(1356,1363), xlab = 'Time Interval, Julian Day Number', ylab = 'Total Solar Irradiance (TSI) at 1-AU, W/m^2')

# Период?

# 27,5 дней
plot(D[,1], D[,2], ylim=c(1360.4,1360.75), xlim=c(2454500,2454750), xlab = 'Time Interval, Julian Day Number', ylab = 'Total Solar Irradiance (TSI) at 1-AU, W/m^2')

for(i in -5:5){
  abline(v = 2454580.5 + i*27, col = 'blue')
}
abline(v = 24545780.5, col = 'red')

abline(v = 2454607, col = "blue")
abline(v = 2454634.5, col = "blue")
abline(v = 2454662, col = 'blue')
abline(v = 2454689.5, col = 'blue')
abline(v = 2454717, col = 'blue')
abline(v = 2454579.5, col = 'blue')
abline(v = 2454552, col = 'blue')
abline(v = 2454524.5, col = 'blue')
abline(v = 2454497, col = 'blue')


# Тест 1

input <- read.table("D:/Paveloom/Курсовая/input")

y = input[,2]
var1 = var(y)

test_A <- read.table("D:/Paveloom/Курсовая/output1-test")

test_v = seq(1:15)
test_delta_t = 1
test_x = (test_v-1)*test_delta_t

plot(test_x, test_A[,2], type='l', xlab = "Time Interval, Julian Day Number", ylab = "Correlogram of Total Solar Irradiance (TSI) at 1-AU, W^2/m^4")


# Тест 2

sinus <- read.table("D:/Paveloom/Курсовая/sinus")

plot(sinus[,1], sinus[,2], type = 'l')

sinus1 <- read.table("D:/Paveloom/Курсовая/sinus1")
sinus2 <- read.table("D:/Paveloom/Курсовая/sinus2")

plot(sinus1[,1], sinus1[,2], type = 'l')
plot(sinus2[,1], sinus2[,2], type = 'l', xlim = c(0,10))

# Тест 3 (проверка пиков)

library(mgcv)

par(mfrow=c(3,2))
plot(D[,1], D[,2], xlim = c(2452696,2452696+5110), ylim = c(1356,1363), main = "TSI Overview for 5110 days", xlab = 'Time Interval, Julian Day Number', ylab = 'TSI at 1-AU, W/m^2')
plot(x, A[,2], type='l', xlim = c(0,5110), ylim = c(0,0.5), main = "SAC with L from 0 to 5110 days", xlab = "Time Lag, Days, k", ylab = "SAC, r(L)")
plot(x, A[,2], type='l', xlim = c(0,5110), ylim = c(-0.1,0.1), main = "SAC with L from 0 to 5110 days, Zoomed", xlab = "Time Lag, Days, k", ylab = "SAC, r(L)")
plot(x, A[,2], type='l', xlim = c(3400,3800), main = "SAC with L from 3400 to 3800 days", xlab = "Time Lag, Days, k", ylab = "SAC, r(L)")
plot(x, A[,2], type='l', xlim = c(1900,2300), main = "SAC with L from 1900 to 2300 days", xlab = "Time Lag, Days, k", ylab = "SAC, r(L)")
plot(x, A[,2], type='l', xlim = c(3800,4200), ylim = c(0,0.5), main = "SAC with L from 3800 to 4200 days", xlab = "Time Lag, Days, k", ylab = "SAC, r(L)")
dev.off()

plot(B[,1]/100, B[,2]/100, xlim = c(0,1),type = 'l')

for (i in 0:199){
  l_b = 1000 + i*10
  r_b = l_b + 10
g<-gam(y~s(t),data=data.frame(t=x[l_b:r_b],y=A[l_b:r_b,2]))
lines(x[l_b:r_b],g$fitted.values,col=3)}

plot(x, A[,2], type='l', xlim = c(1000,2000), ylim = c(0,0.2))
g<-gam(y~s(t),data=data.frame(t=x[1000:2000],y=A[1000:2000,2]), fit = TRUE)
lines(x[1000:2000],g$fitted.values,col=3)

help(gam)

x <- 1:10
y <- c(2,4,6,8,7,8,14,16,18,20)
lo <- loess(y~x)
plot(x,y)
xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
lines(xl, predict(lo,xl), col='red', lwd=2)
library(ggplot2)
ggplot() + xlim(1000,2000) + ylim(0,0.2) + geom_point(aes(x,A[,2])) + geom_smooth(aes(x,A[,2]),span = 0.001)

# x, y: the x and y coordinates of the hull points
# n: the number of points in the curve.
bezierCurve <- function(x, y, n=10)
{
  outx <- NULL
  outy <- NULL
  
  i <- 1
  for (t in seq(0, 1, length.out=n))
  {
    b <- bez(x, y, t)
    outx[i] <- b$x
    outy[i] <- b$y
    
    i <- i+1
  }
  
  return (list(x=outx, y=outy))
}

bez <- function(x, y, t)
{
  outx <- 0
  outy <- 0
  n <- length(x)-1
  for (i in 0:n)
  {
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]
  }
  
  return (list(x=outx, y=outy))
}
help(gam)
# Example usage
x <- c(4,6,4,5,6,7)
y <- 1:6
plot(x, y, "o", pch=20)
points(bezierCurve(x,y,20), type="l", col="red")

# Log Scale Function:

function(xlim=c(1,1000), xlog=TRUE, xbase=10, 
         ylim=c(1,1000), ylog=TRUE, ybase=10, ...)
{
  if(xlog) # rounding the X-axis limits on the log scale
  {
    xlim[1] <- floor(log(xlim[1], base=xbase))
    xlim[2] <- round(log(xlim[2], base=xbase))
    xbreaks <- xlim[1]:xlim[2]
  }
  if(ylog) # rounding the Y-axis limits on the log scale
  { 
    ylim[1] <- floor(log(ylim[1], base=ybase))
    ylim[2] <- round(log(ylim[2], base=ybase))
    ybreaks <- ylim[1]:ylim[2]
  }
  
  ### the empty plot into which the axes will be drawn ###
  plot(xlim, ylim, type="n", axes=FALSE, frame=TRUE, ...)
  
  if(xlog) # plotting the X-axis tickmarks and grids
  {
    for(x in xbase^xbreaks)
    {
      subx <- log(seq(from=x, to=x*xbase, length=xbase) , base=xbase )
      abline(v=subx, col="grey")
    }
    axis(side=1, at=xbreaks, labels=xbase^xbreaks, tck=0.02)
  }
  else  axis(side=1, tck=0.02)
  
  if(ylog) # plotting the Y-axis tickmarks and grids
  {
    for(y in ybase^ybreaks)
    {
      suby <- log(seq(from=y, to=y*ybase, length=ybase), base=ybase)  
      abline(h=suby, col="grey")
    }
    axis(side=2, at=ybreaks, labels=ybase^ybreaks, las=2, tck=0.02)
  }
  else axis(side=2, las=2, tck=0.02)
}


# Тест 4(вау, здесь все время была функция для функции автокорреляции)

library(tseries)
acf(D[,2], lag.max = 100)
lines(x,A[,2])

# (вау, и для периодограммы тоже)

par(mfrow=c(1,2))
TSA::periodogram(D[,2], xlim = c(0,0.01), ylim = c(0,3e+6), lwd = 0.01)
plot(B[,1]/100,B[,2], type = 'l', xlim = c(0,58.6), ylim = c(0,3e+6))
dev.off()

# Тест 5
var(D[,2])
