#  Loading Required Libraries
library(stats)
library(fpp2)
library(ggplot2)
library(seasonal)
library(forecast)

## Loading Data

Data = read.csv('Time Series Data.csv')
Data = data.frame(Data)


# Test for trend estimation of a time series
#Relating Ordering test
Relative_Ordering_Test = function(time_series){
  
  n = length(time_series)
  Q = 0
  
  # Total no. of discordents.....
  for (i in 1:n){
    for (j in (i+1):n)
    {
      if (!is.na(time_series[i]) && !is.na(time_series[j]) && time_series[i] > time_series[j])
        Q = Q + 1
    }
  }
  
  Tau = 1 - 4 * Q / (n * (n - 1))
  Var = 2 * (2 * n + 5) / (9 * n * (n - 1))
  
  Z = Tau/sqrt(Var)
  
  
  #  Let level of significance, alpha=0.05
  alpha = 0.05
  q = qnorm(1-alpha/2,0,1)
  
  if (Z <= q)    # No trend....
    return (1)   
  else 
    return (0)   
}



#  Test for randomness of a time series data.....
Turning_Point_Test = function(time_series)
{
  n = length(time_series)
  P = 0
  for (i in 2:(n-1))
  {
    if (((time_series[i]>time_series[i-1]) && (time_series[i]>time_series[i+1])) || ((time_series[i]<time_series[i-1]) && (time_series[i]<time_series[i+1])))
       P = P + 1
  }
  
  E = 2*(n-2)/3
  Var = (16*n - 29)/90
  
  Z = (P - E)/sqrt(Var)
  
  alpha = 0.05
  q = qnorm(1-alpha/2, 0, 1)
  
  if (Z <= q)
    return (1)
  else
    return (0)
  
}


#  To check the stationarity of a time series....
Dickey_Fuller_Test = function(time_series)
{
  result <- adf.test(time_series,alpha = 0.05)
  return(c(result$statistic,result$p.value))
}



#  Analysis done on the average temperature.....

View(Data)
avg_temp = ts(Data$Avg, start = as.Date("2016-01-01"), end = as.Date("2020-11-30"))
plot.ts(avg_temp)

###average temp vs timeplot

date_sequence <- seq(as.Date("2016-01-01"), by = "1 day", length.out = 1796)
ggplot(avg_temp, aes(x = date_sequence)) +
  geom_line(aes(y = avg_temp, color = "Original Data")) +
  scale_color_manual(values = c("Original Data" = "dodgerblue4", "15-point MA" = "lightslateblue", "30-point MA" = "black")) +
  labs(title = "Time Series Plot",
       x = "Year",
       y = "Average Temperature") +
  theme_minimal()


# ..................................TREND ESTIMATION........................
# We will use a 30-point and a 15-point MA to estimate the trend 
# We will also use differencing to remove the trend 
# And we will be sure by relative-ordering test of trend removal....


# 15-point and 30-point MA
trend_3point = stats::filter(avg_temp, rep(1,15)/15, sides = 2)
lines(trend_3point, col = "red")


trend_3point = stats::filter(avg_temp, rep(1,30)/30, sides = 2)
lines(trend_3point, col = "green")
legend("topleft", legend = c("MA_15", "MA_30", "avg_temp"),col = c("red", "green", "black"),lty = c(2, 3, 1),lwd = 2,box.lty = 0,box.lwd = 0.1)


library(ggplot2)
library(dplyr)

# Assuming 'Data' is your data frame with a column named 'Avg'
D <- data.frame(
  Date = seq(as.Date("2016-01-01"), as.Date("2020-11-30"), by = "days"),
)

# Convert to time series
avg_temp <- ts(Data$Avg, start = as.Date("2016-01-01"), end = as.Date("2020-11-30"))

# Create a data frame for ggplot
plot_data <- data.frame(Date = time(avg_temp), AvgTemp = as.vector(avg_temp))

# Plotting with ggplot
ggplot(plot_data, aes(x = Date, y = AvgTemp)) +
  geom_line(color = "black", size = 1) +
  geom_line(aes(y = stats::filter(AvgTemp, rep(1, 15)/15, sides = 2)), color = "red", linetype = 2, size = 1) +
  geom_line(aes(y = stats::filter(AvgTemp, rep(1, 30)/30, sides = 2)), color = "green", linetype = 3, size = 1) +
  labs(x = "Date", y = "Temperature") +
  scale_linetype_manual(name = "Lines", values = c(1, 2, 3), labels = c("Avg_temp", "MA_15", "MA_30")) +
  theme_minimal() +
  theme(legend.position = "topleft")




# Differencing for trend removal........
a = 0
new_data = avg_temp
b = Turning_Point_Test(diff(new_data))    # A 1st order difference makes the series stationary
print(b)                                  # removing both trend and seasonality....
                          
#  We will use Dickey-Fuller test for to confirm stationarity 
Dickey_Fuller_Test(new_data)


plot.ts(diff(avg_temp), col = "blue")
abline(0,0)
hist(avg_temp, breaks = 500)

hist(diff(avg_temp ), xlim =c(-10,12), breaks = 100)   #Random ~ Nearly N(0.00133, 2.96756)



# .................Seasonality Estimation........................

# creating monthly data
dddd = numeric(60)
ind1 = 0
ind2 = 1
val = 0
while(ind2+ind1*30 <= 1796)
{
  while(ind2 <= 30)
  {
    val = val + avg_temp[ind2 + (ind1*30)]
    ind2 = ind2 + 1
  }
  ind2 = 1
  ind1 = ind1 + 1
  dddd[ind1] = val/30
  val = 0
}
plot.ts(dddd)
val = 0
ind1 = 0
ind2 = 1
yyyy = numeric(5)

while(ind1*12 +ind2 <= 60){
  while(ind2 <= 12 & 12*ind1+ind2 <= 59){
    val = val + dddd[12*ind1+ind2]
    ind2 = ind2 + 1
  }
  ind1 = ind1 + 1
  yyyy[ind1] = val/ind2
  ind2 = 1
  val = 0
}
data <- data.frame(year = seq_along(yyyy), value = yyyy)

# Plot using ggplot2
p <- ggplot(data, aes(x = year, y = value)) +
  geom_line() +
  labs(title = "Time Series Plot", x = "Year", y = "Average Value") +
  theme_minimal()

p + ylim(20,30)
plot.ts(yyyy,ylim = 20:30)


###### as yyyy value for every year is almost same (~23.4) we can say seasonal component is there.

#  From Relative-Ordering test, we confirm that the data has no trend......
plot.ts(avg_temp)

length(avg_temp)


nrow(Data)
train_data = Data[1:1456,]
View(train_data)


avg_temp_train = ts(train_data$Avg, start = as.Date("2016-01-01"), end = as.Date("2019-12-31"))
plot.ts(avg_temp_train)

summary(avg_temp_train)
class(avg_temp_train)

  #Actual Data ACF And PACF
acf(avg_temp, main="ACF For Average Temperature")
pacf(avg_temp, main="PACF For Average Temperature")

acf(diff(avg_temp), main="ACF of First Order Differencing")
pacf(diff(avg_temp), main="PACF of First Order Differencing ")



#  Model Building
arima.model = auto.arima(avg_temp, seasonal = TRUE)
arima.model

summary(arima.model)      # AIC=8580.3,   AICc=8580.36,   BIC=8618.75


arima.forecast = forecast(arima.model, h=365)
autoplot(arima.forecast, avg_temp)
#lines(avg_temp, col = "pink")

residuals <- resid(arima.forecast)
plot(residuals, type = "l", col = "blue", ylab = "Residuals", main = "Residual Plot")
hist(residuals, breaks = 100)

#######################


#######################


