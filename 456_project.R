# Math 456 Project

rm(list=ls())
dev.off()

library(dplyr)
library(zoo)
library(tidyverse)

setwd("~/Downloads/wwu/wwu 25-26/math 456/456 Project")


df <- read_csv("sockeyeSalmon_countRiverConditions_PilgrimRvr_2003-2014_carey.csv")

# cleaning data by removing commas in counts and converting to standard date
# remove 2012 year due to abnormally short observation period
df_clean <- df %>%
  mutate(
    # Convert SockeyeCount to numeric, handling commas (if present)
    SockeyeCount = as.numeric(gsub(",", "", SockeyeCount)),
    # Convert Date and extract Year (assuming 2-digit year is 2000s, e.g., 03 -> 2003)
    Date = as.Date(Date, format = "%m/%d/%y"),
    Year = format(Date, "%Y")
  ) %>%
  # Remove 2012 observations
  filter(Year != "2012")

# Filter and Shorten Remaining Years

# Target length is 63 observations, based on the user's specification for 2014.
TARGET_N <- 63

df_processed <- df_clean %>%
  group_by(Year) %>%
  mutate(
    # Total observations in the current year
    N_current = n(),
    # Total observations to remove
    N_remove = N_current - TARGET_N,
    # Number of observations to remove from the beginning (e.g., for 21 removals, remove 10 from start)
    N_start_remove = floor(N_remove / 2),
    # Index of the row within the current year (1, 2, 3, ...)
    Row_Index = row_number()
  ) %>%
  # Filter step:
  # 1. Row_Index must be greater than the number of rows removed from the start.
  # 2. Row_Index must be less than or equal to the index of the last row to keep.
  filter(
    Row_Index > N_start_remove,
    Row_Index <= N_current - (N_remove - N_start_remove)
  ) %>%
  ungroup() %>%
  # Select the original columns, removing the helper columns used for indexing
  select(Date, SockeyeCount, RiverTemperature, RiverHeight, Year)

# Verification

# This code block confirms that all remaining years now have 63 observations.
print("--- Final Counts by Year ---")
df_processed %>%
  group_by(Year) %>%
  summarise(N_final = n())

sockeye_count <- df_processed$SockeyeCount 
temp_river <- df_processed$RiverTemperature
height_river <- df_processed$RiverHeight

plot.ts(temp_river)

#Dealing with Missing Values
sum(is.na(sockeye_count))
sum(is.na(temp_river))
sum(is.na(height_river))
# perfroming linear interpolation on the missing values (basically jsut averaging
# the two values around the missing point)
sockeye_count <- na.approx(sockeye_count)
temp_river[692:693] = 13
temp_river <- na.approx(temp_river)
height_river <- na.approx(height_river)


# Weekly Averages of Data
N <- length(sockeye_count)
# Calculate the number of full weeks
num_full_weeks <- floor(N / 7)

# Split the vector into chunks of 7
# This approach ignores the remainder (the last partial week)
weekly_averages_sockeye <- tapply(sockeye_count[1:(num_full_weeks * 7)], 
                                 rep(1:num_full_weeks, each = 7), 
                                 mean)
weekly_averages_temp <- tapply(temp_river[1:(num_full_weeks * 7)], 
                                  rep(1:num_full_weeks, each = 7), 
                                  mean)
weekly_averages_height <- tapply(height_river[1:(num_full_weeks * 7)], 
                                  rep(1:num_full_weeks, each = 7), 
                                  mean)

sockeye_count <- as.vector(round(weekly_averages_sockeye, digits = 2))
temp_river <- as.vector(round(weekly_averages_temp, digits = 2))
height_river <- as.vector(round(weekly_averages_height, digits = 2))

plot.ts(sockeye_count, ylab="Sockeye Count")
# shift the data up 1 to deal with 0s, this avoids infinity after log transformation
log_sockeye <- ts(log(sockeye_count + 1))
plot(log_sockeye, ylab="Log Sockeye Count")




#ARIMA Modeling
auto.arima(log_sockeye)

# Single Spectrum Analysis (SSA)
library("Rssa")

# change to weekly average 
n = length(log_sockeye)
log_sockeye_ssa = ssa(log_sockeye, kind='1d-ssa')
# Diagnostics for grouping
plot(log_sockeye_ssa) # Eigenvalues
plot(log_sockeye_ssa, type = "vectors") # Eigenvectors
plot(log_sockeye_ssa, type = "paired") # Pairs of eigenvectors
plot(wcor(log_sockeye_ssa)) # w-correlation matrix plot

#list containing groups
group_log_sockeye = list(c(1,4,5), c(2,3), c(6,7))
#here we have one trend grouping and two periodicity groupings

num_components = length(group_log_sockeye)
recon_log_sockeye = reconstruct(log_sockeye_ssa, groups = group_log_sockeye)

#isolating trend plot:
#group_log_sockeye2 = list(c(1,2,3,4))
#here we have one trend grouping and one periodicity groupings

#par(mfrow=c(1,1))
#num_components2 = length(group_log_sockeye2)
#recon_log_sockeye2 = reconstruct(log_sockeye_ssa, groups = group_log_sockeye2)
#for(i in 1:num_components2){
#  plot(recon_log_sockeye2[[i]], main=main_vector[i])
#}

# Finding periodicities
parestimate(log_sockeye_ssa, list(c(2,3)), method = "pairs")$periods
parestimate(log_sockeye_ssa, list(c(6,7)), method = "pairs")$periods

par(mfrow=c(1,1))
main_vector = c("Trend", "Periodicity #1", "Periodicity #2")
ylab_vector = c("Trend Reconstruction", "Periodicity #1 Reconstruction", 
                "Periodicity #2 Reconstruction")
for(i in 1:num_components){
  plot(recon_log_sockeye[[i]], main=main_vector[i], ylab=ylab_vector[i])
}

d_log_sockeye = rep(0,n)
for(i in 1:num_components){
  d_log_sockeye = d_log_sockeye + recon_log_sockeye[[i]]
}

par(mfrow=c(1,1))
plot(cbind(log_sockeye, d_log_sockeye), plot.type="single", col=c("black","red"), 
     main="log_sockeye Data", lwd=2, xlab="Week", ylab="log_sockeye Concentration")
legend("topright", legend=c("Original", "Deterministic"),
       col=c("black", "red"), lwd=2, cex=1.2, bg="lightblue")

# Forecasting
n_ahead = 9 # Number of forecasts
for_log_sockeye_ssa = rforecast(log_sockeye_ssa, groups = group_log_sockeye, len = n_ahead, only.new = FALSE)
for_log_sockeye = rep(0,(n+n_ahead))
for(i in 1:num_components){
  for_log_sockeye = for_log_sockeye + for_log_sockeye_ssa[[i]]
}
plot(cbind(log_sockeye, for_log_sockeye), plot.type="single", col=c("black","red"), 
     main="Log Sockeye Data", lwd=2, xlab="Week", ylab="Log Sockeye Concentration")
legend("topright", legend=c("Original", "Deterministic"),
       col=c("black", "red"), lwd=2, cex=1.2, bg="lightblue")

# The actual forecast values
for_log_sockeye[(n+1):(n+n_ahead)]

#residuals
resid <- d_log_sockeye - log_sockeye
hist(resid, main="Histogram of Residuals", xlab="Residuals")
qqnorm(resid)
qqline(resid)
shapiro.test(resid)

acf(resid, main="ACF of Residuals")
pacf(resid, main="PACF of Residuals")

### Portmanteau test
pvalues = double(10)
for(lag in 1:10){
  pvalues[lag] = Box.test(resid, lag=lag)$p.value
}
par(mfrow=c(1,1))
plot(pvalues, main="P-Values of the Portmanteau test", xlab="Lag", ylab="P-Value",ylim=c(0,1))
abline(h=0.05, lty=2, col="blue")

auto.arima(resid, )


### RMSE value on the log scale
errors <- (d_log_sockeye - log_sockeye)^2
error_max <- max(errors)
error_min <- min(errors)

RMSE_ssa <- sqrt((sum(errors) / 99))
# Normalized RMSE
norm_RMSE <- RMSE_ssa / (error_max - error_min)



#ARIMA analysis 
plot.ts(log_sockeye)
length(log_sockeye)
#Seasonal differencing with d = 48
sockeye_sdiff = diff(log_sockeye, lag=9)
plot.ts(sockeye_sdiff)
#First differencing of the seasonal-differenced data.
sockeye_sdiff_diff = diff(sockeye_sdiff)
plot.ts(sockeye_sdiff_diff)
#ACF and PACF plot
acf(sockeye_sdiff_diff, lag.max=99)
pacf(sockeye_sdiff_diff, lag.max=99)
#ASK
#Determining an appropriate order p for the AR(p) model using BIC (Bayesian Information Criterion).
library("forecast")
sockeye_bic=auto.arima(log_sockeye, d=1, D=1, max.p=2, max.q=1, 
                       max.P=1, max.Q=1, stationary=TRUE, seasonal=FALSE, ic="bic", 
                       allowmean=FALSE)
sockeye_bic
n_ahead = 9 # Number of forecasts. It is h_max in the description.
#check number of predictions with Jalen 
fcast= predict(sockeye_bic, n.ahead=n_ahead)
fcast$pred

#plot residuals of ARIMA model 

sockeye_residuals <- residuals(sockeye_bic)
qq_sockeye <- qqnorm(sockeye_residuals)
qqline(sockeye_residuals)
acf(sockeye_residuals)
pacf(sockeye_residuals)
plot_residuals<- plot(sockeye_residuals)

#rmse of ARIMA model.
fitted_vals <- fitted(sockeye_bic)
#fits ARIMA model estimates for sockeye count 
#at each time point to compare to real data for RMSE
rmse_ARIMA <- sqrt(mean((log_sockeye - fitted_vals)^2))
rmse_ARIMA











