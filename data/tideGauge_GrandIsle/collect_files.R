# collect_files.R
#
# round up all of the yearly tide gauge files
# They were accessed from https://tidesandcurrents.noaa.gov/waterlevels.html?id=8761724
# on 11 Jan 2017.
# Then calculate MLE estimates of GEV parameters.
# I accounted for leap years; every four years, the end date decreases by 1,
# because the maximum number of days of data that can be downloaded per CSV
# is 365.
#==============================================================================
dat.dir <- '../data/tideGauge_GrandIsle/'
files.tg <- list.files(path=dat.dir,pattern='^[CO-OPS_8761724_hr]')

# all dates/times are relative to first observation (1980-11-11 16:00 GMT)
data <- read.csv(paste(dat.dir,'CO-OPS__8761724__hr.csv',sep=''))
origin <- as.POSIXlt(data$Date.Time, tz='GMT')[1]

# first file to start the data array
data.new <- read.csv(paste(dat.dir,files.tg[1],sep=''))
time.new <- as.POSIXlt(data.new$Date.Time, tz='GMT')
time.new.rel <- difftime(time.new, origin, unit='days')
data.new[,1] <- time.new.rel
data.all <- data.new

for (ff in 2:length(files.tg)) {
  data.new <- read.csv(paste(dat.dir,files.tg[ff],sep=''))
  time.new <- as.POSIXlt(data.new$Date.Time, tz='GMT')
  time.new.rel <- difftime(time.new, origin, unit='days')
  data.new[,1] <- time.new.rel
  data.all <- rbind(data.all,data.new)
}

# sort the tide gauge data array by the time.relative column
# (not necessarily in order, and possibly some overlap between adjacent year files)
data.all.sort <- data.all[order(data.all$Date.Time),]

# Get water levels in mm
data.all.sort$Water.Level <- data.all.sort$Water.Level*1000

# bin by year, or otherwise make the timestamp actually useful
# Note: to get actual datestamp back, just do data.all.sort$Date.TIme + origin
days.year <- 365.25 # number of days in a year, on average
nyear <- floor(as.numeric(max(data.all.sort$Date.Time)/days.year))
data.all.sort$Water.Level.avg <- data.all.sort$Water.Level
data.all.sort$lsl <- data.all.sort$Water.Level
data.all.sort$lsl.avg <- data.all.sort$Water.Level
data.all.sort$lsl.avg.rel <- data.all.sort$Water.Level
lsl.max <- rep(NA,nyear)  # annual block maxima

# One method:
# - Calculate annual means
# - Subtract these off
# - Calculate annual block maxima
# - Fit GEV
# Alternative method:
# - Use NOAA/USACE SLR trend of 9.24 mm/year (https://tidesandcurrents.noaa.gov/est/curves.shtml?stnid=8761724)
# - Subtract this SLR trend from the tidge gauge data.
# - Normalize the result to have mean of zero.
# - Calculate annual block maxima, fit GEV to these, as usual.
# - Yields result which is lower than the annual means method (both are nice ways
# to account for sea-level rise.)
slr.rate <- 9.24/365  # sea level rise per day
slr <- as.numeric(slr.rate*data.all.sort$Date.Time)

data.all.sort$lsl <- data.all.sort$Water.Level - slr
data.all.sort$lsl <- data.all.sort$lsl - mean(data.all.sort$lsl)
for (tt in 1:nyear) {
  iyear <- which(data.all.sort$Date.Time <  days.year*tt &
                 data.all.sort$Date.Time >= days.year*(tt-1))
  data.all.sort$lsl.avg[iyear] <- rep(mean(data.all.sort$Water.Level[iyear]),length(iyear))
  data.all.sort$lsl.avg.rel[iyear] <- data.all.sort$Water.Level[iyear]-data.all.sort$lsl.avg[iyear]
  # Method 1:
  #lsl.max[tt] <- max(data.all.sort$lsl.avg.rel[iyear])
  # Method 2:
  lsl.max[tt] <- max(data.all.sort$lsl[iyear])
}

# install.packages("extRemes")
# install.packages("fExtremes")
# install.packages("ismev")
# install.packages("lubridate")
# install.packages("zoo")
library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)

# Fit GEV, re-make return level plot. Still funky?
#todo
gev.fit <- fevd(coredata(lsl.max))   # extRemes package
gev.fit2 <- gevFit(coredata(lsl.max))   # fExtremes package
gev.fit3 <- gev.fit(coredata(lsl.max), show = FALSE)   # ismev package

print(year.res.max.fit2@fit$par.ests)

storm_surgeL <- 1e5 # desired storm surge level
q = seq(0,1,length.out= storm_surgeL +1)  # quantile array

# Find closed-form solution of GEV fit
fit_q_year = qgev(q, gev.fit2@fit$par.ests[1], gev.fit2@fit$par.ests[2], gev.fit2@fit$par.ests[3])
fit_q_year = fit_q_year[fit_q_year< max(fit_q_year)]
fits <- dgev(lsl.max, gev.fit2@fit$par.ests[1], gev.fit2@fit$par.ests[2], gev.fit2@fit$par.ests[3])

# Find which q is the 100-yr value
num = which(q <= 0.99)
num.max = which.max(num)
year100prob <- num.max +1
# check to make sure the year100prob value is the 100-yr value (10^-2)
print(1-q[year100prob])

plot(fit_q_year/1000,log10(1-q[2:100001]),type='l', xlim=c(0,6))
  points(lsl.max/1000, log10(1-fits))

#==============================================================================
# End
#==============================================================================
if(FALSE){ # scratch work
storm_surgeL <- 1e5 # desired storm surge level
q = seq(0,1,length.out= storm_surgeL +1)  # quantile array
fit1 = qgev(q, xi=gev.mle1$shape, mu=gev.mle1$location, beta=gev.mle1$scale)
fit1 = fit1[fit1< max(fit1)]

# Find which q is the 100-yr value
num = which(q <= 0.99)
num.max = which.max(num)
year100prob <- num.max +1
# check to make sure the year100prob value is the 100-yr value (10^-2)
print(fit1[year100prob])
}
