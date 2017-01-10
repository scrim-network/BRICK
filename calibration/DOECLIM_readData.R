##==============================================================================
#
#  -file = "DOECLIM_readData.R"   Original GSIC code written July 3, 2014
#								   Modified DOECLIM version written 9 June 2016
#  - Original GSIC version author: Kelsey Ruckert (klr324@psu.edu)
#	 - Modified DOECLIM version: Tony Wong (twong@psu.edu)
#
#  -This program loads in temperature and ocean heat data for use in DOECLIM model
#       and uncertainty calculations:
#		+ surface temperature, HadCRUT4 global means (1850-2012), Morice et al (2012)
#		+ ocean heat, 0-3000 m (1953-1996), Gouretski and Koltermann (2007)
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
#   -NOTE: This file contains data that is sourced into the other programs. Information
#       regarding this data can be found below:
#
#       -RCP8.5 is used to create temperature simulations to 2300
#       -RCP8.5 simulates "Business as usual" and is similar to the A2 scenario
#       -The RCP8.5 temperatures are from the CNRM-CM5 model
#           (Centre National de Recherches Meteorologiques)
#       -These RCP8.5 temperatures closly resemble historical temperatures
#           from the National Oceanic and Atmospheric Administration (NOAA)
#
#       -Annual global land and ocean temperature anomalies (C)
#       -Anomalies with respect to the 20th century average
#       -http://www.ncdc.noaa.gov/cag/time-series/global
#
##==============================================================================

# HADCRUT4 annual global mean surface temperature
# Note: all ensemble member files have same time series of ucnertainties, so just
# grabbing the first one.
dat = read.table("../data/HadCRUT.4.4.0.0.annual_ns_avg.txt")
obs.temp = dat[,2]
obs.temp.time = dat[,1]
dat = read.table("../data/HadCRUT.4.4.0.0.annual_ns_avg_realisations/HadCRUT.4.4.0.0.annual_ns_avg.1.txt")
obs.temp.err = dat[,3]

# Normalize temperature anomaly so 1961-1990 mean is 0
ibeg=which(obs.temp.time==1850)
iend=which(obs.temp.time==1870)
obs.temp = obs.temp - mean(obs.temp[ibeg:iend])

# annual global ocean heat content (0-3000 m)
dat = read.table("../data/gouretski_ocean_heat_3000m.txt",skip=1)
obs.ocheat = dat[,2]
obs.ocheat.time = dat[,1]
obs.ocheat.err = dat[,3]

idx = compute_indices(obs.time=obs.temp.time, mod.time=mod.time)
oidx.temp = idx$oidx; midx.temp = idx$midx

idx = compute_indices(obs.time=obs.ocheat.time, mod.time=mod.time)
oidx.ocheat = idx$oidx; midx.ocheat = idx$midx

##==============================================================================
## End
##==============================================================================
