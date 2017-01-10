##==============================================================================
## Paragraph [7] from Wigley and Raper (2005, doi:10.1029/2004GL021238)
## (right above Equation 3) gives beta(t) = d2Gs/dtdT.
##
## This routine uses finite differences to estimate an upper bound on the
## parameter beta0 by calculating the time series beta(t).
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================
## Copyright 2016 Tony Wong, Alexander Bakker
## This file is part of BRICK (Building blocks for Relevant Ice and Climate
## Knowledge). BRICK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## BRICK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

## Read in temperature and sea level data
data = read.csv("../GSIC-MAGICC/NOAA_IPCC_RCPtempsscenarios.csv")

#Historical time frame and temperatures from NOAA
year.temp=data[,1]

## Historical global mean sea level contribution from Glacier and Small Ice Cap melt
GSIC = read.csv("../GSIC-MAGICC/GSICobservations_UPDATED.csv", skip = 1)
year.Gs = GSIC[1:43, 1] #1961-2003
Gs.sle = GSIC[1:43, 5]/1000 # m of melt contribution to sea level rise (Note -- data are in mm)

## Reduce both to only when there are GSIC data (year.Gs)
ibeg = which(year.temp==year.Gs[1])
iend = which(year.temp==year.Gs[length(year.Gs)])

time = year.Gs
temp = data[ibeg:iend,2]
ntime=length(time)-1

##==============================================================================
## Option 1 -- Use a forward difference to approximate dGs/dt
dGsdt = rep(0,ntime-1)
dGsdt = (Gs.sle[2:ntime]-Gs.sle[1:(ntime-1)]) / (time[2:ntime]-time[1:(ntime-1)])

        # then use another forward difference to approximate d(dGs/dt)/dT
d2GsdTdt = rep(0,ntime-2)
d2GsdTdt = (dGsdt[2:(ntime-1)]-dGsdt[1:(ntime-2)]) / (temp[2:(ntime-1)]-temp[1:(ntime-2)])

print(paste('max(abs(beta))=',max(abs(d2GsdTdt))))
##==============================================================================

##==============================================================================
## Option 2 -- same, reverse order
dGsdT = rep(0,ntime-1)
dGsdT = (Gs.sle[2:ntime]-Gs.sle[1:(ntime-1)]) / (temp[2:ntime]-temp[1:(ntime-1)])

d2GsdtdT = rep(0,ntime-2)
d2GsdtdT = (dGsdT[2:(ntime-1)]-dGsdT[1:(ntime-2)]) / (time[2:(ntime-1)]-time[1:(ntime-2)])

print(paste('max(abs(beta))=',max(abs(d2GsdtdT))))
##==============================================================================

##==============================================================================
## Option 3 -- centered differences, with respect to t, then T
dGsdt = rep(NA,ntime-2)
dGsdt = (Gs.sle[3:ntime]-Gs.sle[1:(ntime-2)]) / (time[3:ntime]-time[1:(ntime-2)])
    # note: now the i-th element of dGsdt is the derivative at time t(i+1),
    # so this corresponds to temperature T(i+1),
    # so we have derivatives at times t(2) through t(ntime-1),
    # most transparent way to do this is not fancy indexing, but resetting what
    # our temperature is, to match the derivative dGsdt, before calculating the
    # second derivative
temp2 = temp[2:(ntime-1)]
ntime2= ntime-2

d2GsdTdt = rep(NA,ntime2-2)
d2GsdTdt = (dGsdt[3:ntime2]-dGsdt[1:(ntime2-2)]) / (temp2[3:ntime2]-temp2[1:(ntime2-2)])

    # note: getting an inf because one of temp[t]-temp[t-2] is 0
print(paste('max(abs(beta))=',max(abs(d2GsdTdt[which(is.finite(d2GsdTdt))]))))
##==============================================================================

##==============================================================================
## Option 4 -- centered differences, with respect to T, then t
dGsdT = rep(NA,ntime-2)
dGsdT = (Gs.sle[3:ntime]-Gs.sle[1:(ntime-2)]) / (temp[3:ntime]-temp[1:(ntime-2)])
    # note: now the i-th element of dGsdT is the derivative at time t(i+1),
    # so this corresponds to temperature T(i+1),
    # so we have derivatives at times t(2) through t(ntime-1),
    # most transparent way to do this is not fancy indexing, but resetting what
    # our temperature is, to match the derivative dGsdT, before calculating the
    # second derivative
time2 = time[2:(ntime-1)]
ntime2= ntime-2

d2GsdtdT = rep(NA,ntime2-2)
d2GsdtdT = (dGsdT[3:ntime2]-dGsdT[1:(ntime2-2)]) / (time2[3:ntime2]-time2[1:(ntime2-2)])

    # note: getting an inf because one of temp[t]-temp[t-2] is 0
print(paste('max(abs(beta))=',max(abs(d2GsdtdT[which(is.finite(d2GsdtdT))]))))
##==============================================================================


##==============================================================================
## End
##==============================================================================
