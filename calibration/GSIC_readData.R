##==============================================================================
#
#  -file = "Read_GSIC_data.R"   Code written July 3, 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in temperature and Glacier & Small Ice Caps (GSIC) data for use in model
#       and uncertainty calculations.
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
#       -Glacier & Small Ice Caps (GSIC) Data from National Snow & Ice Data Center (NSIDC)
#       -http://nsidc.org/forms/G10002_or.html?major_version=1
#       -Institute of Arctic and Alpine Research University of Colorado
#       -updated version - Occasional Paper No. 58 (2005)
#
##==============================================================================

## Historical global mean sea level contribution from Glacier and Small Ice Cap melt
dat = read.csv("../data/GSICobservations_UPDATED.csv", skip = 1)
obs.gsic.time = dat[1:43, 1] #1961-2003
obs.gsic = dat[1:43, 5]/1000 # m of melt contribution to sea level rise (Note -- data are in mm)
obs.gsic.err = dat[1:43,9]/1000 # m (standard error) (Note -- data are in mm)

idx = compute_indices(obs.time=obs.gsic.time, mod.time=mod.time)
oidx.gsic = idx$oidx; midx.gsic = idx$midx

##==============================================================================
## End
##==============================================================================
