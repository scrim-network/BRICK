##==============================================================================
#
#  -original file = "Read_GSIC_data.R"   Code written July 3, 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#  - Modified by: Tony Wong (twong@psu.edu) (June 2016)
#
#  -This program loads in temperature and Glacier & Small Ice Caps (GSIC) data for use in model
#       and uncertainty calculations.
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
