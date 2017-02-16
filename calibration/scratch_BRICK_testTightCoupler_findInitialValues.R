##==============================================================================
##	This script is for testing running ensembles of fully-forward coupled
##	simulations by using the results for the 1850 "initial condition" as the
##	uncertain initial condition in:
##		DOECLIM	--	T0, H0
##		GSIC		--	Gs0 (may need to adjust V0 higher)
##		SIMPLE	--	V0
##		TE			--	TE0
##		DAIS		--	n/a (initial ice sheet volume calculated from other parameters)
##
##	Questions? -- Tony Wong <twong@psu.edu
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

rm(list=ls())

setwd('~/codes/BRICK/calibration')

library(ncdf4)

# read BRICK v0.1 calibrated model simulations
filename.brick <- '../output_model/BRICK-model_physical_control_01Nov2016.nc'
ncdata <- nc_open(filename.brick)
  t.hind <- ncvar_get(ncdata, 'time_hind')
  t.proj <- ncvar_get(ncdata, 'time_proj')
  gmsl.ctrl.hind <- ncvar_get(ncdata, 'GlobalSeaLevel_hind')
  temp.ctrl.hind <- ncvar_get(ncdata, 'temp_hind')
  ocheat.ctrl.hind <- ncvar_get(ncdata, 'ocheat_hind')
  ais.ctrl.hind <- ncvar_get(ncdata, 'AIS_hind')
  gis.ctrl.hind <- ncvar_get(ncdata, 'GIS_hind')
  te.ctrl.hind <- ncvar_get(ncdata, 'TE_hind')
  gsic.ctrl.hind <- ncvar_get(ncdata, 'GSIC_hind')
  gmsl.ctrl.rcp26 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gmsl.ctrl.rcp45 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gmsl.ctrl.rcp85 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
nc_close(ncdata)

# read corresponding ensemble parameters
filename.parameters <- '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
ncdata <- nc_open(filename.parameters)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters <- t(parameters)
colnames(parameters) <- parnames


# Determine suitable bounds for the initial condtions in BRICK v0.2 (fully-
# forward integrating, single fortran function call)


# Gs0 is the initial sea-level contribution from GSIC. It is defined to be 0,
# but there is uncertainty. The previous calibration runs give an idea as to
# what this uncertainty should be.
Gs0.gsic <- gsic.ctrl.hind[1,] - mean(gsic.ctrl.hind[1,])

# V0.gsic is the initial GSIC volume. In the previous calibration runs, this was
# defined in 1960 (Wigley and Raper, 2005), but we need it in 1850 (when the
# BRICK-[other model] runs would start). The previous calibration (as described
# in Wong et al 2017, GMDD) starts GSIC in 1961 and integrates forward
# to 2009, and backward to 1850 to obtain the estimates below for V0.gsic.
V0.gsic <- parameters[,'V0.gsic'] +
            (parameters[,'Gs0'] - gsic.ctrl.hind[1,])

# V0.te is the initial TE sea-level contribution. It is defined to be 0, but
# (like Gs0) there is uncertainty. The previous calibratoin runs give an idea as
# to what this uncertainty is.
V0.te <- te.ctrl.hind[1,] - mean(te.ctrl.hind[1,])

# V0.simple is the initial GIS volume (not sea-level contribution). Previous
# calibration has this defined to be relative to 1961-1990.
i1961 <- which(t.hind==1961)
V0.simple <- parameters[,'V0'] +
            (gis.ctrl.hind[i1961,] - gis.ctrl.hind[1,])

print(paste('Gs0.gsic bounds: ',min(Gs0.gsic),' to ',max(Gs0.gsic),sep=''))
print(paste('V0.gsic bounds: ',min(V0.gsic),' to ',max(V0.gsic),sep=''))
print(paste('V0.te bounds: ',min(V0.te),' to ',max(V0.te),sep=''))
print(paste('V0.simple bounds: ',min(V0.simple),' to ',max(V0.simple),sep=''))


##==============================================================================


##==============================================================================
## End
##==============================================================================
