##==============================================================================
# doeclimF.R
#
# Nathan M. Urban (nurban@psu.edu)
# Department of Geosciences, Penn State
#
#	Modified for BRICK framework by Tony Wong (twong@psu.edu)
#
# Implements DOECLIM, a simple climate model
#
# DOECLIM is a 0-dimensional energy balance model (EBM) for the
# atmosphere coupled to a 1-dimensional diffusive ocean.  The
# model outputs temperature and ocean heat content time series
# as well as ocean heat fluxes.  See:
#
#  Elmar Kriegler, "Imprecise probability analysis for integrated
#  assessment of climate change", Ph.D. thesis, Potsdam (2005).
#   http://opus.kobv.de/ubp/volltexte/2005/561/
#  (Fortran port by Marlos Goes and Nathan Urban.)
#
# The model is implemented in Fortran and called from R.  This
# file also contains functions to load and process forcing data
# and model output.  Any pre/post-processing of input/output
# should be done in R or otherwise separate from the main
# model subroutine, which for computational efficiency
# should not perform any file I/O itself.  The Fortran model
# must be implemented as a standalone subroutine.
#
# For further information on R/Fortran coupling, see:
#
#   http://www.stat.umn.edu/~charlie/rc/
#   http://math.acadiau.ca/ACMMaC/howtos/Fortran_R.html
#   http://cran.r-project.org/doc/manuals/R-exts.pdf (chapter 5)
#
##==============================================================================
#
# Input forcing:
#  forcing					CO2 (RCP8.5) and aerosol, non-CO2 radiative forcings
#										(co2,nonco2.land/ocean,aerosol.land/ocean,solar.land/ocean,
#                    volc.land/ocean,tot.land,tot.ocean)
#
# Input parameters:
#  S								climate sensitivity (inc temp from 2xCO2) [deg C]
#  kappa            ocean vertical heat diffusivity [cm2/s]
#
# Output:
#  time             years when model was run [year]
#  temp             surface temperature anomaly [deg C]
#  ocheat           ocean heat uptake [10^22 J]
#  ocheat.mixed     mixed layer ocean heat anomaly [10^22 J]
#  ocheat.interior  interior ocean heat anomaly [10^22 J]
#                   (Note: ocheat = ocheat.mixed + ocheat.interior)
#  ocheatflux.mixed     heat uptake of the ocean mixed layer [W/m2]
#  ocheatflux.interior  heat uptake of the ocean interior [W/m2]
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

# convert annual ocean heat flux (W/m^2) to cumulative ocean heat content anomaly (10^22 J)
flux.to.heat = function(heatflux.mixed, heatflux.interior)
{
	flnd = 0.29 # area land fraction
	fso = 0.95 # ocean area fraction of interior
	secs.per.year = 31556926
	earth.area = 510065600 * 10^6
	ocean.area = (1-flnd)*earth.area
	powtoheat = ocean.area*secs.per.year / 10^22 # in 10^22 J/yr

	heat.mixed = cumsum(heatflux.mixed) * powtoheat
	heat.interior = fso * cumsum(heatflux.interior) * powtoheat
	ocean.heat = heat.mixed + heat.interior

	return(list(ocean.heat=ocean.heat, heat.mixed=heat.mixed, heat.interior=heat.interior))
}

## load DOECLIM model shared library
# dyn.load("../fortran/doeclim.so")
if(.Platform$OS.type == "unix") {
    dyn.load("../fortran/doeclim.so")
} else {
    dyn.load("../fortran/doeclim")
}

# DOECLIM climate model (0D EBM atmosphere + 1D diffusive ocean)
# inputs: climate sensitivity (S), ocean vertical diffusivity (kappa), aerosol forcing scale factor (alpha)
# outputs: annual global temperature (temp, K) and total ocean heat (ocheat, 10^22 J), mixed-layer and interior ocean heat (ocheat.mixed and ocheat.interior, 10^22 J), atmosphere->mixed and mixed->interior heat fluxes (ocheatflux.mixed and ocheatflux.interior, W/m^2)
doeclimF = function(
							S     = 3.1,
							kappa = 3.5,
							forcing.total,
							mod.time
							)
{
	n = length(mod.time)

	# call Fortran DOECLIM
	# doeclim.so must be already dynamically loaded (see above this function)
	fout = .Fortran( "run_doeclim",
			ns = n,
			time_out = as.double(mod.time),
			forcing_in = as.double(forcing.total),
			t2co_in = as.double(S),
			kappa_in = as.double(kappa),
			temp_out = as.double(rep(0,n)),
			heatflux_mixed_out = as.double(rep(0,n)),
			heatflux_interior_out = as.double(rep(0,n))
		)

	ocheat = flux.to.heat(fout$heatflux_mixed, fout$heatflux_interior)

	model.output = list(time=mod.time, temp=fout$temp_out, ocheat=ocheat$ocean.heat,
											ocheat.mixed=ocheat$heat.mixed, ocheat.interior=ocheat$heat.interior,
											ocheatflux.mixed = fout$heatflux_mixed, ocheatflux.interior = fout$heatflux_interior)

	return(model.output)
}
