##==============================================================================
## Interface between differential optimization "DEoptim" and SIMPLE model
##
## Must use standard parameter names as given in SIMPLE_calib_driver.R
##
## Questions? Tony Wong <twong@psu.edu>
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


minimize_residuals_simple = function( parameters.in ,
																			parnames.in ,
																			forcing.temp.in,
																			oidx,
																			midx,
																			obs,
																			obs.err,
																			ind.norm.data
																			){

	a.simple    =parameters.in[match("a.simple"    ,parnames.in)]
  b.simple    =parameters.in[match("b.simple"    ,parnames.in)]
  alpha.simple=parameters.in[match("alpha.simple",parnames.in)]
  beta.simple =parameters.in[match("beta.simple" ,parnames.in)]
  V0       		=parameters.in[match("V0"     	   ,parnames.in)]

  simple.out = simple(a=a.simple, b=b.simple, alpha=alpha.simple,
                       beta=beta.simple, V0=V0, Tg=forcing.temp.in)

	## Subtract off normalization period
	itmp = ind.norm.data[match("gis",ind.norm.data[,1]),2]:ind.norm.data[match("gis",ind.norm.data[,1]),3]
	simple.out$sle.gis = simple.out$sle.gis - mean(simple.out$sle.gis[itmp])

	# optimize based on GIS data
	resid = sum(abs( (simple.out$sle.gis[midx$gis] - obs$gis[oidx$gis])/obs.err$gis ))

	if(is.nan(resid)) resid=Inf

  return(resid)
}

##==============================================================================
## End
##==============================================================================
