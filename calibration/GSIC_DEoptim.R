##==============================================================================
##  -File = "GSIC_minimize_resid_DE.R"
##  -Function of the sum of the absolute value residuals of the
##  -Wigley and Raper 2005 Model to source into DEoptim
##  -The point is to minimize the residuals to find the optimal parameters
##  -For the model
##
##	- Code separately from gsic
##
##  -Original Author: Kelsey Ruckert (klr324@psu.edu)
##  -Modifed by: Tony Wong (twong@psu.edu)
##==============================================================================
##  -Original: May 27, 2014
##  -Modified: June 2016
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

minimize_residuals_gsic = function( parameters.in,
                                    parnames.in,
                                    forcing.temp.in,
                                    oidx,
                                    midx,
                                    obs,
                                    obs.err,
                                    ind.norm.data
                                    ){

		beta0  = parameters.in[match("beta0"  ,parnames.in)]
		V0.gsic= parameters.in[match("V0.gsic",parnames.in)]
		n      = parameters.in[match("n"      ,parnames.in)]
		Gs0    = parameters.in[match("Gs0"    ,parnames.in)]

    gsic.out = gsic_magiccF(beta0=beta0, V0=V0.gsic, n=n, Gs0=Gs0, Tg=forcing.temp.in)

		## Subtract off normalization period model GSIC output as the zero point
		itmp = ind.norm.data[match("gsic",ind.norm.data[,1]),2]:ind.norm.data[match("gsic",ind.norm.data[,1]),3]
		gsic.out.norm = gsic.out - mean(gsic.out[itmp])

    err.sum=sum(abs( (obs$gsic[oidx$gsic] - gsic.out.norm[midx$gsic])/obs.err$gsic[oidx$gsic] ))

    if(is.nan(err.sum)) err.sum=Inf

    return(err.sum)
 }

##==============================================================================
## End
##==============================================================================
