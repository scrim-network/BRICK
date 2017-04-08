## sneasyF.R
##
##==============================================================================
## Copyright 2017 Tony Wong, Alexander Bakker
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
#
# Klaus Keller (klaus@psu.edu)
# Brian Tuttle (btuttle@psu.edu)
# Nathan Urban
# Marlos Goes
# Department of Geosciences, Penn State
#
# Implements SNEASY, Simple Nonlinear EArth SYstem model:
#   Composite DOECLIM, Carbon Cycle model, and MOC box model
#
# Modified for BRICK by Tony Wong (twong@psu.edu), 4 April 2017
# ------------------------------------------------------------------
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
# ------------------------------------------------------------------
# CCM.f90:
#   Nonlinear impulse response carbon/climate model
#   based on the NICCS model from Hooss et al. (2001)
#   original version by DMRicciuto 7/16/2004
#
#   See model details in:
#     Ricciuto, D. M., K. J. Davis, and K. Keller (2008), A Bayesian
#         calibration of a simple carbon cycle model: The role of
#         observations in estimating and reducing uncertainty, Global
#         Biogeochem. Cycles, 22, GB2030, doi:10.1029/2006GB002908.
# ------------------------------------------------------------------

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

# load ocean anomaly table.
if(!exists("oceanomtable")) oceanomtable = t(read.table("../data/sneasy_tmp/anomtable.txt"))
ntab1 = nrow(oceanomtable)
ntab2 = ncol(oceanomtable)

#####      load SNEASY model shared library     #####
if(.Platform$OS.type == "unix") {
    dyn.load("../fortran/sneasy.so")
	#dyn.load("/Users/axw322/codes/sneasy/sneasy/sneasy.so")
} else {
    dyn.load("../fortran/sneasy")
}

setup.sneasy = function()
{
	fout = .Fortran( "init_sneasy",
            tab_ni = as.integer(ntab1),
            tab_nj = as.integer(ntab2),
            oceanomtab = as.matrix(oceanomtable,ntab1,ntab2)
            )
}

cleanup.sneasy = function()
{
	fout = .Fortran("fin_sneasy")
}

# SNEASY coupled climate, carbon cycle, moc model
# inputs:   climate sensitivity (S),
#           ocean vertical diffusivity (kappa),
#           soil respiration (Q10),
#           carbon fertilization (beta),
#           thermocline diffusion (eta),
#           hydrological sensitivity of North Atlantic surface (hydsens)
# outputs:  time [years]
#           MOC strength [Sv]
#           radiative forcing [W/m^2]
#           atmospheric CO2 [ppm]
#           atmosphere-ocean flux [W/m^2]
#           global surface temperature [K]
#           global ocean heat [10^22 J]

sneasy = function(S       = 2.7,
                   kappa   = 2.9,
				   alpha   = 1.0,
				   Q10     = 4.2,
				   beta    = 0.9,
				   eta     = 23,
				   hydsens = 0.03,
				   init.CO2= 280,
				   init.MOC= 20,
				   tstep   = 1,
				   mod.time,
				   forcing.co2,
				   forcing.aero,
				   forcing.other)
{
	n = length(mod.time)

	# call Fortran SNEASY
	# sneasy.so/dll must be already dynamically loaded (see above this function)
	fout = .Fortran( "run_sneasy_model",
			ns = as.integer(n),
            timestep = as.double(tstep),
			Clim_sens = as.double(S),
			Oc_vert_diff = as.double(kappa),
			Aer_ampl = as.double(alpha),
            Soil_resp = as.double(Q10),
            CO2_fert = as.double(beta),
            TC_diff = as.double(eta),
            HS_NA_surf = as.double(hydsens),
            init_MOC = as.double(init.MOC),
            init_CO2 = as.double(init.CO2),
            carbon_emis = as.double(forcing.co2),
			aer_forc = as.double(forcing.aero),
            other_forc = as.double(forcing.other),
			MOC_strength = as.double(rep(0,n)),
			radiative_forc = as.double(rep(0,n)),
			ATM_CO2 = as.double(rep(0,n)),
			atm_oc_flux = as.double(rep(0,n)),
			GL_surface_temp = as.double(rep(0,n)),
			GL_ocean_heat = as.double(rep(0,n))
	)

	model.output = list(
                    time = mod.time,
                    forcing = fout$radiative_forc,
                    temp = fout$GL_surface_temp,
                    ocheat = fout$GL_ocean_heat,
                    co2 = fout$ATM_CO2,
                    ocflux = fout$atm_oc_flux,
                    moc = fout$MOC_strength
                   )

	return(model.output)
}
##==============================================================================
## End
##==============================================================================
