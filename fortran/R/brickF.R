#===============================================================================
# BRICK Fortran90
#
# Call full model, stepped forward together
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
#
# Input forcing:
#
#   forcing             CO2 (RCP8.5) and aerosol, non-CO2 radiative forcings
#                       (co2,nonco2.land/ocean,aerosol.land/ocean,solar.land/ocean,
#                       volc.land/ocean,tot.land,tot.ocean)
#
# Other input:
#
#   slope.Ta2Tg         slope of linear regression between Antarctic temperature
#                       and global surface temperature anomalies
#   intercept.Ta2Tg     intercept of this linear regression (if slope and intercept
#                       are not provided, assumed to be Antarctic temperature already)
#
#===============================================================================
# Input parameters:
#
#	tstep		model time step in years
# -- DOECLIM --
#   S           climate sensitivity (inc temp from 2xCO2) [deg C]
#   kappa       ocean vertical heat diffusivity [cm2/s]
#
# -- GSIC-MAGICC --
#   beta0       initial mass balance sensitivity (how long it takes GSIC to respond to
#               increasing temps) [m/yr/C]
#   V0          initial volume = max(Gs) [meter sle]
#   n           exponent for area-volume scaling [-]
#   Gs0         Gs[1]: the corrected corresponding sea-level rise in 1961 [m]
#   Teq         equilibrium temperature (at which there is no change) [deg C]
#
# -- GIS-SIMPLE --
#   a           sensitivity of equilibrium volume Veq [m sle/degC]
#   b           equilibrium volume Veq [m sle] for temperature Tg = 0
#   alpha       sensitivity of exponential decay rate (1/tau)
#   beta        exponential decay rate [1 / degC] at Tg = 0
#   Vgrl_0      initial ice-sheet volume [m sle]
#
# -- TE --
#   a           sensitivity of TE equilibrium [m/K]
#   b           equilibrium TE [m] for temperature anomaly Tg = 0
#   invtau      1/time-scale of exponential decay (e-folding time) [years^-1]
#   TE_0        initial sea level
#
# -- DAIS --
#   anto.a      Sensitivity Ta [degC/degC]
#   anto.b      Ta(Tg=0)
#   b0          Undisturbed bed height at the continent center [m]
#   slope       Slope of ice sheet bed before loading [-]
#   mu          Profile parameter for parabolic ice sheet surface (related to ice stress) [m0.5]
#   h0          hr(Ta=0): Height of runoff line at Ta = 0 [m]
#   c           Sensitivity of Height of runoff line (hr) [m/degC]
#   P0          P(Ta=0): Annual precipitation for Ta = 0 [m (ice equivalent)]
#   kappa       Coefficient for the exponential dependency of precipitation on Ta [degC-1]
#   nu          Proportionality constant relating runoff decrease with height to precipitation [m^(-1/2) yr^(-1/2)]
#   f0          Proportionality constant for ice flow at grounding line [m/yr]
#   gamma       Power for the relation of ice flow speed to water depth [-]
#   alpha       Partition parameter for effect of ocean subsurface temperature on ice flux [-]
#   Toc_0       Present-day, high latitude ocean subsurface temperature [degC]
#   Rad0        Reference ice sheet radius [m]
#   dSL0        Initial sea level rate
#
#===============================================================================
# Output:
#
#   time                years when model was run [year]
#
# -- DOECLIM --
#   temp                surface temperature anomaly [deg C]
#   ocheat              ocean heat uptake [10^22 J]
#   ocheat.mixed        mixed layer ocean heat anomaly [10^22 J]
#   ocheat.interior     interior ocean heat anomaly [10^22 J]
#                       (Note: ocheat = ocheat.mixed + ocheat.interior)
#   ocheatflux.mixed    heat uptake of the ocean mixed layer [W/m2]
#   ocheatflux.interior heat uptake of the ocean interior [W/m2]
#
# -- GSIC-MAGICC --
#   sl_gsic             cumulative sea-level contribution since t0 (i=1)a [m]
#                       by definition Gs(1) = 0
#
# -- GIS-SIMPLE --
#   sl_gis              Volume of Greenland ice sheet [m sea-level equivalent (sle)]
#
# -- TE --
#   sl_te               cumulative sea-level contribution from thermal expansion [m]
#
# -- DAIS --
#   rad_ais             Antarctic ice sheet radius [m]
#   V_ais               Volume of Antarctic ice sheet [m3]
#   Ta                  Antarctic mean surface temperature [degC] (or global temperature)
#   Toc                 High latitude ocean subsurface temperatures [degC] (calculated by ANTO)
#
#===============================================================================
# Internal variables and constants:
#   TEeq                Equilibrium TE [m]
#   b                   Undisturbed bed profile [m]
#   h                   Ice sheet surface height [m]
#   Btot                Total mass accumulation rate on the ice-sheet [ ??? ]
#   F                   Total ice flux across the grounding line [m3/yr]
#   rc                  Distance from the continent center to where the ice sheets enters the sea [m]
#   hr                  Height of runoff line above which precipitation accumulates as snow [m]
#   P                   Precipitation [m (ice equivalent)]
#   beta                Mass balance gradient [ m-0.5 ]
#   rR                  Distance from the continent center to where the runoff line intersects the ice sheet surface [m]
#   Hw                  Water depth at grounding line
#   r                   Radial coordinate [m]
#   Tf                  Freezing temperature sea water [degC]
#   rho_w               (Sea) water density [kg/m3]
#   rho_i               Ice density [kg/m3]
#   rho_m               Rock density [kg/m3]
#===============================================================================


#===============================================================================
# convert annual ocean heat flux (W/m^2) to cumulative ocean heat content anomaly (10^22 J)
flux.to.heat = function(heatflux.mixed, heatflux.interior)
{
	flnd = 0.29 # area land fraction
	fso = 0.95 # ocean area fraction of interior
	sec.per.year = 31556926
	earth.area = 510065600 * 10^6
	ocean.area = (1-flnd)*earth.area
	powtoheat = ocean.area*sec.per.year / 10^22 # in 10^22 J/yr

	heat.mixed = cumsum(heatflux.mixed) * powtoheat
	heat.interior = fso * cumsum(heatflux.interior) * powtoheat
	ocean.heat = heat.mixed + heat.interior

	return(list(ocean.heat=ocean.heat, heat.mixed=heat.mixed, heat.interior=heat.interior))
}
#===============================================================================


#===============================================================================
# load fortran subroutine
# to check if library is loaded: is.loaded("run_brick") (for example)
# dyn.load("../fortran/brick_te.so")
if(.Platform$OS.type == "unix") {
    #dyn.load("../fortran/brick.so")
    dyn.load("../fortran/dais.so")
    dyn.load("../fortran/doeclim.so")
    dyn.load("../fortran/gsic_magicc.so")
    dyn.load("../fortran/simple.so")
    dyn.load("../fortran/brick_te.so")
} else {
    #dyn.load("../fortran/brick")
    dyn.load("../fortran/dais")
    dyn.load("../fortran/doeclim")
    dyn.load("../fortran/gsic_magicc")
    dyn.load("../fortran/simple")
    dyn.load("../fortran/brick_te")
}

brickF <- function(
	tstep,
    mod.time,
    forcing.total,
    S.doeclim = 3.1,
    kappa.doeclim = 3.5,
    beta0.gsic = 0.000577,
    V0.gsic = 0.4,
    n.gsic = 0.82,
    Gs0.gsic = 0,
    Teq.gsic = -0.15,
    a.simple = -0.827,
    b.simple = 7.242,
    alpha.simple = 1.630e-4,
    beta.simple = 2.845e-05,
    V0.simple = 7.242,
    a.te = 0.5,
    b.te = 0,
    invtau.te = 0.005,
    V0.te = 0,
    a.anto = 0.26,
    b.anto = 0.62,
    slope.Ta2Tg = 1,
    intercept.Ta2Tg = 0,
    b0.dais = 775,
    slope.dais = 6 * 10^(-4),
    mu.dais = 8.7,
    h0.dais = 1471,
    c.dais = 95,
    P0.dais = 0.35,
    kappa.dais = 4 * 10^(-2),
    nu.dais = 1.2 * 10^(-2),
    f0.dais = 1.2,
    gamma.dais = 2.5,
    alpha.dais = 0.5,
    Tf.dais = -1.8,
    rho_w.dais = 1030,
    rho_i.dais = 917,
    rho_m.dais = 4000,
    Toc_0.dais = 0.72,
    Rad0.dais = 1.864 * 10^6,
    Aoc.dais = 3.619e14,
    lf = -1.18,
    includes_dSLais = 0
){

  # determine series length
  ns <- length(forcing.total)

  # wrap up the DAIS parameters
  # TODO - bundle the other models' parameter like this?
  parameters.dais <- c( b0.dais,
                        slope.dais,
                        mu.dais,
                        h0.dais,
                        c.dais,
                        P0.dais,
                        kappa.dais,
                        nu.dais,
                        f0.dais,
                        gamma.dais,
                        alpha.dais,
                        Tf.dais,
                        rho_w.dais,
                        rho_i.dais,
                        rho_m.dais,
                        Toc_0.dais,
                        Rad0.dais,
                        Aoc.dais,
                        lf,
                        includes_dSLais)

# TODO - HERE NOW
# putting in the parameters - need to do DAIS
# might also want to bundle all of the models' parameters like DAIS so it
# is easier to swap them in and out
# TODO - HERE NOW

  # call fortran
  f.output <- .Fortran("run_brick",
        ns = ns,
        tstep = as.double(tstep),
        forcing_in = as.double(forcing.total),
        doeclim_t2co = as.double(S.doeclim),
        doeclim_kappa = as.double(kappa.doeclim),
        gsic_magicc_beta0 = as.double(beta0.gsic),
        gsic_magicc_V0 = as.double(V0.gsic),
        gsic_magicc_n = as.double(n.gsic),
        gsic_magicc_Gs0 = as.double(Gs0.gsic),
        gsic_magicc_Teq = as.double(Teq.gsic),
        simple_a = as.double(a.simple),
        simple_b = as.double(b.simple),
        simple_alpha = as.double(alpha.simple),
        simple_beta = as.double(beta.simple),
        simple_V0 = as.double(V0.simple),
        brick_te_a = as.double(a.te),
        brick_te_b = as.double(b.te),
        brick_te_invtau = as.double(invtau.te),
        brick_te_V0 = as.double(V0.te),
        anto_a = as.double(a.anto),
        anto_b = as.double(b.anto),
        slope_Ta2Tg = as.double(slope.Ta2Tg),
        intercept_Ta2Tg = as.double(intercept.Ta2Tg),
        dais_parameters = as.double(parameters.dais),
        time_out = as.double(mod.time),
        temp_out = as.double(rep(0,ns)),
        heatflux_mixed_out = as.double(rep(0,ns)),
        heatflux_interior_out = as.double(rep(0,ns)),
        sl_te_out = as.double(rep(-999.99,ns)),
        sl_gsic_out = as.double(rep(-999.99,ns)),
        sl_gis_out = as.double(rep(-999.99,ns)),
        GIS_Volume_out = as.double(rep(-999.99,ns)),
        sl_ais_out = as.double(rep(-999.99,ns)),
        AIS_Radius_out = as.double(rep(-999.99,ns)),
        AIS_Volume_out = as.double(rep(-999.99,ns)),
		sl_out = as.double(rep(-999.99,ns))

    )

# TODO - HERE NOW
# calculate the right outputs internally (in the fortran routines)
# and externally (here)
# TODO - HERE NOW

    ocheat = flux.to.heat(f.output$heatflux_mixed, f.output$heatflux_interior)

	model.output = list(time=mod.time, temp=f.output$temp_out, ocheat=ocheat$ocean.heat,
                        ocheat.mixed=ocheat$heat.mixed, ocheat.interior=ocheat$heat.interior,
						ocheatflux.mixed = f.output$heatflux_mixed, ocheatflux.interior = f.output$heatflux_interior)

    return(f.output$brick_out)

}
#===============================================================================
# End
#===============================================================================
