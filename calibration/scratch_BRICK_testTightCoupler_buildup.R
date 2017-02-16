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
n.ensemble <- 10

# read calibrated model simulations
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

#TODO

# draw ensemble parameters
filename.parameters <- '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
ncdata <- nc_open(filename.parameters)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters <- t(parameters)
colnames(parameters) <- parnames

#TODO

# use for initial condition the 1850 hindcast (set up index i0 for this)
#TODO

# call the tightly coupled model
mod.time <- t.hind
l.project = FALSE

begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

source('../R/compute_indices.R')				# function to determine the model and
source('../calibration/DOECLIM_readData.R')		# read DOECLIM calibration data
source('../calibration/GSIC_readData.R')			# read GSIC calibration data
source('../calibration/TE_readData.R')        # read TE data
source('../calibration/SIMPLE_readData.R')		# GIS data, and trends in mass balance
source('../calibration/DAIS_readData.R')			# DAIS forcing data (if at all uncoupled)

# initialize matrix to store model ensemble output
brick.out <- vector("list", n.ensemble)
source('../R/forcing_total.R')					# function to add up the total forcing
forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )

# source the models
source('../R/BRICK_coupledModel.R')
source('../fortran/R/doeclimF.R')       # the DOECLIM model
source('../fortran/R/GSIC_magiccF.R')   # the GSIC model
source('../fortran/R/brick_te_F.R')     # TE (thermosteric expansion) model
source('../fortran/R/simpleF.R')        # GIS (Greenland Ice Sheet) model
source('../fortran/R/daisantoF.R')			# DAIS (Antarctic Ice Sheet) model

# set up the model
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermosteric expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.mymodel  = FALSE   # Example of adding your own model component
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)
source('../calibration/BRICK_parameterSetup.R')

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Set the indices of the initial condition for each sub-model
i0 = vector("list",nrow(ind.norm.data)); names(i0)=as.character(ind.norm.data[,1])

## GSIC initial conditions are actually relative to 1990 (Wigley and Raper 2005)
## Re-set these. The simulation is relative to 1990, but results and comparison
## to data is relative to 1960.
i0$gsic = which(mod.time==1990)

## GIS initial conditions are relative to 1961-1990
i0$gis = which(mod.time==1961)



#debug?
parameters.in <- as.numeric(parameters[1,])
parnames.in <- parnames
forcing.in <- forcing
slope.Ta2Tg.in <- slope.Ta2Tg
intercept.Ta2Tg.in <- intercept.Ta2Tg
ind.norm.sl <- ind.norm



## Run the sample, and enjoy a nice progress bar
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

	brick.out[[i]] <- brick_model(parameters.in			= as.numeric(parameters[i,]),
																parnames.in				= parnames,
																forcing.in				= forcing,
																l.project					= l.project,
																slope.Ta2Tg.in		= slope.Ta2Tg,
																intercept.Ta2Tg.in= intercept.Ta2Tg,
																mod.time					= mod.time,
																ind.norm.data 		= ind.norm.data,
																ind.norm.sl 			= ind.norm,
																luse.brick 				= luse.brick,
																i0=i0
																)

  setTxtProgressBar(pb, i)
}
close(pb)

iH0 = match("H0",parnames)
iT0 = match("T0",parnames)
iGs0= match("Gs0",parnames)

par(mfrow=c(2,2))

plot(mod.time, brick.out[[1]]$doeclim.out$temp+parameters[1,iT0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$doeclim.out$temp+parameters[i,iT0], type='l', col='black')}
points(obs.temp.time, obs.temp, pch=16, col='red')

plot(mod.time, brick.out[[1]]$doeclim.out$ocheat+parameters[1,iH0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$doeclim.out$ocheat+parameters[i,iH0], type='l', col='black')}
points(obs.ocheat.time, obs.ocheat, pch=16, col='red')

plot(mod.time, brick.out[[1]]$gsic.out+parameters[1,iGs0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$gsic.out+parameters[i,iGs0], type='l', col='black')}
points(obs.gsic.time, obs.gsic, pch=16, col='red')

plot(mod.time, brick.out[[1]]$simple.out$sle.gis, type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$simple.out$sle.gis, type='l', col='black')}
points(obs.gis.time, obs.gis, pch=16, col='red')

##==============================================================================
## Okay, it all works. Now screw around with the "fully-forward" BRICK model.
## Want to code this as a two-step (run_BRICK_forward and BRICK_step_forward)
## R-calling-into-Fortran model, similar to the other BRICK sub-models.
##==============================================================================

##==============================================================================
## Test running just hte first ensemble member in "fully-forward" mode versus
## the "old" mode. Trouble-shooting time.
##==============================================================================

# Step 1:   want to be able to run the model in "old" mode, with the original
#           i0 for each submodel, then use as initial conditions the 1850 values
#           (which might be tricky for GIS and AIS, since they are in terms of
#           ice sheet volume), and hopefully get the same result with i0=1.


# Step 2:   get the same result with the "fully-forward"/"step-together" version
#...

source('../fortran/R/brickF.R')

tstep <- 1
parameters.in <- parameters[1,]
parnames.in <- names(parameters.in)
alpha.doeclim <- parameters.in[match("alpha.doeclim",parnames.in)]
forcing.total <- forcing_total(forcing=forcing,
							  alpha.doeclim=alpha.doeclim,
                              l.project=l.project,
							  begyear=mod.time[1],
							  endyear=mod.time[length(mod.time)]
							  )

dyn.load("../fortran/brick.so")
dyn.load("../fortran/dais.so")
dyn.load("../fortran/doeclim.so")
dyn.load("../fortran/gsic_magicc.so")
dyn.load("../fortran/simple.so")
dyn.load("../fortran/brick_te.so")

S.doeclim = 3.1;
kappa.doeclim = 3.5;
T0.doeclim = 0;
beta0.gsic = 0.000577;
V0.gsic = 0.4;
n.gsic = 0.82;
Gs0.gsic = 0;
Teq.gsic = -0.15;
a.simple = -0.827;
b.simple = 7.242;
alpha.simple = 1.630e-4;
beta.simple = 2.845e-05;
V0.simple = 7.242;
a.te = 0.5;
b.te = 0;
invtau.te = 0.005;
V0.te = 0;
a.anto = 0.26;
b.anto = 0.62;
slope.Ta2Tg = 0.8364527;
intercept.Ta2Tg = 15.4235;
b0.dais = 775;
slope.dais = 6 * 10^(-4);
mu.dais = 8.7;
h0.dais = 1471;
c.dais = 95;
P0.dais = 0.35;
kappa.dais = 4 * 10^(-2);
nu.dais = 1.2 * 10^(-2);
f0.dais = 1.2;
gamma.dais = 2.5;
alpha.dais = 0.5;
Tf.dais = -1.8;
rho_w.dais = 1030;
rho_i.dais = 917;
rho_m.dais = 4000;
Toc_0.dais = 0.72;
Rad0.dais = 1.864 * 10^6;
Aoc.dais = 3.619e14;
lf = -1.18;
includes_dSLais = 0;
parameters.dais <- c( b0.dais, slope.dais, mu.dais,
                      h0.dais, c.dais, P0.dais,
                      kappa.dais, nu.dais, f0.dais,
                      gamma.dais, alpha.dais, Tf.dais,
                      rho_w.dais, rho_i.dais, rho_m.dais,
                      Toc_0.dais, Rad0.dais, Aoc.dais,
                      lf, includes_dSLais)

ns <- length(forcing.total)

f.output <- .Fortran("run_brick",
      ns = ns,
      tstep = as.double(tstep),
      forcing_in = as.double(forcing.total),
      doeclim_t2co = as.double(S.doeclim),
      doeclim_kappa = as.double(kappa.doeclim),
      doeclim_T0 = as.double(T0.doeclim),
      time_out = as.double(mod.time),
      temp_out = as.double(rep(0,ns)),
      heatflux_mixed_out = as.double(rep(0,ns)),
      heatflux_interior_out = as.double(rep(0,ns)),
      gsic_magicc_beta0 = as.double(beta0.gsic),
      gsic_magicc_V0 = as.double(V0.gsic),
      gsic_magicc_n = as.double(n.gsic),
      gsic_magicc_Gs0 = as.double(Gs0.gsic),
      gsic_magicc_Teq = as.double(Teq.gsic),
      sl_gsic_out = as.double(rep(-999.99,ns)),
      brick_te_a = as.double(a.te),
      brick_te_b = as.double(b.te),
      brick_te_invtau = as.double(invtau.te),
      brick_te_V0 = as.double(V0.te),
      sl_te_out = as.double(rep(-999.99,ns)),
      simple_a = as.double(a.simple),
      simple_b = as.double(b.simple),
      simple_alpha = as.double(alpha.simple),
      simple_beta = as.double(beta.simple),
      simple_V0 = as.double(V0.simple),
      sl_gis_out = as.double(rep(-999.99,ns)),
      vol_gis_out = as.double(rep(-999.99,ns)),
      anto_a = as.double(a.anto),
      anto_b = as.double(b.anto),
      slope_Ta2Tg = as.double(slope.Ta2Tg),
      intercept_Ta2Tg = as.double(intercept.Ta2Tg),
      dais_parameters = as.double(parameters.dais),
      sl_ais_out = as.double(rep(-999.99,ns)),
      rad_ais_out = as.double(rep(-999.99,ns)),
      vol_ais_out = as.double(rep(-999.99,ns)),
      sl_out = as.double(rep(-999.99,ns))
  )


## Test against the models run sequentially:
doeclim.output <- .Fortran( "run_doeclim",
        ns = ns,
        time_out = as.double(mod.time),
        forcing_in = as.double(forcing.total),
        t2co_in = as.double(S.doeclim),
        kappa_in = as.double(kappa.doeclim),
        temp_out = as.double(rep(0,ns)),
        heatflux_mixed_out = as.double(rep(0,ns)),
        heatflux_interior_out = as.double(rep(0,ns))
        )

gsic.output <- .Fortran("run_gsic_magicc",
                ns                = ns,
                tstep             = as.double(tstep),
                gsic_magicc_beta0 = as.double(beta0.gsic),
                gsic_magicc_V0    = as.double(V0.gsic),
                gsic_magicc_n     = as.double(n.gsic),
                gsic_magicc_Gs0   = as.double(Gs0.gsic),
                gsic_magicc_Teq   = as.double(Teq.gsic),
                Gl_Temp           = as.double(doeclim.output$temp_out),
                gsic_magicc_i0    = as.double(1),
                SL_contribution_out = as.double(rep(-999.99,ns))
                )

gis.output <- .Fortran("run_simple",
                ns            = ns,
                tstep         = as.double(tstep),
                simple_a      = as.double(a.simple),
                simple_b      = as.double(b.simple),
                simple_alpha  = as.double(alpha.simple),
                simple_beta   = as.double(beta.simple),
                simple_V0     = as.double(V0.simple),
                Grl_Temp      = as.double(doeclim.output$temp_out),
                simple_i0     = as.double(1),
                GIS_Volume_out = as.double(rep(-999.99,ns))
)
sle.gis <- V0.simple - gis.output$GIS_Volume_out


te.output <- .Fortran("run_brick_te",
                  ns            = ns,
                  tstep         = as.double(tstep),
                  brick_te_a    = as.double(a.te),
                  brick_te_b    = as.double(b.te),
                  brick_te_invtau = as.double(invtau.te),
                  brick_te_TE_0 = as.double(V0.te),
                  Gl_Temp       = as.double(doeclim.output$temp_out),
                  brick_te_i0   = as.double(1),
                  TE_out        = as.double(rep(-999.99,ns))
)


SL <- f.output$sl_ais_out + sle.gis + te.output$TE_out + gsic.output$SL_contribution_out
for (i in 2:length(SL)) {SL[i]=f.output$sl_ais_out[i-1] + sle.gis[i-1] +
                                te.output$TE_out[i-1] + gsic.output$SL_contribution_out[i-1] +
                                        1.1*diff(sle.gis[(i-1):i]) +
                                        1.0*diff(te.output$TE_out[(i-1):i]) +
                                        1.1*diff(gsic.output$SL_contribution_out[(i-1):i])}
dSL <- c(0, 1.1*diff(sle.gis) +
        1.0*diff(te.output$TE_out) +
        1.1*diff(gsic.output$SL_contribution_out) )
Toc <- anto(a=a.anto, b=b.anto, Tf=Tf.dais, Tg=doeclim.output$temp_out)
Ta.recon = (doeclim.output$temp_out-intercept.Ta2Tg)/slope.Ta2Tg
dais.output <- .Fortran("run_dais",
                ns                 = ns,
                tstep              = as.double(tstep),
                dais_parameters    = as.double(parameters.dais),
                Ant_Temp           = as.double(Ta.recon),
                Ant_Sea_Level      = as.double(SL),
                Ant_Sur_Ocean_Temp = as.double(Toc),
                Ant_SL_rate        = as.double(dSL),
                AIS_Radius_out     = as.double(rep(-999.99,ns)),
                AIS_Volume_out     = as.double(rep(-999.99,ns))
)
Vsle = 57*(1-dais.output$AIS_Volume_out/dais.output$AIS_Volume_out[1]) #Takes steady state present day volume to correspond to 57m SLE


tmp.out <- brickF(  tstep=tstep,
                    mod.time=mod.time,
                    forcing.total = forcing.total,
                    S.doeclim = parameters.in[match("S",parnames.in)],
                    kappa.doeclim = parameters.in[match("kappa.doeclim",parnames.in)],
                    T0.doeclim = parameters.in[match("T0",parnames.in)],
                    H0.doeclim = parameters.in[match("H0",parnames.in)],
                    beta0.gsic = parameters.in[match("beta0",parnames.in)],
                    V0.gsic = parameters.in[match("V0.gsic",parnames.in)],
                    n.gsic = parameters.in[match("n",parnames.in)],
                    Gs0.gsic = parameters.in[match("Gs0",parnames.in)],
                    a.simple = parameters.in[match("a.simple",parnames.in)],
                    b.simple = parameters.in[match("b.simple",parnames.in)],
                    alpha.simple = parameters.in[match("alpha.simple",parnames.in)],
                    beta.simple = parameters.in[match("beta.simple",parnames.in)],
                    V0.simple = parameters.in[match("V0",parnames.in)],
                    a.te = parameters.in[match("a.te",parnames.in)],
                    b.te = parameters.in[match("b.te",parnames.in)],
                    invtau.te = parameters.in[match("invtau.te",parnames.in)],
                    V0.te = parameters.in[match("TE0",parnames.in)],
                    a.anto = parameters.in[match("anto.a",parnames.in)],
                    b.anto = parameters.in[match("anto.b",parnames.in)],
                    slope.Ta2Tg = slope.Ta2Tg,
                    intercept.Ta2Tg = intercept.Ta2Tg,
                    b0.dais = parameters.in[match("b0",parnames.in)],
                    slope.dais = parameters.in[match("slope",parnames.in)],
                    mu.dais = parameters.in[match("mu",parnames.in)],
                    h0.dais = parameters.in[match("h0",parnames.in)],
                    c.dais = parameters.in[match("c",parnames.in)],
                    P0.dais = parameters.in[match("P0",parnames.in)],
                    kappa.dais = parameters.in[match("kappa.dais",parnames.in)],
                    nu.dais = parameters.in[match("nu",parnames.in)],
                    f0.dais = parameters.in[match("f0",parnames.in)],
                    gamma.dais = parameters.in[match("gamma",parnames.in)],
                    alpha.dais = parameters.in[match("alpha.dais",parnames.in)]
                    )

# NOTE: the above run will produce wacky results
# this is because it is applying ~1960 initial conditions to 1850

# NOTE: below, use the i=1 run, but find suitable initial conditions.

##==============================================================================

## Reproduce an old run, but start from 1850
i <- 1

Gs0.gsic <- brick.out[[i]]$gsic.out[1]
V0.gsic <- parameters.in[match("V0.gsic",parnames.in)] +
            (parameters.in[match("Gs0",parnames.in)] - brick.out[[i]]$gsic.out[1])
V0.te <- brick.out[[i]]$te.out[1]
V0.simple <- brick.out[[i]]$simple.out$Vgrl[1]

tmp.out <- brickF(  tstep=tstep,
                    mod.time=mod.time,
                    forcing.total = forcing.total,
                    S.doeclim = parameters.in[match("S",parnames.in)],
                    kappa.doeclim = parameters.in[match("kappa.doeclim",parnames.in)],
                    T0.doeclim = parameters.in[match("T0",parnames.in)],
                    H0.doeclim = parameters.in[match("H0",parnames.in)],
                    beta0.gsic = parameters.in[match("beta0",parnames.in)],
                    V0.gsic = V0.gsic,
                    n.gsic = parameters.in[match("n",parnames.in)],
                    Gs0.gsic = Gs0.gsic,
                    a.simple = parameters.in[match("a.simple",parnames.in)],
                    b.simple = parameters.in[match("b.simple",parnames.in)],
                    alpha.simple = parameters.in[match("alpha.simple",parnames.in)],
                    beta.simple = parameters.in[match("beta.simple",parnames.in)],
                    V0.simple = V0.simple,
                    a.te = parameters.in[match("a.te",parnames.in)],
                    b.te = parameters.in[match("b.te",parnames.in)],
                    invtau.te = parameters.in[match("invtau.te",parnames.in)],
                    V0.te = V0.te,
                    a.anto = parameters.in[match("anto.a",parnames.in)],
                    b.anto = parameters.in[match("anto.b",parnames.in)],
                    slope.Ta2Tg = slope.Ta2Tg,
                    intercept.Ta2Tg = intercept.Ta2Tg,
                    b0.dais = parameters.in[match("b0",parnames.in)],
                    slope.dais = parameters.in[match("slope",parnames.in)],
                    mu.dais = parameters.in[match("mu",parnames.in)],
                    h0.dais = parameters.in[match("h0",parnames.in)],
                    c.dais = parameters.in[match("c",parnames.in)],
                    P0.dais = parameters.in[match("P0",parnames.in)],
                    kappa.dais = parameters.in[match("kappa.dais",parnames.in)],
                    nu.dais = parameters.in[match("nu",parnames.in)],
                    f0.dais = parameters.in[match("f0",parnames.in)],
                    gamma.dais = parameters.in[match("gamma",parnames.in)],
                    alpha.dais = parameters.in[match("alpha.dais",parnames.in)]
                    )

##==============================================================================
## Timing test, old vs new

i <- 1

Gs0.gsic <- brick.out[[i]]$gsic.out[1]
V0.gsic <- parameters.in[match("V0.gsic",parnames.in)] +
            (parameters.in[match("Gs0",parnames.in)] - brick.out[[i]]$gsic.out[1])
V0.te <- brick.out[[i]]$te.out[1]
V0.simple <- brick.out[[i]]$simple.out$Vgrl[1]

niter <- 1e5

t0.old <- proc.time()
for (i in 1:niter){
    # call old brick
    old.out <- brick_model(parameters.in			= as.numeric(parameters.in),
    		  			   parnames.in				= parnames.in,
    					   forcing.in				= forcing,
    					l.project					= l.project,
    							slope.Ta2Tg.in		= slope.Ta2Tg,
    							intercept.Ta2Tg.in= intercept.Ta2Tg,
    							mod.time					= mod.time,
    							ind.norm.data 		= ind.norm.data,
    							ind.norm.sl 			= ind.norm,
    							luse.brick 				= luse.brick,
    							i0=i0
    							)

}
t1.old <- proc.time()

t0.new <- proc.time()
for (i in 1:niter){
    # call new brick
    tmp.out <- brickF(  tstep=tstep,
                        mod.time=mod.time,
                        forcing.total = forcing.total,
                        S.doeclim = parameters.in[match("S",parnames.in)],
                        kappa.doeclim = parameters.in[match("kappa.doeclim",parnames.in)],
                        T0.doeclim = parameters.in[match("T0",parnames.in)],
                        H0.doeclim = parameters.in[match("H0",parnames.in)],
                        beta0.gsic = parameters.in[match("beta0",parnames.in)],
                        V0.gsic = V0.gsic,
                        n.gsic = parameters.in[match("n",parnames.in)],
                        Gs0.gsic = Gs0.gsic,
                        a.simple = parameters.in[match("a.simple",parnames.in)],
                        b.simple = parameters.in[match("b.simple",parnames.in)],
                        alpha.simple = parameters.in[match("alpha.simple",parnames.in)],
                        beta.simple = parameters.in[match("beta.simple",parnames.in)],
                        V0.simple = V0.simple,
                        a.te = parameters.in[match("a.te",parnames.in)],
                        b.te = parameters.in[match("b.te",parnames.in)],
                        invtau.te = parameters.in[match("invtau.te",parnames.in)],
                        V0.te = V0.te,
                        a.anto = parameters.in[match("anto.a",parnames.in)],
                        b.anto = parameters.in[match("anto.b",parnames.in)],
                        slope.Ta2Tg = slope.Ta2Tg,
                        intercept.Ta2Tg = intercept.Ta2Tg,
                        b0.dais = parameters.in[match("b0",parnames.in)],
                        slope.dais = parameters.in[match("slope",parnames.in)],
                        mu.dais = parameters.in[match("mu",parnames.in)],
                        h0.dais = parameters.in[match("h0",parnames.in)],
                        c.dais = parameters.in[match("c",parnames.in)],
                        P0.dais = parameters.in[match("P0",parnames.in)],
                        kappa.dais = parameters.in[match("kappa.dais",parnames.in)],
                        nu.dais = parameters.in[match("nu",parnames.in)],
                        f0.dais = parameters.in[match("f0",parnames.in)],
                        gamma.dais = parameters.in[match("gamma",parnames.in)],
                        alpha.dais = parameters.in[match("alpha.dais",parnames.in)]
                        )
}
t1.new <- proc.time()


##==============================================================================


##==============================================================================
## End
##==============================================================================
