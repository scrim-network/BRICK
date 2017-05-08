for (i in 1:4) {

	setwd('~/codes/BRICK/calibration')
	library(ncdf4)
    today=Sys.Date(); today=format(today,format="%d%b%Y")
	if(i==1) {
        print('processing RCP85, gamma priors...')
        filename.in = "../output_model/BRICK-fastdyn_physical_gamma_07May2017.nc"
	    ncdata <- nc_open(filename.in)
	    sea_level = ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
        sea_level_nofd = ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
        mod.time =ncvar_get(ncdata, 'time_proj')
        nc_close(ncdata)
        filename.vdout = paste('../output_model/vanDantzig_RCP85_gamma_',today,'.nc',sep="")
    } else if(i==2) {
        print('processing RCP45, gamma priors...')
        filename.in = "../output_model/BRICK-fastdyn_physical_gamma_07May2017.nc"
	    ncdata <- nc_open(filename.in)
	    sea_level = ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
        sea_level_nofd = ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP45')
        mod.time =ncvar_get(ncdata, 'time_proj')
        nc_close(ncdata)
        filename.vdout = paste('../output_model/vanDantzig_RCP45_gamma_',today,'.nc',sep="")
    } else if(i==3) {
        print('processing RCP26, gamma priors...')
        filename.in = "../output_model/BRICK-fastdyn_physical_gamma_07May2017.nc"
	    ncdata <- nc_open(filename.in)
	    sea_level = ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
        sea_level_nofd = ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26')
        mod.time =ncvar_get(ncdata, 'time_proj')
        nc_close(ncdata)
        filename.vdout = paste('../output_model/vanDantzig_RCP26_gamma_',today,'.nc',sep="")
    } else if(i==4) {
        print('processing RCP85, uniform priors...')
        filename.in = "../output_model/BRICK-fastdyn_physical_uniform_07May2017.nc"
	    ncdata <- nc_open(filename.in)
	    sea_level = ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
        sea_level_nofd = ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
        mod.time =ncvar_get(ncdata, 'time_proj')
        nc_close(ncdata)
        filename.vdout = paste('../output_model/vanDantzig_RCP85_uniform_',today,'.nc',sep="")
    } else {
        print('ERROR IN RUN INDEX')
    }

## Evaluate van Dantzig model for each of these realizations

## Set up
currentyear		= 2015		# start year
endyear			= 2100		# reevaluation year
lowleveeheight	= 0			# lowest levee heightening considered (m)
highleveeheight	= 10		# largest levee heightening considered (m)
increment		= 0.05		# increment of levee heightening range (m)

## Time horizon until next evaluation of levees (years)
T = endyear - currentyear
time_frame=mod.time

## Range of considered levee heights (meters)
X = seq(from = lowleveeheight, to = highleveeheight, by = increment)

## Create a matrix of sea-level data starting from the current year for each model simulation.
## Note that if you included the NOLA fingerprinting above (in by default), then
## what is called "global_sea_level" below is actually local, but does not include
## subsidence.
global_sea_level = data.matrix(sea_level[match(currentyear, time_frame):match(endyear, time_frame),])
global_sea_level_nofd = data.matrix(sea_level_nofd[match(currentyear, time_frame):match(endyear, time_frame),])

## Set up basic information for the van Dantzig model.

## Number of observations (ensemble members):
n_obs = length(global_sea_level[1, ])

## Conversion from feet to meters:
feet_to_meter = function(ft) { ft * 0.3048 }

## Investment costs for levee heightening. Table # 2 in Jonkman et al. (2009)
design_surge_level_ft = c(9, 11, 13, 15, 17, 21, 25)
design_surge_level_M = feet_to_meter(design_surge_level_ft)

investments_US = c(2.2E09, 2.4E09, 2.6E09, 2.9E09, 3.1E09, 3.6E09, 4.1E09)

I_table = matrix(c(design_surge_level_M, investments_US), ncol=2)
colnames(I_table) = c("Levee heightening", "Associated cost")

p_zero_p = 0.0038                 # Initial flood frequency (1/yr) with zero height increase (Van Dantzig (1956))
alpha_p = 2.6                     # Exponential flood frequency constant (Van Dantzig (1956))
V_p = c(5e+9, 3e+10)              # Value of goods protected by dike (based on estimates in Jonkman et al. (2009)) (US$)
delta_prime_p = c(0.02, 0.06)     # Discount rate (percent/year) (based on estimates in Jonkman et al. (2009))
I_range = c(-0.5,1.0)             # Investment cost uncertainty range (as fraction of investment cost) (Jonkman and Dutch Perspective use -50%, +100%)
subs_rate = 0.0056                # Rate of land subsidence (meter/year) (Dixon et al. (2006))

investment.fit <- lm(I_table[,2]~I_table[,1])
intercept.h2i = investment.fit$coefficients[[1]]
slope.h2i = investment.fit$coefficients[[2]]

source("../R/VD_NOLA.R")

## Sample parameters (using Dutch Perspective as a guide)
## Uses Latin Hypercube to sample Van Dantzig parameters, assigns each ensemble
## member of SLR a set of parameters.
params.vd = parameter_sampling_DP(n_obs, p_zero_p, alpha_p, V_p, delta_prime_p, I_range, subs_rate)

## Test simulation to get the names
vanDantzig.out = VD_NOLA_R(params.vd[1,], I_table, T, X, global_sea_level[,1])

n.ensemble = nrow(params.vd)
vanDantzig.ensemble = vector("list", ncol(vanDantzig.out))
vanDantzig.nofd.ensemble = vector("list", ncol(vanDantzig.out))
names(vanDantzig.ensemble) = names(vanDantzig.out)
names(vanDantzig.nofd.ensemble) = names(vanDantzig.out)
for (i in 1:ncol(vanDantzig.out)) {
  vanDantzig.ensemble[[i]] = mat.or.vec(nrow(vanDantzig.out), n.ensemble)
  vanDantzig.nofd.ensemble[[i]] = mat.or.vec(nrow(vanDantzig.out), n.ensemble)
}

## Run the simulations
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {
  vanDantzig.out = VD_NOLA_R(params.vd[i,], I_table, T, X, global_sea_level[,i])
  for (j in 1:ncol(vanDantzig.out)) {
		vanDantzig.ensemble[[j]][,i] = vanDantzig.out[,j]
	}
  vanDantzig.out = VD_NOLA_R(params.vd[i,], I_table, T, X, global_sea_level_nofd[,i])
  for (j in 1:ncol(vanDantzig.out)) {
		vanDantzig.nofd.ensemble[[j]][,i] = vanDantzig.out[,j]
	}
  setTxtProgressBar(pb, i)
}
close(pb)

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Output results of van Dantzig to netCDF file

library(ncdf4)

dim.heightening <- ncdim_def('H', 'meters', as.double(X))
dim.ensemble <- ncdim_def('ens', 'ensemble member', (1:n.ensemble))

cost <- ncvar_def('ExpectedCost', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost accounting for fast dynamics')
cost.nofd <- ncvar_def('ExpectedCost_nofd', 'US$', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected cost without accounting for fast dynamics')
loss <- ncvar_def('ExpectedLoss', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss accounting for fast dynamics')
loss.nofd <- ncvar_def('ExpectedLoss_nofd', 'US$', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected loss without accounting for fast dynamics')
investment <- ncvar_def('ExpectedInvestment', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment accounting for fast dynamics')
investment.nofd <- ncvar_def('ExpectedInvestment_nofd', 'US$', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected investment without accounting for fast dynamics')
preturn <- ncvar_def('ExpectedPreturn', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period accounting for fast dynamics')
preturn.nofd <- ncvar_def('ExpectedPreturn_nofd', 'years', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected return period without accounting for fast dynamics')

today=Sys.Date(); today=format(today,format="%d%b%Y")

outnc <- nc_create(filename.vdout,
										list(cost, cost.nofd, loss, loss.nofd, investment, investment.nofd, preturn, preturn.nofd),
										force_v4 = TRUE)

ncvar_put(outnc, cost, vanDantzig.ensemble$Expected_costs)
ncvar_put(outnc, cost.nofd, vanDantzig.nofd.ensemble$Expected_costs)
ncvar_put(outnc, loss, vanDantzig.ensemble$Expected_loss)
ncvar_put(outnc, loss.nofd, vanDantzig.nofd.ensemble$Expected_loss)
ncvar_put(outnc, investment, vanDantzig.ensemble$Investment)
ncvar_put(outnc, investment.nofd, vanDantzig.nofd.ensemble$Investment)
ncvar_put(outnc, preturn, 1/vanDantzig.ensemble$Average_p_fail)
ncvar_put(outnc, preturn.nofd, 1/vanDantzig.nofd.ensemble$Average_p_fail)

nc_close(outnc)

}
