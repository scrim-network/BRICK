# =======================================================================================
# Driver script for the flood risk/Van Dantzig analysis. This module serves as
# a coupling link between the physical and socioeconomic impacts modules of
# BRICK.
#
#
#
# Questions? Tony Wong (twong@psu.edu)
# =======================================================================================
#
#   Requires (input variables):
# - sea_level     simulated local sea level [m sle]
# - time          times of model steps [years]
#
#   Simulates (output variables):
# - vanDantzig.ensemble   ensemble of Van Dantzig outputs (heightenings and
#                         corresponding investment costs, losses, and failure
#                         probabilities)
#
#   Parameters:
# - currentyear           year in which we initially evaluate the dikes
# - endyear               year the dikes must hold strong until
# - lowleveeheight        lowest heightening considered [m]
# - highleveeheight       highest dike heightening considered [m]
# - increment             consider dike heightenings from lowleveeheight to
#                         highleveeheight in increments of 'increment'
# - p_zero_p              initial flood frequency (1/yr) if no heightening
# - alpha_p               exponential flood frequency constant
# - V_p                   value of goods protected by the dike
# - delta_prime_p         net discount rate
# - I_range               investment cost uncertainty range
# - subs_rate             land subsidence rate (m/yr)
# - design_surge_level    Table #2 from Jonkman et al (2009)
# - investments_US        Table #2 from Jonkman et al (2009)
#                         These two above give the cost in US$ to build a levee
#                         to withstand the given surge level.
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

brick_vandantzig <- function( sea_level,
                              time
                              )
{

source("../R/VD_NOLA.R")

## Set up
currentyear			= 2015		# start year
endyear					= 2100		# reevaluation year
lowleveeheight	= 0				# lowest levee heightening considered (m)
highleveeheight	= 10			# largest levee heightening considered (m)
increment				= 0.05		# increment of levee heightening range (m)

## Time horizon until next evaluation of levees (years)
T = endyear - currentyear
time_frame=time

## Range of considered levee heights (meters)
X = seq(from = lowleveeheight, to = highleveeheight, by = increment)

## Make sure sea_level has time as the first dimension (number of rows = number of years)
dims = dim(sea_level)
if(dims[1]!=length(time)) {sea_level=t(sea_level)}

## Create a matrix of sea-level data starting from the current year for each model simulation.
## (but does not include subsidence, which is accounted for as uncertain parameter)
local_sea_level = data.matrix(sea_level[match(currentyear, time_frame):match(endyear, time_frame),])

## Set up basic information for the van Dantzig model.

## Number of observations (ensemble members):
n_obs = length(local_sea_level[1, ])

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

## Sample parameters (using Dutch Perspective as a guide)
## Uses Latin Hypercube to sample Van Dantzig parameters, assigns each ensemble
## member of SLR a set of parameters.
params.vd = parameter_sampling_DP(n_obs, p_zero_p, alpha_p, V_p, delta_prime_p, I_range, subs_rate)

## Test simulation to get the names
vanDantzig.out = VD_NOLA_R(params.vd[1,], I_table, T, X, local_sea_level[,1], intercept.h2i, slope.h2i)

n.ensemble = nrow(params.vd)
vanDantzig.ensemble = vector("list", ncol(vanDantzig.out))
names(vanDantzig.ensemble) = names(vanDantzig.out)
for (i in 1:ncol(vanDantzig.out)) {
  vanDantzig.ensemble[[i]] = mat.or.vec(nrow(vanDantzig.out), n.ensemble)
}

## Run the simulations

print(paste('Starting the risk assessment simulations now...',sep=''))

pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {
  vanDantzig.out = VD_NOLA_R(params.vd[i,], I_table, T, X, local_sea_level[,i], intercept.h2i, slope.h2i)
  for (j in 1:ncol(vanDantzig.out)) {
		vanDantzig.ensemble[[j]][,i] = vanDantzig.out[,j]
	}
  setTxtProgressBar(pb, i)
}
close(pb)

print(paste(' ... done with the risk assessment!',sep=''))

  return(vanDantzig.ensemble)
}
