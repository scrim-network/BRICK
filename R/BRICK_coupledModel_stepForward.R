##==============================================================================
## Step fully coupled BRICK model forward one timestep
##
## Input:
##	parameters					input vector of model parameters
##	parnames						vector of parameter names
##	forcing							matrix of radiative forcing input
##	l.project						making projections or hindcasts?
##	slope.Ta2Tg					slope of Antarctic vs global mean temperature regression
##	intercept.Ta2Tg			intercept of Antarctic vs global mean temperature regression
##	mod.time						time (in years) of the model simulation. for stepForward
##											BRICK function, mod.time should have the current [1] and
##											next [2] time steps' times. This corresponds to the two
##											entries in 'forcing'
##	timestep						model timestep [years]
##	luse.brick					which submodels are called (logical list)
##	SL.old							sea level from previous time step (mod.time[1])
##	l.fprint						(logical) fingerprint to sea level at Antarctica for DAIS?
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

brick_model_stepforward = function(
	parameters,
	parnames,
	forcing,
	l.project = FALSE,
	slope.Ta2Tg = 1,
	intercept.Ta2Tg = 0,
	mod.time,
	luse.brick,
	timestep = 1,
	SL.old = 0,
	l.fprint = TRUE
){

	# read parameters
	## Grab the DOECLIM parameters
	S            =parameters[,"S"            ]
	kappa.doeclim=parameters[,"kappa.doeclim"]
	alpha.doeclim=parameters[,"alpha.doeclim"]
	T0           =parameters[,"T0"           ]
	H0           =parameters[,"H0"           ]

	## Grab the GSIC parameters
	beta0  = parameters[,"beta0"  ]
	V0.gsic= parameters[,"V0.gsic"]
	n      = parameters[,"n"      ]
	Gs0    = parameters[,"Gs0"    ]

  ## Grab the BRICK-TE parameters
	a.te  =parameters[,"a.te"  ]
	b.te  =parameters[,"b.te"  ]
	invtau.te=parameters[,"invtau.te"]
	TE0   =parameters[,"TE0"   ]

  ## Grab SIMPLE parameters
	a.simple    =parameters[,"a.simple"    ]
	b.simple    =parameters[,"b.simple"    ]
	alpha.simple=parameters[,"alpha.simple"]
	beta.simple =parameters[,"beta.simple" ]
	V0          =parameters[,"V0"          ]

	## Grab DAIS parameters
	anto.a=parameters[,"anto.a"]
	anto.b=parameters[,"anto.b"]
	gamma =parameters[,"gamma"]
	alpha.dais=parameters[,"alpha.dais"]
	mu=parameters[,"mu"]
	nu=parameters[,"nu"]
	P0 =parameters[,"P0"]
	kappa.dais =parameters[,"kappa.dais"]
	f0 =parameters[,"f0"]
	h0 =parameters[,"h0"]
	c =parameters[,"c"]
	b0 =parameters[,"b0"]
	slope =parameters[,"slope"]

##==============================================================================

	# Initialize the list of output (do NOT grow lists/arrays in R)
	# The +1 is to have total global mean sea level (relative to ind.relative) in
	# the output.
	brick.out <- vector('list',sum(luse.brick)+1)
	outcnt <- 1

##==============================================================================

	# step the forcing forward
	forcing.total <- forcing_total(	forcing=forcing,
																	alpha.doeclim=alpha.doeclim,
																	l.project=l.project,
																	begyear=mod.time[1],
																	endyear=mod.time[length(mod.time)]
																	)

#print(paste(mod.time[1],mod.time[2],forcing.total))

##==============================================================================

if(luse.brick[,'luse.doeclim']){

	# step DOECLIM forward
	doeclim.out <- doeclimF(S=S, kappa=kappa.doeclim, forcing.total=forcing.total, mod.time=mod.time)

	brick.out[[outcnt]] <- doeclim.out; names(brick.out)[outcnt] <- "doeclim.out"; outcnt <- outcnt+1;

	# old and new temperatures for coupling to sea level components
	# NB: the DOECLIM output temperature and ocean heat uptake do *not* have the
	#			offsets T0 and H0 added in to them, respectively; the coupling temperature
	#			temp.couple does.
	temp.couple <- doeclim.out$temp + T0

	# Normalize temperature to match what the sub-model expects (the parameters
	# may assume a particular time period associated with Tg=0, for example)
	# GSIC-MAGICC expects temp.couple relative to late 1800s, which it already is
	# with ind.norm.data for temperature rel to 1850-70 (Wigley and Raper 2005)

# Normalize for what each component expects <-- TODO: what to do about this?

}

##==============================================================================

if(luse.brick[,'luse.gsic']) {

	# step glaciers and small ice caps forward
	gsic.out <- gsic_magiccF(beta0=beta0, V0=V0.gsic, n=n, Gs0=Gs0 , Tg=temp.couple)
	brick.out[[outcnt]] <- gsic.out; names(brick.out)[outcnt] <- "gsic.out"; outcnt <- outcnt+1;

}

##==============================================================================

if(luse.brick[,'luse.te']) {

	# step thermal expansion forward
  te.out <- brick_te_F(a=a.te , b=b.te, invtau=invtau.te, TE_0=TE0, Tg=temp.couple)
	brick.out[[outcnt]] <- te.out; names(brick.out)[outcnt] <- "te.out"; outcnt <- outcnt+1;

}

##==============================================================================

if(luse.brick[,'luse.simple']) {

	# step SIMPLE (Greenland ice sheet) forward
  simple.out <- simpleF(a=a.simple, b=b.simple, alpha=alpha.simple,
										 	 	beta=beta.simple, V0=V0, Tg=temp.couple)
	brick.out[[outcnt]] <- simple.out; names(brick.out)[outcnt] <- "simple.out"; outcnt <- outcnt+1;

}

##==============================================================================

if(luse.brick[,'luse.dais']) {

	# step DAIS (Antarcic ice sheet) forward

##
#TODO
##
## fingerprinting not working - if you only have one time slice, then diff won't work
##
#TODO
##

	if(l.fprint) {
		dSL.gis <- diff(simple.out$sle.gis)
		dSL.gsic <- diff(gsic.out)
		dSL.te <- diff(te.out)
		SL.couple <- c(SL.old, SL.old + timestep*(1.1*dSL.gis + 1.1*dSL.gsic + 1.0*dSL.te))
	}

	# Normalize <-- TODO: what to do about this?
	#SL.couple <- SL.couple - mean(SL.couple[ind.norm.sl])
	dSL.couple <- SL.couple - SL.old

	include_dSLais <- 0		# in coupled model, feeding AIS dSL without AIS contribution

	# Check to make sure output from other models was reasonable
	if(any(is.na(SL.couple))) {
		slr.out <- rep(NA,length(mod.time))
		brick.out[[outcnt]] <- slr.out; names(brick.out)[outcnt] <- "dais.out"; outcnt <- outcnt+1;
	}	else {
		dais.out <- daisantoF(anto.a=anto.a , anto.b=anto.b,
                     			gamma=gamma   , alpha=alpha.dais,
                     			mu=mu         , nu=nu        ,
                     			P0=P0         , kappa=kappa.dais,
                     			f0=f0         , h0=h0        ,
                     			c=c           , b0=b0        ,
                     			slope=slope   ,
										 			slope.Ta2Tg=slope.Ta2Tg, intercept.Ta2Tg=intercept.Ta2Tg,
                     			Tg=temp.couple, SL=SL.couple, dSL=dSL.couple ,
													includes_dSLais = include_dSLais)
		brick.out[[outcnt]] <- dais.out; names(brick.out)[outcnt] <- "dais.out"; outcnt <- outcnt+1;
  }

}

##==============================================================================

	# add up total global mean sea level rise forward, and step forward
	slr.out <- SL.old + dSL.gis + dSL.gsic + dSL.te + timestep*diff(dais.out)
	brick.out[[sum(luse.brick)+1]] <- slr.out; names(brick.out)[sum(luse.brick)+1] <- "slr.out"

##==============================================================================

	# Check to make sure all the output made it
	if(outcnt!=sum(luse.brick)+1) print('ERROR - missing model output!')

##==============================================================================

	return(brick.out)
}

##==============================================================================
## End
##==============================================================================
