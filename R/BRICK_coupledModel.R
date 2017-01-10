##==============================================================================
## Define the coupled BRICK model and its inputs/outputs
## The inputs/outputs will vary based on which components you use, but an overview
## is given below.
##
## Input:
##	parameters.in					input vector of model parameters
##	parnames.in						vector of parameter names
##	forcing.in						matrix of radiative forcing input
##	l.project							making projections or hindcasts?
##	slope.Ta2Tg.in				slope of Antarctic vs global mean temperature regression
##	intercept.Ta2Tg.in		intercept of Antarctic vs global mean temperature regression
##	mod.time							time (in years) of the model simulation
##	obs.temp							mean global surface temperature anomalies, used if no climate module
##	ind.norm.data					indices within the model output for setting zero anomaly of calibration data fields
##	ind.norm.sl						indices within model output for setting zero sea level
##	timestep							model timestep [years]
##
## Requires:
##  luse.brick, includes: luse.doeclim, luse.gsic, luse.te, luse.simple,
##                        luse.dais, and luse.XXX, where XXX
##                        may be replaced with your favorite model component
##
## Questions? Tony Wong <twong@psu.edu>
##==============================================================================

brick_model = function(
												parameters.in,
												parnames.in,
												forcing.in,
												l.project = FALSE,
												slope.Ta2Tg.in = 1,
												intercept.Ta2Tg.in = 0,
												mod.time,
												obs.temp = NULL,
												ind.norm.data = NULL,
												ind.norm.sl = NULL,
												luse.brick,
												timestep = 1
												){

	# Initialize the list of output (do NOT grow lists/arrays in R)
	# The +1 is to have total global mean sea level (relative to ind.relative) in
	# the output.
	brick.out = vector('list',sum(luse.brick)+1)
	slr.out = rep(0,length(mod.time))
	outcnt=1

	# Initialize temperature for coupling
	temp.couple = rep(NA, length(mod.time))

	#=============================================================================
	# DOECLIM - climate and ocean energy balance

  if (luse.brick[,"luse.doeclim"]) {

	  ## Grab the DOECLIM parameters
    S            =parameters.in[match("S"            ,parnames.in)]
		kappa.doeclim=parameters.in[match("kappa.doeclim",parnames.in)]
		alpha.doeclim=parameters.in[match("alpha.doeclim",parnames.in)]
		T0           =parameters.in[match("T0"           ,parnames.in)]
		H0           =parameters.in[match("H0"           ,parnames.in)]

		## Set up the radiative forcing
		forcing.total = forcing_total(forcing=forcing.in,
																	alpha.doeclim=alpha.doeclim,
																	l.project=l.project,
																	begyear=mod.time[1],
																	endyear=mod.time[length(mod.time)]
																	)

		## Run DOECLIM at these parameter values
		doeclim.out = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=forcing.total, mod.time=mod.time)

		## Normalize temperature and ocean heat to match the observations
		itmp = ind.norm.data[match("temp",ind.norm.data[,1]),2]:ind.norm.data[match("temp",ind.norm.data[,1]),3]
		doeclim.out$temp = doeclim.out$temp - mean(doeclim.out$temp[itmp])

		#itmp = ind.norm.data[match("ocheat",ind.norm.data[,1]),2]:ind.norm.data[match("ocheat",ind.norm.data[,1]),3]
		#doeclim.out$ocheat = doeclim.out$ocheat - mean(doeclim.out$ocheat[itmp])

		temp.couple = doeclim.out$temp + T0
		brick.out[[outcnt]] = doeclim.out; names(brick.out)[outcnt]="doeclim.out"; outcnt=outcnt+1;

	} else {temp.couple=obs.temp}

	#=============================================================================
	# GSIC-MAGICC - glaciers and small ice caps

  if (luse.brick[,"luse.gsic"]) {

		## Grab the GSIC parameters
		beta0  = parameters.in[match("beta0"  ,parnames.in)]
		V0.gsic= parameters.in[match("V0.gsic",parnames.in)]
		n      = parameters.in[match("n"      ,parnames.in)]
		Gs0    = parameters.in[match("Gs0"    ,parnames.in)]

		## Run GSIC-MAGICC at these parameter values, using temperature output from DOECLIM
		gsic.out = gsic_magiccF(beta0=beta0, V0=V0.gsic, n=n, Gs0=Gs0 , Tg=temp.couple)

		## Subtract off normalization period model GSIC output as the zero point
		itmp = ind.norm.data[match("gsic",ind.norm.data[,1]),2]:ind.norm.data[match("gsic",ind.norm.data[,1]),3]

		gsic.out.norm = gsic.out - mean(gsic.out[itmp])

		brick.out[[outcnt]] = gsic.out.norm; names(brick.out)[outcnt]="gsic.out"; outcnt=outcnt+1;

		## Add this contribution to the total sea level rise
		slr.out = slr.out + (gsic.out - mean(gsic.out[ind.norm.sl]))

  }

	#=============================================================================
	# TE - thermosteric expansion

	if (luse.brick[,"luse.te"]) {

	  ## Grab the BRICK-TE parameters
		a.te  =parameters.in[match("a.te"  ,parnames.in)]
		b.te  =parameters.in[match("b.te"  ,parnames.in)]
		invtau.te=parameters.in[match("invtau.te",parnames.in)]
		TE0   =parameters.in[match("TE0"   ,parnames.in)]

    ## Run BRICK-TE (thermosteric expansion) model, using temp output from DOECLIM
    te.out = brick_te_F(a=a.te , b=b.te, invtau=invtau.te, TE_0=TE0, Tg=temp.couple)

		## Subtract off normalization period
		itmp = ind.norm.data[match("te",ind.norm.data[,1]),2]:ind.norm.data[match("te",ind.norm.data[,1]),3]
		te.out.norm = te.out - mean(te.out[itmp])

		brick.out[[outcnt]] = te.out.norm; names(brick.out)[outcnt]="te.out"; outcnt=outcnt+1;

		## Add this contribution to the total sea level rise
		slr.out = slr.out + (te.out.norm - mean(te.out.norm[ind.norm.sl]))

  }

	#=============================================================================
	# SIMPLE - Greenland ice sheet

	if (luse.brick[,"luse.simple"]) {

	  ## Grab SIMPLE parameters
		a.simple    =parameters.in[match("a.simple"    ,parnames.in)]
		b.simple    =parameters.in[match("b.simple"    ,parnames.in)]
		alpha.simple=parameters.in[match("alpha.simple",parnames.in)]
		beta.simple =parameters.in[match("beta.simple" ,parnames.in)]
		V0          =parameters.in[match("V0"          ,parnames.in)]

    ## Run SIMPLE (Greenland Ice Sheet model)
    simple.out = simpleF(a=a.simple, b=b.simple, alpha=alpha.simple,
											 	 beta=beta.simple, V0=V0, Tg=temp.couple)

		## Add this contribution to the total sea level rise
		slr.out = slr.out + (simple.out$sle.gis - mean(simple.out$sle.gis[ind.norm.sl]))

		## Subtract off normalization period
		itmp = ind.norm.data[match("gis",ind.norm.data[,1]),2]:ind.norm.data[match("gis",ind.norm.data[,1]),3]
		simple.out$sle.gis = simple.out$sle.gis - mean(simple.out$sle.gis[itmp])

		brick.out[[outcnt]] = simple.out; names(brick.out)[outcnt]="simple.out"; outcnt=outcnt+1;

  }

	#=============================================================================
	# DAIS - Antarctic ice sheet

  if (luse.brick[,"luse.dais"]) {

		## Grab DAIS parameters
		anto.a=parameters.in[match("anto.a",parnames.in)]
		anto.b=parameters.in[match("anto.b",parnames.in)]
		gamma =parameters.in[match("gamma" ,parnames.in)]
		alpha.dais =parameters.in[match("alpha.dais" ,parnames.in)]
		mu =parameters.in[match("mu" ,parnames.in)]
		nu =parameters.in[match("nu" ,parnames.in)]
		P0 =parameters.in[match("P0" ,parnames.in)]
		kappa.dais =parameters.in[match("kappa.dais" ,parnames.in)]
		f0 =parameters.in[match("f0" ,parnames.in)]
		h0 =parameters.in[match("h0" ,parnames.in)]
		c =parameters.in[match("c" ,parnames.in)]
		b0 =parameters.in[match("b0" ,parnames.in)]
		slope =parameters.in[match("slope" ,parnames.in)]

		## Calculate the sea level updated from the other model components'
		## contributions. From Shaffer (2014), SL should be relative to 1961-1990
		## mean. Implement fingerprinting of local sea-level sources on AIS?

		l.fprint=TRUE

		SL.couple = slr.out
		if(l.fprint) {
			dSL.gis = diff(simple.out$sle.gis)
			dSL.gsic= diff(gsic.out)
			dSL.te	= diff(te.out)
			for (i in 2:length(mod.time)) {
				SL.couple[i] = SL.couple[i-1] + timestep*(	1.1*dSL.gis[i] +
																										1.1*dSL.gsic[i] +
																										1.0*dSL.te[i]
																										)
			}
		}

		## Normalize
		SL.couple = SL.couple - mean(SL.couple[ind.norm.sl])

		itmp = ind.norm.data[match("sl",ind.norm.data[,1]),2]:ind.norm.data[match("sl",ind.norm.data[,1]),3]
		SL.couple = slr.out
  	SL.couple = SL.couple - mean(SL.couple[itmp])
		dSL.couple = c(-999,diff(slr.out))
		include_dSLais = 0		# in coupled model, feeding AIS dSL without AIS contribution

		## Check to make sure output from other models was reasonable
		if(any(is.na(SL.couple))) {
			slr.out=rep(NA,length(mod.time))
			brick.out[[outcnt]] = slr.out; names(brick.out)[outcnt]="dais.out"; outcnt=outcnt+1;
		}	else {

  		dais.out = daisantoF(
#			dais.out = daisanto(
                       		anto.a=anto.a , anto.b=anto.b,
                       		gamma=gamma   , alpha=alpha.dais,
                       		mu=mu         , nu=nu        ,
                       		P0=P0         , kappa=kappa.dais,
                       		f0=f0         , h0=h0        ,
                       		c=c           , b0=b0        ,
                       		slope=slope   ,
											 		slope.Ta2Tg=slope.Ta2Tg.in, intercept.Ta2Tg=intercept.Ta2Tg.in,
                       		Tg=temp.couple, SL=SL.couple, dSL=dSL.couple ,
													includes_dSLais = include_dSLais)

		  ## Subtract off normalization period
		  itmp = ind.norm.data[match("ais",ind.norm.data[,1]),2]:ind.norm.data[match("ais",ind.norm.data[,1]),3]
  		dais.out.norm = dais.out - mean(dais.out[itmp])

    	brick.out[[outcnt]] = dais.out.norm; names(brick.out)[outcnt]="dais.out"; outcnt=outcnt+1;

			## Add this contribution to the total sea level rise
			slr.out = slr.out + (dais.out.norm - mean(dais.out.norm[ind.norm.sl]))

  	}
	}

	#=============================================================================
	# Total sea-level rise

	## Add the SLR to the output
	brick.out[[sum(luse.brick)+1]] = slr.out - mean(slr.out[ind.norm.sl])
  names(brick.out)[sum(luse.brick)+1]="slr.out"

	## Check to make sure all the output made it
	if(outcnt!=sum(luse.brick)+1) print('ERROR - missing model output!')

	#=============================================================================

	return(brick.out)
}



##==============================================================================
## End
##==============================================================================
