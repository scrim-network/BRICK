##==============================================================================
##
## Define parameters and their prior ranges
##
## Note1 -- when defining the names, be sure to use the names which are looked up
##			 -- in the BRICK model code that separates the parameters for each model
##       -- Valid parameter names, and their models can be found in the README file
##
## Question? Tony Wong <twong@psu.edu>
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

## DOECLIM (Urban and Keller, 2010, values
parnames.doeclim   =NULL; p0.doeclim       =NULL; bound.lower.doeclim=NULL;
bound.upper.doeclim=NULL; step.mcmc.doeclim=NULL; index.model.doeclim=NULL;
if (luse.doeclim) {
	parnames.doeclim   =c("S" ,"kappa.doeclim","alpha.doeclim","T0"  ,"H0" ,"sigma.T","sigma.H","rho.T","rho.H")	# parameters names
	p0.doeclim		     =c(3.1 , 3.5   , 1.1           , -0.06, -33 , 0.1     , 2       , 0.55  , 0.9   )	# initial parameter guesses
	bound.lower.doeclim=c(0.1 , 0.1   , 0             , -0.3 , -50 , 0.05    , 0.1     , 0     , 0     )	# prior range lower bounds
	bound.upper.doeclim=c(10  , 4     , 2             ,  0.3 ,   0 , 5       , 10      , 0.999 , 0.999 )	# prior range upper bounds
	step.mcmc.doeclim	 =c(0.16, 0.17  ,0.025          ,0.003 , 0.9 , 5e-4    , 0.025   , 0.007 , 0.006 )	# step size for parameters in MCMC (proposals)
	index.model.doeclim=c(1,2,3,4,5)		# which are model parameters? (index within parnames.doeclim)
}

## GSIC-MAGICC
parnames.gsic   =NULL; p0.gsic       =NULL; bound.lower.gsic=NULL;
bound.upper.gsic=NULL; step.mcmc.gsic=NULL; index.model.gsic=NULL;
if (luse.gsic) {
	parnames.gsic   =c("beta0","V0.gsic","n"  ,"Gs0"    , "sigma.gsic", "rho.gsic")	# parameters names
	p0.gsic		      =c(0.00058, 0.41    , 0.82, 0.0     , 0.00045     , 0.5       )	# initial parameter guesses
	bound.lower.gsic=c(0      , 0.3     , 0.55, -0.0041 , 0           , -0.999    )	# prior range lower bounds
	bound.upper.gsic=c(0.041  , 0.5     , 1.0 ,  0.0041 , 0.00150     ,  0.999    )	# prior range upper bounds
	step.mcmc.gsic	=c(0.01   , 0.01    , 0.1 , 0.01    , 0.0001      , 0.01      )	# step size for parameters in MCMC (proposals)
	index.model.gsic=c(1,2,3,4)			# which are model parameters? (index within parnames.gsic)
}

## BRICK-TE
parnames.te   =NULL; p0.te       =NULL; bound.lower.te=NULL;
bound.upper.te=NULL; step.mcmc.te=NULL; index.model.te=NULL;
if (luse.te) {
	parnames.te   =c("a.te","b.te" ,"invtau.te" ,"TE0"  )        # parameters names
	p0.te		      =c(0.3616, 0.5   , 1/200      , 0.0   )        # initial parameter guesses
	bound.lower.te=c(0.0   , 0.0   , 0          ,-0.0484)        # prior range lower bounds
	bound.upper.te=c(0.8595, 2.193 , 1          , 0.0484)        # prior range upper bounds
	step.mcmc.te	=0.05*(bound.upper.te-bound.lower.te) # set size for parameters in MCMC (proposals)
	index.model.te=c(1,2,3,4)			# which are model parameters? (index within parnames.te)
}

## SIMPLE
parnames.simple   =NULL; p0.simple       =NULL; bound.lower.simple=NULL;
bound.upper.simple=NULL; step.mcmc.simple=NULL; index.model.simple=NULL;
if (luse.simple) {
	parnames.simple   =c("a.simple" ,"b.simple" ,"alpha.simple","beta.simple","V0"      ,"sigma.simple") # parameters names
	p0.simple         =c(-0.825     , 7.36      , 1.63e-4      , 2.85e-5     , 7.36     , 5e-4         ) # initial parameter guesses
	bound.lower.simple=c( -4        , 5.888     , 0            , 0           , 7.16     , 0            ) # prior range lower bounds
	bound.upper.simple=c( -1e-3     , 8.832     , 1e-3         , 1e-3        , 7.56     , 0.002        ) # prior range upper bounds
	step.mcmc.simple  =c( 0.2       , 0.05      , 1e-5         , 1e-5        , 0.05     , 0.0001       ) # step sizes for initial MCMC
	index.model.simple=c(1,2,3,4,5)			# which are model parameters? (index within parnames.simple)
}

## DAIS
parnames.dais   =NULL; p0.dais       =NULL; bound.lower.dais=NULL;
bound.upper.dais=NULL; step.mcmc.dais=NULL; index.model.dais=NULL;
if (luse.dais) {
	parnames.dais   =c('anto.a','anto.b','gamma','alpha.dais','mu'  ,'nu'  ,'P0' ,'kappa.dais','f0','h0'  ,'c'  , 'b0','slope','var.dais' ) # parameter names
	p0.dais         =c( 0.1574 , 0.6677 ,  2    , 0.35       , 8.7  , 0.012, 0.35, 0.04       , 1.2, 1471 , 95  , 775 , 0.0006,0.0004656  )
	bound.lower.dais=c( 0.0    , 0      ,  0.5  , 0          , 7.05 , 0.003,0.026, 0.025      , 0.6, 735.5, 47.5, 740 , 4.5e-4,0          )
	bound.upper.dais=c( 1.0    , 2      ,  4.25 , 1          , 13.65, 0.015, 1.5 , 0.085      , 1.8,2206.5,142.5, 820 , 7.5e-4,2          )
	step.mcmc.dais	=0.05*(bound.upper.dais-bound.lower.dais) # set size for parameters in MCMC (proposals)
	index.model.dais=c(1,2,3,4,5,6,7,8,9,10,11,12,13)			# which are model parameters? (index within parnames.dais)
}

##==============================================================================
## Combine for coupled model
parnames    = c(parnames.doeclim, parnames.gsic, parnames.te, parnames.simple, parnames.dais)
p0          = c(p0.doeclim      , p0.gsic      , p0.te      , p0.simple      , p0.dais      )
bound.lower = c(bound.lower.doeclim, bound.lower.gsic, bound.lower.te, bound.lower.simple,
								bound.lower.dais   )
bound.upper = c(bound.upper.doeclim, bound.upper.gsic, bound.upper.te, bound.upper.simple,
								bound.upper.dais   )
step.mcmc   = c(step.mcmc.doeclim  , step.mcmc.gsic  , step.mcmc.te  , step.mcmc.simple  ,
								step.mcmc.dais     )
index.model = c(match(parnames.doeclim[index.model.doeclim],parnames),
								match(parnames.gsic[index.model.gsic]      ,parnames),
								match(parnames.te[index.model.te]          ,parnames),
								match(parnames.simple[index.model.simple]  ,parnames),
								match(parnames.dais[index.model.dais]      ,parnames))
index.all = 1:length(p0); index.stat = index.all[is.na(pmatch(index.all,index.model))]

## Reshape so all the model parameters are first; easier bookkeeping for DEoptim
parnames = c(parnames[index.model],parnames[index.stat])
p0 = c(p0[index.model],p0[index.stat])
bound.lower = c(bound.lower[index.model],bound.lower[index.stat])
bound.upper = c(bound.upper[index.model],bound.upper[index.stat])
step.mcmc   = c(step.mcmc[index.model],step.mcmc[index.stat])
index.model = 1:length(index.model)
index.stat  = (length(index.model)+1):length(p0)

##==============================================================================
## End
##==============================================================================
