##==============================================================================
## Script to use IPCC AR5 (Ch 13, p 1151) trends in land water storage to
## close the assumed global mean sea level budget.
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

## Sea-level obs accounting for land water contributions
## All sea-level (including land water here) relative to 1961-1990 mean 0

#1901-1990: –0.11 [–0.16 to –0.06] (5-95% range) [IPCC AR5, Ch13, p1151]
lw.time.1900 <- 1900:1989
i1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
lw.1900 <- (-0.11/1000)*(lw.time.1900 - 1989)
lw.1900.norm <- lw.1900 # no normalization needed when formed this way
lw.err.1900 <- (0.25*(-0.06--0.16)/1000)*sqrt(abs(lw.time.1900 - 1989))

#1971-2010: 0.12 [0.03 to 0.22]
lw.time.1970 <- 1970:2009
i1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
lw.1970 <- (0.12/1000)*(lw.time.1970 - 2009)
lw.1970.norm <- lw.1970
lw.err.1970 <- (0.25*(0.2-0.03)/1000)*sqrt(abs(lw.time.1970 - 2009))

#1993-2010: 0.38 [0.26 to 0.49]
lw.time.1992 <- 1992:2009
i1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])
lw.1992 <- (0.38/1000)*(lw.time.1992 - 2009)
lw.1992.norm <- lw.1992
lw.err.1992 <- (0.25*(0.49-0.26)/1000)*sqrt(abs(lw.time.1992 - 2009))

# normalize, subtract and add error in quadrature
itmp <- which(obs.sl.time==mod.time[ind.norm.sl[1]]):which(obs.sl.time==mod.time[ind.norm.sl[length(ind.norm.sl)]])
obs.sl <- obs.sl - mean(obs.sl[itmp])

# (... and re-normalize to END of land water period from IPCC trends)
# (do end because lower errors should be closer to normalization period - less ucnertain)
obs.sl.lw.1900 <- obs.sl[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])]
obs.sl.lw.1900 <- obs.sl.lw.1900 - obs.sl[(which(obs.sl.time==lw.time.1900[length(lw.time.1900)]))]
obs.sl.lw.1970 <- obs.sl[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])]
obs.sl.lw.1970 <- obs.sl.lw.1970 - obs.sl[(which(obs.sl.time==lw.time.1970[length(lw.time.1970)]))]
obs.sl.lw.1992 <- obs.sl[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])]
obs.sl.lw.1992 <- obs.sl.lw.1992 - obs.sl[(which(obs.sl.time==lw.time.1992[length(lw.time.1992)]))]

obs.sl.lw.1900 <- obs.sl.lw.1900 - lw.1900.norm
obs.sl.lw.1970 <- obs.sl.lw.1970 - lw.1970.norm
obs.sl.lw.1992 <- obs.sl.lw.1992 - lw.1992.norm

obs.sl.lw.err.1900 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])]^2 + lw.err.1900^2)
obs.sl.lw.err.1970 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])]^2 + lw.err.1970^2)
obs.sl.lw.err.1992 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])]^2 + lw.err.1992^2)

imod.1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
imod.1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
imod.1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])

obs.sl_lw <- vector("list", 3)
names(obs.sl_lw) <- c("r1900","r1970","r1992")
obs.sl_lw$r1900 <- list(obs.sl.lw.1900,obs.sl.lw.err.1900,imod.1900); names(obs.sl_lw$r1900) <- c('sl','err','midx')
obs.sl_lw$r1970 <- list(obs.sl.lw.1970,obs.sl.lw.err.1970,imod.1970); names(obs.sl_lw$r1970) <- c('sl','err','midx')
obs.sl_lw$r1992 <- list(obs.sl.lw.1992,obs.sl.lw.err.1992,imod.1992); names(obs.sl_lw$r1992) <- c('sl','err','midx')

##==============================================================================
## End
##==============================================================================
