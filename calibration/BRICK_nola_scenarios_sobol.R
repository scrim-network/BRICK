##==============================================================================
## Sobol sensitvity analysis for drivers of flood risk.
## Need to set up the design matrix with each of the three levers of RCP
## scenario, AIS fast dynamics contribution to sea level (in 2065), and surge
## factor parameter (0 or between 1.5 and 2), assocaited with the model response
## of average AEP between 2015 and 2065 (or return period, inversely).
##
## Questions? Tony Wong (twong@psu.edu)
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

# how large of a design matrix will you have?
numrow <- n.ensemble[[1]]*length(scen.ais)*length(scen.rcp)*length(scen.ss)

design.matrix <- data.frame(mat.or.vec(numrow,4))

colnames(design.matrix) <- c('surge.rise','rcp','ais','return.period')

# get the sea-levle rise projections, specifically the AIS fast dynamical contributions
ncdata <- nc_open(filename.brick.allslr)
    disint_gamma_rcp26 <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP26') -
                          ncvar_get(ncdata, 'LocalSeaLevel_gamma_nofd_RCP26')
    disint_gamma_rcp45 <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP45') -
                          ncvar_get(ncdata, 'LocalSeaLevel_gamma_nofd_RCP45')
    disint_gamma_rcp85 <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP85') -
                          ncvar_get(ncdata, 'LocalSeaLevel_gamma_nofd_RCP85')
    disint_uniform_rcp26 <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP26') -
                            ncvar_get(ncdata, 'LocalSeaLevel_uniform_nofd_RCP26')
    disint_uniform_rcp45 <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP45') -
                            ncvar_get(ncdata, 'LocalSeaLevel_uniform_nofd_RCP45')
    disint_uniform_rcp85 <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP85') -
                            ncvar_get(ncdata, 'LocalSeaLevel_uniform_nofd_RCP85')
    time_proj <- ncvar_get(ncdata, 'time_proj')
nc_close(ncdata)
inorm <- which(time_proj==2015)
iproj <- which(time_proj==2065)
disint_gamma_rcp26 <- disint_gamma_rcp26[iproj,]-disint_gamma_rcp26[inorm,]
disint_gamma_rcp45 <- disint_gamma_rcp45[iproj,]-disint_gamma_rcp45[inorm,]
disint_gamma_rcp85 <- disint_gamma_rcp85[iproj,]-disint_gamma_rcp85[inorm,]
disint_uniform_rcp26 <- disint_uniform_rcp26[iproj,]-disint_uniform_rcp26[inorm,]
disint_uniform_rcp45 <- disint_uniform_rcp45[iproj,]-disint_uniform_rcp45[inorm,]
disint_uniform_rcp85 <- disint_uniform_rcp85[iproj,]-disint_uniform_rcp85[inorm,]
disint_none_rcp26 <- rep(0, length(disint_gamma_rcp26))
disint_none_rcp45 <- rep(0, length(disint_gamma_rcp26))
disint_none_rcp85 <- rep(0, length(disint_gamma_rcp26))

# assign rcp scenario
design.matrix$rcp <- as.numeric(substring(scen.names[return.period[,1]], 4,6))

# assign return periods
design.matrix$return.period <- return.period[,2]

for (i in 1:length(scen.names)) {

    # indices of all simulations in this scenario
    iscen <- which(return.period[,1]==i)

    # assign storm surge factor
    surge <- substring(scen.names[i], (nchar(scen.names[i])-1),nchar(scen.names[i]) )
    if(surge=='st') {
        design.matrix$surge.rise[iscen] <- rep(0, length(iscen))
    } else if(surge=='ns') {
        design.matrix$surge.rise[iscen] <- surge.factor
    } else {
        print('ERROR - unrecognized storm surge')
    }

    # assign AIS fast dynamics contribution to sea-level (in 2065)
    aisfd <- substring(scen.names[i], 8, 11)
    if(aisfd=='none') {
        design.matrix$ais[iscen] <- rep(0, length(iscen))
    } else if(aisfd=='gamm') {
        if(median(design.matrix$rcp[iscen])==2.6) {
            design.matrix$ais[iscen] <- disint_gamma_rcp26
        } else if(median(design.matrix$rcp[iscen])==4.5) {
            design.matrix$ais[iscen] <- disint_gamma_rcp45
        } else if(median(design.matrix$rcp[iscen])==8.5) {
            design.matrix$ais[iscen] <- disint_gamma_rcp85
        } else {
            print('ERROR - unrecognized AIS-FD')
        }
    } else if(aisfd=='unif') {
        if(median(design.matrix$rcp[iscen])==2.6) {
            design.matrix$ais[iscen] <- disint_uniform_rcp26
        } else if(median(design.matrix$rcp[iscen])==4.5) {
            design.matrix$ais[iscen] <- disint_uniform_rcp45
        } else if(median(design.matrix$rcp[iscen])==8.5) {
            design.matrix$ais[iscen] <- disint_uniform_rcp85
        } else {
            print('ERROR - unrecognized AIS-FD')
        }
    }

}

# install/use the library for Sobol sensitivity analysis
# install.packages('sensitivity')
library(sensitivity)

# Split sample design matrix into two sub-samples
i1 <- sample(x=1:numrow, size=floor(numrow/2), replace=FALSE)
thing1 <- design.matrix[i1,]
thing2 <- design.matrix[-i1,]

# rethink how this is done - sample each of AIS, RCP and Surge from [0,1], then
# map to their respective ranges (RCP discrete in 3 bins, AIS is )
n <- 1000
X1 <- data.frame(matrix(runif(3*n), nrow=n))
X2 <- data.frame(matrix(runif(3*n), nrow=n))

sobol.samples <- sobol2002(model=NULL, thing1, thing2, nboot = 0, conf = 0.95)


##TODO
##TODO
##TODO
# here, we can write a function to take a set of the three inputs (RPC has to be
# discrete), run BRICK, and map to a new output
##TODO
##TODO
##TODO

testfun <- function( x.in ){
    return(as.numeric(x.in[length(x.in)]))
}

##==============================================================================
## End
##==============================================================================
