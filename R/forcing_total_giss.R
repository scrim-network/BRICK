
forcing_total <- function(  forcing       ,
                            alpha.doeclim ,
                            l.project     ,
                            begyear       ,
                            endyear       ,
                            flnd = 0.29
){
  forcing.total = forcing$ghg + forcing$o3 + forcing$sh2o + forcing$stra + forcing$solar + forcing$land + alpha.doeclim*(forcing$refa + forcing$aie + forcing$bc + forcing$snow)
  
  ## Clip forcing at the beginning and end of the model simulation
  ibeg=which(forcing$year==begyear)
  iend=which(forcing$year==endyear)
  if(length(ibeg)==0 | length(iend)==0) print('ERROR - begyear/endyear not within forcing data')
  forcing.total = forcing.total[ibeg:iend]
  
  return(forcing.total)
}
