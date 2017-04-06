##==============================================================================
## Source:
## http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
## http://www.somersault1824.com/wp-content/uploads/2015/02/color-blindness-palette.png
##
## Code by Tony Wong (twong@psu.edu)
##==============================================================================
mycol=rbind(
              c(0,0,0),
              c(0,73,73),
              c(0,146,146),
              c(255,109,182),
              c(255,182,119),
              c(73,0,146),
              c(0,109,219),
              c(182,109,255),
              c(109,182,255),
              c(182,219,255),
              c(146,0,0),
              c(146,73,0),
              c(219,109,0),
              c(36,255,36),
              c(255,255,109)
            )
mycol=mycol/max(mycol)

mycol.rgb <- rep(0,nrow(mycol))
for (i in 1:nrow(mycol)) {
    mycol.rgb[i] <- rgb(mycol[i,1],mycol[i,2],mycol[i,3])
}

##==============================================================================
## End
##==============================================================================
