##==============================================================================
## BRICK_Sobol_functions.R
##
## Original codes by Calvin Whealton, Cornell University
## https://github.com/calvinwhealton/SensitivityAnalysisPlots
## and
## Perry Oddo, Penn State
##
## Modified (condensed - codes are almost entirely unchanged) for brevity by
## Tony Wong (twong@psu.edu). The only code that was changed is the plotting
## routine 'plotRadCon'; Tony added
## (1) inputs for the first-, total- and second-order % values used for legend
## (2) generate labels for first-, total- and second-order in legend
## (3) write legend in black ink instead of gray
## (4) include '%' sign and cut legend labels off at the decimal (whole numbers
##     only)
##
## Tony also modified the 'sig' test for significance to test for confidence
## interval bounds of the same sign (otherwise, 0 is in CI for sensitivity
## index) and greater than 1%.
##==============================================================================

##==============================================================================
###################################
# file: upper.diag.R
###################################
# Author and copyright: Perry Oddo
# Pennsylvania State University
# poddo@psu.edu
###################################
# Function to convert second-order output of SALib analysis
# into upper-triangular matrix for radial convergence plotting
####################################

upper.diag<-function(x){
  m<-(-1+sqrt(1+8*length(x)))/2
  X<-lower.tri(matrix(NA,m,m),diag=TRUE)
  X[X==TRUE]<-x
  X[upper.tri(X, diag = FALSE)] <- NaN
  t(X)
}
##==============================================================================

##==============================================================================
## Original routine: groupAssign.R

# functions for substituting group information into dataframe

# function for assigning group name and color to group
gp_name_col <- function(name_list  # list of variables names for each group
                        ,col_list   # list of colors for each group
                        ,df         # data frame of values with ind as the variable name
                        ){

  # initializing columns
  df$gp_name <- NA # group name
  df$gp_col <- NA  # group colors

  # checking for the same variable names
  n1 <- sort(names(name_list))
  n2 <- sort(names(col_list))
  if(length(setdiff(n1,n2)) != 0){
    print('Group names do not match across the two lists')
  }
  else{
    # loop over the variables in each group and assign group name
    for(i in 1:length(n1)){

      # extracting the variable names for the given element
      var_names <- unlist(name_list[names(name_list)[i]])

      # loop over the values in each list
      for(j in 1:length(var_names)){
        # substituting group name in for the variable
        df$gp_name[which(df$Parameter %in% var_names[j])] <- names(name_list)[i]
      }
    }

    # loop over the group names and assign the color
    for(i in 1:length(names(col_list))){
      df$gp_col[which(df$gp_name %in% unlist(names(name_list))[i])] <- unlist(col_list[names(col_list)[i]])
    }
  }

  # returning data frame with additional columns
  return(df)
}
##==============================================================================

##==============================================================================
## Original routine: sigTests.R

# functions for testing statistical significance
# determines which indices to plot based on values

# function inputs are:
#   df = data frame with sensitivity indices (S1 and ST)
#         includes columns for S1, ST, S1_conf, and ST_conf
#   method = method of testing
#           'sig' is statistically significant with alpha
#           'gtr' is greater than a specified value
#   sigCri = significance criteria
#             'either' parameter is signficant if either S1 or ST is significant
#             'S1' parameter is significant if S1 (or S1 and ST) is significant
#             'ST' parameter is significant if ST (or S1 and ST) is significant

#####################################################
# (Tony-modified) -- function for testing significance of S1 and ST
# functions assume the confidence are for already defined type I error
stat_sig_s1st <- function(df
                     ,greater = 0.01
                     ,method='sig'
                     ,sigCri = 'either'){

  # initializing columns for the statistical significance of indices
  df$s1_sig <- 0
  df$st_sig <- 0
  df$sig <- 0

  # testing for statistical significance
  if(method == 'sig'){
    # testing for statistical significance using the confidence intervals
    df$s1_sig[which(abs(s1st$S1) - s1st$S1_conf > 0)] <- 1
    df$st_sig[which(abs(s1st$ST) - s1st$ST_conf > 0)] <- 1
  }
  else if(method == 'gtr'){
    # finding indicies that are greater than the specified values
    df$s1_sig[which(abs(s1st$S1) > greater)] <- 1
    df$st_sig[which(abs(s1st$ST) > greater)] <- 1
  } else if(method == 'con') {
    df$s1_sig[which(s1st$S1_conf_low * s1st$S1_conf_high > 0)] <- 1
    df$st_sig[which(s1st$ST_conf_low * s1st$ST_conf_high > 0)] <- 1
  } else if(method == 'congtr'){
    df$s1_sig[which(s1st$S1_conf_low * s1st$S1_conf_high > 0 &
                    abs(s1st$S1) > greater)] <- 1
    df$st_sig[which(s1st$ST_conf_low * s1st$ST_conf_high > 0 &
                    abs(s1st$ST) > greater)] <- 1
  } else {
    print('Not a valid parameter for method')
  }

  # determining whether the parameter is significant
  if(sigCri == 'either'){
    for(i in 1:nrow(df)){
      df$sig[i] <- max(df$s1_sig[i],df$st_sig[i])
    }
  }
  else if(sigCri == 'S1'){
    df$sig <- df$s1_sig
  }
  else if(sigCri == 'ST'){
    df$sig <- df$st_sig
  }
  else{
    print('Not a valid parameter for SigCri')
  }

  # returned dataframe will have columns for the test of statistical significance
  return(df)
}

#####################################################
# function to test statistical significane of S2 indices
stat_sig_s2 <- function(dfs2
                          ,dfs2Conf_low
                          ,dfs2Conf_high
                          ,greater = 0.01
                          ,method='sig'){

  # initializing matrix to return values
  s2_sig <- matrix(0,nrow(s2),ncol(s2))

  # testing for statistical significance
  if(method == 'sig'){
    # testing for statistical significance using the confidence intervals
    s2_sig[which(abs(s2) - s2_conf > 0)] <- 1
  }
  else if(method == 'gtr'){
    # finding indicies that are greater than the specified values
    s2_sig[which(abs(s2) > greater)] <- 1
  }
  else if(method == 'con'){
    s2_sig[which(dfs2Conf_low * dfs2Conf_high > 0)] <- 1
  }
  else if(method == 'congtr'){
    s2_sig[which(dfs2Conf_low * dfs2Conf_high > 0 &
                 abs(s2) > greater)] <- 1
  }
  else{
    print('Not a valid parameter for method')
  }

  # returned dataframe will have columns for the test of statistical significance
  return(s2_sig)
}
##==============================================================================


##==============================================================================
## Original routine: plotRadSAinds.R

# legFirLabs, legSecLabs and legTotLabs added by Tony Wong (twong@psu.edu)

# function for plotting the radial convergence plots
# function takes input data frame

plotRadCon <- function(df                   # dataframe with S1 and ST indices
                       ,s2                  # S2 indices
                       ,s2_sig              # S2 significance matrix
                       ,filename = 'plot'   # file name for the saved plot
                       ,plotType = 'EPS'    # plot type
                       ,plotS2 = TRUE       # whether to plot S2 indices
                       ,radSc = 2           # radius scaling of entire plot
                       ,scaling=1           # scaling factor for plot
                       ,widthSc = 0.5       # power used in scaling width, 0.5 is root, 1 is simple multiple
                       ,STthick = 0.05      # value used in determining the width of the ST circle
                       ,RingThick = 0.08    # thickness of ring plotted around variable symbols
                       ,line_col = 'gray48'#"#7A7A7ABF" # color used for lines
                       ,st_col = 'black'    # color for total-order index circle
                       ,s1_col = 'gray48'   # color for first-order index disk(filled circle)
                       ,asp=1               # aspect ratio
                       ,varNameMult = 1.25   # location of variable name with respect to the plot radius
                       ,gpNameMult = 1.45    # location of the group name with respect to the plot radius
                       ,legLoc = 'topleft'  # legend location
                       ,legThick=c(0.1,0.5) # legend thickensses
                       ,legPos=1.9         # legend relative position
                       ,cex = 1
                       ,legFirLabs=NULL     # legend labels for first order
                       ,legSecLabs=NULL     # legend labels for second order
                       ,legTotLabs=NULL     # legend labels for total order
                       ,lBuildRCPhoriz=FALSE # horizontal legends for Emissions and Protection?
                       ,lnoGEVhoriz=FALSE   # horizontal legends for Emissions, Protection, Storm Surge, and Land Subsidence?
                       ,lnoHRhoriz=FALSE   # horizontal legends for Emissions, Protection?
                       ,lsetback=FALSE   # don't rotate, just stagger labels radially?
                       ){

  # Shift plot up
  shift = 1

  # finding number of points to plot
  num_plot <- n_params#sum(df$sig)

  # polar cooridantes angular-values of locations
  angles <- radSc*pi*seq(0,num_plot-1)/num_plot

  # assigning coordinates to varaibles based on groups
  df$rad <- radSc
  df$ang <- NA

  ## coordinates in polar for each variable
  # finding number of groups with a significant variable
  sig_gps <- unique(df$gp_name)#unique(df$gp_name[which(df$sig %in% 1)])#

  # initializing vector to hold the number of significant variables for each group
  num_sig_gp <- rep(0,length(sig_gps))

  # counter used for indexing the values in the angle vector
  counter <- 0

  for(i in 1:length(sig_gps)){
    # indices of variables in group and significant
    sig_in_gp <- intersect(which(df$gp_name %in% sig_gps[i]),which(df$sig >= 0))
      #intersect(which(df$gp_name %in% sig_gps[i]),which(df$sig %in% 1))

    # taking sequential values in the angles vector
    df$ang[sig_in_gp] <- angles[seq(counter+1,counter+length(sig_in_gp),1)]

    # indexing counter
    counter <- counter + length(sig_in_gp)

    # vector for counting number of statistically signficant variables in the applicable groups
    num_sig_gp[i] <- length(sig_in_gp)
  }

  ## converting to Cartesian coordinates
  df$x_val <- df$rad*cos(df$ang)
  df$y_val <- df$rad*sin(df$ang) + shift

  # colors and scales used in plots
  line_scaling <- widthSc*scaling#line_sc_mult*scaling

  ## file set-up storage
  if(plotType == 'EPS'){
    fname <- paste(filename,'.eps',sep='')
    savePlot <- TRUE
    setEPS()
    postscript(fname)
  } else {
    print('Plot not automatically saved')
    savePlot <- FALSE

  }
  ## plotting
  # initial plot is empty
  plot(NA
       , NA
       , xlim = c(-2*radSc,2*radSc)
       , ylim = c(-2*radSc,2*radSc)
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = ''
       , ylab = ''
       #, bty = 'n'
       ,asp=asp)

  # plot circle
  draw.circle(0
              ,0 + shift
              ,radius <- varNameMult * radSc * (1 + RingThick)
              ,nv=200
              ,border=NA
              ,col="gray90"
  )
  draw.circle(0
              ,0 + shift
              ,radius <- varNameMult * radSc * (1 - RingThick)
              ,nv=200
              ,border=NA
              ,col="white"
  )

  # plotting all lines that were significant----
  if(plotS2 == TRUE){
    for(i in 1:ncol(s2_sig)){       # i indexes across the rows
      for(j in 1:nrow(s2_sig)){     # j indexes down the column

        # only plot second order when the two indices are significant
        if(s2_sig[j,i]*(df$sig[i]*df$sig[j]) == 1){
#        if(s2_sig[j,i] == 1){
          # coordinates of the center line
          clx <- c(df$x_val[i],df$x_val[j])
          cly <- c(df$y_val[i],df$y_val[j])

          # calculating the angle of the center line
          # calculating tangent as opposite (difference in y)
          # divided by adjacent (difference in x)
          clAngle1 <- atan((cly[2]-cly[1])/(clx[2]-clx[1]))

          # adding angle when both values are negative to make it on (-pi/2,3*pi/2)
          if(cly[2]-cly[1] < 0){
            clAngle1 <- clAngle1 + pi
          }

          # half width of the line
          line_hw <- scaling*(s2[j,i]^widthSc)/2

          # color of line
          max_col <- max(s2[j,i])

          # creating vector of polygon coordinates
          polyx <- rep(0,4)
          polyy <- rep(0,4)

          polyx[1] <- clx[1] - line_hw*sin(clAngle1)
          polyx[2] <- clx[1] + line_hw*sin(clAngle1)
          polyx[3] <- clx[2] + line_hw*sin(clAngle1)
          polyx[4] <- clx[2] - line_hw*sin(clAngle1)

          polyy[1] <- cly[1] + line_hw*cos(clAngle1)
          polyy[2] <- cly[1] - line_hw*cos(clAngle1)
          polyy[3] <- cly[2] - line_hw*cos(clAngle1)
          polyy[4] <- cly[2] + line_hw*cos(clAngle1)

          # making polygons
          polygon(polyx,polyy
                  ,density=300
                  #,border=NA
                  ,border='black'
                  ,lwd=.3
                  ,col=line_col)
        }
      }
    }
  }

  for(i in 1:nrow(df)){
    if(df$sig[i] == 1){

      # circle for total order index
      draw.circle(df$x_val[i]#*0.95 # for k scaling
                  ,df$y_val[i]#*0.95 # for k scaling
                  ,radius <- scaling*(df$ST[i]^widthSc)/2
                  ,nv=200
                  ,border=NA
                  ,col=st_col
                  )

# Commented out by Tony 9 March 2017 to highlight the importance of total order
#      # white circle to make total-order an outline
#      draw.circle(df$x_val[i]#*0.95
#                  ,df$y_val[i]#*0.95
#                  ,radius <- (1-STthick)*scaling*(df$ST[i]^widthSc)/2
#                  ,nv=200
#                  ,border=NA
#                  ,col="white"
#                  )

      # gray circle for first-order
      draw.circle(df$x_val[i]#*0.95
                  ,df$y_val[i]#*0.95
                  ,radius <- scaling*(df$S1[i]^widthSc)/2
                  ,nv=200
                  ,border=NA
                  ,col=s1_col
                  )
    }
  }

  ## adding text to the plots
  # adding variable names
  for(i in 1:nrow(df)){
    if((df$ang[i]*360/(2*pi)) >= 0 & (df$ang[i]*360/(2*pi)) <= 180){
      if(is.na(df$ang[i]) == FALSE){
        text(varNameMult*df$rad[i]*cos(df$ang[i]), varNameMult*df$rad[i]*sin(df$ang[i]) + shift
             , df$symbols[i]#df$Symbol
             , cex = cex
             , col = df$gp_col[i]
             , adj = 0.5#0
             , font = 1
              #srt = df$ang[i]*360/(2*pi) #- 90 #0
        )
      }
    } else {
      if(is.na(df$ang[i]) == FALSE){
        text(varNameMult*df$rad[i]*cos(df$ang[i]), varNameMult*df$rad[i]*sin(df$ang[i]) + shift
             , df$symbols[i]#df$Symbol
             , col = df$gp_col[i]
             , cex = cex
             , adj = 0.5#0
             , font = 1
              #srt = df$ang[i]*360/(2*pi) #+ 90 #0
        )
      }
    }
  }

  # adding group names
  counter <- 0
  for(i in 1:length(num_sig_gp)){

    angle_gp <- mean(angles[seq(counter+1,counter+num_sig_gp[i],1)])
    #print(angle_gp[i] * 360/(2*pi))
    counter <- counter + num_sig_gp[i]

if(lBuildRCPhoriz) {
  if(sig_gps[i]=='Emissions' | sig_gps[i]=='Protection') {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = 0# angle_gp*360/(2*pi) + 90
           , adj = 0.3 # for centering
           , font = 1
      )
  } else {
    if((angle_gp*360/(2*pi)) >= 0 & (angle_gp*360/(2*pi)) <= 180){
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) - 90
           , adj = 0.5 # for centering
           , font = 1
      )
    } else {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90
           , adj = 0.5 # for centering
           , font = 1
      )
    }
  }
} else if(lnoGEVhoriz) {
  if(sig_gps[i]=='Emissions' | sig_gps[i]=='Protection' | sig_gps[i]=='Storm Surge' | sig_gps[i]=='Land Subsidence') {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = 0# angle_gp*360/(2*pi) + 90
           , adj = 0.15 # for centering
           , font = 1
      )
  } else {
    if((angle_gp*360/(2*pi)) >= 0 & (angle_gp*360/(2*pi)) <= 180){
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) - 90
           , adj = 0.5 # for centering
           , font = 1
      )
    } else {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90
           , adj = 0.5 # for centering
           , font = 1
      )
    }
  }
} else if(lnoHRhoriz) {
  if(sig_gps[i]=='Emissions' | sig_gps[i]=='Protection') {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = 0# angle_gp*360/(2*pi) + 90
           , adj = 0.15 # for centering
           , font = 1
      )
  } else {
    if((angle_gp*360/(2*pi)) >= 0 & (angle_gp*360/(2*pi)) <= 180){
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) - 90
           , adj = 0.5 # for centering
           , font = 1
      )
    } else {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90
           , adj = 0.5 # for centering
           , font = 1
      )
    }
  }
} else if(lsetback) {
  if(sig_gps[i]=='Emissions') {
      text(1.08*gpNameMult*radSc*cos(angle_gp), 1.08*gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90
           , adj = 0.55 # for centering
           , font = 1
      )
  } else if(sig_gps[i]=='Protection') {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90
           , adj = 0.35 # for centering
           , font = 1
      )
  } else if(sig_gps[i]=='Land Subsidence') {
      text(1.15*gpNameMult*radSc*cos(angle_gp), 1.15*gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90 + 5
           , adj = 0.35 # for centering
           , font = 1
      )
  } else if(sig_gps[i]=='Sea Level:\nLand Water Storage') {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90
           , adj = 0.75 # for centering
           , font = 1
      )
  } else {
    if((angle_gp*360/(2*pi)) >= 0 & (angle_gp*360/(2*pi)) <= 180){
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) - 90
           , adj = 0.5 # for centering
           , font = 1
      )
    } else {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt = angle_gp*360/(2*pi) + 90
           , adj = 0.5 # for centering
           , font = 1
      )
    }
  }
} else {
    if((angle_gp*360/(2*pi)) >= 0 & (angle_gp*360/(2*pi)) <= 180){
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt =  angle_gp*360/(2*pi) - 90
           , adj = 0.5 # for centering
           , font = 1
      )
    } else {
      text(gpNameMult*radSc*cos(angle_gp), gpNameMult*radSc*sin(angle_gp) + shift
           , sig_gps[i]
           , col = df$gp_col[which(df$gp_name %in% sig_gps[i])]#[i]]
           , cex = cex
           , srt =  angle_gp*360/(2*pi) + 90
           , adj = 0.5 # for centering
           , font = 1
      )
    }
  }
}

  ## adding legend
  if(legLoc == 'topleft'){
    xloc <- rep(-legPos*radSc ,length(legThick))
    yloc <- seq(legPos*radSc ,1*radSc ,by=-0.3*radSc )
  } else if(legLoc=='topright'){
    xloc <- rep(legPos*radSc ,length(legThick))
    yloc <- seq(legPos*radSc ,1*radSc ,by=-0.3*radSc )
  } else if(legLoc=='bottomleft'){
    xloc <- rep(-legPos*radSc ,length(legThick))
    yloc <- seq(-legPos*radSc ,-1*radSc ,by=0.3*radSc )
  } else if(legLoc=='bottomcenter'){
    xloc <- rep(-legPos*radSc*0.85 ,length(legThick))
    yloc <- seq(-legPos*radSc ,-1*radSc ,by=0.3*radSc )
  } else{
    xloc <- rep(legPos*radSc ,length(legThick))
    yloc <- seq(-legPos*radSc ,-1*radSc ,by=0.3*radSc )
  }

  # gray circle for legend (1st Order)
  for(i in 1:length(xloc)){
    if(is.null(legFirLabs)) {
      legend_max = max(df$S1)#scaling*(max(df$S1)^widthSc)/2
      legend_min = min(df$S1[df$S1>0])#scaling*(min(df$S1[df$S1>0])^widthSc)/2
      if(legend_max==legend_min){
        legend_min <- 0
      }
    } else {
      legend_min = legFirLabs[1]
      legend_max = legFirLabs[2]
    }
    legend_scale = c(legend_max, legend_min)

    draw.circle(-3.25#xloc[i]
                ,yloc[i]
                ,radius <- scaling*(legend_scale[i]^widthSc)/2#scaling*(legThick[i]^widthSc)/2
                ,nv=200
                ,border=NA
                ,col=s1_col
                )
    text(-2#xloc[i]*0.65,
         ,yloc[i]
#         ,as.character(formatC((legend_scale[i]*100), digits = 2, format = "f"))#legThick[i])
         ,paste(as.character(formatC((legend_scale[i]*100), digits = 0, format = "f")),'%',sep='')
         ,col='black'
         ,adj = c(1,0.5))
    }

text(-2.7, -2.7,'First-order', cex=0.9)

  # white circle for legend (Total Order)
  for(i in 1:length(xloc)){
    if(is.null(legTotLabs)) {
      legend_max = max(df$ST)#scaling*(max(df$S1)^widthSc)/2
      legend_min = min(df$ST[df$ST>0])#scaling*(min(df$S1[df$S1>0])^widthSc)/2
      if(legend_max==legend_min){
        legend_min <- 0
      }
    } else {
      legend_min = legTotLabs[1]
      legend_max = legTotLabs[2]
    }
    legend_scale = c(legend_max, legend_min)

    draw.circle(-0.5#xloc[i] + 3
                ,yloc[i]
                ,radius <- scaling*(legend_scale[i]^widthSc)/2#scaling*(legThick[i]^widthSc)/2
                ,nv=200
                ,border=st_col
                #,col="white"
                ,col=st_col
    )
    text(scaling*(legend_max^widthSc)/2 + 0.5
         ,yloc[i]
#         ,as.character(formatC((legend_scale[i]*100), digits = 2, format = "f"))#legThick[i])
         ,paste(as.character(formatC((legend_scale[i]*100), digits = 0, format = "f")),'%',sep='')
         ,col='black'
         ,adj = c(1,0.5))
    }

text(0, -2.7,'Total sensitivity', cex=0.9)

  # Lines (Second Order)
  for(i in 1:length(xloc)){
    if(is.null(legSecLabs)) {
      line_max = abs(max(s2_table[,3]))
      line_min = abs(min(s2_table[,3]))
      if(line_max==line_min){
        line_min <- 0
      }
    } else {
      line_min = legSecLabs[1]
      line_max = legSecLabs[2]
    }
    line_scale = c(line_max, line_min)
    line_hw[i] <- scaling*(line_scale[i]^widthSc)/2

    line_x <- rep(0,4)
    line_y <- rep(0,4)

    line_x[1] <- xloc[i] + 5.5
    line_x[2] <- xloc[i] + 6.25
    line_x[3] <- xloc[i] + 6.25
    line_x[4] <- xloc[i] + 5.5

    line_y[1] <- yloc[i] + line_hw[i]
    line_y[2] <- yloc[i] + line_hw[i]
    line_y[3] <- yloc[i] - line_hw[i]
    line_y[4] <- yloc[i] - line_hw[i]

    polygon(line_x, line_y,
            density=200
            #,border=NA
            ,border='black'
            ,lwd=.3
            ,col=line_col)

    text((xloc[i]*0.65) + 6, yloc[i]
#         ,as.character(formatC((line_scale[i]*100), digits = 2, format = "f"))#legThick[i])
         ,paste(as.character(formatC((line_scale[i]*100), digits = 0, format = "f")),'%',sep='')
         ,col='black'
         ,adj = c(1,0.5))

  }

text(3.0, -2.7,'Second-order', cex=0.9)

  # closing plot if save to external file
  if(savePlot == TRUE){
    dev.off()
    print(paste('Plot saved to ',fname,sep=''))
  }
}
##==============================================================================

##==============================================================================
## End
##==============================================================================
