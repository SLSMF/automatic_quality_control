# SPIKE_DETECT_SPLINE_MEDIAN
#
# Spike detection function using the median method
#
# Usage: flags=spike_detect_spline_median(heights,dates,Nwin,nMAD,minspike,ndegree)
#
# This version is a function, accepts data and parameters for spline
#
# Input parameters:
# heights: vector with sea level height data 
# dates: vector with dates corresponding to samples in vector 'heights'
# Nwin: number of data of the window to fit the spline (recommended value = 60)
# nMAD: number of scaled median absolute deviations to identify an event (recommended value = 6)
# minspike: minimum deviation to identify an event (recommended value = 0)
# ndegree: degree of the spline (polynomial, ndegree=[1,2], recommended value = 2)
# mthd: method to compute fit; 1=use polyfitNoTests.R, 2=use R's fastLm
#
# One problem that remains is that if there are a number of outliers close
# together, they will likely affect the spline. A quadratic fitted using
# least squares (as in polyfit) is not robust against an outlier.
# A robust method for fitting the quadratic would be better, but it is very
# likely to severely impact performance.
#
# Instead, it is recommended to run the spike detection (at least) twice,
# removing the spikes after each pass.
#
# Flags returned are 1 for GOOD data, 4 for BAD data (consistent with EuroGOOS)
# 
# ---------------
# Code history: 
# Based on python code from Begona Perez Gomez, Created on Oct 2017 %  author: EPPE\bpg
# 
# Converted to Matlab June 2018 and optimised by Joanne Williams (NOC)
# 
# Edited 9 Jan 2019 (Joanne Williams): change test from >= to >. Otherwise constant
# heights (or where heights are exactly equal to a spline) will also be
# reported as a SPIKE.
# If the instrument is stuck it should return a STUCK flag instead.
#
# Edited 1st Feb 2019 (Joanne Williams): allow option for median absolute deviation instead
# of RMSE difference. This is a more robust measure of the deviation.
# Also the point being tested is excluded from the spline fit.
#
# Split to two different functions (one using the RMSE and one using the 
# median method) and converted to R November 2023 by Nikos Kalligeris
# --------------------
# 
# Joanne Williams, Feb 2019, joll@noc.ac.uk

spike_detect_spline_median <- function(heights,dates,Nwin,nMAD,minspike,
                                       ndegree,mthd) {

  # library(caTools)
  # library(RcppEigen)
  # 
  # source("qc_scripts/polyfitNoTests.R")
  
  ##########################################################################
  # standard inputs
  erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)
  medianscaling <- -1/(sqrt(2)*erfcinv(3/2))  
  
  # find NA values
  flags <- rep(1,length(heights)) # set flags initially to 1 as default
  isn <- is.na(heights)
  heights <- heights[!isn] # without NA values
  dates <- dates[!isn] # without NA values
  
  # convert dates to numeric (to days to keep numbers smaller)
  numDates <- as.numeric(dates)/ 86400 / 365
  
  #
  ndatamax <- length(heights) # length of data input vectors
  ifl <- rep(1,ndatamax) # set flags initially to 1 as default
  
  ##########################################################################
  # loop through data to fit spline
  
  for (i in 1:ndatamax) {
    #print(i)
    
    ####################################################################
    # find i-dependent indices for Nwin points 
    if (i-1 < Nwin/2) {
      j1<-1
      j2<-Nwin
    } else if (i-1 >= ndatamax - Nwin/2) {
      j1<-ndatamax-Nwin
      j2<-ndatamax
    } else {
      j1<-i-Nwin/2
      j2<-i+Nwin/2
    }
    
    # exclude current value for spline fitting
    if (i-1 < 1) {
      ind<-2:j2
    } else if (i+1 > ndatamax) {
      ind<-j1:ndatamax-1
    } else {
      ind <- j1:(i-1)
      tmp <- (i+1):j2
      ind[length(ind)+1:length(tmp)]<-tmp
    }
    
    # apply indices to get input dates and heights for spline fitting
    datewin <- numDates[ind]
    Datawin <- heights[ind]
    
    ####################################################################
    # fit spline to input dates and heights
    
    if (mthd == 1) {
      
      datemu1 <- (numDates[j2]+numDates[j1])/2
      datemu2 <- (numDates[j2]-numDates[j1])/3
      
      xin <- (datewin-datemu1)/datemu2
      
      p <- polyfitNoTests(xin,Datawin,ndegree)
      
      if (ndegree == 1) {
        yfit <- p[2] + p[1]*xin
        yifit <- p[2] + p[1]*(numDates[i]-datemu1)/datemu2
      } else if (ndegree == 2) {
        yfit <- p[3] + p[2]*xin + p[1]*xin^2
        yifit <- p[3] + p[2]*((numDates[i]-datemu1)/datemu2) +
          p[1]*((numDates[i]-datemu1)/datemu2)^2
      } else {
        print('ndegree > 2 is not supported currently for spline fitting')
        break
      }
    } else if (mthd == 2) {
      
      if (ndegree == 1) {
        xyData <- data.frame(x=Datawin,y=datewin)
        tmp <- data.frame(x=heights[i],y=numDates[i])
        pfOut <- fastLm(formula=x ~  y, data=xyData, method=2)
        yfit <- predict(pfOut, data.frame(x=Datawin,y=datewin))
        yifit <- predict(pfOut, newdata=tmp)} else if (ndegree == 2) {
          xyData <- data.frame(x=Datawin,y=datewin,y2=datewin^2)
          tmp <- data.frame(x=heights[i],y=numDates[i],y2=numDates[i]^2)
          pfOut <- fastLm(formula=x ~  y + y2, data=xyData, method=2)
          yfit <- predict(pfOut, data.frame(x=Datawin,y=datewin,y2=datewin^2))
          yifit <- predict(pfOut, newdata=tmp)} else {
            print('ndegree > 2 is not supported currently for spline fitting')
            break
          }
    }
    
    # case median: median absolute deviation
    A <- yfit-Datawin
    scaledMAD <- medianscaling*median(abs(A-median(A)))
    if (abs(yifit-heights[i])>nMAD*scaledMAD & 
        abs(yifit-heights[i])>minspike) {
      ifl[i] <- 4 # flag for spike detected
      # print(i)
    }
  }
  flags[!isn] <- ifl
  return(flags) 
}
