
calculate_spikes_via_median <- function(data, Nwin=60, ndegree=2, mthd=1, nMAD=6, minspike=3) {
  time = data$time
  slevel = data$slevel
  
  #detect spikes using the NOC (spline-fitting) methodology, median method
  
  # Nwin <- 4 # number of data of the window to fit the spline (recommended value = 60)
  # ndegree <- 2 # degree of the spline (polynomial, ndegree=[1,2], recommended value = 2)
  # mthd <- 1 # method to compute fit; 1=use polyfitNoTests.R, 2=use R's fastLm (recommended value = 1)
  # For fitting the spline, I struggled a bit with R’s polynomial-fitting functions. While R has more and richer functions for polynomial fitting compared to Matlab, none of them runs very fast.
  # The fastest I could find and use is fastLm, however, it is not nearly as fast and doesn’t give the same results as NOC’s polynomial fitting Matlab function.
  # I therefore decided to translate NOC’s polyfitNoTests Matlab script to R and call that instead of R’s built in functions. 
  # Using polyfitNoTests.R, the spike detection routine now runs almost as fast as Matlab, and gives exactly the same results, at least for the Ajaccio data tested here.
  # Nevertheless, I have kept the option of using fastLm for polynomial fitting by setting variable mthd <- 2 in the user inputs. 
  
  # nMAD is number of scaled median absolute deviations of the data over Nwin points to identify an event (recommended value = 6)
  # minspike is minimum number of typical standard deviations to identify an event (NOC default value = 0)
  
  
  # Estimate the usual standard deviation of the data over Nwin points. This is
  # used as a measure of noise for parameter selection.
  # runsd is different to movstd of Matlab in that
  # movstd for the first element computes std(x(1:k/2+i-1))
  # whereas runsd computes std(x(1:k2+i))
  tmp <-runsd(slevel, Nwin, 
              center = runmean(slevel,Nwin),
              endrule=c("func"))
  typicalstd <- median(tmp) # typical standard deviation of input data over Nwin points
  minspike <- minspike*typicalstd
  
  # run spike detection script
  spike_flags1 <- spike_detect_spline_median(slevel,time,Nwin,nMAD,minspike,
                                             ndegree,mthd)
  lg3 <- spike_flags1 == 4 
  # tmp <- length(lg3[lg3==TRUE])
  # cat(tmp," spikes detected after 1 round\n")
  # run spike detection script twice, as recommended in NOC manual
  slevel[lg3] <- NA
  spike_flags2 <- spike_detect_spline_median(slevel,time,Nwin,nMAD,minspike,
                                             ndegree,mthd)
  spikes_via_median <- spike_flags2 == 4 | spike_flags1 == 4
  # tmp <- length(lg3[lg3==TRUE])
  # cat(tmp," spikes detected after 2 rounds\n")
  # spike_flags <- rep(1,length(slevel)) 
  # spike_flags[lg3] <- 4
  
  return(list(result=spikes_via_median, calculation_parameters=list()))
}