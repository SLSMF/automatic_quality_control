
calculate_spikes_via_rmse <- function(data, Nwin=4, ndegree=2, mthd=1, nsigma=4) {
  time = data$time
  slevel = data$slevel
  
  #detect spikes using the NOC (spline-fitting) methodology, RMSE method
  
  # Nwin <- 4 # number of data of the window to fit the spline (recommended value = 60)
  # ndegree <- 2 # degree of the spline (polynomial, ndegree=[1,2], recommended value = 2)
  # mthd <- 1 # method to compute fit; 1=use polyfitNoTests.R, 2=use R's fastLm (recommended value = 1)
  # For fitting the spline, I struggled a bit with R’s polynomial-fitting functions. While R has more and richer functions for polynomial fitting compared to Matlab, none of them runs very fast.
  # The fastest I could find and use is fastLm, however, it is not nearly as fast and doesn’t give the same results as NOC’s polynomial fitting Matlab function.
  # I therefore decided to translate NOC’s polyfitNoTests Matlab script to R and call that instead of R’s built in functions. 
  # Using polyfitNoTests.R, the spike detection routine now runs almost as fast as Matlab, and gives exactly the same results, at least for the Ajaccio data tested here.
  # Nevertheless, I have kept the option of using fastLm for polynomial fitting by setting variable mthd <- 2 in the user inputs. 
  
  # nsigma is number of sigmas to identify an event (recommended value = 4)
  
  
  spike_flags1 <- spike_detect_spline_rmse(slevel,time,Nwin,nsigma,ndegree,mthd)
  lg3 <- spike_flags1 == 4 
  # tmp <- length(lg3[lg3==TRUE])
  # cat(tmp," spikes detected after 1 round\n")
  # run spike detection script twice, as recommended in NOC manual
  slevel[lg3] <- NA
  spike_flags2 <- spike_detect_spline_rmse(slevel,time,Nwin,nsigma,ndegree,mthd)
  spikes_via_rmse <- spike_flags2 == 4 | spike_flags1 == 4
  # tmp <- length(lg3[lg3==TRUE])
  # cat(tmp," spikes detected after 2 rounds\n")
  # spike_flags <- rep(1,length(slevel)) 
  # spike_flags[lg3] <- 4
  
  return(list(result=spikes_via_rmse, calculation_parameters=list()))
}