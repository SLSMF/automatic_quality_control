automatic_qc_noa <- function(dates, values, msl=NA, stdev=NA) {
  data=data.frame(dates,values)
  
  ##sea level spike detection of NOA
  # script accepts minimum 61 points
  # the R script to run the example of post-processing the Ajaccio_data and is based on Petra’s Matlab script shared with me (thank you again Petra!).
  # The script does: i) remove out of range values, ii) remove values that are much higher or lower than their neighbours, and iii) detect spikes using the NOC (spline-fitting) methodology.
  # To run it you need to set the working directory in the first line and select the parameters for spike detection in the user inputs section.
  
  # load(file = "Ajaccio_data.Rdata") # the Rdata file containing the Ajaccio tide gauge data used in Petra’s example.
  
  # source("qc_scripts/spike_detect_spline_median.R") # the R script for detecting spikes using the NOC median method and based on NOC’s Matlab script.
  # source("qc_scripts/spike_detect_spline_rmse.R")  # the R script for detecting spikes using the NOC RMSE method and based on NOC’s Matlab script.
  # 
  # library(caTools)
  
  ##########################################################################
  # user inputs 
  
  # extreme values
  # extr_val <- 100 # remove out of range values above |extr_val| #need to take into account height differences of stations
  
  
  # remove values that are much higher or lower than their neighbours
  # 1) peak higher at least n_val1 from both its neighbours 
  # 2) peak lower at least n_val1 from both its neighbors 
  # 3) peak higher/lower at least n_val2 from one of its neighbor
  n_val1 <- 0.3 
  n_val2 <- 0.5
  
  # spike detection
  spike_mthd <- 1 # spike detection method; 1=median, 2=RMSE
  Nwin <- 4 # number of data of the window to fit the spline (recommended value = 60)
  ndegree <- 2 # degree of the spline (polynomial, ndegree=[1,2], recommended value = 2)
  mthd <- 1 # method to compute fit; 1=use polyfitNoTests.R, 2=use R's fastLm (recommended value = 1)
  # For fitting the spline, I struggled a bit with R’s polynomial-fitting functions. While R has more and richer functions for polynomial fitting compared to Matlab, none of them runs very fast.
  # The fastest I could find and use is fastLm, however, it is not nearly as fast and doesn’t give the same results as NOC’s polynomial fitting Matlab function.
  # I therefore decided to translate NOC’s polyfitNoTests Matlab script to R and call that instead of R’s built in functions. 
  # Using polyfitNoTests.R, the spike detection routine now runs almost as fast as Matlab, and gives exactly the same results, at least for the Ajaccio data tested here.
  # Nevertheless, I have kept the option of using fastLm for polynomial fitting by setting variable mthd <- 2 in the user inputs. 
  
  
  # input only relevant to median method
  nMAD <- 6 # number of scaled median absolute deviations of the data over Nwin points to identify an event (recommended value = 6)
  minspike <- 3 # minimum number of typical standard deviations to identify an event (NOC default value = 0)
  # input only relevant to RMSE method
  nsigma <- 4 # number of sigmas to identify an event (recommended value = 4)
  
  
  #speed: 1982.43 records /sec
  
  
  # data$values[1]=90
  
  # shorten the input vectors for testing, otherwise comment out first line below
  # data<-data[-c(200001:length(data$values)),]
  dates <- data$dates
  heights <- data$values
  ndatamax <- length(heights)
  
  ##########################################################################
  # remove out of range sea level values
  
  outofrangeflags <- rep(1,ndatamax) # set flags initially to 1 as default
  
  #take into account height difference of stations: https://www.ioc-sealevelmonitoring.org/service.php#outlier
  # All values X where (abs(X) – abs(median)) > tolerance are hidden.
  # With tolerance = 3*abs(percentile90 - median)
  # The statistics (median and percentile90) are calculated on the data being plotted (12h, 1d, 7d, 30d)
  
  # lg1 <- abs(heights) > extr_val
  
  if (is.na(msl)) msl=median(heights)
  if (is.na(stdev)) {
    tolerance=3*abs(quantile(heights, probs = 0.9) - msl)
  } else {
    tolerance=3*abs(stdev)
  }
  
  lg1=(abs(heights) - abs(msl)) > tolerance
  
  outofrangeflags[lg1] <- 4 # set flags to 4 for out of range values
  
  ##########################################################################
  # remove values that are much higher or lower than their neighbors
  
  neighbourflags <- rep(1,ndatamax) # set flags initially to 1 as default
  for (j in 2:(ndatamax-1)) {
    
    if ( (heights[j]>(heights[j-1]+n_val1) && heights[j]>(heights[j+1]+n_val1)) ||
         (heights[j-1]>(heights[j]+n_val1) && heights[j+1]>(heights[j]+n_val1)) ||
         abs(heights[j]-heights[j-1])>n_val2 || 
         abs(heights[j]-heights[j+1])>n_val2 ) {
      neighbourflags[j] <- 4
    }
  }
  lg2 <- neighbourflags > 1 
  
  ##########################################################################
  # spike detection 
  
  
  if (spike_mthd == 1) {
    
    # median method
    
    # Estimate the usual standard deviation of the data over Nwin points. This is
    # used as a measure of noise for parameter selection.
    # runsd is different to movstd of Matlab in that
    # movstd for the first element computes std(x(1:k/2+i-1))
    # whereas runsd computes std(x(1:k2+i))
    tmp <-runsd(heights, Nwin, 
                center = runmean(heights,Nwin),
                endrule=c("func"))
    typicalstd <- median(tmp) # typical standard deviation of input data over Nwin points
    minspike <- minspike*typicalstd
    
    # run spike detection script
    spike_flags1 <- spike_detect_spline_median(heights,dates,Nwin,nMAD,minspike,
                                               ndegree,mthd)
    lg3 <- spike_flags1 == 4 
    tmp <- length(lg3[lg3==TRUE])
    # cat(tmp," spikes detected after 1 round\n")
    # run spike detection script twice, as recommended in NOC manual
    heights[lg3] <- NA
    spike_flags2 <- spike_detect_spline_median(heights,dates,Nwin,nMAD,minspike,
                                               ndegree,mthd)
    lg3 <- spike_flags2 == 4 | spike_flags1 == 4
    tmp <- length(lg3[lg3==TRUE])
    # cat(tmp," spikes detected after 2 rounds\n")
    spike_flags <- rep(1,length(heights)) 
    spike_flags[lg3] <- 4
  } else {
    
    # RMSE method
    
    spike_flags1 <- spike_detect_spline_rmse(heights,dates,Nwin,nsigma,ndegree,mthd)
    lg3 <- spike_flags1 == 4 
    tmp <- length(lg3[lg3==TRUE])
    # cat(tmp," spikes detected after 1 round\n")
    # run spike detection script twice, as recommended in NOC manual
    heights[lg3] <- NA
    spike_flags2 <- spike_detect_spline_rmse(heights,dates,Nwin,nsigma,ndegree,mthd)
    lg3 <- spike_flags2 == 4 | spike_flags1 == 4
    tmp <- length(lg3[lg3==TRUE])
    # cat(tmp," spikes detected after 2 rounds\n")
    spike_flags <- rep(1,length(heights)) 
    spike_flags[lg3] <- 4
  }
  ##########################################################################
  # plot detected spikes
  
  
  data$out_of_range=lg1
  data$exceeded_neighbours=lg2
  data$spike=lg3
  
  data=dplyr::rename(data, Time=dates, Slevel=values)
  return(data)
  # # plot
  # # quartz()
  # plot(data$dates,data$values,type="l", ylim = c(0,3), xlim= c(1700646819-(60*60*24*365*5.481),1700646819-(60*60*24*365*5.48)))
  # # out of range values plotted as red
  # points(data$dates[lg1],data$values[lg1],type="p",pch=19,col="red")
  # # values much heigher or lower than their neighbors plotted as blue
  # points(data$dates[lg2],data$values[lg2],type="p",pch=19,col="blue")
  # # spike flags plotted as green
  # points(data$dates[lg3],data$values[lg3],type="p",pch=19,col="green")
  
  ##########################################################################
  # write flags to .txt file
  #write.table(spike_flags, file = "Rflags.txt", sep = "\t",
  #            row.names = FALSE,col.names = FALSE)
  
}

# load(file = "Ajaccio_data.Rdata")
# data=Ajaccio_tmsr[1:61,]
# automatic_qc(data$DATE, data$ETA)
