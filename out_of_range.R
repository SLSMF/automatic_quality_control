calculate_out_of_range <- function(data) {
  # extreme values
  # extr_val <- 100 # remove out of range values above |extr_val| #need to take into account height differences of stations
  
  # ndatamax <- length(data$slevel)
  # outofrangeflags <- rep(1,ndatamax) # set flags initially to 1 as default
  
  #take into account height difference of stations: https://www.ioc-sealevelmonitoring.org/service.php#outlier
  # All values X where (abs(X) – abs(median)) > tolerance are hidden.
  # With tolerance = 3*abs(percentile90 - median)
  # The statistics (median and percentile90) are calculated on the data being plotted (12h, 1d, 7d, 30d)
  
  # lg1 <- abs(data$slevel) > extr_val
  
  msl=median(data$slevel)
  tolerance=3*abs(quantile(data$slevel, probs = 0.9) - msl)
  out_of_range=abs(data$slevel - msl) > tolerance
  
  return(list(result=out_of_range, calculation_parameters=list(tolerance=tolerance)))
}
