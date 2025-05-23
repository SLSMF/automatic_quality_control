
calculate_shift <- function(data, data_moonmonth) {
  #detect break by comparing msl daily simple mean with q10 and q90 of previous month without current day
  
  
  if (nrow(data_moonmonth)==0) {
    return(list(result=FALSE, calculation_parameters=list()))
  }
  
  
  #remove current day
  x = dplyr::mutate(data_moonmonth, date = lubridate::date(time)) |>
    dplyr::filter(date != as.Date(data$time[1]))
  if (nrow(x) != 0) {
    data_moonmonth=x
  }
  
  #to calculate msl_daily_simple_mean, it is better to first remove out of range and neighbours, but out of range doesn't work well if shift happened, so here is a temporary mean of all raw values of day
  msl_daily_mean = mean(data$slevel)
  
  cumulative_distribution_function <- ecdf(data_moonmonth$slevel)
  q_msl <- cumulative_distribution_function(msl_daily_mean)
  
  shift=(q_msl < 0.1  || q_msl > 0.9)
  return(list(result=shift, calculation_parameters=list(q_msl=q_msl)))
}

