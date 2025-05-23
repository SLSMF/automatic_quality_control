
calculate_msl_daily_doodson <- function(data_of_three_days) {
  
  #[39 hourly values needed]:
  # 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 [17 18 19 20 21 22 23 (day-1)
  # 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16  17 18 19 20 21 22 23 (day of which msl is calculated)
  # 0  1  2  3  4  5  6  7] 8  9 10 11 12 13 14 15 16  17 18 19 20 21 22 23 (day+1)
  
  doodson_coefficients=c(1,0,1,0,0,1,0,1,1,0,2,0,1,1,0,2,1,1,2,0,2,1,1,2,0,1,1,0,2,0,1,1,0,1,0,0,1,0,1)
  
  mslh=calculate_msl_hourly(data_of_three_days)$result |> unlist()
  mslh=mslh[18:56]
  
  #max 1 empty value allowed
  if (sum(is.na(mslh)) > 1) {
    message("can't calculate as more than one empty hour mean sea level")
    return(NA)
  }
  
  #doodson formula
  msl=sum(na.omit(mslh*doodson_coefficients))/30
  
  return(list(result = msl, calculation_parameters = list()))
}