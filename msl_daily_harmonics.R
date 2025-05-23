calculate_msl_daily_harmonics <- function(data) {
  msl=TideHarmonics::ftide(x = data$slevel, dto = data$time, TideHarmonics::hc60)$msl
  return(list(result = msl, calculation_parameters = list()))
}