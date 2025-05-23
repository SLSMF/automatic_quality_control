calculate_distinctness <- function(data) {
  data$day = lubridate::date(data$time)
  distinctness=dplyr::summarise(group_by(data, day), distinctness=n_distinct(slevel)/n()) |> 
    dplyr::pull(distinctness)
  
  return(list(result=distinctness, calculation_parameters=list()))
}
