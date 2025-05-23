calculate_flat_line <- function(data, sample_rate=NA, min_value_frequency=40, min_flat_line_size=40, max_gap_size=2) {
  
  output = data |>
    dplyr::mutate(date=as.Date(time)) |>
    dplyr::group_by(date) |>
    dplyr::reframe(flat_line = calculate_flat_lines_in_day(
      data = data.frame(slevel = slevel, timeslot_id = timeslot_id), 
      sample_rate, min_value_frequency, min_flat_line_size, max_gap_size))
  
  return(list(result=output$flat_line, calculation_parameters=list()))
}

calculate_flat_lines_in_day <- function(data, sample_rate, min_value_frequency, min_flat_line_size, max_gap_size) {
  #processes one day of one sensorid
  
  #TODO: increase or decrease min values in case sample rate is too high or too low
  values_too_frequent=which(table(data$slevel) >= min_value_frequency) |> 
    names() |>
    as.numeric()
  
  flat_line_data_index=lapply(values_too_frequent, function(v) {
    index=which(data$slevel==v)
    #max_gap_size + 1 used as minimum difference between points before splitting
    index_consecutive=split_at(index, (which(diff(index) > (max_gap_size+1))+1))
    
    index_size=sapply(index_consecutive, length)
    index_consecutive[which(index_size >= min_flat_line_size)]
  }) |> 
    unlist() |>
    sort()
  
  flat_line_index=data$timeslot_id[flat_line_data_index]
  flat_line=data$timeslot_id %in% flat_line_index
  
  return(flat_line)
}