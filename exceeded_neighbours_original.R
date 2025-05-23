
calculate_exceeded_neighbours <- function(time, slevel, timeslot_id, max_diff_both_neighbours=0.3, max_diff_either_neighbour=0.5) {
  # remove values that are much higher or lower than their neighbours
  # 1) peak higher at least max_diff_both_neighbours from both its neighbours
  # 2) peak lower at least max_diff_both_neighbours from both its neighbors
  # 3) peak higher/lower at least max_diff_both_neighbours from one of its neighbor

  # remove values that are much higher or lower than their neighbors

  ##detect gaps, because exceeded neighbours frequently marks start and end of gap
  #maybe there is a tolerance for small gaps to achieve this timeslot_id_diff can be larger than 1
  ##this will detect less gaps which will result with more exceeded neighbours flagged
  data = data.frame(time, slevel, timeslot_id)
  
  data=data |>
    dplyr::mutate(date=as.Date(time)) |>
    dplyr::group_by(date) |>
    dplyr::mutate(
      timeslot_id_diff=c(1, diff(timeslot_id)),
      gap_end=timeslot_id_diff > 1,
      gap_start=lead(gap_end),
      gap=(gap_start | gap_end))

  slevel=data$slevel
  ndatamax <- length(slevel)
  neighbourflags <- rep(1,ndatamax) # set flags initially to 1 as default
  for (j in 1:(ndatamax)) {

    if ( (slevel[j]>(slevel[j-1]+max_diff_both_neighbours) && slevel[j]>(slevel[j+1]+max_diff_both_neighbours)) ||
         (slevel[j-1]>(slevel[j]+max_diff_both_neighbours) && slevel[j+1]>(slevel[j]+max_diff_both_neighbours)) ||
         abs(slevel[j]-slevel[j-1])>max_diff_either_neighbour ||
         abs(slevel[j]-slevel[j+1])>max_diff_either_neighbour
    ) {
      neighbourflags[j] <- 4
    }
  }
  data$exceeded_neighbours <- neighbourflags > 1
  # data$exceeded_neighbours=(data$gap==FALSE & data$exceeded_neighbours==TRUE)

  return(list(result=data$exceeded_neighbours, calculation_parameters=list()))
}