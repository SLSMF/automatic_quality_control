



#min_gap_size=0 -> selects everything
#min_gap_size=1 -> . x .    -> . (gap) x . (gap)  
#min_gap_size=2 -> . x x .  -> . (gap) x x . (gap)  
calculate_gaps <- function(data, sample_rate, min_gap_size) {
  data=data |>
    dplyr::mutate(date=as.Date(time)) |>
    dplyr::group_by(date) |>
    dplyr::mutate(
      timeslot_id_diff=c(1, diff(timeslot_id)),
      gap_end=timeslot_id_diff > min_gap_size,
      gap_start=lead(gap_end))
  
  
  #after lead() last value of each group is NA
  #since using timeslot_id to detect gaps, can't know if last point is start of a gap
  #need to compare with total values
  i=which(is.na(data$gap_start))
  total_values=60*60*24/sample_rate #if sample_rate wrong, won't detect last gap 
  data$gap_start[i]=(total_values - data$timeslot_id[i]) > min_gap_size
  
  i=i[-length(i)]
  data$gap_end[i[which(data$gap_start[i])] + 1]=TRUE
  
  gap = (data$gap_start | data$gap_end)
  
  return(list(result=gap, calculation_parameters=list()))
}


calculate_exceeded_neighbours <- function(data, sample_rate, max_diff_either_neighbour=0.3, max_distance_to_fill=3, min_gap_size=1) {
  
  #ignore start and end point of gaps
  data$gap=calculate_gaps(data, sample_rate, min_gap_size)$result
  
  
  data$exceeded_neighbours=F
  
  
  #get diff from left to right
  exceeded_neighbours=c(F, abs(diff(data$slevel)) > max_diff_either_neighbour)
  
  i=which(exceeded_neighbours)
  if (length(i)==0) {
    return(list(result=exceeded_neighbours, calculation_parameters=list(gap=data$gap)))
  }
  
  #previous point is also suspicious
  exceeded_neighbours[i-1]=TRUE
  
  
  
  ##from each group of exceeded neighbours select right points
  i=which(exceeded_neighbours)
  group=rep(F, length(i))
  group[which(diff(i) > 1)+1]=T
  group=cumsum(group)+1
  
  groups=data.frame(group, slevel=data$slevel[i], index=i) |>
    dplyr::mutate(n=1:n()) |>
    dplyr::group_by(group) |>
    dplyr::mutate(
      group_size=n(),
      group_index=1:n()) 
  
  #assumption: with 3 it's always the middle one
  index_spike_of_triplets=groups |>
    dplyr::filter(
      group_size==3,
      group_index==2) |>
    dplyr::pull(index)
  
  #add points from the left and right of each group
  #abs(slevel - median of added points) > half of max_diff_either_neighbour are the outlier spikes
  index_spikes_of_group=dplyr::filter(groups, group_size != 3)
  if (nrow(index_spikes_of_group) != 0) {  
    index_spikes_of_group=index_spikes_of_group |>
      dplyr::slice(1) |>
      dplyr::mutate(p=list((index-10):(index+group_size-1+10))) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        p=list(data$slevel[ifelse(p < 1, NA, p)]),
        pl=length(p),
        pi=list(c(1:10,(pl-9):pl)),
        sp=list((abs(unlist(p) - median(unlist(p)[pi], na.rm=T)) > (max_diff_either_neighbour/2))[-pi]),
        index=list(c(index:(index+(group_size-1)))[sp])
      ) |>
      dplyr::pull(index)
  }
  
  i=c(index_spikes_of_group, index_spike_of_triplets) |> unlist()
  data$exceeded_neighbours[i]=TRUE
  
  
  
  #detect blunt spikes by filling in points between nearby exceeding neighbours
  i=which(data$exceeded_neighbours)
  blunt_spike_index=which(c(F, diff(i) <= (max_distance_to_fill+1)))
  start_pos=i[(blunt_spike_index-1)]
  end_pos=i[blunt_spike_index]
  i=Map(':',start_pos, end_pos) |> unlist()
  data$exceeded_neighbours[i]=TRUE
  
  data$exceeded_neighbours=(data$gap==FALSE & data$exceeded_neighbours==TRUE)
  
  return(list(result=data$exceeded_neighbours, calculation_parameters=list(gap=data$gap)))
}
