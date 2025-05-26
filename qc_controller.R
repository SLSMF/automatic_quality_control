
source("src/qc_scripts/distinctness.R")
source("src/qc_scripts/exceeded_neighbours.R")
source("src/qc_scripts/flat_line.R")
source("src/qc_scripts/msl_daily_doodson.R")
source("src/qc_scripts/msl_daily_harmonics.R")
source("src/qc_scripts/msl_hourly.R")
source("src/qc_scripts/out_of_range.R")
source("src/qc_scripts/polyfitNoTests.R")
source("src/qc_scripts/shift.R")
source("src/qc_scripts/spikes_via_median.R")
source("src/qc_scripts/spikes_via_rmse.R")

create_indexes <- function(dbs, data, boolean_parameters, suffix="_index") {
  #empty index column means all were FALSE, unless for that day all flat_lines=T or less than 24 values, then that index was not calculated
  
  index_columns=paste0(boolean_parameters, suffix)
  
  #if data not present, give empty index columns
  if (!is.data.frame(data) || nrow(data)==0) {
    dbs[index_columns]=list(rep(list(integer(0)), nrow(dbs)))
    return(dbs)
  }
  
  
  missing_columns=setdiff(boolean_parameters, colnames(data))
  data[missing_columns]=F
  
  
  data=dplyr::mutate(data, "date"=as.Date(time)) |> 
      dplyr::group_by(date)
  
  indexes=list()
  for (n in 1:length(boolean_parameters)) {
    indexes[[n]]=dplyr::summarise(data, !!index_columns[n] := list(timeslot_id[which(get(boolean_parameters[n]))])) |>
      dplyr::select(-date)
  }
  
  #add new columns, but if already exists, replace it
  dbs=dplyr::mutate(dbs, !!!indexes)
  
  return(dbs)
}

create_values_after_qc <- function(dbs, index_parameters_to_remove) {
  
  dbs=dplyr::rowwise(dbs)
  
  if (any(index_parameters_to_remove=="remove_all")) {
    dbs=dplyr::mutate(dbs,
      "values_after_qc"=list(rep(NA, list(length(unlist(`_values`))))),
      "n_values_after_qc"=0)
  } else {
    dbs=dbs |>
      dplyr::mutate(
        "indexes_to_remove" = list(unname(unlist(across(all_of(index_parameters_to_remove))))), 
        "values_after_qc" = list(replace(`_values`, indexes_to_remove, NA)),
        "n_values_after_qc" = na.omit(values_after_qc) |> length(),
        "completeness_after_qc" = n_values_after_qc / `_n_values`) |>
      dplyr::ungroup() |>
      dplyr::select(-any_of(c("indexes_to_remove","remove_all")))
  }
  
  return(dbs)
}

split_at <- function(vector, pos) {
  out <- list()
  pos2 <- c(1, pos, length(vector)+1)
  for (i in seq_along(pos2[-1])) {
    out[[i]] <- vector[pos2[i]:(pos2[i+1]-1)]
  }
  return(out)
}

