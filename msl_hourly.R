calculate_msl_hourly <- function(data) {
  
  mslh=dplyr::mutate(data, date=as.Date(time), hour=as.integer(format(time,"%H"))) |> 
    dplyr::group_by(date, hour) |>
    dplyr::filter(!is.na(date)) |>
    dplyr::mutate(msl_hourly=median(slevel, na.rm=TRUE)) |>
    dplyr::slice(1) 
  
  mslh=expand.grid(date=unique(mslh$date),hour=0:23) |>
    dplyr::left_join(mslh, by=c("date","hour")) |>
    dplyr::group_by(date) |>
    dplyr::mutate(msl_hourly=list(msl_hourly)) |>
    dplyr::slice(1)
  
  mslh=sapply(mslh$msl_hourly,function(m) {
    append(rep(NA,24-length(m)),m) #add NAs from the back, to get hour 0-23
  })
  mslh=lapply(seq_len(ncol(mslh)), function(i) mslh[,i]) #columns to list
  
  return(list(result = mslh, calculation_parameters = list()))
}
