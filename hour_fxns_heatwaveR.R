require(heatwaveR)
require(data.table)
require(tidyverse)
require(ggplot2)
#ts2clm##############################################
ts2clm_posix <- function(data,
                   x = t,
                   y = temp,
                   climatologyPeriod,
                   robust = FALSE,
                   maxPadLength = FALSE,
                   windowHalfWidth = 5,
                   pctile = 90,
                   smoothPercentile = TRUE,
                   smoothPercentileWidth = 31,
                   clmOnly = FALSE,
                   var = FALSE,
                   roundClm = 4) {
  
  if (missing(climatologyPeriod))
    stop("Oops! Please provide a period (two dates) for calculating the climatology.")
  if (length(climatologyPeriod) != 2)
    stop("Bummer! Please provide BOTH start and end dates for the climatology period.")
  if (!(is.logical(robust)))
    stop("Please ensure that 'robust' is either TRUE or FALSE.")
  if (robust)
    message("The 'robust' argument has been deprecated and will be removed from future versions.")
  if (maxPadLength != FALSE & !is.numeric(maxPadLength))
    stop("Please ensure that 'maxPadLength' is either FALSE or a numeric/integer value.")
  if (!(is.numeric(pctile)))
    stop("Please ensure that 'pctile' is a numeric/integer value.")
  if (!(is.numeric(windowHalfWidth)))
    stop("Please ensure that 'windowHalfWidth' is a numeric/integer value.")
  if (!(is.logical(smoothPercentile)))
    stop("Please ensure that 'smoothPercentile' is either TRUE or FALSE.")
  if (!(is.numeric(smoothPercentileWidth)))
    stop("Please ensure that 'smoothPercentileWidth' is a numeric/integer value.")
  if (!(is.logical(clmOnly)))
    stop("Please ensure that 'clmOnly' is either TRUE or FALSE.")
  if (!(is.numeric(roundClm))) {
    if (!roundClm == FALSE) {
      stop("Please ensure that 'roundClm' is either a numeric value or FALSE.")
    }
  }
  
  clim_start <- climatologyPeriod[1]
  clim_end <- climatologyPeriod[2]
  temp <- doy <- .SD <-  NULL
  
  ts_x <- eval(substitute(x), data)
  if (is.null(ts_x) | is.function(ts_x))
    stop("Please ensure that a column named 't' is present in your data.frame or that you have assigned a column to the 'x' argument.")
  ts_y <- eval(substitute(y), data)
  if (is.null(ts_y) | is.function(ts_y))
    stop("Please ensure that a column named 'temp' is present in your data.frame or that you have assigned a column to the 'y' argument.")
  # rm(data) # Need to keep this for the end
  
  if (!inherits(ts_x[1], "POSIXct"))
    stop("Please ensure your date values are type 'POSIXct'. This may be done with 'as.POSIXct()'.")
  if (!is.numeric(ts_y[1]))
    stop("Please ensure the temperature values you are providing are type 'num' for numeric.")
  
  ts_xy <- data.table::data.table(ts_x = ts_x, ts_y = ts_y)[base::order(ts_x)]
  rm(list = c("ts_x", "ts_y"))
  
  ts_whole <- make_whole_fast_posix(ts_xy)
  
  if (length(stats::na.omit(ts_whole$ts_y)) < length(ts_whole$ts_y) & is.numeric(maxPadLength)) {
    ts_whole <- na_interp(doy = ts_whole$doy,
                          hoy = ts_whole$hoy,
                          x = ts_whole$ts_x,
                          y = ts_whole$ts_y,
                          maxPadLength = maxPadLength)
  }
  
  if (ts_whole$ts_x[1] > clim_start)
    stop(paste("The specified start date precedes the first day of series, which is",
               ts_whole$ts_x[1]))
  
  #if (clim_end > utils::tail(ts_whole$ts_x, 1))
   # stop(paste("The specified end date follows the last day of series, which is",
    #           ts_whole$ts_x[nrow(ts_whole)]))
  
  if (as.Date(clim_end) - as.Date(clim_start) < 1095)
    stop("The climatologyPeriod must be at least three years to calculate thresholds")
  
 
  
  ts_wide <- clim_spread_hour(ts_whole, clim_start, clim_end, windowHalfWidth)
  
  if (nrow(stats::na.omit(ts_wide)) < nrow(ts_wide) | var) {
    ts_mat <- clim_calc_hour(ts_wide, windowHalfWidth, pctile)
    ts_mat[is.nan(ts_mat)] <- NA
  } else {
    ts_mat <- clim_calc_hour(ts_wide, windowHalfWidth, pctile)
  }
  rm(ts_wide)
  
  if (smoothPercentile) {
    ts_clim <- smooth_percentile_hour(ts_mat, smoothPercentileWidth, var)
  } else {
    ts_clim <- data.table::data.table(ts_mat)
  }
  
  cols <- names(ts_clim)
  if (is.numeric(roundClm)) {
    ts_clim[,(cols) := round(.SD, roundClm), .SDcols = cols]
  }
  rm(ts_mat)
  
  if (clmOnly) {
    
    return(ts_clim)
    
  } else {
    
    data.table::setkey(ts_whole, hoy)
    data.table::setkey(ts_clim, hoy)
    ts_res <- merge(ts_whole, ts_clim, all = TRUE)
    rm(ts_whole); rm(ts_clim)
    data.table::setorder(ts_res, ts_x)
    names(ts_res)[3] <- paste(substitute(x))
    names(ts_res)[4] <- paste(substitute(y))
    
    if (ncol(data) > 2) {
      # It would be better to order the columns
      ts_res1 <- merge(data, ts_res, all = TRUE)
    }
    
    ts_res <- tibble::as_tibble(ts_res)
    return(ts_res)
    
  }
}


#make_whole_fast#####################

make_whole_fast_posix <- function(data) {
  
  feb28 <- 59
  
  # create full, complete time series for joining against
  date_start <- as.POSIXct(utils::head(data$ts_x, 1))
  date_end <- as.POSIXct(utils::tail(data$ts_x, 1))
  ts_full <- data.table::data.table(ts_x = seq.POSIXt(date_start, date_end, "hour"))
  ts_full2<-lubridate::ymd_hms(ts_full$ts_x)
  
  ts_merged <- merge(ts_full, data, all.x = TRUE)
  rm(ts_full)
  
  v_date <- ts_merged$ts_x
  v_doy <- lubridate::yday(v_date)
  v_doy <- as.integer(ifelse(
    lubridate::leap_year(lubridate::year(ts_merged$ts_x)) == FALSE,
    ifelse(v_doy > feb28, v_doy + 1, v_doy),
    v_doy)
  )
  v_hoy<-  hour(v_date) + (yday(v_date) - 1) * 24
  v_ts_y <- as.numeric(ts_merged$ts_y)
  
  t_series <- data.table::data.table(doy = v_doy,
                                     ts_x = v_date,
                                     ts_y = v_ts_y,
                                     hoy=v_hoy)
  rm(list = c("v_date", "v_doy", "v_ts_y","v_hoy"))
  
  return(t_series)
}

#clim_spread###################
clim_spread_hour <- function(data2, clim_start, clim_end, windowHalfWidth) {
  
  .NA2mean <- function(x) {
    z <- round(mean(x, na.rm = TRUE), 2)
    x[is.na(x)] <- z
    return(x)
  }
  
  ts_x <- ts_y <- NULL
  
  ts_clim <- data.table::as.data.table(data2)[ts_x %between% c(clim_start, clim_end)]
  
  rm(data2)
  
  data.table::setDT(ts_clim)[, ts_x := format(as.POSIXct(ts_x), "%Y") ]
  ts_spread <- data.table::dcast(ts_clim, hoy ~ ts_x, value.var = "ts_y", mean)
  rm(ts_clim)
  
  ts_spread_filled <- data.table::data.table((sapply(ts_spread[59:61, ],
                                                     function(x) .NA2mean(x))))
  ts_spread[60, ] <- ts_spread_filled[2, ]
  rm(ts_spread_filled)
  
  begin_pad <- utils::tail(ts_spread, windowHalfWidth)
  end_pad <- utils::head(ts_spread, windowHalfWidth)
  l <- list(begin_pad, ts_spread, end_pad)
  rm(list = c("begin_pad", "end_pad"))
  
  ts_spread <- data.table::rbindlist(l)
  rm(l)
  
  len_yr <- length(lubridate::year(clim_start):lubridate::year(clim_end))
  
  # clim_calc_cpp needs a matrix...
  ts_mat <- as.matrix(ts_spread)[, 2:(len_yr + 1)]
  
  if (nrow(stats::na.omit(ts_mat)) < nrow(ts_mat)) {
    plugs <- which(is.na(ts_mat), arr.ind = TRUE)
    ts_mat[plugs] <- rowMeans(ts_mat, na.rm = TRUE)[plugs[,1]]
  }
  
  return(ts_mat)
}

#clim_calc##############
clim_calc_hour <- function(data, windowHalfWidth, pctile) {
  
  seas <- rep(NA, nrow(data))
  thresh <- rep(NA, nrow(data))
  var <- rep(NA, nrow(data))
  
  for (i in (windowHalfWidth + 1):((nrow(data) - windowHalfWidth))) {
    seas[i] <-
      mean(
        c(t(data[(i - (windowHalfWidth)):(i + windowHalfWidth), seq_len(ncol(data))])),
        na.rm = TRUE)
    thresh[i] <-
      stats::quantile(
        c(t(data[(i - (windowHalfWidth)):(i + windowHalfWidth), seq_len(ncol(data))])),
        probs = pctile/100,
        type = 7,
        na.rm = TRUE,
        names = FALSE
      )
    var[i] <-
      stats::sd(
        c(t(data[(i - (windowHalfWidth)):(i + windowHalfWidth), seq_len(ncol(data))])),
        na.rm = TRUE
      )
  }
  
  len_clim_year <- 366*24
  hoy <- seq_len(366*24)
  
  seas <- seas[(windowHalfWidth + 1):((windowHalfWidth) + len_clim_year)]
  thresh <- thresh[(windowHalfWidth + 1):((windowHalfWidth) + len_clim_year)]
  var <- var[(windowHalfWidth + 1):((windowHalfWidth) + len_clim_year)]
  
  clim <- matrix(c(hoy, seas, thresh, var), ncol = 4, byrow = FALSE,
                 dimnames = list(NULL, c("hoy", "seas", "thresh", "var")))
  
  return(clim)
}

#smooth_percentil########
smooth_percentile_hour <- function(data, smoothPercentileWidth, var_calc) {
  
  seas <- thresh <- NULL
  
  prep <- rbind(utils::tail(data[,-1], smoothPercentileWidth),
                data[,-1],
                utils::head(data[,-1], smoothPercentileWidth))
  rm(data)
  
  len_clim_year <- 366*24
  
  seas <- RcppRoll::roll_mean(as.numeric(prep[,1]), n = smoothPercentileWidth, na.rm = FALSE)
  thresh <- RcppRoll::roll_mean(as.numeric(prep[,2]), n = smoothPercentileWidth, na.rm = FALSE)
  
  clim <- data.table::data.table(hoy = seq_len(len_clim_year),
                                 seas = seas[(smoothPercentileWidth/2 + 2):((smoothPercentileWidth/2 + 1) + len_clim_year)],
                                 thresh = thresh[(smoothPercentileWidth/2 + 2):((smoothPercentileWidth/2 + 1) + len_clim_year)])
  
  if (var_calc) {
    var <- NULL
    
    var <- RcppRoll::roll_mean(as.numeric(prep[,3]), n = smoothPercentileWidth, na.rm = FALSE)
    
    clim$var <- var[(smoothPercentileWidth/2 + 2):((smoothPercentileWidth/2 + 1) + len_clim_year)]
  }
  rm(prep)
  
  return(clim)
}

#detect_event################
detect_event_hour<-function (data, x = t, y = temp, seasClim = seas, threshClim = thresh, 
                             threshClim2 = NA, minDuration = 0, minDuration2 = minDuration, 
                             joinAcrossGaps = TRUE, maxGap = 2*24, maxGap2 = maxGap, coldSpells = FALSE, 
                             protoEvents = FALSE, categories = FALSE, roundRes = 4, ...) 
{
  if (!(is.numeric(minDuration))) 
    stop("Please ensure that 'minDuration' is a numeric/integer value.")
  if (!(is.logical(joinAcrossGaps))) 
    stop("Please ensure that 'joinAcrossGaps' is either TRUE or FALSE.")
  if (!(is.numeric(maxGap))) 
    stop("Please ensure that 'maxGap' is a numeric/integer value.")
  if (!(is.numeric(roundRes))) {
    if (!roundRes == FALSE) {
      stop("Please ensure that 'roundRes' is either a numeric value or FALSE.")
    }
  }
  temp <- seas <- thresh <- threshCriterion <- durationCriterion <- event <- NULL
  ts_x <- eval(substitute(x), data)
  if (is.null(ts_x) | is.function(ts_x)) 
    stop("Please ensure that a column named 't' is present in your data.frame or that you have assigned a column to the 'x' argument.")
  ts_y <- eval(substitute(y), data)
  if (is.null(ts_y) | is.function(ts_y)) 
    stop("Please ensure that a column named 'temp' is present in your data.frame or that you have assigned a column to the 'y' argument.")
  ts_seas <- eval(substitute(seasClim), data)
  if (is.null(ts_seas) | is.function(ts_seas)) 
    stop("Please ensure that a column named 'seas' is present in your data.frame or that you have assigned a column to the 'seasClim' argument.")
  ts_thresh <- eval(substitute(threshClim), data)
  if (is.null(ts_thresh) | is.function(ts_thresh)) 
    stop("Please ensure that a column named 'thresh' is present in your data.frame or that you have assigned a column to the 'threshClim' argument.")
  t_series <- data.frame(ts_x, ts_y, ts_seas, ts_thresh)
  rm(ts_x)
  rm(ts_y)
  rm(ts_seas)
  rm(ts_thresh)
  if (coldSpells) {
    t_series$ts_y <- -t_series$ts_y
    t_series$ts_seas <- -t_series$ts_seas
    t_series$ts_thresh <- -t_series$ts_thresh
  }
  t_series$ts_y[is.na(t_series$ts_y)] <- t_series$ts_seas[is.na(t_series$ts_y)]
  t_series$threshCriterion <- t_series$ts_y > t_series$ts_thresh
  t_series$threshCriterion[is.na(t_series$threshCriterion)] <- FALSE
  events_clim <- proto_event(t_series, criterion_column = t_series$threshCriterion, 
                             minDuration = minDuration, joinAcrossGaps = joinAcrossGaps, 
                             maxGap = maxGap)
  if (!is.na(threshClim2[1])) {
    if (!is.logical(threshClim2[1])) 
      stop("Please ensure 'threshClim2' contains logical values (e.g. TRUE and/or FALSE)")
    events_clim <- proto_event(t_series, criterion_column = events_clim$event & 
                                 threshClim2, minDuration = minDuration2, joinAcrossGaps = joinAcrossGaps, 
                               maxGap = maxGap2)
  }
  if (protoEvents) {
    events_clim <- data.frame(data, threshCriterion = events_clim$threshCriterion, 
                              durationCriterion = events_clim$durationCriterion, 
                              event = events_clim$event, event_no = events_clim$event_no)
    return(events_clim)
  }else {
    intensity_mean <- intensity_max <- intensity_cumulative <- intensity_mean_relThresh <- intensity_max_relThresh <- intensity_cumulative_relThresh <- intensity_mean_abs <- intensity_max_abs <- intensity_cumulative_abs <- rate_onset <- rate_decline <- mhw_rel_thresh <- mhw_rel_seas <- event_no <- row_index <- index_start <- index_peak <- index_end <- NULL
    if (nrow(stats::na.omit(events_clim)) > 0) {
      events <- data.frame(events_clim, row_index = base::seq_len(nrow(events_clim)), 
                           mhw_rel_seas = events_clim$ts_y - events_clim$ts_seas, 
                           mhw_rel_thresh = events_clim$ts_y - events_clim$ts_thresh)
      events <- events[stats::complete.cases(events$event_no), 
      ]
      events <- plyr::ddply(events, c("event_no"), .fun = plyr::summarise, 
                            index_start = min(row_index), 
                            index_peak = row_index[mhw_rel_seas == 
                              max(mhw_rel_seas)][1], 
                            index_end = max(row_index), 
                            duration = index_end - index_start + 1, date_start = min(ts_x), 
                            date_peak = ts_x[mhw_rel_seas == max(mhw_rel_seas)][1], 
                            date_end = max(ts_x), intensity_mean = mean(mhw_rel_seas), 
                            intensity_max = max(mhw_rel_seas), intensity_var = sqrt(stats::var(mhw_rel_seas)), 
                            intensity_cumulative = sum(mhw_rel_seas), intensity_mean_relThresh = mean(mhw_rel_thresh), 
                            intensity_max_relThresh = max(mhw_rel_thresh), 
                            intensity_var_relThresh = sqrt(stats::var(mhw_rel_thresh)), 
                            intensity_cumulative_relThresh = sum(mhw_rel_thresh), 
                            intensity_mean_abs = mean(ts_y), intensity_max_abs = max(ts_y), 
                            intensity_var_abs = sqrt(stats::var(ts_y)), intensity_cumulative_abs = sum(ts_y))
      events <- tibble::as_tibble(events)
      mhw_rel_seas <- t_series$ts_y - t_series$ts_seas
      A <- mhw_rel_seas[events$index_start]
      B <- t_series$ts_y[events$index_start - 1]
      C <- t_series$ts_seas[events$index_start - 1]
      if (length(B) + 1 == length(A)) {
        B <- c(NA, B)
        C <- c(NA, C)
      }
      mhw_rel_seas_start <- 0.5 * (A + B - C)
      events$rate_onset <- ifelse(events$index_start > 
                                    1, (events$intensity_max - mhw_rel_seas_start)/(as.numeric(difftime(events$date_peak, 
                                                                                                        events$date_start, units = "hours")) + 0.5), NA)
      D <- mhw_rel_seas[events$index_end]
      E <- t_series$ts_y[events$index_end + 1]
      G <- t_series$ts_seas[events$index_end + 1]
      mhw_rel_seas_end <- 0.5 * (D + E - G)
      events$rate_decline <- ifelse(events$index_end < 
                                      nrow(t_series), (events$intensity_max - mhw_rel_seas_end)/(as.numeric(difftime(events$date_end, 
                                                                                                                     events$date_peak, units = "hours")) + 0.5), NA)
      if (coldSpells) {
        events$intensity_mean <- -events$intensity_mean
        events$intensity_max <- -events$intensity_max
        events$intensity_cumulative <- -events$intensity_cumulative
        events$intensity_mean_relThresh <- -events$intensity_mean_relThresh
        events$intensity_max_relThresh <- -events$intensity_max_relThresh
        events$intensity_cumulative_relThresh <- -events$intensity_cumulative_relThresh
        events$intensity_mean_abs <- -events$intensity_mean_abs
        events$intensity_max_abs <- -events$intensity_max_abs
        events$intensity_cumulative_abs <- -events$intensity_cumulative_abs
        events$rate_onset <- -events$rate_onset
        events$rate_decline <- -events$rate_decline
      }
    }else {
      events <- data.frame(event_no = NA, index_start = NA, 
                           index_peak = NA, index_end = NA, duration = NA, 
                           date_start = NA, date_peak = NA, date_end = NA, 
                           intensity_mean = NA, intensity_max = NA, intensity_var = NA, 
                           intensity_cumulative = NA, intensity_mean_relThresh = NA, 
                           intensity_max_relThresh = NA, intensity_var_relThresh = NA, 
                           intensity_cumulative_relThresh = NA, intensity_mean_abs = NA, 
                           intensity_max_abs = NA, intensity_var_abs = NA, 
                           intensity_cumulative_abs = NA, rate_onset = NA, 
                           rate_decline = NA)
      events <- tibble::as_tibble(events)
    }
    event_cols <- names(events)[9:22]
    clim_cols <- names(events_clim)[2:4]
    if (nrow(events) == 1) {
      if (is.na(events$rate_onset)) {
        event_cols <- event_cols[-grep(pattern = "rate_onset", 
                                       x = event_cols, value = FALSE)]
      }
      if (is.na(events$rate_decline)) {
        event_cols <- event_cols[-grep(pattern = "rate_decline", 
                                       x = event_cols, value = FALSE)]
      }
    }
    if (is.numeric(roundRes)) {
      if (nrow(events) > 0) {
        events <- dplyr::mutate(events, dplyr::across(dplyr::all_of(event_cols), 
                                                      round, roundRes))
        events_clim <- dplyr::mutate(events_clim, dplyr::across(dplyr::all_of(clim_cols), 
                                                                round, roundRes))
      }
    }
    data_clim <- tibble::as_tibble(cbind(data, events_clim[, 
                                                           5:8]))
    data_res <- list(climatology = data_clim, event = events)
    if (categories) {
      data_cat <- category(data_res, ...)
      if (is.data.frame(data_cat)) {
        data_res <- dplyr::left_join(events, data_cat, 
                                     by = c("event_no", "duration", intensity_max = "i_max", 
                                            date_peak = "peak_date"))
      }else {
        data_res <- list(climatology = dplyr::left_join(data_res$climatology, 
                                                        data_cat$climatology, by = c("t", "event_no")), 
                         event = dplyr::left_join(data_res$event, data_cat$event, 
                                                  by = c("event_no", "duration", intensity_max = "i_max", 
                                                         date_peak = "peak_date")))
      }
    }
    return(data_res)
  }
}
#event_line hour####################
event_line_hour<-function (data, x = t, y = temp, min_duration = 0, spread = 150, 
                           metric = "intensity_cumulative", start_date = NULL, end_date = NULL, 
                           category = FALSE, x_axis_title = NULL, x_axis_text_angle = NULL, 
                           y_axis_title = NULL, y_axis_range = NULL) 
{
  date_end <- date_start <- duration <- temp <- NULL
  if (!(exists("event", data)) | !(exists("climatology", data))) 
    stop("Please ensure you are running this function on the output of 'heatwaveR::detect_event()'")
  ts_x <- eval(substitute(x), data$climatology)
  data$climatology$ts_x <- ts_x
  ts_y <- eval(substitute(y), data$climatology)
  data$climatology$ts_y <- ts_y
  if (is.null(start_date)) 
    start_date <- min(data$climatology$ts_x)
  if (is.null(end_date)) 
    end_date <- max(data$climatology$ts_x)
  event <- data$event %>% dplyr::filter(date_end >= start_date & 
                                          date_start <= end_date) %>% data.frame()
  if (nrow(event) == 0) 
    stop("No events detected! Consider changing the 'start_date' or 'end_date' values.")
  if (!(metric %in% c("intensity_mean", "intensity_max", "intensity_var", 
                      "intensity_cumulative", "intensity_mean_relThresh", "intensity_max_relThresh", 
                      "intensity_var_relThresh", "intensity_cumulative_relThresh", 
                      "intensity_mean_abs", "intensity_max_abs", "intensity_var_abs", 
                      "intensity_cumulative_abs", "rate_onset", "rate_decline"))) {
    stop("Please ensure you have spelled the desired metric correctly.")
  }
  index_start <- index_end <- event_idx <- NULL
  event_idx <- as.vector(-abs(event[colnames(event) == metric])[, 
                                                                1])
  event <- event[base::order(event_idx), ]
  event <- event %>% dplyr::filter(duration >= min_duration) %>% 
    dplyr::mutate(index_start_fix = index_start - 1, index_end_fix = index_end + 
                    1)
  event_top <- event[1, ]
  date_spread <- seq((event_top$date_start - spread), (event_top$date_end + 
                                                         spread), by = "hour")
  event_sub <- event %>% dplyr::filter(date_start >= min(date_spread), 
                                       date_end <= max(date_spread))
  thresh_2x <- thresh_3x <- thresh_4x <- NULL
  clim_diff <- data$climatology %>% dplyr::mutate(diff = thresh - 
                                                    seas, thresh_2x = thresh + diff, thresh_3x = thresh_2x + 
                                                    diff, thresh_4x = thresh_3x + diff)
  clim_events <- data.frame()
  for (i in seq_len(nrow(event_sub))) {
    clim_sub <- clim_diff[(event_sub$index_start_fix[i]):(event_sub$index_end_fix[i]),]
    clim_events <- rbind(clim_events, clim_sub)
  }
  clim_top <- clim_diff[event_top$index_start_fix:event_top$index_end_fix, 
  ]
  clim_spread <- clim_diff %>% dplyr::filter(ts_x %in% date_spread)
  thresh <- seas <- y1 <- y2 <- NULL
  if (event_top$intensity_mean > 0) {
    fillCol <- c(events = "salmon", `peak event` = "red")
    clim_events$y1 <- clim_events$ts_y
    clim_events$y2 <- clim_events$thresh
    clim_top$y1 <- clim_top$ts_y
    clim_top$y2 <- clim_top$thresh
  }else {
    fillCol <- c(events = "steelblue3", `peak event` = "navy")
    clim_events$y1 <- clim_events$thresh
    clim_events$y2 <- clim_events$ts_y
    clim_top$y1 <- clim_top$thresh
    clim_top$y2 <- clim_top$ts_y
  }
  if (!is.null(y_axis_title)) {
    if (!is.character(y_axis_title)) 
      stop("Please ensure that the argument provided to 'y_axis_title' is a character string.")
    ylabel <- y_axis_title
  }else {
    ylabel <- expression(paste("Temperature [", degree, "C]"))
  }
  if (!is.null(x_axis_title)) {
    if (!is.character(x_axis_title)) 
      stop("Please ensure that the argument provided to 'x_axis_title' is a character string.")
    xlabel <- x_axis_title
  }else {
    xlabel <- NULL
  }
  if (!is.null(x_axis_text_angle)) {
    if (!is.numeric(x_axis_text_angle)) 
      stop("Please ensure that the argument provided to 'x_axis_text_angle' is a number.")
    xtangle <- x_axis_text_angle
  }else {
    xtangle <- 0
  }
  ep <- ggplot(data = clim_spread, aes(x = ts_x, y = ts_y)) + 
    scale_x_datetime(expand = c(0, 0), date_labels = "%b %Y") + 
    labs(x = xlabel, y = ylabel) + theme(plot.background = element_blank(), 
                                         panel.background = element_rect(fill = "white"), panel.border = element_rect(colour = "black", 
                                                                                                                      fill = NA, size = 0.75), panel.grid.minor = element_line(colour = NA), 
                                         panel.grid.major = element_line(colour = "black", size = 0.2, 
                                                                         linetype = "dotted"), axis.text = element_text(colour = "black"), 
                                         axis.text.x = element_text(margin = unit(c(0.5, 0.5, 
                                                                                    0.5, 0.5), "cm"), angle = xtangle), axis.text.y = element_text(margin = unit(c(0.5, 
                                                                                                                                                                   0.5, 0.5, 0.5), "cm")), axis.ticks.length = unit(-0.25, 
                                                                                                                                                                                                                    "cm"), legend.background = element_rect(colour = "black"), 
                                         legend.direction = "horizontal", legend.justification = c(0, 
                                                                                                   0), legend.position = c(0.005, 0.015), legend.key = element_blank())
  if (category) {
    lineColCat <- c(Temperature = "black", Climatology = "gray20", 
                    Threshold = "darkgreen", `2x Threshold` = "darkgreen", 
                    `3x Threshold` = "darkgreen", `4x Threshold` = "darkgreen")
    if (event_top$intensity_mean < 0) {
      fillColCat <- c(Moderate = "#C7ECF2", Strong = "#85B7CC", 
                      Severe = "#4A6A94", Extreme = "#111433")
      ep <- ep + geom_flame(data = clim_events, size = 0.5, 
                            aes(x = ts_x, y = thresh, y2 = ts_y, fill = "Moderate")) + 
        geom_flame(data = clim_events, size = 0.5, aes(x = ts_x, 
                                                       y = thresh_2x, y2 = ts_y, fill = "Strong")) + 
        geom_flame(data = clim_events, size = 0.5, aes(x = ts_x, 
                                                       y = thresh_3x, y2 = ts_y, fill = "Severe")) + 
        geom_flame(data = clim_events, size = 0.5, aes(x = ts_x, 
                                                       y = thresh_4x, y2 = ts_y, fill = "Extreme"))
    }else {
      fillColCat <- c(Moderate = "#ffc866", Strong = "#ff6900", 
                      Severe = "#9e0000", Extreme = "#2d0000")
      ep <- ep + geom_flame(data = clim_events, size = 0.5, 
                            aes(x = ts_x, y = y1, y2 = y2, fill = "Moderate")) + 
        geom_flame(data = clim_events, size = 0.5, aes(x = ts_x, 
                                                       y = y1, y2 = thresh_2x, fill = "Strong")) + 
        geom_flame(data = clim_events, size = 0.5, aes(x = ts_x, 
                                                       y = y1, y2 = thresh_3x, fill = "Severe")) + 
        geom_flame(data = clim_events, size = 0.5, aes(x = ts_x, 
                                                       y = y1, y2 = thresh_4x, fill = "Extreme"))
    }
    ep <- ep + geom_line(aes(y = thresh_2x, col = "2x Threshold"), 
                         size = 0.7, linetype = "dashed") + geom_line(aes(y = thresh_3x, 
                                                                          col = "3x Threshold"), size = 0.7, linetype = "dotdash") + 
      geom_line(aes(y = thresh_4x, col = "4x Threshold"), 
                size = 0.7, linetype = "dotted") + geom_line(aes(y = seas, 
                                                                 col = "Climatology"), size = 0.7, alpha = 1) + geom_line(aes(y = thresh, 
                                                                                                                              col = "Threshold"), size = 0.7, alpha = 1) + geom_line(aes(y = ts_y, 
                                                                                                                                                                                         col = "Temperature"), size = 0.6) + scale_colour_manual(name = NULL, 
                                                                                                                                                                                                                                                 values = lineColCat, breaks = c("Temperature", "Climatology", 
                                                                                                                                                                                                                                                                                 "Threshold", "2x Threshold", "3x Threshold", 
                                                                                                                                                                                                                                                                                 "4x Threshold")) + scale_fill_manual(name = NULL, 
                                                                                                                                                                                                                                                                                                                      values = fillColCat, guide = "none") + guides(colour = guide_legend(override.aes = list(linetype = c("solid", 
                                                                                                                                                                                                                                                                                                                                                                                                                           "solid", "solid", "dashed", "dotdash", "dotted"), 
                                                                                                                                                                                                                                                                                                                                                                                                              size = c(0.6, 0.7, 0.7, 0.7, 0.7, 0.7)))) + theme(legend.direction = "vertical")
    ep
  }
  else {
    lineCol <- c(Temperature = "black", Climatology = "blue", 
                 Threshold = "darkgreen")
    ep <- ep + geom_flame(data = clim_events, size = 0.5, 
                          aes(x = ts_x, y = y1, y2 = y2, fill = "events")) + 
      geom_flame(data = clim_top, size = 0.5, aes(x = ts_x, 
                                                  y = y1, y2 = y2, fill = "peak event")) + geom_line(aes(y = seas, 
                                                                                                         col = "Climatology"), size = 0.7, alpha = 1) + geom_line(aes(y = thresh, 
                                                                                                                                                                      col = "Threshold"), size = 0.7, alpha = 1) + geom_line(aes(y = ts_y, 
                                                                                                                                                                                                                                 col = "Temperature"), size = 0.6) + scale_colour_manual(name = NULL, 
                                                                                                                                                                                                                                                                                         values = lineCol, breaks = c("Temperature", "Climatology", 
                                                                                                                                                                                                                                                                                                                      "Threshold")) + scale_fill_manual(name = NULL, 
                                                                                                                                                                                                                                                                                                                                                        values = fillCol, guide = FALSE)
    if (!is.null(y_axis_range)) {
      if (length(y_axis_range) != 2) 
        stop("Please ensure that exactly two numbers are provided to 'y_axis_range' (e.g. c(10, 20)).")
      if (!is.numeric(y_axis_range[1]) | !is.numeric(y_axis_range[2])) 
        stop("Please ensure that only numeric values are provided to 'y_axis_range'.")
      ep <- ep + coord_cartesian(ylim = c(y_axis_range[1], 
                                          y_axis_range[2]))
    }
    ep
  }
}

#proto_event#########
proto_event <- function(t_series,
                       criterion_column,
                       minDuration,
                       joinAcrossGaps,
                       maxGap) {
  
  index_start <- index_end <- duration <- NULL
  
  ex1 <- rle(criterion_column)
  ind1 <- rep(seq_along(ex1$lengths), ex1$lengths)
  s1 <- base::split(base::seq_len(nrow(t_series)), ind1)
  
  proto_events <- do.call(rbind, lapply(s1[ex1$values == TRUE], function(x)
    data.frame(index_start = min(x), index_end = max(x))))
  
  duration <- proto_events$index_end - proto_events$index_start + 1
  
  suppressWarnings(
    if (is.null(proto_events) | max(duration) < minDuration){
      res <- data.frame(t_series,
                        durationCriterion = FALSE,
                        event = FALSE,
                        event_no = NA)
      return(res)
    } else {
      proto_events$duration <- duration
    }
  )
  
  proto_events <- proto_events[proto_events$duration >= minDuration, ]
  
  # NB: Apparently using for loops on pre-allocated memory size vectors is faster than using *apply()
  # https://johanndejong.wordpress.com/2016/07/07/r-are-apply-loops-faster-than-for-loops/
  # Meaning this is arguably one of the faster ways to do this
  durationCriterion <- rep(FALSE, nrow(t_series))
  for (i in base::seq_len(nrow(proto_events))) {
    durationCriterion[proto_events$index_start[i]:proto_events$index_end[i]] <-
      rep(TRUE, length = proto_events$duration[i])
  }
  
  if (joinAcrossGaps) {
    ex2 <- rle(durationCriterion)
    ind2 <- rep(seq_along(ex2$lengths), ex2$lengths)
    s2 <- base::split(base::seq_len(nrow(t_series)), ind2)
    
    proto_gaps <- do.call(rbind, lapply(s2[ex2$values == FALSE], function(x)
      data.frame(index_start = min(x), index_end = max(x))))
    proto_gaps$duration <- proto_gaps$index_end - proto_gaps$index_start + 1
    proto_gaps <- proto_gaps[proto_gaps$index_end > proto_events$index_start[1], ]
    
    if (any(proto_gaps$duration >= 1 & proto_gaps$duration <= maxGap)) {
      proto_gaps <- proto_gaps[proto_gaps$duration >= 1 & proto_gaps$duration <= maxGap, ]
      
      event <- durationCriterion
      for (i in base::seq_len(nrow(proto_gaps))) {
        event[proto_gaps$index_start[i]:proto_gaps$index_end[i]] <-
          rep(TRUE, length = proto_gaps$duration[i])
      }
      
    } else {
      
      event <- durationCriterion
      
    }
    
  } else {
    
    event <- durationCriterion
    
  }
  
  ex3 <- rle(event)
  ind3 <- rep(seq_along(ex3$lengths), ex3$lengths)
  s3 <- base::split(base::seq_len(nrow(t_series)), ind3)
  
  proto_final <- do.call(rbind, lapply(s3[ex3$values == TRUE], function(x)
    data.frame(index_start = min(x), index_end = max(x))))
  proto_final$duration <- proto_final$index_end - proto_final$index_start + 1
  
  event_no <- rep(NA, nrow(t_series))
  for (i in base::seq_len(nrow(proto_final))) {
    event_no[proto_final$index_start[i]:proto_final$index_end[i]] <-
      rep(i, length = proto_final$duration[i])
  }
  
  res <- cbind(t_series, durationCriterion, event, event_no)
  
  return(res)
  
}


#
#category_fixed
category_fixed<-function (data, y = temp, S = TRUE, name = "Event", climatology = FALSE, 
          MCScorrect = F, season = "range", roundVal = 4) 
{
  temp <- NULL
  ts_y <- eval(substitute(y), data$climatology)
  data$climatology$ts_y <- ts_y
  rm(ts_y)
  event_no <- event_name <- peak_date <- category <- duration <- NULL
  if (nrow(stats::na.omit(data$event)) == 0) {
    cat_res <- tibble::as_tibble(data.frame(event_no = NA, 
                                            event_name = NA, peak_date = NA, category = NA, i_max = NA, 
                                            duration = NA, p_moderate = NA, p_strong = NA, p_severe = NA, 
                                            p_extreme = NA, season = NA))
    if (climatology) {
      clim_res <- tibble::as_tibble(data.frame(t = NA, 
                                               event_no = NA, intensity = NA, category = NA))
      res <- list(climatology = clim_res, event = cat_res)
      return(res)
    }
    else {
      return(cat_res)
    }
  }
  cat_frame <- data.frame(event_no = data$event$event_no, event_name = paste0(as.character(name), 
                                                                              " ", lubridate::year(data$event$date_peak)), peak_date = data$event$date_peak, 
                          category = NA, i_max = round(data$event$intensity_max, 
                                                       roundVal), duration = data$event$duration)
  
  #cat_frame<-cat_frame%>%drop_na(peak_date)
  seasons <- data.frame(event_no = data$event$event_no, date_start = data$event$date_start, 
                        date_peak = data$event$date_peak, date_end = data$event$date_end, 
                        duration = data$event$duration, season = NA)
  ss <- as.POSIXlt(data$event$date_start)
  ss$day <- 1
  ss$mo <- ss$mo + 1
  sp <- as.POSIXlt(data$event$date_peak)
  sp$day <- 1
  sp$mo <- sp$mo + 1
  se <- as.POSIXlt(data$event$date_end)
  se$day <- 1
  se$mo <- se$mo + 1
  start_season <- peak_season <- end_season <- NULL
  if (S) {
    seasons$start_season <- factor(quarters(ss, abbreviate = F), 
                                   levels = c("Q1", "Q2", "Q3", "Q4"), labels = c("Summer", 
                                                                                  "Fall", "Winter", "Spring"))
    seasons$peak_season <- factor(quarters(sp, abbreviate = F), 
                                  levels = c("Q1", "Q2", "Q3", "Q4"), labels = c("Summer", 
                                                                                 "Fall", "Winter", "Spring"))
    seasons$end_season <- factor(quarters(se, abbreviate = F), 
                                 levels = c("Q1", "Q2", "Q3", "Q4"), labels = c("Summer", 
                                                                                "Fall", "Winter", "Spring"))
  }else {
    seasons$start_season <- factor(quarters(ss, abbreviate = F), 
                                   levels = c("Q1", "Q2", "Q3", "Q4"), labels = c("Winter", 
                                                                                  "Spring", "Summer", "Fall"))
    seasons$peak_season <- factor(quarters(sp, abbreviate = F), 
                                  levels = c("Q1", "Q2", "Q3", "Q4"), labels = c("Winter", 
                                                                                 "Spring", "Summer", "Fall"))
    seasons$end_season <- factor(quarters(se, abbreviate = F), 
                                 levels = c("Q1", "Q2", "Q3", "Q4"), labels = c("Winter", 
                                                                                "Spring", "Summer", "Fall"))
  }
  if (season == "range") {
    seasons <- seasons %>% dplyr::mutate(diff_season = as.integer(start_season) - 
                                           as.integer(end_season))
    for (i in 1:nrow(seasons)) {
      if (seasons$diff_season[i] == 0 & seasons$duration[i] < 
          100) {
        seasons$season[i] <- paste0(seasons$start_season[i])
      }
      else if (seasons$diff_season[i] %in% c(-1, 3) & seasons$duration[i] < 
               180) {
        seasons$season[i] <- paste0(seasons$start_season[i], 
                                    "/", seasons$end_season[i])
      }
      else if (seasons$diff_season[i] %in% c(-3, -2, 2)) {
        seasons$season[i] <- paste0(seasons$start_season[i], 
                                    "-", seasons$end_season[i])
      }
      if (seasons$duration[i] > 270) {
        seasons$season[i] <- "Year-round"
      }
    }
  }else if (season == "start") {
    seasons$season <- as.character(seasons$start_season)
  }else if (season == "peak") {
    seasons$season <- as.character(seasons$peak_season)
  }else if (season == "end") {
    seasons$season <- as.character(seasons$end_season)
  }else {
    stop("Please provide one of the following to the `season` argument: 'range', 'start', 'peak', 'end'.")
  }
  seas <- thresh <- thresh_2x <- thresh_3x <- thresh_4x <- NULL
  clim_diff <- data$climatology %>% dplyr::filter(!is.na(event_no)) %>% 
    dplyr::mutate(diff = round(thresh - seas, roundVal), 
                  thresh_2x = round(thresh + diff, roundVal), thresh_3x = round(thresh_2x + 
                                                                                  diff, roundVal), thresh_4x = round(thresh_3x + 
                                                                                                                       diff, roundVal))
  if (MCScorrect) {
    clim_diff <- clim_diff %>% dplyr::mutate(diff = round(dplyr::case_when(thresh_4x + 
                                                                             diff <= -1.8 ~ -(thresh + 1.8)/4, TRUE ~ diff), roundVal), 
                                             thresh_2x = round(thresh + diff, roundVal), thresh_3x = round(thresh_2x + 
                                                                                                             diff, roundVal), thresh_4x = round(thresh_3x + 
                                                                                                                                                  diff, roundVal))
  }
  moderate <- strong <- severe <- extreme <- NULL
  if (max(cat_frame$i_max) < 0) {
    clim_diff$ts_y <- -clim_diff$ts_y
    clim_diff$seas <- -clim_diff$seas
    clim_diff$thresh <- -clim_diff$thresh
    clim_diff$thresh_2x <- -clim_diff$thresh_2x
    clim_diff$thresh_3x <- -clim_diff$thresh_3x
    clim_diff$thresh_4x <- -clim_diff$thresh_4x
  }
  moderate_n <- clim_diff %>% dplyr::filter(ts_y > thresh) %>% 
    dplyr::group_by(event_no) %>% dplyr::summarise(moderate = dplyr::n(), 
                                                   .groups = "drop")
  strong_n <- clim_diff %>% dplyr::filter(ts_y > thresh_2x) %>% 
    dplyr::group_by(event_no) %>% dplyr::summarise(strong = dplyr::n(), 
                                                   .groups = "drop")
  severe_n <- clim_diff %>% dplyr::filter(ts_y > thresh_3x) %>% 
    dplyr::group_by(event_no) %>% dplyr::summarise(severe = dplyr::n(), 
                                                   .groups = "drop")
  extreme_n <- clim_diff %>% dplyr::filter(ts_y > thresh_4x) %>% 
    dplyr::group_by(event_no) %>% dplyr::summarise(extreme = dplyr::n(), 
                                                   .groups = "drop")
  cat_n <- dplyr::left_join(moderate_n, strong_n, by = "event_no") %>% 
    dplyr::left_join(severe_n, by = "event_no") %>% dplyr::left_join(extreme_n, 
                                                                     by = "event_no")
  cat_n[is.na(cat_n)] <- 0
  p_moderate <- p_strong <- p_severe <- p_extreme <- event_count <- event_name_letter <- NULL
  cat_join <- dplyr::left_join(cat_frame, cat_n, by = "event_no") %>% 
    dplyr::mutate(p_moderate = round(((moderate - strong)/duration * 
                                        100), 0), p_strong = round(((strong - severe)/duration * 
                                                                      100), 0), p_severe = round(((severe - extreme)/duration * 
                                                                                                    100), 0), p_extreme = round((extreme/duration * 100), 
                                                                                                                                0), category = ifelse(p_extreme > 0, "IV Extreme", 
                                                                                                                                                      ifelse(p_severe > 0, "III Severe", ifelse(p_strong > 
                                                                                                                                                                                                  0, "II Strong", "I Moderate"))), event_name = replace(event_name, 
                                                                                                                                                                                                                                                        which(category == "I Moderate"), NA), event_name = as.character(event_name)) %>% 
    dplyr::arrange(event_no) %>% dplyr::left_join(seasons[, 
                                                          c("event_no", "season")], by = "event_no") %>% dplyr::select(event_no:duration, 
                                                                                                                       p_moderate:season) %>% droplevels() %>% dplyr::group_by(event_name) %>% 
    dplyr::mutate(event_count = dplyr::row_number()) %>% 
    dplyr::mutate(event_name_letter = dplyr::case_when(max(event_count) > 
                                                         1 & !is.na(event_name) ~ letters[event_count])) %>% 
    dplyr::ungroup() %>% dplyr::mutate(event_name = dplyr::case_when(!is.na(event_name_letter) ~ 
                                                                       paste0(event_name, event_name_letter), TRUE ~ event_name), 
                                       event_name = as.factor(event_name)) %>% dplyr::select(-event_count, 
                                                                                             -event_name_letter)
  cat_res <- tibble::as_tibble(cat_join) %>% dplyr::arrange(-p_moderate, 
                                                            -p_strong, -p_severe, -p_extreme)
  if (climatology) {
    doy <- intensity <- NULL
    clim_res <- clim_diff %>% dplyr::mutate(category = ifelse(ts_y > 
                                                                thresh_4x, "IV Extreme", ifelse(ts_y > thresh_3x, 
                                                                                                "III Severe", ifelse(ts_y > thresh_2x, "II Strong", 
                                                                                                                     ifelse(ts_y > thresh, "I Moderate", NA)))), intensity = round(ts_y - 
                                                                                                                                                                                     seas, roundVal)) %>% dplyr::select(t, event_no, intensity, 
                                                                                                                                                                                                                        category)
    if (max(cat_frame$i_max) < 0) 
      clim_res$intensity <- -clim_res$intensity
    list(climatology = tibble::as_tibble(clim_res), event = cat_res)
  }else {
    return(cat_res)
  }
}
