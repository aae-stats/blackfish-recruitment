# create a system lookup to match CPUE and longnames
waterbody_lu <- list(
  "broken_creek_r4" = list(waterbody = "Broken Creek", reach_no = 4), 
  "glenelg_river_r1" = list(waterbody = "Glenelg River", reach_no = 1),
  "glenelg_river_r2" = list(waterbody = "Glenelg River", reach_no = 2), 
  "glenelg_river_r3" = list(waterbody = "Glenelg River", reach_no = 3), 
  "goulburn_river_r4" = list(waterbody = "Goulburn River", reach_no = 4), 
  "holland_creek_r1" = list(waterbody = "Holland Creek", reach_no = 1), 
  "hughes_creek_r1" = list(waterbody = "Hughes Creek", reach_no = 1), 
  "kiewa_river_r1" = list(waterbody = "Kiewa River", reach_no = 1), 
  "kiewaeast_river_r1" = list(waterbody = "Kiewa River East Branch", reach_no = 1),
  "kingparrot_creek_r1" = list(waterbody = "King Parrot Creek", reach_no = 1),
  "moorabool_river_r3" = list(waterbody = "Moorabool River", reach_no = 3), 
  "moorabool_river_r4" = list(waterbody = "Moorabool River", reach_no = 4),
  "murray_river_r6" = list(waterbody = "Murray River", reach_no = 6),
  "ovens_river_r0" = list(waterbody = "Ovens River", reach_no = 0), 
  "ovens_river_r4" = list(waterbody = "Ovens River", reach_no = 4),
  "ovens_river_r5" = list(waterbody = "Ovens River", reach_no = 5),
  "seven_creeks_r1" = list(waterbody = "Seven Creeks", reach_no = 1),
  "seven_creeks_r2" = list(waterbody = "Seven Creeks", reach_no = 2),
  "thomson_river_r2" = list(waterbody = "Thomson River", reach_no = 2), 
  "thomson_river_r3" = list(waterbody = "Thomson River", reach_no = 3)
)

# wrapper function to calculate all metrics
calculate_metrics <- function(x, system, baseline = NULL, level = 0.1) {
  
  # set baseline equal to x if not provided
  if (is.null(baseline))
    baseline <- x$stream_discharge_mld
  
  # calculate quantile threshold over all baseline data
  threshold <- quantile(baseline, probs = level)
  
  # and calculate rescale value of all baseline data
  lt_median <- median(baseline)
  
  # calculate and return all flow metrics
  data.frame(
    waterbody = waterbody_lu[[system]]$waterbody,
    reach_no = waterbody_lu[[system]]$reach_no,
    water_year = calculate_qx(
      x$stream_discharge_mld,
      x$date_formatted,
      season = 7:18, 
      threshold = threshold,
      output = "date"
    ),
    lowflow_wateryear = calculate_qx(
      x$stream_discharge_mld,
      x$date_formatted,
      season = 7:18, 
      threshold = threshold
    ),
    dailyflow_spring = calculate_std(
      x$stream_discharge_mld,
      x$date_formatted, 
      season = 9:11, 
      fun = median
    ) / lt_median,
    dailyflow_summer = calculate_std(
      x$stream_discharge_mld,
      x$date_formatted, 
      season = 12:15, 
      fun = median
    ) / lt_median,
    dailyflow_winter = calculate_std(
      x$stream_discharge_mld,
      x$date_formatted, 
      season = 6:8, 
      fun = median
    ) / lt_median,
    maxantecedent_wateryear = calculate_std(
      x$stream_discharge_mld, 
      x$date_formatted,
      season = 7:18,
      lag = 1, 
      fun = max
    ) / lt_median,
    cvflow_spawning = calculate_std(
      x$stream_discharge_mld,
      x$date_formatted,
      season = 10:12, 
      fun = \(x) sd(x) / mean(x)
    ),
    spawning_temperature = calculate_std(
      x$water_temperature_c, 
      x$date_formatted, 
      season = 10:12, 
      fun = median
    )
  )
  
}

# function to calculate number of days below a given quantile
calculate_qx <- function(
    value, date, threshold, quantile = 0.1, season = 7:18, output = "metric"
) {
  
  # and count days below that threshold per season
  out <- calculate(
    value = value,
    date = date,
    resolution = survey(season = season),
    fun = days_below,
    threshold = threshold
  )
  
  # return date if needed
  if (output == "date") {
    out <- out$date
  } else {
    out <- out$metric
  }
  
  # and return 
  out
  
}

# function to calculate median-standardised flow metrics
calculate_std <- function(value, date, fun = median, season = 7:18, lag = 0, ...) {
  
  # and count days below that threshold per season
  out <- calculate(
    value = value,
    date = date,
    resolution = survey(season = season, lag = lag),
    fun = fun,
    standardise = NULL
  )
  
  # and return 
  out$metric
  
}

# function to load flow data or download it from WMIS if required
fetch_flow <- function(start, end, recompile = FALSE) {
  
  # load flow data already there, download otherwise
  flow_exists <- grepl("flow-daily-compiled", dir("data"))
  if (any(flow_exists) & !recompile) {
    
    # list all and sort to newest
    flow_files <- dir("data")[flow_exists]
    flow_files <- sort(flow_files, decreasing = TRUE)[1]
    flow <- qread(paste0("data/", flow_files))
    
  } else {
    
    # check if the WMIS data have been downloaded already, load saved
    #   version if so
    flow_list_exists <- grepl("flow-list", dir("data"))
    if (any(flow_list_exists) & !recompile) {
      
      # list all and load newest      
      flow_list_file <- dir("data")[flow_list_exists]
      flow_list_file <- sort(flow_list_file, decreasing = TRUE)[1]
      flow <- qread(paste0("data/", flow_list_file))
      
    } else {

      # list target reaches and associated flow/temp gauges
      targets <- c(
        "broken_creek_r4" = 404210,
        "glenelg_river_r1" = 238210,
        "glenelg_river_r2" = 238211,
        "glenelg_river_r3" = 238206,
        "goulburn_river_r4" = 405200,
        "holland_creek_r1" = 404207,
        "hughes_creek_r1" = 405228,
        "kiewa_river_r1" = 402203,
        "kiewawest_river_r1" = 402223,
        "kingparrot_creek_r1" = 405231,
        "moorabool_river_r3" = 232204,
        "moorabool_river_r4" = 232202,
        "moorabool_river_r4_temp" = 232242,
        "murray_river_r6" = 409025,  # Mulwala to Tocumwal (Yarrawonga Gauge?)
        "ovens_river_r0" = 403250,  # above myrtleford
        "ovens_river_r4" = 403230,  # Myrtleford to trib above ?? (check VEWH docs)
        "ovens_river_r5" = 403241,
        "seven_creeks_r1" = 405307,
        "seven_creeks_r2" = 405237,
        "thomson_river_r2" = 225231,
        "thomson_river_r2_temp" = 225212 ## WANDOCKA, DOWNSTREAM OF TARGET REACHES
        # "thomson_river_r3" = NA,  ## COOPERS CK FLOW IN data/
      )
      
      # specify location for NSW gauges (Murray River)
      gauge_location <- rep("vic", length(targets))
      nsw_gauges <- c(
        "murray_river_r6"
      )
      gauge_location[names(targets) %in% nsw_gauges] <- "nsw"
      
      # set the variable codes for water temp data (not the same for
      #   Vic vs the other two statesa)
      water_temp_code <- c("vic" = "450.00", "nsw" = "2080.00")
      
      # download data from WMIS
      flow <- vector("list", length = length(targets))
      names(flow) <- names(targets)
      for (i in seq_along(flow)) {
        
        # set variable codes by site
        varfrom <- c("100.00", "450.00")
        varto <- c("141.00", "450.00")
        
        # update for specific sites
        if (targets[i] == 404210)
          varfrom[1] <- "141.00"
        
        # grab data for each variable
        flow[[i]] <- fetch_hydro(
          sites = targets[i],
          state = gauge_location[i],
          start = dmy(paste0("01-01-", start)),
          end = dmy(paste0("31-12-", end)),
          options = list(
            varfrom = c(varfrom[1], water_temp_code[gauge_location[i]]),
            varto = c(varto[1], water_temp_code[gauge_location[i]])
          ),
          include_missing = TRUE
        )
        
      }
      
      # widen to put each variable in its own column
      flow <- lapply(
        flow,
        \(x) x |> pivot_wider(
          id_cols = c("date_formatted"),
          names_from = variable_name,
          values_from = c(value)
        )
      )
      
      # save to file
      qsave(flow, file = paste0("data/flow-list-", Sys.Date(), ".qs"))
    }
    
    # repeat for telemetry data source
    flow_telem_list_exists <- grepl("flow-telem-list", dir("data"))
    if (any(flow_telem_list_exists) & !recompile) {
      
      # load most recent version of saved file if it exists
      flow_telem_list_file <- dir("data")[flow_telem_list_exists]
      flow_telem_list_file <- sort(flow_telem_list_file, decreasing = TRUE)[1]
      flow_telem <- qread(paste0("data/", flow_telem_list_file))
      
    } else {
      
      # grab extra temperature data
      telem_targets <- c(
        "hughes_creek_r1" = 405228, 
        "seven_creeks_r1" = 405307
      )
      flow_telem <- vector("list", length = length(telem_targets))
      for (i in seq_along(telem_targets)) {
        
        flow_telem[[i]] <- fetch_hydro(
          sites = telem_targets[i],
          data_source = "TELEM",
          start = dmy(paste0("01-01-", start)),
          end = dmy(paste0("31-12-", end)),
          options = list(
            varfrom = c("450.00"),
            varto = c("450.00")
          ),
          include_missing = TRUE
        )
        
      }
      
      # add some names to this list
      names(flow_telem) <- names(telem_targets)
      
      # widen to put each variable in its own column
      flow_telem <- lapply(
        flow_telem,
        \(x) x |> pivot_wider(
          id_cols = c("date_formatted"),
          names_from = variable_name,
          values_from = c(value)
        )
      )
      
      # save to file
      qsave(flow_telem, file = paste0("data/flow-telem-list-", Sys.Date(), ".qs"))
      
    }
    
    # rename Murray R6 fields because NSW use different terms
    flow$murray_river_r6 <- flow$murray_river_r6 |>
      rename(
        stream_discharge_mld = discharge_rate, 
        water_temperature_c = water_temperature
      )
    
    # merge temperature data from backup gauges into main file
    flow$moorabool_river_r4$water_temperature_c <-
      flow$moorabool_river_r4_temp$water_temperature_c
    flow$moorabool_river_r4$`NA` <- NULL
    flow$thomson_river_r2$water_temperature_c <-
      flow$thomson_river_r2_temp$water_temperature_c
    flow$thomson_river_r2$`NA` <- NULL
    
    # and remove the _temp information that is now not needed
    flow <- flow[!grepl("_temp", names(flow))]

    # estimate flows in east branch of Kiewa by subtracting 
    #    west branch values from the DS (combined) flows
    flow$kiewaeast_river_r1 <- flow$kiewa_river_r1
    flow$kiewaeast_river_r1$stream_discharge_mld <-
      flow$kiewaeast_river_r1$stream_discharge_mld -
      flow$kiewawest_river_r1$stream_discharge_mld
    
    # and fill water temps with values from the west branch
    flow$kiewa_river_r1$water_temperature_c <- 
      flow$kiewawest_river_r1$water_temperature_c
    flow$kiewa_river_r1$`NA` <- NULL
    flow$kiewaeast_river_r1$water_temperature_c <- 
      flow$kiewawest_river_r1$water_temperature_c
    flow$kiewaeast_river_r1$`NA` <- NULL
    flow$kiewawest_river_r1 <- NULL
    
    # fill Thomson reach 3 from Coopers Ck flows, with other info
    #    copied from reach 2 gauge
    thomson_tmp <- read_xlsx("data/flow-raw/CoopersCk daily data_20050701-20221130.xlsx")
    thomson_tmp <- thomson_tmp |>
      mutate(
        date_formatted = parse_date_time(Date, orders = c("ymd_HMS", "ymd_HM", "ymd")),
        stream_discharge_mld = `ML/d`
      ) |> 
      select(date_formatted, stream_discharge_mld)
    flow$thomson_river_r3 <- flow$thomson_river_r2 |>
      left_join(thomson_tmp, by = "date_formatted", suffix = c("", "_manual")) |>
      mutate(stream_discharge_mld = stream_discharge_mld_manual) |>
      select(-stream_discharge_mld_manual)
    
    # and some more for the Ovens R0, Mullaroo, and Sevens R1
    #  (check whether these can be improved, particularly Sevens
    #   which might be filled with PMQ weir)
    flow$ovens_river_r0 <- flow$ovens_river_r0 |>
      left_join(
        flow$ovens_river_r4 |> select(date_formatted, stream_discharge_mld),
        by = "date_formatted",
        suffix = c("", "_backup")
      ) |>
      mutate(stream_discharge_mld = stream_discharge_mld_backup) |>
      select(-stream_discharge_mld_backup)
    flow$seven_creeks_r1 <- flow$seven_creeks_r1 |>
      left_join(
        flow$seven_creeks_r2 |> select(date_formatted, stream_discharge_mld),
        by = "date_formatted",
        suffix = c("", "_backup")
      ) |>
      mutate(stream_discharge_mld = stream_discharge_mld_backup) |>
      select(-stream_discharge_mld_backup)
    
    # fix up remaining gauges that have NA column instead of temperature
    flow$ovens_river_r4$water_temperature_c <-
      flow$ovens_river_r5$water_temperature_c
    flow$ovens_river_r4$`NA` <- NULL
    flow$hughes_creek_r1$water_temperature_c <-
      flow_telem$hughes_creek_r1$water_temperature_c
    flow$hughes_creek_r1$`NA` <- NULL
    flow$seven_creeks_r1$water_temperature_c <-
      flow_telem$seven_creeks_r1$water_temperature_c
    flow$seven_creeks_r1$`NA` <- NULL
    flow$seven_creeks_r2$water_temperature_c <-
      flow$seven_creeks_r1$water_temperature_c
    flow$seven_creeks_r2$`NA` <- NULL
    flow$kingparrot_creek_r1$water_temperature_c <-
      flow$hughes_creek_r1$water_temperature_c
    flow$kingparrot_creek_r1$`NA` <- NULL
    flow$holland_creek_r1$water_temperature_c <-
      flow$hughes_creek_r1$water_temperature_c
    
    # remove any NA dates
    flow <- lapply(flow, \(x) x[!is.na(x$date_formatted), ])
    
    # and fill gaps and cut off extreme values
    for (i in seq_along(flow)) {
      flow[[i]]$stream_discharge_mld <- impute_year(
        flow[[i]]$stream_discharge_mld,
        date = flow[[i]]$date_formatted,
        threshold = 100
      )
      flow[[i]]$stream_discharge_mld <- impute_rolling(
        flow[[i]]$stream_discharge_mld, 
        recursive = TRUE
      )
      flow[[i]]$water_temperature_c <- ifelse(
        flow[[i]]$water_temperature_c == 0, 
        NA, 
        flow[[i]]$water_temperature_c
      )
      flow[[i]]$water_temperature_c <- impute_year(
        flow[[i]]$water_temperature_c,
        date = flow[[i]]$date_formatted,
        threshold = 100
      )
      flow[[i]]$water_temperature_c <- impute_rolling(
        flow[[i]]$water_temperature_c, 
        n = 20,
        recursive = TRUE
      )
      flow[[i]]$water_temperature_c <- rm_extremes(
        flow[[i]]$water_temperature_c,
        range = c(3, 32)
      )
    }
    
    # save to file
    qsave(flow, file = paste0("data/flow-daily-compiled-", Sys.Date(), ".qs"))
    
  }
  
  # and return
  flow
  
}

# function to remove extreme values
rm_extremes <- function(x, range) {
  x <- ifelse(x < range[1], range[1], x)
  x <- ifelse(x > range[2], range[2], x)
  x
}
