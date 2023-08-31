# internal function to download fish data from AAEDB
fetch_fish <- function(recompile = FALSE) {
  
  # list target ARI projects
  .project_list <- c(1, 2, 4, 6, 9, 10:13, 15:50)
  
  # list target species and waterbodies
  .waterbody_list <- c(
    "Broken Creek",
    "Glenelg River", 
    "Goulburn River",
    "Holland Creek",
    "Hughes Creek",
    "Kiewa River", 
    "Kiewa River East Branch",
    "King Parrot Creek", 
    "Moorabool River",
    "Murray River",
    "Ovens River",
    "Seven Creeks",
    "Thomson River"
  )
  .species_list <- c(
    "Gadopsis bispinosus",
    "Gadopsis marmoratus"
  )
  
  # list of reaches to exclude
  .exclusions <- c(
    "Broken Creek R5",
    "Goulburn River R5",
    "Moorabool River R1",
    "Murray River R8",
    "Thomson River R4",
    "Thomson River R5",
    "Thomson River R6"
  )
  
  # check if data exist
  fish_exists <- any(grepl("fish-compiled.qs", dir("data/")))
  
  # if data exist and !recompile, load saved version. Re-extract otherwise
  if (fish_exists & !recompile) {
    
    # load data
    cpue <- qread("data/fish-compiled.qs")
    
  } else {

    # grab CPUE values for 0+ and 1+ fish based on length intervals
    bf_0plus <- fetch_cpue(.project_list, criterion = list(var = "length_cm", lower = 0, upper = 8)) |>
      filter(
        scientific_name %in% !!.species_list,
        waterbody %in% !!.waterbody_list
      ) |>
      collect() |>
      mutate(estimated_age = 0)
    bf_1plus <- fetch_cpue(.project_list, criterion = list(var = "length_cm", lower = 8, upper = 16.5)) |>
      filter(
        scientific_name %in% !!.species_list,
        waterbody %in% !!.waterbody_list
      ) |>
      collect() |>
      mutate(estimated_age = 1)
    
    # we want to work with a single CPUE data set
    cpue <- bind_rows(bf_0plus, bf_1plus)
        
    # add some site info
    site_info <- cpue |> fetch_site_info() |> collect()
    st_geometry(site_info) <- st_as_sfc(site_info$geom_pnt, crs = 4283)
    
    # can use VEWH reach table to add reach info    
    vewh_reaches <- fetch_table("eflow_reaches_20171214", "spatial") |>
      collect()
    st_geometry(vewh_reaches) <- st_as_sfc(vewh_reaches$geom, crs = 4283)
    site_info <- site_info |>
      st_join(vewh_reaches |> select(vewh_reach), join = st_within) |>
      mutate(reach_no = ifelse(is.na(reach_no), vewh_reach, reach_no)) |>
      select(-vewh_reach)
    
    # grab a few Ovens sites and duplicate for id_project = 9
    ovens_sub <- site_info |>
      filter(id_site %in% c(3160, 3162, 3163, 3183, 3185)) |>
      mutate(id_project = 9)
    site_info <- bind_rows(site_info, ovens_sub)
    
    # add reach info for unknown reaches and then append to CPUE data
    site_info <- site_info |>
      mutate(
        id_site = as.numeric(id_site),
        reach_no = ifelse(waterbody == "Boggy Creek", 1, reach_no),
        reach_no = ifelse(waterbody == "Yea River", 1, reach_no),
        reach_no = ifelse(waterbody == "Gellibrand River", 1, reach_no),
        reach_no = ifelse(waterbody == "Love Creek", 1, reach_no),
        reach_no = ifelse(waterbody == "Murray River", 6, reach_no),
        reach_no = ifelse(waterbody == "Gunbower Creek", 1, reach_no),
        reach_no = ifelse(waterbody == "Holland Creek", 1, reach_no),
        reach_no = ifelse(waterbody == "Hughes Creek", 1, reach_no),
        reach_no = ifelse(waterbody == "Kiewa River", 1, reach_no),
        reach_no = ifelse(waterbody == "Kiewa River East Branch", 1, reach_no),
        reach_no = ifelse(waterbody == "King Parrot Creek", 1, reach_no),
        reach_no = ifelse(waterbody == "Lindsay River", 1, reach_no),
        reach_no = ifelse(waterbody == "Mullaroo Creek", 1, reach_no),
        reach_no = ifelse(id_site %in% c(2058, 4464), 3, reach_no),
        reach_no = ifelse(id_site %in% c(4466, 4467, 4470, 4471), 4, reach_no),
        reach_no = ifelse(id_site %in% c(4468, 4066, 4068, 4069, 4073), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4067, 4070:4072), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3193), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3194, 4266), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4060, 4061), 1, reach_no),
        reach_no = ifelse(id_site %in% c(4382, 4383), 3, reach_no),
        reach_no = ifelse(id_site %in% c(4109, 4110, 4115), 4, reach_no),
        reach_no = ifelse(id_site %in% c(4116), 5, reach_no),
        reach_no = ifelse(id_site %in% c(3133), 0, reach_no),
        reach_no = ifelse(id_site %in% c(3134), 5, reach_no),
        reach_no = ifelse(id_site %in% c(1643:1646, 3164, 4225, 4229, 4231), 0, reach_no),
        reach_no = ifelse(id_site %in% c(3160:3163, 3182, 4172:4185, 4194:4197, 4199:4204, 4208:4212, 4217:4224, 4228, 4232:4241), 4, reach_no),
        reach_no = ifelse(waterbody == "Seven Creeks", 1, reach_no),
        reach_no = ifelse(id_site %in% c(3322, 3324, 3325), 2, reach_no),
        reach_no = ifelse(id_site %in% c(3167), 3, reach_no),
        reach_no = ifelse(id_site %in% c(3112, 3177:3180, 4451), 5, reach_no),
        reach_no = ifelse(id_site %in% c(3113, 3181), 4, reach_no),
        reach_no = ifelse(id_site %in% c(3171:3176), 6, reach_no),
        reach_no = ifelse(id_site %in% c(4275:4280), 8, reach_no),  # Murray Tocumwal to Gunbower ## kEEP THESE
        reach_no = ifelse(id_site %in% c(4260:4265), 9, reach_no),  # Murray below Gunbower  # DROP THESE
        reach_no = ifelse(waterbody == "Broken Creek" & is.na(reach_no), 5, reach_no),
        reach_no = ifelse(waterbody == "Ovens River" & is.na(reach_no), 5, reach_no)
      )
    
    # add reach and lat/long info, filter to remove Buffalo sites and add reach info
    #   for Ovens sites in id_project 9
    cpue <- cpue |>
      left_join(
        site_info |>
          select(id_project, id_site, reach_no, latitude, longitude),
        by = c("id_project", "id_site")
      ) |>
      filter(!(id_site %in% c(4226, 4227, 4230))) |>
      mutate(waterbody = ifelse(id_site %in% c(3134), "Thomson River", waterbody))
    
    # generate a list of Murray River sites and filter to "Lower" and 
    #   id_site %in% c(4275:4280)
    murray_site_list <- cpue |>
      filter(waterbody == "Murray River") |>
      pull(id_site) |>
      unique()
    murray_sites <- fetch_table("site") |>
      filter(id_site %in% !!murray_site_list) |>
      collect()
    murray_to_keep <- murray_sites |> 
      filter(site_desc == "Lower") |>
      pull(id_site) |>
      unique()
    others_to_keep <- cpue |>
      filter(waterbody != "Murray River") |>
      pull(id_site) |>
      unique()
    cpue <- cpue |>
      filter(id_site %in% !!(c(murray_to_keep, others_to_keep, 4275:4280)))
    
    # filter out some upper reaches of the Moorabool and Macalister (0, 2, 0)
    cpue <- cpue |>
      filter(
        !(waterbody == "Macalister River" & reach_no == 0),
        !(waterbody == "Moorabool River" & reach_no %in% c(0, 2))
      )
    
    # filter out sites without geom information
    cpue <- cpue |>
      filter(
        id_site %in% !!(site_info |> filter(!is.na(geom_pnt)) |> pull(id_site) |> unique())
      )
    
    # and remove a final list of excluded reaches based on low catch
    cpue <- cpue |>
      mutate(waterbody_reach = paste(waterbody, reach_no, sep = " R")) |>
      filter(!waterbody_reach %in% .exclusions) |>
      select(-waterbody_reach)

    # save this
    qsave(cpue, file = "data/fish-compiled.qs")
    
  }
  
  # return
  cpue
  
}
