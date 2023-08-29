# internal function to download fish data from AAEDB
fetch_fish <- function(recompile = FALSE) {
  
  # list target species and waterbodies
  .waterbody_list <- c(
    "Broken Creek",
    "Broken River", 
    "Campaspe River",
    "Gellibrand River",
    "Glenelg River", 
    "Goulburn River",
    "Gunbower Creek",
    "Holland Creek",
    "Hughes Creek",
    "Kiewa River", 
    "Kiewa River East Branch",
    "King Parrot Creek", 
    "Latrobe River",
    "Lindsay River",
    "Little Murray River",
    "Loddon River", 
    "Love Creek",
    "Macalister River",
    "Moorabool River", 
    "Mullaroo Creek", 
    "Murray River",
    "Ovens River",
    "Pyramid Creek",
    "Seven Creeks",
    "Thomson River",
    "Wimmera River", 
    "Yarra River"
  )
  .species_list <- c(
    "Anguilla australis", 
    "Anguilla reinhardtii", 
    "Anguilla spp.", 
    "Bidyanus bidyanus", 
    "Carassius auratus",
    "Craterocephalus stercusmuscarum", 
    "Craterocephalus stercusmuscarum fulvus", 
    "Cyprinus carpio", 
    "fam. Eleotridae gen. Philypnodon",
    "Gadopsis bispinosus",
    "Gadopsis marmoratus",
    "Galaxias maculatus", 
    "Galaxias truttaceus", 
    "Maccullochella macquariensis",
    "Maccullochella peelii",
    "Macquaria ambigua", 
    "Macquaria australasica",
    "Melanotaenia fluviatilis",
    "Nannoperca australis", 
    "Nannoperca variegata", 
    "Nematalosa erebi", 
    "Perca fluviatilis", 
    "Percalates colonorum",
    "Percalates novemaculeatus",
    "Philypnodon grandiceps", 
    "Philypnodon macrostomus", 
    "Prototroctes maraena", 
    "Pseudaphritis urvillii",
    "Retropinna semoni", 
    "Salmo trutta", 
    "Tandanus tandanus", 
    "unidentified eel"
  )
  
  # check if data exist
  fish_exists <- any(grepl("fish-compiled.qs", dir("data/")))
  
  # if data exist and !recompile, load saved version. Re-extract otherwise
  if (fish_exists & !recompile) {
    
    # load data
    cpue <- qread("data/fish-compiled.qs")
    
  } else {
    
    # grab cpue data from AAEDB, filtered to targets
    cpue <- fetch_cpue(c(1, 2, 4, 6, 9, 10:13, 15:50)) |>
      filter(
        scientific_name %in% !!.species_list,
        waterbody %in% !!.waterbody_list
      ) |>
      collect()
    
    # group some species together
    cpue <- cpue |>
      mutate(
        species = scientific_name,
        species = ifelse(grepl("Philypnodon", species), "Philypnodon spp.", species),
        species = ifelse(grepl("unidentified eel|Anguilla", species), "Anguilla spp.", species),
        species = ifelse(grepl("Craterocephalus", species), "Craterocephalus spp.", species),
        species = ifelse(grepl("Carassius|Cyprinus", species), "Carp Goldfish", species)
      )
    
    # need to sum up multiple records of species groups at a site
    cpue <- cpue |>
      group_by(id_project, id_survey, species) |>
      summarise(catch = sum(catch)) |>
      ungroup() |>
      left_join(
        cpue |> distinct(id_project, id_survey, species, waterbody, id_site, survey_year, gear_type, effort_h),
        by = c("id_project", "id_survey", "species")
      )
    
    # add some site info
    site_info <- cpue |> fetch_site_info() |> collect()
    st_geometry(site_info) <- st_as_sfc(site_info$geom_pnt, crs = 4283)
    
    # ignored for now, can use VEWH reach table to add reach info    
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
    
    # save this
    qsave(cpue, file = "data/fish-compiled.qs")
    
  }
  
  # return
  cpue
  
}

# function to collate (coarse) fyke net data for a specific 
#   project
fetch_coarse_fykes <- function(recompile = FALSE, ...) {
  
  # list target species and waterbodies
  .waterbody_list <- c(
    "Gellibrand River",
    "King Parrot Creek", 
    "Love Creek"
  )
  .species_list <- c(
    "Anguilla australis", 
    "Anguilla reinhardtii", 
    "Anguilla spp.", 
    "Bidyanus bidyanus", 
    "Carassius auratus",
    "Craterocephalus stercusmuscarum", 
    "Craterocephalus stercusmuscarum fulvus", 
    "Cyprinus carpio", 
    "fam. Eleotridae gen. Philypnodon",
    "Gadopsis bispinosus",
    "Gadopsis marmoratus",
    "Galaxias maculatus", 
    "Galaxias truttaceus", 
    "Maccullochella macquariensis",
    "Maccullochella peelii",
    "Macquaria ambigua", 
    "Macquaria australasica",
    "Melanotaenia fluviatilis",
    "Nannoperca australis", 
    "Nannoperca variegata", 
    "Nematalosa erebi", 
    "Perca fluviatilis", 
    "Percalates colonorum",
    "Percalates novemaculeatus",
    "Philypnodon grandiceps", 
    "Philypnodon macrostomus", 
    "Prototroctes maraena", 
    "Pseudaphritis urvillii",
    "Retropinna semoni", 
    "Salmo trutta", 
    "Tandanus tandanus", 
    "unidentified eel"
  )
  # check if data exist
  fyke_exists <- any(grepl("fykes-compiled.qs", dir("data/")))
  
  # if data exist and !recompile, load saved version. Re-extract otherwise
  if (fyke_exists & !recompile) {
    
    # load data
    cpue <- qread("data/fykes-compiled.qs")
    
  } else {
    
    # download all catch data for the project (electro and netting
    #   to get all species)
    catch <- fetch_project(2:20) |> 
      collect()
    
    # calculate total catch by species and survey
    catch <- catch |> 
      filter(condition == "FISHABLE") |>
      mutate(
        collected = ifelse(is.na(collected), 0, collected),
        observed = ifelse(is.na(observed), 0, observed),
        catch = collected + observed
      ) |>
      group_by(
        id_project, waterbody, id_site, site_name, site_desc, id_survey,
        survey_date, survey_year, gear_type, scientific_name
      ) |>
      summarise(catch = sum(catch, na.rm = TRUE)) |>
      ungroup()
    
    # download the survey table for the project (all survey visits)
    survey_table <- fetch_netting_table(2:20) |>
      filter(gear_type == "FYKE") |>
      collect() |>
      left_join(catch |> distinct(id_survey, scientific_name), by = "id_survey") |>
      complete(nesting(id_project, waterbody, id_site, id_survey, gear_type, gear_count), scientific_name) |>
      filter(!is.na(scientific_name))
    
    # reduce catch to coarse fykes only
    catch <- catch |> filter(gear_type == "FYKE")
    
    # and use the survey table to pull out catch for every survey and species
    #   and calculate CPUE
    cpue <- survey_table |>
      left_join(
        catch |> 
          distinct(
            id_project, waterbody, id_site, site_name,
            site_desc, id_survey, survey_date, survey_year, 
            gear_type
          ), 
        by = c("id_project", "waterbody", "id_site", "id_survey", "gear_type")
      ) |>
      rename(effort_s = gear_count) |>
      left_join(
        catch |> select(id_survey, scientific_name, catch), 
        by = c("id_survey", "scientific_name")
      ) |>
      mutate(
        catch = ifelse(is.na(catch), 0, catch),
        effort_h = effort_s,
        cpue = catch / effort_h
      )
    
    # add timezone and pull out target columns
    cpue <- cpue |>
      mutate(extracted_ts = sql("timezone('Australia/Melbourne'::text, now())")) |>
      select(all_of(aae.db:::cpue_return_cols)) |>
      filter(
        !is.na(scientific_name),
        !is.na(effort_h)
      )
    
    # collect
    cpue <- cpue |> collect()
    
    # filter to target species and waterbodies
    cpue <- cpue |> 
      filter(
        scientific_name %in% !!.species_list,
        waterbody %in% !!.waterbody_list
      )
    
    # group some species together
    cpue <- cpue |>
      mutate(
        species = scientific_name,
        species = ifelse(grepl("Philypnodon", species), "Philypnodon spp.", species),
        species = ifelse(grepl("unidentified eel|Anguilla", species), "Anguilla spp.", species),
        species = ifelse(grepl("Craterocephalus", species), "Craterocephalus spp.", species),
        species = ifelse(grepl("Carassius|Cyprinus", species), "Carp Goldfish", species)
      )
    
    # need to sum up multiple records of species groups at a site
    cpue <- cpue |>
      group_by(id_project, id_survey, species) |>
      summarise(catch = sum(catch)) |>
      ungroup() |>
      left_join(
        cpue |> distinct(id_project, id_survey, species, waterbody, id_site, survey_year, gear_type, effort_h),
        by = c("id_project", "id_survey", "species")
      )
    
    # add site info
    site_info <- cpue |> fetch_site_info() |> collect()
    st_geometry(site_info) <- st_as_sfc(site_info$geom_pnt, crs = 4283)
    cpue <- cpue |>
      left_join(
        site_info |>
          select(id_project, id_site, reach_no, latitude, longitude),
        by = c("id_project", "id_site")
      )
    
    
    # assign reach info
    cpue <- cpue |> 
      mutate(
        reach_no = ifelse(is.na(reach_no), 1, reach_no),
        reach_no = as.character(reach_no)
      )
    
    # filter out sites without geom information
    cpue <- cpue |>
      filter(
        id_site %in% !!(site_info |> filter(!is.na(geom_pnt)) |> pull(id_site) |> unique())
      )
    
    # and cast id_site to int
    cpue <- cpue |> mutate(id_site = as.integer(id_site))
    
    # save this
    qsave(cpue, file = "data/fykes-compiled.qs")
    
  }
  
  # return
  cpue
  
}

# function to collate info on netting effort by survey
fetch_netting_table <- function (project_id, ...) {
  
  # grab full survey table
  survey_table <- fetch_table("site") |>
    select(id_site, waterbody) |>
    left_join(
      fetch_table("survey") |> select(id_site, id_survey, id_project, gear_type, released),
      by = "id_site"
    ) |>
    left_join(
      fetch_table("survey_event") |>
        select(id_survey, id_surveyevent, condition),
      by = "id_survey"
    )
  
  # add netting info for each survey event
  survey_table <- survey_table |> 
    left_join(
      fetch_table("netting") |> select(id_surveyevent, id_netting, gear_count),
      by = "id_surveyevent"
    )
  
  # return just coarse fykes at fishable locations
  survey_table <- survey_table |>
    filter(
      id_project %in% project_id, 
      gear_type == "FYKE",
      released,
      condition == !!"FISHABLE"
    ) |>
    select(-released, -condition)
  
  # want total number of nets set per survey
  survey_table <- survey_table %>%
    distinct(id_project, id_site, waterbody, id_survey, id_surveyevent, gear_type, gear_count) |>
    group_by(id_project, id_site, waterbody, id_survey, gear_type) |>
    summarise(gear_count = sum(gear_count, na.rm = TRUE)) |>
    ungroup()
  
  # tidy up classes
  survey_table <- survey_table |>
    mutate(
      id_project = as.integer(id_project), 
      id_site = as.integer(id_site),
      id_survey = as.integer(id_survey), 
      gear_count = as.integer(gear_count)
    )
  
  # return
  survey_table
  
}
