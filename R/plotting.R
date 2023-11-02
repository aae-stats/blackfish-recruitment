# extract YCS/recruitment estimates (work out which this is)
plot_flow_effects <- function(x, dat, xname = NULL, rescale = NULL) {

  # reformat and tidy up names
  tmp <- x |>
    rename(
      Waterbody = facet,
      Species = group
    ) |>
    select(-conf.low, -conf.high)
  
  # rescale predictor if possible
  if (!is.null(rescale)) {
    tmp <- tmp |>
      mutate(x = rescale[1] + rescale[2] * x)
  }
  
  # set up a filter based on predicted counts being largely greater than zero
  include <- dat |> 
    group_by(species, waterbody) |> 
    summarise(catch = sum(catch)) |>
    mutate(facet = paste(waterbody, species, sep = "_sp"))
  tmp <- tmp |>
    left_join(include, by = c("Species" = "species", "Waterbody" = "waterbody")) |>
    filter(catch > 0)
  
  # average over waterbodies if required
  ave <- tmp |>
    group_by(x, Species) |>
    summarise(predicted = mean(predicted)) |>
    mutate(Waterbody = "Averaged")
  tmp <- bind_rows(tmp, ave)
  tmp <- tmp |>
    mutate(
      label = ifelse(x == max(x), Waterbody, NA),
      label = ifelse(label == "Averaged", "", label)
    )
  
  # set up a line width based on whether it's an average or waterbody value
  tmp <- tmp |>
    mutate(line_width = ifelse(Waterbody == "Averaged", 1.05, 1.03))
  
  # set up base plot, with wtaerbodies if needed
  p <-   tmp |>
    ggplot(aes(y = predicted, x = x, col = Waterbody))
  
  # group by species
  p <- p +
    geom_line(aes(linewidth = line_width)) +
    facet_wrap( ~ Species, scales = "free") +
    ylab("Young-of-year CPUE") +
    ylim(0, NA)
  
  # add xlab is known
  if (!is.null(xname)) 
    p <- p + xlab(xname)
  
  # add legend if required
  p <- p + 
    geom_label_repel(
      aes(label = label),
      nudge_x = 1,
      na.rm = TRUE
    ) +
    theme(legend.position = "none")
  
  # return
  p
  
}

# function to plot system-level associations for each variable
plot_flow_effects_system <- function(x, dat, metric_range, lt_vals, xname = NULL, rescale = NULL) {
  
  tmp <- x
  
  # rescale predictor if possible
  if (!is.null(rescale)) {
    tmp <- tmp |>
      mutate(x = rescale[1] + rescale[2] * x)
  }
  
  # restrict to range observed in that system
  tmp <- tmp |>
    left_join(metric_range, by = c("Waterbody", "Reach")) |>
    filter(x >= min_value, x <= max_value)
  
  # and rescale appropriately for the metric
  tmp <- tmp |>
    left_join(lt_vals, by = c("Waterbody", "Reach")) |>
    mutate(x = x * scale_factor)
  
  # set up a filter based on predicted counts being largely greater than zero
  include <- dat |> 
    group_by(species, waterbody, reach_no) |> 
    summarise(catch = sum(catch))
  tmp <- tmp |>
    left_join(include, by = c("Species" = "species", "Waterbody" = "waterbody", "Reach" = "reach_no")) |>
    filter(catch > 0)
  
  # plot each waterbody in turn
  waterbody_list <- tmp |>
    distinct(Waterbody, Reach)
  p <- vector("list", length = length(waterbody_list))
  for (i in seq_len(nrow(waterbody_list))) {
    
    # set up base plot, with wtaerbodies if needed
    p[[i]] <- tmp |>
      filter(Waterbody == waterbody_list$Waterbody[i], Reach == waterbody_list$Reach[i]) |>
      ggplot(aes(y = predicted, x = x))
    
    # group by species
    p[[i]] <- p[[i]] +
      geom_line() +
      facet_wrap( ~ Species, scales = "free") +
      ylab("Young-of-year CPUE") +
      ylim(0, NA)
    
    # add xlab is known
    if (!is.null(xname)) 
      p[[i]] <- p[[i]] + xlab(xname)
    
  }
  
  # add some names so we know which is which
  names(p) <- tolower(gsub(" ", "_", paste(waterbody_list$Waterbody, waterbody_list$Reach, sep = "_r")))
  
  # return
  p
  
}
