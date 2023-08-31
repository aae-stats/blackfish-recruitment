
# extract YCS/recruitment estimates (work out which this is)
plot_flow_effects <- function(x, dat, average = FALSE, xname = NULL, rescale = NULL) {
  
  # reformat and tidy up names
  tmp <- x |>
    mutate(
      Waterbody = gsub("Gadopsis marmoratus|Gadopsis bispinosus|_sp", "", facet),
      Species = sapply(strsplit(as.character(facet) ,"_sp"), \(x) x[2])
    ) |>
    filter(Species == group)
  
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
  if (average) {
    tmp <- tmp |>
      group_by(x, Species) |>
      summarise(predicted = mean(predicted))
  } else {
    tmp <- tmp |>
      mutate(label = ifelse(x == max(x), Waterbody, NA))
  }
  
  # set up base plot, with wtaerbodies if needed
  if (average) {
    p <- tmp |>
      ggplot(aes(y = predicted, x = x))
  } else {
    p <-   tmp |>
      ggplot(aes(y = predicted, x = x, col = Waterbody))
  }
  
  # group by species
  p <- p +
    geom_line(linewidth = 1.25) +
    facet_wrap( ~ Species, scales = "free") +
    ylab("Young-of-year CPUE") +
    ylim(0, NA)
  
  # add xlab is known
  if (!is.null(xname)) 
    p <- p + xlab(xname)
  
  # add legend if required
  if (!average) {
    p <- p + 
      geom_label_repel(
        aes(label = label),
        nudge_x = 1,
        na.rm = TRUE
      ) +
      theme(legend.position = "none")
  }
  
  # return
  p
  
}
