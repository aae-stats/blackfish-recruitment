# analysis of blackfish recruitment based on age-length estimates

# Author: Jian Yen (jian.yen [at] deeca.vic.gov.au)
# Created: August 2023
# Last updated: 20 March 2024

# load packages
library(qs)
library(readxl)
library(aae.db)
library(aae.hydro)
library(lubridate)
library(dplyr)
library(tidyr)
library(sf)
library(brms)
library(ggplot2)
library(bayesplot)
library(ragg)
library(ggeffects)
library(ggrepel)
library(patchwork)

# load helpers
source("R/fish.R")
source("R/flow.R")
source("R/plotting.R")

# flags to refit model
refit_model <- FALSE
rerun_diagnostics <- FALSE

# check correlations among different levels
cpue_corr <- check_correlations(recompile = FALSE)
threshold_list <- expand.grid(
  yoy = c(7, 8, 9),
  oneplus = c(15, 16, 18, 20)
)
cpue_corr <- mapply(
  \(x, y) x |> mutate(yoy_threshold = y[1, 1], oneplus_threshold = y[1, 2]),
  cpue_corr,
  lapply(seq_len(nrow(threshold_list)), \(i, .x) .x[i, ], .x = threshold_list),
  SIMPLIFY = FALSE
)
cpue_corr <- bind_rows(cpue_corr) |>
  select(waterbody, id_site, id_survey, gear_type, scientific_name, estimated_age, catch, yoy_threshold, oneplus_threshold) |>
  mutate(threshold = paste0("threshold_", paste(yoy_threshold, oneplus_threshold, sep = "_"))) |>
  select(-yoy_threshold, -oneplus_threshold) |>
  pivot_wider(
    id_cols = c(waterbody, id_site, id_survey, gear_type, scientific_name, estimated_age),
    values_from = catch,
    names_from = threshold
  )
groupings <- cpue_corr |> distinct(scientific_name, estimated_age)
corr <- vector("list", length = nrow(groupings))
for (i in seq_len(nrow(groupings))) {
  corr[[i]] <- cpue_corr |>
    filter(
      scientific_name == groupings$scientific_name[i], 
      estimated_age == groupings$estimated_age[i]
    ) |>
    select(contains("threshold")) |>
    cor()
  write.csv(
    corr[[i]],
    paste0(
      "outputs/tables/corr-", 
      paste(
        tolower(gsub(" ", "_", groupings$scientific_name[i])),
        groupings$estimated_age[i], sep = "-age"
      ),
      ".csv"
    )
  )
}


# need fish catches with age information (filtered to 0+ and 1+)
cpue <- fetch_fish(recompile = FALSE)

# separate adults and add recuit_year info to juves
cpue_ad <- cpue |> 
  filter(estimated_age > 2)
cpue <- cpue |> 
  filter(estimated_age <= 2)

# work out waterbody lengths
bf_spatial <- st_read("data/blackfish-spatial/blackfish_reaches_231222.shp")
bf_spatial <- bf_spatial |>
  mutate(segment_length = st_length(bf_spatial)) |>
  group_by(waterbody) |>
  summarise(surveyed_length_m = sum(segment_length)) |>
  ungroup() |>
  mutate(surveyed_length_km = surveyed_length_m / 1000)
bf_spatial |>
  as_tibble() |>
  select(waterbody, surveyed_length_m, surveyed_length_km) |>
  write.csv("outputs/tables/surveyed-waterbody-lengths.csv")

# want flow data for all target systems and reaches
flow <- fetch_flow(start = 1997, end = 2023, recompile = FALSE)
sapply(
  flow,
  \(x) x |>
    group_by(year(date_formatted)) |> 
    summarise(annual_discharge = sum(stream_discharge_mld)) |> 
    pull(annual_discharge) |> 
    mean(na.rm = TRUE)
) |>
  as_tibble() |>
  rename(mean_annual_discharge_ml = value) |>
  mutate(
    waterbody_reach = names(flow),
    mean_annual_discharge_gl = mean_annual_discharge_ml * 0.001
  ) |>
  write.csv("outputs/tables/average-annual-discharge.csv")

# and need metrics based on these flow data
metrics <- mapply(calculate_metrics, flow, names(flow), SIMPLIFY = FALSE)
metrics <- bind_rows(metrics)

# pull out daily flow for each survey as a detection covariate
flow <- mapply(
  \(x, y) x |> mutate(
    waterbody = waterbody_lu[[y]]$waterbody, 
    reach_no = as.character(waterbody_lu[[y]]$reach_no)
  ),
  x = flow,
  y = names(flow),
  SIMPLIFY = FALSE
)
flow <- bind_rows(flow)
cpue <- cpue |>
  left_join(
    flow |> select(date_formatted, waterbody, reach_no, stream_discharge_mld),
    by = c("survey_date" = "date_formatted", "waterbody", "reach_no")
  ) |>
  rename(survey_flow = stream_discharge_mld)

# summary info for reporting
fish_summary <- cpue |>
  group_by(waterbody, reach_no)|>
  left_join(
    cpue |> 
      group_by(waterbody, reach_no, survey_year) |> 
      summarise(no_sites = length(unique(id_site))),
    by = c("waterbody", "reach_no", "survey_year")
  ) |>
  summarise(
    average_no_sites_per_year = mean(no_sites),
    min_no_sites_per_year = min(no_sites),
    max_no_sites_per_year = max(no_sites),
    years_surveyed = paste(sort(unique(survey_year)), collapse = ", "),
    months_surveyed = paste(sort(unique(month(survey_date))), collapse = ", "),
    gear_used = paste(sort(unique(gear_type)), collapse = ", "),
    electro_seconds_range = paste(min(effort_s), max(effort_s), sep = "-"),
    electro_seconds_mean = mean(effort_s)
  )
fish_summary_simple <- cpue |>
  group_by(waterbody)|>
  left_join(
    cpue |> 
      group_by(waterbody, survey_year) |> 
      summarise(no_sites = length(unique(id_site))),
    by = c("waterbody", "survey_year")
  ) |>
  summarise(
    average_no_sites_per_year = mean(no_sites),
    min_no_sites_per_year = min(no_sites),
    max_no_sites_per_year = max(no_sites),
    years_surveyed = paste(sort(unique(survey_year)), collapse = ", "),
    months_surveyed = paste(sort(unique(month(survey_date))), collapse = ", "),
    gear_used = paste(sort(unique(gear_type)), collapse = ", "),
    electro_seconds_range = paste(min(effort_s), max(effort_s), sep = "-"),
    electro_seconds_mean = mean(effort_s)
  )

#   summarise flow
flow_summary <- flow |>
  mutate(water_year = ifelse(month(date_formatted) > 6, year(date_formatted) + 1, year(date_formatted))) |>
  group_by(waterbody, reach_no, water_year) |>
  summarise(
    average_discharge_mld = mean(stream_discharge_mld, na.rm = TRUE),
    total_discharge_mld = sum(stream_discharge_mld, na.rm = TRUE),
    average_water_temperature_c = mean(water_temperature_c, na.rm = TRUE),
    min_water_temperature_c = min(water_temperature_c, na.rm = TRUE),
    max_water_temperature_c = max(water_temperature_c, na.rm = TRUE)
  )
flow_summary_simple <- flow_summary |>
  group_by(waterbody, reach_no) |>
  summarise(
    mean_annual_discharge_mld = mean(total_discharge_mld),
    mean_water_temperature_c = mean(average_water_temperature_c),
    mean_min_water_temperature_c = mean(min_water_temperature_c),
    mean_max_water_temperature_c = mean(max_water_temperature_c)
  )

# write to file
write.csv(fish_summary, file = "outputs/tables/blackfish-survey-info-by-reach.csv")
write.csv(fish_summary_simple, file = "outputs/tables/blackfish-survey-info.csv")
write.csv(flow_summary, file = "outputs/tables/blackfish-flow-info-by-year.csv")
write.csv(flow_summary_simple, file = "outputs/tables/blackfish-flow-info.csv")

# flatten to a single observation per year
cpue <- cpue |>
  group_by(
    id_site,
    waterbody,
    reach_no,
    survey_year, 
    gear_type,
    scientific_name,
    estimated_age
  ) |>
  summarise(
    survey_flow = median(survey_flow),
    effort_h = sum(effort_h),
    catch = sum(catch)
  ) |>
  ungroup() |>
  mutate(survey_flow_std = scale(survey_flow))
cpue_ad <- cpue_ad |>
  group_by(
    id_site,
    survey_year, 
    gear_type,
    scientific_name
  ) |>
  summarise(catch = sum(catch)) |>
  ungroup()

# include adult_ym1 catch as a predictor (species-specific)
cpue <- cpue |> 
  mutate(
    recruit_year = survey_year - estimated_age,
    recruit_year_m1 = recruit_year - 1
  ) |>
  left_join(
    cpue_ad |> 
      select(id_site, survey_year, gear_type, scientific_name, catch) |>
      rename(adult_catch_m1 = catch),
    by = c("id_site", "recruit_year_m1" = "survey_year", "gear_type", "scientific_name")
  ) |>
  mutate(
    adult_catch_m1 = ifelse(is.na(adult_catch_m1), 0, adult_catch_m1),
    log_adult_catch_m1 = log(adult_catch_m1 + 1)
  )

# attach flow to cpue data, aligning 1+ fish with flow in the previous year
#    (the year in which they were spawned)
cpue <- cpue |>
  mutate(reach_no = as.numeric(reach_no)) |>
  left_join(
    metrics, 
    by = c("waterbody", "reach_no", "recruit_year" = "water_year")
  ) |>
  mutate(log_effort_h = log(effort_h))

# standardise flow metrics
cpue <- cpue |>
  mutate(
    species = scientific_name,
    lowflow_wateryear_std = scale(lowflow_wateryear),
    dailyflow_spring_std = scale(dailyflow_spring),
    dailyflow_summer_std = scale(dailyflow_summer),
    dailyflow_winter_std = scale(dailyflow_winter),
    maxantecedent_wateryear_std = scale(maxantecedent_wateryear),
    cvflow_spawning_std = scale(cvflow_spawning),
    spawning_temperature_std = scale(spawning_temperature)
  )

# save a copy of raw data for mapping
write.csv(cpue |> distinct(id_site, waterbody, reach_no), file = "outputs/tables/bf-cpue.csv")

# plot flow conditions
var_lookup <- c(
  "cvflow_spawning" = "Flow variability",
  "dailyflow_spring" = "Spring flow",
  "dailyflow_summer" = "Summer flow",
  "dailyflow_winter" = "Winter flow",
  "lowflow_wateryear" = "Days of low flow",
  "maxantecedent_wateryear" = "Maximum antecedent flow",
  "spawning_temperature" = "Water temp. during spawning"
)
metrics_to_plot <- metrics |>
  select(
    waterbody, reach_no, water_year,
    lowflow_wateryear, dailyflow_spring,
    dailyflow_summer, dailyflow_winter,
    maxantecedent_wateryear, cvflow_spawning,
    spawning_temperature
  ) |>
  pivot_longer(
    cols = c(contains("flow"), contains("temperature"), contains("wateryear")),
    names_to = "metric"
  ) |>
  mutate(
    waterbody_reach = paste(waterbody, reach_no, sep = ": Reach "),
    metric = var_lookup[metric]
  )

p_coastal <- metrics_to_plot |>
  filter(grepl("Glenelg|Moorabool|Thomson", waterbody)) |>
  ggplot(aes(y = value, x = water_year, col = waterbody_reach)) +
  geom_line() +
  facet_wrap( ~ metric, scales = "free") +
  scale_color_brewer(palette = "Set2", name = "") +
  xlab("") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 3))
p_inland <- metrics_to_plot |>
  filter(!grepl("Glenelg|Moorabool|Thomson", waterbody)) |>
  ggplot(aes(y = value, x = water_year, col = waterbody_reach)) +
  geom_line() +
  facet_wrap( ~ metric, scales = "free") +
  scale_color_discrete(name = "") +
  xlab("") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 3))
p_var_coastal <- metrics_to_plot |>
  filter(
    grepl("Glenelg|Moorabool|Thomson", waterbody),
    metric == "Flow variability"
  ) |>
  ggplot(aes(y = value, x = water_year)) +
  geom_line() +
  facet_wrap( ~ waterbody_reach) +
  xlab("") +
  theme(legend.position = "bottom")
p_var_inland <- metrics_to_plot |>
  filter(
    !grepl("Glenelg|Moorabool|Thomson", waterbody),
    metric == "Flow variability"
  ) |>
  ggplot(aes(y = value, x = water_year)) +
  geom_line() +
  facet_wrap( ~ waterbody_reach) +
  xlab("") +
  theme(legend.position = "bottom")

# save to file
ggsave(
  file = "outputs/figures/flow-metrics-inland.png",
  plot = p_inland,
  device = ragg::agg_png,
  width = 7,
  height = 7,
  units = "in",
  res = 600
)
ggsave(
  file = "outputs/figures/flow-metrics-coastal.png",
  plot = p_coastal,
  device = ragg::agg_png,
  width = 7,
  height = 7,
  units = "in",
  res = 600
)
ggsave(
  file = "outputs/figures/flow-variability-coastal.png",
  plot = p_var_coastal,
  device = ragg::agg_png,
  width = 7,
  height = 7,
  units = "in",
  res = 600
)
ggsave(
  file = "outputs/figures/flow-variability-inland.png",
  plot = p_var_inland,
  device = ragg::agg_png,
  width = 7,
  height = 7,
  units = "in",
  res = 600
)

# fit model or load saved model (set at top)
if (refit_model) {
  
  # fit a model
  stan_seed <- 2023-08-29
  mod <- brm(
    bf(
      catch ~ estimated_age +
        log_adult_catch_m1 +
        survey_flow_std + 
        lowflow_wateryear_std +
        dailyflow_spring_std + 
        I(dailyflow_spring_std ^ 2) +
        dailyflow_summer_std + 
        I(dailyflow_summer_std ^ 2) +
        dailyflow_winter_std +
        maxantecedent_wateryear_std +
        cvflow_spawning_std +
        spawning_temperature_std + 
        I(spawning_temperature_std ^ 2) +
        (estimated_age +
           lowflow_wateryear_std +
           dailyflow_spring_std + 
           I(dailyflow_spring_std ^ 2) +
           dailyflow_summer_std + 
           I(dailyflow_summer_std ^ 2) +
           dailyflow_winter_std +
           maxantecedent_wateryear_std +
           cvflow_spawning_std +
           spawning_temperature_std + 
           I(spawning_temperature_std ^ 2) | species) +
        (estimated_age +
           lowflow_wateryear_std +
           dailyflow_spring_std +
           dailyflow_summer_std +
           dailyflow_winter_std +
           maxantecedent_wateryear_std +
           cvflow_spawning_std +
           spawning_temperature_std | species:waterbody) +
        (1 | waterbody) +
        (1 | waterbody:reach_no) +
        (1 | species:waterbody:reach_no) +
        (1 | survey_year) +
        (1 | species:survey_year) +
        (1 | species:survey_year:waterbody) +
        (1 | species:survey_year:waterbody:reach_no) +
        (1 | recruit_year) +
        (1 | recruit_year:species) +
        (1 | id_site) +
        (1 | id_site:species) +
        (1 | gear_type) +
        (1 | gear_type:species) +
        offset(log_effort_h),
      shape ~ (1 | species) +
        (1 | waterbody) +
        (1 | species:waterbody) +
        (1 | waterbody:reach_no) +
        (1 | species:waterbody:reach_no)
    ),
    data = cpue,
    family = negbinomial(),
    chains = 4,
    cores = 4,
    seed = stan_seed,
    iter = 3000,
    warmup = 2000,
    control = list(adapt_delta = 0.9, max_treedepth = 15),
    backend = "rstan",
    refresh = 100,
    silent = 0,
    threads = threading(3)
  )
  
  # save fitted
  qsave(mod, file = "outputs/fitted/bf-mod-quad.qs")
  
} else {
  
  # load saved model  
  mod <- qread("outputs/fitted/bf-mod-quad.qs")
  
}

# run diagnostics
if (rerun_diagnostics) {
  pp_vals <- pp_check(mod, type = "hist")
  pp_max <- pp_check(mod, group = "species", type = "stat", stat = "max")
  pp_pzero_grouped <- pp_check(mod, group = "species", type = "stat_grouped", stat = \(x) mean(x == 0))
  rhat <- brms::rhat(mod)
  neff <- brms::neff_ratio(mod)
  standard_diag <- tibble(
    stat = c(rep("Rhat", length(rhat)), rep("Neff ratio", length(neff))),
    par = c(names(rhat), names(neff)),
    value = c(rhat, neff)
  ) |>
    ggplot(aes(x = value)) +
    geom_histogram() +
    xlab("Value") +
    ylab("Count") +
    facet_wrap( ~ stat, scales = "free")
  ggsave(
    filename = "outputs/figures/pp-hist.png",
    plot = pp_vals + scale_x_log10(),
    device = agg_png,
    width = 6,
    height = 6,
    dpi = 600,
    units = "in"
  )
  ggsave(
    filename = "outputs/figures/pp-max.png",
    plot = pp_max,
    device = agg_png,
    width = 6,
    height = 6,
    dpi = 600,
    units = "in"
  )
  ggsave(
    filename = "outputs/figures/pp-pzero.png",
    plot = pp_pzero_grouped,
    device = agg_png,
    width = 6,
    height = 6,
    dpi = 600,
    units = "in"
  )
  ggsave(
    filename = "outputs/figures/diagnostics.png",
    plot = standard_diag,
    device = agg_png,
    width = 6,
    height = 6,
    dpi = 600,
    units = "in"
  )
}

# plot fitted vs observed by year
fitted_plot <- posterior_epred(mod)
fitted_obs_plot <- tibble(
  catch = c(mod$data$catch, apply(fitted_plot, 2, median)),
  lower = c(
    rep(NA, nrow(mod$data)), 
    apply(fitted_plot, 2, quantile, probs = 0.1)
  ),
  upper = c(rep(NA, nrow(mod$data)),
            apply(fitted_plot, 2, quantile, probs = 0.9)
  ),
  estimated_age = rep(mod$data$estimated_age, 2), 
  species = rep(mod$data$species, 2),
  waterbody = rep(mod$data$waterbody, 2),
  reach_no = rep(mod$data$reach_no, 2),
  survey_year = rep(mod$data$survey_year, 2),
  id_site = rep(mod$data$id_site, 2)
)
include_subset <- fitted_obs_plot |>
  group_by(species, waterbody, reach_no) |>
  summarise(sum_catch = sum(catch))
fitted_0plus_ts <- fitted_obs_plot |>
  left_join(include_subset, by = c("species", "waterbody", "reach_no")) |>
  mutate(category = rep(c("Observed", "Modelled"), each = nrow(mod$data))) |>
  filter(estimated_age == 0, species == "Gadopsis bispinosus") |>
  filter(sum_catch > 1) |>
  group_by(waterbody, reach_no, category, survey_year) |>
  summarise(
    catch = sum(catch),
    lower = sum(lower),
    upper = sum(upper)
  ) |>
  ungroup() |>
  mutate(waterbody = paste(waterbody, reach_no, sep = ": Reach ")) |>
  ggplot(
    aes(y = catch, x = survey_year, ymin = lower, ymax = upper, col = category)
  ) +
  geom_point(position = position_dodge(0.1)) +
  geom_errorbar(position = position_dodge(0.1)) +
  scale_color_brewer(palette = "Set2", name = "") +
  xlab("Survey year") +
  ylab("Catch") +
  facet_wrap( ~ waterbody, scales = "free_y") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 8)
  )
fitted_1plus_ts <- fitted_obs_plot |>
  left_join(include_subset, by = c("species", "waterbody", "reach_no")) |>
  mutate(category = rep(c("Observed", "Modelled"), each = nrow(mod$data))) |>
  filter(estimated_age == 1, species == "Gadopsis bispinosus") |>
  filter(sum_catch > 1) |>
  group_by(waterbody, reach_no, category, survey_year) |>
  summarise(
    catch = sum(catch),
    lower = sum(lower),
    upper = sum(upper)
  ) |>
  ungroup() |>
  mutate(waterbody = paste(waterbody, reach_no, sep = ": Reach ")) |>
  ggplot(
    aes(y = catch, x = survey_year, ymin = lower, ymax = upper, col = category)
  ) +
  geom_point(position = position_dodge(0.1)) +
  geom_errorbar(position = position_dodge(0.1)) +
  scale_color_brewer(palette = "Set2", name = "") +
  xlab("Survey year") +
  ylab("Catch") +
  facet_wrap( ~ waterbody, scales = "free_y") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 8)
  )
fitted_0plus_rb <- fitted_obs_plot |>
  left_join(include_subset, by = c("species", "waterbody", "reach_no")) |>
  mutate(category = rep(c("Observed", "Modelled"), each = nrow(mod$data))) |>
  filter(estimated_age == 0, species == "Gadopsis marmoratus") |>
  filter(sum_catch > 1) |>
  group_by(waterbody, reach_no, category, survey_year) |>
  summarise(
    catch = sum(catch),
    lower = sum(lower),
    upper = sum(upper)
  ) |>
  ungroup() |>
  mutate(waterbody = paste(waterbody, reach_no, sep = ": Reach ")) |>
  ggplot(
    aes(y = catch, x = survey_year, ymin = lower, ymax = upper, col = category)
  ) +
  geom_point(position = position_dodge(0.1)) +
  geom_errorbar(position = position_dodge(0.1)) +
  scale_color_brewer(palette = "Set2", name = "") +
  xlab("Survey year") +
  ylab("Catch") +
  facet_wrap( ~ waterbody, scales = "free_y") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 8)
  )
fitted_1plus_rb <- fitted_obs_plot |>
  left_join(include_subset, by = c("species", "waterbody", "reach_no")) |>
  mutate(category = rep(c("Observed", "Modelled"), each = nrow(mod$data))) |>
  filter(estimated_age == 1, species == "Gadopsis marmoratus") |>
  filter(sum_catch > 1) |>
  group_by(waterbody, reach_no, category, survey_year) |>
  summarise(
    catch = sum(catch),
    lower = sum(lower),
    upper = sum(upper)
  ) |>
  ungroup() |>
  mutate(waterbody = paste(waterbody, reach_no, sep = ": Reach ")) |>
  ggplot(
    aes(y = catch, x = survey_year, ymin = lower, ymax = upper, col = category)
  ) +
  geom_point(position = position_dodge(0.1)) +
  geom_errorbar(position = position_dodge(0.1)) +
  scale_color_brewer(palette = "Set2", name = "") +
  xlab("Survey year") +
  ylab("Catch") +
  facet_wrap( ~ waterbody, scales = "free_y") +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 8)
  )
ggsave(
  file = "outputs/figures/fitted-observed-0plus-bispinosus.png",
  plot = fitted_0plus_ts,
  device = ragg::agg_png,
  width = (5/4) * 9,
  height = 9,
  units = "in",
  dpi = 600
)
ggsave(
  file = "outputs/figures/fitted-observed-1plus-bispinosus.png",
  plot = fitted_1plus_ts,
  device = ragg::agg_png,
  width = (5/4) * 9,
  height = 9,
  units = "in",
  dpi = 600
)
ggsave(
  file = "outputs/figures/fitted-observed-0plus-marmoratus.png",
  plot = fitted_0plus_rb,
  device = ragg::agg_png,
  width = (5/4) * 9,
  height = 9,
  units = "in",
  dpi = 600
)
ggsave(
  file = "outputs/figures/fitted-observed-1plus-marmoratus.png",
  plot = fitted_1plus_rb,
  device = ragg::agg_png,
  width = (5/4) * 9,
  height = 9,
  units = "in",
  dpi = 600
)

# calculate model fit by species
fitted_values <- tibble(
  species = mod$data$species,
  waterbody = mod$data$waterbody,
  fitted = fitted(mod, robust = TRUE)[, "Estimate"],
  observed = mod$data$catch
)
r2 <- fitted_values |>
  group_by(species, waterbody) |>
  summarise(
    r = cor(fitted, observed),
    r2 = r ^ 2, 
    rsp = cor(fitted, observed, method = "spearman"),
    r2sp = rsp ^ 2
  )
write.csv(r2, file = "outputs/tables/r2.csv")
fitted_plot <- fitted_values |>
  ggplot(aes(y = fitted, x = observed, col = waterbody)) +
  geom_point() +
  facet_wrap( ~ species, scales = "free")
ggsave(
  filename = "outputs/figures/fitted-vs-observed.png",
  plot = fitted_plot,
  device = ragg::agg_png,
  width = 10,
  height = 10,
  dpi = 600,
  units = "in"
)

# summarise flow effects
flow_vars <- c(
  "lowflow_wateryear_std",
  "dailyflow_spring_std",
  "dailyflow_summer_std",
  "dailyflow_winter_std",
  "maxantecedent_wateryear_std",
  "cvflow_spawning_std",
  "spawning_temperature_std [-8:2.5 by=0.1]"
)
flow_effects <- vector("list", length = length(flow_vars))
flow_effects <- lapply(
  flow_vars,
  \(x) ggpredict(mod, terms = c(x, "species", "waterbody"), type = "random", ci.lvl = 0.5, allow_new_levels = TRUE) |> as_tibble()
)

# rename to remove the range specification from below
flow_vars[7] <- "spawning_temperature_std"

# plot these
flow_vars_clean <- c(
  "Number of days of low flow",
  "Relative spring flow",
  "Relative summer flow",
  "Relative winter flow",
  "Relative antecedent flow",
  "Flow variability during spawning season",
  "Water temperature during spawning season"
)
rescale_values <- list(
  low_flow = c(mean(cpue$lowflow_wateryear), sd(cpue$lowflow_wateryear)),
  spring_flow = c(mean(cpue$dailyflow_spring), sd(cpue$dailyflow_spring)),
  summer_flow = c(mean(cpue$dailyflow_summer), sd(cpue$dailyflow_summer)),
  winter_flow = c(mean(cpue$dailyflow_winter), sd(cpue$dailyflow_winter)),
  antecedent_flow = c(mean(cpue$maxantecedent_wateryear), sd(cpue$maxantecedent_wateryear)),
  cv_flow = c(mean(cpue$cvflow_spawning), sd(cpue$cvflow_spawning)),
  spawning_temp = c(mean(cpue$spawning_temperature), sd(cpue$spawning_temperature))
)
for (i in seq_along(flow_effects)) {
  
  # quick check to drop out really low values of temps (included to plot Kiewa outputs only)
  flow_effects_tmp <- flow_effects[[i]]
  if (i == 7)
    flow_effects_tmp <- flow_effects_tmp |> filter(x > -5)
  
  # plot it
  p <- plot_flow_effects(
    flow_effects_tmp, 
    dat = cpue,
    xname = flow_vars_clean[i],
    rescale = rescale_values[[i]]
  )
  ggsave(
    filename = paste0("outputs/figures/flow-", flow_vars[i], ".png"),
    plot = p,
    device = ragg::agg_png,
    width = 7,
    height = 5,
    dpi = 600,
    units = "in"
  )
}

# Plots for each system with ml/d as axes.
min_vals <- metrics |>
  select(-water_year) |>
  group_by(waterbody, reach_no) |>
  summarise(across(everything(), .fns = list(min = min), .names = "{.col}")) |>
  mutate(
    type = "min",
    spawning_temperature = ifelse(grepl("Kiewa", waterbody), 8, spawning_temperature)
  )
max_vals <- metrics |>
  select(-water_year) |>
  group_by(waterbody, reach_no) |>
  summarise(across(everything(), .fns = list(max = max), .names = "{.col}")) |>
  mutate(
    type = "max",
    lowflow_wateryear = ifelse(waterbody == "Holland Creek", 10, lowflow_wateryear),
    lowflow_wateryear = ifelse(waterbody == "Glenelg River" & reach_no == 1, 10, lowflow_wateryear),
    lowflow_wateryear = ifelse(waterbody == "Moorabool River" & reach_no == 4, 10, lowflow_wateryear),
    spawning_temperature = ifelse(grepl("Kiewa", waterbody), 13, spawning_temperature)
  )
metrics_range <- bind_rows(min_vals, max_vals)

# and calculate long-term median to rescale flow values
lt_values <- flow |>
  group_by(waterbody, reach_no) |>
  summarise(
    median = median(stream_discharge_mld),
    mean = mean(stream_discharge_mld)
  )

# plot flow effects rescaled to their observation scale (e.g. ML/d)
flow_vars_clean2 <- c(
  "Number of days of flow below Q90",
  "Average spring flow (ML/d)",
  "Average summer flow (ML/d)",
  "Average winter flow (ML/d)",
  "Maximum antecedent flow (ML/d)",
  "Std. dev. of flow during spawning season (ML/d)",
  "Water temperature during spawning season (C)"
)
p <- vector("list", length = length(flow_effects))
for (i in seq_along(flow_effects)) {
  
  # pull out range as select(Waterbody, min_value, max_value)
  metrics_sub <- metrics_range |>
    select(waterbody, reach_no, type, all_of(gsub("_std", "", flow_vars[i]))) |>
    pivot_wider(id_cols = c(waterbody, reach_no), names_from = type, values_from = gsub("_std", "", flow_vars[i])) |>
    rename(Waterbody = waterbody, Reach = reach_no, min_value = min, max_value = max) |>
    mutate(Reach = as.character(Reach))
  
  # and lt_values as select(Waterbody, scale_factor)
  lt_sub <- lt_values |>
    mutate(
      scale_factor = median,
      scale_factor = ifelse(grepl("temperature|lowflow", flow_vars[i]), 1, scale_factor),
      scale_factor = ifelse(grepl("cvflow", flow_vars[i]), mean, scale_factor)
    ) |>
    rename(Waterbody = waterbody, Reach = reach_no)
  
  # exand and include reaches in this output
  flow_sub <- full_join(
    flow_effects[[i]] |>
      select(-conf.low, -conf.high) |>
      rename(Species = group, Waterbody = facet),
    lt_sub |> distinct(Waterbody, Reach),
    by = c("Waterbody"),
    relationship = "many-to-many"
  )
  
  # and plot this
  p[[i]] <- plot_flow_effects_system(
    x = flow_sub, 
    dat = cpue |> mutate(reach_no = as.character(reach_no)),
    metric_range = metrics_sub,
    lt_vals = lt_sub,
    xname = flow_vars_clean2[i],
    rescale = rescale_values[[i]]
  )
  
}

# make up a plot for each system (all variables)
all_systems <- lt_values |> 
  mutate(waterbody_reach = tolower(gsub(" ", "_", paste(waterbody, reach_no, sep = "_r")))) |>
  pull(waterbody_reach) |>
  unique()
for (i in seq_along(all_systems)) {
  
  flow_plots <- lapply(p, \(x, sys) x[[sys]], sys = all_systems[i])
  sys_plot <- (flow_plots[[1]] + flow_plots[[2]]) /
    (flow_plots[[3]] + flow_plots[[4]]) /
    (flow_plots[[5]] + flow_plots[[6]]) /
    (flow_plots[[7]] + plot_spacer())
  ggsave(
    filename = paste0("outputs/figures/flow-", all_systems[i], ".png"),
    plot = sys_plot,
    device = ragg::agg_png,
    width = 8,
    height = 8,
    dpi = 600,
    units = "in"
  )
  
}
