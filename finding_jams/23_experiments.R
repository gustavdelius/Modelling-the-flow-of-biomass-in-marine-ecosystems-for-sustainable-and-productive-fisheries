library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)
library(glue)
library(ggplot2)
library(colorRamps)
library(future)
library(scales) 

make_second_order_params_kr <- function(lambda = 2.05, resource_decrease = 0.001,
                                        capacity_mult = 1, second_order = TRUE,
                                        ext_diff = 0.00) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)
  
  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  
  given_species_params(params)$D_ext <- ext_diff
  params <- setBevertonHolt(params)
  
  # Independent capacity and rate knobs, both as multipliers on mizer's own
  # kappa/lambda-derived defaults -- capacity_mult = 1 reproduces the
  # unswept default carrying capacity, same convention resource_decrease
  # already uses for rate. balance = FALSE so neither gets recomputed from
  # the other.
  default_capacity <- getResourceCapacity(params)
  r  <- getResourceRate(params) * resource_decrease
  cc <- default_capacity * capacity_mult
  
  params <- setResource(params, resource_rate = r, resource_capacity = cc,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE)
  
  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }
  
  params
}

score_oscillation_kr <- function(resource_decrease, capacity_mult,
                                 lambda = 2.05, t_run = 600) {
  p   <- make_second_order_params_kr(lambda = lambda,
                                     resource_decrease = resource_decrease,
                                     capacity_mult     = capacity_mult)
  sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                 progress_bar = FALSE, effort = 0, method = "tr_bdf2")
  bm   <- getBiomass(sim)[, "Anchovy"]
  tv   <- as.numeric(names(bm))
  late <- bm[tv > t_run * 0.6]  # last 40% of the run -- the "settled down" window
  
  data.frame(resource_decrease = resource_decrease,
             capacity_mult     = capacity_mult,
             max_bm            = max(late),
             min_bm            = min(late),
             mean_bm           = mean(late),
             rel_amplitude     = (max(late) - min(late)) / mean(late))
}


rd_seq <- exp(seq(log(0.0001), log(0.5), length.out = 20))

# live_plot = TRUE passes mizerExperimental's biomass_callback to project(),
# same as run_rd_sweep() above -- see the note there. Open x11() (or leave
# the RStudio plot pane visible) before calling this with live_plot = TRUE so
# you can watch each of the 20+20 steps in the sweep settle in real time.
run_rd_sweep_second_order <- function(rd_seq, init_n = NULL, init_n_pp = NULL,
                                      t_run = 60, lambda = 2.05,
                                      second_order = TRUE, label = "",
                                      live_plot = TRUE) {
  out       <- data.frame(rd = rd_seq, max_bm = NA_real_, min_bm = NA_real_,
                          bm_mean = NA_real_)
  state_n   <- init_n
  state_npp <- init_n_pp
  
  for (i in seq_along(rd_seq)) {
    # A fresh params object is still built at each step (resource_decrease
    # changes the resource_rate vector itself, not just a scaling applied
    # after the fact) -- but state_n/state_npp carry the PREVIOUS run's
    # settled abundances into it as the initial condition, so each step
    # starts from where the last one left off rather than from scratch.
    p <- make_second_order_params(lambda = lambda, resource_decrease = rd_seq[i],
                                  second_order = second_order)
    if (!is.null(state_n)) {
      p@initial_n[]    <- state_n
      p@initial_n_pp[] <- state_npp
    }
    if (live_plot) {
      sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                     progress_bar = FALSE, effort = 0, method = "tr_bdf2",
                     callback = biomass_callback, species = "Anchovy")
    } else {
      sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                     progress_bar = FALSE, effort = 0, method = "tr_bdf2")
    }
    bm   <- getBiomass(sim)[, "Anchovy"]
    tv   <- as.numeric(names(bm))
    late <- bm[tv > t_run * 0.6]  # last 40% of the run -- the "settled down" window
    last <- dim(sim@n)[1]
    
    state_n   <- sim@n[last, , ]
    state_npp <- sim@n_pp[last, ]
    
    out$max_bm[i]  <- max(late)
    out$min_bm[i]  <- min(late)
    out$bm_mean[i] <- mean(late)
  }
  list(df = out, n_final = state_n, npp_final = state_npp)
}

rd_seq_so <- rd_seq  # same log-spaced grid as Part 1 (0.0001-0.5, 20 points) -- keeps the two directly comparable



################################################################################
# Phase diagram: carrying capacity x resource rate, traced as a snake
#
# Straight from the supercooling/superheating discussion: a hysteresis loop
# only shows up if the system carries memory across the sweep, and the
# *path* taken through a 2-parameter space can change the observed
# metastable region (sweeping field at fixed temperature vs temperature at
# fixed field gives different loop widths in that literature). So far every
# sweep here has threaded state along ONE axis at a time (resource_decrease),
# restarting fresh at each capacity_mult row. This instead traces the WHOLE
# 2D grid as a single continuous path -- a boustrophedon ("snake"): sweep
# resource_decrease ascending across one capacity_mult row, then descending
# across the next, so consecutive points -- including the row-to-row
# turn -- are always adjacent in parameter space with a carried-over state
# between them, the way an experimentalist sweeps field back and forth while
# slowly stepping temperature rather than re-preparing a fresh sample at
# every (T, H) pair.
#
# This can't be parallelised with future.apply like the earlier sweeps --
# each point depends on the previous one's settled state -- so it's a plain
# sequential loop, and will take longer wall-clock than a same-size grid run
# in parallel.
################################################################################

# Same threshold as Day 22's baseline sweep -- filters out numerical noise
# masquerading as a viable population before classifying a grid point as
# Collapsed/Oscillating/Fixed point below. Defined here (not just inherited
# from 22_experiments.R's session state) so this file runs standalone.
MIN_VIABLE_BIOMASS <- 1e-2

run_snake_grid <- function(rd_grid, cc_grid, lambda = 2.05, t_run = 600) {
  out       <- vector("list", length(rd_grid) * length(cc_grid))
  state_n   <- NULL
  state_npp <- NULL
  k <- 1
  
  for (j in seq_along(cc_grid)) {
    # Odd rows sweep resource_decrease ascending, even rows descending --
    # the snake turn at each row boundary keeps the path continuous in
    # resource_decrease across the capacity_mult step too, rather than
    # jumping from one end of the grid back to the other.
    rd_order <- if (j %% 2 == 1) rd_grid else rev(rd_grid)
    
    for (rd in rd_order) {
      p <- make_second_order_params_kr(lambda = lambda, resource_decrease = rd,
                                       capacity_mult = cc_grid[j])
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }
      sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                     progress_bar = FALSE, effort = 0, method = "tr_bdf2")
      bm   <- getBiomass(sim)[, "Anchovy"]
      tv   <- as.numeric(names(bm))
      late <- bm[tv > t_run * 0.6]
      last <- dim(sim@n)[1]
      
      state_n   <- sim@n[last, , ]
      state_npp <- sim@n_pp[last, ]
      
      out[[k]] <- data.frame(resource_decrease = rd, capacity_mult = cc_grid[j],
                             max_bm        = max(late), min_bm = min(late),
                             mean_bm       = mean(late),
                             rel_amplitude = (max(late) - min(late)) / mean(late),
                             step          = k)
      k <- k + 1
    }
  }
  bind_rows(out)
}

# 5 resource_decrease points x 4 capacity_mult rows = 20 runs, same budget as
# the earlier sweeps, just threaded as one continuous path instead of an
# independent 4x5 grid.
rd_grid_snake <- exp(seq(log(0.0001), log(0.5), length.out = 5))
cc_grid_snake <- exp(seq(log(3), log(100), length.out = 4))

snake_df <- run_snake_grid(rd_grid_snake, cc_grid_snake) %>%
  arrange(step) %>%
  mutate(collapsed = mean_bm < MIN_VIABLE_BIOMASS,
         phase     = case_when(
           collapsed            ~ "Collapsed",
           rel_amplitude > 1e-2 ~ "Oscillating",
           TRUE                 ~ "Fixed point"
         ))

# Phase-diagram-style plot: discrete regions (like a solid/liquid/gas
# diagram) rather than a continuous colour scale, with the snake path
# overlaid as a thin guide line so the traversal order -- and therefore
# where state was carried across a boundary versus where it wasn't -- stays
# visible on the same plot. Note this single path is only one route through
# the 2D space; per the path-dependence point above, a different route
# (e.g. capacity_mult swept at fixed resource_decrease) could trace out
# different boundaries for the same underlying regions.
snake_plot_fwd <- ggplot(snake_df, aes(x = resource_decrease, y = capacity_mult)) +
  geom_tile(aes(fill = phase), color = "white", linewidth = 0.3) +
  geom_path(color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c("Fixed point" = "#4C72B0",
                               "Oscillating" = "#DD8452",
                               "Collapsed"   = "grey70")) +
  labs(x = "resource_decrease (resource rate)",
       y = "capacity_mult (carrying capacity)",
       title = "Phase diagram: fixed point vs. oscillating regime",
       subtitle = "Traced as a snake path -- state carried continuously between neighbouring points",
       fill = "Regime") +
  theme_minimal()
snake_plot_fwd

################################################################################
# Saving plots overnight
#
# This is set to run unattended for a while, so every plot below gets written
# to disk as soon as it's produced instead of only sitting in the graphics
# device (or getting silently overwritten by the next plot before anyone's
# looked at it). interesting_plots/ already exists at the project root with
# a handful of PNGs saved by hand earlier -- dir.create() below just makes
# sure it exists if this is ever run from a fresh checkout; it's silent
# (showWarnings = FALSE) if the folder is already there.
################################################################################

dir.create("interesting_plots", showWarnings = FALSE)

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 150) {
  ggsave(file.path("interesting_plots", filename), plot = plot,
         width = width, height = height, dpi = dpi)
}

save_plot(snake_plot_fwd, "Phase diagram - snake, capacity rows.png")

################################################################################
# The other intuitive routes: reversed, and axis order swapped
#
# Two more paths through the same (resource_decrease, capacity_mult) grid,
# both motivated by the path-dependence point from the hysteresis research:
#
#   1. The exact reverse of the snake above (same points, same order, just
#      walked backwards) -- the direct analogue of every fwd/bwd pair used
#      throughout this project (rev(rd_seq) etc.), applied to a 2D path
#      instead of a 1D one.
#   2. The transposed route: capacity_mult swept back and forth within each
#      resource_decrease "row" instead of the other way around -- literally
#      the "field at fixed T vs T at fixed field" comparison from the
#      literature, translated into this grid.
#
# build_path_*() below just returns the sequence of (resource_decrease,
# capacity_mult) points to visit, in order; run_along_path() is the same
# project()-and-carry-state loop as run_snake_grid() above, but takes that
# sequence directly instead of reconstructing it from rd_grid/cc_grid, so the
# exact-reverse case can be expressed as literally reversing a data frame.
################################################################################

build_path_capacity_rows <- function(rd_grid, cc_grid) {
  path <- vector("list", length(rd_grid) * length(cc_grid))
  k <- 1
  for (j in seq_along(cc_grid)) {
    rd_order <- if (j %% 2 == 1) rd_grid else rev(rd_grid)
    for (rd in rd_order) {
      path[[k]] <- data.frame(resource_decrease = rd, capacity_mult = cc_grid[j])
      k <- k + 1
    }
  }
  bind_rows(path)
}

# The transposed route: resource_decrease values are the "rows", and
# capacity_mult is swept back and forth within each one.
build_path_rd_rows <- function(rd_grid, cc_grid) {
  path <- vector("list", length(rd_grid) * length(cc_grid))
  k <- 1
  for (i in seq_along(rd_grid)) {
    cc_order <- if (i %% 2 == 1) cc_grid else rev(cc_grid)
    for (cc in cc_order) {
      path[[k]] <- data.frame(resource_decrease = rd_grid[i], capacity_mult = cc)
      k <- k + 1
    }
  }
  bind_rows(path)
}

run_along_path <- function(path_df, lambda = 2.05, t_run = 600, label = "") {
  state_n   <- NULL
  state_npp <- NULL
  out <- vector("list", nrow(path_df))
  
  for (i in seq_len(nrow(path_df))) {
    rd <- path_df$resource_decrease[i]
    cc <- path_df$capacity_mult[i]
    
    p <- make_second_order_params_kr(lambda = lambda, resource_decrease = rd,
                                     capacity_mult = cc)
    if (!is.null(state_n)) {
      p@initial_n[]    <- state_n
      p@initial_n_pp[] <- state_npp
    }
    sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                   progress_bar = FALSE, effort = 0, method = "tr_bdf2")
    bm   <- getBiomass(sim)[, "Anchovy"]
    tv   <- as.numeric(names(bm))
    late <- bm[tv > t_run * 0.6]
    last <- dim(sim@n)[1]
    
    state_n   <- sim@n[last, , ]
    state_npp <- sim@n_pp[last, ]
    
    out[[i]] <- data.frame(resource_decrease = rd, capacity_mult = cc,
                           max_bm        = max(late), min_bm = min(late),
                           mean_bm       = mean(late),
                           rel_amplitude = (max(late) - min(late)) / mean(late),
                           step = i, route = label)
  }
  bind_rows(out)
}

path_fwd         <- build_path_capacity_rows(rd_grid_snake, cc_grid_snake)
path_bwd         <- path_fwd[rev(seq_len(nrow(path_fwd))), ]  # exact reverse retrace
path_transposed  <- build_path_rd_rows(rd_grid_snake, cc_grid_snake)

routes_df <- bind_rows(
  run_along_path(path_fwd,        label = "Forward (capacity rows)"),
  run_along_path(path_bwd,        label = "Reversed (exact retrace)"),
  run_along_path(path_transposed, label = "Transposed (resource-rate rows)")
) %>%
  mutate(collapsed = mean_bm < MIN_VIABLE_BIOMASS,
         phase     = case_when(
           collapsed            ~ "Collapsed",
           rel_amplitude > 1e-2 ~ "Oscillating",
           TRUE                 ~ "Fixed point"
         ))

# Saved immediately in case the overnight run doesn't get all the way
# through -- whatever routes finished are recoverable without rerunning.
saveRDS(routes_df, file.path("interesting_plots", "routes_df.rds"))

routes_plot <- ggplot(routes_df, aes(x = resource_decrease, y = capacity_mult)) +
  geom_tile(aes(fill = phase), color = "white", linewidth = 0.3) +
  geom_path(color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c("Fixed point" = "#4C72B0",
                               "Oscillating" = "#DD8452",
                               "Collapsed"   = "grey70")) +
  facet_wrap(~route) +
  labs(x = "resource_decrease (resource rate)",
       y = "capacity_mult (carrying capacity)",
       title = "Phase diagram by route: same grid, three different paths",
       subtitle = "If these three panels disagree, the boundary is path-dependent, not just a fixed property of the parameters",
       fill = "Regime") +
  theme_minimal()
routes_plot

save_plot(routes_plot, "Phase diagram - three routes compared.png", width = 12, height = 5)

################################################################################
# Follow-up 1: scan-rate dependence
#
# The three routes above agree exactly at every one of the 20 grid points --
# but t_run = 600 per step is long enough for each point to fully settle
# before being sampled, which is also exactly the condition under which any
# route-dependent memory effect would be erased before it could be measured.
# Rerunning the same three routes at shorter t_run per step is the direct
# test: if disagreement shows up (and grows) as t_run shrinks, that confirms
# a real but scan-rate-suppressed hysteresis; if the routes keep agreeing
# even at short t_run, the "no path dependence" read above holds up rather
# than being an artefact of over-relaxing every step.
#
# t_run = 600 doesn't need rerunning -- routes_df above already has it.
################################################################################

classify_phase <- function(mean_bm, rel_amplitude) {
  case_when(
    mean_bm < MIN_VIABLE_BIOMASS ~ "Collapsed",
    rel_amplitude > 1e-2         ~ "Oscillating",
    TRUE                         ~ "Fixed point"
  )
}

run_three_routes <- function(t_run, rd_grid = rd_grid_snake, cc_grid = cc_grid_snake) {
  fwd        <- build_path_capacity_rows(rd_grid, cc_grid)
  bwd        <- fwd[rev(seq_len(nrow(fwd))), ]
  transposed <- build_path_rd_rows(rd_grid, cc_grid)

  bind_rows(
    run_along_path(fwd,        t_run = t_run, label = "Forward (capacity rows)"),
    run_along_path(bwd,        t_run = t_run, label = "Reversed (exact retrace)"),
    run_along_path(transposed, t_run = t_run, label = "Transposed (resource-rate rows)")
  ) %>%
    mutate(t_run = t_run, phase = classify_phase(mean_bm, rel_amplitude))
}

# t_run = 50 and 150 run fresh; t_run = 600 is reused from routes_df above --
# same grid, same routes, same MIN_VIABLE_BIOMASS floor, directly comparable
# without spending another ~45 minutes reproducing it.
scan_rate_new <- bind_rows(
  run_three_routes(t_run = 50),
  run_three_routes(t_run = 150)
)

scan_rate_df <- bind_rows(scan_rate_new, routes_df %>% mutate(t_run = 600))
saveRDS(scan_rate_df, file.path("interesting_plots", "scan_rate_df.rds"))

# Disagreement count at each scan rate -- the actual test, not just the plot.
scan_rate_agreement <- scan_rate_df %>%
  select(resource_decrease, capacity_mult, t_run, route, phase) %>%
  tidyr::pivot_wider(names_from = route, values_from = phase) %>%
  mutate(agree = `Forward (capacity rows)` == `Reversed (exact retrace)` &
                 `Reversed (exact retrace)` == `Transposed (resource-rate rows)`) %>%
  group_by(t_run) %>%
  summarise(n_disagree = sum(!agree), n_total = n(), .groups = "drop")

print(scan_rate_agreement)
write.csv(scan_rate_agreement, file.path("interesting_plots", "scan_rate_agreement.csv"),
          row.names = FALSE)

scan_rate_plot <- ggplot(scan_rate_df, aes(x = resource_decrease, y = capacity_mult)) +
  geom_tile(aes(fill = phase), color = "white", linewidth = 0.3) +
  geom_path(color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c("Fixed point" = "#4C72B0",
                               "Oscillating" = "#DD8452",
                               "Collapsed"   = "grey70")) +
  facet_grid(t_run ~ route) +
  labs(x = "resource_decrease (resource rate)",
       y = "capacity_mult (carrying capacity)",
       title = "Scan-rate dependence: same three routes at t_run = 50, 150, 600",
       subtitle = "Rows = per-step settle time. If low-t_run rows disagree across routes but t_run=600 doesn't, relaxation time was hiding real path dependence.",
       fill = "Regime") +
  theme_minimal()
scan_rate_plot

save_plot(scan_rate_plot, "Scan-rate dependence, three routes.png", width = 12, height = 10)

################################################################################
# Follow-up 2: genuine perturbation retest of the Fixed-point cells
#
# Every sweep since Day 20 (including everything above) starts each project()
# call from mizer's own default initial state, or from whatever the previous
# grid point settled to. Day 17/18's make_limit_cycle_sim() found the limit
# cycle in the first place by actively kicking the system away from the fixed
# point -- depleting the resource to 10% of carrying capacity, a short
# predictor-corrector burn-in, then dividing the mature stock (w >= w_mat) by
# a large factor -- before continuing. A cell that reads "Fixed point" in
# every sweep so far might just mean nothing tried has escaped it, not that
# no limit cycle exists there. This reruns that exact two-stage kick (Day 18's
# version, adapted to make_second_order_params_kr / tr_bdf2) on every cell the
# t_run = 600 sweep above classified as a Fixed point.
################################################################################

make_limit_cycle_sim_kr <- function(params, t_total = 600, effort = 0, perturbation = 1e3) {
  params@initial_n_pp[] <- params@cc_pp * 0.1

  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / perturbation

  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "tr_bdf2")
}

fixed_point_cells <- routes_df %>%
  filter(phase == "Fixed point") %>%
  distinct(resource_decrease, capacity_mult)

print(fixed_point_cells)

retest_with_kick <- function(rd, cc, lambda = 2.05, t_total = 600, perturbation = 1e3) {
  p   <- make_second_order_params_kr(lambda = lambda, resource_decrease = rd,
                                     capacity_mult = cc)
  sim <- make_limit_cycle_sim_kr(p, t_total = t_total, perturbation = perturbation)
  bm   <- getBiomass(sim)[, "Anchovy"]
  tv   <- as.numeric(names(bm))
  late <- bm[tv > (t_total - 10) * 0.6]

  data.frame(resource_decrease = rd, capacity_mult = cc,
             mean_bm       = mean(late),
             rel_amplitude = (max(late) - min(late)) / mean(late))
}

kick_retest_df <- purrr::map2_dfr(
  fixed_point_cells$resource_decrease, fixed_point_cells$capacity_mult,
  retest_with_kick
) %>%
  mutate(phase_after_kick = classify_phase(mean_bm, rel_amplitude)) %>%
  left_join(routes_df %>% select(resource_decrease, capacity_mult, phase) %>%
              distinct() %>%
              rename(phase_before_kick = phase),
            by = c("resource_decrease", "capacity_mult"))

print(kick_retest_df)
write.csv(kick_retest_df, file.path("interesting_plots", "kick_retest_df.csv"),
          row.names = FALSE)

################################################################################
# Follow-up 3: densified grid around the phase boundaries
#
# The 5x4 grid used above is coarse enough that a narrow metastable band
# could sit entirely inside one tile and never show up. This doubles the
# resolution on both axes (5x4 -> 9x7, still log-spaced over the same range)
# and reruns just the Forward route -- Follow-up 1 above already covers
# route-to-route agreement, so this is purely about whether a finer grid
# reveals structure the coarse one smoothed over. t_run is dropped to 300
# (vs 600 above) to keep the 63-point run to a reasonable length; per
# Follow-up 1's own logic that trades some settle time for resolution, so if
# this run's boundaries look different from the grid above, check t_run
# before concluding it's purely a resolution effect.
################################################################################

rd_grid_fine <- exp(seq(log(0.0001), log(0.5), length.out = 9))
cc_grid_fine <- exp(seq(log(3), log(100), length.out = 7))

path_fwd_fine <- build_path_capacity_rows(rd_grid_fine, cc_grid_fine)

fine_grid_df <- run_along_path(path_fwd_fine, t_run = 300,
                               label = "Forward (capacity rows), fine grid") %>%
  mutate(phase = classify_phase(mean_bm, rel_amplitude))

saveRDS(fine_grid_df, file.path("interesting_plots", "fine_grid_df.rds"))

fine_grid_plot <- ggplot(fine_grid_df, aes(x = resource_decrease, y = capacity_mult)) +
  geom_tile(aes(fill = phase), color = "white", linewidth = 0.2) +
  geom_path(color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c("Fixed point" = "#4C72B0",
                               "Oscillating" = "#DD8452",
                               "Collapsed"   = "grey70")) +
  labs(x = "resource_decrease (resource rate)",
       y = "capacity_mult (carrying capacity)",
       title = "Densified phase diagram: 9x7 grid, Forward route only",
       subtitle = "t_run = 300 per step -- double the resolution of the 5x4 grid above",
       fill = "Regime") +
  theme_minimal()
fine_grid_plot

save_plot(fine_grid_plot, "Phase diagram - densified 9x7 grid.png", width = 9, height = 7)

################################################################################
# Follow-up 4: refining capacity_mult in [3, 10]
#
# The 9x7 densified grid above found a wide Fixed-point plateau at
# capacity_mult = 5.38 that the original 5x4 grid's cap=3 / cap=9.65 rows
# didn't show -- the coarse grid made the cap=3 -> cap=9.65 transition look
# like a fast, near-immediate switch from mostly-fixed-point to
# mostly-oscillating, when there's actually a broad fixed-point band in
# between. This narrows in on capacity_mult in [3, 10] specifically -- same
# window, denser -- to see how wide that plateau actually is and where its
# edges sit, using the same 9-point resource_decrease grid as the 9x7 run
# for direct comparability. Forward route only, t_run = 300 again.
################################################################################

cc_grid_refine <- exp(seq(log(3), log(10), length.out = 9))

path_fwd_refine <- build_path_capacity_rows(rd_grid_fine, cc_grid_refine)

refine_grid_df <- run_along_path(path_fwd_refine, t_run = 300,
                                 label = "Forward (capacity rows), cap in [3,10]") %>%
  mutate(phase = classify_phase(mean_bm, rel_amplitude))

saveRDS(refine_grid_df, file.path("interesting_plots", "refine_grid_df.rds"))

refine_grid_plot <- ggplot(refine_grid_df, aes(x = resource_decrease, y = capacity_mult)) +
  geom_tile(aes(fill = phase), color = "white", linewidth = 0.2) +
  geom_path(color = "black", linewidth = 0.3, alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c("Fixed point" = "#4C72B0",
                               "Oscillating" = "#DD8452",
                               "Collapsed"   = "grey70")) +
  labs(x = "resource_decrease (resource rate)",
       y = "capacity_mult (carrying capacity)",
       title = "Refined phase diagram: capacity_mult in [3, 10], 9x9 grid",
       subtitle = "Zooming in on the Fixed-point plateau found at capacity_mult ~ 5.4 in the 9x7 grid",
       fill = "Regime") +
  theme_minimal()
refine_grid_plot

save_plot(refine_grid_plot, "Phase diagram - refined capacity 3-10.png", width = 8, height = 6)

################################################################################
# Follow-up 5: bifurcation sweeps revisited -- was Day 22's "no hysteresis"
# result just a diffusion artefact?
#
# Day 22 found the forward/backward hysteresis gap that Days 18-21 built the
# whole limit-cycle story around disappeared once the resource setup was
# fixed (balance = FALSE, resource_capacity/resource_rate set directly) AND
# ext_diff was defaulted to 0 rather than Day 21's 0.01. Which of those two
# changes actually killed the hysteresis was explicitly left unresolved at
# the time ("disentangling the two... [is a] separate problem for later").
# This picks that apart along two lines:
#
#   1. Classic 1D bifurcation diagrams (forward + backward branches) on the
#      *current* balance = FALSE setup, once sweeping resource_decrease with
#      capacity_mult held fixed, once sweeping capacity_mult with
#      resource_decrease held fixed -- the direct capacity_mult analogue of
#      every resource_decrease bifurcation diagram run since Day 18.
#   2. The *old* balance = TRUE, resource_level = 1 setup (Day 21/22's
#      original, pre-fix resource handling) with resource_decrease swept and
#      ext_diff varied across 4 equally-spaced values from 0.001 to 0.5 --
#      directly testing whether hysteresis reappears at some diffusion
#      strength under the old setup, which would mean the earlier "no
#      hysteresis" result really was a diffusion-strength artefact rather
#      than a consequence of the balance/resource_level fix itself.
################################################################################

make_second_order_params_balanced <- function(lambda = 2.05, resource_decrease = 0.001,
                                              resource_level = 1, second_order = TRUE,
                                              ext_diff = 0.00) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )

  given_species_params(params)$D_ext <- ext_diff
  params <- setBevertonHolt(params)

  # setResource() rejects resource_rate and resource_level together under
  # balance = TRUE ("you should only provide either... because the other is
  # determined by the requirement that the resource replenishes at the same
  # rate at which it is consumed") -- confirmed directly against the
  # installed mizer build, not assumed. resource_level is kept as an
  # argument here (default 1, i.e. the resource left at its full,
  # unsuppressed level) purely to document intent; only resource_rate is
  # actually passed, which is exactly Day 21's original balance = TRUE call
  # before Day 22 moved to balance = FALSE + explicit resource_level = 0.5.
  # That's the version this function is meant to reproduce.
  r <- getResourceRate(params) * resource_decrease
  params <- setResource(params, resource_rate = r,
                        resource_dynamics = "resource_semichemostat",
                        balance = TRUE)

  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }

  params
}

# Generic forward/backward bifurcation sweep over any one parameter of
# params_fn, holding everything else in fixed_params constant -- the same
# forward-then-backward, state-carried-between-steps convention every
# bifurcation diagram since Day 18 has used, just generalised to sweep
# capacity_mult as easily as resource_decrease.
run_bifurcation_sweep <- function(param_seq, param_name, fixed_params = list(),
                                  params_fn = make_second_order_params_kr,
                                  t_run = 600, lambda = 2.05) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_bm = NA_real_, min_bm = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      args <- fixed_params
      args[[param_name]] <- seq_vals[i]
      args$lambda <- lambda
      p <- do.call(params_fn, args)
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }
      sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                     progress_bar = FALSE, effort = 0, method = "tr_bdf2")
      bm   <- getBiomass(sim)[, "Anchovy"]
      tv   <- as.numeric(names(bm))
      late <- bm[tv > t_run * 0.6]
      last <- dim(sim@n)[1]

      state_n   <- sim@n[last, , ]
      state_npp <- sim@n_pp[last, ]

      out$max_bm[i] <- max(late)
      out$min_bm[i] <- min(late)
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(param_seq)
  bwd    <- run_one_direction(rev(param_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, biomass = fwd$df$max_bm, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, biomass = fwd$df$min_bm, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, biomass = bwd_df$max_bm, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, biomass = bwd_df$min_bm, direction = "Backward", branch = "min")
  )
}

plot_bifurcation <- function(df, x_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = biomass, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_x_log10() +
    labs(x = x_label, y = "Biomass", title = title, subtitle = subtitle) +
    theme_minimal()
}

################################################################################
# 1. balance = FALSE (current setup): resource_decrease swept vs.
#    capacity_mult swept, each holding the other fixed at its own baseline
#    default (capacity_mult = 1 reproduces the unswept carrying capacity;
#    resource_decrease = 0.001 is the value hardcoded throughout Days 17-19).
################################################################################

rd_seq_bif <- exp(seq(log(0.0001), log(0.5), length.out = 20))
cc_seq_bif <- exp(seq(log(1), log(100), length.out = 20))

bif_rd_df <- run_bifurcation_sweep(rd_seq_bif, "resource_decrease",
                                   fixed_params = list(capacity_mult = 1),
                                   params_fn = make_second_order_params_kr)

bif_cc_df <- run_bifurcation_sweep(cc_seq_bif, "capacity_mult",
                                   fixed_params = list(resource_decrease = 0.001),
                                   params_fn = make_second_order_params_kr)

saveRDS(list(rd = bif_rd_df, cc = bif_cc_df),
        file.path("interesting_plots", "bifurcation_kr_df.rds"))

bif_rd_plot <- plot_bifurcation(bif_rd_df, "resource_decrease",
                                "Bifurcation: resource_decrease swept",
                                "capacity_mult held at 1 (balance = FALSE)")
bif_cc_plot <- plot_bifurcation(bif_cc_df, "capacity_mult",
                                "Bifurcation: capacity_mult swept",
                                "resource_decrease held at 0.001 (balance = FALSE)")

bif_kr_combined <- bif_rd_plot | bif_cc_plot
bif_kr_combined

save_plot(bif_kr_combined, "Bifurcation - resource_decrease vs capacity_mult (balance=FALSE).png",
         width = 12, height = 5)

################################################################################
# 2. balance = TRUE, resource_level = 1 (the old, pre-Day-22-fix setup):
#    resource_decrease swept only -- resource_level isn't an independent
#    capacity control (Day 22's own finding), so there's no capacity_mult
#    analogue to sweep here -- at 4 equally-spaced ext_diff values between
#    0.001 and 0.5. If hysteresis reappears at some diffusion strength here,
#    Day 22's "no hysteresis" result was a diffusion artefact of that
#    specific setup, not evidence against balance = TRUE / resource_level
#    itself.
################################################################################

ext_diff_values_bif <- seq(0.001, 0.5, length.out = 4)

bif_balanced_df <- bind_rows(lapply(ext_diff_values_bif, function(ed) {
  run_bifurcation_sweep(rd_seq_bif, "resource_decrease",
                        fixed_params = list(resource_level = 1, ext_diff = ed),
                        params_fn = make_second_order_params_balanced) %>%
    mutate(ext_diff = ed)
})) %>%
  mutate(ext_diff_label = factor(sprintf("ext_diff = %.3f", ext_diff),
                                 levels = sprintf("ext_diff = %.3f", sort(unique(ext_diff)))))

saveRDS(bif_balanced_df, file.path("interesting_plots", "bifurcation_balanced_df.rds"))

bif_balanced_plot <- ggplot(bif_balanced_df,
                            aes(x = value, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  scale_x_log10() +
  facet_wrap(~ext_diff_label, nrow = 1) +
  labs(x = "resource_decrease", y = "Biomass",
       title = "Bifurcation: resource_decrease swept at 4 diffusion strengths",
       subtitle = "balance = TRUE, resource_level = 1 -- the pre-Day-22-fix resource setup") +
  theme_minimal()
bif_balanced_plot

save_plot(bif_balanced_plot, "Bifurcation - balanced setup across diffusion strengths.png",
         width = 16, height = 5)
