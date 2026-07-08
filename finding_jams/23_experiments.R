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
