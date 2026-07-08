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
#Will be playing around with different ext_diff
make_second_order_params <- function(lambda = 2.05, resource_decrease = 0.001,
                                     second_order = TRUE, ext_diff = 0) {
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
  
  r      <- getResourceRate(params) * resource_decrease
  params <- setResource(params, resource_rate = r,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE,resource_level = 1)
  
  # Mizer 3.1's second-order finite-volume scheme -- see NOTE above. Both
  # sub-options (bin_average + centred flux) enabled together for the fully
  # second-order variant, rather than just one of the two independently.
  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }
  
  params
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


bwd_so <- run_rd_sweep_second_order(rev(rd_seq_so),
                                    label     = "BWD (2nd order)", live_plot = TRUE)

fwd_so <- run_rd_sweep_second_order(rd_seq_so,
                                    init_n    = bwd_so$n_final,
                                    init_n_pp = bwd_so$npp_final,
                                    label     = "FWD (2nd order)", live_plot = TRUE)


fwd_so_df <- fwd_so$df
bwd_so_df <- bwd_so$df[order(bwd_so$df$rd), ]

# Same max/min-separated overlay as Part 1, now on the second-order +
# ext_diffusion params.
plot(fwd_so_df$rd, fwd_so_df$max_bm, type = "l", col = "steelblue", lwd = 2, lty = 1,
     log = "x", xlab = "resource_decrease", ylab = "Biomass",
     main = "Hysteresis (second-order scheme + ext_diffusion): max/min separately",
     ylim = range(c(fwd_so_df$max_bm, fwd_so_df$min_bm, bwd_so_df$max_bm, bwd_so_df$min_bm)))
lines(fwd_so_df$rd, fwd_so_df$min_bm, col = "steelblue", lwd = 2, lty = 2)
lines(bwd_so_df$rd, bwd_so_df$max_bm, col = "firebrick", lwd = 2, lty = 1)
lines(bwd_so_df$rd, bwd_so_df$min_bm, col = "firebrick", lwd = 2, lty = 2)
legend("topright",
       legend = c("Forward max", "Forward min", "Backward max", "Backward min"),
       col = c("steelblue", "steelblue", "firebrick", "firebrick"),
       lty = c(1, 2, 1, 2), lwd = 2, cex = 0.8)

# Same data, tidier bifurcation-diagram-style ggplot version as Part 1.
bifurcation_df_so <- bind_rows(
  data.frame(rd = fwd_so_df$rd, biomass = fwd_so_df$max_bm, direction = "Forward", branch = "max"),
  data.frame(rd = fwd_so_df$rd, biomass = fwd_so_df$min_bm, direction = "Forward", branch = "min"),
  data.frame(rd = bwd_so_df$rd, biomass = bwd_so_df$max_bm, direction = "Backward", branch = "max"),
  data.frame(rd = bwd_so_df$rd, biomass = bwd_so_df$min_bm, direction = "Backward", branch = "min")
)

ggplot(bifurcation_df_so, aes(x = rd, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_log10() +
  labs(x = "resource_decrease", y = "Biomass",
       title = "Bifurcation diagram: second-order scheme + ext_diffusion of 0.001 ",
       subtitle = "Branches collapsing onto one curve = fixed point; fanning apart = limit cycle") +
  theme_minimal()


################################################################################
# Reintroducing oscillations: sweeping resource_level x resource_decrease
#
# The balance = FALSE, resource_level = 0.5 fix above killed the limit-cycle
# behaviour along with the backwards extinction -- fixing the resource at a
# static fraction of its steady state removes the resource-consumer feedback
# the cycles were probably riding on. resource_level is effectively a
# carrying-capacity knob, and resource_decrease is effectively a resource
# regeneration-rate knob; the classic paradox-of-enrichment intuition is that
# cycles show up at some intermediate combination of the two rather than at
# either extreme, so this scans both together instead of one at a time.
################################################################################


# Same params function as above, but with resource_level exposed as an
# argument instead of hardcoded to 0.5, so it can be swept like
# resource_decrease already is.
make_second_order_params_rl <- function(lambda = 2.05, resource_decrease = 0.001,
                                        resource_level = 0.5, second_order = TRUE,
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

  r      <- getResourceRate(params) * resource_decrease
  params <- setResource(params, resource_rate = r,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE, resource_level = resource_level)

  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }

  params
}

# Oscillation score for one (resource_decrease, resource_level) combination:
# relative amplitude (max - min) / mean over the settled window, rather than
# raw (max - min), so a combination that just sits at a higher biomass level
# doesn't automatically score as more "oscillatory" than one that cycles
# around a smaller mean.
score_oscillation <- function(resource_decrease, resource_level,
                              lambda = 2.05, t_run = 600) {
  p   <- make_second_order_params_rl(lambda = lambda,
                                     resource_decrease = resource_decrease,
                                     resource_level    = resource_level)
  sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                 progress_bar = FALSE, effort = 0, method = "tr_bdf2")
  bm   <- getBiomass(sim)[, "Anchovy"]
  tv   <- as.numeric(names(bm))
  late <- bm[tv > t_run * 0.6]  # last 40% of the run -- the "settled down" window

  data.frame(resource_decrease = resource_decrease,
             resource_level    = resource_level,
             max_bm            = max(late),
             min_bm            = min(late),
             mean_bm           = mean(late),
             rel_amplitude     = (max(late) - min(late)) / mean(late))
}

rd_grid    <- exp(seq(log(0.0001), log(0.5), length.out = 8))
rl_grid    <- seq(0.1, 0.9, by = 0.2)
param_grid <- expand.grid(resource_decrease = rd_grid, resource_level = rl_grid)

# 48 x 5 = 40 independent project() calls, each t_max = 600 -- still more
# expensive than the single-parameter sweeps above, so this runs across cores
# via future.apply rather than one combination at a time.
plan(multisession)
osc_results <- future_lapply(seq_len(nrow(param_grid)), function(i) {
  score_oscillation(param_grid$resource_decrease[i], param_grid$resource_level[i])
}, future.seed = TRUE)
osc_df <- bind_rows(osc_results)
plan(sequential)

# rel_amplitude blows up on floating-point noise once a population has
# collapsed to (numerical) extinction: max_bm/min_bm sitting near
# double-precision zero can differ from each other by many orders of
# magnitude without meaning anything biologically -- that's what was
# happening at resource_level = 0.9/0.7 above (mean_bm as low as 1e-55,
# versus the ~0.01-0.65 scale Day 22's baseline sweep actually sits in).
# MIN_VIABLE_BIOMASS filters those combinations out before ranking, rather
# than letting numerical noise win the "most oscillatory" search.
MIN_VIABLE_BIOMASS <- 1e-2
osc_df <- osc_df %>%
  mutate(collapsed            = mean_bm < MIN_VIABLE_BIOMASS,
         rel_amplitude_viable = ifelse(collapsed, NA_real_, rel_amplitude))

# Heatmap: resource_decrease (log x) against resource_level (y), coloured by
# relative amplitude -- collapsed (functionally extinct) combinations are
# greyed out via NA rather than left in to dominate the colour scale with
# floating-point noise.
ggplot(osc_df, aes(x = resource_decrease, y = resource_level, fill = rel_amplitude_viable)) +
  geom_tile() +
  scale_x_log10() +
  scale_fill_viridis_c(name = "Relative\namplitude", na.value = "grey85") +
  labs(x = "resource_decrease", y = "resource_level",
       title = "Searching for the oscillatory regime",
       subtitle = "Brighter = bigger (max - min) / mean; grey = population collapsed (mean_bm < 1e-2)") +
  theme_minimal()

# Numeric readout of the top viable (non-collapsed) candidates, to pair with
# the heatmap -- filters out the numerically-extinct rows rather than
# re-ranking noise.
osc_df %>%
  filter(!collapsed) %>%
  arrange(desc(rel_amplitude)) %>%
  head(10)


p <- make_second_order_params_rl(resource_level = 0.7,resource_decrease=0.5)
sim <- project(p,t_max=600,method="predictor-corrector")
plotHover(getBiomass(sim,tlim=c(550,600)))

################################################################################
# Reintroducing oscillations, take 2: sweeping carrying capacity x resource rate
#
# The resource_level x resource_decrease sweep above didn't turn up an
# oscillatory regime. resource_level (with balance = FALSE) only
# reparametrises capacity from the CURRENT resource density at the moment
# setResource() runs -- it doesn't give independent control over capacity and
# rate the way a genuine enrichment-vs-regeneration-rate scan needs. This
# sweeps resource_capacity and resource_rate directly and independently
# instead (still balance = FALSE, so neither gets silently recomputed from
# the other): capacity as a multiplier on the default carrying capacity mizer
# derives from kappa/lambda, rate the same way resource_decrease already
# does. capacity_mult = 1, resource_decrease = 0.001 is the closest point in
# this grid to Day 20/21's original (oscillatory, balance = TRUE) setup, so
# it's a useful sanity-check corner to look at first.
#
# NOTE: this assumes getResourceCapacity() exists in your installed mizer
# version, mirroring the getResourceRate() already used above. If it errors,
# swap it for the accessor resource_capacity(params) instead -- same value,
# different name across mizer versions.
################################################################################

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

# Same resource_decrease grid as the resource_level sweep above (4 points),
# capacity swept 10x below to 10x above the default, log-spaced -- 4 x 5 = 20
# combinations, same total budget as that sweep.
rd_grid_kr    <- exp(seq(log(0.0001), log(0.5), length.out = 4))
cc_grid       <- exp(seq(log(0.1), log(10), length.out = 5))
param_grid_kr <- expand.grid(resource_decrease = rd_grid_kr, capacity_mult = cc_grid)

plan(multisession)
osc_results_kr <- future_lapply(seq_len(nrow(param_grid_kr)), function(i) {
  score_oscillation_kr(param_grid_kr$resource_decrease[i], param_grid_kr$capacity_mult[i])
}, future.seed = TRUE)
osc_df_kr <- bind_rows(osc_results_kr)
plan(sequential)

# Same collapsed-population guard as the resource_level sweep: rel_amplitude
# is meaningless once mean_bm has crashed near double-precision zero, so
# those combinations get greyed out / filtered rather than winning the
# ranking on floating-point noise.
osc_df_kr <- osc_df_kr %>%
  mutate(collapsed            = mean_bm < MIN_VIABLE_BIOMASS,
         rel_amplitude_viable = ifelse(collapsed, NA_real_, rel_amplitude))

ggplot(osc_df_kr, aes(x = resource_decrease, y = capacity_mult, fill = rel_amplitude_viable)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(name = "Relative\namplitude", na.value = "grey85") +
  labs(x = "resource_decrease (resource rate multiplier)",
       y = "capacity_mult (carrying capacity multiplier)",
       title = "Searching for the oscillatory regime: carrying capacity x resource rate",
       subtitle = "Brighter = bigger (max - min) / mean; grey = population collapsed (mean_bm < 1e-2)") +
  theme_minimal()

osc_df_kr %>%
  filter(!collapsed) %>%
  arrange(desc(rel_amplitude)) %>%
  head(10)

################################################################################
# Reintroducing oscillations, take 3: bracketing the enrichment threshold
#
# Take 2 found real (non-collapsed, non-numerical-noise) oscillations, but
# only at capacity_mult = 10 -- the top edge of that grid -- across three
# very different resource_decrease values (0.0017, 0.029, 0.5), while
# capacity_mult = 3.16 gave clean fixed points (max_bm == min_bm to machine
# precision) at the same resource_decrease values. That's a paradox-of-
# enrichment signature: raising carrying capacity past some threshold
# destabilises the fixed point into a limit cycle, the opposite direction
# from the original "lower the carrying capacity" hypothesis. Since the
# strongest signal sat at the edge of what was tested, this brackets the
# 3.16-10 transition more finely and pushes further out to 100x to see
# whether amplitude keeps growing or the system restabilises again further
# into enrichment. resource_decrease didn't seem to matter much for whether
# oscillation appeared, so it gets fewer points here in favour of resolving
# capacity_mult, keeping the same ~20-run budget as takes 1 and 2.
################################################################################

rd_grid_kr2 <- exp(seq(log(0.0001), log(0.5), length.out = 4))
cc_grid2    <- exp(seq(log(3), log(100), length.out = 5))
param_grid_kr2 <- expand.grid(resource_decrease = rd_grid_kr2, capacity_mult = cc_grid2)

plan(multisession)
osc_results_kr2 <- future_lapply(seq_len(nrow(param_grid_kr2)), function(i) {
  score_oscillation_kr(param_grid_kr2$resource_decrease[i], param_grid_kr2$capacity_mult[i])
}, future.seed = TRUE)
osc_df_kr2 <- bind_rows(osc_results_kr2)
plan(sequential)

osc_df_kr2 <- osc_df_kr2 %>%
  mutate(collapsed            = mean_bm < MIN_VIABLE_BIOMASS,
         rel_amplitude_viable = ifelse(collapsed, NA_real_, rel_amplitude))

ggplot(osc_df_kr2, aes(x = resource_decrease, y = capacity_mult, fill = rel_amplitude_viable)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(name = "Relative\namplitude", na.value = "grey85") +
  labs(x = "resource_decrease (resource rate multiplier)",
       y = "capacity_mult (carrying capacity multiplier)",
       title = "Bracketing the enrichment threshold: capacity_mult 3-100",
       subtitle = "Brighter = bigger (max - min) / mean; grey = population collapsed (mean_bm < 1e-2)") +
  theme_minimal()

osc_df_kr2 %>%
  filter(!collapsed) %>%
  arrange(desc(rel_amplitude)) %>%
  head(10)

################################################################################
# Sanity check: does rel_amplitude actually mean "oscillating"?
#
# rel_amplitude is a single settled-window (max - min) / mean number -- it
# can't distinguish a genuine limit cycle from, say, a slow monotonic drift
# that hasn't fully settled by t = 600, or a one-off transient spike. This
# plots the raw biomass trace over t = 550-600 for a handful of combinations,
# the same late-window overlay Day 18 used to actually confirm the limit
# cycle was a real attractor rather than trusting a summary statistic at
# face value -- a genuine cycle should show a clean repeating waveform here,
# not just a nonzero number in the table.
################################################################################

plot_late_biomass <- function(combos, t_run = 600, t_window = c(550, 600)) {
  traces <- lapply(seq_len(nrow(combos)), function(i) {
    p   <- make_second_order_params_kr(resource_decrease = combos$resource_decrease[i],
                                       capacity_mult     = combos$capacity_mult[i])
    sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                   progress_bar = FALSE, effort = 0, method = "tr_bdf2")
    bm  <- getBiomass(sim)[, "Anchovy"]
    tv  <- as.numeric(names(bm))
    data.frame(time    = tv, biomass = as.numeric(bm),
               combo   = sprintf("rd = %.4g, cap_mult = %.3g",
                                 combos$resource_decrease[i], combos$capacity_mult[i]))
  })
  traces_df <- bind_rows(traces) %>%
    filter(time >= t_window[1], time <= t_window[2])

  ggplot(traces_df, aes(x = time, y = biomass, color = combo)) +
    geom_line(linewidth = 0.8) +
    labs(x = "Time", y = "Biomass",
         title = paste0("Biomass, t = ", t_window[1], "-", t_window[2]),
         subtitle = "Flat line = fixed point; repeating waveform = real limit cycle",
         color = "Combination") +
    theme_minimal()
}

# The three real (non-collapsed) hits reported from take 2, plus one
# capacity_mult = 3.162 case for contrast -- that one scored rel_amplitude
# ~5e-15 (a fixed point to machine precision), so it should show up here as a
# flat line while the other three should show a repeating waveform if the
# rel_amplitude numbers are trustworthy.
check_combos <- data.frame(
  resource_decrease = c(0.029240177, 0.001709976, 0.500000000, 0.029240177),
  capacity_mult      = c(10,          10,          10,          3.162278)
)
plot_late_biomass(check_combos)

# To check the take-3 (3-100x) sweep's own top candidates instead, once
# osc_df_kr2 has been computed:
# check_combos_kr2 <- osc_df_kr2 %>%
#   filter(!collapsed) %>%
#   arrange(desc(rel_amplitude)) %>%
#   head(4) %>%
#   select(resource_decrease, capacity_mult)
# plot_late_biomass(check_combos_kr2)

################################################################################
# Size spectrum for resource_decrease = 0.5, capacity_mult = 10 -- one of take
# 2's three real hits (max_bm = 10.33, min_bm = 6.70, mean_bm = 8.54,
# rel_amplitude = 0.42). Uses mizer's own plotSpectra() rather than a custom
# plot.
################################################################################

p_rd05_cap10   <- make_second_order_params_kr(resource_decrease = 0.5, capacity_mult = 10)
sim_rd05_cap10 <- project(p_rd05_cap10, t_max = 600, dt = 0.1, t_save = 0.5,
                          progress_bar = FALSE, effort = 0, method = "tr_bdf2")

plotSpectra(sim_rd05_cap10, power = 2)
animate(sim_rd05_cap10,tlim=c(550,600))

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
