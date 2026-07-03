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
library(future.apply)

# Day 20 picks up two threads left open by Day 18/19:
#   Part 1-2: the forward/backward resource_decrease sweep left a gap at
#     rd ~ 0.008 that couldn't be written off as a sweep-start detection
#     artefact -- plot max and min separately (not collapsed to
#     max - min) to actually see the shape of that gap, then try to pin down
#     whether it's real hysteresis or a run-length/resolution artefact.
#   Part 3: recommendations (diagnostic and experimental) for investigating
#     why so much biomass piles up at juvenile sizes, continuing the
#     unfinished survivorship-through-growth thread at the bottom of
#     19_experiments.R.

# Helper functions, carried over from Day 18/19. balance defaults to TRUE here
# (unlike 19_experiments.R's FALSE) because every experiment below reuses
# Day 18's original oscillating regime, where balance = TRUE.
make_params <- function(lambda = 2.05, resource_decrease = 0.001,
                        beta = NULL, sigma = NULL, balance = TRUE) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  args <- list(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  if (!is.null(beta))  args$beta  <- beta
  if (!is.null(sigma)) args$sigma <- sigma

  params <- do.call(newSingleSpeciesParams, args)
  r      <- getResourceRate(params) * resource_decrease
  setResource(params, resource_rate = r,
             resource_dynamics = "resource_semichemostat", balance = balance)
}

make_limit_cycle_sim <- function(params, t_total = 600, effort = 0, perturbation = 1e3) {
  params@initial_n_pp[] <- params@cc_pp * 0.1

  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "tr_bdf2")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / perturbation

  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "tr_bdf2")
}

################################################################################
# Part 1: Max and min plotted separately, forward vs backward, instead of
# collapsing to a single amp = max - min number. Day 18's amp-only plot could
# not tell apart "the whole cycle sits differently depending on direction"
# (real hysteresis) from "one branch drifts while the other tracks tightly"
# (an asymmetric artefact, e.g. critical slowing down dragging the mean level
# around). run_rd_sweep is extended to record max_bm/min_bm, not just amp.
################################################################################

# live_plot = TRUE passes mizerExperimental's biomass_callback to project(),
# which redraws a live biomass plot into the currently active graphics
# device at every saved timestep (see
# https://sizespectrum.org/mizerExperimental/reference/biomass_callback.html)
# -- open a device that stays open across all the sequential project() calls
# in this loop first (x11() in RStudio, as the callback's own doc example
# recommends; the default RStudio plot pane also works but redraws slower),
# otherwise each call opens/reuses whatever the current device is and you
# only see the last one settle.
#
# biomass_callback is a stateful closure (it caches which species to plot
# the first time it's called, in its own enclosing environment, then reuses
# that cache on every later call) -- it works out which species to draw from
# params@linecolour by default, and for a hand-built single-species params
# object like the one from make_params(), that inference can come back
# empty. When it does, the callback just silently no-ops on every timestep:
# the device opens, but nothing ever gets drawn onto it. Passing species =
# explicitly (forwarded through project()'s ... into the callback call) is
# what actually fixes that, rather than relying on the default inference.
run_rd_sweep <- function(rd_seq, init_n = NULL, init_n_pp = NULL,
                         t_run = 300, lambda = 2.05, label = "",
                         live_plot = FALSE) {
  out       <- data.frame(rd = rd_seq, max_bm = NA_real_, min_bm = NA_real_,
                          bm_mean = NA_real_)
  state_n   <- init_n
  state_npp <- init_n_pp

  for (i in seq_along(rd_seq)) {
    p <- make_params(lambda = lambda, resource_decrease = rd_seq[i])
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
    late <- bm[tv > t_run * 0.6]
    last <- dim(sim@n)[1]

    state_n   <- sim@n[last, , ]
    state_npp <- sim@n_pp[last, ]

    out$max_bm[i]  <- max(late)
    out$min_bm[i]  <- min(late)
    out$bm_mean[i] <- mean(late)
  }
  list(df = out, n_final = state_n, npp_final = state_npp)
}

# Log-spaced rather than linear: spans 0.0001-0.5 (3.7 orders of magnitude,
# vs the previous 0.001-0.02) while keeping roughly the same per-decade
# resolution as before (40 points / 3.7 decades =~ 11 points per decade)
# instead of thinning out badly at the low end the way a linear seq() would
# over this much wider a range.
rd_seq <- exp(seq(log(0.0001), log(0.5), length.out = 40))

fwd <- run_rd_sweep(rd_seq, label = "FWD")
bwd <- run_rd_sweep(rev(rd_seq),
                    init_n    = fwd$n_final,
                    init_n_pp = fwd$npp_final,
                    label     = "BWD")

fwd_df <- fwd$df
bwd_df <- bwd$df[order(bwd$df$rd), ]

# Base-R overlay: four curves (fwd max, fwd min, bwd max, bwd min) on the same
# axes. Below the transition, max and min should collapse onto each other
# (the stable fixed point); above it they fan out into two branches (the
# cycle's peak and trough). Where forward and backward branches fail to
# overlap is the real hysteresis region -- not just a difference in a single
# collapsed "amp" number.
# log = "x" because rd_seq is now log-spaced over 3.7 orders of magnitude
# (0.0001-0.5) -- on a linear x-axis the low end would bunch up into an
# unreadable cluster of points.
plot(fwd_df$rd, fwd_df$max_bm, type = "l", col = "steelblue", lwd = 2, lty = 1,
     log = "x", xlab = "resource_decrease", ylab = "Biomass",
     main = "Hysteresis: max and min plotted separately",
     ylim = range(c(fwd_df$max_bm, fwd_df$min_bm, bwd_df$max_bm, bwd_df$min_bm)))
lines(fwd_df$rd, fwd_df$min_bm, col = "steelblue", lwd = 2, lty = 2)
lines(bwd_df$rd, bwd_df$max_bm, col = "firebrick", lwd = 2, lty = 1)
lines(bwd_df$rd, bwd_df$min_bm, col = "firebrick", lwd = 2, lty = 2)
legend("topright",
       legend = c("Forward max", "Forward min", "Backward max", "Backward min"),
       col = c("steelblue", "steelblue", "firebrick", "firebrick"),
       lty = c(1, 2, 1, 2), lwd = 2, cex = 0.8)

# Same data, tidier bifurcation-diagram-style ggplot version -- long format,
# one line per (direction, branch) combination.
bifurcation_df <- bind_rows(
  data.frame(rd = fwd_df$rd, biomass = fwd_df$max_bm, direction = "Forward", branch = "max"),
  data.frame(rd = fwd_df$rd, biomass = fwd_df$min_bm, direction = "Forward", branch = "min"),
  data.frame(rd = bwd_df$rd, biomass = bwd_df$max_bm, direction = "Backward", branch = "max"),
  data.frame(rd = bwd_df$rd, biomass = bwd_df$min_bm, direction = "Backward", branch = "min")
)

ggplot(bifurcation_df, aes(x = rd, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_log10() +
  labs(x = "resource_decrease", y = "Biomass",
       title = "Bifurcation diagram: max/min branches, forward vs backward",
       subtitle = "Branches collapsing onto one curve = fixed point; fanning apart = limit cycle") +
  theme_minimal()

################################################################################
# Part 2: Is the gap real hysteresis, or a run-length/resolution artefact?
# Two follow-ups, both flagged in Day 18's revised What's Next.
################################################################################

# 2a. Finer rd grid, zoomed on the region where the backward curve was still
# oscillating (rd ~ 0.008) while the forward curve looked already flat. If the
# gap between forward and backward transition points shrinks to near-nothing
# at this resolution, the original 10-point sweep was just undersampling a
# sharp but continuous transition. If a clear, stable gap survives, that's
# much harder to explain away as sweep resolution.
rd_fine <- seq(0.005, 0.012, length.out = 15)

fwd_fine <- run_rd_sweep(rd_fine, t_run = 300, label = "FWD fine")
bwd_fine <- run_rd_sweep(rev(rd_fine),
                         init_n    = fwd_fine$n_final,
                         init_n_pp = fwd_fine$npp_final,
                         t_run     = 300, label = "BWD fine")

fwd_fine_df <- fwd_fine$df
bwd_fine_df <- bwd_fine$df[order(bwd_fine$df$rd), ]

plot(fwd_fine_df$rd, fwd_fine_df$max_bm - fwd_fine_df$min_bm, type = "b",
     col = "steelblue", lwd = 2, xlab = "resource_decrease",
     ylab = "Amplitude (max - min)",
     main = "Fine-grid hysteresis check (rd = 0.005-0.012)",
     ylim = c(0, max(fwd_fine_df$max_bm - fwd_fine_df$min_bm,
                     bwd_fine_df$max_bm - bwd_fine_df$min_bm) * 1.1))
lines(bwd_fine_df$rd, bwd_fine_df$max_bm - bwd_fine_df$min_bm, type = "b",
      col = "firebrick", lwd = 2)
legend("topright", legend = c("Forward", "Backward"),
       col = c("steelblue", "firebrick"), lwd = 2)

################################################################################
# Part 3: mizer 3.1's second-order finite-volume scheme (second_order_w()),
# plus a fixed background external diffusion (setExtDiffusion()), as a params
# constructor to reuse in Part 4 below.
#
# setExtDiffusion()'s ext_diffusion argument has to be a species x size
# array, not a bare scalar ("ext_diffusion is not an array" is that shape
# check failing) -- getDiffusion(params) is already a validly-shaped species
# x size array at that point (all zero, since no ext diffusion has been set
# yet), so it's used as a template for dim/dimnames and filled with the
# scalar ext_diff.
#
# NOTE: second_order_w() is new/inferred from mizer 3.1's announcement
# (https://blog.mizer.sizespectrum.org/posts/2026-06-26-mizer-3-1-announcement/)
# rather than tested against this project's installed mizer version --
# smoke-test make_second_order_params() on a short run before trusting Part 4.
################################################################################

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
  r      <- getResourceRate(params) * resource_decrease
  params <- setResource(params, resource_rate = r,
                        resource_dynamics = "resource_semichemostat",
                        balance = TRUE)
  
  # Mizer 3.1's second-order finite-volume scheme -- see NOTE above. Both
  # sub-options (bin_average + centred flux) enabled together for the fully
  # second-order variant, rather than just one of the two independently.
  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }
  
  # setExtDiffusion() needs ext_diffusion as an array (species x size) --
  # "ext_diffusion is not an array" is that shape check failing on a bare
  # scalar. getDiffusion(params) at this point is already a validly-shaped
  # species x size array (all zero, since no ext diffusion is set yet), so
  # it's used as a template for dim/dimnames and filled with the scalar.
  
  
  diffusion_template  <- getDiffusion(params)
  ext_diffusion_array <- array(ext_diff, dim = dim(diffusion_template),
                               dimnames = dimnames(diffusion_template))
  
  params <- setExtDiffusion(params, ext_diffusion = ext_diffusion_array)
  params
}
lambda <- 2.05
a0    <- 100
kappa <- a0 * exp(-6.9 * (lambda - 1))
no_w  <- round(log(66.5 / 0.0003) / 0.1)

params <- newSingleSpeciesParams(
  species_name = "Anchovy",
  w_min = 0.0003, w_max = 66.5, w_mat = 10,
  no_w = no_w, lambda = lambda, kappa = kappa,
  alpha = 0.1, gamma = 750, ks = 0
)
plotSpectra(params,power=2,log="y",resource=FALSE)
plot(getBiomass(params))
################################################################################
# Part 4: the same forward/backward, max/min-separated bifurcation plot as
# Part 1, rebuilt on make_second_order_params() instead of make_params(). The
# sweep mechanics are unchanged from Part 1: multiple project() calls in a
# row, each one carrying the previous run's final state (n, n_pp) forward as
# the next run's initial condition, with resource_decrease stepped a bit
# lower each time going forward and a bit higher each time going backward --
# still discrete, sequential projections, NOT a continuous in-simulation
# ramp. second_order and ext_diff are held fixed across the whole sweep; only
# resource_decrease varies step to step, same as Part 1.
################################################################################

# live_plot = TRUE passes mizerExperimental's biomass_callback to project(),
# same as run_rd_sweep() above -- see the note there. Open x11() (or leave
# the RStudio plot pane visible) before calling this with live_plot = TRUE so
# you can watch each of the 20+20 steps in the sweep settle in real time.
run_rd_sweep_second_order <- function(rd_seq, init_n = NULL, init_n_pp = NULL,
                                      t_run = 300, lambda = 2.05,
                                      second_order = TRUE, ext_diff = 0, label = "",
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
                                  second_order = second_order, ext_diff = ext_diff)
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

rd_seq_so <- rd_seq  # same log-spaced grid as Part 1 (0.0001-0.5, 40 points) -- keeps the two directly comparable

# x11() is Unix/X11 only -- it errors out (or silently does nothing useful)
# on native Windows R, which is almost certainly why the live plot didn't
# show anything before. windows() is the Windows equivalent graphics device;
# quartz() is macOS's. This opens a persistent external window so
# biomass_callback's live plot has somewhere to keep drawing into across all
# 40 project() calls below, rather than each one fighting over whatever
# device happens to be active. Comment this out (and set live_plot = FALSE
# below) if running non-interactively, e.g. Rscript on a machine with no
# display.
if (.Platform$OS.type == "windows") {
  windows()
} else if (Sys.info()[["sysname"]] == "Darwin") {
  quartz()
} else {
  x11()
}
fwd_so <- run_rd_sweep_second_order(rd_seq_so, label = "FWD (2nd order)", live_plot = TRUE)
bwd_so <- run_rd_sweep_second_order(rev(rd_seq_so),
                                    init_n    = fwd_so$n_final,
                                    init_n_pp = fwd_so$npp_final,
                                    label     = "BWD (2nd order)", live_plot = TRUE)

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
       title = "Bifurcation diagram: second-order scheme + ext_diffusion = 0.001",
       subtitle = "Branches collapsing onto one curve = fixed point; fanning apart = limit cycle") +
  theme_minimal()
# Compare directly against Part 1's bifurcation_df/plot (plain make_params,
# no second-order scheme, no ext_diffusion): a narrower gap between the
# forward and backward branches here would point at the Part 1 gap being at
# least partly a numerical artefact of the default first-order scheme and/or
# zero ext_diffusion; a similar or wider gap argues the opposite.
