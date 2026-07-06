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

make_second_order_params <- function(lambda = 2.05, resource_decrease = 0.001,
                                     second_order = TRUE, ext_diff = 0.01) {
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
  
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)
  
  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  params <- setBevertonHolt(params)
  
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
  
  params
}

################################################################################

################################################################################

# live_plot = TRUE passes mizerExperimental's biomass_callback to project(),
# same as run_rd_sweep() above -- see the note there. Open x11() (or leave
# the RStudio plot pane visible) before calling this with live_plot = TRUE so
# you can watch each of the 20+20 steps in the sweep settle in real time.
run_rd_sweep_second_order <- function(rd_seq, init_n = NULL, init_n_pp = NULL,
                                      t_run = 600, lambda = 2.05,
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
       title = "Bifurcation diagram: second-order scheme + ext_diffusion ",
       subtitle = "Branches collapsing onto one curve = fixed point; fanning apart = limit cycle") +
  theme_minimal()






################################################################################
# Time series: how biomass actually fluctuates over time, rather than just the
# settled max/min summary the sweep above collapses each run down to.
# Picks a few resource_decrease values spanning the hysteresis region from the
# sweep (deep in the oscillating regime, through the transition, out past
# where the backward branch collapsed) and plots getBiomass(sim) across the
# whole run for each, on one time axis, so the shrinking oscillation
# amplitude/eventual flat line is visible directly rather than inferred from
# max - min.
################################################################################

rd_time_series <- c(0.001, 0.01, 0.05, 0.15, 0.3)

biomass_over_time <- lapply(rd_time_series, function(rd) {
  p   <- make_second_order_params(resource_decrease = rd)
  sim <- project(p, t_max = 300, dt = 0.1, t_save = 0.5,
                progress_bar = FALSE, effort = 0, method = "tr_bdf2")
  bm  <- getBiomass(sim)[, "Anchovy"]
  data.frame(time = as.numeric(names(bm)), biomass = as.numeric(bm),
            resource_decrease = factor(rd))
})
biomass_over_time_df <- bind_rows(biomass_over_time)

ggplot(biomass_over_time_df, aes(x = time, y = biomass, color = resource_decrease)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Time", y = "Biomass",
       title = "Biomass over time for a few resource_decrease values",
       subtitle = "Flat line = settled fixed point; oscillating = limit cycle",
       color = "resource_decrease") +
  theme_minimal()

################################################################################
# Same rd_time_series values, but now run as a proper forward/backward pair --
# each step carries the previous step's settled n/n_pp forward as its initial
# condition (like fwd_so/bwd_so above), rather than every rd starting fresh
# the way biomass_over_time above does. Backward starts from the forward run's
# final state, same convention as bwd_so <- ...(init_n = fwd_so$n_final, ...).
# Time is offset cumulatively across steps so each direction reads as one
# continuous trajectory, and the two get overlaid on the same time axis to
# compare the dynamics of approaching each resource_decrease from below vs
# from above -- not just the settled max/min, which is all the rd sweep above
# captures.
################################################################################

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

make_second_order_params <- function(lambda = 2.05, resource_decrease = 0.001,
                                     second_order = TRUE, ext_diff = 0.01) {
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
  
  given_species_params(params)$D_ext <- ext_diff
  params
}

################################################################################

################################################################################

# live_plot = TRUE passes mizerExperimental's biomass_callback to project(),
# same as run_rd_sweep() above -- see the note there. Open x11() (or leave
# the RStudio plot pane visible) before calling this with live_plot = TRUE so
# you can watch each of the 20+20 steps in the sweep settle in real time.
run_rd_sweep_second_order <- function(rd_seq, init_n = NULL, init_n_pp = NULL,
                                      t_run = 600, lambda = 2.05,
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

bwd_so <- run_rd_sweep_second_order(rev(rd_seq_so),
                                    label     = "BWD (2nd order)", live_plot = TRUE)
fwd_so <- run_rd_sweep_second_order(rd_seq_so,,
                                    init_n    = bwd_so$n_final,
                                    init_n_pp = bwd_so$npp_final, label = "FWD (2nd order)", live_plot = TRUE)


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
       title = "Bifurcation diagram: second-order scheme + ext_diffusion ",
       subtitle = "Branches collapsing onto one curve = fixed point; fanning apart = limit cycle") +
  theme_minimal()


p <- make_second_order_params(lambda = lambda, resource_decrease = 0.5,
                              second_order = second_order)
sim <- project(p,t_max = 300, dt = 0.1, t_save = 0.5,
               progress_bar = FALSE, effort = 0, method = "tr_bdf2")
plotHover(getBiomass(sim))
