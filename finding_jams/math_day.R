library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)
library(glue)
library(ggplot2)
#Useful Plots
nice_animation <- function(sim,t){
  nf <- melt(sim@n)
  n_ppf <- melt(sim@n_pp)
  n_ppf$sp <- "Plankton"
  nf <- rbind(nf, n_ppf)
  
  plot_ly(nf) %>%
    # show only part of plankton spectrum
    filter(w > 10^-5) %>%
    # start at time 50
    filter(time >= t) %>%
    # calculate biomass density with respect to log size
    mutate(b = value * w^2) %>%
    # Plot lines
    add_lines(
      x = ~w, y = ~b,
      color = ~sp,
      frame = ~time,
      line = list(simplify = FALSE)
    ) %>%
    # Use logarithmic axes
    layout(xaxis = list(type = "log", exponentformat = "power",
                        title_text = "body mass (g)"),
           yaxis = list(type = "log", exponentformat = "power",
                        title_text = "biomass (g/m^3)",
                        range = c(-8, 0)))
}

nice_biomass_plot <- function(sim,t){
  params <- getParams(sim,c(550,600))
  abm <- melt(getBiomass(sim))
  abmr <- melt(getBiomass(sim, min_w = 0.01, max_w = 0.4))
  abmr$sp <- "small Anchovy"
  pbm <- sim@n_pp %*% (params@w_full * params@dw_full)
  pbm <- melt(pbm)
  names(pbm)[names(pbm) == "Var1"] <- "time"
  pbm$Var2 <- NULL
  pbm$sp <- "Plankton"
  bm <- rbind(pbm, abm, abmr)
  plot_ly(bm) %>%
    filter(time >= t) %>%
    add_lines(x = ~time, y = ~value, color = ~sp) %>%
    # Use logarithmic axes
    layout(yaxis = list(type = "log", exponentformat = "power",
                        title_text = "biomass (g/m^3)",
                        range = c(-7, 2)),
           xaxis = list(title_text = "time (year)"))
}

################## Code from yesterday that I'll experiment on #################

modelling_steady_state <- function(resource_decrease=1, f_effort=0.3){
  p <- list(
    dt     = 0.001,
    dx     = 0.1,
    w_min  = 0.0003,
    w_max  = 66.5,     # was w_inf
    w_mat  = 10,
    alpha  = 0.1,      # assimilation efficiency (was p$K)
    gamma  = 750,
    lambda = 2.05,
    a0     = 100
  )
  
  kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))
  no_w  <- round(log(p$w_max / p$w_min) / p$dx)
  
  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min  = p$w_min,
    w_max  = p$w_max,
    w_mat  = p$w_mat,
    no_w   = no_w,
    lambda = p$lambda,
    kappa  = kappa,
    alpha  = p$alpha,
    gamma  = p$gamma,
    ks     = 0
  )
  r <- getResourceRate(params)
  r <- r*resource_decrease
  params <- setResource(params, resource_rate = r, resource_dynamics = "resource_semichemostat")
  
  params@initial_n_pp[] <- params@cc_pp * 0.1
  
  test_sizes <- list(
    "small"  = c(0.1, 1),
    "medium" = c(1, 10),
    "large"  = c(10, 100)
  )
  
  # --- Shared phase: perturbation + settling, no fishing ---
  sim <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                 effort = 0, progress_bar = FALSE, method = "predictor-corrector")
  
  # rng <- test_sizes[["large"]]
  # idx <- params@w >= rng[1] & params@w <= rng[2]
  # last <- dim(sim@n)[1]
  # sim@n[last, , idx] <- sim@n[last, , idx] / 10^3
  
  sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
                 effort = 0, progress_bar = FALSE, method = "predictor-corrector")
  
  sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                     effort = 0, progress_bar = FALSE, method = "predictor-corrector")

  # --- Constant fishing from t = 300 onwards: realistic gear (knife-edge at
  # maturity, no juveniles caught), with a constant effort thereafter ---
  params_fished <- sim_300@params

  gp <- params_fished@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- p$w_mat   # minimum landing size = maturity weight, nothing below w_mat is caught
  gp$catchability    <- 1
  gp$w_low           <- NULL
  gp$w_high          <- NULL
  gear_params(params_fished) <- gp     # <-- triggers setFishing() and rebuilds @selectivity

  sim_300_fished        <- sim_300
  sim_300_fished@params <- params_fished

  sim_600 <- project(sim_300_fished, t_max = 300, dt = 0.1, t_save = 0.2,
                     effort = f_effort, progress_bar = FALSE, method = "predictor-corrector")

}

# Custom FMort rate function for "fish only after a biomass peak".
# Mizer resolves a rate function set via setRateFunction() by looking up its
# name with get() from the global environment (or a package) -- a closure
# defined inside modelling_fishing_at_peaks() would never be found there, so
# this factory and the name it's assigned to must live at top level. Each
# call to modelling_fishing_at_peaks() builds a fresh closure and reassigns
# it under the same name, so separate runs don't share peak-tracking state.

make_dynamic_peak_FMort <- function(window) {
  state <- new.env()
  state$t_prev     <- NA_real_   # time of the last distinct recorded point
  state$bm_prev1   <- NA_real_   # biomass one recorded step back
  state$bm_prev2   <- NA_real_   # biomass two recorded steps back
  state$fish_until <- -Inf       # fishing stays on until this time

  function(params, n, n_pp, n_other, t, effort, e_growth, pred_mort, ...) {
    bm <- sum(n[1, ] * params@w * params@dw)

    # predictor-corrector calls FMort twice per dt at the same t; only
    # advance the rolling history once per distinct time point
    if (is.na(state$t_prev) || t > state$t_prev) {
      if (!is.na(state$bm_prev2) &&
          state$bm_prev1 > state$bm_prev2 && state$bm_prev1 > bm) {
        # bm_prev1, recorded at state$t_prev, was a local peak in the
        # actual (fished) trajectory -> open/extend the fishing window
        state$fish_until <- max(state$fish_until, state$t_prev + window)
      }
      state$bm_prev2 <- state$bm_prev1
      state$bm_prev1 <- bm
      state$t_prev   <- t
    }

    dynamic_effort <- effort
    dynamic_effort[] <- if (t <= state$fish_until) effort else 0

    mizerFMort(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
               effort = dynamic_effort, e_growth = e_growth, pred_mort = pred_mort, ...)
  }
}

modelling_fishing_at_peaks <- function(resource_decrease=0.001,f_effort=0.3,peak_window=1){
  p <- list(
    dt     = 0.001,
    dx     = 0.1,
    w_min  = 0.0003,
    w_max  = 66.5,     # was w_inf
    w_mat  = 10,
    alpha  = 0.1,      # assimilation efficiency (was p$K)
    gamma  = 750,
    lambda = 2.05,
    a0     = 100
  )
  
  kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))
  no_w  <- round(log(p$w_max / p$w_min) / p$dx)
  
  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min  = p$w_min,
    w_max  = p$w_max,
    w_mat  = p$w_mat,
    no_w   = no_w,
    lambda = p$lambda,
    kappa  = kappa,
    alpha  = p$alpha,
    gamma  = p$gamma,
    ks     = 0
  )
  r <- getResourceRate(params)
  r <- r*resource_decrease
  params <- setResource(params, resource_rate = r, resource_dynamics = "resource_semichemostat")
  
  params@initial_n_pp[] <- params@cc_pp * 0.1
  
  test_sizes <- list(
    "small"  = c(0.1, 1),
    "medium" = c(1, 10),
    "large"  = c(10, 100)
  )
  
  # --- Shared phase: perturbation + settling, no fishing ---
  sim <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                 effort = 0, progress_bar = FALSE, method = "predictor-corrector")
  
  sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
                 effort = 0, progress_bar = FALSE, method = "predictor-corrector")
  
  sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                     effort = 0, progress_bar = FALSE, method = "predictor-corrector")
  
  # --- Branch 1: no fishing (this is your existing pipeline) ---
  sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                     effort = 0, progress_bar = FALSE, method = "predictor-corrector")
  
  #----Branch 2: fishing always, realistic gear (knife-edge at maturity, no juveniles caught) ------
  params_fished <- sim_300@params

  gp <- params_fished@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- p$w_mat   # minimum landing size = maturity weight, nothing below w_mat is caught
  gp$catchability    <- 1
  gp$w_low           <- NULL
  gp$w_high          <- NULL
  gear_params(params_fished) <- gp     # <-- triggers setFishing() and rebuilds @selectivity
  
  # --- Branch 3: fish only at peaks, never at troughs (dynamic, online) ---
  # Peaks are detected live, from the fished trajectory's own biomass as it
  # is simulated, via the custom FMort rate function above -- not pre-computed
  # once from the unfished sim_600. Only effort is gated; selectivity (set on
  # params_fished above) is never touched here.
  window     <- peak_window
  fish_level <- f_effort

  assign("dynamic_peak_FMort", make_dynamic_peak_FMort(window = window), envir = .GlobalEnv)

  params_fished_peaks <- setRateFunction(params_fished, "FMort", "dynamic_peak_FMort")

  sim_300_fished_peaks        <- sim_300
  sim_300_fished_peaks@params <- params_fished_peaks

  sim_peaks_only <- project(sim_300_fished_peaks, t_max = 300, dt = 0.1, t_save = 0.2,
                            effort = fish_level, progress_bar = FALSE,
                            method = "predictor-corrector")
}
sim_const_fish <- modelling_steady_state(0.001)
sim_peaks      <- modelling_fishing_at_peaks()

y_fished <- getYield(sim_peaks)
y_base   <- getYield(sim_const_fish)
t_f      <- as.numeric(rownames(y_fished))
t_b      <- as.numeric(rownames(y_base))

mean(y_fished[t_f > 550, ])   # average yield, fishing the bottleneck (peaks-only)
mean(y_base[t_b > 550, ])     # average yield, fishing constantly (sub-adults only)

plotYield(sim_const_fish,sim_peaks,log="x",tlim=c(550,600))

nice_biomass_plot(sim_peaks,550)
nice_biomass_plot(sim_const_fish,550)
plotRelative(getBiomass(sim_const_fish),getBiomass(sim_peaks),tlim=c(550,600))

peak_yield_summary <- function(f_effort, peak_window) {
  sim <- modelling_fishing_at_peaks(f_effort = f_effort, peak_window = peak_window)
  yt  <- getYield(sim)
  t   <- as.numeric(rownames(yt))
  mean(yt[t > 300, ])
}

effort_grid <- seq(0.1, 1, length.out = 5)
window_grid <- seq(0.5, 5, length.out = 5)

grid <- expand.grid(f_effort = effort_grid, peak_window = window_grid)
grid$mean_yield <- mapply(peak_yield_summary, grid$f_effort, grid$peak_window)

best <- grid[which.max(grid$mean_yield), ]
print(grid[order(-grid$mean_yield), ])
print(best)

ggplot(grid, aes(x = f_effort, y = peak_window, fill = mean_yield)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Fishing effort", y = "Peak window (years)",
       fill = "Mean yield", title = "Yield surface, peaks-only fishing")

effort_grid <- c(1, 1.5, 2, 2.5, 3, 4, 5, 7, 10)
window_grid <- c(0.1, 0.25, 0.5, 0.75, 1)

grid2 <- expand.grid(f_effort = effort_grid, peak_window = window_grid)
grid2$mean_yield <- mapply(peak_yield_summary, grid2$f_effort, grid2$peak_window)

print(grid2[order(-grid2$mean_yield), ][1:15, ])
best2 <- grid2[which.max(grid2$mean_yield), ]
print(best2)

ggplot(grid2, aes(x = f_effort, y = peak_window, fill = mean_yield)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(x = "Fishing effort", y = "Peak window (years)",
       fill = "Mean yield", title = "Refined yield surface, peaks-only fishing")

effort_grid <- seq(0.1, 2, length.out = 8)
window_grid <- seq(0.5, 3, length.out = 6)

grid3 <- expand.grid(f_effort = effort_grid, peak_window = window_grid)
grid3$mean_yield <- mapply(peak_yield_summary, grid3$f_effort, grid3$peak_window)

print(grid3[order(-grid3$mean_yield), ][1:15, ])
best3 <- grid3[which.max(grid3$mean_yield), ]
print(best3)

ggplot(grid3, aes(x = f_effort, y = peak_window, fill = mean_yield)) +
  geom_tile() + scale_fill_viridis_c() +
  labs(x = "Fishing effort", y = "Peak window (years)",
       fill = "Mean yield", title = "Bounded yield surface, peaks-only fishing")
