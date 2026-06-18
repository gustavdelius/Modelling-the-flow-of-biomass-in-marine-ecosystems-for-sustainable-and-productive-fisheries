library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)
library(glue)
library(ggplot2)

################## Code from yesterday that I'll experiment on #################

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

r <- resource_rate(params)
r <- r * 0.001
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

rng <- test_sizes[["large"]]
idx <- params@w >= rng[1] & params@w <= rng[2]
last <- dim(sim@n)[1]
sim@n[last, , idx] <- sim@n[last, , idx] / 10^3

sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               effort = 0, progress_bar = FALSE, method = "predictor-corrector")

sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                   effort = 0, progress_bar = FALSE, method = "predictor-corrector")

# --- Branch 1: no fishing (this is your existing pipeline) ---
sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                   effort = 0, progress_bar = FALSE, method = "predictor-corrector")

# --- Branch 2: fishing targeted at the bottleneck size range ---
# Bottleneck = juveniles below maturation (w_mat = 10g), where the
# Day 8 phase analysis put the food-limitation queue (amplitude peak ~0.01-0.05g).

bottleneck_sel <- function(w, w_low, w_high) {
  as.numeric(w >= w_low & w <= w_high)
}

params_fished <- sim_300@params
params_fished@gear_params$sel_func     <- "bottleneck_sel"
params_fished@gear_params$w_low        <- 0.01
params_fished@gear_params$w_high       <- p$w_mat
params_fished@gear_params$catchability <- 1

sim_300_fished        <- sim_300
sim_300_fished@params <- params_fished

sim_fished <- project(sim_300_fished, t_max = 300, dt = 0.1, t_save = 0.2,
                      effort = 0.3, progress_bar = FALSE,
                      method = "predictor-corrector")

# --- Compare yield ---
plotYield(sim_fished, sim_600, tlim = c(550, 600))

y_fished <- getYield(sim_fished)
y_base   <- getYield(sim_600)
t_f      <- as.numeric(rownames(y_fished))
t_b      <- as.numeric(rownames(y_base))

mean(y_fished[t_f > 550, ])   # average yield, fishing the bottleneck
mean(y_base[t_b > 550, ])     # average yield, no fishing (will be 0)


# --- Branch 3: fish only at peaks, never at troughs ---

# 1. Find peak times from the settled, unfished trajectory (sim_600),
#    reusing sim_600's own dimnames directly to avoid floating-point
#    drift from re-parsing/re-formatting times, and including the
#    fork point t = 300 itself (not just times strictly greater than it)
all_times <- as.numeric(dimnames(sim_600@n)[[1]])
keep      <- all_times >= 300
t_p       <- all_times[keep]
time_lbls <- dimnames(sim_600@n)[[1]][keep]

bm   <- getBiomass(sim_600)[, "Anchovy"]
bm_p <- bm[time_lbls]

is_peak    <- c(FALSE, diff(sign(diff(bm_p))) == -2, FALSE)
peak_times <- t_p[is_peak]

# 2. Build a 0/fish_level effort schedule: "on" only within `window`
#    years either side of each peak, "off" everywhere else (troughs included)
window     <- 1     # years either side of peak to fish
fish_level <- 0.3

effort_vec <- rep(0, length(t_p))
for (pt in peak_times) {
  effort_vec[t_p >= pt - window / 2 & t_p <= pt + window / 2] <- fish_level
}

gear_name  <- params_fished@gear_params$gear[1]
effort_arr <- matrix(effort_vec, ncol = 1,
                     dimnames = list(time_lbls, gear_name))

# 3. Run the bottleneck-gear branch with this schedule instead of constant effort
sim_peaks_only <- project(sim_300_fished, t_max = 300, dt = 0.1, t_save = 0.2,
                          effort = effort_arr, progress_bar = FALSE,
                          method = "predictor-corrector")

# Sanity check before plotting — should be TRUE
identical(dimnames(sim_peaks_only@n)[[1]], dimnames(sim_600@n)[[1]])

plotYield(sim_peaks_only, sim_600, tlim = c(550, 600))
y_fished <- getYield(sim_fished)
y_base   <- getYield(sim_600)
t_f      <- as.numeric(rownames(y_fished))
t_b      <- as.numeric(rownames(y_base))

mean(y_fished[t_f > 550, ])   # average yield, fishing the bottleneck
mean(y_base[t_b > 550, ])     # average yield, no fishing (will be 0)
plotHover(getBiomass(sim_600),tlim=c(550,600))
plotHover(getBiomass(sim_peaks_only),tlim=c(550,600))
plotHover(getFlux(sim_600,power=2),tlim=c(550,600))
plotHover(getFlux(sim_peaks_only,power=2),tlim=c(550,600),wlim=c(5e-3,5e-2))
