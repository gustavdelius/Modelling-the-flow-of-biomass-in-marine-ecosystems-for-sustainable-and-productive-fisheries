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
r <- r*0.001
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

# --- Branch 1: no fishing (this is your existing pipeline) ---
sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                   effort = 0, progress_bar = FALSE, method = "predictor-corrector")

plotHover(getBiomass(sim_600,power=2),tlim=c(550,600))
nice_biomass_plot(sim_600,550)
#nice_animation(sim_600,550)
plotHover(getFlux(sim_600,power=2),log="xy",tlim=c(550,600))
#----Branch 2: fishing always ------
#----Branch 2: fishing always, sub-adults only ------
params_fished <- sim_300@params

gp <- params_fished@gear_params
gp$sel_func        <- "box_selectivity_func"
gp$w_low           <- 0.5 * p$w_mat
gp$w_high          <- p$w_mat
gp$catchability    <- 1
gp$knife_edge_size <- NULL
gear_params(params_fished) <- gp     # <-- triggers setFishing() and rebuilds @selectivity

sim_300_fished        <- sim_300
sim_300_fished@params <- params_fished

# ===== juvenile-inclusive comparison variant =====
params_fished_juv <- sim_300@params

gp_juv <- params_fished_juv@gear_params
gp_juv$sel_func        <- "box_selectivity_func"
gp_juv$w_low           <- 0.01
gp_juv$w_high          <- p$w_mat
gp_juv$catchability    <- 1
gp_juv$knife_edge_size <- NULL
gear_params(params_fished_juv) <- gp_juv

sim_300_fished_juv        <- sim_300
sim_300_fished_juv@params <- params_fished_juv

sim_fished <- project(sim_300_fished, t_max = 300, dt = 0.1, t_save = 0.2,
                      effort = 0.3, progress_bar = FALSE,
                      method = "predictor-corrector")
plotRelative(getFlux(sim_600,power=2),getFlux(sim_fished,power=2),tlim=c(550,600))
plotRelative(getBiomass(sim_600),getBiomass(sim_fished),tlim=c(550,600))
nice_biomass_plot(sim_fished,550)
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

n <- length(bm_p)
is_peak <- logical(n)
is_peak[2:(n-1)] <- bm_p[2:(n-1)] > bm_p[1:(n-2)] & bm_p[2:(n-1)] > bm_p[3:n]
peak_times <- t_p[is_peak]

# 2. Build a 0/fish_level effort schedule: "on" only within `window`
#    years either side of each peak, "off" everywhere else (troughs included)
window     <- 3    # years either side of peak to fish
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

plotRelative(getFlux(sim_600,power=2),getFlux(sim_peaks_only,power=2),tlim=c(550,600))
plotRelative(getFlux(sim_fished,power=2),getFlux(sim_peaks_only,power=2),tlim=c(550,600))

plotYield(sim_peaks_only,sim_fished,tlim=c(550,600))

plotYield(sim_peaks_only, sim_fished, tlim = c(550, 600))
y_fished <- getYield(sim_peaks_only)
y_base   <- getYield(sim_fished)
t_f      <- as.numeric(rownames(y_fished))
t_b      <- as.numeric(rownames(y_base))

mean(y_fished[t_f > 550, ])   # average yield, fishing the bottleneck
mean(y_base[t_b > 550, ])     # average yield, no fishing (will be 0)


# ================== Confirm: does excluding juveniles explain the improvement? ==================
# Build a "juvenile-inclusive" gear variant (your original setup) to compare
# against the "sub-adults only" variant in params_fished -- now with the
# selectivity fix applied, so this is a real comparison, not a no-op.

params_fished_juv <- sim_300@params
params_fished_juv@gear_params$sel_func        <- "box_selectivity_func"
params_fished_juv@gear_params$w_low           <- 0.01   # juveniles included
params_fished_juv@gear_params$w_high          <- p$w_mat
params_fished_juv@gear_params$catchability    <- 1
params_fished_juv@gear_params$knife_edge_size <- NULL

sim_300_fished_juv        <- sim_300
sim_300_fished_juv@params <- params_fished_juv

# --- Constant effort, juveniles included ---
sim_fished_juv <- project(sim_300_fished_juv, t_max = 300, dt = 0.1, t_save = 0.2,
                          effort = 0.3, progress_bar = FALSE,
                          method = "predictor-corrector")

# --- Peaks-only, juveniles included (same effort_arr schedule, different gear) ---
sim_peaks_only_juv <- project(sim_300_fished_juv, t_max = 300, dt = 0.1, t_save = 0.2,
                              effort = effort_arr, progress_bar = FALSE,
                              method = "predictor-corrector")

# --- Compare juvenile-inclusive vs sub-adult-only, within each fishing strategy ---

# Constant effort: does excluding juveniles change the bottleneck?
plotRelative(getFlux(sim_fished_juv, power = 2), getFlux(sim_fished, power = 2),
             tlim = c(550, 600))
plotRelative(getBiomass(sim_fished_juv), getBiomass(sim_fished), tlim = c(550, 600))

# Peaks-only: does excluding juveniles change the bottleneck?
plotRelative(getFlux(sim_peaks_only_juv, power = 2), getFlux(sim_peaks_only, power = 2),
             tlim = c(550, 600))
plotRelative(getBiomass(sim_peaks_only_juv), getBiomass(sim_peaks_only), tlim = c(550, 600))

# Biomass ranges, sanity check that the gear change actually did something this time
range(getBiomass(sim_fished)[, "Anchovy"])
range(getBiomass(sim_peaks_only)[, "Anchovy"])
range(getBiomass(sim_fished_juv)[, "Anchovy"])
range(getBiomass(sim_peaks_only_juv)[, "Anchovy"])

plotBiomass(sim_fished,         species = "Anchovy", start_time = 550, end_time = 600)
plotBiomass(sim_peaks_only,     species = "Anchovy", start_time = 550, end_time = 600)
plotBiomass(sim_fished_juv,     species = "Anchovy", start_time = 550, end_time = 600)
plotBiomass(sim_peaks_only_juv, species = "Anchovy", start_time = 550, end_time = 600)

# Yield trade-off: are you gaining bottleneck relief at the cost of yield?
y_fished_juv     <- getYield(sim_fished_juv)
y_peaks_only_juv <- getYield(sim_peaks_only_juv)
t_fj <- as.numeric(rownames(y_fished_juv))
t_pj <- as.numeric(rownames(y_peaks_only_juv))

mean(y_fished_juv[t_fj > 550, "Anchovy"])       # constant, juveniles included
mean(y_base[t_b > 550, "Anchovy"])              # constant, sub-adults only
mean(y_peaks_only_juv[t_pj > 550, "Anchovy"])   # peaks-only, juveniles included
mean(y_fished[t_f > 550, "Anchovy"])            # peaks-only, sub-adults only
