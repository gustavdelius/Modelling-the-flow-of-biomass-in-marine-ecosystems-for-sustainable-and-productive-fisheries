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
  lambda = 2,
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

# Baseline: same Anchovy params, no perturbation, same dt/t_save/t_max as sim_600
sim_base <- project(params, t_max = 600, dt = 0.1, t_save = 0.2,
                    progress_bar = FALSE, method = "predictor-corrector")

r <- resource_rate(params)
r <- r * 0.001
params <- setResource(params, resource_rate = r, resource_dynamics = "resource_semichemostat")

params@initial_n_pp[] <- params@cc_pp * 0.1



test_sizes <- list(
  "small"  = c(0.1, 1),
  "medium" = c(1, 10),
  "large"  = c(10, 100)
)

# Run up to T1
sim <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

# Select the size range to perturb
rng <- test_sizes[["large"]]
idx <- params@w >= rng[1] & params@w <= rng[2]

# Divide that size range, at the current end of sim, by 10^3
last <- dim(sim@n)[1]
sim@n[last, , idx] <- sim@n[last, , idx] / 10^3

# Continue for the rest of the run
sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

plotHover(getBiomass(sim_600), tlim = c(550, 600))
plotHover(getBiomass(sim_base), tlim = c(550, 600))
plotYield(sim_600, sim_base, tlim = c(550, 600))