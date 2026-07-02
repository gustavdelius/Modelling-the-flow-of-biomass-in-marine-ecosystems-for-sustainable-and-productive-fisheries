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

# Helper Functions
# Builds the standard single-species Anchovy params.
make_params <- function(lambda = 2.05, resource_decrease = 0.001) {
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
                        resource_dynamics = "resource_semichemostat")
  params
}

# Runs the two-stage perturbation needed to kick the system onto the limit cycle.
# Stage 1 (t=0–10): run with depleted resource to destabilise.
# Stage 2 (t=10): knock down mature fish by 1000×, then run to t_total.
make_limit_cycle_sim <- function(params, t_total = 600, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1
  
  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 1e3
  
  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "predictor-corrector")
}


############# Experiment 1: Finding where the species goes extinct ############

test_constant_fishing <- function(effort, lambda = 2.05, t_total = 600) {
  p  <- make_params(lambda = lambda)
  gp <- p@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- p@species_params$w_mat
  gp$catchability    <- 1
  gear_params(p)     <- gp
  
  sim <- make_limit_cycle_sim(p, t_total = t_total, effort = effort)
  bm  <- getBiomass(sim)[, "Anchovy"]
  t   <- as.numeric(names(bm))
  list(sim = sim, bm = bm, t = t)
}

efforts_fishing <- c(6,7, 8, 9, 10, 15,20)
results_fishing <- lapply(efforts_fishing, test_constant_fishing)

# Fixed point reference for lambda = 2.05
p_ref    <- make_params(lambda = 2.05)
ss_ref   <- steadyNewton(p_ref)
bm_fp_ref <- sum(getBiomass(ss_ref))

# Trajectory panels
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (i in seq_along(efforts_fishing)) {
  r       <- results_fishing[[i]]
  bm_mean <- mean(r$bm[r$t > 400])
  plot(r$t, r$bm, type = "l",
       xlab = "Time (years)", ylab = "Biomass",
       main = paste0("Effort = ", efforts_fishing[i]))
  abline(h = bm_fp_ref, col = "firebrick", lty = 2)
  abline(h = bm_mean,   col = "steelblue", lty = 2)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Summary table: mean biomass, peak, trough, yield
cat("\nExp 7 — Constant fishing summary (lambda = 2.05):\n")
for (i in seq_along(efforts_fishing)) {
  r       <- results_fishing[[i]]
  settled <- r$bm[r$t > 400]
  y       <- getYield(r$sim)[, "Anchovy"]
  t_y     <- as.numeric(rownames(getYield(r$sim)))
  cat(sprintf("  Effort %.2f:  mean = %.5f  max = %.5f  min = %.5f  yield = %.6f\n",
              efforts_fishing[i],
              mean(settled), max(settled), min(settled),
              mean(y[t_y > 400])))
}

# Experiment 2: Seeing whether resource alone can cause oscillations

#make_params <- function(resource_decrease = 0.001) {
  params <- newSingleSpeciesParams(
    species_name = "Anchovy"
  )
  r      <- getResourceRate(params) * resource_decrease
  params <- setResource(params, resource_rate = r,
                        resource_dynamics = "resource_semichemostat")
  params
}

# Runs the two-stage perturbation needed to kick the system onto the limit cycle.
# Stage 1 (t=0–10): run with depleted resource to destabilise.
# Stage 2 (t=10): knock down mature fish by 1000×, then run to t_total.
#make_limit_cycle_sim <- function(params, t_total = 600, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1
  
  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 2
  
  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "predictor-corrector")
}
params_only_resource_different <- make_params()
lim_cycle_only_resource <- make_limit_cycle_sim(params_only_resource_different)
plotHover(getBiomass(lim_cycle_only_resource),tlim=c(550,600))
#Makes sense - no perturbation, no getting away from the steady state solution.


# Experiment 3: Building the model up

# =============================================================================
# Helpers
# =============================================================================

make_params <- function(lambda = 2.05, resource_decrease = 0.001) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)
  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  r <- getResourceRate(params) * resource_decrease
  setResource(params, resource_rate = r,
              resource_dynamics = "resource_semichemostat")
}

# steady() → projectToSteady() → project_simple() has a bug in mizer 3.1.0.9000:
# r$rdd is missing from the simplified rates list, crashing project_n_no_diffusion.
# project() uses the full rates path and works; this wrapper extracts the final state.
get_steady <- function(params, t_max = 300, dt = 0.1) {
  sim  <- project(params, t_max = t_max, dt = dt, t_save = t_max,
                  method = "predictor-corrector", progress_bar = FALSE)
  last <- dim(sim@n)[1]
  n_new               <- params@initial_n
  n_new[]             <- sim@n[last, , ]
  params@initial_n    <- n_new
  params@initial_n_pp <- sim@n_pp[last, ]
  params
}

add_bump <- function(params_ss, w_c = 0.01, amplitude = 1.0, sigma_log = 0.5) {
  w          <- params_ss@w
  bump       <- amplitude * exp(-(log(w / w_c))^2 / sigma_log^2)
  n_new      <- params_ss@initial_n
  n_new[1, ] <- n_new[1, ] * (1 + bump)
  params_ss@initial_n <- n_new
  params_ss
}



# =============================================================================
# Baseline
# =============================================================================
p_base    <- make_params(resource_decrease = 1)
ss_base   <- get_steady(p_base)          # project() for 300 yr; take final state
n_ss_base <- ss_base@initial_n[1, ]

# RDI: the density-independent recruitment rate at steady state (eggs produced
# before any density-dependent competition is applied)
rdi_ss    <- getRDI(ss_base, n = ss_base@initial_n,
                    n_pp = ss_base@initial_n_pp,
                    n_other = ss_base@initial_n_other)["Anchovy"]

# R_max: the half-saturation constant in the Beverton-Holt function;
# when RDI = R_max, exactly half the maximum recruitment survives competition
r_max_std <- species_params(ss_base)$R_max

# RDD: the actual (density-dependent) recruitment at steady state, computed via
# Beverton-Holt: RDD = R_max * RDI / (R_max + RDI).
# Setting a model's R_max to this value forces BH to always return rdd_ss
# regardless of RDI, making recruitment effectively constant (Level 1).
rdd_ss    <- r_max_std * rdi_ss / (r_max_std + rdi_ss)

# =============================================================================
# Level 0: static steady-state spectrum
# =============================================================================
plot(p_base@w, n_ss_base * p_base@w, log = "xy", type = "l", lwd = 2,
     col = "steelblue", xlab = "Body mass w (g)",
     ylab = "Biomass density B(w)",
     main = "Level 0: steady-state biomass spectrum")
abline(v = 10, col = "darkgreen", lty = 3)

# =============================================================================
# Level 1: pure wave transport  (R_max = rdd_ss → BH saturated → constant recruitment)
# =============================================================================
p1 <- ss_base
species_params(p1)$R_max <- rdd_ss

sim1 <- project(add_bump(p1), t_max = 20, dt = 0.01, t_save = 0.1,
                method = "predictor-corrector", progress_bar = FALSE)

animateSpectra(sim1,power=2,log="xy",tlim=c(0,20))

# =============================================================================
# Level 2: reproduction feedback  (standard R_max restored → BH active)
# =============================================================================
p2 <- ss_base
species_params(p2)$R_max <- r_max_std

sim2 <- project(add_bump(p2), t_max = 60, dt = 0.01, t_save = 0.1,
                method = "predictor-corrector", progress_bar = FALSE)

animateSpectra(sim2,power=2,log="xy",tlim=c(0,20))

# =============================================================================
# Level 3: resource depletion → limit cycle
# (steadyNewton uses a Newton solver, not project_simple, so it works)
# =============================================================================
p3 <- make_params(resource_decrease = 0.001)
species_params(p3)$R_max <- r_max_std

ss3   <- steadyNewton(p3)
n_ss3 <- ss3@initial_n[1, ]

sim3 <- project(add_bump(ss3), t_max = 200, dt = 0.01, t_save = 0.1,
                method = "predictor-corrector", progress_bar = FALSE)

animateSpectra(sim3,power=2,log="xy",tlim=c(150,200),resource=FALSE)

