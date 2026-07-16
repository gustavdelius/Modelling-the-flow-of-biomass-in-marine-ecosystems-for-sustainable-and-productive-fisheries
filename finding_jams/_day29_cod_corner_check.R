suppressPackageStartupMessages({library(mizer); library(mizerExperimental); library(dplyr)})

cod_params <- readRDS("cod_params.rds")

resource_limitation_2d_cod <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}
build_cod <- function(resource_decrease, capacity_mult) resource_limitation_2d_cod(cod_params, resource_decrease, capacity_mult)

make_limit_cycle_sim <- function(params, t_total = 600, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1
  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0, method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 1e3
  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort, method = "predictor-corrector")
}

rel_amplitude_at <- function(sim, window_frac = 0.4) {
  bm    <- rowSums(getBiomass(sim))
  times <- as.numeric(names(bm))
  late  <- bm[times >= max(times) * (1 - window_frac)]
  (max(late) - min(late)) / mean(late)
}

# The hottest cod corner point from the full-grid run.
p <- build_cod(resource_decrease = 0.00187382, capacity_mult = 1)

t_max_seq <- c(600, 1200, 2400, 4800, 9600)
cat("Rerunning cod's hottest corner point (resource_decrease=0.00187, capacity_mult=1) at increasing t_max:\n")
for (tm in t_max_seq) {
  sim <- make_limit_cycle_sim(p, t_total = tm)
  ra  <- rel_amplitude_at(sim)
  cat(sprintf("  t_max=%5d: rel_amplitude = %.6g\n", tm, ra))
}

# Also check the endpoint distance to steady state directly, like Day 26's
# check -- pull the actual biomass steady state via steadyNewton() for
# comparison against where the long run ends up.
cat("\nFor reference, steadyNewton()'s own verdict at this point:\n")
st <- steadyNewton(p, stability = TRUE)
stab <- attr(st, "stability")
cat(sprintf("  spectral_radius = %.6g, stable = %s\n", stab$spectral_radius, stab$stable))
