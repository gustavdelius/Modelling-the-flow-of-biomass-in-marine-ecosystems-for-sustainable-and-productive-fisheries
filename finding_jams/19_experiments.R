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
# Helper Functions
# Builds the standard single-species Anchovy params.
# beta/sigma default to NULL so, when left unspecified, newSingleSpeciesParams()
# falls back to its own built-in defaults (i.e. behaviour is unchanged from
# before beta/sigma were added here).
make_params <- function(lambda = 2.05, resource_decrease = 0.005,
                        beta = NULL, sigma = NULL, balance = FALSE) {
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

# Runs the two-stage perturbation needed to kick the system onto the limit cycle.
# Stage 1 (t=0–10): run with depleted resource to destabilise.
# Stage 2 (t=10): knock down mature fish by 1000×, then run to t_total.
make_limit_cycle_sim <- function(params, t_total = 600, effort = 0,perturbation=1e3) {
  params@initial_n_pp[] <- params@cc_pp * 0.1
  
  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / perturbation
  
  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "tr_bdf2")
}

#Verification: setResource()'s docs say that with balance = TRUE (the default),
# supplying only resource_rate makes mizer silently back-solve resource_capacity
# (cc_pp) so that "the resource replenishes at the same rate at which it is
# consumed". This checks whether cc_pp actually differs between balance = TRUE
# and balance = FALSE for the same resource_decrease, i.e. whether balance = TRUE
# is depleting the resource ceiling as well as its replenishment rate.

check_balance_effect <- function(resource_decrease, lambda = 2.05) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  base <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  r <- getResourceRate(base) * resource_decrease

  p_balanced   <- setResource(base, resource_rate = r,
                              resource_dynamics = "resource_semichemostat", balance = TRUE)
  p_unbalanced <- setResource(base, resource_rate = r,
                              resource_dynamics = "resource_semichemostat", balance = FALSE)

  data.frame(
    resource_decrease    = resource_decrease,
    cc_pp_default_max    = max(base@cc_pp),
    cc_pp_balanced_max   = max(p_balanced@cc_pp),
    cc_pp_unbalanced_max = max(p_unbalanced@cc_pp)
  )
}

balance_check <- bind_rows(lapply(c(0.001, 0.01, 0.1, 0.3), check_balance_effect))
balance_check
# If cc_pp_balanced_max << cc_pp_default_max (and shrinks further as
# resource_decrease shrinks) while cc_pp_unbalanced_max stays at the default,
# that confirms balance = TRUE is compounding the depletion by lowering the
# resource ceiling, not just slowing its replenishment.

#Experiment: Testing resource decrease to learn more about what leads to oscillations when balance=False


test_effort <- function(r_d, lambda = 2.05, t_total = 600) {
  p   <- make_params(lambda = lambda,resource_decrease = r_d)   # resource_decrease stays at 0.001
  sim <- make_limit_cycle_sim(p, t_total = t_total)
  bm  <- getBiomass(sim)[, "Anchovy"]
  t   <- as.numeric(names(bm))
  list(sim = sim, bm = bm, t = t)
}

r_d <- exp(seq(log(0.001), log(0.3), length.out = 10))  # log-spaced, matches Day 18's oscillatory window

plan(multisession)  # spawns background R workers; uses all available cores by default
results_fishing <- future_lapply(r_d, test_effort, future.seed = TRUE)

# Trajectory panels
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (i in seq_along(results_fishing)) {  # fixed variable name
  r <- results_fishing[[i]]
  plot(r$t, r$bm, type = "l",
       xlab = "Time (years)", ylab = "Biomass",
       main = paste0("r_d = ", r_d[i]))
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Filter plot for t > 550
plot(NULL, xlim = c(550, 600), ylim = range(sapply(results_fishing, function(r) range(r$bm[r$t > 550]))),
     xlab = "Time (years)", ylab = "Biomass", main = "Biomass trajectories (t > 550)")
cols <- colorRampPalette(c("blue", "red"))(length(r_d))
for (i in seq_along(results_fishing)) {
  r    <- results_fishing[[i]]
  mask <- r$t > 550
  lines(r$t[mask], r$bm[mask], col = cols[i])
}
legend("topright", legend = paste0("r_d = ", r_d), col = cols, lty = 1, cex = 0.8)

#Control experiment: the cc_pp check on a fresh, never-projected params object
# came back a no-op (balance = TRUE and balance = FALSE gave identical cc_pp),
# which means it doesn't actually explain the flat-vs-oscillating difference
# seen so far. This isolates *only* balance, holding every other setting
# (solvers, perturbation, r_d grid) fixed, to check whether balance is really
# the causal variable or whether something else (e.g. the predictor-corrector/
# tr_bdf2 solver mismatch in make_limit_cycle_sim) is responsible instead.

test_effort <- function(r_d, lambda = 2.05, t_total = 600, balance = FALSE) {
  p   <- make_params(lambda = lambda, resource_decrease = r_d, balance = balance)
  sim <- make_limit_cycle_sim(p, t_total = t_total)
  bm  <- getBiomass(sim)[, "Anchovy"]
  t   <- as.numeric(names(bm))
  list(sim = sim, bm = bm, t = t)
}

plan(multisession)
results_balance_true  <- future_lapply(r_d, test_effort, balance = TRUE,  future.seed = TRUE)
results_balance_false <- future_lapply(r_d, test_effort, balance = FALSE, future.seed = TRUE)

amp_by_balance <- data.frame(
  r_d          = r_d,
  amp_balance_true  = sapply(results_balance_true,  function(r) {
    late <- r$bm[r$t > 550]; max(late) - min(late)
  }),
  amp_balance_false = sapply(results_balance_false, function(r) {
    late <- r$bm[r$t > 550]; max(late) - min(late)
  })
)
amp_by_balance
# If amp_balance_true is much larger than amp_balance_false at the same r_d,
# balance is confirmed causal. If both columns are equally flat, balance
# isn't the cause and the solver mismatch (or something else) is.

# amp_balance_false climbs smoothly across ~16 orders of magnitude with r_d,
# which looks like numerical noise scaling with biomass magnitude rather than
# a flat "no oscillation" floor. Check the actual biomass level (not just
# amplitude) to see whether balance = FALSE is driving the population toward
# near-extinction at small r_d rather than a healthy non-oscillating steady
# state.
bm_final_balance <- data.frame(
  r_d           = r_d,
  bm_final_true  = sapply(results_balance_true,  function(r) tail(r$bm, 1)),
  bm_final_false = sapply(results_balance_false, function(r) tail(r$bm, 1))
)
bm_final_balance
# If bm_final_false is near zero at small r_d and grows toward bm_final_true's
# scale as r_d increases, that confirms extinction rather than a stable
# non-oscillating equilibrium.
#
# CONFIRMED: at r_d = 0.001, bm_final_false ~ 1.6e-17 (extinct) vs
# bm_final_true ~ 0.025 (real biomass). Even at r_d = 0.3, the top of the
# tested range, bm_final_false ~ 1.1e-3 is still ~76x below bm_final_true's
# stable-fixed-point plateau (~0.087). So balance = FALSE chronically
# underfeeds the population across the whole grid, and drives it toward
# extinction at small r_d, instead of reaching either a healthy stable fixed
# point or the oscillatory regime. This is why the beta/sigma sweep above
# (run at r_d = 0.001, i.e. deep in the near-extinct regime) found nothing:
# no feeding-kernel shape matters when there is essentially no population left
# to feed.

#Experiment: Testing whether varying the feeding kernel (beta, sigma) at a
# fixed, relatively low resource_decrease can bring the oscillations back
# under balance = FALSE. beta is the preferred predator:prey mass ratio,
# sigma is the width of the (log-normal) selectivity around that ratio.

test_beta_sigma <- function(beta, sigma, r_d = 0.001, lambda = 2.05, t_total = 600) {
  p   <- make_params(lambda = lambda, resource_decrease = r_d, beta = beta, sigma = sigma)
  sim <- make_limit_cycle_sim(p, t_total = t_total)
  bm  <- getBiomass(sim)[, "Anchovy"]
  t   <- as.numeric(names(bm))
  list(sim = sim, bm = bm, t = t, beta = beta, sigma = sigma)
}

beta_seq  <- c(25,50, 100, 200, 400,1000)
sigma_seq <- c(0.1,0.3,0.5,0.8, 1.3, 1.8, 2.3)
grid      <- expand.grid(beta = beta_seq, sigma = sigma_seq)

plan(multisession)
results_beta_sigma <- future_lapply(seq_len(nrow(grid)), function(i) {
  test_beta_sigma(beta = grid$beta[i], sigma = grid$sigma[i], r_d = 0.001)
}, future.seed = TRUE)

# Trajectory panels, one per (beta, sigma) combination
par(mfrow = c(4, 4), mar = c(2, 2, 2, 1))
for (i in seq_along(results_beta_sigma)) {
  r <- results_beta_sigma[[i]]
  plot(r$t, r$bm, type = "l", xlab = "", ylab = "",
       main = paste0("beta=", r$beta, ", sigma=", r$sigma), cex.main = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Amplitude heatmap: does any (beta, sigma) combination restore oscillation
# amplitude at this low, previously-oscillatory resource_decrease?
amp_df <- data.frame(
  beta  = grid$beta,
  sigma = grid$sigma,
  amp   = sapply(results_beta_sigma, function(r) {
    late <- r$bm[r$t > 550]
    max(late) - min(late)
  })
)

ggplot(amp_df, aes(x = factor(beta), y = factor(sigma), fill = amp)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "beta", y = "sigma",
       fill = "Amplitude\n(max - min, t > 550)",
       title = "Oscillation amplitude across beta/sigma grid (balance = FALSE, r_d = 0.01)") +
  theme_minimal()

#Experiment: the beta/sigma sweep above failed because it was run at a
# resource_rate/capacity combination that was already starving the
# population (see bm_final_balance above). This time, resource_rate and
# resource_capacity are pulled directly from the run that produced the
# strongest oscillation (r_d = 0.001, balance = TRUE, amplitude ~0.159 in
# amp_by_balance), fixed explicitly, and beta/sigma are varied on top of
# that known-good baseline instead.

best_params   <- make_params(resource_decrease = r_d[1], balance = TRUE)
best_rate     <- best_params@rr_pp
best_capacity <- best_params@cc_pp

# Supplying both resource_rate and resource_capacity leaves nothing for
# balance to back-solve, so balance = FALSE here just sets them directly.
make_params_fixed_resource <- function(lambda = 2.05, resource_rate, resource_capacity,
                                       beta = NULL, sigma = NULL) {
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
  setResource(params, resource_rate = resource_rate, resource_capacity = resource_capacity,
             resource_dynamics = "resource_semichemostat", balance = FALSE)
}

test_beta_sigma_fixed <- function(beta, sigma, resource_rate, resource_capacity,
                                  lambda = 2.05, t_total = 600) {
  p   <- make_params_fixed_resource(lambda = lambda, resource_rate = resource_rate,
                                    resource_capacity = resource_capacity,
                                    beta = beta, sigma = sigma)
  sim <- make_limit_cycle_sim(p, t_total = t_total)
  bm  <- getBiomass(sim)[, "Anchovy"]
  t   <- as.numeric(names(bm))
  list(sim = sim, bm = bm, t = t, beta = beta, sigma = sigma)
}

# newSingleSpeciesParams() defaults to beta = 100, sigma = 1.3, so the grid
# is centred on those defaults with room to move in both directions.
beta_seq_fixed  <- c(50, 75, 100, 150, 200, 300)
sigma_seq_fixed <- c(0.8, 1.0, 1.3, 1.6, 2.0)
grid_fixed      <- expand.grid(beta = beta_seq_fixed, sigma = sigma_seq_fixed)

plan(multisession)
results_beta_sigma_fixed <- future_lapply(seq_len(nrow(grid_fixed)), function(i) {
  test_beta_sigma_fixed(beta = grid_fixed$beta[i], sigma = grid_fixed$sigma[i],
                        resource_rate = best_rate, resource_capacity = best_capacity)
}, future.seed = TRUE)

amp_df_fixed <- data.frame(
  beta  = grid_fixed$beta,
  sigma = grid_fixed$sigma,
  amp   = sapply(results_beta_sigma_fixed, function(r) {
    late <- r$bm[r$t > 550]
    max(late) - min(late)
  })
)

ggplot(amp_df_fixed, aes(x = factor(beta), y = factor(sigma), fill = amp)) +
  geom_tile() +
  geom_point(data = subset(amp_df_fixed, beta == 100 & sigma == 1.3),
             aes(x = factor(beta), y = factor(sigma)),
             shape = 4, size = 4, color = "red") +
  scale_fill_viridis_c() +
  labs(x = "beta", y = "sigma",
       fill = "Amplitude\n(max - min, t > 550)",
       title = "Oscillation amplitude across beta/sigma grid",
       subtitle = "resource_rate/capacity fixed at the best-oscillating balance = TRUE values; red x = defaults (beta=100, sigma=1.3)") +
  theme_minimal()
# Cells brighter than the red-marked default cell identify a feeding-kernel
# shape that sustains a larger-amplitude oscillation than the default kernel,
# at the same resource baseline that is already known to support oscillation.
amp_df_fixed

#Follow-up 1: look at the actual biomass trajectories, not just the amplitude
# summary, for every (beta, sigma) combination in the grid.
par(mfrow = c(5, 6), mar = c(2, 2, 2, 1))
for (i in seq_along(results_beta_sigma_fixed)) {
  r    <- results_beta_sigma_fixed[[i]]
  mask <- r$t > 500
  plot(r$t[mask], r$bm[mask], type = "l", xlab = "", ylab = "",
       main = paste0(r$beta, ", ", r$sigma), cex.main = 0.7)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

#Follow-up 2: overlay the highest-amplitude cell, the default, and the
# lowest-amplitude cell directly, to compare waveform shape and period, not
# just amplitude, across the grid.
compare_spec <- list(
  list(beta = 50,  sigma = 0.8, label = "beta=50, sigma=0.8 (highest amp)",  col = "firebrick"),
  list(beta = 100, sigma = 1.3, label = "beta=100, sigma=1.3 (default)",     col = "black"),
  list(beta = 300, sigma = 2.0, label = "beta=300, sigma=2.0 (lowest amp)",  col = "steelblue")
)
compare_idx <- sapply(compare_spec, function(s) {
  which(grid_fixed$beta == s$beta & grid_fixed$sigma == s$sigma)
})
compare_runs <- results_beta_sigma_fixed[compare_idx]

plot(NULL, xlim = c(550, 600),
     ylim = range(sapply(compare_runs, function(r) range(r$bm[r$t > 550]))),
     xlab = "Time (years)", ylab = "Biomass",
     main = "Biomass trajectories: highest-amp vs default vs lowest-amp (t > 550)")
for (j in seq_along(compare_runs)) {
  r    <- compare_runs[[j]]
  mask <- r$t > 550
  lines(r$t[mask], r$bm[mask], col = compare_spec[[j]]$col, lwd = 2)
}
legend("topright", legend = sapply(compare_spec, `[[`, "label"),
       col = sapply(compare_spec, `[[`, "col"), lwd = 2, cex = 0.8)

#Follow-up 3: check biomass LEVEL (not just amplitude) across the grid, to
# tell apart a genuine stable fixed point (healthy nonzero biomass, no
# oscillation) from starvation (near-zero biomass) for the low-amplitude
# cells -- the same distinction that mattered for bm_final_balance earlier.
bm_level_fixed <- data.frame(
  beta         = grid_fixed$beta,
  sigma        = grid_fixed$sigma,
  bm_mean_late = sapply(results_beta_sigma_fixed, function(r) mean(r$bm[r$t > 550])),
  bm_final     = sapply(results_beta_sigma_fixed, function(r) tail(r$bm, 1))
)

ggplot(bm_level_fixed, aes(x = factor(beta), y = factor(sigma), fill = bm_mean_late)) +
  geom_tile() +
  geom_point(data = subset(bm_level_fixed, beta == 100 & sigma == 1.3),
             aes(x = factor(beta), y = factor(sigma)),
             shape = 4, size = 4, color = "red") +
  scale_fill_viridis_c() +
  labs(x = "beta", y = "sigma", fill = "Mean biomass\n(t > 550)",
       title = "Mean late-time biomass across beta/sigma grid",
       subtitle = "Low here + low amplitude = starvation; healthy biomass + low amplitude = genuine stable fixed point") +
  theme_minimal()
# If the low-amplitude cells in the top-right of the amplitude heatmap have
# bm_mean_late comparable to bm_final_true's ~0.087 plateau (not collapsed
# toward zero), that confirms they're real stable fixed points, not another
# starvation case.

