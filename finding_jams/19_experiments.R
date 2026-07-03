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

#Follow-up 4: the trajectory panels showed several "flat-looking" cells
# (75,0.8 plus most of beta = 150/200/300) whose y-axis range is absurdly
# tiny (~1e-6), which is the signature of a run that hasn't finished
# converging, not a small real oscillation. Plot the FULL trajectory from
# t = 0 (not just the t > 500 tail) for a few of these alongside a
# genuinely oscillating cell, to see whether they're damped spirals
# (decaying oscillation settling to a fixed point -- a real dynamical
# signature) or were never oscillating at all (pure monotonic relaxation).
full_trace_spec <- list(
  list(beta = 50,  sigma = 0.8, label = "beta=50, sigma=0.8 (oscillating)"),
  list(beta = 75,  sigma = 0.8, label = "beta=75, sigma=0.8 (the odd one out)"),
  list(beta = 150, sigma = 0.8, label = "beta=150, sigma=0.8 (looks flat)"),
  list(beta = 300, sigma = 1.3, label = "beta=300, sigma=1.3 (looks flat)")
)
full_trace_idx <- sapply(full_trace_spec, function(s) {
  which(grid_fixed$beta == s$beta & grid_fixed$sigma == s$sigma)
})

par(mfrow = c(2, 2), mar = c(3, 3, 2, 1))
for (j in seq_along(full_trace_idx)) {
  r <- results_beta_sigma_fixed[[full_trace_idx[j]]]
  plot(r$t, r$bm, type = "l", xlab = "Time (years)", ylab = "Biomass",
       main = full_trace_spec[[j]]$label, cex.main = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
# Damped spiral (decaying oscillation settling to a constant) => these cells
# are near a genuine bifurcation boundary. Pure monotonic decay/rise the
# whole way => they were never oscillating, and the earlier "amplitude"
# number was always just leftover convergence slope, not a real signal.

#Follow-up 5: re-run the same flat-looking cells with a much longer t_total
# to check whether they eventually start oscillating, settle to a genuinely
# flat constant, or are still moving even after a much longer window (i.e.
# whether 600 years was simply too short for this region of the grid).
plan(multisession)
results_long_run <- future_lapply(full_trace_spec, function(s) {
  test_beta_sigma_fixed(beta = s$beta, sigma = s$sigma,
                        resource_rate = best_rate, resource_capacity = best_capacity,
                        t_total = 3000)
}, future.seed = TRUE)

par(mfrow = c(2, 2), mar = c(3, 3, 2, 1))
for (j in seq_along(results_long_run)) {
  r <- results_long_run[[j]]
  plot(r$t, r$bm, type = "l", xlab = "Time (years)", ylab = "Biomass",
       main = paste0(full_trace_spec[[j]]$label, " (t_total = 3000)"), cex.main = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

par(mfrow = c(2, 2), mar = c(3, 3, 2, 1))
for (j in seq_along(results_long_run)) {
  r    <- results_long_run[[j]]
  mask <- r$t > 2900
  plot(r$t[mask], r$bm[mask], type = "l", xlab = "Time (years)", ylab = "Biomass",
       main = paste0(full_trace_spec[[j]]$label, " (t > 2900)"), cex.main = 0.8)
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
# If the tail (t > 2900) is still flat, these are genuine stable fixed
# points and 600 years was already enough time. If it's now oscillating or
# still visibly moving, 600 years was too short a window to trust the
# earlier amplitude heatmap for this part of the grid.

#Experiment: Apply fishing effort within the two confirmed genuinely-
# oscillating kernels (default beta=100/sigma=1.3, and the highest-amplitude
# beta=50/sigma=0.8), on the same fixed best-oscillating resource baseline
# used throughout today, to see how yield and biomass respond to effort.
# effort = 0 (no fishing) is the baseline the fished runs are compared
# against.

test_effort_fixed <- function(effort, beta, sigma, resource_rate, resource_capacity,
                              lambda = 2.05, t_total = 600) {
  p     <- make_params_fixed_resource(lambda = lambda, resource_rate = resource_rate,
                                      resource_capacity = resource_capacity,
                                      beta = beta, sigma = sigma)
  sim   <- make_limit_cycle_sim(p, t_total = t_total, effort = effort)
  bm    <- getBiomass(sim)[, "Anchovy"]
  yield <- getYield(sim)[, "Anchovy"]
  t     <- as.numeric(names(bm))
  list(sim = sim, bm = bm, yield = yield, t = t,
       effort = effort, beta = beta, sigma = sigma)
}

effort_seq <- seq(0, 6, length.out = 13)  # extended from 0-2: both curves were still climbing at effort = 2

fishing_kernels <- list(
  list(beta = 100, sigma = 1.3, label = "default (beta=100, sigma=1.3)"),
  list(beta = 50,  sigma = 0.8, label = "highest-amplitude (beta=50, sigma=0.8)")
)

plan(multisession)
fishing_results <- lapply(fishing_kernels, function(k) {
  runs <- future_lapply(effort_seq, function(e) {
    test_effort_fixed(effort = e, beta = k$beta, sigma = k$sigma,
                      resource_rate = best_rate, resource_capacity = best_capacity)
  }, future.seed = TRUE)
  list(kernel = k, runs = runs)
})

fishing_summary <- bind_rows(lapply(fishing_results, function(fr) {
  data.frame(
    kernel     = fr$kernel$label,
    effort     = effort_seq,
    mean_bm    = sapply(fr$runs, function(r) mean(r$bm[r$t > 550])),
    mean_yield = sapply(fr$runs, function(r) mean(r$yield[r$t > 550])),
    amp        = sapply(fr$runs, function(r) {
      late <- r$bm[r$t > 550]; max(late) - min(late)
    })
  )
}))
fishing_summary

# Yield vs effort. effort = 0 is the unfished baseline (yield = 0 by
# definition); the shape of the curve from there shows whether yield keeps
# climbing with effort or peaks and declines (a maximum sustainable yield).
ggplot(fishing_summary, aes(x = effort, y = mean_yield, color = kernel)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(x = "Fishing effort", y = "Mean yield (t > 550)",
       title = "Yield vs fishing effort in the oscillating regime",
       subtitle = "effort = 0 (leftmost point) is the unfished baseline") +
  theme_minimal()

# Biomass vs effort, with each kernel's own unfished (effort = 0) biomass
# drawn in as a dashed reference line -- the explicit baseline to compare
# fished biomass against.
baseline_bm <- fishing_summary[fishing_summary$effort == 0,
                               c("kernel", "mean_bm")]
names(baseline_bm)[2] <- "baseline_bm"

ggplot(fishing_summary, aes(x = effort, y = mean_bm, color = kernel)) +
  geom_line(linewidth = 1) + geom_point() +
  geom_hline(data = baseline_bm, aes(yintercept = baseline_bm, color = kernel),
             linetype = "dashed") +
  labs(x = "Fishing effort", y = "Mean biomass (t > 550)",
       title = "Biomass vs fishing effort in the oscillating regime",
       subtitle = "Dashed lines = each kernel's unfished (effort = 0) baseline biomass") +
  theme_minimal()

# Bonus: does fishing effort damp the oscillation itself, as the effort
# scans in Days 16-17 suggested? Reuses the same runs, no extra simulation.
ggplot(fishing_summary, aes(x = effort, y = amp, color = kernel)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(x = "Fishing effort", y = "Oscillation amplitude (max - min, t > 550)",
       title = "Does fishing effort dampen the oscillation?") +
  theme_minimal()

#Follow-up: the highest-amplitude kernel has more than double the default
# kernel's unfished biomass but produces roughly 8x LESS yield at every
# effort level. Gear selectivity (knife-edge, by default at w_mat) only
# catches fish above a size cutoff, so if that kernel's biomass sits mostly
# below the cutoff, high standing biomass wouldn't translate into catch.
# Reuses the sim objects already computed above at effort = 0 -- no new
# simulation needed.
unfished_default     <- fishing_results[[1]]$runs[[which(effort_seq == 0)]]$sim
unfished_highest_amp <- fishing_results[[2]]$runs[[which(effort_seq == 0)]]$sim

# Both plotSpectra() and plotSpectra2() hit internal errors with this mizer
# version, so the spectrum is built directly from the sim objects instead,
# using the same @n / @params@w accessors already used in
# make_limit_cycle_sim above. Biomass density = N(w) * w (power = 1
# convention), taken at the final saved timestep of each unfished run.
last_default <- dim(unfished_default@n)[1]
last_highamp <- dim(unfished_highest_amp@n)[1]

spectrum_df <- bind_rows(
  data.frame(
    w               = unfished_default@params@w,
    biomass_density = unfished_default@n[last_default, 1, ] * unfished_default@params@w,
    kernel          = "default (beta=100, sigma=1.3)"
  ),
  data.frame(
    w               = unfished_highest_amp@params@w,
    biomass_density = unfished_highest_amp@n[last_highamp, 1, ] * unfished_highest_amp@params@w,
    kernel          = "highest-amplitude (beta=50, sigma=0.8)"
  )
)

ggplot(spectrum_df, aes(x = w, y = biomass_density, color = kernel)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 10, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Body weight (g)", y = "Biomass density (N x w)",
       title = "Size spectra at effort = 0",
       subtitle = "Dashed line = w_mat = 10, the default gear knife-edge selectivity cutoff") +
  theme_minimal()
# If the highest-amplitude kernel's spectrum sits mostly to the left of the
# dashed line (below the selectivity cutoff), that's the biomass that's
# present but not catchable, and it explains the yield/biomass mismatch.

# The two curves cross a couple of times around the cutoff, so it's hard to
# tell by eye which kernel actually has more catchable (w >= 10) biomass.
# Quantify it directly instead of reading it off the log-log plot.
#
# n(w) is a density (individuals per unit weight), and the size grid is
# log-spaced, so summing n(w) * w across bins WITHOUT weighting by bin width
# dw is not a valid integral -- bins at large w are far wider than bins at
# small w in linear units. The correct biomass integral is
# sum(n(w) * w * dw), using params@dw (mizer's own bin-width slot, the same
# one getBiomass() uses internally).
catchable_fraction <- function(sim, cutoff = 10) {
  last         <- dim(sim@n)[1]
  w            <- sim@params@w
  dw           <- sim@params@dw
  bm_density   <- sim@n[last, 1, ] * w * dw
  total_bm     <- sum(bm_density)
  catchable_bm <- sum(bm_density[w >= cutoff])
  data.frame(total_bm = total_bm, catchable_bm = catchable_bm,
             catchable_fraction = catchable_bm / total_bm)
}

catchable_summary <- rbind(
  data.frame(kernel = "default (beta=100, sigma=1.3)",
             catchable_fraction(unfished_default)),
  data.frame(kernel = "highest-amplitude (beta=50, sigma=0.8)",
             catchable_fraction(unfished_highest_amp))
)
catchable_summary
getBiomass(finalParams(unfished_default),min_w = 10)

# If catchable_fraction is much lower for the highest-amplitude kernel, that
# confirms and quantifies the mechanism: most of its extra standing biomass
# is juveniles below the gear's selectivity cutoff, not catchable stock.
getParams(unfished_default)@species_params
plot(getResourceMort(getParams(unfished_default)))

#Experiment: Investigate why juveniles appear more impacted than adults, in
# the original Day 18-style regime (balance = TRUE, resource_decrease =
# r_d[1] = 0.001, default perturbation). Reuses the sim already computed in
# results_balance_true[[1]] above -- no new simulation needed.
# getBiomass() takes min_w/max_w directly, so there's no need to hand-roll
# the dw-weighted sum (which is exactly where the catchable_fraction bug was
# earlier) -- it does the correct integration internally.
track_life_stages <- function(sim, w_mat = 10) {
  bm_juv   <- getBiomass(sim, max_w = w_mat)[, "Anchovy"]
  bm_adult <- getBiomass(sim, min_w = w_mat)[, "Anchovy"]
  t        <- as.numeric(names(bm_juv))
  data.frame(t = t, bm_juv = bm_juv, bm_adult = bm_adult)
}

day18_sim     <- results_balance_true[[1]]$sim  # r_d = 0.001, balance = TRUE
life_stage_bm <- track_life_stages(day18_sim)
mask          <- life_stage_bm$t > 550

plot(NULL, xlim = c(550, 600),
     ylim = range(c(life_stage_bm$bm_juv[mask], life_stage_bm$bm_adult[mask])),
     xlab = "Time (years)", ylab = "Biomass",
     main = "Juvenile vs adult biomass (balance = TRUE, r_d = 0.001)")
lines(life_stage_bm$t[mask], life_stage_bm$bm_juv[mask],   col = "darkgreen", lwd = 2)
lines(life_stage_bm$t[mask], life_stage_bm$bm_adult[mask], col = "purple",    lwd = 2)
legend("topright", legend = c("Juveniles (w < 10)", "Adults (w >= 10)"),
       col = c("darkgreen", "purple"), lwd = 2)

# Absolute amplitude isn't a fair comparison -- juveniles and adults sit at
# very different biomass scales. Use relative amplitude (amp / mean) instead.
rel_amp <- function(x) (max(x) - min(x)) / mean(x)

life_stage_summary <- data.frame(
  stage   = c("Juvenile (w < 10)", "Adult (w >= 10)"),
  mean_bm = c(mean(life_stage_bm$bm_juv[mask]),  mean(life_stage_bm$bm_adult[mask])),
  amp     = c(max(life_stage_bm$bm_juv[mask]) - min(life_stage_bm$bm_juv[mask]),
             max(life_stage_bm$bm_adult[mask]) - min(life_stage_bm$bm_adult[mask])),
  rel_amp = c(rel_amp(life_stage_bm$bm_juv[mask]), rel_amp(life_stage_bm$bm_adult[mask]))
)
life_stage_summary
# If rel_amp is much larger for juveniles than adults, that confirms and
# quantifies "juveniles are more impacted" as a real feature of the Day
# 18-style oscillation, not an artefact of today's balance = FALSE detour.

#Recommended follow-up test A: deplete the resource selectively by size,
# rather than uniformly. resource_rate/resource_capacity are vectors over
# the resource's OWN size grid (juveniles and adults prefer different prey
# sizes -- an individual of weight w prefers prey around w/beta), but
# resource_decrease currently scales that whole vector by one constant.
# Depleting only the portion of the resource spectrum that juveniles feed on
# (small prey) vs. only the portion adults feed on (large prey) tests
# whether WHICH food source gets depleted drives the juvenile/adult
# asymmetry above, rather than just the overall amount of depletion.
#
# NOTE: params@w_full is this script's best guess at the resource's own
# weight-grid slot name -- verify against your installed mizer version
# (e.g. slotNames(p) or ?MizerParams) and adjust if it errors. getDiet() or
# getResourceMort() could also be used to empirically identify which prey
# sizes juveniles vs adults actually draw down, rather than the w/beta
# approximation used here.
make_params_localized_resource <- function(lambda = 2.05, resource_decrease = 0.001,
                                           deplete_w_range = c(0, Inf), balance = TRUE) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  r_full <- getResourceRate(params)   # vector over the resource's own size grid
  w_res  <- params@w_full             # verify this slot name -- see note above
  idx    <- w_res >= deplete_w_range[1] & w_res <= deplete_w_range[2]
  r_full[idx] <- r_full[idx] * resource_decrease

  setResource(params, resource_rate = r_full,
             resource_dynamics = "resource_semichemostat", balance = balance)
}

# beta = 100 by default, so juveniles (w < 10) prefer prey roughly w_full <
# 0.1, and adults (up to w_max = 66.5) prefer prey up to roughly w_full <
# 0.665. Depleting only the low end targets juvenile food; only the high end
# targets adult food.
p_deplete_juv_food   <- make_params_localized_resource(deplete_w_range = c(0, 0.1))
p_deplete_adult_food <- make_params_localized_resource(deplete_w_range = c(0.1, Inf))

sim_deplete_juv_food   <- make_limit_cycle_sim(p_deplete_juv_food,   t_total = 600)
sim_deplete_adult_food <- make_limit_cycle_sim(p_deplete_adult_food, t_total = 600)

summarize_life_stages <- function(sim, w_mat = 10) {
  ls   <- track_life_stages(sim, w_mat)
  m    <- ls$t > 550
  data.frame(juv_rel_amp = rel_amp(ls$bm_juv[m]), adult_rel_amp = rel_amp(ls$bm_adult[m]))
}

rbind(
  data.frame(scenario = "deplete juvenile food only", summarize_life_stages(sim_deplete_juv_food)),
  data.frame(scenario = "deplete adult food only",    summarize_life_stages(sim_deplete_adult_food))
)

#Recommended follow-up test B: perturb a different life stage. The two-stage
# perturbation in make_limit_cycle_sim always knocks down MATURE fish (w
# 10-100) at t = 10. Generalizing which size range gets shocked tests
# whether "juveniles more impacted" is a property of the cohort-wave
# mechanism itself (i.e. shows up regardless of where the shock lands), or
# specific to always shocking adults first.
make_limit_cycle_sim_custom <- function(params, t_total = 600, effort = 0,
                                        perturbation = 1e3, perturb_range = c(10, 100)) {
  params@initial_n_pp[] <- params@cc_pp * 0.1

  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= perturb_range[1] & params@w <= perturb_range[2]
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / perturbation

  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "tr_bdf2")
}

p_day18style <- make_params(resource_decrease = r_d[1], balance = TRUE)

sim_perturb_adult <- make_limit_cycle_sim_custom(p_day18style, perturb_range = c(10, 100))    # default, for comparison
sim_perturb_juv   <- make_limit_cycle_sim_custom(p_day18style, perturb_range = c(0.0003, 1))  # juveniles instead

rbind(
  data.frame(perturbed = "adults (default)", summarize_life_stages(sim_perturb_adult)),
  data.frame(perturbed = "juveniles",        summarize_life_stages(sim_perturb_juv))
)
# If juveniles still show a much larger rel_amp even when the shock is
# applied to THEM directly rather than to adults, that's strong evidence the
# cohort-wave mechanism itself amplifies juvenile fluctuations, rather than
# it being an artefact of always perturbing adults.
plotHover(getEGrowth(unfished_default,time_range=c(550,600)),log="xy")

#Reframing: rel_amp shows adults oscillate MORE in relative terms, not less
# -- so "juveniles are more impacted" isn't about oscillation amplitude at
# all. It's that there are several orders of magnitude more juveniles than
# adults, standing-stock-wise, in the first place. That's a demographic
# question (how much cumulative mortality does an individual experience
# while growing from w_min to w_mat?), not a dynamical one, and it's
# checkable on the UNPERTURBED steady state, before the two-stage kick that
# starts the oscillation at all.
#
# Standard size-spectrum survivorship-through-growth calculation: an
# individual spends dt = dw / g(w) time growing through each size bin, and
# accumulates hazard mu(w) * dt in that bin. Survivorship to size w is
# exp(-cumulative hazard). g() is somatic growth rate (getEGrowth), mu() is
# total mortality rate (getMort). There is no cannibalism in this sim (single
# species, no self-predation), so mu(w) here is background/senescence
# mortality plus fishing (zero in the unfished baseline), not predation from
# conspecifics -- worth decomposing directly (see getMort() components
# below) rather than assumed.
survivorship_curve <- function(params) {
  w      <- params@w
  dw     <- params@dw
  g      <- getEGrowth(params)[1, ]
  mu     <- getMort(params)[1, ]
  dt_bin <- dw / g

  cum_hazard   <- cumsum(mu * dt_bin)
  survivorship <- exp(-cum_hazard)

  data.frame(w = w, g = g, mu = mu, dt_bin = dt_bin,
             cum_hazard = cum_hazard, survivorship = survivorship)
}

p_day18style <- make_params(resource_decrease = r_d[1], balance = TRUE)
surv_depleted <- survivorship_curve(p_day18style)

# Comparison case: essentially undepleted resource (resource_decrease = 1),
# to check whether the steep survivorship drop is specific to the severe
# depletion used to trigger oscillation, or present regardless -- i.e.
# whether this is baked into growth/mortality calibration (lambda, alpha,
# gamma) alone.
p_undepleted <- make_params(resource_decrease = 1, balance = TRUE)
surv_undepleted <- survivorship_curve(p_undepleted)

plot(surv_depleted$w, surv_depleted$survivorship, type = "l", log = "xy",
     col = "firebrick", lwd = 2,
     xlab = "Body weight (g)", ylab = "Survivorship (fraction of w_min cohort reaching w)",
     main = "Survivorship-through-growth (unperturbed steady state)")
lines(surv_undepleted$w, surv_undepleted$survivorship, col = "steelblue", lwd = 2)
abline(v = 10, lty = 2)
legend("bottomleft",
       legend = c("resource_decrease = 0.001 (oscillating baseline)",
                 "resource_decrease = 1 (undepleted)"),
       col = c("firebrick", "steelblue"), lwd = 2, cex = 0.8)

surv_at_mat <- data.frame(
  scenario = c("resource_decrease = 0.001", "resource_decrease = 1"),
  survivorship_at_w_mat = c(
    surv_depleted$survivorship[which.min(abs(surv_depleted$w - 10))],
    surv_undepleted$survivorship[which.min(abs(surv_undepleted$w - 10))]
  )
)
surv_at_mat
# BUG: survivorship_at_w_mat came back bit-identical (1.944762e-10) for both
# scenarios. That's not a real finding -- getEGrowth(params)/getMort(params)
# on a freshly built params object use params@initial_n_pp, which
# setResource() never touches (it only changes the rr_pp/cc_pp parameters,
# not the current abundance state). initial_n_pp is set once by
# newSingleSpeciesParams() before setResource() runs, so both scenarios are
# being evaluated at the SAME default resource abundance regardless of
# resource_decrease -- the same root cause as the earlier
# check_balance_effect() dead end.
#
# Fix: getEGrowth(sim, time_range=)/getMort(sim, time_range=) don't return
# the same shape as getEGrowth(params)/getMort(params) -- indexing [1, ]
# on their output errors with "incorrect number of dimensions". Rather than
# guess at that return format, build a modified params object with
# initial_n/initial_n_pp set directly from the sim's own state (averaged
# over t > 550), and reuse the original survivorship_curve(params) function
# unchanged -- it already works correctly, using only the @n/@n_pp array
# accessors already used elsewhere in this script (e.g. in
# make_limit_cycle_sim's perturbation step).
params_at_sim_state <- function(sim, t_min = 550) {
  p    <- sim@params
  tvec <- as.numeric(dimnames(sim@n)[[1]])
  idx  <- which(tvec > t_min)

  p@initial_n[]    <- apply(sim@n[idx, , , drop = FALSE], c(2, 3), mean)
  p@initial_n_pp[] <- apply(sim@n_pp[idx, , drop = FALSE], 2, mean)
  p
}

sim_undepleted <- make_limit_cycle_sim(p_undepleted, t_total = 600)

p_depleted_dynamic   <- params_at_sim_state(day18_sim,      t_min = 550)
p_undepleted_dynamic <- params_at_sim_state(sim_undepleted, t_min = 550)

surv_depleted_dynamic   <- survivorship_curve(p_depleted_dynamic)
surv_undepleted_dynamic <- survivorship_curve(p_undepleted_dynamic)

plot(surv_depleted_dynamic$w, surv_depleted_dynamic$survivorship, type = "l", log = "xy",
     col = "firebrick", lwd = 2,
     xlab = "Body weight (g)", ylab = "Survivorship (fraction of w_min cohort reaching w)",
     main = "Survivorship-through-growth (dynamic resource state, t > 550)")
lines(surv_undepleted_dynamic$w, surv_undepleted_dynamic$survivorship, col = "steelblue", lwd = 2)
abline(v = 10, lty = 2)
legend("bottomleft",
       legend = c("resource_decrease = 0.001 (oscillating baseline)",
                 "resource_decrease = 1 (undepleted)"),
       col = c("firebrick", "steelblue"), lwd = 2, cex = 0.8)

surv_at_mat_dynamic <- data.frame(
  scenario = c("resource_decrease = 0.001 (dynamic)", "resource_decrease = 1 (dynamic)"),
  survivorship_at_w_mat = c(
    surv_depleted_dynamic$survivorship[which.min(abs(surv_depleted_dynamic$w - 10))],
    surv_undepleted_dynamic$survivorship[which.min(abs(surv_undepleted_dynamic$w - 10))]
  )
)
surv_at_mat_dynamic
# If survivorship_at_w_mat is now different between the two scenarios (and
# much lower for the depleted/oscillating case), that confirms resource
# depletion IS compounding an already-steep baseline decline. If they're
# still both extremely low but now genuinely differ, growth/mortality
# calibration sets the floor and depletion makes it worse; if they end up
# close again even with dynamic rates, the imbalance really is baked into
# growth/mortality calibration alone (background mortality + slow growth,
# not cannibalism -- see decomposition below).

#Follow-up: decompose getMort() into its components to identify what's
# actually driving the steep survivorship decline, now that cannibalism is
# ruled out (no self-predation in this single-species sim). Candidates:
# background/senescence mortality (mu_b, often allometric in w) vs fishing
# (zero here) vs how much dt_bin = dw/g(w) itself is inflated by slow growth
# (low alpha = 0.1 assimilation efficiency), which would mean individuals
# simply spend a very long time exposed at small sizes rather than being
# killed off unusually fast.
getMort(p_depleted_dynamic)[1, 1:10]     # inspect the first few size classes directly
getMort(p_undepleted_dynamic)[1, 1:10]

# Compare how much of cum_hazard's growth comes from mu vs from dt_bin
# (slow growth) by looking at which term dominates near w_mat.
idx_near_mat <- which.min(abs(surv_depleted_dynamic$w - 10))
surv_depleted_dynamic[max(1, idx_near_mat - 5):idx_near_mat, c("w", "g", "mu", "dt_bin")]
