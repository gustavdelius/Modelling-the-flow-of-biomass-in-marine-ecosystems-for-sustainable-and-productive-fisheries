library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)
library(glue)
library(ggplot2)
# ── Shared helper ──────────────────────────────────────────────────────────────
# Builds the standard single-species Anchovy params.
# kappa is scaled with lambda to keep total resource biomass constant.
# resource_decrease scales down the replenishment rate to create oscillations.
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


# ══ Experiment 1: steadyNewton() validation in the stable regime ══════════════
# At lambda = 2.05, the system converges to a fixed point under time-integration.
# steadyNewton() should agree with steady() to machine precision.

p_stable <- make_params(lambda = 2.05)
ss_old   <- steady(p_stable, t_max = 500)
ss_new   <- steadyNewton(p_stable)

mask     <- ss_old@initial_n[1, ] > 1e-20 & ss_new@initial_n[1, ] > 1e-20
rel_diff <- abs(ss_new@initial_n[1, mask] - ss_old@initial_n[1, mask]) /
              ss_old@initial_n[1, mask]
cat("Exp 1 — Max relative difference (non-zero bins):", max(rel_diff), "\n")

plot(ss_old@w, ss_old@initial_n[1, ], log = "xy", type = "l",
     col = "steelblue", lwd = 2,
     xlab = "Body mass (g)", ylab = "Number density",
     main = "Stable regime (lambda = 2.05): steady() vs steadyNewton()")
lines(ss_new@w, ss_new@initial_n[1, ], col = "firebrick", lwd = 2, lty = 2)
legend("topright", c("steady()", "steadyNewton()"),
       col = c("steelblue", "firebrick"), lwd = 2, lty = 1:2)


# ══ Experiment 2: Finding the unstable fixed point ════════════════════════════
# At lambda = 2.05 with low resource replenishment, the fixed point is unstable
# and the system oscillates. steadyNewton() finds it directly; time-integration
# would never converge to it.

p_osc       <- make_params(lambda = 2.05)
ss_unstable <- steadyNewton(p_osc)
bm_fp       <- sum(getBiomass(ss_unstable))

# Verify it is a genuine fixed point: project 1 year and check how much it moves.
sim_check  <- project(ss_unstable, t_max = 1, dt = 0.1, t_save = 1,
                      progress_bar = FALSE)
n_start    <- ss_unstable@initial_n[1, ]
n_end      <- sim_check@n[2, 1, ]
mask2      <- n_start > 1e-20
rel_change <- max(abs(n_end[mask2] - n_start[mask2]) / n_start[mask2], na.rm = TRUE)
cat("Exp 2 — Max relative change after 1 year from Newton state:", rel_change, "\n")
# Small -> genuine fixed point. It diverges from here because it is *unstable*.

# Confirm instability: project 300 years from the Newton state.
# Tiny numerical noise should grow and the trajectory should drift to the limit cycle.
sim_from_fp <- project(ss_unstable, t_max = 300, dt = 0.1, t_save = 1,
                       progress_bar = FALSE)
bm_from_fp  <- getBiomass(sim_from_fp)[, "Anchovy"]
t_from_fp   <- as.numeric(names(bm_from_fp))

plot(t_from_fp, bm_from_fp, type = "l",
     xlab = "Time (years)", ylab = "Biomass",
     main = "From Newton fixed point — unstable: drifts to limit cycle")
abline(h = bm_fp, col = "firebrick", lty = 2)


# ══ Experiment 3: Limit cycle dynamics and comparison to the fixed point ══════
# The perturbation (depleted resource + mature fish knocked down at t=10)
# kicks the system onto the limit cycle. We then compare the time-averaged
# biomass on the cycle against the unstable fixed point.

sim_long <- make_limit_cycle_sim(p_osc, t_total = 600)
bm2      <- getBiomass(sim_long)[, "Anchovy"]
t2       <- as.numeric(names(bm2))
bm_avg   <- exp(mean(log(bm2[t2 > 400])))

cat("Exp 3 — Fixed point biomass:   ", round(bm_fp,  6), "\n")
cat("Exp 3 — Time-averaged biomass: ", round(bm_avg, 6), "\n")
cat("Exp 3 — Ratio (avg / fp):      ", round(bm_avg / bm_fp, 4), "\n")

plot(t2, bm2, type = "l",
     xlab = "Time (years)", ylab = "Biomass",
     main = "Limit cycle vs unstable fixed point")
abline(h = bm_fp,  col = "firebrick", lty = 2, lwd = 2)
abline(h = bm_avg, col = "steelblue", lty = 2, lwd = 2)
legend("topright",
       c("Trajectory", "Unstable fixed point", "Time average"),
       col = c("black", "firebrick", "steelblue"),
       lty = c(1, 2, 2), lwd = 2)

# Power spectrum of settled oscillation to extract the dominant period.
signal <- bm2[t2 > 400]
n_sig  <- length(signal)
dt_sig <- 0.2
spec   <- Mod(fft(signal - mean(signal)))^2
half   <- seq_len(floor(n_sig / 2))
ps     <- data.frame(period = (n_sig * dt_sig) / half, power = spec[half + 1])
cat("Exp 3 — Dominant period (years):", ps$period[which.max(ps$power)], "\n")


# ══ Experiment 4: Lambda scan — is the ~2× ratio robust? ═════════════════════
# Scan lambda from 1.9 to 2.2 and compute the fixed-point biomass and
# time-averaged limit-cycle biomass at each value.

scan_lambda <- function(lam) {
  p <- make_params(lambda = lam)

  sim  <- make_limit_cycle_sim(p, t_total = 600)
  bm   <- getBiomass(sim)[, "Anchovy"]
  t    <- as.numeric(names(bm))

  ss    <- tryCatch(steadyNewton(p), error = function(e) NULL)
  bm_fp_l <- if (!is.null(ss)) sum(getBiomass(ss)) else NA

  data.frame(lambda      = lam,
             fixed_point = bm_fp_l,
             time_avg    = mean(bm[t > 400]),
             ratio       = mean(bm[t > 400]) / bm_fp_l)
}

lambdas    <- seq(1.9, 2.2, by = 0.025)
df_lambdas <- do.call(rbind, lapply(lambdas, scan_lambda))
n_avg_spec <- apply(sim_long@n[t_idx, 1, ], 2, mean)
n_min_spec <- apply(sim_long@n[t_idx, 1, ], 2, min)
n_max_spec <- apply(sim_long@n[t_idx, 1, ], 2, max)

plot(df_lambdas$lambda, df_lambdas$ratio, type = "b", pch = 19,
     xlab = "lambda", ylab = "Time avg / fixed point",
     main = "Limit-cycle advantage across lambda")
abline(h = 1, lty = 2, col = "grey50")

# 1. Print the actual fixed-point and time-avg numbers from the scan
print(df_lambdas)


# ══ Experiment 5: Duty-cycle decomposition of the ~2× ratio ══════════════════
# The time average above the fixed point can be decomposed exactly into the
# fraction of time spent above fp and the conditional means above/below.
# This accounts for the full ratio with no residual.

settled    <- bm2[t2 > 400]
f_above    <- mean(settled > bm_fp)
mean_above <- mean(settled[settled > bm_fp])
mean_below <- mean(settled[settled <= bm_fp])
ratio_pred <- (f_above * mean_above + (1 - f_above) * mean_below) / bm_fp

cat("Exp 5 — Fraction of time above fp:", round(f_above,    3), "\n")
cat("Exp 5 — Mean above fp:            ", round(mean_above, 4), "\n")
cat("Exp 5 — Mean below fp:            ", round(mean_below, 4), "\n")
cat("Exp 5 — Predicted ratio:          ", round(ratio_pred, 3), "\n")
cat("Exp 5 — Actual ratio:             ", round(mean(settled) / bm_fp, 3), "\n")

hist(settled, breaks = 50, col = "steelblue",
     xlab = "Biomass", main = "Time spent at each biomass level")
abline(v = bm_fp,         col = "firebrick", lwd = 2, lty = 2)
abline(v = mean(settled), col = "darkblue",  lwd = 2, lty = 2)
legend("topright", c("Fixed point", "Time average"),
       col = c("firebrick", "darkblue"), lwd = 2, lty = 2)
#Question this plots usefulness

# ══ Experiment 6: Size spectrum — fixed point vs limit cycle ══════════════════
# Compare n(w) from steadyNewton() against the time-averaged n(w,t) from the
# limit cycle simulation. The pointwise ratio reveals whether the cycle is a
# uniform scaling or has size-dependent structure.

w          <- p_osc@w
n_fp_spec  <- ss_unstable@initial_n[1, ]
t_all      <- as.numeric(dimnames(sim_long@n)[[1]])
t_idx      <- which(t_all > 400)
n_avg_spec <- apply(sim_long@n[t_idx, 1, ], 2, mean)
n_min_spec <- apply(sim_long@n[t_idx, 1, ], 2, min)
n_max_spec <- apply(sim_long@n[t_idx, 1, ], 2, max)

# Plot 1: spectra on log-log
plot(w, n_fp_spec, type = "l", log = "xy", col = "firebrick", lwd = 2,
     xlab = "Body mass w (g)", ylab = "Number density n(w)",
     main = "Size spectrum: fixed point vs time-averaged limit cycle")
lines(w, n_avg_spec, col = "steelblue", lwd = 2)
legend("topright", c("Fixed point (steadyNewton)", "Time-avg limit cycle"),
       col = c("firebrick", "steelblue"), lwd = 2)

# Plot 2: pointwise ratio — flat means uniform scaling, slope means size-dependent
plot(w, n_avg_spec / n_fp_spec, type = "l", log = "x",
     col = "darkgreen", lwd = 2,
     xlab = "Body mass w (g)", ylab = expression(bar(n)(w) / n[fp](w)),
     main = "Ratio: limit-cycle average / fixed point")
abline(h = bm_avg / bm_fp, lty = 2, col = "steelblue")
abline(h = 1,               lty = 2, col = "grey50")
legend("bottomright",
       c("Pointwise ratio",
         paste0("Total biomass ratio (~", round(bm_avg / bm_fp, 2), "x)"),
         "Ratio = 1"),
       col = c("darkgreen", "steelblue", "grey50"), lty = c(1, 2, 2), lwd = 2)

# Plot 3: min–max envelope across the cycle
plot(w, n_fp_spec, type = "n", log = "xy",
     xlab = "Body mass w (g)", ylab = "Number density n(w)",
     main = "Spectrum envelope across the limit cycle",
     ylim = range(c(n_min_spec[n_min_spec > 0], n_max_spec), na.rm = TRUE))
polygon(c(w, rev(w)), c(n_max_spec, rev(n_min_spec)),
        col = adjustcolor("steelblue", alpha.f = 0.25), border = NA)
lines(w, n_avg_spec, col = "steelblue", lwd = 2)
lines(w, n_fp_spec,  col = "firebrick", lwd = 2, lty = 2)
legend("topright",
       c("Fixed point", "Time average", "Min-max envelope"),
       col = c("firebrick", "steelblue", "steelblue"),
       lty = c(2, 1, 1), lwd = c(2, 2, 8))


# ══ Experiment 7: Constant fishing scan ═══════════════════════════════════════
# Does fishing destroy the limit cycle or reduce mean biomass significantly?
# Uses the same limit-cycle initial condition throughout so that differences
# are due to fishing mortality, not initialisation.

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

efforts_fishing <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
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
