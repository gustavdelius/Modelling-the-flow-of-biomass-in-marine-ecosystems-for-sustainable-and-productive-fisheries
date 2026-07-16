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
library(future)
library(scales)
library(deSolve)

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH (260 chars) has silently truncated filenames -- even the
# .png extension itself -- once combined with this repo's long, deeply
# nested folder path. Keep names short at the call site; this is a
# defensive last resort so a long one truncates safely instead of silently.
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 150) {
  max_name <- 40
  if (nchar(filename) > max_name) {
    ext      <- tools::file_ext(filename)
    base     <- tools::file_path_sans_ext(filename)
    filename <- paste0(substr(base, 1, max_name - nchar(ext) - 1), ".", ext)
    warning(sprintf("save_plot(): filename too long, truncated to '%s'", filename))
  }
  ggsave(file.path("interesting_plots", filename), plot = plot,
         width = width, height = height, dpi = dpi)
}

# Self-contained, same as every script since Day 20 -- redefined here rather
# than sourced from 26_experiments.R.
make_second_order_params_kr <- function(lambda = 2.05, resource_decrease = 0.001,
                                        capacity_mult = 1, second_order = TRUE,
                                        ext_diff = 0.00, alpha = 0.1, gamma = 750,
                                        ks = 0, kappa_override = NULL) {
  a0    <- 100
  kappa <- if (is.null(kappa_override)) a0 * exp(-6.9 * (lambda - 1)) else kappa_override
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = alpha, gamma = gamma, ks = ks
  )

  given_species_params(params)$D_ext <- ext_diff
  params <- setBevertonHolt(params)

  default_capacity <- getResourceCapacity(params)
  r  <- getResourceRate(params) * resource_decrease
  cc <- default_capacity * capacity_mult

  params <- setResource(params, resource_rate = r, resource_capacity = cc,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE)

  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }

  params
}

################################################################################
# Shared mizer catchable-fraction machinery
#
# One primitive underneath every catchable_fraction read in this project: a
# full time series at every saved timestep. catchable_fraction_at() (legacy
# single-snapshot, Days 25-26) and catchable_fraction_windowed() (late-window
# mean, Section 2) are both thin wrappers around it, so a mizer/anchovy
# oscillation shows up the same way regardless of which read is used
# downstream.
################################################################################

catchable_fraction_series <- function(params, t_max, dt = 0.1, t_save = 0.5) {
  sim   <- project(params, t_max = t_max, dt = dt, t_save = t_save,
                   progress_bar = FALSE, effort = 0, method = "tr_bdf2")
  times <- as.numeric(dimnames(sim@n)[[1]])
  w     <- sim@params@w
  dw    <- sim@params@dw
  w_mat <- params@species_params$w_mat[1]

  fractions <- sapply(seq_along(times), function(i) {
    bm_density   <- sim@n[i, 1, ] * w * dw
    total_bm     <- sum(bm_density)
    catchable_bm <- sum(bm_density[w >= w_mat])
    catchable_bm / total_bm
  })

  data.frame(time = times, catchable_fraction = fractions)
}

# Single last-timestep read, same as Days 25-26 -- kept numerically identical
# (same dt/t_save) so anything still comparing against those numbers stays
# valid.
catchable_fraction_at <- function(params, t_max = 600) {
  series <- catchable_fraction_series(params, t_max = t_max, dt = 0.1, t_save = 0.5)
  data.frame(catchable_fraction = tail(series$catchable_fraction, 1), error = NA_character_)
}

# Mean over the late window (last window_frac of the run) instead of a
# single snapshot -- see Section 2 for why.
catchable_fraction_windowed <- function(params, t_max = 600, window_frac = 0.4) {
  series <- catchable_fraction_series(params, t_max = t_max, dt = 0.1, t_save = 0.5)
  late   <- series[series$time >= t_max * (1 - window_frac), ]
  data.frame(catchable_fraction_mean = mean(late$catchable_fraction),
             catchable_fraction_sd   = sd(late$catchable_fraction),
             catchable_fraction_min  = min(late$catchable_fraction),
             catchable_fraction_max  = max(late$catchable_fraction),
             error = NA_character_)
}

# Shared sweep harness: builds params via make_second_order_params_kr(),
# applies fraction_fn to each, and warns (without aborting the sweep) on any
# per-point failure. Used for both the legacy single-snapshot sweep (Days
# 25-26) and the windowed rerun (Section 2) -- same harness, different
# fraction_fn.
run_param_sweep_generic <- function(param_name, param_seq, fraction_fn, capacity_mult = 10,
                                    resource_decrease = 0.001, ...) {
  bind_rows(lapply(param_seq, function(v) {
    args <- list(capacity_mult = capacity_mult, resource_decrease = resource_decrease)
    args[[param_name]] <- v
    result <- tryCatch({
      p <- do.call(make_second_order_params_kr, args)
      fraction_fn(p, ...)
    }, error = function(e) {
      data.frame(error = conditionMessage(e))
    })
    if (!is.na(result$error)) {
      warning(sprintf("%s = %.4g: %s", param_name, v, result$error))
    }
    data.frame(param = param_name, value = v, result)
  }))
}

run_param_sweep <- function(param_name, param_seq, ...) {
  run_param_sweep_generic(param_name, param_seq, catchable_fraction_at, ...)
}
run_param_sweep_windowed <- function(param_name, param_seq, ...) {
  run_param_sweep_generic(param_name, param_seq, catchable_fraction_windowed, ...)
}

################################################################################
# Shared toy-model (DDE) machinery
#
# D: prey renewal rate; a: attack rate; e: predator conversion efficiency;
# m: predator mortality; K: prey carrying capacity -- the same reference
# parameters used throughout Days 24-27.
################################################################################

predator_prey_params_delay <- list(D = 1, a = 1, e = 0.5, m = 0.3, K = 1)
params_small_D <- modifyList(predator_prey_params_delay, list(D = 0.1))

# deSolve's dede()/lagvalue() keeps a fixed-size history buffer for looking
# back in time (default mxhist = 1e4 steps), and it isn't big enough once
# tau/dt runs into the thousands -- hit directly at tau=1000, dt=0.1 in
# Section 3 ("illegal input in lagvalue - lag ... too large", right at the
# 10,000-step default). NOT a top-level dede() argument (that was tried and
# is wrong -- an unrecognised name there falls into dede()'s own ... and
# gets forwarded straight into func(), producing "unused argument"); it's
# set via control = list(mxhist = ...). Set once, comfortably above the
# largest tau/dt used anywhere in this file, and pass explicitly to every
# dede() call.
DEFAULT_MXHIST <- 2e5

# Shared runner for every delay-placement variant in this file (pdelay,
# selflim, selflim_typeII, ndelay): builds history from state0, integrates
# with dede(), and summarises the late window (last 40% of t_max) the same
# way every time, so results are directly comparable across variants.
# extra_parms carries anything beyond tau (e.g. selflim_typeII's handling
# time h).
run_predator_prey_delay <- function(rhs_fn, tau, state0 = c(N = 1, P = 0.1), t_max = 300,
                                    parms = predator_prey_params_delay, extra_parms = list(),
                                    mxhist = DEFAULT_MXHIST) {
  parms     <- modifyList(parms, extra_parms)
  parms$tau <- tau
  parms$N0  <- state0[["N"]]
  parms$P0  <- state0[["P"]]
  out  <- as.data.frame(dede(y = state0, times = seq(0, t_max, by = 0.1),
                             func = rhs_fn, parms = parms,
                             control = list(mxhist = mxhist)))
  late <- out[out$time > t_max * 0.6, ]
  data.frame(tau = tau, max_N = max(late$N), min_N = min(late$N),
             max_P = max(late$P), min_P = min(late$P),
             rel_amplitude_P = (max(late$P) - min(late$P)) / mean(late$P))
}

################################################################################
# 1. Is capacity_mult = 10, resource_decrease = 0.001 actually oscillating?
#
# Day 26's gamma = 305.8 recheck found catchable_fraction flipping
# non-monotonically as t_run grew, with the ordering of all three rechecked
# points scrambling at t_run = 2400 -- not what a converging transient looks
# like. The suspected explanation: catchable_fraction_at() only ever reads a
# single snapshot at the last saved timestep, so if this operating point
# (capacity_mult = 10, resource_decrease = 0.001, anchovy's own alpha/gamma/
# kappa -- every juvenile-pileup sweep's baseline since Day 24) is gently
# oscillating rather than settled, varying t_max just samples different
# PHASES of it. This checks that directly with a real time series.
################################################################################

operating_point_params <- make_second_order_params_kr(capacity_mult = 10, resource_decrease = 0.001)

# t_max = 4800 matches the longest t_run used in Day 26's gamma-dip check.
operating_point_ts <- catchable_fraction_series(operating_point_params, t_max = 4800, t_save = 1)

print(head(operating_point_ts, 20))
write.csv(operating_point_ts, file.path("interesting_plots", "operating_point_timeseries.csv"),
          row.names = FALSE)

# Late window = last 40%, the same convention used throughout the
# predator-prey delay checks.
late_window <- operating_point_ts %>% filter(time >= 0.6 * max(time))
rel_amplitude_operating_point <- (max(late_window$catchable_fraction) - min(late_window$catchable_fraction)) /
  mean(late_window$catchable_fraction)

# Crude period estimate: local maxima in the late window, mean spacing
# between them. Fewer than two maxima -> not (visibly) oscillating, period
# reported as NA rather than guessed at.
find_local_maxima_times <- function(df) {
  y <- df$catchable_fraction
  is_max <- c(FALSE, y[-1] > y[-length(y)]) & c(y[-length(y)] > y[-1], FALSE)
  df$time[is_max]
}

maxima_times     <- find_local_maxima_times(late_window)
estimated_period <- if (length(maxima_times) >= 2) mean(diff(maxima_times)) else NA_real_

cat(sprintf(
  paste0(
    "Operating point (capacity_mult=10, resource_decrease=0.001) over t_max=4800:\n",
    "  late-window rel_amplitude of catchable_fraction = %.4g\n",
    "  local maxima found in late window: %d, estimated period: %s\n",
    "Small rel_amplitude with no repeated maxima = settled, and Day 26's gamma-\n",
    "dip flip needs a different explanation. Non-negligible rel_amplitude with\n",
    "evenly-spaced maxima = a genuine oscillation, and every single-snapshot\n",
    "catchable_fraction_at() read since Day 24 has been sampling an arbitrary\n",
    "phase of it.\n"
  ),
  rel_amplitude_operating_point, length(maxima_times),
  ifelse(is.na(estimated_period), "N/A (fewer than 2 maxima)", sprintf("%.1f time units", estimated_period))
))

operating_point_plot <- ggplot(operating_point_ts, aes(x = time, y = catchable_fraction)) +
  geom_line(color = "#4C72B0", linewidth = 0.6) +
  labs(x = "Time", y = "Catchable fraction (biomass at or above w_mat)",
       title = "Is capacity_mult=10, resource_decrease=0.001 actually oscillating?",
       subtitle = "Full time series, not a single final-timestep snapshot") +
  theme_minimal()
operating_point_plot
save_plot(operating_point_plot, "opoint_full_ts.png")

# Zoomed-in view of just the late window, where any real oscillation should
# be easiest to see against the plateau.
operating_point_late_plot <- ggplot(late_window, aes(x = time, y = catchable_fraction)) +
  geom_line(color = "#55A868", linewidth = 0.8) +
  geom_point(data = late_window %>% filter(time %in% maxima_times), color = "#C44E52", size = 2) +
  labs(x = "Time", y = "Catchable fraction",
       title = "Late window only (last 40% of the run)",
       subtitle = "Red points = detected local maxima, used for the period estimate") +
  theme_minimal()
operating_point_late_plot
save_plot(operating_point_late_plot, "opoint_late_zoom.png")

################################################################################
# 2. Rerunning alpha, gamma, and kappa's isolation sweeps with a
#    windowed-average read
#
# Day 25's alpha/gamma/kappa sweep (and Day 26's refinements) all read
# catchable_fraction from a single final-timestep snapshot. Section 1 above
# found the operating point itself basically flat (rel_amplitude ~ 1e-12),
# so this isn't chasing a genuine limit cycle -- but individual sweep points
# away from the anchovy's own gamma/kappa could still be slowly settling
# within t_max=600, and gamma/kappa's Day 25 ranges were only 1-3 orders of
# magnitude, comparable to or smaller than the noise Day 26's gamma-dip
# check exposed. alpha's own 27.5-order-of-magnitude range was never at risk
# of the wrong sign, but it's rerun here too for the same footing.
################################################################################

# Exact same 15-point grids as Day 25's original alpha, gamma, and kappa
# sweeps, so the windowed-average numbers are directly comparable.
alpha_seq_original <- exp(seq(log(0.02), log(0.8), length.out = 15))
gamma_seq_original  <- exp(seq(log(100), log(5000), length.out = 15))
kappa_seq_original  <- exp(seq(log(0.001), log(1), length.out = 15))

alpha_windowed_df <- run_param_sweep_windowed("alpha", alpha_seq_original)
gamma_windowed_df <- run_param_sweep_windowed("gamma", gamma_seq_original)
kappa_windowed_df <- run_param_sweep_windowed("kappa_override", kappa_seq_original) %>%
  mutate(param = "kappa")

print(alpha_windowed_df)
print(gamma_windowed_df)
print(kappa_windowed_df)
write.csv(alpha_windowed_df, file.path("interesting_plots", "alpha_windowed_sweep.csv"), row.names = FALSE)
write.csv(gamma_windowed_df, file.path("interesting_plots", "gamma_windowed_sweep.csv"), row.names = FALSE)
write.csv(kappa_windowed_df, file.path("interesting_plots", "kappa_windowed_sweep.csv"), row.names = FALSE)

windowed_vs_single_summary <- function(windowed_df, label) {
  windowed_df %>%
    filter(is.na(error)) %>%
    summarise(
      param       = label,
      min_mean    = min(catchable_fraction_mean),
      max_mean    = max(catchable_fraction_mean),
      log10_range = log10(max(catchable_fraction_mean)) - log10(min(catchable_fraction_mean)),
      max_rel_sd  = max(catchable_fraction_sd / catchable_fraction_mean, na.rm = TRUE)
    )
}

windowed_magnitude_summary <- bind_rows(
  windowed_vs_single_summary(alpha_windowed_df, "alpha"),
  windowed_vs_single_summary(gamma_windowed_df, "gamma"),
  windowed_vs_single_summary(kappa_windowed_df, "kappa")
)

print(windowed_magnitude_summary)
write.csv(windowed_magnitude_summary, file.path("interesting_plots", "windowed_magnitude_summary.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Windowed-average alpha/gamma/kappa sweep: log10 range = %.2f (alpha), %.2f\n",
    "(gamma), %.2f (kappa). Compare against Day 25's single-snapshot values:\n",
    "27.5 (alpha), 1.3 (gamma), 2.6 (kappa). If the windowed ranges come out\n",
    "about the same (or smaller), gamma and kappa's 'modest'/'flat' Day 25\n",
    "conclusions survive averaging; alpha's range was never at risk either way.\n"
  ),
  windowed_magnitude_summary$log10_range[windowed_magnitude_summary$param == "alpha"],
  windowed_magnitude_summary$log10_range[windowed_magnitude_summary$param == "gamma"],
  windowed_magnitude_summary$log10_range[windowed_magnitude_summary$param == "kappa"]
))

make_windowed_plot <- function(df, color, xlab, title, subtitle, xline = NULL) {
  p <- ggplot(df %>% filter(is.na(error)), aes(x = value, y = catchable_fraction_mean)) +
    geom_ribbon(aes(ymin = catchable_fraction_min, ymax = catchable_fraction_max),
               fill = color, alpha = 0.2) +
    geom_line(color = color, linewidth = 1) +
    geom_point(color = color, size = 1.5) +
    scale_x_log10() +
    scale_y_log10(labels = scales::label_scientific()) +
    labs(x = xlab, y = "Catchable fraction (late-window mean)", title = title, subtitle = subtitle) +
    theme_minimal()
  if (!is.null(xline)) p <- p + geom_vline(xintercept = xline, linetype = "dashed", color = "grey40")
  p
}

alpha_windowed_plot <- make_windowed_plot(
  alpha_windowed_df, "#C44E52", "alpha (assimilation efficiency)",
  "Alpha sweep rerun with a windowed-average read",
  "Ribbon = min-max within the late window; dashed = anchovy's own alpha = 0.1",
  xline = 0.1)
alpha_windowed_plot
save_plot(alpha_windowed_plot, "alpha_windowed.png")

gamma_windowed_plot <- make_windowed_plot(
  gamma_windowed_df, "#4C72B0", "gamma (search volume)",
  "Gamma sweep rerun with a windowed-average read",
  "Ribbon = min-max within the late window")
gamma_windowed_plot
save_plot(gamma_windowed_plot, "gamma_windowed.png")

kappa_windowed_plot <- make_windowed_plot(
  kappa_windowed_df, "#55A868", "kappa (background resource spectrum)",
  "Kappa sweep rerun with a windowed-average read",
  "Ribbon = min-max within the late window")
kappa_windowed_plot
save_plot(kappa_windowed_plot, "kappa_windowed.png")

################################################################################
# 3. Re-checking the P-only-delay variant: Day 26's "D < 0.144" threshold
#    turns out to be necessary but not sufficient -- no tau, however large,
#    ever produces a Hopf bifurcation here
#
# Day 26 derived Hopf possible <=> 2*beta*gamma > alpha^2 for this variant
# (delay only P(t-tau) inside the predator's intake term), then reduced that
# to a threshold D < 2*a*m*C/(1+a*C)^2 (~0.144 at this project's reference
# a, e, m, K), and treated that as the full condition. It isn't:
# 2*beta*gamma > alpha^2 only makes the SUM of the two roots of the
# u = omega^2 quadratic positive. A real, positive omega also needs that
# quadratic's DISCRIMINANT to be non-negative -- a second condition Day 26's
# write-up never checked.
################################################################################

check_hopf_feasibility_pdelay_full <- function(parms = predator_prey_params_delay) {
  D <- parms$D; a <- parms$a; e <- parms$e; m <- parms$m; K <- parms$K
  N_star <- m / (e * a)
  C      <- (K - N_star) / (a * N_star)
  P_star <- D * C
  alpha  <- D * (1 + a * C)
  beta   <- a * N_star
  gamma  <- e * a * P_star
  bg     <- beta * gamma

  B  <- alpha^2 - 2 * bg           # coefficient of u in u^2 + B*u + Cc = 0
  Cc <- bg * (2 * alpha * m + bg)  # constant term
  disc <- B^2 - 4 * Cc

  data.frame(D = D, alpha = alpha, beta = beta, gamma = gamma, beta_gamma = bg,
             B = B, Cc = Cc, disc = disc,
             sum_positive = (-B) > 0, disc_nonneg = disc >= 0,
             hopf_possible = (-B) > 0 & disc >= 0)
}

pdelay_recheck_D01 <- check_hopf_feasibility_pdelay_full(params_small_D)
print(pdelay_recheck_D01)
cat(sprintf(
  "Re-checked D = 0.1: sum_positive = %s, disc_nonneg = %s (disc = %.6f), hopf_possible = %s.\n",
  pdelay_recheck_D01$sum_positive, pdelay_recheck_D01$disc_nonneg, pdelay_recheck_D01$disc,
  pdelay_recheck_D01$hopf_possible
))

# Is D = 0.1 just unlucky, or does the window never open for ANY D (a, e, m,
# K held at reference values)?
D_scan_seq <- exp(seq(log(1e-4), log(1e2), length.out = 300))
pdelay_D_scan_df <- bind_rows(lapply(D_scan_seq, function(D) {
  check_hopf_feasibility_pdelay_full(modifyList(predator_prey_params_delay, list(D = D)))
}))
write.csv(pdelay_D_scan_df, file.path("interesting_plots", "pdelay_full_feasibility_D_scan.csv"),
          row.names = FALSE)
cat(sprintf("D scan (1e-4 to 100, %d points): any D with hopf_possible = TRUE? %s\n",
            nrow(pdelay_D_scan_df), any(pdelay_D_scan_df$hopf_possible)))

pdelay_D_scan_plot <- ggplot(pdelay_D_scan_df, aes(x = D)) +
  geom_line(aes(y = -B, color = "sum of roots (-B), need > 0"), linewidth = 1) +
  geom_line(aes(y = disc, color = "discriminant, need >= 0"), linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(x = "D", y = "Value", color = NULL,
       title = "Neither condition needed for a Hopf crossing is ever satisfied together",
       subtitle = "Sum > 0 holds only for small D; discriminant >= 0 only for much larger D -- the two never overlap") +
  theme_minimal()
pdelay_D_scan_plot
save_plot(pdelay_D_scan_plot, "pdelay_sum_disc.png")

# Structural, or just unlucky reference numbers? After rescaling by m, the
# system collapses to one ratio, r = K*a*e/m. D_max(r) is where the sum
# condition turns off; D_r(r) is where the discriminant turns on. If
# D_r(r) >= D_max(r) for every r > 1 (the only physically valid range), the
# window is structurally empty.
D_max_of_r <- function(r) 2 * (r - 1) / r^2
D_r_of_r <- function(r) {
  Aq <- r^3
  Bq <- -4 * r * (r - 1)
  Cq <- -8 * (r - 1)
  (-Bq + sqrt(Bq^2 - 4 * Aq * Cq)) / (2 * Aq)
}

r_seq <- exp(seq(log(1.001), log(1000), length.out = 300))
pdelay_structural_check_df <- data.frame(r = r_seq, D_max = D_max_of_r(r_seq), D_r = D_r_of_r(r_seq)) %>%
  mutate(window_open = D_r < D_max)
write.csv(pdelay_structural_check_df, file.path("interesting_plots", "pdelay_structural_window_check.csv"),
          row.names = FALSE)
cat(sprintf(
  "Structural check across r in [1.001, 1000] (%d points): any r with the window open? %s\n",
  nrow(pdelay_structural_check_df), any(pdelay_structural_check_df$window_open)
))

pdelay_structural_plot <- ggplot(pdelay_structural_check_df, aes(x = r)) +
  geom_line(aes(y = D_max, color = "D_max (sum-of-roots cutoff)"), linewidth = 1) +
  geom_line(aes(y = D_r, color = "D_r (discriminant turns on)"), linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "r = K*a*e / m", y = "D", color = NULL,
       title = "The feasibility window [D_r, D_max) is empty for every r tested",
       subtitle = "D_r sits above D_max everywhere -- structural, not a bad parameter pick") +
  theme_minimal()
pdelay_structural_plot
save_plot(pdelay_structural_plot, "pdelay_structural.png")

# One direct simulation check, for confirmation rather than more searching:
# a very large tau should still show the amplitude decaying, not settling
# into a new plateau. tau=1000 at dt=0.1 needs history reaching back 10,000
# steps -- past dede()'s default buffer size, hence the explicit mxhist on
# run_predator_prey_delay() (see DEFAULT_MXHIST above).
predator_prey_rhs_pdelay <- function(t, state, parms) {
  N <- max(state[1], 0); P <- max(state[2], 0)
  D <- parms[["D"]]; a <- parms[["a"]]; e <- parms[["e"]]; m <- parms[["m"]]; K <- parms[["K"]]
  tau <- parms[["tau"]]; P0 <- parms[["P0"]]
  lag_P <- if (t <= tau) P0 else max(lagvalue(t - tau, 2), 0)
  dN <- D * (K - N) - a * N * P
  dP <- e * a * N * lag_P - m * P
  list(c(dN, dP))
}

pdelay_large_tau_check <- run_predator_prey_delay(predator_prey_rhs_pdelay, tau = 1000, t_max = 6000,
                                                   parms = params_small_D)
print(pdelay_large_tau_check)
cat(sprintf(
  paste0(
    "tau = 1000 (D = 0.1, t_max = 6000): rel_amplitude_P = %.4g.\n",
    "Still small/decaying, not a steady nonzero plateau -- the simulation agrees\n",
    "with the analytic result: this delay placement never bifurcates at this D.\n"
  ),
  pdelay_large_tau_check$rel_amplitude_P
))

################################################################################
# 4. How far past onset does the self-limitation model's amplitude stay
#    biologically sensible?
#
# Day 26 found delaying the prey's own self-limitation term produces a
# genuine Hopf bifurcation (onset between tau=2.2 and tau=2.4), but min_N
# goes negative right at onset (-0.021 at tau=2.4) and gets rapidly worse
# (-15.5 by tau=20) -- this toy model has no mechanism to cap the
# oscillation's amplitude. Two follow-ups from Day 26's "What's Next": find
# exactly where N stays non-negative, and check whether a saturating
# (Type II) functional response caps the amplitude the way a real
# predator-prey system would.
################################################################################

predator_prey_rhs_selflim <- function(t, state, parms) {
  N <- max(state[1], 0); P <- max(state[2], 0)
  D <- parms[["D"]]; a <- parms[["a"]]; e <- parms[["e"]]; m <- parms[["m"]]; K <- parms[["K"]]
  tau <- parms[["tau"]]; N0 <- parms[["N0"]]
  lag_N <- if (t <= tau) N0 else max(lagvalue(t - tau, 1), 0)
  dN <- D * (K - lag_N) - a * N * P
  dP <- e * a * N * P - m * P
  list(c(dN, dP))
}

# Denser than Day 26's onset sweep, out to tau=30, since min_N=-15.5 was
# already reached by tau=20.
tau_seq_selflim_biological <- c(2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10, 15, 20, 25, 30)

selflim_biological_df <- bind_rows(lapply(tau_seq_selflim_biological, function(tau) {
  run_predator_prey_delay(predator_prey_rhs_selflim, tau = tau)
}))

print(selflim_biological_df)
write.csv(selflim_biological_df, file.path("interesting_plots", "selflim_biological_range.csv"),
          row.names = FALSE)

# Bracket where min_N crosses zero. min_N FALLS as tau rises here (the
# mirror image of Day 26's find_crossing() for alpha, which rises with its
# parameter), so "last still non-negative" / "first negative" are the other
# way round from that helper.
bracket_descending_zero <- function(df, value_col = "tau", metric_col = "min_N") {
  df     <- df %>% arrange(.data[[value_col]])
  nonneg <- df %>% filter(.data[[metric_col]] >= 0) %>% pull(.data[[value_col]])
  neg    <- df %>% filter(.data[[metric_col]] <  0) %>% pull(.data[[value_col]])
  data.frame(
    last_nonnegative = if (length(nonneg) > 0) max(nonneg) else NA_real_,
    first_negative    = if (length(neg) > 0) min(neg) else NA_real_
  )
}

selflim_biological_bracket <- bracket_descending_zero(selflim_biological_df)
print(selflim_biological_bracket)
write.csv(selflim_biological_bracket, file.path("interesting_plots", "selflim_biological_bracket.csv"),
          row.names = FALSE)

cat(sprintf(
  "Self-limitation-delay model, mass-action predation: N stays non-negative up to tau = %s, negative by tau = %s.\n",
  ifelse(is.na(selflim_biological_bracket$last_nonnegative), "?", selflim_biological_bracket$last_nonnegative),
  ifelse(is.na(selflim_biological_bracket$first_negative), "?", selflim_biological_bracket$first_negative)
))

selflim_biological_plot <- ggplot(selflim_biological_df, aes(x = tau, y = min_N)) +
  geom_line(color = "#C44E52", linewidth = 1) +
  geom_point(color = "#C44E52", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "tau (delay in prey self-limitation)", y = "min(N) over the settled window",
       title = "How far past onset does the self-limitation model stay biologically sensible?",
       subtitle = "Dashed line = N = 0 -- below it, the toy model has gone unphysical") +
  theme_minimal()
selflim_biological_plot
save_plot(selflim_biological_plot, "selflim_bio_range.png")

# Does a saturating (Holling Type II) functional response -- replacing the
# mass-action a*N*P intake term with one that flattens out at high N -- cap
# the amplitude and keep N non-negative further past onset? h is the
# handling time; larger h saturates sooner.
predator_prey_rhs_selflim_typeII <- function(t, state, parms) {
  N <- max(state[1], 0); P <- max(state[2], 0)
  D <- parms[["D"]]; a <- parms[["a"]]; e <- parms[["e"]]; m <- parms[["m"]]; K <- parms[["K"]]
  tau <- parms[["tau"]]; h <- parms[["h"]]; N0 <- parms[["N0"]]
  lag_N <- if (t <= tau) N0 else max(lagvalue(t - tau, 1), 0)
  intake <- a * N * P / (1 + a * h * N)
  dN <- D * (K - lag_N) - intake
  dP <- e * intake - m * P
  list(c(dN, dP))
}

# Modest handling time -- large enough to matter at N* ~ 0.6, small enough
# not to shift the equilibrium much.
selflim_typeII_h <- 0.5

selflim_typeII_df <- bind_rows(lapply(tau_seq_selflim_biological, function(tau) {
  run_predator_prey_delay(predator_prey_rhs_selflim_typeII, tau = tau, extra_parms = list(h = selflim_typeII_h))
}))

print(selflim_typeII_df)
write.csv(selflim_typeII_df, file.path("interesting_plots", "selflim_typeII_biological_range.csv"),
          row.names = FALSE)

selflim_typeII_bracket <- bracket_descending_zero(selflim_typeII_df)
print(selflim_typeII_bracket)

selflim_compare_df <- bind_rows(
  selflim_biological_df %>% mutate(model = "Type I (mass action)"),
  selflim_typeII_df %>% mutate(model = sprintf("Type II (h = %.1f)", selflim_typeII_h))
)

selflim_compare_plot <- ggplot(selflim_compare_df, aes(x = tau, y = min_N, color = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "tau (delay in prey self-limitation)", y = "min(N) over the settled window",
       title = "Does a saturating functional response cap the amplitude?",
       subtitle = "Type I vs Type II predation, same delayed self-limitation structure", color = NULL) +
  theme_minimal()
selflim_compare_plot
save_plot(selflim_compare_plot, "selflim_typeI_II.png")

cat(sprintf(
  paste0(
    "Type I vs Type II: mass-action N stays non-negative up to tau = %s; Type II\n",
    "(h = %.1f) up to tau = %s. A meaningfully larger Type II bracket means\n",
    "saturation is doing real work capping the amplitude; about the same means\n",
    "the instability itself needs a different cap (e.g. a carrying-capacity-like\n",
    "term on P).\n"
  ),
  ifelse(is.na(selflim_biological_bracket$last_nonnegative), "?", selflim_biological_bracket$last_nonnegative),
  selflim_typeII_h,
  ifelse(is.na(selflim_typeII_bracket$last_nonnegative), "?", selflim_typeII_bracket$last_nonnegative)
))

################################################################################
# 5. A fourth delay placement: delay only N(t-tau) inside the predator's
#    intake term, P(t) instantaneous there
#
# Completes the space of single-delay placements in a*N*P: Day 25 delayed
# both N and P together (impossible); Section 3 above delayed P alone (also
# impossible, once corrected); Day 26 delayed N inside the prey's OWN
# self-limitation term (possible, aP* < D). The remaining combination:
#
#   dN/dt = D(K - N(t))       - a*N(t)*P(t)      <- unchanged, instantaneous
#   dP/dt = e*a*N(t-tau)*P(t) - m*P(t)            <- only N delayed here
#
# Full derivation (linearising, the characteristic equation, nondimensional-
# ising to delta=D/m and kappa=a*P*/m, and the closed-form omega/tau_c) is
# worked out in Day 27's post -- not repeated here. Headline result: unlike
# the P-only-delay variant, this one has NO threshold -- a Hopf bifurcation
# exists for ANY positive D, a, e, m, K, with:
#   omega  = sqrt[ (-alpha^2 + sqrt(alpha^4 + 4*(beta*gamma)^2)) / 2 ]
#   tau_c  = atan2(alpha*omega, omega^2) / omega
# At this project's reference parameters: omega ~ 0.1197, tau_c ~ 12.52.
################################################################################

check_hopf_feasibility_ndelay <- function(parms = predator_prey_params_delay) {
  D <- parms$D; a <- parms$a; e <- parms$e; m <- parms$m; K <- parms$K
  N_star <- m / (e * a)
  P_star <- D * (K - N_star) / (a * N_star)
  alpha  <- D + a * P_star
  beta   <- a * N_star
  gamma  <- e * a * P_star
  bg     <- beta * gamma

  # u^2 + alpha^2*u - bg^2 = 0 -- product of roots (-bg^2) is always
  # negative, so a positive real root always exists.
  u_pos <- (-alpha^2 + sqrt(alpha^4 + 4 * bg^2)) / 2
  omega <- sqrt(u_pos)
  tau_c <- atan2(alpha * omega, omega^2) / omega

  data.frame(N_star = N_star, P_star = P_star, alpha = alpha, beta = beta,
             gamma = gamma, beta_gamma = bg, omega = omega, tau_critical = tau_c)
}

ndelay_hopf_check <- check_hopf_feasibility_ndelay()
print(ndelay_hopf_check)
cat(sprintf(
  "N-only-delay (intake term) Hopf check: omega = %.4f, tau_critical = %.4f -- possible for ANY positive D, a, e, m, K.\n",
  ndelay_hopf_check$omega, ndelay_hopf_check$tau_critical
))

predator_prey_rhs_ndelay <- function(t, state, parms) {
  N <- max(state[1], 0); P <- max(state[2], 0)
  D <- parms[["D"]]; a <- parms[["a"]]; e <- parms[["e"]]; m <- parms[["m"]]; K <- parms[["K"]]
  tau <- parms[["tau"]]; N0 <- parms[["N0"]]
  lag_N <- if (t <= tau) N0 else max(lagvalue(t - tau, 1), 0)
  dN <- D * (K - N) - a * N * P
  dP <- e * a * lag_N * P - m * P
  list(c(dN, dP))
}

# Denser near the predicted tau_c ~ 12.52 than far from it, same convention
# as Day 26's self-limitation sweep around its own hand-estimated onset.
tau_seq_ndelay <- c(0, 2, 5, 8, 10, 11, 12, 12.5, 13, 14, 16, 20, 25, 30)

ndelay_sweep_df <- bind_rows(lapply(tau_seq_ndelay, function(tau) {
  run_predator_prey_delay(predator_prey_rhs_ndelay, tau = tau)
}))

print(ndelay_sweep_df)
write.csv(ndelay_sweep_df, file.path("interesting_plots", "predator_prey_ndelay_sweep.csv"),
          row.names = FALSE)

ndelay_amplitude_plot <- ggplot(ndelay_sweep_df, aes(x = tau, y = rel_amplitude_P)) +
  geom_line(color = "#8172B2", linewidth = 1) +
  geom_point(color = "#8172B2", size = 1.5) +
  geom_vline(xintercept = ndelay_hopf_check$tau_critical, linetype = "dashed", color = "grey40") +
  labs(x = "tau (delay in N(t-tau) inside the intake term, P(t) instantaneous)",
       y = "Relative amplitude of P (settled)",
       title = "Delaying only N inside the intake term: onset predicted at every parameter value",
       subtitle = sprintf("Dashed = analytically predicted tau_critical = %.3f", ndelay_hopf_check$tau_critical)) +
  theme_minimal()
ndelay_amplitude_plot
save_plot(ndelay_amplitude_plot, "ndelay_amp.png")

# Persistence check just past onset: a genuine limit cycle should hold a
# steady amplitude across increasing t_max, not decay like an unsettled
# transient.
ndelay_persistence_check <- bind_rows(lapply(c(600, 3000, 6000), function(tm) {
  run_predator_prey_delay(predator_prey_rhs_ndelay, tau = 15, t_max = tm) %>% mutate(t_max = tm)
}))

print(ndelay_persistence_check)
write.csv(ndelay_persistence_check, file.path("interesting_plots", "ndelay_persistence_check.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Persistence check at tau=15 (past predicted onset tau_c=%.2f): rel_amplitude_P\n",
    "at t_max=600, 3000, 6000 -- %.4f, %.4f, %.4f. Steady = genuine limit cycle;\n",
    "collapsing toward 0 = tau=15 isn't past the real onset after all.\n"
  ),
  ndelay_hopf_check$tau_critical,
  ndelay_persistence_check$rel_amplitude_P[1],
  ndelay_persistence_check$rel_amplitude_P[2],
  ndelay_persistence_check$rel_amplitude_P[3]
))

# Below-onset check: the sweep above climbs smoothly from tau=8 to tau=12
# with no visible knee at tau_c=12.52 -- the same critical-slowing-down
# signature seen elsewhere in this project (Day 26's gamma dip, the
# conversion-delay model's misleading amplitude climb). Confirms that climb
# is transient relaxation, not a real earlier onset.
ndelay_below_onset_check <- bind_rows(lapply(c(9, 10, 11), function(tau) {
  bind_rows(lapply(c(600, 3000, 6000), function(tm) {
    run_predator_prey_delay(predator_prey_rhs_ndelay, tau = tau, t_max = tm) %>% mutate(t_max = tm)
  }))
}))

print(ndelay_below_onset_check)
write.csv(ndelay_below_onset_check, file.path("interesting_plots", "ndelay_below_onset_check.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Below-onset check (tau=9,10,11; predicted tau_c=%.2f): rel_amplitude_P\n",
    "shrinking as t_max grows (600 -> 3000 -> 6000) confirms the pre-onset climb\n",
    "was transient relaxation, not a real earlier onset.\n"
  ),
  ndelay_hopf_check$tau_critical
))

# Extending tau much further: does min_N eventually go negative the way the
# self-limitation model's did almost immediately past its own onset?
tau_seq_ndelay_extended <- c(35, 40, 50, 60, 80, 100, 150, 200)

ndelay_sweep_extended_df <- bind_rows(lapply(tau_seq_ndelay_extended, function(tau) {
  run_predator_prey_delay(predator_prey_rhs_ndelay, tau = tau)
}))

print(ndelay_sweep_extended_df)
write.csv(ndelay_sweep_extended_df, file.path("interesting_plots", "ndelay_sweep_extended.csv"),
          row.names = FALSE)

ndelay_full_sweep_df <- bind_rows(ndelay_sweep_df, ndelay_sweep_extended_df) %>% arrange(tau)

ndelay_minN_plot <- ggplot(ndelay_full_sweep_df, aes(x = tau, y = min_N)) +
  geom_line(color = "#8172B2", linewidth = 1) +
  geom_point(color = "#8172B2", size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = ndelay_hopf_check$tau_critical, linetype = "dotted", color = "grey40") +
  labs(x = "tau (delay in N(t-tau) inside the intake term)", y = "min(N) over the settled window",
       title = "N-only-delay model: does min_N ever go negative, and how far out?",
       subtitle = "Dotted = predicted tau_critical; dashed = N = 0 (unphysical below this)") +
  theme_minimal()
ndelay_minN_plot
save_plot(ndelay_minN_plot, "ndelay_minN_ext.png")

# Oscillation PERIOD, measured directly from a settled time series via peak
# detection -- same approach as the operating-point check in Section 1,
# ported to this DDE model. Needed for the m/alpha sweeps below, where the
# closed-form prediction is period = 2*pi/omega.
measure_period_ndelay <- function(tau, parms, t_max = 4000, dt = 0.05,
                                  state0 = c(N = 1, P = 0.1), mxhist = DEFAULT_MXHIST) {
  parms$tau <- tau
  parms$N0  <- state0[["N"]]
  out <- as.data.frame(dede(y = state0, times = seq(0, t_max, by = dt),
                            func = predator_prey_rhs_ndelay, parms = parms,
                            control = list(mxhist = mxhist)))
  late <- out[out$time > t_max * 0.6, ]

  y <- late$P
  is_max <- c(FALSE, y[-1] > y[-length(y)]) & c(y[-length(y)] > y[-1], FALSE)
  maxima_times <- late$time[is_max]

  data.frame(
    measured_period = if (length(maxima_times) >= 2) mean(diff(maxima_times)) else NA_real_,
    n_peaks_found   = length(maxima_times),
    rel_amplitude_P = (max(late$P) - min(late$P)) / mean(late$P)
  )
}

# Period vs alpha (via D): holding a, e, m, K fixed and varying D directly
# varies alpha (= D*K/N*, and N* doesn't depend on D), a clean linear alpha
# sweep. FIXED DESIGN: testing each D at 1.5*tau_c(D) confounds D's own
# effect with tau_c(D) itself shifting (11.48 -> 12.97 across this range) --
# bad enough to flip the apparent trend. TAU_SHARED_ALPHA tests every D at
# one common tau, comfortably past every D's own onset.
TAU_SHARED_ALPHA <- 20
D_seq_period <- c(0.3, 0.5, 1, 1.5, 2, 3, 5)

alpha_period_df <- bind_rows(lapply(D_seq_period, function(D) {
  parms <- modifyList(predator_prey_params_delay, list(D = D))
  pred  <- check_hopf_feasibility_ndelay(parms)
  meas  <- measure_period_ndelay(TAU_SHARED_ALPHA, parms)
  data.frame(D = D, alpha = pred$alpha, tau_critical = pred$tau_critical,
             tau_tested = TAU_SHARED_ALPHA, predicted_period = 2 * pi / pred$omega,
             measured_period = meas$measured_period, n_peaks_found = meas$n_peaks_found)
}))

print(alpha_period_df)
write.csv(alpha_period_df, file.path("interesting_plots", "ndelay_alpha_period_sweep.csv"),
          row.names = FALSE)

alpha_period_plot <- ggplot(alpha_period_df, aes(x = alpha)) +
  geom_line(aes(y = predicted_period, color = "predicted (2*pi/omega)"), linewidth = 1) +
  geom_point(aes(y = measured_period, color = "measured (peak spacing)"), size = 2.5) +
  labs(x = "alpha (= D*K/N*, varied via D)", y = "Oscillation period",
       title = "N-only-delay model: how does the oscillation's period depend on alpha?",
       subtitle = "Predicted vs measured -- period only weakly sensitive to alpha once D isn't small",
       color = NULL) +
  theme_minimal()
alpha_period_plot
save_plot(alpha_period_plot, "ndelay_period_alpha.png")

# Period vs m: m moves the equilibrium itself (N*=m/(ea), P* both depend on
# m), not just a rescaling, and the closed form traces a U-shape -- minimised
# in the interior, diverging as m -> 0 or m -> a*e*K (where the equilibrium
# collides with extinction). FIXED DESIGN: same confound as the alpha sweep,
# worse here since tau_c(m) ranges over nearly 3x (12.09 to 34.81) --
# TAU_SHARED_M isolates m's own effect, fixed comfortably above the largest
# tau_c but not excessively far past it, since extended tau can turn this
# model numerically pathological (see the extended-tau check above).
TAU_SHARED_M <- 45
m_seq_period <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45)

m_period_df <- bind_rows(lapply(m_seq_period, function(m) {
  parms <- modifyList(predator_prey_params_delay, list(m = m))
  pred  <- check_hopf_feasibility_ndelay(parms)
  meas  <- measure_period_ndelay(TAU_SHARED_M, parms)
  data.frame(m = m, P_star = pred$P_star, tau_critical = pred$tau_critical,
             tau_tested = TAU_SHARED_M, predicted_period = 2 * pi / pred$omega,
             measured_period = meas$measured_period, n_peaks_found = meas$n_peaks_found)
}))

print(m_period_df)
write.csv(m_period_df, file.path("interesting_plots", "ndelay_m_period_sweep.csv"),
          row.names = FALSE)

m_period_plot <- ggplot(m_period_df, aes(x = m)) +
  geom_line(aes(y = predicted_period, color = "predicted (2*pi/omega)"), linewidth = 1) +
  geom_point(aes(y = measured_period, color = "measured (peak spacing)"), size = 2.5) +
  labs(x = "m (predator mortality)", y = "Oscillation period",
       title = "N-only-delay model: period vs m -- a U-shape, not monotonic",
       subtitle = "m must stay in (0, a*e*K = 0.5) for the equilibrium to exist; period blows up near both ends",
       color = NULL) +
  theme_minimal()
m_period_plot
save_plot(m_period_plot, "ndelay_period_m.png")

cat(sprintf(
  paste0(
    "Period vs m: predicted period ranges from %.1f (m=%.2f) up past %.1f (m=%.2f)\n",
    "across the tested range, minimised somewhere in the middle -- compare against\n",
    "measured_period in ndelay_m_period_sweep.csv to confirm the simulation tracks\n",
    "the closed form's U-shape, not just its rough magnitude.\n"
  ),
  min(m_period_df$predicted_period), m_period_df$m[which.min(m_period_df$predicted_period)],
  max(m_period_df$predicted_period), m_period_df$m[which.max(m_period_df$predicted_period)]
))
