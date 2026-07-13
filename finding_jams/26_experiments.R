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

dir.create("interesting_plots", showWarnings = FALSE)

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 150) {
  ggsave(file.path("interesting_plots", filename), plot = plot,
         width = width, height = height, dpi = dpi)
}

# Self-contained, same as every script since Day 20 -- redefined here rather
# than sourced from 25_experiments.R.
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

# Same catchable_fraction_at()/run_param_sweep() as 25_experiments.R, carried
# over unchanged so today's refinements are directly comparable to the sweep
# that found the transition and the dip in the first place.
catchable_fraction_at <- function(params, t_max = 600) {
  sim <- project(params, t_max = t_max, dt = 0.1, t_save = 0.5,
                 progress_bar = FALSE, effort = 0, method = "tr_bdf2")
  last  <- dim(sim@n)[1]
  w_mat <- params@species_params$w_mat[1]
  w     <- sim@params@w
  dw    <- sim@params@dw

  bm_density   <- sim@n[last, 1, ] * w * dw
  total_bm     <- sum(bm_density)
  catchable_bm <- sum(bm_density[w >= w_mat])

  data.frame(catchable_fraction = catchable_bm / total_bm, error = NA_character_)
}

run_param_sweep <- function(param_name, param_seq, capacity_mult = 10,
                            resource_decrease = 0.001, t_max = 600) {
  bind_rows(lapply(param_seq, function(v) {
    args <- list(capacity_mult = capacity_mult, resource_decrease = resource_decrease)
    args[[param_name]] <- v
    result <- tryCatch({
      p <- do.call(make_second_order_params_kr, args)
      catchable_fraction_at(p, t_max = t_max)
    }, error = function(e) {
      data.frame(catchable_fraction = NA_real_, error = conditionMessage(e))
    })
    if (!is.na(result$error)) {
      warning(sprintf("%s = %.4g: %s", param_name, v, result$error))
    }
    data.frame(param = param_name, value = v,
               catchable_fraction = result$catchable_fraction,
               error = result$error)
  }))
}

################################################################################
# 1. Pinning down alpha's transition window
#
# Day 25's 15-point alpha sweep (0.02-0.8) found the catchable fraction
# jumps roughly 8 orders of magnitude between alpha = 0.057 (1.1e-9) and
# alpha = 0.28 (25.5%) -- a genuine viability threshold, not a gradual
# slope. That grid was never dense enough to say exactly where the steep
# part starts and ends, only that it's somewhere in that range. Same move
# as Day 24's capacity_mult refinement: rerun at much higher density,
# zoomed into just the window where the climb happens, at the same t_run
# (600, catchable_fraction_at()'s own default -- not overridden here, same
# as Day 25's original sweep) so the two are directly comparable.
#
# capacity_mult = 10, resource_decrease = 0.001, gamma/ks/kappa pinned at
# the anchovy's own values throughout -- the same operating point every
# isolation check since Day 24 has used.
#
# Unlike capacity_mult's collapse bifurcation, there's no natural
# alive/dead split here -- catchable_fraction rises smoothly (if steeply)
# from near-zero to a real plateau, not a bistable jump at one point. So
# the "transition window" is defined the same way Day 24 pinned down the
# capacity_mult edges, but with two threshold crossings instead of one:
# where the climb clearly starts (still under LOW_THRESHOLD) and where it's
# clearly finished (already over HIGH_THRESHOLD). Both thresholds sit
# comfortably inside [0.05, 0.3] going by Day 25's coarse table, so the
# refined grid below should bracket both edges without needing to guess.
################################################################################

alpha_seq_refine <- exp(seq(log(0.05), log(0.3), length.out = 25))

alpha_refine_df <- run_param_sweep("alpha", alpha_seq_refine)

saveRDS(alpha_refine_df, file.path("interesting_plots", "alpha_refine_df.rds"))
write.csv(alpha_refine_df, file.path("interesting_plots", "alpha_refine_df.csv"), row.names = FALSE)

# below/above the given threshold, walking the swept value ascending. The
# true crossing sits somewhere in (below_last, above_first]. Same shape as
# Day 24's find_transition(), generalised to an arbitrary threshold rather
# than a fixed alive/dead split, since alpha's climb is monotonic rather
# than bistable.
find_crossing <- function(df, threshold, value_col = "value",
                          metric_col = "catchable_fraction") {
  df    <- df %>% arrange(.data[[value_col]])
  below <- df %>% filter(.data[[metric_col]] < threshold)  %>% pull(.data[[value_col]])
  above <- df %>% filter(.data[[metric_col]] >= threshold) %>% pull(.data[[value_col]])
  data.frame(
    threshold   = threshold,
    below_last  = if (length(below) > 0) max(below) else NA_real_,
    above_first = if (length(above) > 0) min(above) else NA_real_
  )
}

LOW_THRESHOLD  <- 1e-3  # 0.1% -- still clearly pre-transition
HIGH_THRESHOLD <- 0.10  # 10%  -- clearly past the steep part

alpha_crossings <- bind_rows(
  find_crossing(alpha_refine_df %>% filter(is.na(error)), LOW_THRESHOLD)  %>% mutate(edge = "climb starts"),
  find_crossing(alpha_refine_df %>% filter(is.na(error)), HIGH_THRESHOLD) %>% mutate(edge = "climb ends")
)

print(alpha_crossings)
write.csv(alpha_crossings, file.path("interesting_plots", "alpha_transition_window.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "alpha transition window (capacity_mult = 10, resource_decrease = 0.001, t_run = 600):\n",
    "  Climb starts (crosses %.3g%%) somewhere in (%.4f, %.4f]\n",
    "  Climb ends   (crosses %.3g%%) somewhere in (%.4f, %.4f]\n"
  ),
  LOW_THRESHOLD * 100,
  alpha_crossings$below_last[alpha_crossings$edge == "climb starts"],
  alpha_crossings$above_first[alpha_crossings$edge == "climb starts"],
  HIGH_THRESHOLD * 100,
  alpha_crossings$below_last[alpha_crossings$edge == "climb ends"],
  alpha_crossings$above_first[alpha_crossings$edge == "climb ends"]
))

alpha_refine_plot <- ggplot(alpha_refine_df %>% filter(is.na(error)),
                            aes(x = value, y = catchable_fraction)) +
  geom_line(color = "#4C72B0", linewidth = 1) +
  geom_point(color = "#4C72B0", size = 1.5) +
  geom_hline(yintercept = c(LOW_THRESHOLD, HIGH_THRESHOLD),
            linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey40") +
  scale_y_log10(labels = scales::label_scientific()) +
  scale_x_log10() +
  labs(x = "alpha (assimilation efficiency)",
       y = "Catchable fraction (biomass at or above w_mat)",
       title = "Refined alpha sweep: pinning down the transition window",
       subtitle = "25 points over [0.05, 0.3] -- dotted lines = the two threshold crossings bracketed above; dashed = anchovy's own alpha = 0.1") +
  theme_minimal()
alpha_refine_plot

save_plot(alpha_refine_plot, "Alpha transition window refined.png")

################################################################################
# 2. Is the gamma = 305.8 dip real, or an under-settled transient?
#
# Day 25's 15-point gamma sweep found one point that breaks monotonicity:
# gamma = 305.8 read 9.48e-5, dipping below both its neighbours (231 ->
# 2.45e-4, 404 -> 4.03e-4), at t_run = 600. Flagged at the time as possibly
# critical slowing down near a smaller feature -- the same explanation Day
# 24 used for its own default-species Forward-branch spike, where
# t_run = 600 wasn't necessarily long enough to fully settle right at a
# transition. This reruns exactly that point and its immediate neighbours
# (located dynamically off the original 15-point grid, not hardcoded, so
# this stays correct if the grid ever changes) at increasingly long t_run
# to see whether the dip shrinks away (transient) or holds up (real).
################################################################################

gamma_full_seq <- exp(seq(log(100), log(5000), length.out = 15))  # Day 25's exact grid
dip_idx        <- which.min(abs(gamma_full_seq - 305.8))
neighbour_idx  <- unique(pmax(1, pmin(length(gamma_full_seq), c(dip_idx - 1, dip_idx, dip_idx + 1))))
gamma_check_vals <- gamma_full_seq[neighbour_idx]

cat(sprintf("gamma points being rechecked: %s\n",
           paste(sprintf("%.1f", gamma_check_vals), collapse = ", ")))

# Doubling from the sweep's own default (600) rather than an arbitrary
# jump, so the point where any convergence happens is itself informative.
t_run_check_seq <- c(600, 1200, 2400, 4800)

gamma_dip_check_df <- bind_rows(lapply(t_run_check_seq, function(tr) {
  bind_rows(lapply(gamma_check_vals, function(g) {
    p   <- make_second_order_params_kr(gamma = g, capacity_mult = 10,
                                       resource_decrease = 0.001)
    res <- catchable_fraction_at(p, t_max = tr)
    data.frame(gamma = g, t_run = tr, catchable_fraction = res$catchable_fraction)
  }))
}))

print(gamma_dip_check_df)
write.csv(gamma_dip_check_df, file.path("interesting_plots", "gamma_dip_check.csv"),
          row.names = FALSE)

# At each t_run: is the middle point still below both neighbours? If that
# stops being true as t_run grows, the dip was a settling artefact: if it's
# still true at t_run = 4800, the dip is a real, if small, feature.
dip_persistence <- gamma_dip_check_df %>%
  group_by(t_run) %>%
  arrange(gamma, .by_group = TRUE) %>%
  summarise(
    dip_still_present = catchable_fraction[2] < catchable_fraction[1] &
                        catchable_fraction[2] < catchable_fraction[3],
    .groups = "drop"
  )

print(dip_persistence)
write.csv(dip_persistence, file.path("interesting_plots", "gamma_dip_persistence.csv"),
          row.names = FALSE)

gamma_dip_plot <- ggplot(gamma_dip_check_df,
                         aes(x = t_run, y = catchable_fraction,
                             color = factor(round(gamma, 1)))) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_scientific()) +
  labs(x = "t_run", y = "Catchable fraction",
       title = "Is the gamma = 305.8 dip real or an under-settled transient?",
       subtitle = "Same point and its two neighbours, rerun at increasing t_run",
       color = "gamma") +
  theme_minimal()
gamma_dip_plot

save_plot(gamma_dip_plot, "Gamma dip check - t_run convergence.png")

cat(sprintf(
  paste0(
    "gamma = %.1f dip check: still below both neighbours at t_run = 600? %s. ",
    "Still below both neighbours at t_run = %d (longest tested)? %s.\n",
    "If that flips to FALSE before the longest t_run, the dip was critical ",
    "slowing down, not a real feature; if it's still TRUE, treat it as real ",
    "but small.\n"
  ),
  gamma_check_vals[2],
  dip_persistence$dip_still_present[dip_persistence$t_run == min(t_run_check_seq)],
  max(t_run_check_seq),
  dip_persistence$dip_still_present[dip_persistence$t_run == max(t_run_check_seq)]
))

################################################################################
# 3a. Reproducing Day 25's time-delay predator-prey model, unchanged
#
# This is the exact model, functions, and sweep from 25_experiments.R,
# copied here verbatim (not re-derived) so it can be run and checked
# standalone in this file, side by side with the follow-up in 3b below.
# Running this block alone should reproduce Day 25's reported numbers:
# hopf_possible = FALSE, rel_amplitude_P climbing from ~1.8e-10 at tau = 0
# to ~0.25 at tau = 40, and a monotonic (non-oscillating) staircase in the
# tau = 40 time series.
#
#   dN/dt = D*(K - N(t)) - a*N(t)*P(t)
#   dP/dt = e*a*N(t-tau)*P(t-tau) - m*P(t)
#
# Predation and prey loss are instantaneous; only the resulting predator
# biomass gain is delayed by tau (a maturation/gestation lag), following
# the classical Wangersky-Cunningham (1957) delayed-conversion form.
# Predator mortality stays instantaneous.
################################################################################

library(deSolve)

# D: prey renewal rate; a: attack rate; e: predator conversion efficiency;
# m: predator mortality; K: prey carrying capacity -- same as the Type I/II
# toy models in 24_experiments.R.
predator_prey_params_delay <- list(D = 1, a = 1, e = 0.5, m = 0.3, K = 1)

# Numeric check of the analytic Hopf-feasibility condition beta*gamma >
# 2*m*alpha derived in the write-up: if this ever prints hopf_possible =
# TRUE, the algebra is wrong and the simulation below becomes the more
# trustworthy half of this section, not the other way around.
check_hopf_feasibility <- function(parms = predator_prey_params_delay) {
  D <- parms$D; a <- parms$a; e <- parms$e; m <- parms$m; K <- parms$K
  N_star <- m / (e * a)
  P_star <- D * (K - N_star) / (a * N_star)
  alpha  <- D + a * P_star
  beta   <- a * N_star
  gamma  <- e * a * P_star
  data.frame(N_star = N_star, P_star = P_star, alpha = alpha, beta = beta,
             gamma = gamma, beta_gamma = beta * gamma, two_m_alpha = 2 * m * alpha,
             hopf_possible = (beta * gamma) > (2 * m * alpha))
}

hopf_check <- check_hopf_feasibility()
print(hopf_check)
cat(sprintf(
  paste0(
    "Hopf feasibility check: beta*gamma = %.4f, 2*m*alpha = %.4f -- %s\n",
    "(algebra predicts this is ALWAYS false, for any positive parameters)\n"
  ),
  hopf_check$beta_gamma, hopf_check$two_m_alpha,
  ifelse(hopf_check$hopf_possible, "HOPF POSSIBLE (unexpected!)", "no Hopf possible, as predicted")
))

# Delay differential equation via deSolve::dede(), using lagvalue() for the
# lagged state. For t <= tau there's no history yet, so the lag falls back
# to the initial condition (constant-history assumption). N, P clamped to
# >= 0 since a transiently negative lagged value would flip the sign of the
# a*N*P-style term.
predator_prey_rhs_delay <- function(t, state, parms) {
  N <- max(state[1], 0)
  P <- max(state[2], 0)
  D   <- parms[["D"]];   a <- parms[["a"]]; e <- parms[["e"]]
  m   <- parms[["m"]];   K <- parms[["K"]]; tau <- parms[["tau"]]
  N0  <- parms[["N0"]]; P0 <- parms[["P0"]]

  if (t <= tau) {
    lag_N <- N0
    lag_P <- P0
  } else {
    lag   <- lagvalue(t - tau)
    lag_N <- max(lag[1], 0)
    lag_P <- max(lag[2], 0)
  }

  dN <- D * (K - N) - a * N * P
  dP <- e * a * lag_N * lag_P - m * P
  list(c(dN, dP))
}

# Original short-form runner, exactly as in 25_experiments.R: one fresh
# integration to t_max, amplitude read off the last 40% ("settled window").
run_predator_prey_delay <- function(tau, state0 = c(N = 1, P = 0.1), t_max = 300,
                                    parms = predator_prey_params_delay) {
  parms$tau <- tau
  parms$N0  <- state0[["N"]]
  parms$P0  <- state0[["P"]]
  out  <- as.data.frame(dede(y = state0, times = seq(0, t_max, by = 0.1),
                             func = predator_prey_rhs_delay, parms = parms))
  late <- out[out$time > t_max * 0.6, ]
  data.frame(tau = tau, max_N = max(late$N), min_N = min(late$N),
             max_P = max(late$P), min_P = min(late$P),
             rel_amplitude_P = (max(late$P) - min(late$P)) / mean(late$P))
}

# Perturbed well away from equilibrium (N=1 vs N*=0.6, P=0.1 vs P*=0.667) --
# a system started exactly at a stable fixed point never oscillates
# regardless of tau, so this is a genuine test of whether delay lets a real
# perturbation grow into sustained oscillation rather than damping out.
# tau up to 40, two orders of magnitude past the model's own natural
# timescales (1/D, 1/m).
tau_seq <- c(0, 0.5, 1, 2, 5, 10, 20, 40)

delay_sweep_df <- bind_rows(lapply(tau_seq, run_predator_prey_delay))

print(delay_sweep_df)
write.csv(delay_sweep_df, file.path("interesting_plots", "predator_prey_delay_sweep_check.csv"),
          row.names = FALSE)

delay_amplitude_plot <- ggplot(delay_sweep_df, aes(x = tau, y = rel_amplitude_P)) +
  geom_line(color = "#55A868", linewidth = 1) +
  geom_point(color = "#55A868", size = 1.5) +
  labs(x = "tau (delay in predator conversion)", y = "Relative amplitude of P (settled)",
       title = "Does delay reintroduce the oscillation the instantaneous model doesn't have?",
       subtitle = "Perturbed start (N=1, P=0.1), D=1, a=1, e=0.5, m=0.3, K=1 -- analytic prediction: no, for any tau") +
  theme_minimal()
delay_amplitude_plot

save_plot(delay_amplitude_plot, "Predator-prey toy model - delay sweep amplitude - check.png")

# Dedicated time-series run at the largest tested delay -- if delay were
# going to produce a growing oscillation instead of damping to the fixed
# point, it should be visible here, not just in the settled-window amplitude
# summary above.
longest_delay_ts <- as.data.frame(dede(
  y = c(N = 1, P = 0.1),
  times = seq(0, 300, by = 0.1),
  func = predator_prey_rhs_delay,
  parms = modifyList(predator_prey_params_delay, list(tau = 40, N0 = 1, P0 = 0.1))
))

delay_timeseries_plot <- ggplot(longest_delay_ts, aes(x = time)) +
  geom_line(aes(y = N, color = "N (prey)")) +
  geom_line(aes(y = P, color = "P (predator)")) +
  geom_hline(yintercept = hopf_check$N_star, linetype = "dotted", color = "grey40") +
  geom_hline(yintercept = hopf_check$P_star, linetype = "dotted", color = "grey40") +
  labs(x = "Time", y = "Biomass", color = NULL,
       title = "Time series at tau = 40 (the longest delay tested)",
       subtitle = "Dotted lines = analytic equilibrium N*, P* -- both unchanged by tau, as derived above") +
  theme_minimal()
delay_timeseries_plot

save_plot(delay_timeseries_plot, "Predator-prey toy model - delay timeseries tau40 - check.png")

cat(paste0(
  "Time-delay predator-prey check: if rel_amplitude_P stays small and flat across\n",
  "every tau tested, and the tau=40 time series damps back to N*, P* rather than\n",
  "sustaining or growing an oscillation, that confirms the analytic prediction --\n",
  "a Hopf bifurcation is impossible for this delayed-conversion model regardless of\n",
  "tau, because beta*gamma - 2*m*alpha = -a*m*P* - 2*m*D is always negative for\n",
  "positive parameters. Chemostat-style prey growth blocks this route to\n",
  "hysteresis too, not just the Type I/II functional-response route Day 24 ruled out.\n"
))

################################################################################
# 3b. Does rel_amplitude_P actually decay to ~0 at large tau, or plateau?
#
# The sweep above (3a, reproducing Day 25) shows rel_amplitude_P climbing
# with tau -- up to ~0.25 at tau = 40 -- even though the algebra says a
# Hopf bifurcation is impossible and the tau = 40 time series is
# monotonic. The claim is that this is critical-slowing-down, not real
# growth: larger tau relaxes more slowly, so a fixed t_max = 300 catches
# the system still mid-approach. That predicts rel_amplitude_P should
# shrink toward 0 as t_max grows -- worth checking directly rather than
# assuming, exactly as Day 25's own "What's Next" flagged.
#
# Rather than rerunning tau = 40/80 separately at every t_max (which would
# just re-integrate the same trajectory's early portion over and over),
# each tau is integrated ONCE out to the longest t_max tested, and
# rel_amplitude_P is recomputed by restricting that single trajectory to
# progressively later cutoffs. That's mathematically equivalent to
# rerunning at each shorter t_max -- this is a deterministic DDE with a
# fixed history, so the trajectory up to any given time doesn't depend on
# how much further it's eventually integrated -- and it's far cheaper than
# 2 taus x 6 t_max values worth of independent long integrations. Reuses
# predator_prey_rhs_delay(), predator_prey_params_delay, and hopf_check's
# N_star/P_star from 3a above rather than redefining them.
################################################################################

run_predator_prey_delay_long <- function(tau, t_max_long, state0 = c(N = 1, P = 0.1),
                                         parms = predator_prey_params_delay) {
  parms$tau <- tau
  parms$N0  <- state0[["N"]]
  parms$P0  <- state0[["P"]]
  as.data.frame(dede(y = state0, times = seq(0, t_max_long, by = 0.1),
                     func = predator_prey_rhs_delay, parms = parms))
}

# Same "late window = last 40%" convention as run_predator_prey_delay() in
# 25_experiments.R, just applied to a slice of an already-computed long run
# instead of a fresh integration.
amplitude_at_cutoff <- function(ts_df, t_max_effective) {
  late <- ts_df[ts_df$time > t_max_effective * 0.6 & ts_df$time <= t_max_effective, ]
  data.frame(t_max = t_max_effective,
             rel_amplitude_P = (max(late$P) - min(late$P)) / mean(late$P))
}

tau_check_seq   <- c(40, 80)  # "and perhaps larger" -- one further past 40
t_max_check_seq <- c(300, 600, 1200, 2400, 4800, 9600)

# One long integration per tau, out to the longest t_max tested -- reused
# below for both the windowed-amplitude check and the endpoint-distance
# check, rather than integrating twice per tau for two different metrics.
ts_long_by_tau <- setNames(
  lapply(tau_check_seq, function(tau) run_predator_prey_delay_long(tau, max(t_max_check_seq))),
  tau_check_seq
)

decay_check_df <- bind_rows(lapply(tau_check_seq, function(tau) {
  ts_long <- ts_long_by_tau[[as.character(tau)]]
  bind_rows(lapply(t_max_check_seq, function(tm) {
    amplitude_at_cutoff(ts_long, tm) %>% mutate(tau = tau)
  }))
}))

print(decay_check_df)
write.csv(decay_check_df, file.path("interesting_plots", "delay_amplitude_decay_check.csv"),
          row.names = FALSE)

decay_check_plot <- ggplot(decay_check_df, aes(x = t_max, y = rel_amplitude_P,
                                               color = factor(tau))) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_scientific()) +
  labs(x = "t_max (how long the run is allowed to settle)",
       y = "Relative amplitude of P, late window",
       title = "Does rel_amplitude_P decay to ~0 given enough time, or plateau?",
       subtitle = "Same tau = 40 run as Day 25, plus tau = 80, re-windowed at longer cutoffs -- analytic prediction: decays to 0",
       color = "tau") +
  theme_minimal()
decay_check_plot

save_plot(decay_check_plot, "Predator-prey toy model - amplitude decay check.png")

# Direct endpoint check, not just the windowed amplitude metric: how close
# is the trajectory to the analytic equilibrium at the very end of the
# longest run tested, for each tau. If this keeps shrinking as t_max grows
# rather than levelling off at some nonzero floor, that's a second,
# independent confirmation of decay toward the fixed point rather than a
# sustained oscillation.
endpoint_check_df <- bind_rows(lapply(tau_check_seq, function(tau) {
  ts_long   <- ts_long_by_tau[[as.character(tau)]]
  final_row <- tail(ts_long, 1)
  data.frame(tau = tau, t_max = max(t_max_check_seq),
             N_end = final_row$N, P_end = final_row$P,
             dist_to_equilibrium = sqrt((final_row$N - hopf_check$N_star)^2 +
                                        (final_row$P - hopf_check$P_star)^2))
}))

print(endpoint_check_df)
write.csv(endpoint_check_df, file.path("interesting_plots", "delay_endpoint_check.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Amplitude decay check: if rel_amplitude_P in delay_amplitude_decay_check.csv\n",
    "keeps shrinking as t_max grows (rather than levelling off at a nonzero value)\n",
    "for both tau = 40 and tau = 80, and dist_to_equilibrium in\n",
    "delay_endpoint_check.csv is very small by t_max = %d, that confirms the\n",
    "algebra's prediction directly: the earlier climbing amplitude sweep was\n",
    "measuring an unfinished monotonic relaxation, not a real or growing\n",
    "oscillation, and rel_amplitude_P genuinely decays to ~0 given enough time.\n"
  ),
  max(t_max_check_seq)
))

################################################################################
# 4. A model that DOES support a Hopf bifurcation: delay the prey's own
#    self-limitation instead of the predator's conversion
#
# Section 3 proved that with the delay on the predator's conversion term
# (dP/dt = e*a*N(t-tau)*P(t-tau) - m*P(t)) and instantaneous chemostat prey
# growth, beta*gamma - 2*m*alpha = -a*m*P* - 2*m*D is always negative, for
# ANY positive D, a, e, m, K -- a Hopf bifurcation is structurally
# impossible for that model, not just for the specific numbers tested.
#
# Moving the delay onto the prey's OWN self-limitation term instead --
# dN/dt = D*(K - N(t-tau)) - a*N(t)*P(t), predation/conversion now
# instantaneous -- is the classical route to delay-induced instability
# (the two-species analogue of Hutchinson's single-species delayed
# logistic equation). Re-deriving via the same method as section 3:
#
# Equilibrium is UNCHANGED (N* = m/(e*a), P* = D*(K-N*)/(a*N*) -- same
# reasoning as before, the delayed term equals its instantaneous value at
# steady state).
#
# Linearising gives:
#   dn/dt = -alpha*n(t) - beta*p(t) - D*n(t-tau)   [alpha = a*P*, NOT D+a*P* --
#                                                     D now sits on its own,
#                                                     as the coefficient of
#                                                     the delayed term]
#   dp/dt = gamma*n(t)                              [gamma = e*a*P*, no delay
#                                                     in this equation at all]
#
# Characteristic equation: lambda^2 + alpha*lambda + beta*gamma +
#                           D*lambda*exp(-lambda*tau) = 0
# At tau=0: lambda^2 + (alpha+D)*lambda + beta*gamma = 0 -- alpha+D = a*P*+D
# is exactly the OLD alpha (D+a*P*), so this matches Day 24's non-delayed
# result exactly, as it must.
#
# Setting lambda=i*omega and matching moduli (same method as section 3)
# gives a quadratic in u=omega^2:
#   u^2 + (alpha^2 - D^2 - 2*beta*gamma)*u + (beta*gamma)^2 = 0
#
# The product of roots here is (beta*gamma)^2 -- always a perfect square,
# always >= 0 -- which rules out the "opposite signs" case from section 3
# entirely. Both roots always share the same sign; the sum decides which.
# Requiring both real (discriminant >= 0) AND positive (sum > 0) and
# simplifying collapses to a single clean condition:
#
#   Hopf possible  <=>  alpha < D  <=>  a*P* < D
#
# Checked against this project's own reference parameters (D=1, a=1,
# e=0.5, m=0.3, K=1 -> P*=0.667): a*P* = 0.667 < D = 1 -- the condition
# HOLDS, with the exact same numbers used throughout Days 25-26.
################################################################################

# Only N is delayed here (inside the self-limitation term); P depends only
# on instantaneous N(t), P(t) -- lagvalue(, 1) pulls just the N component
# rather than the whole lagged state vector.
predator_prey_rhs_selflim <- function(t, state, parms) {
  N <- max(state[1], 0)
  P <- max(state[2], 0)
  D   <- parms[["D"]];  a <- parms[["a"]]; e <- parms[["e"]]
  m   <- parms[["m"]];  K <- parms[["K"]]; tau <- parms[["tau"]]
  N0  <- parms[["N0"]]

  if (t <= tau) {
    lag_N <- N0
  } else {
    lag_N <- max(lagvalue(t - tau, 1), 0)
  }

  dN <- D * (K - lag_N) - a * N * P
  dP <- e * a * N * P - m * P
  list(c(dN, dP))
}

# Numeric check of the derived condition alpha < D (i.e. a*P* < D).
check_hopf_feasibility_selflim <- function(parms = predator_prey_params_delay) {
  D <- parms$D; a <- parms$a; e <- parms$e; m <- parms$m; K <- parms$K
  N_star  <- m / (e * a)
  P_star  <- D * (K - N_star) / (a * N_star)
  alpha_s <- a * P_star
  beta_s  <- a * N_star
  gamma_s <- e * a * P_star
  data.frame(N_star = N_star, P_star = P_star, alpha_s = alpha_s, D = D,
             beta_gamma = beta_s * gamma_s,
             hopf_possible = alpha_s < D)
}

hopf_check_selflim <- check_hopf_feasibility_selflim()
print(hopf_check_selflim)
cat(sprintf(
  paste0(
    "Self-limitation-delay Hopf check: a*P* = %.4f vs D = %.4f -- %s\n",
    "(derived condition: Hopf possible iff a*P* < D; contrast with section 3,\n",
    "where the analogous condition could never hold for any parameters)\n"
  ),
  hopf_check_selflim$alpha_s, hopf_check_selflim$D,
  ifelse(hopf_check_selflim$hopf_possible, "HOPF POSSIBLE", "not possible")
))

run_predator_prey_selflim <- function(tau, state0 = c(N = 1, P = 0.1), t_max = 300,
                                      parms = predator_prey_params_delay) {
  parms$tau <- tau
  parms$N0  <- state0[["N"]]
  out  <- as.data.frame(dede(y = state0, times = seq(0, t_max, by = 0.1),
                             func = predator_prey_rhs_selflim, parms = parms))
  late <- out[out$time > t_max * 0.6, ]
  data.frame(tau = tau, max_N = max(late$N), min_N = min(late$N),
             max_P = max(late$P), min_P = min(late$P),
             rel_amplitude_P = (max(late$P) - min(late$P)) / mean(late$P))
}

# Denser near the analytically-estimated critical tau (~2.4, worked out by
# hand from the phase-matching condition -- treat that as a rough guide for
# where to look, not a certified number; the simulation below is what
# actually pins down where the onset sits).
tau_seq_selflim <- c(0, 0.5, 1, 1.5, 2, 2.2, 2.4, 2.6, 2.8, 3, 4, 5, 7, 10, 15, 20)

selflim_sweep_df <- bind_rows(lapply(tau_seq_selflim, run_predator_prey_selflim))

print(selflim_sweep_df)
write.csv(selflim_sweep_df, file.path("interesting_plots", "predator_prey_selflim_delay_sweep.csv"),
          row.names = FALSE)

selflim_amplitude_plot <- ggplot(selflim_sweep_df, aes(x = tau, y = rel_amplitude_P)) +
  geom_line(color = "#C44E52", linewidth = 1) +
  geom_point(color = "#C44E52", size = 1.5) +
  labs(x = "tau (delay in prey self-limitation)", y = "Relative amplitude of P (settled)",
       title = "Delaying the prey's self-limitation: does a genuine Hopf bifurcation appear?",
       subtitle = "Same D=1,a=1,e=0.5,m=0.3,K=1 as section 3 -- analytic condition a*P*=0.667 < D=1 predicts yes") +
  theme_minimal()
selflim_amplitude_plot

save_plot(selflim_amplitude_plot, "Predator-prey toy model - selflim delay sweep amplitude.png")

################################################################################
# Confirming it's a genuine sustained oscillation, not another relaxation
# artefact: rerun a tau clearly above the apparent onset at a much longer
# t_max. If this is a real limit cycle, the amplitude should hold roughly
# steady rather than decaying away the way section 3's did.
################################################################################

selflim_persistence_check <- bind_rows(
  run_predator_prey_selflim(tau = 10, t_max = 300)  %>% mutate(t_max = 300),
  run_predator_prey_selflim(tau = 10, t_max = 1500) %>% mutate(t_max = 1500),
  run_predator_prey_selflim(tau = 10, t_max = 3000) %>% mutate(t_max = 3000)
)

print(selflim_persistence_check)
write.csv(selflim_persistence_check, file.path("interesting_plots", "selflim_persistence_check.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Persistence check at tau=10: rel_amplitude_P at t_max=300, 1500, 3000 -- ",
    "%.4f, %.4f, %.4f.\n",
    "If these hold roughly steady (unlike section 3's collapse toward 0), ",
    "that's a genuine limit cycle, not an unsettled transient.\n"
  ),
  selflim_persistence_check$rel_amplitude_P[1],
  selflim_persistence_check$rel_amplitude_P[2],
  selflim_persistence_check$rel_amplitude_P[3]
))

################################################################################
# The actual "here's what a Hopf bifurcation looks like" plot: time series
# below vs. above the onset, side by side.
################################################################################

selflim_ts_below <- as.data.frame(dede(
  y = c(N = 1, P = 0.1), times = seq(0, 300, by = 0.1),
  func = predator_prey_rhs_selflim,
  parms = modifyList(predator_prey_params_delay, list(tau = 1, N0 = 1))
)) %>% mutate(regime = "tau = 1 (below onset)")

selflim_ts_above <- as.data.frame(dede(
  y = c(N = 1, P = 0.1), times = seq(0, 300, by = 0.1),
  func = predator_prey_rhs_selflim,
  parms = modifyList(predator_prey_params_delay, list(tau = 10, N0 = 1))
)) %>% mutate(regime = "tau = 10 (above onset)")

selflim_ts_both <- bind_rows(selflim_ts_below, selflim_ts_above) %>%
  mutate(regime = factor(regime, levels = c("tau = 1 (below onset)", "tau = 10 (above onset)")))

selflim_ts_plot <- ggplot(selflim_ts_both, aes(x = time, y = P)) +
  geom_line(color = "#4C72B0", linewidth = 0.7) +
  facet_wrap(~regime, ncol = 1, scales = "free_y") +
  labs(x = "Time", y = "Predator biomass (P)",
       title = "Self-limitation-delay model: decay below onset vs. sustained oscillation above it",
       subtitle = "Same perturbed start (N=1, P=0.1) both panels -- this is what an actual Hopf bifurcation looks like") +
  theme_minimal()
selflim_ts_plot

save_plot(selflim_ts_plot, "Predator-prey toy model - selflim timeseries below vs above onset.png",
         width = 8, height = 8)

################################################################################
# 5. A third variant: delay only P(t-tau) inside the predator's growth term,
#    N(t) instantaneous there
#
# dN/dt = D*(K - N(t))      - a*N(t)*P(t)          <- unchanged, instantaneous
# dP/dt = e*a*N(t)*P(t-tau) - m*P(t)               <- only P delayed here
#
# Structurally between the other two: alpha = D + a*P* (D included, as in
# section 3), but the characteristic equation's beta*gamma term no longer
# carries exp(-lambda*tau) (same effect as section 4 -- both roots of the
# resulting quadratic end up forced to the same sign).
#
# Linearising (N* = m/(e*a), P* = D*(K-N*)/(a*N*), unchanged):
#   dn/dt = -alpha*n(t) - beta*p(t)                [alpha = D+a*P*, as section 3]
#   dp/dt = gamma*n(t) + m*p(t-tau) - m*p(t)        [gamma*n(t) now instantaneous;
#                                                      m*p(t-tau) still delayed,
#                                                      inherited from N* times
#                                                      the still-delayed P(t-tau)]
#
# Characteristic equation: (lambda+alpha)*(lambda+m-m*exp(-lambda*tau)) +
#                           beta*gamma = 0
# (beta*gamma has NO exp(-lambda*tau) factor here -- that's the whole
# structural difference from delaying n as well, as in section 3.)
#
# Setting lambda=i*omega and matching moduli gives a quadratic in
# u = omega^2:
#   u^2 + (alpha^2 - 2*beta*gamma)*u + beta*gamma*(2*alpha*m + beta*gamma) = 0
#
# The product of roots, beta*gamma*(2*alpha*m+beta*gamma), is always
# positive (every factor positive) -- same situation as section 4, so both
# roots share a sign and the sum decides which:
#
#   Hopf possible  <=>  2*beta*gamma > alpha^2
#
# Writing P* = D*C (C = (K-N*)/(a*N*), independent of D, same trick as
# section 4) gives alpha = D*(1+a*C) and beta*gamma = a*m*C*D, and
# substituting these in reduces the condition to a genuine THRESHOLD on D
# (unlike section 3's absolute impossibility, and unlike section 4's
# condition which held outright at the reference parameters):
#
#   D < 2*a*m*C / (1+a*C)^2
#
# At this project's reference parameters (D=1, a=1, e=0.5, m=0.3, K=1 ->
# N*=0.6, C=0.667) that threshold works out to ~0.144. D=1 is far above
# it, so no Hopf at the reference D -- tested below alongside D=0.1
# (comfortably under the threshold) to confirm the derived formula
# actually discriminates correctly on both sides, not just at one point.
################################################################################

# Only P is delayed here; N(t) is instantaneous, so lagvalue(, 2) pulls just
# the P component (state order is c(N, P), so index 2 is P).
predator_prey_rhs_pdelay <- function(t, state, parms) {
  N <- max(state[1], 0)
  P <- max(state[2], 0)
  D   <- parms[["D"]];  a <- parms[["a"]]; e <- parms[["e"]]
  m   <- parms[["m"]];  K <- parms[["K"]]; tau <- parms[["tau"]]
  P0  <- parms[["P0"]]

  if (t <= tau) {
    lag_P <- P0
  } else {
    lag_P <- max(lagvalue(t - tau, 2), 0)
  }

  dN <- D * (K - N) - a * N * P
  dP <- e * a * N * lag_P - m * P
  list(c(dN, dP))
}

# Numeric check of the derived condition, in both forms: 2*beta*gamma >
# alpha^2, and the equivalent D-threshold reduction.
check_hopf_feasibility_pdelay <- function(parms = predator_prey_params_delay) {
  D <- parms$D; a <- parms$a; e <- parms$e; m <- parms$m; K <- parms$K
  N_star  <- m / (e * a)
  C       <- (K - N_star) / (a * N_star)
  P_star  <- D * C
  alpha_p <- D * (1 + a * C)
  beta_p  <- a * N_star
  gamma_p <- e * a * P_star
  D_threshold <- 2 * a * m * C / (1 + a * C)^2
  data.frame(N_star = N_star, P_star = P_star, C = C,
             alpha_p = alpha_p, beta_gamma = beta_p * gamma_p,
             D = D, D_threshold = D_threshold,
             hopf_possible = D < D_threshold)
}

hopf_check_pdelay <- check_hopf_feasibility_pdelay()
print(hopf_check_pdelay)
cat(sprintf(
  paste0(
    "P-delay Hopf check at reference D=%.4f: threshold D < %.4f -- %s\n",
    "(derived condition: Hopf possible iff D < 2*a*m*C/(1+a*C)^2)\n"
  ),
  hopf_check_pdelay$D, hopf_check_pdelay$D_threshold,
  ifelse(hopf_check_pdelay$hopf_possible, "HOPF POSSIBLE", "not possible")
))

# Rerun the same check at a smaller D, comfortably under the threshold, so
# the derived formula is tested on both sides rather than asserted at only
# the one reference value.
params_small_D <- modifyList(predator_prey_params_delay, list(D = 0.1))
hopf_check_pdelay_smallD <- check_hopf_feasibility_pdelay(params_small_D)
print(hopf_check_pdelay_smallD)
cat(sprintf(
  "At D=%.4f (below threshold %.4f): %s\n",
  hopf_check_pdelay_smallD$D, hopf_check_pdelay_smallD$D_threshold,
  ifelse(hopf_check_pdelay_smallD$hopf_possible, "HOPF POSSIBLE", "not possible")
))

run_predator_prey_pdelay <- function(tau, state0 = c(N = 1, P = 0.1), t_max = 300,
                                     parms = predator_prey_params_delay) {
  parms$tau <- tau
  parms$P0  <- state0[["P"]]
  out  <- as.data.frame(dede(y = state0, times = seq(0, t_max, by = 0.1),
                             func = predator_prey_rhs_pdelay, parms = parms))
  late <- out[out$time > t_max * 0.6, ]
  data.frame(tau = tau, max_N = max(late$N), min_N = min(late$N),
             max_P = max(late$P), min_P = min(late$P),
             rel_amplitude_P = (max(late$P) - min(late$P)) / mean(late$P))
}

tau_seq_pdelay <- c(0, 0.5, 1, 2, 3, 5, 7, 10, 15, 20)

# Reference D=1: above threshold, so no sustained oscillation predicted --
# amplitude should decay toward 0 at every tau, same signature as section 3.
pdelay_sweep_refD <- bind_rows(lapply(tau_seq_pdelay, run_predator_prey_pdelay)) %>%
  mutate(D_label = "D = 1 (above threshold, no Hopf predicted)")

# D=0.1: below threshold, so a genuine onset is predicted somewhere in this
# tau range.
pdelay_sweep_smallD <- bind_rows(lapply(tau_seq_pdelay, run_predator_prey_pdelay,
                                        parms = params_small_D)) %>%
  mutate(D_label = "D = 0.1 (below threshold, Hopf predicted)")

pdelay_sweep_df <- bind_rows(pdelay_sweep_refD, pdelay_sweep_smallD)

print(pdelay_sweep_df)
write.csv(pdelay_sweep_df, file.path("interesting_plots", "predator_prey_pdelay_sweep.csv"),
          row.names = FALSE)

pdelay_amplitude_plot <- ggplot(pdelay_sweep_df, aes(x = tau, y = rel_amplitude_P, color = D_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(x = "tau (delay in P(t-tau) term only, N instantaneous)",
       y = "Relative amplitude of P (settled)",
       title = "Delaying only P inside the conversion term: a threshold on D, not an absolute block",
       subtitle = "Same a=1,e=0.5,m=0.3,K=1 -- derived threshold D < 0.144 tested on both sides",
       color = NULL) +
  theme_minimal()
pdelay_amplitude_plot

save_plot(pdelay_amplitude_plot, "Predator-prey toy model - pdelay sweep amplitude, D above and below threshold.png")

################################################################################
# Persistence check at D=0.1 (below threshold): confirm any oscillation
# found there is a genuine limit cycle, not another relaxation artefact --
# same style of check as sections 3 and 4.
################################################################################

pdelay_persistence_check <- bind_rows(
  run_predator_prey_pdelay(tau = 10, t_max = 300,  parms = params_small_D) %>% mutate(t_max = 300),
  run_predator_prey_pdelay(tau = 10, t_max = 1500, parms = params_small_D) %>% mutate(t_max = 1500),
  run_predator_prey_pdelay(tau = 10, t_max = 3000, parms = params_small_D) %>% mutate(t_max = 3000)
)

print(pdelay_persistence_check)
write.csv(pdelay_persistence_check, file.path("interesting_plots", "pdelay_persistence_check.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "P-delay persistence check at D=0.1, tau=10: rel_amplitude_P at ",
    "t_max=300, 1500, 3000 -- %.4f, %.4f, %.4f.\n",
    "Steady across these -> genuine limit cycle; collapsing toward 0 -> just ",
    "another unsettled transient.\n"
  ),
  pdelay_persistence_check$rel_amplitude_P[1],
  pdelay_persistence_check$rel_amplitude_P[2],
  pdelay_persistence_check$rel_amplitude_P[3]
))
