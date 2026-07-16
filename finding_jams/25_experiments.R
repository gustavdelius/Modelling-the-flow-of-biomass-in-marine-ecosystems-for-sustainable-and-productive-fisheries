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

################################################################################
# Rerunning the alpha/gamma/kappa juvenile-pileup isolation as a real sweep
#
# Day 24's version (24_experiments.R) only ever swapped each of alpha,
# gamma, ks, and kappa to ONE alternate value -- the other species' own
# default -- not a real range, so the resulting ranking (alpha dominant,
# gamma sufficient-but-not-necessary, ks/kappa smaller, not-yet-disentangled
# contributors) wasn't solid enough to trust. It also surfaced a genuine
# confound: mizer's own default-parameter calibration is interdependent, so
# leaving a parameter unspecified in a "Default params, anchovy X" call let
# mizer silently re-derive OTHER parameters too (caught directly: a
# "Default params, anchovy kappa" call came back reporting gamma = 145, not
# the default species' own 2066, even though nothing in that call touched
# gamma).
#
# This version sidesteps that confound entirely rather than working around
# it: every sweep starts from the anchovy's own fully-specified parameter
# set (every biological argument always hardcoded, nothing left for mizer to
# fill in), and only the ONE parameter under test is varied across a real
# range. Every other parameter stays explicitly pinned at the anchovy's own
# value throughout, so there is no re-derivation to worry about, and every
# sweep point is directly comparable to every other.
#
# ks was dropped after a first pass: swept 0-6, it collapses the population
# by ~50 orders of magnitude on its own (a straightforward metabolic-cost
# effect, not the same kind of juvenile-pileup question the other three are
# testing), which swamped the shared y-axis on the combined plot and made
# alpha/gamma/kappa's own, much more interesting variation look flat by
# comparison. ks stays pinned at the anchovy's own value (0) in every sweep
# below, same as before -- just no longer swept in its own right.
################################################################################

# Self-contained, same as every script since Day 20 -- redefined here rather
# than sourced from 24_experiments.R. Extended with alpha/gamma/ks/kappa
# arguments (all defaulting to the anchovy's own values), so the same single
# builder can sweep any one of the four while pinning the rest.
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

# The default species' own alpha/gamma/kappa, purely as reference lines on
# the sweep plots below -- not used to build anything, so there's no risk of
# reintroducing the re-derivation confound this rerun is designed to avoid.
default_species_raw <- newSingleSpeciesParams(species_name = "Anchovy")
default_alpha <- default_species_raw@species_params$alpha[1]
default_gamma <- default_species_raw@species_params$gamma[1]
default_kappa <- tryCatch(resource_params(default_species_raw)$kappa,
                          error = function(e) NA_real_)
anchovy_kappa <- 100 * exp(-6.9 * (2.05 - 1))  # same formula make_second_order_params_kr() uses

cat(sprintf(
  paste0(
    "Reference values (anchovy | default species):\n",
    "  alpha: %.4f | %.4f\n  gamma: %.1f | %.1f\n  kappa: %.4f | %s\n"
  ),
  0.1, default_alpha, 750, default_gamma, anchovy_kappa,
  ifelse(is.na(default_kappa), "NOT FOUND", sprintf("%.4f", default_kappa))
))

# Same catchable_fraction() logic as 24_experiments.R's run_spectrum_check(),
# pulled out on its own since this version doesn't need the "which builder,
# which label" machinery that existed to compare anchovy-vs-default params --
# every call here uses the same anchovy-based builder, varying one argument.
# Takes an already-built params object, so it can only ever fail inside
# project() -- the viability check below happens earlier, at construction
# time, which is why run_param_sweep() (not this function) is what needs to
# wrap the whole thing in tryCatch.
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

# capacity_mult = 10, resource_decrease = 0.001: same operating point every
# isolation check in 24_experiments.R used -- comfortably inside the
# coexistence region so a marginal population never confounds the pileup
# reading, and directly comparable to Day 24's single-swap numbers.
#
# The tryCatch has to wrap BOTH the params construction and the simulation,
# not just the simulation -- newSingleSpeciesParams() (inside
# make_second_order_params_kr()) does its own viability check at
# construction time ("the feeding level is not sufficient to maintain the
# fish") and errors out hard before project() is ever reached. A first
# version of this sweep only wrapped catchable_fraction_at()'s call to
# project(), which left that earlier construction-time error completely
# unprotected -- exactly the failure mode Day 24's run_spectrum_check()
# guarded against, and this rerun is specifically testing wider parameter
# ranges than Day 24 did, so hitting it here is expected, not a fluke.
run_param_sweep <- function(param_name, param_seq, capacity_mult = 10,
                            resource_decrease = 0.001) {
  bind_rows(lapply(param_seq, function(v) {
    args <- list(capacity_mult = capacity_mult, resource_decrease = resource_decrease)
    args[[param_name]] <- v
    result <- tryCatch({
      p <- do.call(make_second_order_params_kr, args)
      catchable_fraction_at(p)
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

plot_param_sweep <- function(df, param_label, ref_values, log_x = TRUE) {
  p <- ggplot(df, aes(x = value, y = catchable_fraction)) +
    geom_line(color = "#4C72B0", linewidth = 1) +
    geom_point(color = "#4C72B0", size = 1.5) +
    # percent_format() with a fixed decimal accuracy collapses every break
    # below its rounding threshold to the same "0.0000%" label -- exactly
    # what happens here, since catchable_fraction spans from near-extinct up
    # to ~11% across a single sweep. Scientific notation on the raw fraction
    # has no such threshold, so every decade of the log scale stays legible.
    scale_y_log10(labels = scales::label_scientific()) +
    labs(x = param_label, y = "Catchable fraction (biomass at or above w_mat)",
         title = sprintf("Juvenile-pileup sweep: %s", param_label),
         subtitle = "capacity_mult = 10, resource_decrease = 0.001 -- every other parameter pinned at the anchovy's own value") +
    theme_minimal()
  if (log_x) p <- p + scale_x_log10()
  for (rv in ref_values) {
    p <- p + geom_vline(xintercept = rv$value, linetype = rv$linetype, color = "grey40")
  }
  p
}

################################################################################
# alpha (assimilation efficiency): anchovy = 0.1, default species' own value
# printed above. Day 24's single swap found this the dominant lever (5+
# orders of magnitude either way), so this is the one most worth resolving
# into a real curve rather than two endpoints.
################################################################################

alpha_seq <- exp(seq(log(0.02), log(0.8), length.out = 15))

alpha_sweep_df <- run_param_sweep("alpha", alpha_seq)

saveRDS(alpha_sweep_df, file.path("interesting_plots", "alpha_sweep_df.rds"))
write.csv(alpha_sweep_df, file.path("interesting_plots", "alpha_sweep_df.csv"), row.names = FALSE)

alpha_sweep_plot <- plot_param_sweep(
  alpha_sweep_df, "alpha (assimilation efficiency)",
  ref_values = list(list(value = 0.1, linetype = "dashed"),
                    list(value = default_alpha, linetype = "dotted"))
)
alpha_sweep_plot

save_plot(alpha_sweep_plot, "juv_alpha.png")

################################################################################
# gamma (search volume / intake rate): anchovy = 750, default species' own
# value printed above. Day 24's single swap found this sufficient but not
# necessary -- moved the default species a lot, barely moved the anchovy.
################################################################################

gamma_seq <- exp(seq(log(100), log(5000), length.out = 15))

gamma_sweep_df <- run_param_sweep("gamma", gamma_seq)

saveRDS(gamma_sweep_df, file.path("interesting_plots", "gamma_sweep_df.rds"))
write.csv(gamma_sweep_df, file.path("interesting_plots", "gamma_sweep_df.csv"), row.names = FALSE)

gamma_sweep_plot <- plot_param_sweep(
  gamma_sweep_df, "gamma (search volume)",
  ref_values = list(list(value = 750, linetype = "dashed"),
                    list(value = default_gamma, linetype = "dotted"))
)
gamma_sweep_plot

save_plot(gamma_sweep_plot, "juv_gamma.png")

################################################################################
# kappa (background resource-spectrum coefficient): anchovy's own value
# (anchovy_kappa, computed above) vs. the default species' own value.
################################################################################

kappa_seq <- exp(seq(log(0.001), log(1), length.out = 15))

kappa_sweep_df <- run_param_sweep("kappa_override", kappa_seq) %>%
  mutate(param = "kappa")

saveRDS(kappa_sweep_df, file.path("interesting_plots", "kappa_sweep_df.rds"))
write.csv(kappa_sweep_df, file.path("interesting_plots", "kappa_sweep_df.csv"), row.names = FALSE)

kappa_sweep_plot <- plot_param_sweep(
  kappa_sweep_df, "kappa (background resource spectrum)",
  ref_values = list(list(value = anchovy_kappa, linetype = "dashed"),
                    list(value = default_kappa, linetype = "dotted"))
)
kappa_sweep_plot

save_plot(kappa_sweep_plot, "juv_kappa.png")

################################################################################
# Putting the three sweeps side by side, and ranking them the same way Day
# 24's single-swap version tried to: how many orders of magnitude does each
# parameter move the catchable fraction on its own, across a real range,
# with every other parameter pinned.
################################################################################

all_sweeps_df <- bind_rows(alpha_sweep_df, gamma_sweep_df, kappa_sweep_df)
write.csv(all_sweeps_df, file.path("interesting_plots", "all_param_sweeps.csv"), row.names = FALSE)

combined_sweep_plot <- ggplot(all_sweeps_df %>% filter(is.na(error)),
                              aes(x = value, y = catchable_fraction)) +
  geom_line(color = "#4C72B0", linewidth = 1) +
  geom_point(color = "#4C72B0", size = 1) +
  facet_wrap(~param, scales = "free_x", nrow = 1) +
  scale_y_log10(labels = scales::label_scientific()) +
  labs(x = "Parameter value", y = "Catchable fraction",
       title = "Juvenile-pileup sweeps: alpha, gamma, kappa",
       subtitle = "Each panel varies one parameter with the other two pinned at the anchovy's own value") +
  theme_minimal()
combined_sweep_plot

save_plot(combined_sweep_plot, "juv_all3.png", width = 14, height = 5)

sweep_magnitude <- all_sweeps_df %>%
  filter(is.na(error), catchable_fraction > 0) %>%
  group_by(param) %>%
  summarise(
    min_fraction = min(catchable_fraction),
    max_fraction = max(catchable_fraction),
    log10_range  = log10(max_fraction) - log10(min_fraction),
    .groups = "drop"
  ) %>%
  arrange(desc(log10_range))

print(sweep_magnitude)
write.csv(sweep_magnitude, file.path("interesting_plots", "param_sweep_magnitude.csv"), row.names = FALSE)

cat(paste0(
  "Parameter sweep ranking (orders of magnitude the catchable fraction moves\n",
  "across each parameter's own tested range, every other parameter pinned):\n",
  "  Compares directly against Day 24's single-swap ranking (alpha dominant,\n",
  "  gamma sufficient-but-not-necessary, kappa smaller and not disentangled --\n",
  "  ks dropped from this rerun, see the header comment above) -- see\n",
  "  param_sweep_magnitude.csv for the numbers, and the three individual sweep\n",
  "  plots for whether each relationship is smooth/monotonic or has a sharp\n",
  "  threshold the single-swap version could never have shown.\n"
))

################################################################################
# Does a time delay in the predator's growth term unlock hysteresis, where
# instantaneous Type I and Type II versions (24_experiments.R) couldn't?
#
# Day 24 found that neither a linear (Type I) nor a saturating (Type II)
# functional response can destabilise the toy predator-prey model, as long
# as prey growth stays chemostat-style -- proven analytically (the Jacobian
# trace stays negative regardless of f's shape) and confirmed by
# simulation. The next candidate mechanism flagged in Day 24's own "What's
# Next": a time delay between a predator eating and that intake showing up
# as predator growth (a maturation/gestation lag), a well-known route to
# oscillation in delay differential equation (DDE) models that an
# instantaneous ODE structurally can't produce.
#
# The specific form used here follows the classical Wangersky-Cunningham
# (1957) delayed-conversion model: predation removes prey from the system
# instantaneously (a predator eats now, so N(t) drops now), but the
# resulting predator biomass gain only shows up tau time units later,
# reflecting a successful "encounter -> offspring" event that started back
# at t - tau. Predator mortality stays instantaneous, since dying isn't
# subject to the same conversion lag as being born:
#
#   dN/dt = D*(K - N(t)) - a*N(t)*P(t)
#   dP/dt = e*a*N(t-tau)*P(t-tau) - m*P(t)
#
# WORKED OUT ANALYTICALLY FIRST, the same way K_crit and the Type II
# no-hysteresis result were derived on Day 24.
#
# Equilibrium: setting both derivatives to zero with N, P held constant (so
# every delayed term equals the un-delayed one) reproduces EXACTLY the same
# N* = m/(e*a) and P* = D*(K-N*)/(a*N*) as the non-delayed model -- the
# delay cannot move where the equilibrium sits, only whether it's stable.
#
# Linearising N = N*+n, P = P*+p around that equilibrium (dropping
# second-order terms, and using e*a*N* = m from the equilibrium condition
# to simplify) gives:
#
#   dn/dt = -alpha*n(t) - beta*p(t)
#   dp/dt = gamma*n(t-tau) + m*p(t-tau) - m*p(t)
#
#   where alpha = D + a*P*, beta = a*N*, gamma = e*a*P*
#
# Substituting n, p ~ exp(lambda*t) gives the characteristic equation
#   (lambda+alpha)*(lambda+m-m*exp(-lambda*tau)) + beta*gamma*exp(-lambda*tau) = 0
# which at tau = 0 reduces to lambda^2 + alpha*lambda + beta*gamma = 0 --
# exactly Day 24's trace = -alpha, det = beta*gamma result, a direct check
# that this reduces correctly to the known non-delayed case.
#
# For tau > 0, a Hopf bifurcation (a pair of roots crossing the imaginary
# axis, lambda = i*omega, as tau increases from 0) needs a REAL, POSITIVE
# omega^2 solving omega^4 + alpha^2*omega^2 + beta*gamma*(2*m*alpha -
# beta*gamma) = 0, which has one only when
#   beta*gamma > 2*m*alpha                                            (*)
#
# Whether (*) can ever hold turns out not to depend on tau, K, or the
# specific parameter values at all: substituting beta = a*N*,
# gamma = e*a*P*, and using e*a*N* = m again,
#   beta*gamma - 2*m*alpha = a^2*e*N**P* - 2*m*(D + a*P*)
#                          = a*P*(a*e*N* - 2*m) - 2*m*D
#                          = a*P*(m - 2*m) - 2*m*D
#                          = -a*m*P* - 2*m*D
# which is strictly NEGATIVE for every positive D, a, m, P* -- so (*) can
# NEVER hold. A Hopf bifurcation is analytically impossible for this
# delayed-conversion model, for ANY delay tau and ANY parameter choice,
# exactly mirroring the Type II saturation result: chemostat-style prey
# growth blocks this route to instability too, not just the Type I/II
# functional-response route. Simulated below anyway, since Day 24's own
# standard is that confirming a null result empirically is worth more than
# trusting the algebra alone.
################################################################################

library(deSolve)

# Same parameters as the Type I/Type II toy models in 24_experiments.R --
# D: prey renewal rate; a: attack rate; e: predator conversion efficiency;
# m: predator mortality; K: prey carrying capacity.
predator_prey_params_delay <- list(D = 1, a = 1, e = 0.5, m = 0.3, K = 1)

# alpha/beta/gamma and both sides of the Hopf-feasibility condition (*)
# above, as a numeric double-check of the algebra -- if this ever printed
# TRUE, the analytic proof would be wrong and the simulation below would be
# the more trustworthy half of this section, not the other way around.
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

# Delay differential equation via deSolve::dede(), which supports lagged
# state access through lagvalue(). For t <= tau there's no history yet, so
# the lagged state falls back to the same initial condition used to start
# the integration -- a constant-history assumption, the standard choice
# when nothing else is known about the system's past. Same clamping lesson
# as the Type I/II models: a transiently negative lagged N or P would flip
# the sign of the a*N*P-style term.
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

# Perturbed well away from equilibrium (N=1 vs N*=0.6, P=0.1 vs P*=0.667)
# rather than starting AT equilibrium -- a system exactly at a stable fixed
# point never oscillates regardless of tau, so this is a genuine test of
# whether delay lets a real perturbation grow into a sustained oscillation
# instead of damping back out. tau up to 40 -- two orders of magnitude past
# any of the model's own natural timescales (1/D, 1/m) -- to stress-test
# the "no Hopf for any tau" prediction rather than only checking modest
# delays.
tau_seq <- c(0, 0.5, 1, 2, 5, 10, 20, 40)

delay_sweep_df <- bind_rows(lapply(tau_seq, run_predator_prey_delay))

print(delay_sweep_df)
write.csv(delay_sweep_df, file.path("interesting_plots", "predator_prey_delay_sweep.csv"),
          row.names = FALSE)

delay_amplitude_plot <- ggplot(delay_sweep_df, aes(x = tau, y = rel_amplitude_P)) +
  geom_line(color = "#55A868", linewidth = 1) +
  geom_point(color = "#55A868", size = 1.5) +
  labs(x = "tau (delay in predator conversion)", y = "Relative amplitude of P (settled)",
       title = "Does delay reintroduce the oscillation the instantaneous model doesn't have?",
       subtitle = "Perturbed start (N=1, P=0.1), D=1, a=1, e=0.5, m=0.3, K=1 -- analytic prediction: no, for any tau") +
  theme_minimal()
delay_amplitude_plot

save_plot(delay_amplitude_plot, "pp_delay_amp.png")

# A dedicated time-series run at the largest tested delay, to see the
# transient directly rather than only reading off the settled-window
# amplitude above -- if delay were going to produce a growing oscillation
# instead of damping to the fixed point, it should be visible here.
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

save_plot(delay_timeseries_plot, "pp_delay_ts_tau40.png")

cat(paste0(
  "Time-delay predator-prey check: if rel_amplitude_P stays small and flat across\n",
  "every tau tested, and the tau=40 time series damps back to N*, P* rather than\n",
  "sustaining or growing an oscillation, that confirms the analytic prediction --\n",
  "a Hopf bifurcation is impossible for this delayed-conversion model regardless of\n",
  "tau, because beta*gamma - 2*m*alpha = -a*m*P* - 2*m*D is always negative for\n",
  "positive parameters. Chemostat-style prey growth blocks this route to\n",
  "hysteresis too, not just the Type I/II functional-response route Day 24 ruled out.\n"
))



