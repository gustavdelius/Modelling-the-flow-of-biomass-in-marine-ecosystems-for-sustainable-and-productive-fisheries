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

# Self-contained, same as every script since Day 20 -- redefined here rather
# than sourced from 23_experiments.R.
make_second_order_params_kr <- function(lambda = 2.05, resource_decrease = 0.001,
                                        capacity_mult = 1, second_order = TRUE,
                                        ext_diff = 0.00) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
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

# Same generic forward/backward bifurcation sweep as 23_experiments.R --
# carried over unchanged so this refinement is directly comparable to the
# sweep that found the loop in the first place.
run_bifurcation_sweep <- function(param_seq, param_name, fixed_params = list(),
                                  params_fn = make_second_order_params_kr,
                                  t_run = 600, lambda = 2.05) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_bm = NA_real_, min_bm = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      args <- fixed_params
      args[[param_name]] <- seq_vals[i]
      args$lambda <- lambda
      p <- do.call(params_fn, args)
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }
      sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                     progress_bar = FALSE, effort = 0, method = "tr_bdf2")
      bm   <- getBiomass(sim)[, "Anchovy"]
      tv   <- as.numeric(names(bm))
      late <- bm[tv > t_run * 0.6]
      last <- dim(sim@n)[1]

      state_n   <- sim@n[last, , ]
      state_npp <- sim@n_pp[last, ]

      out$max_bm[i] <- max(late)
      out$min_bm[i] <- min(late)
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(param_seq)
  bwd    <- run_one_direction(rev(param_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, biomass = fwd$df$max_bm, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, biomass = fwd$df$min_bm, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, biomass = bwd_df$max_bm, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, biomass = bwd_df$min_bm, direction = "Backward", branch = "min")
  )
}

plot_bifurcation <- function(df, x_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = biomass, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_x_log10() +
    labs(x = x_label, y = "Biomass", title = title, subtitle = subtitle) +
    theme_minimal()
}

dir.create("interesting_plots", showWarnings = FALSE)

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 150) {
  ggsave(file.path("interesting_plots", filename), plot = plot,
         width = width, height = height, dpi = dpi)
}

################################################################################
# Refining the capacity_mult ~ 1.6-2.1 hysteresis window
#
# Day 23's 20-point capacity_mult bifurcation sweep (cc_seq_bif, log-spaced
# 1-100, resource_decrease held at 0.001) turned up a genuine forward/backward
# hysteresis loop, not an artefact: at capacity_mult = 1.62 Forward (arriving
# from nothing) was still dead (7.71e-24) while Backward (arriving from an
# established population at capacity_mult = 100) was already alive (6.66e-4)
# -- about 20 orders of magnitude apart at the identical parameter value.
# That grid only bracketed the two edges of the loop: Backward's collapse
# point somewhere in (1.27, 1.62], Forward's escape point somewhere in
# (1.62, 2.07). This is the same move as Follow-up 4's [3, 10] plateau
# refinement in 23_experiments.R -- rerun the same kind of sweep at much
# higher density, zoomed into just the bracketed window, instead of guessing
# where the two edges actually sit.
#
# t_run is kept at the sweep's own default (600, matching the run that found
# the loop) rather than cut down the way the phase-diagram grid refinements
# were. This is specifically a hysteresis-width measurement, so trading
# settle time for point density here would confound a genuine scan-rate
# effect (Day 23's own Follow-up 1) with the resolution increase this is
# meant to isolate.
################################################################################

# Brackets both suspected edges (1.27 and 2.07) with margin on either side,
# at roughly 4x the points-per-decade of the original 20-point/1-100 sweep.
cc_seq_refine <- exp(seq(log(1.2), log(2.3), length.out = 25))

bif_cc_refine_df <- run_bifurcation_sweep(cc_seq_refine, "capacity_mult",
                                          fixed_params = list(resource_decrease = 0.001),
                                          params_fn = make_second_order_params_kr)

saveRDS(bif_cc_refine_df, file.path("interesting_plots", "bif_cc_refine_df.rds"))

bif_cc_refine_plot <- plot_bifurcation(
  bif_cc_refine_df, "capacity_mult",
  "Refined bifurcation: capacity_mult in [1.2, 2.3]",
  "resource_decrease held at 0.001 (balance = FALSE) -- 25 points vs. the original 20 over [1, 100]"
)
bif_cc_refine_plot

save_plot(bif_cc_refine_plot, "Bifurcation - capacity_mult refined 1.2-2.3.png")

################################################################################
# Pin down both edges numerically, not just by eye on the plot
#
# "Dead" here means numerically-collapsed to floating-point-underflow levels
# (Day 23's dead values ranged 1e-24 to 1e-57), not just "small" -- a real,
# viable population right at the edge of this loop can itself be small
# (Backward's 6.66e-4 and 1.46e-3 near-boundary points were both called
# "alive" in Day 23). MIN_VIABLE_BIOMASS (1e-2) is tuned for a different job
# -- filtering numerical noise out of the relative-amplitude oscillation
# score in the phase-diagram sweeps -- and is too strict to use here, since it
# would wrongly call some genuinely-alive near-boundary points dead. A
# threshold of 1e-6 sits in the wide, empty gap between every dead value and
# every alive value seen in Day 23's table, so it cleanly separates the two
# without needing to be precisely tuned.
################################################################################

ALIVE_THRESHOLD <- 1e-6

edge_calls <- bif_cc_refine_df %>%
  filter(branch == "max") %>%
  transmute(capacity_mult = value, direction, biomass, alive = biomass > ALIVE_THRESHOLD) %>%
  arrange(direction, capacity_mult)

print(edge_calls)
write.csv(edge_calls, file.path("interesting_plots", "capacity_mult_edge_calls.csv"),
          row.names = FALSE)

# For each direction, walking the swept value ascending: the last point
# still called dead, and the first point called alive. The true transition
# sits somewhere in (dead_below, alive_above]. Returns NA for either side
# rather than erroring if a group is entirely dead or entirely alive.
# value_col defaults to "capacity_mult" (every call above uses that); the
# non-size-structured model below sweeps "K" instead, hence the parameter
# rather than a third hardcoded copy of this logic.
find_transition <- function(df, value_col = "capacity_mult") {
  df         <- df %>% arrange(.data[[value_col]])
  dead_vals  <- df %>% filter(!alive) %>% pull(.data[[value_col]])
  alive_vals <- df %>% filter(alive)  %>% pull(.data[[value_col]])
  data.frame(
    dead_below  = if (length(dead_vals)  > 0) max(dead_vals)  else NA_real_,
    alive_above = if (length(alive_vals) > 0) min(alive_vals) else NA_real_
  )
}

transitions <- edge_calls %>%
  group_by(direction) %>%
  group_modify(~ find_transition(.x)) %>%
  ungroup()

print(transitions)
write.csv(transitions, file.path("interesting_plots", "capacity_mult_transitions.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Refined window (capacity_mult, resource_decrease = 0.001):\n",
    "  Backward collapses somewhere in (%.4f, %.4f]\n",
    "  Forward escapes    somewhere in (%.4f, %.4f]\n"
  ),
  transitions$dead_below[transitions$direction == "Backward"],
  transitions$alive_above[transitions$direction == "Backward"],
  transitions$dead_below[transitions$direction == "Forward"],
  transitions$alive_above[transitions$direction == "Forward"]
))

################################################################################
# Does the capacity_mult bifurcation point move with resource rate?
#
# Gustav's prediction (from the shape of Day 23's original snake-sweep plot,
# capacity_mult x resource_decrease traced as a phase diagram): a higher
# resource_decrease -- more of the default resource rate surviving, i.e. a
# faster-replenishing resource -- should need a SMALLER carrying capacity to
# sustain the population, moving the capacity_mult bifurcation point left.
# Everything so far (Day 23's original discovery, the refinement above) only
# tested this at a single resource_decrease = 0.001, so there's no way yet to
# tell whether that's a general relationship or one cherrypicked rate.
#
# Reruns the capacity_mult bifurcation sweep at the same five
# resource_decrease values as Day 23's original snake grid (rd_grid_snake in
# 23_experiments.R) rather than picking new ones -- those are literally the
# points behind the plot Gustav was reading when he made the prediction.
################################################################################

rd_grid_multi <- exp(seq(log(0.0001), log(0.5), length.out = 5))

# Capacity range widened well past [1, 100] on both sides relative to the
# single-rate sweeps above: if Gustav's prediction holds, the lowest
# resource_decrease (0.0001, more depleted than the 0.001 baseline) should
# need MORE capacity than the (1.30, 1.41) window found above, and the
# highest (0.5) should need LESS -- possibly below capacity_mult = 1
# entirely. Density is coarser per decade than the single-rate sweeps -- this
# is a first look at the trend across rates, not a precision pass on any one
# of them; a refinement in the style of the one above can follow up on
# whichever rate needs it once the rough location is known.
cc_seq_multi <- exp(seq(log(0.05), log(200), length.out = 12))

# t_run cut from the single-rate sweeps' 600 down to 150: Day 23's own
# scan-rate check (23_experiments.R) found t_run = 150 already agreed with
# t_run = 600 everywhere on a similar grid, so this trades some settle-time
# margin for a ~4x speedup rather than guessing at a safe shortcut. 12
# points instead of 20 for a further ~1.7x on top of that -- lower accuracy
# is fine here since this pass is only meant to locate the rough
# neighbourhood, not pin down an exact edge.
bif_multi_df <- bind_rows(lapply(rd_grid_multi, function(rd) {
  run_bifurcation_sweep(cc_seq_multi, "capacity_mult",
                        fixed_params = list(resource_decrease = rd),
                        params_fn = make_second_order_params_kr,
                        t_run = 150) %>%
    mutate(resource_decrease = rd)
})) %>%
  mutate(rd_label = factor(sprintf("resource_decrease = %.4g", resource_decrease),
                           levels = sprintf("resource_decrease = %.4g", sort(unique(resource_decrease)))))

saveRDS(bif_multi_df, file.path("interesting_plots", "bif_multi_df.rds"))

bif_multi_plot <- ggplot(bif_multi_df, aes(x = value, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  scale_x_log10() +
  facet_wrap(~rd_label, nrow = 1, scales = "free_y") +
  labs(x = "capacity_mult", y = "Biomass",
       title = "capacity_mult bifurcation at five different resource rates",
       subtitle = "If Gustav's prediction holds, the collapse boundary should shift left as resource_decrease increases") +
  theme_minimal()
bif_multi_plot

save_plot(bif_multi_plot, "Bifurcation - capacity_mult across resource rates.png", width = 18, height = 5)

################################################################################
# Reading the bifurcation point off each panel numerically
#
# Same alive/dead classification and find_transition() as the refinement
# above, generalised across resource_decrease and direction instead of just
# direction.
################################################################################

edge_calls_multi <- bif_multi_df %>%
  filter(branch == "max") %>%
  transmute(resource_decrease, capacity_mult = value, direction,
           biomass, alive = biomass > ALIVE_THRESHOLD) %>%
  arrange(resource_decrease, direction, capacity_mult)

transitions_multi <- edge_calls_multi %>%
  group_by(resource_decrease, direction) %>%
  group_modify(~ find_transition(.x)) %>%
  ungroup() %>%
  arrange(resource_decrease, direction)

print(transitions_multi)
write.csv(transitions_multi, file.path("interesting_plots", "capacity_mult_transitions_by_rate.csv"),
          row.names = FALSE)

# The actual test: does the bifurcation location -- geometric mean of each
# direction's dead_below/alive_above bracket, consistent with the log-spaced
# grid -- decrease as resource_decrease increases? Both directions are kept
# side by side so a widening/narrowing hysteresis gap across rates is visible
# too, not just the overall leftward/rightward shift.
bifurcation_trend <- transitions_multi %>%
  mutate(bifurcation_estimate = sqrt(dead_below * alive_above)) %>%
  select(resource_decrease, direction, bifurcation_estimate) %>%
  tidyr::pivot_wider(names_from = direction, values_from = bifurcation_estimate) %>%
  arrange(resource_decrease)

cat("Bifurcation estimate vs. resource_decrease (expect DECREASING if Gustav's prediction holds):\n")
print(bifurcation_trend)
write.csv(bifurcation_trend, file.path("interesting_plots", "capacity_mult_bifurcation_trend.csv"),
          row.names = FALSE)

################################################################################
# Zooming in: capacity_mult bifurcation across resource rates, refined to
# [1, 15]
#
# The original [0.05, 200], 20-point, t_run = 600 grid gave the SAME
# bifurcation estimate for every resource_decrease value (Backward = 1.32,
# Forward = 4.89, no variation at all) -- but the grid step there is ~55%
# (the ratio between adjacent log-spaced points over 3.6 decades in 20
# steps), while the actual window at a single rate (resource_decrease =
# 0.001) was already found above to be just (1.30, 1.41), an 8% span. A
# shift narrower than the grid's own resolution would look exactly like no
# shift at all -- so that flat result doesn't yet distinguish "the boundary
# really doesn't move" from "the grid is too coarse to see it move". Same
# fix as every other refinement in this file: narrow the window to where
# the previous coarse pass says the action is, and raise the point count.
################################################################################

cc_seq_multi_refine <- exp(seq(log(1), log(15), length.out = 20))

# Same t_run = 150 shortcut as the coarse pass above -- Day 23's own
# scan-rate check found it agreed with t_run = 600 on a similar grid, so
# this is a deliberate accuracy-for-speed trade, not an unvalidated guess.
# Point count halved from the original 40 to 20 for a further speedup on top
# of that; still 4x the coarse pass's density over a much narrower range.
bif_multi_refine_df <- bind_rows(lapply(rd_grid_multi, function(rd) {
  run_bifurcation_sweep(cc_seq_multi_refine, "capacity_mult",
                        fixed_params = list(resource_decrease = rd),
                        params_fn = make_second_order_params_kr,
                        t_run = 150) %>%
    mutate(resource_decrease = rd)
})) %>%
  mutate(rd_label = factor(sprintf("resource_decrease = %.4g", resource_decrease),
                           levels = sprintf("resource_decrease = %.4g", sort(unique(resource_decrease)))))

saveRDS(bif_multi_refine_df, file.path("interesting_plots", "bif_multi_refine_df.rds"))

bif_multi_refine_plot <- ggplot(bif_multi_refine_df, aes(x = value, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1) +
  scale_x_log10() +
  facet_wrap(~rd_label, nrow = 1, scales = "free_y") +
  labs(x = "capacity_mult", y = "Biomass",
       title = "capacity_mult bifurcation across resource rates, refined to [1, 15]",
       subtitle = "20 points, t_run = 150 vs. the coarse pass's 12 points over [0.05, 200], t_run = 150") +
  theme_minimal()
bif_multi_refine_plot

save_plot(bif_multi_refine_plot, "Bifurcation - capacity_mult across resource rates, refined 1-15.png",
         width = 18, height = 5)

################################################################################
# Reading the refined bifurcation point off each panel numerically
################################################################################

edge_calls_multi_refine <- bif_multi_refine_df %>%
  filter(branch == "max") %>%
  transmute(resource_decrease, capacity_mult = value, direction,
           biomass, alive = biomass > ALIVE_THRESHOLD) %>%
  arrange(resource_decrease, direction, capacity_mult)

transitions_multi_refine <- edge_calls_multi_refine %>%
  group_by(resource_decrease, direction) %>%
  group_modify(~ find_transition(.x)) %>%
  ungroup() %>%
  arrange(resource_decrease, direction)

print(transitions_multi_refine)
write.csv(transitions_multi_refine,
          file.path("interesting_plots", "capacity_mult_transitions_by_rate_refined.csv"),
          row.names = FALSE)

bifurcation_trend_refine <- transitions_multi_refine %>%
  mutate(bifurcation_estimate = sqrt(dead_below * alive_above)) %>%
  select(resource_decrease, direction, bifurcation_estimate) %>%
  tidyr::pivot_wider(names_from = direction, values_from = bifurcation_estimate) %>%
  arrange(resource_decrease)

cat(paste0(
  "Refined bifurcation estimate vs. resource_decrease (expect DECREASING if Gustav's prediction holds --\n",
  "if this is STILL flat at 20 points over [1, 15], that's evidence for genuine invariance, not just\n",
  "an under-resolved grid):\n"
))
print(bifurcation_trend_refine)
write.csv(bifurcation_trend_refine,
          file.path("interesting_plots", "capacity_mult_bifurcation_trend_refined.csv"),
          row.names = FALSE)

################################################################################
# A non-size-structured two-species predator-prey model: does the collapse
# bifurcation need size resolution at all?
#
# Gustav's question from the 2026-07-08 meeting: the mizer model's collapse
# as capacity_mult decreases below ~1.30-1.41 (refined above) is basically a
# predator-prey system with size structure layered on top. Does the collapse
# bifurcation itself need that size structure, or would the same qualitative
# behaviour already show up in the simplest possible unstructured version --
# a 2-ODE prey/predator system with no age, size, or stage distinction at
# all?
#
# The model, exactly as discussed: prey grows chemostat-style (a
# replenishment rate toward a carrying capacity, the same functional form as
# mizer's resource_semichemostat), predator's intake is proportional to the
# product of prey and predator abundance (mass-action encounter, no
# saturation), and predator mortality is proportional to its own abundance:
#
#   dN/dt = D * (K - N) - a * N * P
#   dP/dt = e * a * N * P - m * P
#
# Setting dP/dt = 0 (excluding the trivial P = 0 branch) gives the
# predator's break-even prey level, N* = m / (e * a) -- independent of K.
# Coexistence needs the prey-alone equilibrium (N = K) to exceed that level,
# i.e. K > m / (e * a). Crossing that threshold from below is a
# transcritical bifurcation: the prey-only state and the coexistence state
# exchange stability, exactly the "fixed point joining the unstable fixed
# point at zero" Gustav described watching the earlier capacity_mult sweep.
# This is standard consumer-resource theory (Tilman's R*; see Smith &
# Waltman, "The Theory of the Chemostat") -- the analytic K_crit below is a
# check that the simulation behaves the way the textbook math says it
# should, not just a guess pulled from the sweep.
#
# Caveat worth flagging up front: this is a Type I (linear, unsaturated)
# functional response, which is what "proportional to prey x predator"
# means literally. Classical theory (Rosenzweig 1971, the "paradox of
# enrichment") says a Type I model like this one has NO Hopf bifurcation for
# any K -- it can show the collapse-at-low-K transcritical bifurcation
# tested here, but it structurally cannot reproduce the separate
# oscillation-onset-at-high-capacity_mult behaviour from Day 22. That's a
# different question (whether oscillation needs size structure) needing a
# saturating (Type II) intake term instead -- not attempted here, since
# Gustav's specific question was about the collapse bifurcation, not the
# oscillation.
################################################################################

library(deSolve)

predator_prey_rhs <- function(t, state, parms) {
  # Positional, not named, indexing into state: confirmed directly (the
  # state[["N"]] version above raised "subscript out of bounds") that
  # deSolve's default lsoda strips state's names during internal numerical
  # Jacobian estimation, so name-based lookup can't be relied on inside this
  # function at all. Position is safe -- deSolve guarantees state's order
  # always matches the original y = c(N = ..., P = ...) vector, even on the
  # calls where the names themselves are gone. parms isn't perturbed by the
  # solver the way state is, so name-based indexing into it stays safe.
  # Clamped to >= 0: the model is bounded and only ever spirals into the
  # coexistence equilibrium (the Jacobian's trace, -D - a*P*, is always
  # negative once that equilibrium exists), so it can't blow up on its own.
  # But an under-damped oscillation sampled mid-swing at the end of one
  # sweep step can hand the next step a transiently negative N or P; a
  # negative N times a positive P flips the sign of the a*N*P term, and that
  # sign flip is what produced the NaN TOLSF failure -- lsoda's step-size
  # adaptation choking on a fast, badly-conditioned transient the model was
  # never meant to produce. N, P are biological abundances, so clamping to
  # zero isn't discarding real dynamics, just refusing to integrate a
  # negative population that shouldn't exist in the first place.
  N <- max(state[1], 0)
  P <- max(state[2], 0)
  D <- parms[["D"]]
  a <- parms[["a"]]
  e <- parms[["e"]]
  m <- parms[["m"]]
  K <- parms[["K"]]

  dN <- D * (K - N) - a * N * P
  dP <- e * a * N * P - m * P
  list(c(dN, dP))
}

# D: prey renewal rate; a: attack rate; e: predator conversion efficiency;
# m: predator mortality. Kept as simple, round default numbers rather than
# reverse-engineered from the anchovy parameters -- the point of this model
# is to test for the qualitative bifurcation shape, not to match mizer's
# numbers exactly.
predator_prey_params <- list(D = 1, a = 1, e = 0.5, m = 0.3, K = 1)
N_star <- with(predator_prey_params, m / (e * a))
K_crit <- N_star  # coexistence needs K > N_star exactly, since N* doesn't depend on K
cat(sprintf("Analytic prediction: N* = %.4f, coexistence requires K > %.4f\n", N_star, K_crit))

# Falls back to the last available point (with a warning) rather than
# silently returning a zero-length settle window if ode() ever terminates
# before reaching t_max -- an empty `late` here would otherwise produce a
# zero-length n_final/p_final, which is exactly what caused the "initial
# conditions vector (0)" failure when that got fed into the next ode() call
# downstream in run_bifurcation_sweep_pp().
run_predator_prey <- function(K, state0 = c(N = 1, P = 0.1), t_max = 500, parms = predator_prey_params) {
  parms$K <- K
  out  <- as.data.frame(ode(y = state0, times = seq(0, t_max, by = 0.5),
                            func = predator_prey_rhs, parms = parms))
  late <- out[out$time > t_max * 0.6, ]
  if (nrow(late) == 0) {
    warning(sprintf(
      "K = %.4f: integration stopped at t = %.2f of %.2f -- falling back to the last available point",
      K, max(out$time), t_max
    ))
    late <- tail(out, 1)
  }
  data.frame(K = K, max_P = max(late$P), min_P = min(late$P), mean_P = mean(late$P),
             n_final = tail(late$N, 1), p_final = tail(late$P, 1))
}

# Quick sanity check before the full sweep: does a single run at K below and
# above K_crit behave the way the analytic threshold predicts? Analytic
# coexistence equilibrium (only valid for K > K_crit): N* as above, P* from
# setting dN/dt = 0 at that N*.
P_star_at <- function(K, parms = predator_prey_params) {
  # Explicit indexing here too -- with(parms, ...) would evaluate against
  # parms$K (always 1, predator_prey_params's hardcoded default), silently
  # shadowing this function's own K argument instead of using whatever K
  # was actually passed in.
  D <- parms[["D"]]
  a <- parms[["a"]]
  D * (K - N_star) / (a * N_star)
}

sanity_below <- run_predator_prey(K = K_crit * 0.5)
sanity_above <- run_predator_prey(K = K_crit * 2)

cat(sprintf(
  paste0(
    "Sanity check:\n",
    "  K = %.3f (below K_crit): simulated P settles to %.6f (expect ~0)\n",
    "  K = %.3f (above K_crit): simulated P settles to %.6f, analytic P* = %.6f\n"
  ),
  K_crit * 0.5, sanity_below$p_final,
  K_crit * 2, sanity_above$p_final, P_star_at(K_crit * 2)
))

# NOT a state-carried forward/backward sweep, unlike every capacity_mult
# sweep above -- deliberately so. A transcritical bifurcation has no
# hysteresis at all (Forward and Backward are only ever expected to
# disagree exactly at K_crit itself), so state-carrying was never going to
# reveal anything real here. Worse, it actively breaks: this is a
# single-scalar Lotka-Volterra-style model with no immigration term, so
# dP/dt = 0 identically once P underflows to exact floating-point zero.
# Starting the sweep well below K_crit (as K_seq_pp does, for margin)
# guarantees several early points where P decays toward zero as expected --
# but once it hits TRUE zero there, carrying that state forward means P
# stays zero for every later K too, no matter how favourable, since there's
# nothing left to regrow it from. A first attempt at this sweep hit exactly
# that: Forward read as ~1e-31 (numerically extinct) all the way out to
# K = 3, deep in territory the analytic K_crit = 0.6 says should sustain a
# real, order-1 population, and Backward then inherited that same stuck-at-
# zero state and stayed there for its entire run too.
#
# The fix is to test each K independently, fresh from the same default
# initial condition every time -- exactly what run_predator_prey() already
# does, and what the sanity check above already validated at two points.
K_seq_pp <- exp(seq(log(K_crit * 0.2), log(K_crit * 5), length.out = 25))

pp_fresh_df <- bind_rows(lapply(K_seq_pp, run_predator_prey))

saveRDS(pp_fresh_df, file.path("interesting_plots", "predator_prey_fresh_df.rds"))

# Analytic P*(K): only a real, meaningful solution for K > K_crit (below
# that, the formula gives a negative number, which isn't a population --
# the analytic answer there is just P = 0).
pp_analytic_df <- data.frame(K = K_seq_pp) %>%
  mutate(P_star = ifelse(K > K_crit, vapply(K, P_star_at, numeric(1)), NA_real_))

pp_fresh_plot <- ggplot(pp_fresh_df, aes(x = K, y = p_final)) +
  geom_line(color = "#4C72B0", linewidth = 1) +
  geom_point(color = "#4C72B0", size = 1.5) +
  geom_line(data = pp_analytic_df, aes(x = K, y = P_star),
           color = "grey40", linetype = "dashed", inherit.aes = FALSE) +
  geom_vline(xintercept = K_crit, linetype = "dotted", color = "grey40") +
  scale_x_log10() +
  labs(x = "K (prey carrying capacity)", y = "Predator biomass (P), settled",
       title = "Non-size-structured predator-prey model: fresh run at each K",
       subtitle = "Solid = simulated; dashed = analytic P*(K); dotted = analytic K_crit -- no state-carrying, so no extinction-trapping artefact") +
  theme_minimal()
pp_fresh_plot

save_plot(pp_fresh_plot, "Predator-prey toy model - fresh runs vs K.png")

# Self-contained transition finder, not reused from the capacity_mult
# section above -- that's exactly the kind of cross-section dependency that
# broke this section once already (find_transition()'s signature had
# changed between sessions). Cheap enough to duplicate; this section
# shouldn't silently break again if that one changes.
find_transition_pp <- function(df) {
  df         <- df %>% arrange(K)
  dead_vals  <- df %>% filter(!alive) %>% pull(K)
  alive_vals <- df %>% filter(alive)  %>% pull(K)
  data.frame(
    dead_below  = if (length(dead_vals)  > 0) max(dead_vals)  else NA_real_,
    alive_above = if (length(alive_vals) > 0) min(alive_vals) else NA_real_
  )
}

pp_transition <- pp_fresh_df %>%
  transmute(K, alive = p_final > ALIVE_THRESHOLD) %>%
  find_transition_pp()

write.csv(pp_transition, file.path("interesting_plots", "predator_prey_transition.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Non-size-structured toy model bifurcation check (fresh runs, no state-carrying):\n",
    "  Analytic prediction:         K_crit = %.4f\n",
    "  Simulated transition bracket: (%.4f, %.4f]\n"
  ),
  K_crit, pp_transition$dead_below, pp_transition$alive_above
))
