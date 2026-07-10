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

################################################################################
# Is the capacity_mult bifurcation SHAPE a consequence of our own choice of
# species parameters, or does it show up with mizer's own defaults too?
#
# Every capacity_mult sweep so far (Day 22 onward) has used the same
# anchovy-like calibration -- w_min = 0.0003, w_max = 66.5, w_mat = 10,
# lambda = 2.05, alpha = 0.1, gamma = 750, ks = 0 -- carried over unchanged
# from Day 17. That's one specific point in a large parameter space; nothing
# so far has checked whether the collapse-at-low-capacity shape survives
# swapping in a completely different species. newSingleSpeciesParams()
# called with none of those arguments gives mizer's own built-in defaults
# instead -- a different, independently-chosen parameterisation -- so
# rerunning the same kind of sweep against it is the direct test.
#
# make_default_species_params_kr() mirrors make_second_order_params_kr()
# exactly (same D_ext/setBevertonHolt() order, same independent
# capacity_mult/resource_decrease knobs, same balance = FALSE), with the
# species itself as the only thing that changes. It still accepts a lambda
# argument purely so it stays call-compatible with run_bifurcation_sweep(),
# which always injects one -- the argument is unused here since nothing
# about the default species is being overridden.
################################################################################

make_default_species_params_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                            second_order = TRUE, ext_diff = 0.00,
                                            lambda = NULL) {
  # species_name is a label, not a biological parameter -- kept as "Anchovy"
  # purely so run_bifurcation_sweep()'s shared harness (which hardcodes
  # getBiomass(sim)[, "Anchovy"] for extracting the biomass series) doesn't
  # need touching. Every actual biological argument (w_min, w_max, w_mat,
  # lambda, kappa, alpha, gamma, ks, ...) is left unspecified, so mizer
  # fills in its own defaults instead of this project's anchovy calibration.
  params <- newSingleSpeciesParams(species_name = "Anchovy")

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

# Wide and sparse rather than a multi-stage refinement: the point here is
# just to see whether a collapse-then-rise shape shows up at all for a
# different species, not to pin its edges down precisely the way the
# anchovy version was. t_run = 150 reuses the same validated shortcut as
# the multi-rate sweeps above; capacity_mult spans 5 orders of magnitude
# since there's no prior information about where a default species'
# threshold would fall.
cc_seq_default <- exp(seq(log(0.01), log(1000), length.out = 15))

bif_default_df <- run_bifurcation_sweep(cc_seq_default, "capacity_mult",
                                        fixed_params = list(resource_decrease = 0.001),
                                        params_fn = make_default_species_params_kr,
                                        t_run = 150)

saveRDS(bif_default_df, file.path("interesting_plots", "bif_default_species_df.rds"))

bif_default_plot <- plot_bifurcation(
  bif_default_df, "capacity_mult",
  "capacity_mult bifurcation with mizer's own default species parameters",
  "resource_decrease held at 0.001, t_run = 150 -- same sweep style as the anchovy version, different species"
)
bif_default_plot

save_plot(bif_default_plot, "Bifurcation - capacity_mult, default species params.png")

edge_calls_default <- bif_default_df %>%
  filter(branch == "max") %>%
  transmute(capacity_mult = value, direction, biomass, alive = biomass > ALIVE_THRESHOLD) %>%
  arrange(direction, capacity_mult)

transitions_default <- edge_calls_default %>%
  group_by(direction) %>%
  group_modify(~ find_transition(.x)) %>%
  ungroup()

print(transitions_default)
write.csv(transitions_default, file.path("interesting_plots", "capacity_mult_transitions_default_species.csv"),
          row.names = FALSE)

cat(paste0(
  "Default-species capacity_mult check: if this shows the same collapse-then-rise\n",
  "shape (with or without a hysteresis gap) as the anchovy version, that's evidence\n",
  "the shape is a general mizer/consumer-resource feature, not an artefact of this\n",
  "project's specific species calibration. If it looks qualitatively different (no\n",
  "collapse, a different-shaped transition, or no hysteresis at all), the anchovy\n",
  "parameters are doing more work than assumed.\n"
))

################################################################################
# Zooming into the default-species hysteresis region
#
# The coarse 15-point/[0.01, 1000] pass above already resolved real,
# non-overlapping brackets -- Backward (0.61, 1.39], Forward (3.16, 7.20] --
# a genuine ~2.3x gap, not grid-coarseness noise. Same move as every other
# refinement in this file: narrow the window to bracket both edges with
# margin, and raise the point count.
#
# t_run is put back to the sweep's own default (600), not the t_run = 150
# shortcut used for the coarse pass above. That shortcut was validated
# specifically for the anchovy-calibrated model (Day 23's scan-rate check);
# it hasn't been separately checked against this different default species,
# so reusing it here would be stacking one unvalidated assumption on
# another. The window is narrow enough now that going back to full accuracy
# doesn't cost much.
################################################################################

cc_seq_default_refine <- exp(seq(log(0.5), log(8), length.out = 25))

bif_default_refine_df <- run_bifurcation_sweep(cc_seq_default_refine, "capacity_mult",
                                               fixed_params = list(resource_decrease = 0.001),
                                               params_fn = make_default_species_params_kr,
                                               t_run = 600)

saveRDS(bif_default_refine_df, file.path("interesting_plots", "bif_default_refine_df.rds"))

bif_default_refine_plot <- plot_bifurcation(
  bif_default_refine_df, "capacity_mult",
  "Refined capacity_mult bifurcation, default species: [0.5, 8]",
  "resource_decrease held at 0.001, t_run = 600 -- 25 points vs. the coarse pass's 15 over [0.01, 1000]"
)
bif_default_refine_plot

save_plot(bif_default_refine_plot, "Bifurcation - capacity_mult, default species, refined 0.5-8.png")

edge_calls_default_refine <- bif_default_refine_df %>%
  filter(branch == "max") %>%
  transmute(capacity_mult = value, direction, biomass, alive = biomass > ALIVE_THRESHOLD) %>%
  arrange(direction, capacity_mult)

transitions_default_refine <- edge_calls_default_refine %>%
  group_by(direction) %>%
  group_modify(~ find_transition(.x)) %>%
  ungroup()

print(transitions_default_refine)
write.csv(transitions_default_refine,
          file.path("interesting_plots", "capacity_mult_transitions_default_refined.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Refined default-species hysteresis window:\n",
    "  Backward collapses somewhere in (%.4f, %.4f]\n",
    "  Forward escapes    somewhere in (%.4f, %.4f]\n"
  ),
  transitions_default_refine$dead_below[transitions_default_refine$direction == "Backward"],
  transitions_default_refine$alive_above[transitions_default_refine$direction == "Backward"],
  transitions_default_refine$dead_below[transitions_default_refine$direction == "Forward"],
  transitions_default_refine$alive_above[transitions_default_refine$direction == "Forward"]
))

################################################################################
# Is the juvenile-biomass pileup general, or specific to our anchovy
# calibration?
#
# Day 19 found the anchovy size spectrum drops several orders of magnitude
# right at w_mat = 10: only ~0.031% (default kernel) to ~0.0033%
# (highest-amplitude kernel) of total biomass was ever "catchable" (mature).
# That was flagged as an open demographic question (survivorship through
# growth -- not yet explained) and carried forward unresolved through Days
# 20-24. The default-species capacity_mult check earlier in this file only
# looked at AGGREGATE biomass and found the same collapse-then-rise shape;
# it never looked at how that biomass is distributed across sizes. This is
# the actual pileup check: same idea (mizer's own default species instead
# of the anchovy calibration), but looking at the size spectrum itself.
#
# capacity_mult = 10 rather than 1: the refined sweep above found
# capacity_mult = 1 sits right at this species' own collapse boundary
# (Backward barely alive there, Forward dead until ~1.26) -- a marginal,
# not-fully-established population would be a poor basis for a clean
# size-spectrum comparison. capacity_mult = 10 is comfortably inside the
# coexistence region for this species, matching the kind of healthy,
# established regime Day 19's original anchovy spectrum was taken from.
################################################################################

default_spectrum_params <- make_default_species_params_kr(resource_decrease = 0.001,
                                                           capacity_mult = 10, lambda = 2.05)

default_spectrum_sim <- project(default_spectrum_params, t_max = 600, dt = 0.1, t_save = 0.5,
                                progress_bar = FALSE, effort = 0, method = "tr_bdf2")

last_default_spectrum <- dim(default_spectrum_sim@n)[1]
default_w_mat         <- default_spectrum_params@species_params$w_mat[1]

default_spectrum_df <- data.frame(
  w               = default_spectrum_sim@params@w,
  biomass_density = default_spectrum_sim@n[last_default_spectrum, 1, ] * default_spectrum_sim@params@w
)

default_spectrum_plot <- ggplot(default_spectrum_df, aes(x = w, y = biomass_density)) +
  geom_line(linewidth = 1, color = "#4C72B0") +
  geom_vline(xintercept = default_w_mat, linetype = "dashed", color = "grey40") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Body weight (g)", y = "Biomass density (N x w)",
       title = "Size spectrum, mizer's own default species parameters",
       subtitle = sprintf(
         "capacity_mult = 10, resource_decrease = 0.001 -- dashed line = this species' own w_mat = %.3g",
         default_w_mat
       )) +
  theme_minimal()
default_spectrum_plot

save_plot(default_spectrum_plot, "Size spectrum - default species params.png")

# Same catchable_fraction() quantification Day 19 used for the anchovy, so
# the two are directly comparable numbers, not just similar-looking plots.
catchable_fraction_default <- function(sim, cutoff) {
  last         <- dim(sim@n)[1]
  w            <- sim@params@w
  dw           <- sim@params@dw
  bm_density   <- sim@n[last, 1, ] * w * dw
  total_bm     <- sum(bm_density)
  catchable_bm <- sum(bm_density[w >= cutoff])
  data.frame(total_bm = total_bm, catchable_bm = catchable_bm,
             catchable_fraction = catchable_bm / total_bm)
}

default_catchable <- catchable_fraction_default(default_spectrum_sim, cutoff = default_w_mat)

print(default_catchable)
write.csv(default_catchable, file.path("interesting_plots", "default_species_catchable_fraction.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Default-species juvenile-pileup check (capacity_mult = 10, resource_decrease = 0.001):\n",
    "  w_mat = %.4g\n",
    "  Fraction of total biomass at or above w_mat (\"catchable\"): %.5f%%\n",
    "Day 19's anchovy figures for comparison: 0.031%% (default kernel), 0.0033%% (highest-amplitude kernel).\n",
    "If this default species' fraction is similarly tiny (well under 1%%), that's evidence the juvenile\n",
    "pileup is a general mizer/growth-mortality feature, not an artefact of this project's anchovy\n",
    "calibration. A fraction much closer to 50%% (or otherwise much larger) would suggest the anchovy's\n",
    "specific w_min/w_max/mortality choices are what's driving the extreme skew.\n"
  ),
  default_w_mat, default_catchable$catchable_fraction * 100
))

################################################################################
# Does a saturating (Type II) functional response unlock the mizer-style
# hysteresis, or does the chemostat prey growth rule that out regardless?
#
# Worked out analytically before touching the sweep, the same way K_crit
# was derived for the Type I model. General 2-species system:
#   dN/dt = g(N) - f(N)*P,  dP/dt = e*f(N)*P - m*P
# At the coexistence equilibrium, f(N*) = m/e always, which makes the
# dP/dt-vs-P Jacobian entry (e*f(N*) - m) exactly zero regardless of f's
# shape. That leaves:
#   trace = g'(N*) - f'(N*)*P*
#   det   = m * f'(N*) * P*          (always positive if f is increasing)
# so stability depends entirely on trace's sign. For CHEMOSTAT prey growth,
# g(N) = D*(K - N), g'(N) = -D -- a constant, independent of K and of
# whether f is Type I (f' = a) or Type II (f' = a/(1+a*h*N)^2, still > 0).
# trace = -D - f'(N*)*P* is negative either way, for every K where the
# predator survives at all. The paradox of enrichment specifically needs
# the prey's OWN growth curve to weaken (or reverse) as K grows at fixed
# N* -- true for logistic growth r*N*(1-N/K), structurally impossible for
# chemostat growth. So the analytic prediction here is: Type II saturation
# alone, with prey growth kept chemostat-style as specified, should still
# give a single transcritical bifurcation -- a different K_crit formula,
# but no oscillation and no hysteresis. Confirming that empirically is more
# convincing than trusting the algebra, and if it holds, it sharpens the
# open question to "needs logistic prey growth too", not just "needs
# Type II" -- not a wasted run either way.
################################################################################

predator_prey_rhs_typeII <- function(t, state, parms) {
  # Same positional-indexing and clamping lessons as the Type I model above.
  N <- max(state[1], 0)
  P <- max(state[2], 0)
  D <- parms[["D"]]
  a <- parms[["a"]]
  e <- parms[["e"]]
  m <- parms[["m"]]
  h <- parms[["h"]]
  K <- parms[["K"]]

  intake <- a * N / (1 + a * h * N)  # Holling Type II -- reduces to a*N (Type I) as h -> 0

  dN <- D * (K - N) - intake * P
  dP <- e * intake * P - m * P
  list(c(dN, dP))
}

# h = 1: comfortably below the e/m = 1.667 ceiling above which the predator
# can never break even even at infinite prey -- as N -> infinity, intake
# maxes out at 1/h, so viability at any N at all needs e/h > m, i.e.
# h < e/m.
predator_prey_params_typeII <- list(D = 1, a = 1, e = 0.5, m = 0.3, h = 1, K = 1)

# f(N*) = m/e => a*N* / (1 + a*h*N*) = m/e => N* = m / (a*(e - h*m)),
# valid only while e > h*m (checked above). K_crit = N*, same logic as the
# Type I model: coexistence needs the prey-alone equilibrium (N = K) to
# exceed the predator's break-even level.
N_star_typeII <- with(predator_prey_params_typeII, m / (a * (e - h * m)))
K_crit_typeII <- N_star_typeII
cat(sprintf("Type II analytic prediction: N* = %.4f, coexistence requires K > %.4f\n",
           N_star_typeII, K_crit_typeII))

run_predator_prey_typeII <- function(K, state0 = c(N = 1, P = 0.1), t_max = 500,
                                     parms = predator_prey_params_typeII) {
  parms$K <- K
  out  <- as.data.frame(ode(y = state0, times = seq(0, t_max, by = 0.5),
                            func = predator_prey_rhs_typeII, parms = parms))
  late <- out[out$time > t_max * 0.6, ]
  if (nrow(late) == 0) late <- tail(out, 1)
  data.frame(K = K, max_P = max(late$P), min_P = min(late$P), mean_P = mean(late$P),
             n_final = tail(late$N, 1), p_final = tail(late$P, 1))
}

# Same sanity check pattern as the Type I model: does a single run below
# and above K_crit_typeII behave the way the analytic threshold predicts?
sanity_below_typeII <- run_predator_prey_typeII(K = K_crit_typeII * 0.5)
sanity_above_typeII  <- run_predator_prey_typeII(K = K_crit_typeII * 2)

cat(sprintf(
  paste0(
    "Type II sanity check:\n",
    "  K = %.3f (below K_crit): simulated P settles to %.6f (expect ~0)\n",
    "  K = %.3f (above K_crit): simulated P settles to %.6f\n"
  ),
  K_crit_typeII * 0.5, sanity_below_typeII$p_final,
  K_crit_typeII * 2, sanity_above_typeII$p_final
))

# Same "fresh run at each K, no state-carrying" fix as the Type I model --
# a transcritical bifurcation (which is what the analytic derivation above
# predicts this still is, even with Type II) has no hysteresis to detect,
# and state-carrying can only break here the same way it did before: once P
# underflows to exact zero anywhere in the low-K region, there's nothing
# left to regrow it from, no matter how favourable K becomes later.
K_seq_typeII <- exp(seq(log(K_crit_typeII * 0.2), log(K_crit_typeII * 5), length.out = 25))

pp_typeII_df <- bind_rows(lapply(K_seq_typeII, run_predator_prey_typeII))

saveRDS(pp_typeII_df, file.path("interesting_plots", "predator_prey_typeII_df.rds"))

pp_typeII_plot <- ggplot(pp_typeII_df, aes(x = K, y = p_final)) +
  geom_line(color = "#DD8452", linewidth = 1) +
  geom_point(color = "#DD8452", size = 1.5) +
  geom_vline(xintercept = K_crit_typeII, linetype = "dotted", color = "grey40") +
  scale_x_log10() +
  labs(x = "K (prey carrying capacity)", y = "Predator biomass (P), settled",
       title = "Type II (saturating) predator-prey model: fresh run at each K",
       subtitle = sprintf("h = %.2f -- dotted line = analytic K_crit = %.3f",
                          predator_prey_params_typeII$h, K_crit_typeII)) +
  theme_minimal()
pp_typeII_plot

save_plot(pp_typeII_plot, "Predator-prey toy model - Type II, fresh runs vs K.png")

find_transition_pp_typeII <- function(df) {
  df         <- df %>% arrange(K)
  dead_vals  <- df %>% filter(!alive) %>% pull(K)
  alive_vals <- df %>% filter(alive)  %>% pull(K)
  data.frame(
    dead_below  = if (length(dead_vals)  > 0) max(dead_vals)  else NA_real_,
    alive_above = if (length(alive_vals) > 0) min(alive_vals) else NA_real_
  )
}

pp_typeII_transition <- pp_typeII_df %>%
  transmute(K, alive = p_final > ALIVE_THRESHOLD) %>%
  find_transition_pp_typeII()

write.csv(pp_typeII_transition, file.path("interesting_plots", "predator_prey_typeII_transition.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Type II toy model bifurcation check (fresh runs, no state-carrying):\n",
    "  Analytic prediction:         K_crit = %.4f\n",
    "  Simulated transition bracket: (%.4f, %.4f]\n"
  ),
  K_crit_typeII, pp_typeII_transition$dead_below, pp_typeII_transition$alive_above
))

################################################################################
# The direct hysteresis test: does state-carrying reveal any real
# forward/backward disagreement here, or does it just confirm the analytic
# prediction of none?
#
# Run a Forward (ascending K, state-carried) and Backward (descending K,
# state-carried from Forward's endpoint) pass, exactly like the mizer
# sweeps -- but starting comfortably ABOVE K_crit_typeII rather than deep in
# extinction territory, so the absorbing-zero failure mode that broke the
# Type I state-carried attempt can't trigger here either. If Forward and
# Backward agree everywhere, that's the "no hysteresis" prediction holding
# up under direct test, not just algebra. Any visible gap would mean the
# analytic argument missed something, which would be the more interesting
# result.
################################################################################

run_state_carried_typeII <- function(K_seq, parms = predator_prey_params_typeII, t_max = 500) {
  run_one_direction <- function(seq_vals, init_state = c(N = 1, P = 0.1)) {
    out   <- vector("list", length(seq_vals))
    state <- init_state
    for (i in seq_along(seq_vals)) {
      p <- parms
      p$K <- seq_vals[i]
      sim   <- as.data.frame(ode(y = state, times = seq(0, t_max, by = 0.5),
                                 func = predator_prey_rhs_typeII, parms = p))
      late  <- sim[sim$time > t_max * 0.6, ]
      if (nrow(late) == 0) {
        warning(sprintf(
          "K = %.4f: integration stopped at t = %.2f of %.2f -- falling back to the last available point",
          seq_vals[i], max(sim$time), t_max
        ))
        late <- tail(sim, 1)
      }
      state <- pmax(c(N = tail(late$N, 1), P = tail(late$P, 1)), 0)
      out[[i]] <- data.frame(K = seq_vals[i], max_P = max(late$P), min_P = min(late$P))
    }
    list(df = bind_rows(out), final_state = state)
  }

  fwd    <- run_one_direction(K_seq)
  bwd    <- run_one_direction(rev(K_seq), init_state = fwd$final_state)
  bwd_df <- bwd$df[order(bwd$df$K), ]

  bind_rows(
    data.frame(K = fwd$df$K, biomass = fwd$df$max_P, direction = "Forward",  branch = "max"),
    data.frame(K = fwd$df$K, biomass = fwd$df$min_P, direction = "Forward",  branch = "min"),
    data.frame(K = bwd_df$K, biomass = bwd_df$max_P, direction = "Backward", branch = "max"),
    data.frame(K = bwd_df$K, biomass = bwd_df$min_P, direction = "Backward", branch = "min")
  )
}

# Starting at 0.8x K_crit_typeII rather than 0.2x like the fresh sweep above
# -- deliberately avoiding the low-K graveyard that broke the Type I
# state-carried attempt, since the point here is whether Forward and
# Backward disagree ANYWHERE in a range both can actually survive, not to
# re-trigger the same absorbing-zero failure mode.
K_seq_typeII_carried <- exp(seq(log(K_crit_typeII * 0.8), log(K_crit_typeII * 5), length.out = 20))

pp_typeII_carried_df <- run_state_carried_typeII(K_seq_typeII_carried)

saveRDS(pp_typeII_carried_df, file.path("interesting_plots", "predator_prey_typeII_carried_df.rds"))

pp_typeII_carried_plot <- ggplot(pp_typeII_carried_df,
                                 aes(x = K, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  geom_vline(xintercept = K_crit_typeII, linetype = "dotted", color = "grey40") +
  scale_x_log10() +
  labs(x = "K (prey carrying capacity)", y = "Predator biomass (P)",
       title = "Type II predator-prey model: state-carried Forward/Backward",
       subtitle = "Started at K = 0.8 * K_crit_typeII, comfortably clear of the low-K absorbing-zero region") +
  theme_minimal()
pp_typeII_carried_plot

save_plot(pp_typeII_carried_plot, "Predator-prey toy model - Type II, forward-backward.png")

cat(paste0(
  "Type II forward/backward check: if Forward and Backward overlap everywhere in\n",
  "this plot, that confirms the analytic prediction directly -- chemostat prey growth\n",
  "keeps this model's trace negative regardless of functional response, so Type II\n",
  "saturation alone (without also making prey growth density-dependent/logistic) isn't\n",
  "enough to reproduce mizer's hysteresis. Any visible gap between the two lines would\n",
  "be a real, unexpected finding worth digging into further.\n"
))

################################################################################
# Is gamma specifically what's driving the juvenile-biomass pileup?
#
# The default-species check found ~7.1% of biomass above w_mat, versus the
# anchovy's 0.031%/0.0033% -- 200-2000x less extreme. That points at this
# project's own alpha/gamma/ks/w_min/w_max calibration rather than mizer's
# growth-mortality structure in general, but doesn't say WHICH of those is
# doing the work. gamma is the natural first suspect: it scales search
# volume/intake rate directly, so a much higher gamma should mean faster
# growth through the juvenile size range, which could plausibly explain why
# so little biomass survives to reach w_mat.
#
# A clean 2x2 isolates it: species parameters (anchovy vs. default) crossed
# with gamma (anchovy's 750 vs. the default species' own value), holding
# everything else about each species fixed. If gamma alone explains the gap,
# "anchovy params + default gamma" should land close to the default
# baseline (~7%), and "default params + anchovy gamma" should land close to
# the anchovy's extreme end (well under 1%). If gamma only explains part of
# the gap, both swapped variants will land somewhere in between -- still
# informative, just not the whole story.
################################################################################

# The default species' own gamma, pulled from the params object already
# built for the earlier juvenile-pileup check, not re-derived -- same
# species build both times, so this is the exact value that check used.
default_gamma <- default_spectrum_params@species_params$gamma[1]
cat(sprintf("Default species' own gamma: %.4f (anchovy's is 750)\n", default_gamma))

# Anchovy species parameters, but gamma is swappable instead of hardcoded --
# every other anchovy-specific argument (w_min, w_max, w_mat, alpha, ks,
# kappa) is unchanged from make_second_order_params_kr().
make_anchovy_gamma_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                       second_order = TRUE, ext_diff = 0.00,
                                       lambda = 2.05, gamma_override) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = gamma_override, ks = 0
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

# Default species parameters (mizer's own w_min/w_max/w_mat/alpha/ks/kappa),
# but gamma is swappable instead of left to mizer's own default.
make_default_gamma_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                       second_order = TRUE, ext_diff = 0.00,
                                       lambda = NULL, gamma_override) {
  params <- newSingleSpeciesParams(species_name = "Anchovy", gamma = gamma_override)

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

# Same spectrum-and-catchable_fraction logic as the juvenile-pileup check
# above, generalised into one function so all four variants are computed
# identically rather than by hand four times.
run_spectrum_check <- function(params_fn, label, capacity_mult = 10,
                               resource_decrease = 0.001, t_max = 600, ...) {
  # Wrapped in tryCatch: newSingleSpeciesParams() does its own viability
  # check at construction time ("the feeding level is not sufficient to
  # maintain the fish") and errors out hard if a swapped parameter makes the
  # species infeasible before any dynamics even run. That's a real result
  # for an isolation test like this one -- it means the swap pushes the
  # species past viability entirely -- not just a technical failure, so
  # it's recorded as a row (catchable_fraction = NA, error message kept)
  # rather than halting the whole isolation table.
  tryCatch({
    p   <- params_fn(resource_decrease = resource_decrease, capacity_mult = capacity_mult, ...)
    sim <- project(p, t_max = t_max, dt = 0.1, t_save = 0.5,
                   progress_bar = FALSE, effort = 0, method = "tr_bdf2")
    last  <- dim(sim@n)[1]
    w_mat <- p@species_params$w_mat[1]
    w     <- sim@params@w
    dw    <- sim@params@dw

    bm_density   <- sim@n[last, 1, ] * w * dw
    total_bm     <- sum(bm_density)
    catchable_bm <- sum(bm_density[w >= w_mat])

    # Reports gamma, alpha, ks, and kappa regardless of which one a given
    # call is actually swapping, so this one function serves every
    # isolation test without needing a near-duplicate each time -- and so a
    # swap that didn't actually take effect would be visible in the output,
    # not silently assumed.
    #
    # kappa's correct retrieval location wasn't actually verified against a
    # live mizer session (two guesses -- species_params$kappa, then
    # resource_params()$kappa -- both silently returned NULL rather than
    # erroring, which is why tryCatch alone didn't catch it: there was
    # never an error to catch, just an empty result breaking data.frame()'s
    # row-count check). safe_scalar() below defends against BOTH failure
    # modes -- a thrown error AND a silent NULL/zero-length result -- so
    # this can degrade to a clean NA instead of breaking every row's
    # construction a third time.
    safe_scalar <- function(expr) {
      v <- tryCatch(expr, error = function(e) NULL)
      if (is.null(v) || length(v) == 0) NA_real_ else v[1]
    }
    kappa_val <- safe_scalar(resource_params(p)$kappa)
    data.frame(label = label, gamma = p@species_params$gamma[1],
               alpha = p@species_params$alpha[1], ks = p@species_params$ks[1],
               kappa = kappa_val, w_mat = w_mat,
               catchable_fraction = catchable_bm / total_bm, error = NA_character_)
  }, error = function(e) {
    warning(sprintf("%s: %s", label, conditionMessage(e)))
    data.frame(label = label, gamma = NA_real_, alpha = NA_real_, ks = NA_real_,
               kappa = NA_real_, w_mat = NA_real_, catchable_fraction = NA_real_,
               error = conditionMessage(e))
  })
}

# capacity_mult = 10 for all four: already established as comfortably
# inside the coexistence region for both the pure anchovy and pure default
# species, so reused here for consistency rather than re-deriving each
# swapped variant's own collapse threshold just for this targeted check.
gamma_isolation_df <- bind_rows(
  run_spectrum_check(make_second_order_params_kr,
                     label = "Anchovy params, anchovy gamma (750)"),
  run_spectrum_check(make_anchovy_gamma_swap_kr,
                     label = "Anchovy params, default gamma",
                     gamma_override = default_gamma),
  run_spectrum_check(make_default_species_params_kr,
                     label = "Default params, default gamma"),
  run_spectrum_check(make_default_gamma_swap_kr,
                     label = "Default params, anchovy gamma (750)",
                     gamma_override = 750)
)

print(gamma_isolation_df)
write.csv(gamma_isolation_df, file.path("interesting_plots", "gamma_isolation_df.csv"),
          row.names = FALSE)

cat(paste0(
  "Gamma isolation check (capacity_mult = 10, resource_decrease = 0.001 throughout):\n",
  "  If gamma alone explains the pileup gap: 'Anchovy params, default gamma' should land\n",
  "  close to the default baseline (~7%), and 'Default params, anchovy gamma' should land\n",
  "  close to the anchovy's extreme end (well under 1%). If both swapped variants instead\n",
  "  land somewhere in between the two baselines, gamma is part of the story but not all\n",
  "  of it -- alpha, ks, or w_min/w_max would still need checking.\n"
))

################################################################################
# Is alpha the other piece of the puzzle?
#
# The gamma isolation test found gamma is SUFFICIENT to reproduce the
# extreme pileup starting from the default species (7.13% -> ~0.04% just by
# lowering gamma to the anchovy's 750), but NOT sufficient to undo it
# starting from the anchovy (raising gamma to the default's own value
# barely moved the anchovy's own ~0.034% -> ~0.042%). Something else in the
# anchovy's parameter set is locking in the extreme skew regardless of
# gamma. alpha (assimilation efficiency) is the natural next suspect: it
# scales NET growth the same way gamma scales intake, and the anchovy's
# alpha = 0.1 is fairly low if mizer's own default is meaningfully higher.
# Same 2x2 isolation design as the gamma check, alpha instead of gamma, and
# reusing run_spectrum_check() unchanged since it already reports both.
################################################################################

default_alpha <- default_spectrum_params@species_params$alpha[1]
cat(sprintf("Default species' own alpha: %.4f (anchovy's is 0.1)\n", default_alpha))

# Anchovy species parameters, but alpha is swappable instead of hardcoded --
# gamma stays fixed at the anchovy's own 750 throughout, so this isolates
# alpha specifically rather than re-opening the gamma question.
make_anchovy_alpha_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                       second_order = TRUE, ext_diff = 0.00,
                                       lambda = 2.05, alpha_override) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = alpha_override, gamma = 750, ks = 0
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

# Default species parameters, but alpha is swappable instead of left to
# mizer's own default. gamma stays at the default species' own value
# throughout, for the same reason.
make_default_alpha_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                       second_order = TRUE, ext_diff = 0.00,
                                       lambda = NULL, alpha_override) {
  params <- newSingleSpeciesParams(species_name = "Anchovy", alpha = alpha_override)

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

alpha_isolation_df <- bind_rows(
  run_spectrum_check(make_second_order_params_kr,
                     label = "Anchovy params, anchovy alpha (0.1)"),
  run_spectrum_check(make_anchovy_alpha_swap_kr,
                     label = "Anchovy params, default alpha",
                     alpha_override = default_alpha),
  run_spectrum_check(make_default_species_params_kr,
                     label = "Default params, default alpha"),
  run_spectrum_check(make_default_alpha_swap_kr,
                     label = "Default params, anchovy alpha (0.1)",
                     alpha_override = 0.1)
)

print(alpha_isolation_df)
write.csv(alpha_isolation_df, file.path("interesting_plots", "alpha_isolation_df.csv"),
          row.names = FALSE)

cat(paste0(
  "Alpha isolation check (capacity_mult = 10, resource_decrease = 0.001, gamma held at\n",
  "each species' own value throughout): if alpha is the other piece of the puzzle, 'Anchovy\n",
  "params, default alpha' should now pull the anchovy's extreme fraction back up toward the\n",
  "default baseline (~7%), rescuing what the gamma swap alone couldn't. 'Default params,\n",
  "anchovy alpha' should push the default baseline back down toward the anchovy's extreme end.\n",
  "If neither swap moves much, alpha isn't it either, and ks or w_min/w_max are the remaining\n",
  "suspects.\n"
))

################################################################################
# ks and kappa: closing the remaining gap
#
# Alpha turned out to be a dominant lever (0.034% <-> 70.3%, 7.13% <->
# 8.0e-8%, 5+ orders of magnitude either way) -- much bigger than gamma's
# effect. But a genuine loose thread remains: two rows, one from the gamma
# test and one from the alpha test, share the EXACT same alpha = 0.1,
# gamma = 2066 combination, differing only in which species' remaining
# parameters (w_min, w_max, w_mat, kappa, ks) were used --
#   anchovy's other params + alpha=0.1, gamma=2066: 0.042%  (gamma test)
#   default's other params + alpha=0.1, gamma=2066: 8.0e-8% (alpha test)
# still a ~500,000x gap with alpha AND gamma pinned identical. Something in
# that remaining bucket is also doing real work. ks (standard metabolism)
# and kappa (background resource-spectrum coefficient) are the two
# parameters left that plausibly affect the growth/mortality rate balance
# directly -- w_min/w_max mostly set the size RANGE rather than the rate
# balance, so lower priority to test first.
################################################################################

default_ks <- default_spectrum_params@species_params$ks[1]

# kappa's storage location after construction has now been guessed wrong
# twice (species_params$kappa, then resource_params()$kappa) -- both
# silently returned NULL rather than erroring. But the likely reason isn't
# a third wrong accessor name: setResource(..., resource_dynamics =
# "resource_semichemostat", ...) runs AFTER newSingleSpeciesParams() in
# make_default_species_params_kr(), and very plausibly replaces the
# power-law (kappa/lambda) resource setup with the semichemostat one
# (resource_rate/resource_capacity) entirely -- by the time
# default_spectrum_params exists, kappa may simply no longer be part of the
# object, not just stored somewhere unfound. Querying a FRESH default
# species object instead, before that override ever runs, is the fix that
# actually addresses this rather than trying a third accessor name blind.
get_kappa <- function(p) {
  v <- tryCatch(resource_params(p)$kappa, error = function(e) NULL)
  if (is.null(v) || length(v) == 0) NA_real_ else v[1]
}

default_species_raw <- newSingleSpeciesParams(species_name = "Anchovy")
default_kappa <- get_kappa(default_species_raw)
if (is.na(default_kappa)) {
  warning("Still couldn't retrieve kappa, even from a fresh pre-setResource() default species object -- the accessor name itself needs finding by hand this time, e.g. via `str(default_species_raw, max.level = 2)` or `slotNames(default_species_raw)`. The kappa isolation test's 'default kappa' rows will fail cleanly rather than silently testing the wrong value.")
}
anchovy_kappa <- 100 * exp(-6.9 * (2.05 - 1))  # same formula make_second_order_params_kr() uses
cat(sprintf(
  "Default species' own ks: %.4f (anchovy's is 0), kappa: %s (anchovy's is %.4f)\n",
  default_ks, ifelse(is.na(default_kappa), "NOT FOUND", sprintf("%.4f", default_kappa)), anchovy_kappa
))

make_anchovy_ks_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                    second_order = TRUE, ext_diff = 0.00,
                                    lambda = 2.05, ks_override) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = ks_override
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

make_default_ks_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                    second_order = TRUE, ext_diff = 0.00,
                                    lambda = NULL, ks_override) {
  params <- newSingleSpeciesParams(species_name = "Anchovy", ks = ks_override)

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

ks_isolation_df <- bind_rows(
  run_spectrum_check(make_second_order_params_kr,
                     label = "Anchovy params, anchovy ks (0)"),
  run_spectrum_check(make_anchovy_ks_swap_kr,
                     label = "Anchovy params, default ks",
                     ks_override = default_ks),
  run_spectrum_check(make_default_species_params_kr,
                     label = "Default params, default ks"),
  run_spectrum_check(make_default_ks_swap_kr,
                     label = "Default params, anchovy ks (0)",
                     ks_override = 0)
)

print(ks_isolation_df)
write.csv(ks_isolation_df, file.path("interesting_plots", "ks_isolation_df.csv"),
          row.names = FALSE)

make_anchovy_kappa_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                       second_order = TRUE, ext_diff = 0.00,
                                       lambda = 2.05, kappa_override) {
  no_w <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa_override,
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

make_default_kappa_swap_kr <- function(resource_decrease = 0.001, capacity_mult = 1,
                                       second_order = TRUE, ext_diff = 0.00,
                                       lambda = NULL, kappa_override) {
  params <- newSingleSpeciesParams(species_name = "Anchovy", kappa = kappa_override)

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

kappa_isolation_df <- bind_rows(
  run_spectrum_check(make_second_order_params_kr,
                     label = "Anchovy params, anchovy kappa"),
  run_spectrum_check(make_anchovy_kappa_swap_kr,
                     label = "Anchovy params, default kappa",
                     kappa_override = default_kappa),
  run_spectrum_check(make_default_species_params_kr,
                     label = "Default params, default kappa"),
  run_spectrum_check(make_default_kappa_swap_kr,
                     label = "Default params, anchovy kappa",
                     kappa_override = anchovy_kappa)
)

print(kappa_isolation_df)
write.csv(kappa_isolation_df, file.path("interesting_plots", "kappa_isolation_df.csv"),
          row.names = FALSE)

cat(paste0(
  "ks/kappa isolation check (capacity_mult = 10, resource_decrease = 0.001, alpha and gamma\n",
  "held at each species' own natural value throughout): if either ks or kappa alone closes\n",
  "the remaining ~500,000x gap found between the two alpha=0.1,gamma=2066 rows above, that's\n",
  "the last piece. If neither moves the needle much, the w_min/w_max size-range span itself\n",
  "is the remaining suspect.\n"
))
