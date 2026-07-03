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

# Day 20 picks up two threads left open by Day 18/19:
#   Part 1-2: the forward/backward resource_decrease sweep left a gap at
#     rd ~ 0.008 that couldn't be written off as a sweep-start detection
#     artefact -- plot max and min separately (not collapsed to
#     max - min) to actually see the shape of that gap, then try to pin down
#     whether it's real hysteresis or a run-length/resolution artefact.
#   Part 3: recommendations (diagnostic and experimental) for investigating
#     why so much biomass piles up at juvenile sizes, continuing the
#     unfinished survivorship-through-growth thread at the bottom of
#     19_experiments.R.

# Helper functions, carried over from Day 18/19. balance defaults to TRUE here
# (unlike 19_experiments.R's FALSE) because every experiment below reuses
# Day 18's original oscillating regime, where balance = TRUE.
make_params <- function(lambda = 2.05, resource_decrease = 0.001,
                        beta = NULL, sigma = NULL, balance = TRUE) {
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

make_limit_cycle_sim <- function(params, t_total = 600, effort = 0, perturbation = 1e3) {
  params@initial_n_pp[] <- params@cc_pp * 0.1

  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "tr_bdf2")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / perturbation

  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "tr_bdf2")
}

################################################################################
# Part 1: Max and min plotted separately, forward vs backward, instead of
# collapsing to a single amp = max - min number. Day 18's amp-only plot could
# not tell apart "the whole cycle sits differently depending on direction"
# (real hysteresis) from "one branch drifts while the other tracks tightly"
# (an asymmetric artefact, e.g. critical slowing down dragging the mean level
# around). run_rd_sweep is extended to record max_bm/min_bm, not just amp.
################################################################################

run_rd_sweep <- function(rd_seq, init_n = NULL, init_n_pp = NULL,
                         t_run = 300, lambda = 2.05, label = "") {
  out       <- data.frame(rd = rd_seq, max_bm = NA_real_, min_bm = NA_real_,
                          bm_mean = NA_real_)
  state_n   <- init_n
  state_npp <- init_n_pp

  for (i in seq_along(rd_seq)) {
    p <- make_params(lambda = lambda, resource_decrease = rd_seq[i])
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

    out$max_bm[i]  <- max(late)
    out$min_bm[i]  <- min(late)
    out$bm_mean[i] <- mean(late)
  }
  list(df = out, n_final = state_n, npp_final = state_npp)
}

rd_seq <- seq(0.001, 0.02, length.out = 20)  # same grid as Day 18

fwd <- run_rd_sweep(rd_seq, label = "FWD")
bwd <- run_rd_sweep(rev(rd_seq),
                    init_n    = fwd$n_final,
                    init_n_pp = fwd$npp_final,
                    label     = "BWD")

fwd_df <- fwd$df
bwd_df <- bwd$df[order(bwd$df$rd), ]

# Base-R overlay: four curves (fwd max, fwd min, bwd max, bwd min) on the same
# axes. Below the transition, max and min should collapse onto each other
# (the stable fixed point); above it they fan out into two branches (the
# cycle's peak and trough). Where forward and backward branches fail to
# overlap is the real hysteresis region -- not just a difference in a single
# collapsed "amp" number.
plot(fwd_df$rd, fwd_df$max_bm, type = "l", col = "steelblue", lwd = 2, lty = 1,
     xlab = "resource_decrease", ylab = "Biomass",
     main = "Hysteresis: max and min plotted separately",
     ylim = range(c(fwd_df$max_bm, fwd_df$min_bm, bwd_df$max_bm, bwd_df$min_bm)))
lines(fwd_df$rd, fwd_df$min_bm, col = "steelblue", lwd = 2, lty = 2)
lines(bwd_df$rd, bwd_df$max_bm, col = "firebrick", lwd = 2, lty = 1)
lines(bwd_df$rd, bwd_df$min_bm, col = "firebrick", lwd = 2, lty = 2)
legend("topright",
       legend = c("Forward max", "Forward min", "Backward max", "Backward min"),
       col = c("steelblue", "steelblue", "firebrick", "firebrick"),
       lty = c(1, 2, 1, 2), lwd = 2, cex = 0.8)

# Same data, tidier bifurcation-diagram-style ggplot version -- long format,
# one line per (direction, branch) combination.
bifurcation_df <- bind_rows(
  data.frame(rd = fwd_df$rd, biomass = fwd_df$max_bm, direction = "Forward", branch = "max"),
  data.frame(rd = fwd_df$rd, biomass = fwd_df$min_bm, direction = "Forward", branch = "min"),
  data.frame(rd = bwd_df$rd, biomass = bwd_df$max_bm, direction = "Backward", branch = "max"),
  data.frame(rd = bwd_df$rd, biomass = bwd_df$min_bm, direction = "Backward", branch = "min")
)

ggplot(bifurcation_df, aes(x = rd, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(x = "resource_decrease", y = "Biomass",
       title = "Bifurcation diagram: max/min branches, forward vs backward",
       subtitle = "Branches collapsing onto one curve = fixed point; fanning apart = limit cycle") +
  theme_minimal()

################################################################################
# Part 2: Is the gap real hysteresis, or a run-length/resolution artefact?
# Two follow-ups, both flagged in Day 18's revised What's Next.
################################################################################

# 2a. Finer rd grid, zoomed on the region where the backward curve was still
# oscillating (rd ~ 0.008) while the forward curve looked already flat. If the
# gap between forward and backward transition points shrinks to near-nothing
# at this resolution, the original 10-point sweep was just undersampling a
# sharp but continuous transition. If a clear, stable gap survives, that's
# much harder to explain away as sweep resolution.
rd_fine <- seq(0.005, 0.012, length.out = 15)

fwd_fine <- run_rd_sweep(rd_fine, t_run = 300, label = "FWD fine")
bwd_fine <- run_rd_sweep(rev(rd_fine),
                         init_n    = fwd_fine$n_final,
                         init_n_pp = fwd_fine$npp_final,
                         t_run     = 300, label = "BWD fine")

fwd_fine_df <- fwd_fine$df
bwd_fine_df <- bwd_fine$df[order(bwd_fine$df$rd), ]

plot(fwd_fine_df$rd, fwd_fine_df$max_bm - fwd_fine_df$min_bm, type = "b",
     col = "steelblue", lwd = 2, xlab = "resource_decrease",
     ylab = "Amplitude (max - min)",
     main = "Fine-grid hysteresis check (rd = 0.005-0.012)",
     ylim = c(0, max(fwd_fine_df$max_bm - fwd_fine_df$min_bm,
                     bwd_fine_df$max_bm - bwd_fine_df$min_bm) * 1.1))
lines(bwd_fine_df$rd, bwd_fine_df$max_bm - bwd_fine_df$min_bm, type = "b",
      col = "firebrick", lwd = 2)
legend("topright", legend = c("Forward", "Backward"),
       col = c("steelblue", "firebrick"), lwd = 2)

# 2b. Much longer per-step runs at the original grid, to let critical slowing
# down fully resolve before amplitude is measured -- tests Day 18's other
# candidate explanation (run length, not sweep resolution) directly. If the
# gap closes up at 5x the run length, critical slowing down was the whole
# explanation. If it survives, that's stronger evidence for a genuinely
# subcritical bifurcation with a real, if narrow, bistable window.
fwd_long <- run_rd_sweep(rd_seq, t_run = 1500, label = "FWD long")
bwd_long <- run_rd_sweep(rev(rd_seq),
                         init_n    = fwd_long$n_final,
                         init_n_pp = fwd_long$npp_final,
                         t_run     = 1500, label = "BWD long")

fwd_long_df <- fwd_long$df
bwd_long_df <- bwd_long$df[order(bwd_long$df$rd), ]

plot(fwd_long_df$rd, fwd_long_df$max_bm - fwd_long_df$min_bm, type = "b",
     col = "steelblue", lwd = 2, xlab = "resource_decrease",
     ylab = "Amplitude (max - min)",
     main = "Hysteresis check with t_run = 1500 (vs Day 18's 300)")
lines(bwd_long_df$rd, bwd_long_df$max_bm - bwd_long_df$min_bm, type = "b",
      col = "firebrick", lwd = 2)
legend("topright", legend = c("Forward", "Backward"),
       col = c("steelblue", "firebrick"), lwd = 2)

################################################################################
# Part 3: Why so much biomass at juvenile sizes? Continuing the unfinished
# thread at the bottom of 19_experiments.R (the getMort() decomposition that
# was started but never interpreted). Rebuilds the same states so this script
# runs standalone.
################################################################################

r_d           <- exp(seq(log(0.001), log(0.3), length.out = 10))
day18_sim     <- make_limit_cycle_sim(make_params(resource_decrease = r_d[1], balance = TRUE))
sim_undepleted <- make_limit_cycle_sim(make_params(resource_decrease = 1, balance = TRUE))

# Evaluate growth/mortality at the model's own dynamic resource state (not a
# freshly built params object's default initial_n_pp -- that was the bug that
# made the first survivorship attempt in 19_experiments.R come back
# bit-identical for both scenarios).
params_at_sim_state <- function(sim, t_min = 550) {
  p    <- sim@params
  tvec <- as.numeric(dimnames(sim@n)[[1]])
  idx  <- which(tvec > t_min)

  p@initial_n[]    <- apply(sim@n[idx, , , drop = FALSE], c(2, 3), mean)
  p@initial_n_pp[] <- apply(sim@n_pp[idx, , drop = FALSE], 2, mean)
  p
}

p_depleted_dynamic   <- params_at_sim_state(day18_sim,      t_min = 550)
p_undepleted_dynamic <- params_at_sim_state(sim_undepleted, t_min = 550)

# --- Diagnostic recommendation 1: finish the getMort() decomposition ---
# There's no cannibalism in this single-species sim, so total mortality mu(w)
# is background/senescence mortality plus fishing (zero in the unfished
# baseline) plus predation mortality from... nothing, since there's no other
# species or size class preying on this one. mizer stores background
# mortality on params@mu_b and exposes predation mortality via
# getPredMort() -- VERIFY these against the installed mizer version
# (slotNames(p), ?MizerParams) before trusting the split below; if mu_b
# doesn't exist, getMort() minus getPredMort() minus getFMort() should give
# the same background component.
mortality_decomposition <- function(params, label) {
  w         <- params@w
  mu_total  <- getMort(params)[1, ]
  mu_pred   <- tryCatch(getPredMort(params)[1, ], error = function(e) rep(NA_real_, length(w)))
  mu_back   <- if (!is.null(params@mu_b)) params@mu_b[1, ] else rep(NA_real_, length(w))
  data.frame(w = w, label = label, mu_total = mu_total, mu_pred = mu_pred, mu_back = mu_back)
}

mort_decomp <- bind_rows(
  mortality_decomposition(p_depleted_dynamic,   "resource_decrease = 0.001 (dynamic)"),
  mortality_decomposition(p_undepleted_dynamic, "resource_decrease = 1 (dynamic)")
)

ggplot(mort_decomp, aes(x = w, y = mu_total, color = label)) +
  geom_line(linewidth = 1) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Body weight (g)", y = "Total mortality rate",
       title = "Total mortality mu(w): depleted vs undepleted resource state") +
  theme_minimal()
# If mu_total is nearly identical between the two scenarios, mortality isn't
# what's driving any survivorship difference -- it has to be growth rate
# (via dt_bin = dw/g below) instead. If mu_total is visibly higher for the
# depleted scenario, especially at small w, that's a direct, decomposed
# confirmation of a mortality-driven mechanism rather than an inference from
# the collapsed survivorship curve alone.

# --- Diagnostic recommendation 2: time-to-maturity, not just survivorship ---
# Survivorship folds growth and mortality together into one curve. This
# isolates the growth side alone: how much calendar time does an individual
# spend growing from w_min to w_mat? A demographic "residence time" effect
# (slow growth -> long residence -> lots of standing biomass sitting at small
# sizes) would inflate the juvenile:adult ratio even if survivorship were
# perfect (mu = 0 everywhere), which the survivorship curve alone can't show.
time_to_maturity <- function(params, w_mat = 10) {
  w      <- params@w
  dw     <- params@dw
  g      <- getEGrowth(params)[1, ]
  dt_bin <- dw / g
  idx    <- w <= w_mat
  sum(dt_bin[idx])
}

data.frame(
  scenario         = c("resource_decrease = 0.001 (dynamic)", "resource_decrease = 1 (dynamic)"),
  years_to_w_mat   = c(time_to_maturity(p_depleted_dynamic), time_to_maturity(p_undepleted_dynamic))
)
# If years_to_w_mat is much larger for the depleted/oscillating scenario,
# that's direct evidence individuals simply spend far longer exposed as
# juveniles under resource depletion, independent of whether mortality also
# differs -- a "slow clock" effect on top of (or instead of) a "high filter"
# mortality effect.

# --- Further diagnostic recommendations (not run here) ---
# 3. Plot getEGrowth(p_depleted_dynamic) and getEGrowth(p_undepleted_dynamic)
#    directly on log-log axes, overlaid. A growth bottleneck at small w (a
#    dip or flattening in g(w) well before w_mat) would show up directly here
#    and would explain an inflated dt_bin in that region rather than requiring
#    it to be inferred from the survivorship curve's slope.
# 4. Check the maturity ogive / reproductive allocation (species_params$w_mat,
#    w_mat25, or the psi vector via allocation functions such as
#    getReproductionAllocation()) -- if energy is diverted from somatic growth
#    into reproduction gradually below w_mat rather than as a hard step at it,
#    that alone would slow the approach to w_mat and inflate the
#    juvenile:adult standing-stock ratio, independent of any mortality effect.
# 5. Recompute the juvenile:adult biomass ratio from a genuinely steady,
#    non-oscillating run (resource_decrease = 1, already computed above as
#    sim_undepleted) and compare its order of magnitude to the oscillating
#    run's ratio. If the undepleted/steady case shows a similarly lopsided
#    ratio, the juvenile pile-up is a generic feature of this size-spectrum
#    parameterisation (alpha, gamma, kappa, lambda), not something specific to
#    the resource-depletion/oscillation regime being studied here.

# --- Experimental recommendations (not run here; each needs new simulation) ---
# 6. Vary alpha (assimilation efficiency, currently fixed at 0.1) across a
#    small grid (e.g. 0.1, 0.2, 0.4) at fixed resource_decrease, and track the
#    juvenile:adult biomass ratio (via track_life_stages, from
#    19_experiments.R) at each. alpha directly sets growth rate g(w); if the
#    ratio is highly sensitive to alpha, slow growth (not mortality) is
#    confirmed as the dominant lever.
# 7. Vary background mortality (params@mu_b, or whichever species_params
#    column controls it in the installed mizer version) holding alpha fixed,
#    as the mortality-side counterpart to recommendation 6. A full alpha x
#    mu_b factorial would cleanly separate "slow growth" from "high background
#    death rate" as the dominant driver, rather than inferring it post hoc
#    from the mu_total/dt_bin decomposition above.
# 8. Track a single cohort explicitly instead of reading the standing size
#    spectrum: seed initial_n with all biomass concentrated in the smallest
#    size bin at t = 0, project with no further perturbation, and watch how
#    that pulse spreads and thins as it moves up the size axis (a "cohort in
#    a jar" experiment). This is the most direct experimental test of the
#    survivorship-through-growth mechanism, and would show visually where the
#    die-off or slowdown actually happens rather than relying on rate
#    functions evaluated at a single point in time.
# 9. Compare against one of mizer's built-in multi-species community models
#    (e.g. the North Sea example params) to check whether a juvenile-heavy
#    size spectrum is a generic feature of size-spectrum theory (the
#    allometric scaling of growth/mortality that produces roughly a -2 to -1
#    biomass-density slope) rather than something peculiar to this specific
#    single-species anchovy calibration (alpha = 0.1, gamma = 750,
#    kappa/lambda choice). If a well-validated, independently parameterised
#    model shows the same lopsidedness, that points at a structural feature of
#    size-spectrum models generally rather than a calibration artefact here.
