library(mizer)
library(mizerExperimental)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)

# Day 39 takes Day 38's cannibalism model -- the Figure 2e plankton-anchovy
# model (Canales, Delius & Law 2020), cannibalism only, mu_l=0 -- and asks
# this of it:
#
#   Section 1: cut its dependence on cannibalism, to make it MORE sensitive
#   to fishing.
#   Section 2: cut cannibalism the same way again, but this time raise its
#   interaction with the resource too, because Section 1 alone risks
#   dropping growth rates too far. Once the script actually ran (see
#   Section 2's own header), this turned out to be badly insufficient --
#   growth at w>=w_mat stayed crushed to roughly a seventh of its original
#   value, not "restored".
#   Section 3: an actual sweep of theta x interaction_resource together
#   (added after reading Section 1/2's own saved results back), to see
#   what combination, if any, actually restores growth without giving back
#   all of Section 1's fishing sensitivity.
#   Section 4: test whether the result actually depends on fishing effort,
#   sweeping at the default knife-edge size (10), a slightly lower one (9),
#   and a much lower one (5). A first pass at this compared "Constant"
#   fishing against nothing at all, which conflates timing with average
#   effort -- fixed here (see Section 4's own header) into three schedules
#   that all share the same Constant floor, differing only in whether extra
#   effort gets layered on near cycle peaks, near troughs, or not at all.
#   Section 5: summary.
#
# Self-contained convention since Day 20: helpers redefined here, not
# sourced from 38_experiments.R.

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH truncation guard, carried over from Day 30 onward.
save_plot <- function(plot, filename, width = 9, height = 6, dpi = 150) {
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
# Section 0: rebuilding Day 38's Figure 2e plankton-anchovy model, plus
# fishing -- verbatim from 38_experiments.R (p2 list, setAnchovyMort(),
# plankton_logistic(), norm_box_pred_kernel(), setAnchovyModel(),
# make_anchovy_fishing_params()). Day 38 already checked the kernel and
# mortality curve this builds match the paper (day37_anchovy_kernel_check.png,
# day37_anchovy_mort_check.png) and that this exact config sustains a limit
# cycle where Day 37's other two candidates don't -- not repeated here.
################################################################################

p2 <- list(
  dt = 0.001, dx = 0.1, w_min = 0.0003, w_inf = 66.5,
  ppmr_min = 100, ppmr_max = 30000, gamma = 750, alpha = 0.85, K = 0.1,
  mu_l = 0, w_l = 0.03, rho_l = 5,          # mu_l=0 fixed: Figure 2e config
  mu_0 = 1, rho_b = -0.25,
  w_s = 0.5, rho_s = 1,
  w_mat = 10, rho_m = 15, rho_inf = 0.2, epsilon_R = 0.1,
  w_pp_cutoff = 0.1, r0 = 10, a0 = 100, i0 = 100, rho = 0.85, lambda = 2
)

setAnchovyMort <- function(params, p) {
  mu_b <- rep(0, length(params@w))
  mu_b[params@w <= p$w_s] <-
    (p$mu_0 * (params@w / p$w_min)^p$rho_b)[params@w < p$w_s]
  if (p$mu_0 > 0) {
    mu_s <- min(mu_b[params@w <= p$w_s])
  } else {
    mu_s <- p$mu_s
  }
  mu_b[params@w >= p$w_s] <-
    (mu_s * (params@w / p$w_s)^p$rho_s)[params@w >= p$w_s]
  mu_b <- mu_b + p$mu_l / (1 + (params@w / p$w_l)^p$rho_l)

  params@mu_b[] <- mu_b
  return(params)
}

plankton_state <- new.env(parent = emptyenv())
plankton_state$time   <- 0
plankton_state$factor <- 1
plankton_state$random <- FALSE   # random plankton forcing kept off, as Day 37/38

plankton_logistic <- function(params, n, n_pp, n_other, rates, dt = 0.1, ...) {
  plankton_state$time <- plankton_state$time + dt
  f <- params@rr_pp * n_pp * (1 - n_pp / params@cc_pp / plankton_state$factor) +
    anchovy_immigration - rates$resource_mort * n_pp
  f[is.na(f)] <- 0
  return(n_pp + dt * f)
}

norm_box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
  phi <- rep(1, length(ppmr))
  phi[ppmr > ppmr_max] <- 0
  phi[ppmr < ppmr_min] <- 0
  phi[1] <- 0
  logppmr <- log(ppmr)
  dl <- logppmr[2] - logppmr[1]
  N <- sum(phi) * dl
  phi / N
}

setAnchovyModel <- function(p) {
  kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))

  species_params <- data.frame(
    species          = "Anchovy",
    w_min            = p$w_min,
    w_mat            = p$w_mat,
    m                = p$rho_inf + 2 / 3,
    w_inf            = p$w_inf,
    erepro           = p$epsilon_R,
    alpha            = p$K,
    ks               = 0,
    gamma            = p$gamma,
    q                = p$alpha,
    ppmr_min         = p$ppmr_min,
    ppmr_max         = p$ppmr_max,
    pred_kernel_type = "norm_box",
    h                = Inf,
    R_max            = Inf,
    linecolour       = "brown",
    stringsAsFactors = FALSE
  )

  no_w <- round(log(p$w_inf / p$w_min) / p$dx)

  params <- newMultispeciesParams(
    species_params,
    no_w = no_w,
    lambda = p$lambda,
    kappa = kappa,
    w_pp_cutoff = p$w_pp_cutoff,
    resource_dynamics = "plankton_logistic"
  )

  params@rr_pp[] <- p$r0 * params@w_full^(p$rho - 1)
  params
}

# Figure 2e config + knife-edge gear at w_mat, catchability=1 -- Day 36/38's
# own make_anchovy_fishing_params() convention, carried over.
make_anchovy_fishing_params <- function(p = p2) {
  params <- setAnchovyModel(p)
  params <- setAnchovyMort(params, p)
  params@interaction[] <- 1

  gp                 <- params@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- params@species_params$w_mat
  gp$catchability    <- 1
  gear_params(params) <- gp

  params
}

p_scan <- make_anchovy_fishing_params()
anchovy_immigration <- p2$i0 * p_scan@w_full^(-p2$lambda) * exp(-6.9 * (p2$lambda - 1))

cat(sprintf("Figure 2e plankton-anchovy erepro: %.4g (should sit in [0,1])\n",
           species_params(p_scan)$erepro))

################################################################################
# Section 0b: building and calibrating the fork
#
# Day 37/38's own settle+kick recipe (10yr settle from a power-law abundance,
# 10^7 knockdown across the whole grid at t=10), continued unfished up to
# t_fork=20 -- Day 37 found the cycle is already at full amplitude well
# before then. Confirmed below rather than assumed.
################################################################################

t_fork               <- 20
scan_post_fork_years <- 24   # ~4-6 periods of the paper's own ~6yr cycle
scan_summary_window  <- 12   # ~2 periods, matching Day 37/38's own sampling window
scan_t_cut           <- t_fork + scan_post_fork_years - scan_summary_window

make_anchovy_fork_sim <- function(params, t_fork) {
  params@initial_n[]    <- 0.001 * params@w^(-1.8)
  params@initial_n_pp[] <- params@cc_pp
  sim <- project(params, t_max = 10, dt = p2$dt, progress_bar = FALSE)
  sim@n[11, , ] <- sim@n[11, , ] / 10^7
  project(sim, t_max = t_fork - 10, dt = p2$dt, t_save = 0.2, progress_bar = FALSE, effort = 0)
}

compute_selected_biomass_series <- function(sim, params, t_cut) {
  tv    <- as.numeric(dimnames(sim@n)[[1]])
  keep  <- which(tv > t_cut)
  f_ref <- mizerFMortGear(params, effort = 1)

  vapply(keep, function(i) {
    n_i <- array(sim@n[i, , , drop = FALSE], dim = dim(sim@n)[-1])
    biomass_density <- sweep(n_i, 2, params@w * params@dw, "*")
    sum(colSums(f_ref) * biomass_density)
  }, numeric(1))
}

sim_fork <- make_anchovy_fork_sim(p_scan, t_fork)
last_n   <- array(sim_fork@n[dim(sim_fork@n)[1], , , drop = FALSE], dim = dim(sim_fork@n)[-1])
last_npp <- sim_fork@n_pp[dim(sim_fork@n_pp)[1], ]

fork_bp <- compute_selected_biomass_series(sim_fork, p_scan, t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "Fork check (t in [%.0f,%.0f], last %.0fyr of the %.0fyr fork): selected biomass min=%.4g max=%.4g -- should already span close to the full cycle.\n",
  t_fork - scan_summary_window, t_fork, scan_summary_window, t_fork, min(fork_bp), max(fork_bp)
))

################################################################################
# Section 1: decreasing the model's dependence on cannibalism, to make it
# more sensitive to fishing
#
# getEncounter()'s own formula (checked via ?getEncounter this session) is
#   E_i(w) = gamma_i(w) * integral[ theta_ip*N_R(w_p) + sum_j theta_ij*N_j(w_p) ] * phi(w,w_p) w_p dw_p
# so for this single-species model theta_11 (the one entry of @interaction,
# set to 1 by make_anchovy_fishing_params()) IS the cannibalism knob -- how
# strongly an Anchovy's own encounter rate (both food, and via getPredMort()
# the mortality it inflicts on smaller conspecifics) counts other Anchovy as
# prey. theta_1p (species_params$interaction_resource, default 1, untouched
# here) is the separate resource-interaction knob Section 2 uses.
# setInteraction()'s own help page confirms these are two distinct levers,
# not one -- interaction_resource is not part of the interaction matrix.
#
# Why lowering theta should make the model MORE fishing sensitive: Day 38
# Section 3c found this model's whole resistance to fishing (the
# "cultivation effect") comes from cannibalistic adults suppressing
# juveniles -- fish the adults, juveniles gain more than the fishery
# removes. That relief only exists because cannibalism is doing real work in
# the first place. Turn cannibalism down and there is less of that relief
# for fishing to trigger, so the population's response to effort should owe
# less to the cultivation effect and more to straightforward biomass
# removal.
#
# theta_low=0.3 (70% cut): an early interactive check, done BEFORE this file
# first ran end-to-end, evaluated growth at theta=0 vs. theta=1 while
# holding the population state fixed at the theta=1 baseline's own cycling
# abundance -- a snapshot of the DIRECT diet-share effect only. It found
# conspecifics are under 0.1% of an Anchovy's own diet at every size
# checked (kappa's own resource spectrum dwarfs the single-species fish
# spectrum here), and concluded growth "barely moves with theta". That
# conclusion was WRONG, or at least badly incomplete -- caught only once
# the script actually ran end-to-end and growth_compare_df below could be
# read back against real numbers, not predicted:
#
#   - at the two SMALL sizes (w=0.01, 0.1g), growth is fine, even improved
#     (+10%, +19%) -- consistent with the direct diet-share check above.
#   - at w=1g growth is UP substantially (+70%).
#   - but at the two LARGE sizes (w=10.1g, 28.6g -- at and above w_mat),
#     growth CRASHES: -88% and -89% versus the theta=1 baseline.
#
# The direct diet-share check wasn't wrong about the mechanism it measured
# -- it's genuinely true that cannibalism is a negligible direct food
# source here. What it missed is the INDIRECT, population-mediated effect:
# theta=0.3 removes most of the predation mortality that normally caps the
# small/juvenile size classes, so that population booms (Section 1's own
# effort-probe below shows total biomass roughly 50% higher than baseline
# at every fish_level tested). A much bigger population of small-to-mid
# fish grazes the SAME shared resource spectrum that large fish depend on
# for nearly all of their own food -- resource competition between size
# classes, mediated through n_pp, not a direct cannibalism/growth link.
# Fixing this by holding the population fixed and only perturbing theta
# (the first check's own method) structurally cannot see this, because the
# whole effect runs through letting the population size-structure actually
# change.
#
# predation MORTALITY (getPredMort()) is similarly not the clean
# linear-in-theta story a fixed-population check would suggest: at
# w=0.01g it roughly halves (51.9 -> 27.3/yr, not the ~3.3x drop naive
# theta-scaling alone would predict), because the much larger population at
# theta=0.3 partly offsets the lower per-capita theta. At w=0.1g it goes
# the OTHER way -- UP by ~9x (0.49 -> 4.54/yr) -- the population-size effect
# dominates the per-capita theta cut entirely at that size. Turning
# cannibalism's per-capita strength down does not mean turning its
# population-level consequences down uniformly.
#
# That growth check matters directly for Section 2: the "too-low growth"
# risk is not a small, mostly-theoretical concern here -- it is the
# dominant effect at exactly the sizes (w>=w_mat) this whole project's
# fishing analysis is built around. Section 2's fix needs to be judged
# against an ~88% growth deficit at w=10g, not against "barely moves".
################################################################################

make_anchovy_fishing_params_theta <- function(interaction_val, interaction_resource_val,
                                              knife_edge_size = p2$w_mat, p = p2) {
  params <- setAnchovyModel(p)
  params <- setAnchovyMort(params, p)
  params@interaction[] <- interaction_val

  sp <- species_params(params)
  sp$interaction_resource <- interaction_resource_val
  species_params(params) <- sp

  gp                 <- params@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- knife_edge_size
  gp$catchability    <- 1
  gear_params(params) <- gp

  params
}

theta_low <- 0.3

p_theta_low <- make_anchovy_fishing_params_theta(theta_low, 1, knife_edge_size = p2$w_mat)
sim_fork_theta_low <- make_anchovy_fork_sim(p_theta_low, t_fork)
last_n_theta_low   <- array(sim_fork_theta_low@n[dim(sim_fork_theta_low@n)[1], , , drop = FALSE],
                            dim = dim(sim_fork_theta_low@n)[-1])
last_npp_theta_low <- sim_fork_theta_low@n_pp[dim(sim_fork_theta_low@n_pp)[1], ]

bp_theta_low <- compute_selected_biomass_series(sim_fork_theta_low, p_theta_low,
                                                t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "Section 1 (theta=%.2g, interaction_resource=1), unfished fork, last %.0fyr of %.0fyr: selected biomass min=%.4g max=%.4g -- still cycling.\n",
  theta_low, scan_summary_window, t_fork, min(bp_theta_low), max(bp_theta_low)
))

# Growth vs. predation-mortality decomposition, baseline (theta=1) vs.
# theta_low, both evaluated on their OWN cycling fork state.
w_check   <- c(0.01, 0.1, 1, 10, 30)
idx_check <- vapply(w_check, function(x) which.min(abs(p_scan@w - x)), integer(1))

growth_compare_df <- data.frame(
  w                  = round(p_scan@w[idx_check], 4),
  growth_baseline    = getEGrowth(p_scan, n = last_n, n_pp = last_npp)[1, idx_check],
  growth_theta_low   = getEGrowth(p_theta_low, n = last_n_theta_low, n_pp = last_npp_theta_low)[1, idx_check],
  predmort_baseline  = getPredMort(p_scan, n = last_n, n_pp = last_npp)[1, idx_check],
  predmort_theta_low = getPredMort(p_theta_low, n = last_n_theta_low, n_pp = last_npp_theta_low)[1, idx_check]
)
write.csv(growth_compare_df, file.path("interesting_plots", "day39_theta_low_growth_mort.csv"), row.names = FALSE)
cat("Section 1 growth/mortality comparison (theta=1 vs. theta=0.3), by size -- see day39_theta_low_growth_mort.csv:\n")
print(growth_compare_df)

# Constant-effort probe. Originally just fish_level=[1,3,5,7,9] (this
# project's own established range, Day 34/36/38 Sections 2-5) -- extended
# here up to 100 to match Section 4's own scale, since [1,9] alone showed
# almost no curvature (Day 38's own finding: this whole model class needs
# far more effort than [1,9] to register anything). Gear left at
# knife_edge=10 (w_mat, the default) throughout, so this isolates the
# effect of theta/interaction_resource alone, not the gear-cutoff effect
# Section 4 tests separately.
fish_level_seq_theta <- c(1, 3, 5, 7, 9, 20, 30, 50, 75, 100)

cannibalism_fishing_probe <- function(params, seed_n, seed_npp, fish_level_seq,
                                      years = scan_summary_window) {
  params@initial_n[]    <- seed_n
  params@initial_n_pp[] <- seed_npp
  bind_rows(lapply(fish_level_seq, function(fl) {
    sim <- project(params, t_max = years, dt = p2$dt, t_save = 0.2,
                   t_start = t_fork, progress_bar = FALSE, effort = fl)
    tv   <- as.numeric(dimnames(sim@n)[[1]])
    keep <- tv > t_fork + years / 2
    total    <- unname(getBiomass(sim)[, "Anchovy"])[keep]
    selected <- unname(getBiomass(sim, min_w = p2$w_mat)[, "Anchovy"])[keep]
    yield    <- unname(getYield(sim)[, "Anchovy"])[keep]
    data.frame(fish_level = fl, mean_total = mean(total), min_total = min(total),
              mean_selected = mean(selected), mean_juvenile = mean(total - selected),
              mean_yield = mean(yield))
  }))
}

baseline_probe_s1  <- cannibalism_fishing_probe(p_scan, last_n, last_npp, fish_level_seq_theta) %>%
  mutate(model = "Baseline (theta=1)")
theta_low_probe_s1 <- cannibalism_fishing_probe(p_theta_low, last_n_theta_low, last_npp_theta_low,
                                                fish_level_seq_theta) %>%
  mutate(model = "theta=0.3")

theta_probe_df <- bind_rows(baseline_probe_s1, theta_low_probe_s1)
write.csv(theta_probe_df, file.path("interesting_plots", "day39_theta_low_effort_probe.csv"), row.names = FALSE)

theta_probe_plot <- ggplot(theta_probe_df, aes(x = fish_level, y = mean_total, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=10)", y = "Mean total biomass",
       title = "Does cutting cannibalism (theta) flatten the model's own fishing resistance?",
       subtitle = "Gear unchanged (knife_edge=10) -- isolates theta's effect from Section 4's gear-cutoff effect") +
  theme_minimal()
theta_probe_plot
save_plot(theta_probe_plot, "day39_theta_low_effort_probe.png", width = 9, height = 6)

theta_probe_yield_plot <- ggplot(theta_probe_df, aes(x = fish_level, y = mean_yield, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=10)", y = "Mean yield",
       title = "Yield vs. effort -- does cutting cannibalism change where yield peaks or turns over?",
       subtitle = "Gear unchanged (knife_edge=10) -- isolates theta's effect from Section 4's gear-cutoff effect") +
  theme_minimal()
theta_probe_yield_plot
save_plot(theta_probe_yield_plot, "day39_theta_low_effort_yield.png", width = 9, height = 6)

cat("Section 1 (theta=0.3 vs. baseline effort probe, knife_edge=10): see day39_theta_low_effort_probe.png / day39_theta_low_effort_yield.png / .csv.\n")
print(theta_probe_df)

################################################################################
# Section 2: same cannibalism cut, but with interaction_resource raised to
# guard against Section 1's growth-rate risk
#
# Section 1 found growth barely moved with theta in this particular
# parameterisation, but the risk was real in principle, and the fix this
# task asks for is a direct one: species_params$interaction_resource
# (theta_1p in getEncounter()'s own formula, see Section 1's header) scales
# exactly the resource half of the encounter integral, independent of
# theta_11's cannibalism half. Raising it gives the fish more food from the
# resource to offset less food from cannibalism, without touching mortality
# at all -- getPredMort() only involves the species interaction matrix,
# never interaction_resource, so this section's mortality profile is
# identical to Section 1's own theta_low case (not re-shown here).
#
# interaction_resource=1.3 (30% more resource food) was an early, hand-picked
# guess -- checked here against growth_resboost_df, read back once the
# script actually ran, not assumed to have worked:
#   - unfished, the settle+kick fork still cycles.
#   - at the small sizes that were ALREADY fine or improved under theta=0.3
#     alone, +1.3x resource pushes growth further above the ORIGINAL
#     baseline: +31% at w=0.01g, +28% at w=0.1g, +15% at w=1g.
#   - at the two LARGE sizes that actually needed rescuing, it barely moves
#     the needle: w=10.1g goes from -88% (theta=0.3 alone) to -86% versus
#     baseline; w=28.6g goes from -89% to -87%. A 30% resource boost against
#     an ~88% growth deficit was never going to close that gap, and it
#     doesn't -- growth at w>=w_mat is STILL crushed to roughly a seventh to
#     an eighth of its original value.
#
# interaction_resource=1.3 is not a working fix for Section 1's growth
# problem, just a first guess that turned out insufficient by a wide
# margin -- it over-corrects the sizes that didn't need it and
# under-corrects the sizes that did. Section 3 sweeps both knobs together
# to find out what, if anything, actually closes that gap.
################################################################################

interaction_resource_boost <- 1.3

p_theta_low_resboost <- make_anchovy_fishing_params_theta(theta_low, interaction_resource_boost,
                                                           knife_edge_size = p2$w_mat)
sim_fork_resboost <- make_anchovy_fork_sim(p_theta_low_resboost, t_fork)
last_n_resboost   <- array(sim_fork_resboost@n[dim(sim_fork_resboost@n)[1], , , drop = FALSE],
                          dim = dim(sim_fork_resboost@n)[-1])
last_npp_resboost <- sim_fork_resboost@n_pp[dim(sim_fork_resboost@n_pp)[1], ]

bp_resboost <- compute_selected_biomass_series(sim_fork_resboost, p_theta_low_resboost,
                                              t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "Section 2 (theta=%.2g, interaction_resource=%.2g), unfished fork, last %.0fyr of %.0fyr: selected biomass min=%.4g max=%.4g -- still cycling.\n",
  theta_low, interaction_resource_boost, scan_summary_window, t_fork, min(bp_resboost), max(bp_resboost)
))

growth_resboost_df <- data.frame(
  w                         = round(p_scan@w[idx_check], 4),
  growth_baseline           = getEGrowth(p_scan, n = last_n, n_pp = last_npp)[1, idx_check],
  growth_theta_low          = getEGrowth(p_theta_low, n = last_n_theta_low, n_pp = last_npp_theta_low)[1, idx_check],
  growth_theta_low_resboost = getEGrowth(p_theta_low_resboost, n = last_n_resboost, n_pp = last_npp_resboost)[1, idx_check]
)
growth_resboost_df$pct_vs_baseline <- 100 * (growth_resboost_df$growth_theta_low_resboost /
                                              growth_resboost_df$growth_baseline - 1)
write.csv(growth_resboost_df, file.path("interesting_plots", "day39_theta_resboost_growth.csv"), row.names = FALSE)
cat("Section 2 growth comparison (baseline vs. theta=0.3 vs. theta=0.3+interaction_resource=1.3) -- see day39_theta_resboost_growth.csv:\n")
print(growth_resboost_df)

resboost_probe_s2 <- cannibalism_fishing_probe(p_theta_low_resboost, last_n_resboost, last_npp_resboost,
                                               fish_level_seq_theta) %>%
  mutate(model = "theta=0.3, interaction_resource=1.3")

theta_probe_df_s2 <- bind_rows(theta_probe_df, resboost_probe_s2)
write.csv(theta_probe_df_s2, file.path("interesting_plots", "day39_theta_resboost_effort_probe.csv"), row.names = FALSE)

resboost_probe_plot <- ggplot(theta_probe_df_s2, aes(x = fish_level, y = mean_total, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=10)", y = "Mean total biomass",
       title = "Does the resource-interaction bump change the fishing-resistance pattern?",
       subtitle = "All three models share the same gear (knife_edge=10) -- only theta/interaction_resource differ") +
  theme_minimal()
resboost_probe_plot
save_plot(resboost_probe_plot, "day39_theta_resboost_effort_probe.png", width = 9, height = 6)

resboost_probe_yield_plot <- ggplot(theta_probe_df_s2, aes(x = fish_level, y = mean_yield, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=10)", y = "Mean yield",
       title = "Yield vs. effort, all three variants -- does the resource bump change where yield peaks?",
       subtitle = "All three models share the same gear (knife_edge=10) -- only theta/interaction_resource differ") +
  theme_minimal()
resboost_probe_yield_plot
save_plot(resboost_probe_yield_plot, "day39_theta_resboost_effort_yield.png", width = 9, height = 6)

cat("Section 2 (theta=0.3+interaction_resource=1.3 effort probe, knife_edge=10): see day39_theta_resboost_effort_probe.png / day39_theta_resboost_effort_yield.png / .csv.\n")
print(theta_probe_df_s2)

################################################################################
# Section 3: an actual sweep of theta x interaction_resource, to find a
# combination that's actually good, rather than hand-picking one pair and
# hoping (which Sections 1/2 did, and which Section 2's own growth_
# resboost_df shows failed -- see Section 2's header)
#
# Two things Sections 1/2 were separately trying to buy pull in OPPOSITE
# directions, which is exactly why hand-picking one point in a 2D space
# didn't work:
#   - MORE fishing sensitivity (Section 1's goal) wants LESS cannibalism
#     (theta down), because that's what weakens the cultivation-effect
#     relief fishing currently gets rewarded with.
#   - RESTORED growth (Section 2's goal) wants MORE resource food
#     (interaction_resource up) to offset the population boom that theta
#     down causes -- but more resource food also means a bigger, better-fed
#     population overall, which Section 4's own knife_edge=10/9 numbers
#     below show makes the model MORE resistant to fishing, not less. The
#     resource fix risks undoing the cannibalism fix.
#
# So this is a genuine two-objective trade-off, not a single number to
# solve for. Sweep both knobs on a grid, read off three things per cell,
# and look at the trade-off directly rather than asserting a winner:
#   - unfished_cycle_ratio: is it even still a working oscillator?
#   - growth_recovery_w10: growth at w=10.1g (the size Section 1 crashed by
#     -88%) as a fraction of the ORIGINAL theta=1/interaction_resource=1
#     baseline's own growth there. 1.0 = fully restored.
#   - fishing_sensitivity: fraction of unfished total biomass lost at a
#     single moderate probe effort (fish_level=10, Day 38 Section 7's own
#     "healthy" checkpoint for knife_edge=5) -- higher = more responsive to
#     fishing, which is what Section 1 was chasing in the first place.
#
# Kept to knife_edge=10 (the default gear) throughout and a single probe
# effort rather than a full effort sweep per cell -- a 5x5 grid at
# Section 4's own per-cell cost (fork + multi-point sweep) would be far too
# expensive; this section is about finding a better (theta,
# interaction_resource) pair, not about characterising it fully once found.
# The (theta=1, interaction_resource=1) cell IS the original baseline,
# included in the grid rather than computed separately, so it doubles as
# the reference point for both growth_recovery (=1 there by construction)
# and fishing_sensitivity (whatever that baseline's own value comes out
# to) -- read as directly comparable numbers, not two different
# computations.
################################################################################

theta_grid_seq     <- c(1, 0.7, 0.5, 0.3, 0.1)
resource_grid_seq  <- c(1, 1.3, 2, 3, 5)
grid_probe_effort  <- 10

w_grid_check         <- c(0.01, 1, 10, 30)
idx_grid_check       <- vapply(w_grid_check, function(x) which.min(abs(p_scan@w - x)), integer(1))
growth_baseline_grid <- getEGrowth(p_scan, n = last_n, n_pp = last_npp)[1, idx_grid_check]

run_interaction_grid_case <- function(theta_val, ir_val) {
  params <- make_anchovy_fishing_params_theta(theta_val, ir_val, knife_edge_size = p2$w_mat)
  sim_fork_grid <- make_anchovy_fork_sim(params, t_fork)
  last_n_grid   <- array(sim_fork_grid@n[dim(sim_fork_grid@n)[1], , , drop = FALSE],
                        dim = dim(sim_fork_grid@n)[-1])
  last_npp_grid <- sim_fork_grid@n_pp[dim(sim_fork_grid@n_pp)[1], ]

  bp_grid        <- compute_selected_biomass_series(sim_fork_grid, params, t_cut = t_fork - scan_summary_window)
  unfished_ratio <- max(bp_grid) / max(min(bp_grid), 1e-12)

  growth_grid      <- getEGrowth(params, n = last_n_grid, n_pp = last_npp_grid)[1, idx_grid_check]
  growth_recovery  <- growth_grid / growth_baseline_grid   # 1 = fully recovered vs. the ORIGINAL baseline

  tv_fork_grid       <- as.numeric(dimnames(sim_fork_grid@n)[[1]])
  total_unfished     <- mean(unname(getBiomass(sim_fork_grid)[, "Anchovy"])[
    tv_fork_grid > t_fork - scan_summary_window])

  sim_fished <- project(params, t_max = scan_summary_window, dt = p2$dt, t_save = 0.2,
                        t_start = t_fork, progress_bar = FALSE, effort = grid_probe_effort)
  tv_fished           <- as.numeric(dimnames(sim_fished@n)[[1]])
  keep_fished         <- tv_fished > t_fork + scan_summary_window / 2
  total_at_effort     <- mean(unname(getBiomass(sim_fished)[, "Anchovy"])[keep_fished])
  yield_at_effort      <- mean(unname(getYield(sim_fished)[, "Anchovy"])[keep_fished])
  fishing_sensitivity  <- 1 - total_at_effort / total_unfished   # fraction of unfished biomass lost

  data.frame(theta = theta_val, interaction_resource = ir_val,
            unfished_cycle_ratio = unfished_ratio,
            growth_recovery_w0.01 = growth_recovery[1], growth_recovery_w1 = growth_recovery[2],
            growth_recovery_w10 = growth_recovery[3], growth_recovery_w30 = growth_recovery[4],
            total_unfished = total_unfished, total_at_effort10 = total_at_effort,
            mean_yield_at_effort10 = yield_at_effort, fishing_sensitivity = fishing_sensitivity)
}

interaction_grid_df <- bind_rows(lapply(theta_grid_seq, function(th) {
  bind_rows(lapply(resource_grid_seq, function(ir) run_interaction_grid_case(th, ir)))
}))
write.csv(interaction_grid_df, file.path("interesting_plots", "day39_interaction_grid.csv"), row.names = FALSE)
cat("Section 3 (theta x interaction_resource grid): see day39_interaction_grid.csv.\n")
print(interaction_grid_df)

grid_growth_heatmap <- ggplot(interaction_grid_df,
                              aes(x = factor(theta), y = factor(interaction_resource), fill = growth_recovery_w10)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", growth_recovery_w10)), size = 3, color = "white") +
  scale_fill_viridis_c(limits = c(0, NA)) +
  labs(x = "theta (cannibalism)", y = "interaction_resource", fill = "recovery",
       title = "Growth recovery at w=10.1g",
       subtitle = "1.0 = fully back to the original theta=1/interaction_resource=1 baseline's own growth there") +
  theme_minimal()
grid_growth_heatmap
save_plot(grid_growth_heatmap, "day39_grid_growth_heatmap.png", width = 8, height = 6)

grid_sensitivity_heatmap <- ggplot(interaction_grid_df,
                                   aes(x = factor(theta), y = factor(interaction_resource), fill = fishing_sensitivity)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", fishing_sensitivity)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "magma", limits = c(0, NA)) +
  labs(x = "theta (cannibalism)", y = "interaction_resource", fill = "sensitivity",
       title = "Fishing sensitivity: fraction of unfished biomass lost at effort=10",
       subtitle = "Higher = more responsive to fishing -- Section 1's original goal") +
  theme_minimal()
grid_sensitivity_heatmap
save_plot(grid_sensitivity_heatmap, "day39_grid_sensitivity_heatmap.png", width = 8, height = 6)

grid_tradeoff_plot <- ggplot(interaction_grid_df,
                             aes(x = fishing_sensitivity, y = growth_recovery_w10,
                                 color = factor(theta), shape = factor(interaction_resource))) +
  geom_point(size = 3) +
  labs(x = "Fishing sensitivity (fraction of unfished biomass lost at effort=10)",
       y = "Growth recovery at w=10.1g (1.0 = original baseline)",
       color = "theta", shape = "interaction_resource",
       title = "The trade-off this section is actually navigating",
       subtitle = "Top-right is the ideal quadrant: growth restored AND still fishing-sensitive") +
  theme_minimal()
grid_tradeoff_plot
save_plot(grid_tradeoff_plot, "day39_grid_tradeoff.png", width = 9, height = 6)

baseline_grid_row <- interaction_grid_df %>% filter(theta == 1, interaction_resource == 1)
cat(sprintf(
  "Section 3 baseline cell (theta=1, interaction_resource=1): fishing_sensitivity=%.3f at effort=10 -- the bar Section 1 was trying to clear.\n",
  baseline_grid_row$fishing_sensitivity
))

grid_candidates <- interaction_grid_df %>%
  filter(growth_recovery_w10 > 0.5, fishing_sensitivity > baseline_grid_row$fishing_sensitivity) %>%
  arrange(desc(growth_recovery_w10))
write.csv(grid_candidates, file.path("interesting_plots", "day39_grid_candidates.csv"), row.names = FALSE)
cat(sprintf(
  "Candidates beating the baseline's own fishing_sensitivity (%.3f) while recovering over half of its growth at w=10.1g: %d found -- see day39_grid_candidates.csv.\n",
  baseline_grid_row$fishing_sensitivity, nrow(grid_candidates)
))
print(grid_candidates)

################################################################################
# Section 4: is the reduced-cannibalism model actually dependent on fishing
# effort? Sweep at the default knife-edge size (10, = w_mat) and two lower
# ones (9, 5)
#
# FIRST PASS AT THIS SECTION WAS FLAWED, fixed here rather than patched: it
# compared "Constant" (fishing the nominal fish_level all the time) against
# nothing at all otherwise (a schedule concept that was never built into
# this section in the first place -- it was Constant-only). That's not a
# comparison of WHEN to fish, only of HOW MUCH, and it silently dropped the
# one comparison this project actually cares about -- whether *timing*
# fishing around the cycle matters. Day 38's own What's Next (item 1)
# flagged the general fix: "keep Constant fishing running as a baseline all
# the time, and layer extra effort on top specifically during peaks/
# troughs" -- so the question becomes whether topping up an existing
# fishery in good years beats fishing that same floor flat, not whether an
# intermittent fishery beats a continuous one from a standing start.
#
# Built here with thresholdFMort()/attach_threshold_rule() (Day 36/38's own
# rate-function machinery, ported in below, self-contained per this file's
# own convention) using background_level = baseline_effort_s4 for EVERY
# schedule, including "Constant" (which is just the degenerate case of a
# threshold rule that never switches). "Threshold (peaks)"/"Threshold
# (troughs)" add a boost on top of that same floor only near cycle peaks/
# troughs. All three schedules now share the same floor, so any difference
# between them is genuinely about timing, not about average effort.
#
# Sections 1/2's own probes kept the gear fixed at knife_edge=10 to isolate
# theta's effect -- and at that gear, Day 38 Section 3c/7's structural
# juvenile refuge (the gear never touches w<w_mat, so >95% of biomass is
# untouchable regardless of theta) dominates over anything theta does, so
# both variants stayed resistant. Day 38 Section 7 found that refuge is
# exactly what a lower knife_edge_size cuts into. This section reruns that
# same lever -- knife_edge=10 (default), 9 (a small cut just below w_mat),
# 5 (Day 38 Section 7's own collapse-capable value) -- on Section 2's model
# (theta=0.3, interaction_resource=1.3), NOT on whatever Section 3's grid
# sweep turns up as a better candidate. Deliberate, not an oversight: this
# section's own machinery (per-knife_edge fork, threshold calibration
# against that fork) is expensive enough that re-running it against a
# second (theta, interaction_resource) pair roughly doubles this section's
# cost, and Section 3's grid was evaluated at a single knife_edge=10 probe
# point, not validated across knife_edge=9/5 -- rerunning THIS section on a
# grid candidate before checking that candidate still cycles and still
# collapses sensibly at knife_edge=5 would be building on an unverified
# choice. Natural next step, not done today.
#
# baseline_effort_s4=10 is Day 38 Section 7's own "healthy" checkpoint for
# knife_edge=5 (mean_total=0.381, ~unfished). boost_fish_level_seq_s4 pushes
# the ON-period effort from there up into Section 7's own stressed/extinct
# range (30 -> stressed, 100 -> extinct at theta=1). See the blog post for
# why these numbers -- effort in the tens, sometimes 100 -- are not
# remotely realistic fishing efforts; that's flagged there as unresolved,
# not fixed here.
#
# Reference point, not rerun here: Day 38 Section 7's OWN numbers for the
# original (theta=1, interaction_resource=1) model, Constant schedule only,
# at these same three knife-edge sizes are already in
# day38_knife_edge_summary.csv -- knife_edge=10 stayed resistant to
# effort=30 (mean_total still rising), knife_edge=5 was healthy at
# effort=10 (0.381), stressed at 30 (0.083), extinct by 100 (3.5e-10).
# knife_edge=9 wasn't in Day 38 Section 7's own sweep (c(10,5,1,0.1,w_l)) --
# new here, and Day 38 never tested layered schedules at all (only
# Constant) -- so this section has no direct Day 38 precedent to reuse for
# the threshold schedules, only for the Constant floor.
################################################################################

# Day 36/38's own threshold-rule machinery, ported in verbatim (self-
# contained convention).
thresholdFMort <- function(params, n, n_pp, n_other, t, effort, e_growth, pred_mort, ...) {
  p <- other_params(params)

  f_ref            <- mizerFMortGear(params, effort = 1)
  biomass_density  <- sweep(n, 2, params@w * params@dw, "*")
  selected_biomass <- sum(colSums(f_ref) * biomass_density)

  direction <- if (p$mode == "above") 1 else -1
  on_frac   <- if (isTRUE(p$hard_step)) {
    as.numeric(direction * (selected_biomass - p$threshold) > 0)
  } else {
    plogis(direction * (selected_biomass - p$threshold) / max(p$sharpness, 1e-12))
  }

  fmort_on  <- colSums(mizerFMortGear(params, effort = p$fish_level))
  fmort_off <- colSums(mizerFMortGear(params, effort = p$background_level))
  result <- on_frac * fmort_on + (1 - on_frac) * fmort_off
  dim(result)      <- dim(n)
  dimnames(result) <- dimnames(n)
  result
}

attach_threshold_rule <- function(params, threshold, fish_level, background_level = 0,
                                  mode = c("above", "below"), sharpness, hard_step = FALSE) {
  mode <- match.arg(mode)
  other_params(params) <- list(threshold = threshold, fish_level = fish_level,
                               background_level = background_level, mode = mode,
                               sharpness = sharpness, hard_step = hard_step)
  setRateFunction(params, "FMort", "thresholdFMort")
}

compute_on_frac_series <- function(sim, params, threshold, mode, sharpness, t_cut, hard_step = FALSE) {
  tv    <- as.numeric(dimnames(sim@n)[[1]])
  keep  <- which(tv > t_cut)
  f_ref <- mizerFMortGear(params, effort = 1)
  direction <- if (mode == "above") 1 else -1

  vapply(keep, function(i) {
    n_i <- array(sim@n[i, , , drop = FALSE], dim = dim(sim@n)[-1])
    biomass_density  <- sweep(n_i, 2, params@w * params@dw, "*")
    selected_biomass <- sum(colSums(f_ref) * biomass_density)
    if (isTRUE(hard_step)) {
      as.numeric(direction * (selected_biomass - threshold) > 0)
    } else {
      plogis(direction * (selected_biomass - threshold) / max(sharpness, 1e-12))
    }
  }, numeric(1))
}

threshold_diagnostics <- function(sim, threshold, fish_level, background_level, mode, sharpness, t_cut, hard_step = FALSE) {
  tv      <- as.numeric(dimnames(sim@n)[[1]])
  keep    <- tv > t_cut
  on_frac <- compute_on_frac_series(sim, sim@params, threshold, mode, sharpness, t_cut, hard_step)
  is_on   <- on_frac > 0.5

  runs   <- rle(is_on)
  n_runs <- length(runs$values)

  # First/last run in the window are censored -- see Day 36's own note.
  interior <- rep(TRUE, n_runs)
  if (n_runs > 0) { interior[1] <- FALSE; interior[n_runs] <- FALSE }
  on_lengths <- runs$lengths[runs$values & interior]
  dt_save    <- mean(diff(tv[keep]))

  list(
    mean_effort      = mean(on_frac * fish_level + (1 - on_frac) * background_level),
    effective_window = if (length(on_lengths) > 0) mean(on_lengths) * dt_save else NA_real_,
    n_bursts         = length(on_lengths)
  )
}

baseline_effort_s4      <- 10
# Extended from c(20,30,50,75,100) -- the original range never showed
# knife_edge=10/9 turning over (they plateau rather than collapse, see
# below), so it was never clear whether that's a genuine plateau or just
# not enough range. 150/200 settle that.
boost_fish_level_seq_s4 <- c(20, 30, 50, 75, 100, 150, 200)
knife_edge_seq_s4       <- c(10, 9, 5)

scan_metrics_layered <- function(sim) {
  tv   <- as.numeric(dimnames(sim@n)[[1]])
  keep <- tv > scan_t_cut
  total_bm    <- unname(getBiomass(sim)[, "Anchovy"])[keep]
  selected_bm <- unname(getBiomass(sim, min_w = p2$w_mat)[, "Anchovy"])[keep]
  yield_bm    <- unname(getYield(sim)[, "Anchovy"])[keep]
  data.frame(min_total = min(total_bm), mean_total = mean(total_bm),
            mean_selected = mean(selected_bm), mean_juvenile = mean(total_bm - selected_bm),
            mean_yield = mean(yield_bm))
}

run_layered_ke_case <- function(knife_edge_size) {
  params <- make_anchovy_fishing_params_theta(theta_low, interaction_resource_boost,
                                              knife_edge_size = knife_edge_size)
  params@initial_n[]    <- last_n_resboost   # Section 2's own fork state -- same
  params@initial_n_pp[] <- last_npp_resboost # starting point for every knife-edge value

  # "Constant": the shared floor every schedule below is layered on top of.
  # Run once per knife_edge_size and reused twice -- as the Constant result
  # itself, and as the reference series the threshold/sharpness below are
  # calibrated against (the floor's OWN cycle, not a hypothetical unfished
  # one, since the "off" state for every schedule here is "fishing at the
  # floor", not "not fishing at all").
  sim_const <- project(params, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                       t_start = t_fork, progress_bar = FALSE, effort = baseline_effort_s4)
  bp_const     <- compute_selected_biomass_series(sim_const, params, scan_t_cut)
  sharpness_ke <- 0.02 * (max(bp_const) - min(bp_const))
  threshold_ke <- unname(quantile(bp_const, probs = 0.5))

  const_row <- cbind(scan_metrics_layered(sim_const),
                     knife_edge_size = knife_edge_size, schedule = "Constant",
                     boost_level = NA_real_, mean_effort = baseline_effort_s4,
                     effective_window = NA_real_, n_bursts = NA_integer_)

  boosted_rows <- bind_rows(lapply(boost_fish_level_seq_s4, function(boost) {
    bind_rows(lapply(c("above", "below"), function(mode) {
      schedule_name <- if (mode == "above") "Threshold (peaks)" else "Threshold (troughs)"
      p_rule <- attach_threshold_rule(params, threshold = threshold_ke, fish_level = boost,
                                      background_level = baseline_effort_s4, mode = mode,
                                      sharpness = sharpness_ke)
      p_rule@initial_n[]    <- last_n_resboost
      p_rule@initial_n_pp[] <- last_npp_resboost
      sim <- project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                     t_start = t_fork, progress_bar = FALSE, effort = 1)
      diag <- threshold_diagnostics(sim, threshold_ke, boost, baseline_effort_s4, mode,
                                    sharpness_ke, scan_t_cut)
      cbind(scan_metrics_layered(sim),
           knife_edge_size = knife_edge_size, schedule = schedule_name,
           boost_level = boost, mean_effort = diag$mean_effort,
           effective_window = diag$effective_window, n_bursts = diag$n_bursts)
    }))
  }))

  bind_rows(const_row, boosted_rows)
}

theta_ke_layered_summary <- bind_rows(lapply(knife_edge_seq_s4, run_layered_ke_case))
write.csv(theta_ke_layered_summary, file.path("interesting_plots", "day39_theta_ke_layered_summary.csv"), row.names = FALSE)

theta_ke_layered_plot <- ggplot(theta_ke_layered_summary,
                                aes(x = mean_effort, y = pmax(mean_total, 1e-12), color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = baseline_effort_s4, linetype = "dashed", color = "grey50") +
  facet_wrap(~knife_edge_size, labeller = label_both) +
  scale_y_log10() +
  labs(x = "Realised mean effort actually applied", y = "Mean total biomass (log scale)",
       color = "Schedule",
       title = "Constant floor vs. boosting at peaks vs. boosting at troughs -- same floor throughout",
       subtitle = sprintf("baseline_effort=%.2g (dashed line) -- every schedule fishes at least this much, always. Reduced-cannibalism model (theta=0.3, interaction_resource=1.3).",
                          baseline_effort_s4)) +
  theme_minimal()
theta_ke_layered_plot
save_plot(theta_ke_layered_plot, "day39_theta_ke_layered_collapse.png", width = 11, height = 5)

theta_ke_layered_yield_plot <- ggplot(theta_ke_layered_summary,
                                      aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = baseline_effort_s4, linetype = "dashed", color = "grey50") +
  facet_wrap(~knife_edge_size, labeller = label_both, scales = "free_y") +
  labs(x = "Realised mean effort actually applied", y = "Mean yield",
       color = "Schedule",
       title = "Yield vs. realised effort -- same floor, boost only at peaks/troughs",
       subtitle = sprintf("baseline_effort=%.2g (dashed line). Reduced-cannibalism model (theta=0.3, interaction_resource=1.3).",
                          baseline_effort_s4)) +
  theme_minimal()
theta_ke_layered_yield_plot
save_plot(theta_ke_layered_yield_plot, "day39_theta_ke_layered_yield.png", width = 11, height = 5)

# Nominal boost vs. what actually got applied -- Day 38's own headline
# finding was that "Threshold (peaks)" self-limits (a bigger nominal burst
# depletes the peak faster, so effective_window shrinks and mean_effort
# saturates well below the nominal fish_level). Worth checking that holds
# here too, now that there's a nonzero floor under it.
boost_realised_plot <- ggplot(theta_ke_layered_summary %>% filter(schedule != "Constant"),
                              aes(x = boost_level, y = mean_effort, color = schedule)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~knife_edge_size, labeller = label_both) +
  labs(x = "Nominal boost (on-period fish_level)", y = "Realised mean effort",
       color = "Schedule",
       title = "How much of the nominal boost actually gets applied, on average?",
       subtitle = "Dotted line = boost applied in full, all the time (upper bound); floor=10 is the lower bound") +
  theme_minimal()
boost_realised_plot
save_plot(boost_realised_plot, "day39_theta_ke_boost_realised.png", width = 11, height = 5)

cat(sprintf(
  "Section 4 (layered schedules -- Constant floor=%.2g, boosted at peaks/troughs, knife_edge=10/9/5): see day39_theta_ke_layered_collapse.png / day39_theta_ke_layered_yield.png / day39_theta_ke_boost_realised.png / day39_theta_ke_layered_summary.csv.\n",
  baseline_effort_s4
))
print(theta_ke_layered_summary)

theta_ke_layered_trend <- theta_ke_layered_summary %>%
  group_by(knife_edge_size, schedule) %>%
  summarise(total_trend = cor(mean_effort, mean_total), .groups = "drop") %>%
  arrange(desc(knife_edge_size), schedule)
cat("Correlation of realised mean_effort with mean_total, by knife_edge_size and schedule (negative = fishing-dependent/collapsing, positive/NA = still resistant or too few points):\n")
print(theta_ke_layered_trend)

################################################################################
# Section 5: summary
################################################################################

cat("\n===== Day 39 summary =====\n")
cat("Day 38's Figure 2e plankton-anchovy cannibalism model, with its dependence on cannibalism (the species interaction matrix, theta_11) turned down and its dependence on the resource (interaction_resource) turned up.\n")
cat(sprintf(
  "Section 1 (theta=%.2g cannibalism cut, gear unchanged at knife_edge=10): day39_theta_low_growth_mort.csv / day39_theta_low_effort_probe.png / day39_theta_low_effort_yield.png.\n",
  theta_low
))
cat(sprintf(
  "Section 2 (theta=%.2g + interaction_resource=%.2g -- shown NOT to fix Section 1's growth deficit at w>=w_mat, see Section 2's own header): day39_theta_resboost_growth.csv / day39_theta_resboost_effort_probe.png / day39_theta_resboost_effort_yield.png.\n",
  theta_low, interaction_resource_boost
))
cat(sprintf(
  "Section 3 (theta x interaction_resource grid, %dx%d cells): day39_interaction_grid.csv / day39_grid_growth_heatmap.png / day39_grid_sensitivity_heatmap.png / day39_grid_tradeoff.png -- %d candidate(s) beat the original baseline's own fishing_sensitivity while recovering over half its growth at w=10.1g, see day39_grid_candidates.csv.\n",
  length(theta_grid_seq), length(resource_grid_seq), nrow(grid_candidates)
))
cat(sprintf(
  "Section 4 (layered schedules -- Constant floor=%.2g with boosts at peaks/troughs, knife_edge=10/9/5 on Section 2's theta=0.3/interaction_resource=1.3 model, NOT the grid's own best candidate -- see Section 4's header for why): day39_theta_ke_layered_collapse.png / day39_theta_ke_layered_yield.png / day39_theta_ke_boost_realised.png / day39_theta_ke_layered_summary.csv.\n",
  baseline_effort_s4
))
cat(sprintf(
  "NOTE: the effort values Section 4 needed to see any collapse at all (floor=%.2g, boosts up to %.2g) are nowhere near realistic fishing efforts (F is typically O(0.1-1)/yr for real fisheries) -- this is inherited from Day 38 Section 7's own finding that knife_edge=5 doesn't visibly stress until effort~30 and doesn't collapse until ~100. Still unresolved, not fixed today -- see the blog post.\n",
  baseline_effort_s4, max(boost_fish_level_seq_s4)
))
