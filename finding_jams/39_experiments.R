library(mizer)
library(mizerExperimental)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)

# Day 39: cut Day 38's cannibalism model's dependence on cannibalism (S1),
# compensate via resource interaction (S2), grid-sweep both (S3), test fishing-schedule dependence at 3 knife-edges on the S2 model and the original model (S4-5),
# then zoom in on theta=0.3 alone at knife_edge=5 over a realistic effort range (S6),
# and scan the most fishing-sensitive model (theta=0.3+ir=1.3) at the one knife-edge it hadn't covered yet (S7).
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
# Section 0: rebuild Day 38's Figure 2e model, verbatim from 38_experiments.R.
# Kernel/mortality already checked against the paper in Day 38 -- not repeated here.
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
# Section 0b: settle+kick fork (Day 37/38's recipe) to t_fork=20 -- cycle is
# already full-amplitude well before then (confirmed below).
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
# Section 1: decrease dependence on cannibalism, to make the model more
# fishing-sensitive.
#
# theta_11 (@interaction) is the cannibalism knob; interaction_resource (Section 2) is separate.
# Cutting theta to 0.3 crashes growth -88%/-89% at w=10.1/28.6g via indirect resource competition
# (juvenile population boom grazes the shared resource), not direct diet loss -- see growth_compare_df.
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

# Growth vs. predation-mortality, baseline (theta=1) vs. theta_low, both on their OWN fork state.
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

# Extended from [1,3,5,7,9] up to 100, matching Section 4's scale -- [1,9] alone
# showed almost no curvature. Gear fixed at knife_edge=10 to isolate theta's effect.
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
              mean_yield = mean(yield),
              rel_amplitude = (max(total) - min(total)) / ((max(total) + min(total)) / 2))
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

# Relative amplitude alongside yield -- is the population still genuinely cycling
# at each effort level, or has yield kept climbing while the cycle itself died out?
theta_probe_amplitude_plot <- ggplot(theta_probe_df, aes(x = fish_level, y = rel_amplitude, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=10)", y = "Relative amplitude of biomass",
       title = "Relative amplitude vs. effort -- alongside the yield plot above",
       subtitle = "Gear unchanged (knife_edge=10)") +
  theme_minimal()
theta_probe_amplitude_plot
save_plot(theta_probe_amplitude_plot, "day39_theta_low_effort_amplitude.png", width = 9, height = 6)

cat("Section 1 (theta=0.3 vs. baseline effort probe, knife_edge=10): see day39_theta_low_effort_probe.png / day39_theta_low_effort_yield.png / day39_theta_low_effort_amplitude.png / .csv.\n")
print(theta_probe_df)

################################################################################
# Section 2: same cannibalism cut, with interaction_resource raised to guard
# against Section 1's growth-rate risk.
#
# interaction_resource=1.3 does NOT fix it: growth at w=10.1/28.6g only recovers to -86%/-87%
# vs. baseline (was -88%/-89%), while over-correcting the small sizes that didn't need it.
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

resboost_probe_amplitude_plot <- ggplot(theta_probe_df_s2, aes(x = fish_level, y = rel_amplitude, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=10)", y = "Relative amplitude of biomass",
       title = "Relative amplitude vs. effort, all three variants -- alongside the yield plot above",
       subtitle = "All three models share the same gear (knife_edge=10) -- only theta/interaction_resource differ") +
  theme_minimal()
resboost_probe_amplitude_plot
save_plot(resboost_probe_amplitude_plot, "day39_theta_resboost_effort_amplitude.png", width = 9, height = 6)

cat("Section 2 (theta=0.3+interaction_resource=1.3 effort probe, knife_edge=10): see day39_theta_resboost_effort_probe.png / day39_theta_resboost_effort_yield.png / day39_theta_resboost_effort_amplitude.png / .csv.\n")
print(theta_probe_df_s2)

# Same three-way comparison, knife_edge=5 -- Day 38's own collapse-capable gear,
# vs. the knife_edge=10 isolation above. Unfished trajectory is gear-independent
# (Day 38 Section 7), so last_n/last_n_theta_low/last_n_resboost are reused.
p_scan_ke5              <- make_anchovy_fishing_params_theta(1, 1, knife_edge_size = 5)
p_theta_low_ke5         <- make_anchovy_fishing_params_theta(theta_low, 1, knife_edge_size = 5)
p_theta_low_resboost_ke5 <- make_anchovy_fishing_params_theta(theta_low, interaction_resource_boost, knife_edge_size = 5)

theta_probe_df_ke5 <- bind_rows(
  cannibalism_fishing_probe(p_scan_ke5, last_n, last_npp, fish_level_seq_theta) %>%
    mutate(model = "Baseline (theta=1)"),
  cannibalism_fishing_probe(p_theta_low_ke5, last_n_theta_low, last_npp_theta_low, fish_level_seq_theta) %>%
    mutate(model = "theta=0.3"),
  cannibalism_fishing_probe(p_theta_low_resboost_ke5, last_n_resboost, last_npp_resboost, fish_level_seq_theta) %>%
    mutate(model = "theta=0.3, interaction_resource=1.3")
)
write.csv(theta_probe_df_ke5, file.path("interesting_plots", "day39_theta_resboost_effort_yield_ke5.csv"), row.names = FALSE)

theta_probe_yield_plot_ke5 <- ggplot(theta_probe_df_ke5, aes(x = fish_level, y = mean_yield, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=5)", y = "Mean yield",
       title = "Yield vs. effort at knife_edge=5, all three variants -- does gear that actually reaches the collapse zone change the picture?",
       subtitle = "Day 38's own collapse-capable gear setting, vs. Sections 1/2's knife_edge=10 isolation above") +
  theme_minimal()
theta_probe_yield_plot_ke5
save_plot(theta_probe_yield_plot_ke5, "day39_theta_resboost_effort_yield_ke5.png", width = 9, height = 6)

theta_probe_amplitude_plot_ke5 <- ggplot(theta_probe_df_ke5, aes(x = fish_level, y = rel_amplitude, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=5)", y = "Relative amplitude of biomass",
       title = "Relative amplitude vs. effort at knife_edge=5, all three variants -- alongside the yield plot above",
       subtitle = "Day 38's own collapse-capable gear setting") +
  theme_minimal()
theta_probe_amplitude_plot_ke5
save_plot(theta_probe_amplitude_plot_ke5, "day39_theta_resboost_effort_amplitude_ke5.png", width = 9, height = 6)

cat("Section 2 (yield vs. effort at knife_edge=5, all three variants): see day39_theta_resboost_effort_yield_ke5.png / day39_theta_resboost_effort_amplitude_ke5.png / .csv.\n")
print(theta_probe_df_ke5)

################################################################################
# Section 3: actual grid sweep of theta x interaction_resource -- Sections
# 1/2 pull in opposite directions (less theta = more fishing-sensitive but
# less growth; more interaction_resource = more growth but less fishing-sensitive),
# so this looks for a combination that beats the baseline on both, not just one hand-picked pair.
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
# Section 4: is the reduced-cannibalism model dependent on fishing effort?
# Sweep at knife_edge=10/9/5.
#
# First pass compared Constant fishing against nothing at all -- fixed here: every
# schedule (Constant/Threshold peaks/Threshold troughs) shares the same Constant
# floor (background_level=baseline_effort_s4), so only the boost timing differs.
# Runs on Section 2's model (theta=0.3, ir=1.3); Section 5 reruns it on the original.
################################################################################

# Day 36/38's own threshold-rule machinery, ported in verbatim.
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

# Sanity check before trusting the scan below, Day 36/38's own convention: one
# representative case, all three schedules, biomass plotted against the
# calibration threshold with "on" periods shaded -- does peaks/troughs actually
# fire where it should, or does it barely engage?
sanity_ke    <- 5
sanity_boost <- 50

p_sanity <- make_anchovy_fishing_params_theta(theta_low, interaction_resource_boost, knife_edge_size = sanity_ke)
p_sanity@initial_n[]    <- last_n_resboost
p_sanity@initial_n_pp[] <- last_npp_resboost
sim_const_sanity <- project(p_sanity, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                            t_start = t_fork, progress_bar = FALSE, effort = baseline_effort_s4)
bp_const_sanity  <- compute_selected_biomass_series(sim_const_sanity, p_sanity, scan_t_cut)
sharpness_sanity <- 0.02 * (max(bp_const_sanity) - min(bp_const_sanity))
threshold_sanity <- unname(quantile(bp_const_sanity, probs = 0.5))

run_sanity_schedule <- function(schedule_name) {
  if (schedule_name == "Constant") {
    return(list(sim = sim_const_sanity, mode = NA_character_))
  }
  mode <- if (schedule_name == "Threshold (peaks)") "above" else "below"
  p_rule <- attach_threshold_rule(p_sanity, threshold = threshold_sanity, fish_level = sanity_boost,
                                  background_level = baseline_effort_s4, mode = mode, sharpness = sharpness_sanity)
  p_rule@initial_n[]    <- last_n_resboost
  p_rule@initial_n_pp[] <- last_npp_resboost
  sim <- project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 t_start = t_fork, progress_bar = FALSE, effort = 1)
  list(sim = sim, mode = mode)
}

sanity_schedules <- c("Constant", "Threshold (peaks)", "Threshold (troughs)")
sanity_cases     <- setNames(lapply(sanity_schedules, run_sanity_schedule), sanity_schedules)

sanity_series_df <- bind_rows(lapply(names(sanity_cases), function(nm) {
  case <- sanity_cases[[nm]]
  tv   <- as.numeric(dimnames(case$sim@n)[[1]])
  keep <- which(tv > t_fork)
  bp   <- compute_selected_biomass_series(case$sim, p_sanity, t_fork)
  on_frac <- if (is.na(case$mode)) {
    rep(NA_real_, length(bp))
  } else {
    compute_on_frac_series(case$sim, p_sanity, threshold_sanity, case$mode, sharpness_sanity, t_fork)
  }
  data.frame(t = tv[keep], selected_biomass = bp, on_frac = on_frac, schedule = nm)
}))
write.csv(sanity_series_df, file.path("interesting_plots", "day39_schedule_sanity_series.csv"), row.names = FALSE)

sanity_check_plot <- ggplot(sanity_series_df, aes(x = t, y = selected_biomass)) +
  geom_rect(data = sanity_series_df %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = threshold_sanity, linetype = "dashed", color = "grey40") +
  facet_wrap(~schedule, ncol = 1, scales = "free_y") +
  labs(x = "Time (years)", y = "Selected biomass",
       title = "Sanity check: does thresholdFMort() actually fire where it should?",
       subtitle = sprintf("knife_edge=%.0f, floor=%.2g, boost=%.2g -- dashed = calibration threshold, red = boosted 'on'",
                          sanity_ke, baseline_effort_s4, sanity_boost)) +
  theme_minimal()
sanity_check_plot
save_plot(sanity_check_plot, "day39_schedule_sanity_check.png", width = 9, height = 8)

cat(sprintf(
  "Schedule sanity check (knife_edge=%.0f, boost=%.2g): see day39_schedule_sanity_check.png / day39_schedule_sanity_series.csv.\n",
  sanity_ke, sanity_boost
))

baseline_effort_s4      <- 10
# Extended to 150/200 -- knife_edge=10/9 plateaued rather than collapsed at the old max of 100.
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
            mean_yield = mean(yield_bm),
            rel_amplitude = (max(total_bm) - min(total_bm)) / ((max(total_bm) + min(total_bm)) / 2))
}

# Generalised over (theta_val, ir_val, seed_n, seed_npp) so Section 5 can rerun
# this verbatim on the original model -- only the model differs, not the schedule code.
# baseline_effort/boost_seq default to Sections 4/5's own floor+boost range,
# but are now parameters (not globals) so Section 6 can reuse this with a
# different, much lower range instead of duplicating the whole function.
run_layered_ke_case <- function(knife_edge_size, theta_val, ir_val, seed_n, seed_npp,
                                baseline_effort = baseline_effort_s4,
                                boost_seq = boost_fish_level_seq_s4) {
  params <- make_anchovy_fishing_params_theta(theta_val, ir_val,
                                              knife_edge_size = knife_edge_size)
  params@initial_n[]    <- seed_n
  params@initial_n_pp[] <- seed_npp

  # Constant = the shared floor every other schedule is layered on top of; its own
  # cycle also calibrates the threshold/sharpness used below.
  sim_const <- project(params, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                       t_start = t_fork, progress_bar = FALSE, effort = baseline_effort)
  bp_const     <- compute_selected_biomass_series(sim_const, params, scan_t_cut)
  sharpness_ke <- 0.02 * (max(bp_const) - min(bp_const))
  threshold_ke <- unname(quantile(bp_const, probs = 0.5))

  const_row <- cbind(scan_metrics_layered(sim_const),
                     knife_edge_size = knife_edge_size, schedule = "Constant",
                     boost_level = NA_real_, mean_effort = baseline_effort,
                     effective_window = NA_real_, n_bursts = NA_integer_)

  boosted_rows <- bind_rows(lapply(boost_seq, function(boost) {
    bind_rows(lapply(c("above", "below"), function(mode) {
      schedule_name <- if (mode == "above") "Threshold (peaks)" else "Threshold (troughs)"
      p_rule <- attach_threshold_rule(params, threshold = threshold_ke, fish_level = boost,
                                      background_level = baseline_effort, mode = mode,
                                      sharpness = sharpness_ke)
      p_rule@initial_n[]    <- seed_n
      p_rule@initial_n_pp[] <- seed_npp
      sim <- project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                     t_start = t_fork, progress_bar = FALSE, effort = 1)
      diag <- threshold_diagnostics(sim, threshold_ke, boost, baseline_effort, mode,
                                    sharpness_ke, scan_t_cut)
      cbind(scan_metrics_layered(sim),
           knife_edge_size = knife_edge_size, schedule = schedule_name,
           boost_level = boost, mean_effort = diag$mean_effort,
           effective_window = diag$effective_window, n_bursts = diag$n_bursts)
    }))
  }))

  bind_rows(const_row, boosted_rows)
}

theta_ke_layered_summary <- bind_rows(lapply(knife_edge_seq_s4, function(ke) {
  run_layered_ke_case(ke, theta_low, interaction_resource_boost, last_n_resboost, last_npp_resboost)
}))
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

theta_ke_layered_amplitude_plot <- ggplot(theta_ke_layered_summary,
                                          aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = baseline_effort_s4, linetype = "dashed", color = "grey50") +
  facet_wrap(~knife_edge_size, labeller = label_both) +
  labs(x = "Realised mean effort actually applied", y = "Relative amplitude of biomass",
       color = "Schedule",
       title = "Relative amplitude vs. realised effort -- alongside the yield plot above",
       subtitle = sprintf("baseline_effort=%.2g (dashed line). Reduced-cannibalism model (theta=0.3, interaction_resource=1.3).",
                          baseline_effort_s4)) +
  theme_minimal()
theta_ke_layered_amplitude_plot
save_plot(theta_ke_layered_amplitude_plot, "day39_theta_ke_layered_amplitude.png", width = 11, height = 5)

# Does "Threshold (peaks)" still self-limit (Day 38's finding) now there's a nonzero floor?
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
  "Section 4 (layered schedules -- Constant floor=%.2g, boosted at peaks/troughs, knife_edge=10/9/5): see day39_theta_ke_layered_collapse.png / day39_theta_ke_layered_yield.png / day39_theta_ke_layered_amplitude.png / day39_theta_ke_boost_realised.png / day39_theta_ke_layered_summary.csv.\n",
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
# Section 5: same layered schedules, on the ORIGINAL paper model (theta=1,
# interaction_resource=1) -- Section 4 only tested Section 2's modified model,
# and Sections 1-3 found theta/interaction_resource change fishing response a lot,
# so the original model's own schedule comparison can't just be inferred.
################################################################################

original_ke_layered_summary <- bind_rows(lapply(knife_edge_seq_s4, function(ke) {
  run_layered_ke_case(ke, 1, 1, last_n, last_npp)
}))
write.csv(original_ke_layered_summary, file.path("interesting_plots", "day39_original_ke_layered_summary.csv"), row.names = FALSE)

original_ke_layered_plot <- ggplot(original_ke_layered_summary,
                                   aes(x = mean_effort, y = pmax(mean_total, 1e-12), color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = baseline_effort_s4, linetype = "dashed", color = "grey50") +
  facet_wrap(~knife_edge_size, labeller = label_both) +
  scale_y_log10() +
  labs(x = "Realised mean effort actually applied", y = "Mean total biomass (log scale)",
       color = "Schedule",
       title = "Original paper model (theta=1, interaction_resource=1): Constant floor vs. boosting at peaks vs. troughs",
       subtitle = sprintf("baseline_effort=%.2g (dashed line) -- every schedule fishes at least this much, always.",
                          baseline_effort_s4)) +
  theme_minimal()
original_ke_layered_plot
save_plot(original_ke_layered_plot, "day39_original_ke_layered_collapse.png", width = 11, height = 5)

original_ke_layered_yield_plot <- ggplot(original_ke_layered_summary,
                                         aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = baseline_effort_s4, linetype = "dashed", color = "grey50") +
  facet_wrap(~knife_edge_size, labeller = label_both, scales = "free_y") +
  labs(x = "Realised mean effort actually applied", y = "Mean yield",
       color = "Schedule",
       title = "Original paper model: yield vs. realised effort, same floor throughout",
       subtitle = sprintf("baseline_effort=%.2g (dashed line).", baseline_effort_s4)) +
  theme_minimal()
original_ke_layered_yield_plot
save_plot(original_ke_layered_yield_plot, "day39_original_ke_layered_yield.png", width = 11, height = 5)

original_ke_layered_amplitude_plot <- ggplot(original_ke_layered_summary,
                                             aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = baseline_effort_s4, linetype = "dashed", color = "grey50") +
  facet_wrap(~knife_edge_size, labeller = label_both) +
  labs(x = "Realised mean effort actually applied", y = "Relative amplitude of biomass",
       color = "Schedule",
       title = "Original paper model: relative amplitude vs. realised effort -- alongside the yield plot above",
       subtitle = sprintf("baseline_effort=%.2g (dashed line).", baseline_effort_s4)) +
  theme_minimal()
original_ke_layered_amplitude_plot
save_plot(original_ke_layered_amplitude_plot, "day39_original_ke_layered_amplitude.png", width = 11, height = 5)

original_boost_realised_plot <- ggplot(original_ke_layered_summary %>% filter(schedule != "Constant"),
                                       aes(x = boost_level, y = mean_effort, color = schedule)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~knife_edge_size, labeller = label_both) +
  labs(x = "Nominal boost (on-period fish_level)", y = "Realised mean effort",
       color = "Schedule",
       title = "Original paper model: how much of the nominal boost actually gets applied?",
       subtitle = "Dotted line = boost applied in full, all the time (upper bound); floor=10 is the lower bound") +
  theme_minimal()
original_boost_realised_plot
save_plot(original_boost_realised_plot, "day39_original_ke_boost_realised.png", width = 11, height = 5)

cat(sprintf(
  "Section 5 (layered schedules on the ORIGINAL model, theta=1/interaction_resource=1): see day39_original_ke_layered_collapse.png / day39_original_ke_layered_yield.png / day39_original_ke_layered_amplitude.png / day39_original_ke_boost_realised.png / day39_original_ke_layered_summary.csv.\n"
))
print(original_ke_layered_summary)

original_ke_layered_trend <- original_ke_layered_summary %>%
  group_by(knife_edge_size, schedule) %>%
  summarise(total_trend = cor(mean_effort, mean_total), .groups = "drop") %>%
  arrange(desc(knife_edge_size), schedule)
cat("Correlation of realised mean_effort with mean_total, ORIGINAL model, by knife_edge_size and schedule:\n")
print(original_ke_layered_trend)

# Direct side-by-side: does cutting cannibalism change which schedule wins, or just shift the picture?
model_comparison_df <- bind_rows(
  original_ke_layered_summary %>% mutate(model = "Original (theta=1, ir=1)"),
  theta_ke_layered_summary %>% mutate(model = "Reduced-cannibalism (theta=0.3, ir=1.3)")
)
write.csv(model_comparison_df, file.path("interesting_plots", "day39_schedule_model_comparison.csv"), row.names = FALSE)

schedule_model_comparison_plot <- ggplot(model_comparison_df,
                                         aes(x = mean_effort, y = pmax(mean_total, 1e-12), color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 1.5) +
  facet_grid(model ~ knife_edge_size, labeller = label_both, scales = "free_y") +
  scale_y_log10() +
  labs(x = "Realised mean effort actually applied", y = "Mean total biomass (log scale)",
       color = "Schedule",
       title = "Does cutting cannibalism change which fishing schedule wins?",
       subtitle = "Same schedules, same knife-edge sweep, same effort range -- original paper model (top) vs. Section 2's reduced-cannibalism model (bottom)") +
  theme_minimal()
schedule_model_comparison_plot
save_plot(schedule_model_comparison_plot, "day39_schedule_model_comparison.png", width = 12, height = 7)

cat("Section 5 model comparison (original vs. reduced-cannibalism, same schedules): see day39_schedule_model_comparison.png / .csv.\n")

################################################################################
# Section 6: zoomed-in sweep, theta=0.3 alone (NOT Section 2's resource-boosted
# variant -- that one reversed theta=0.3's own good yield-curve shape, see
# Section 2's header), knife_edge=5, effort in [0,7] rather than Section 4/5's
# collapse-hunting [10,200] range -- the actual promising region.
################################################################################

zoom_ke              <- 5
zoom_baseline_effort <- 1
zoom_boost_seq       <- c(2, 3, 4, 5, 6, 7)

zoom_summary <- run_layered_ke_case(zoom_ke, theta_low, 1, last_n_theta_low, last_npp_theta_low,
                                    baseline_effort = zoom_baseline_effort, boost_seq = zoom_boost_seq)
write.csv(zoom_summary, file.path("interesting_plots", "day39_zoom_ke5_summary.csv"), row.names = FALSE)

zoom_collapse_plot <- ggplot(zoom_summary, aes(x = mean_effort, y = mean_total, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = zoom_baseline_effort, linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort", y = "Mean total biomass",
       title = "theta=0.3, knife_edge=5, effort in [0,7]: biomass vs. realised effort",
       subtitle = sprintf("floor=%.2g (dashed line)", zoom_baseline_effort)) +
  theme_minimal()
zoom_collapse_plot
save_plot(zoom_collapse_plot, "day39_zoom_ke5_collapse.png", width = 9, height = 6)

zoom_yield_plot <- ggplot(zoom_summary, aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = zoom_baseline_effort, linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort", y = "Mean yield",
       title = "theta=0.3, knife_edge=5, effort in [0,7]: yield vs. realised effort",
       subtitle = sprintf("floor=%.2g (dashed line)", zoom_baseline_effort)) +
  theme_minimal()
zoom_yield_plot
save_plot(zoom_yield_plot, "day39_zoom_ke5_yield.png", width = 9, height = 6)

zoom_amplitude_plot <- ggplot(zoom_summary, aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_vline(xintercept = zoom_baseline_effort, linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort", y = "Relative amplitude of biomass",
       title = "theta=0.3, knife_edge=5, effort in [0,7]: relative amplitude vs. realised effort",
       subtitle = sprintf("floor=%.2g (dashed line)", zoom_baseline_effort)) +
  theme_minimal()
zoom_amplitude_plot
save_plot(zoom_amplitude_plot, "day39_zoom_ke5_amplitude.png", width = 9, height = 6)

cat("Section 6 (theta=0.3, knife_edge=5, effort in [0,7]): see day39_zoom_ke5_collapse.png / day39_zoom_ke5_yield.png / day39_zoom_ke5_amplitude.png / day39_zoom_ke5_summary.csv.\n")
print(zoom_summary)

################################################################################
# Section 7: yield/relative-amplitude scan, the most fishing-sensitive model
# built so far -- theta=0.3+interaction_resource=1.3 (Section 2's own finding:
# the only variant whose total biomass declines monotonically across the
# whole [1,100] range rather than plateauing). Section 2 already scanned it
# at knife_edge=10 and 5; this fills the one gap, knife_edge=9.
################################################################################

sensitive_ke9 <- make_anchovy_fishing_params_theta(theta_low, interaction_resource_boost, knife_edge_size = 9)

sensitive_probe_df <- cannibalism_fishing_probe(sensitive_ke9, last_n_resboost, last_npp_resboost,
                                                fish_level_seq_theta) %>%
  mutate(model = "theta=0.3, interaction_resource=1.3")
write.csv(sensitive_probe_df, file.path("interesting_plots", "day39_sensitive_ke9_scan.csv"), row.names = FALSE)

sensitive_yield_plot <- ggplot(sensitive_probe_df, aes(x = fish_level, y = mean_yield)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=9)", y = "Mean yield",
       title = "Most fishing-sensitive model (theta=0.3, interaction_resource=1.3): yield vs. effort",
       subtitle = "knife_edge=9 -- the one gear setting Section 2 didn't already cover for this model") +
  theme_minimal()
sensitive_yield_plot
save_plot(sensitive_yield_plot, "day39_sensitive_ke9_yield.png", width = 9, height = 6)

sensitive_amplitude_plot <- ggplot(sensitive_probe_df, aes(x = fish_level, y = rel_amplitude)) +
  geom_line() +
  geom_point(size = 2) +
  labs(x = "Fishing effort (Constant schedule, knife_edge=9)", y = "Relative amplitude of biomass",
       title = "Most fishing-sensitive model (theta=0.3, interaction_resource=1.3): relative amplitude vs. effort",
       subtitle = "knife_edge=9") +
  theme_minimal()
sensitive_amplitude_plot
save_plot(sensitive_amplitude_plot, "day39_sensitive_ke9_amplitude.png", width = 9, height = 6)

cat("Section 7 (theta=0.3+interaction_resource=1.3, knife_edge=9): see day39_sensitive_ke9_yield.png / day39_sensitive_ke9_amplitude.png / day39_sensitive_ke9_scan.csv.\n")
print(sensitive_probe_df)

################################################################################
# Section 8: summary
################################################################################

cat("\n===== Day 39 summary =====\n")
cat("Day 38's Figure 2e plankton-anchovy cannibalism model, with its dependence on cannibalism (the species interaction matrix, theta_11) turned down and its dependence on the resource (interaction_resource) turned up.\n")
cat(sprintf(
  "Section 1 (theta=%.2g cannibalism cut, gear unchanged at knife_edge=10): day39_theta_low_growth_mort.csv / day39_theta_low_effort_probe.png / day39_theta_low_effort_yield.png / day39_theta_low_effort_amplitude.png.\n",
  theta_low
))
cat(sprintf(
  "Section 2 (theta=%.2g + interaction_resource=%.2g -- shown NOT to fix Section 1's growth deficit at w>=w_mat, see Section 2's own header): day39_theta_resboost_growth.csv / day39_theta_resboost_effort_probe.png / day39_theta_resboost_effort_yield.png / day39_theta_resboost_effort_amplitude.png / (knife_edge=5 versions also saved as *_ke5.png).\n",
  theta_low, interaction_resource_boost
))
cat(sprintf(
  "Section 3 (theta x interaction_resource grid, %dx%d cells): day39_interaction_grid.csv / day39_grid_growth_heatmap.png / day39_grid_sensitivity_heatmap.png / day39_grid_tradeoff.png -- %d candidate(s) beat the original baseline's own fishing_sensitivity while recovering over half its growth at w=10.1g, see day39_grid_candidates.csv.\n",
  length(theta_grid_seq), length(resource_grid_seq), nrow(grid_candidates)
))
cat(sprintf(
  "Section 4 (layered schedules -- Constant floor=%.2g with boosts at peaks/troughs, knife_edge=10/9/5 on Section 2's theta=0.3/interaction_resource=1.3 model, NOT the grid's own best candidate -- see Section 4's header for why): day39_theta_ke_layered_collapse.png / day39_theta_ke_layered_yield.png / day39_theta_ke_layered_amplitude.png / day39_theta_ke_boost_realised.png / day39_theta_ke_layered_summary.csv.\n",
  baseline_effort_s4
))
cat(sprintf(
  "Section 5 (the same layered schedules, on the ORIGINAL paper model theta=1/interaction_resource=1): day39_original_ke_layered_collapse.png / day39_original_ke_layered_yield.png / day39_original_ke_layered_amplitude.png / day39_original_ke_boost_realised.png / day39_original_ke_layered_summary.csv -- plus a direct side-by-side against Section 4's model, day39_schedule_model_comparison.png / .csv.\n"
))
cat(sprintf(
  "NOTE: the effort values Sections 4/5 needed to see any collapse at all (floor=%.2g, boosts up to %.2g) are nowhere near realistic fishing efforts (F is typically O(0.1-1)/yr for real fisheries) -- this is inherited from Day 38 Section 7's own finding that knife_edge=5 doesn't visibly stress until effort~30 and doesn't collapse until ~100. Section 6 zooms into a lower, more realistic range instead.\n",
  baseline_effort_s4, max(boost_fish_level_seq_s4)
))
cat(sprintf(
  "Section 6 (theta=%.2g alone, knife_edge=%.0f, effort in [0,%.0f] -- the most promising combination): day39_zoom_ke5_collapse.png / day39_zoom_ke5_yield.png / day39_zoom_ke5_amplitude.png / day39_zoom_ke5_summary.csv.\n",
  theta_low, zoom_ke, max(zoom_boost_seq)
))
cat(sprintf(
  "Section 7 (theta=%.2g+interaction_resource=%.2g -- the most fishing-sensitive model, knife_edge=9): day39_sensitive_ke9_yield.png / day39_sensitive_ke9_amplitude.png / day39_sensitive_ke9_scan.csv.\n",
  theta_low, interaction_resource_boost
))
