library(mizer)
library(dplyr)
library(ggplot2)

# Day 39: cut Day 38's cannibalism model's dependence on cannibalism (S1),
# compensate via resource interaction (S2), grid-sweep both (S3), test fishing-schedule dependence at 3 knife-edges on the S2 model and the original model (S4-5),
# then zoom in on theta=0.3 alone at knife_edge=5 over a realistic effort range (S6),
# and scan the most fishing-sensitive model (theta=0.3+ir=1.3) at the one knife-edge it hadn't covered yet (S7).
#
# Self-contained convention since Day 20: helpers redefined here, not
# sourced from 38_experiments.R.

plot_dir <- "interesting_plots"
dir.create(plot_dir, showWarnings = FALSE)

# Windows MAX_PATH truncation guard, carried over from Day 30 onward.
save_plot <- function(plot, filename, width = 9, height = 6, dpi = 150) {
  max_name <- 40
  if (nchar(filename) > max_name) {
    ext      <- tools::file_ext(filename)
    base     <- tools::file_path_sans_ext(filename)
    filename <- paste0(substr(base, 1, max_name - nchar(ext) - 1), ".", ext)
    warning(sprintf("save_plot(): filename too long, truncated to '%s'", filename))
  }
  print(plot)
  ggsave(file.path(plot_dir, filename), plot = plot,
         width = width, height = height, dpi = dpi)
}

save_table <- function(df, filename) {
  write.csv(df, file.path(plot_dir, filename), row.names = FALSE)
  df
}

# Every scan in this script draws the same picture: a summary metric against
# effort, one line per model or schedule, optionally facetted by gear.
effort_plot <- function(df, y, ylab, title, subtitle,
                        x = "mean_effort",
                        xlab = "Realised mean effort actually applied",
                        colour = "schedule", vline = NULL, facet = NULL,
                        facet_scales = "fixed", log_y = FALSE, floor_y = NULL) {
  if (!is.null(floor_y)) df[[y]] <- pmax(df[[y]], floor_y)
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour]])) +
    geom_line() +
    geom_point(size = 2) +
    labs(x = xlab, y = ylab, colour = tools::toTitleCase(colour),
         title = title, subtitle = subtitle) +
    theme_minimal()
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = "dashed", colour = "grey50")
  }
  if (!is.null(facet)) {
    p <- p + if (length(facet) == 3L) {
      facet_grid(facet, labeller = label_both, scales = facet_scales)
    } else {
      facet_wrap(facet, labeller = label_both, scales = facet_scales)
    }
  }
  if (log_y) p <- p + scale_y_log10()
  p
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
  w <- w(params)
  mu_b <- rep(0, length(w))
  mu_b[w <= p$w_s] <- (p$mu_0 * (w / p$w_min)^p$rho_b)[w < p$w_s]
  mu_s <- if (p$mu_0 > 0) min(mu_b[w <= p$w_s]) else p$mu_s
  mu_b[w >= p$w_s] <- (mu_s * (w / p$w_s)^p$rho_s)[w >= p$w_s]
  mu_b <- mu_b + p$mu_l / (1 + (w / p$w_l)^p$rho_l)

  # ext_mort<-() rather than a direct write to @mu_b: the setter marks the array
  # as manually set, so the later species_params<-() calls (which run setParams()
  # and hence setExtMort()) leave it alone instead of silently recomputing mu_b
  # from the z0 defaults. It needs the value shaped like the slot it replaces.
  mort <- ext_mort(params)
  mort[] <- mu_b
  ext_mort(params) <- mort
  params
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

# Not mizer's own box_pred_kernel(): the paper's kernel is normalised to unit
# integral in log(ppmr) and drops the self-bin, so it stays a custom function.
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

# Figure 2e config + knife-edge gear at `knife_edge_size`, catchability=1 --
# Day 36/38's own fishing convention, carried over. `interaction_val` is the
# cannibalism knob theta_11; `interaction_resource_val` is the separate
# resource-interaction knob (Section 2).
anchovy_params <- function(interaction_val = 1, interaction_resource_val = 1,
                           knife_edge_size = p2$w_mat, p = p2) {
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

  params <- newMultispeciesParams(
    species_params,
    no_w = round(log(p$w_inf / p$w_min) / p$dx),
    lambda = p$lambda,
    kappa = kappa,
    w_pp_cutoff = p$w_pp_cutoff,
    resource_dynamics = "plankton_logistic"
  )
  # plankton_logistic has no balancing function, so setResource() leaves the
  # capacity alone and just records the rate as manually set.
  resource_rate(params) <- p$r0 * w_full(params)^(p$rho - 1)

  interaction_matrix(params) <- interaction_val
  species_params(params)$interaction_resource <- interaction_resource_val

  gp <- gear_params(params)
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- knife_edge_size
  gp$catchability    <- 1
  gear_params(params) <- gp

  # Mortality last: it is the one array set by hand, and setting it through
  # ext_mort<-() protects it from everything above.
  setAnchovyMort(params, p)
}

p_scan <- anchovy_params()
anchovy_immigration <- p2$i0 * w_full(p_scan)^(-p2$lambda) * exp(-6.9 * (p2$lambda - 1))

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

# initialN<-() / initialNResource<-() want values shaped exactly like the slots
# they replace, so recycle through the existing arrays.
set_state <- function(params, n, n_pp) {
  n0 <- initialN(params)
  n0[] <- n
  initialN(params) <- n0
  npp0 <- initialNResource(params)
  npp0[] <- n_pp
  initialNResource(params) <- npp0
  params
}

seed_from <- function(params, sim) {
  set_state(params, finalN(sim), finalNResource(sim))
}

make_anchovy_fork_sim <- function(params, t_fork) {
  params <- set_state(params, 0.001 * w(params)^(-1.8), resource_capacity(params))
  settled <- project(params, t_max = 10, dt = p2$dt, progress_bar = FALSE)
  # The kick: knock the settled spectrum down by 1e7 and let the cycle grow
  # back out. t_start keeps the fork on the same clock as the settling leg.
  params <- set_state(params, finalN(settled) / 1e7, finalNResource(settled))
  project(params, t_max = t_fork - 10, t_start = 10, dt = p2$dt, t_save = 0.2,
          progress_bar = FALSE, effort = 0)
}

# The model has a single species, so every ArrayTimeBySpecies collapses to one
# time series; `t_cut` drops the transient at the front of the window.
after_cut <- function(x, sim, t_cut) unname(x[getTimes(sim) > t_cut, 1])

# Biomass over the sizes the gear actually selects. mizer's knife_edge
# selectivity is 1 above knife_edge_size and 0 below, and catchability is 1
# throughout this script, so getBiomass(min_w = ) is exactly the
# selectivity-weighted sum(F(w) B(w)) that earlier days assembled by hand.
gear_knife_edge <- function(params) gear_params(params)$knife_edge_size[[1]]

selected_biomass <- function(sim, t_cut = -Inf) {
  after_cut(getBiomass(sim, min_w = gear_knife_edge(sim@params)), sim, t_cut)
}

sim_fork <- make_anchovy_fork_sim(p_scan, t_fork)

fork_bp <- selected_biomass(sim_fork, t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "Fork check (t in [%.0f,%.0f], last %.0fyr of the %.0fyr fork): selected biomass min=%.4g max=%.4g -- should already span close to the full cycle.\n",
  t_fork - scan_summary_window, t_fork, scan_summary_window, t_fork, min(fork_bp), max(fork_bp)
))

################################################################################
# Section 1: decrease dependence on cannibalism, to make the model more
# fishing-sensitive.
#
# theta_11 (the interaction matrix) is the cannibalism knob; interaction_resource (Section 2) is separate.
# Cutting theta to 0.3 crashes growth at w=10.1/28.6g via indirect resource competition
# (juvenile population boom grazes the shared resource), not direct diet loss -- see growth_compare_df.
################################################################################

theta_low <- 0.3

p_theta_low <- anchovy_params(theta_low, 1, knife_edge_size = p2$w_mat)
sim_fork_theta_low <- make_anchovy_fork_sim(p_theta_low, t_fork)

bp_theta_low <- selected_biomass(sim_fork_theta_low, t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "Section 1 (theta=%.2g, interaction_resource=1), unfished fork, last %.0fyr of %.0fyr: selected biomass min=%.4g max=%.4g -- still cycling.\n",
  theta_low, scan_summary_window, t_fork, min(bp_theta_low), max(bp_theta_low)
))

# Growth vs. predation-mortality, baseline (theta=1) vs. theta_low, both on their OWN fork state.
w_check   <- c(0.01, 0.1, 1, 10, 30)
idx_check <- vapply(w_check, function(x) which.min(abs(w(p_scan) - x)), integer(1))

# Rates evaluated on a model's own fork state. Seeding the params with the fork
# and reading the rate off it keeps state and model together in one object.
fork_rate <- function(getter, params, sim, idx) {
  getter(seed_from(params, sim))[1, idx]
}

growth_compare_df <- data.frame(
  w                  = round(w(p_scan)[idx_check], 4),
  growth_baseline    = fork_rate(getEGrowth,  p_scan,      sim_fork,           idx_check),
  growth_theta_low   = fork_rate(getEGrowth,  p_theta_low, sim_fork_theta_low, idx_check),
  predmort_baseline  = fork_rate(getPredMort, p_scan,      sim_fork,           idx_check),
  predmort_theta_low = fork_rate(getPredMort, p_theta_low, sim_fork_theta_low, idx_check)
)
save_table(growth_compare_df, "day39_theta_low_growth_mort.csv")
cat("Section 1 growth/mortality comparison (theta=1 vs. theta=0.3), by size -- see day39_theta_low_growth_mort.csv:\n")
print(growth_compare_df)

# Extended from [1,3,5,7,9] up to 100, matching Section 4's scale -- [1,9] alone
# showed almost no curvature. Gear fixed at knife_edge=10 to isolate theta's effect.
fish_level_seq_theta <- c(1, 3, 5, 7, 9, 20, 30, 50, 75, 100)

# The one summary every scan in this script reports. `mean_selected` is biomass
# above w_mat (maturity), which is not the same thing as the gear-selected
# biomass used to drive the threshold rules below.
sim_metrics <- function(sim, t_cut) {
  total  <- after_cut(getBiomass(sim), sim, t_cut)
  mature <- after_cut(getBiomass(sim, min_w = p2$w_mat), sim, t_cut)
  yield  <- after_cut(getYield(sim), sim, t_cut)
  data.frame(mean_total = mean(total), min_total = min(total),
             mean_selected = mean(mature), mean_juvenile = mean(total - mature),
             mean_yield = mean(yield),
             rel_amplitude = (max(total) - min(total)) / ((max(total) + min(total)) / 2))
}

cannibalism_fishing_probe <- function(params, fork_sim, fish_level_seq,
                                      years = scan_summary_window) {
  params <- seed_from(params, fork_sim)
  bind_rows(lapply(fish_level_seq, function(fl) {
    sim <- project(params, t_max = years, dt = p2$dt, t_save = 0.2,
                   t_start = t_fork, progress_bar = FALSE, effort = fl)
    cbind(fish_level = fl, sim_metrics(sim, t_fork + years / 2))
  }))
}

theta_probe_df <- bind_rows(
  cannibalism_fishing_probe(p_scan, sim_fork, fish_level_seq_theta) %>%
    mutate(model = "Baseline (theta=1)"),
  cannibalism_fishing_probe(p_theta_low, sim_fork_theta_low, fish_level_seq_theta) %>%
    mutate(model = "theta=0.3")
)
save_table(theta_probe_df, "day39_theta_low_effort_probe.csv")

probe_x    <- "fish_level"
probe_xlab <- "Fishing effort (Constant schedule, knife_edge=10)"
probe_sub  <- "Gear unchanged (knife_edge=10) -- isolates theta's effect from Section 4's gear-cutoff effect"

save_plot(effort_plot(
  theta_probe_df, "mean_total", "Mean total biomass",
  "Does cutting cannibalism (theta) flatten the model's own fishing resistance?",
  probe_sub, x = probe_x, xlab = probe_xlab, colour = "model"
), "day39_theta_low_effort_probe.png")

save_plot(effort_plot(
  theta_probe_df, "mean_yield", "Mean yield",
  "Yield vs. effort -- does cutting cannibalism change where yield peaks or turns over?",
  probe_sub, x = probe_x, xlab = probe_xlab, colour = "model"
), "day39_theta_low_effort_yield.png")

# Relative amplitude alongside yield -- is the population still genuinely cycling
# at each effort level, or has yield kept climbing while the cycle itself died out?
save_plot(effort_plot(
  theta_probe_df, "rel_amplitude", "Relative amplitude of biomass",
  "Relative amplitude vs. effort -- alongside the yield plot above",
  "Gear unchanged (knife_edge=10)", x = probe_x, xlab = probe_xlab, colour = "model"
), "day39_theta_low_effort_amplitude.png")

cat("Section 1 (theta=0.3 vs. baseline effort probe, knife_edge=10): see day39_theta_low_effort_probe.png / day39_theta_low_effort_yield.png / day39_theta_low_effort_amplitude.png / .csv.\n")
print(theta_probe_df)

################################################################################
# Section 2: same cannibalism cut, with interaction_resource raised to guard
# against Section 1's growth-rate risk.
#
# interaction_resource=1.3 does NOT fix it: growth at w=10.1/28.6g barely recovers
# vs. baseline, while over-correcting the small sizes that didn't need it.
################################################################################

interaction_resource_boost <- 1.3

p_theta_low_resboost <- anchovy_params(theta_low, interaction_resource_boost,
                                       knife_edge_size = p2$w_mat)
sim_fork_resboost <- make_anchovy_fork_sim(p_theta_low_resboost, t_fork)

bp_resboost <- selected_biomass(sim_fork_resboost, t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "Section 2 (theta=%.2g, interaction_resource=%.2g), unfished fork, last %.0fyr of %.0fyr: selected biomass min=%.4g max=%.4g -- still cycling.\n",
  theta_low, interaction_resource_boost, scan_summary_window, t_fork, min(bp_resboost), max(bp_resboost)
))

growth_resboost_df <- data.frame(
  w                         = round(w(p_scan)[idx_check], 4),
  growth_baseline           = fork_rate(getEGrowth, p_scan,               sim_fork,           idx_check),
  growth_theta_low          = fork_rate(getEGrowth, p_theta_low,          sim_fork_theta_low, idx_check),
  growth_theta_low_resboost = fork_rate(getEGrowth, p_theta_low_resboost, sim_fork_resboost,  idx_check)
)
growth_resboost_df$pct_vs_baseline <- 100 * (growth_resboost_df$growth_theta_low_resboost /
                                              growth_resboost_df$growth_baseline - 1)
save_table(growth_resboost_df, "day39_theta_resboost_growth.csv")
cat("Section 2 growth comparison (baseline vs. theta=0.3 vs. theta=0.3+interaction_resource=1.3) -- see day39_theta_resboost_growth.csv:\n")
print(growth_resboost_df)

theta_probe_df_s2 <- bind_rows(
  theta_probe_df,
  cannibalism_fishing_probe(p_theta_low_resboost, sim_fork_resboost, fish_level_seq_theta) %>%
    mutate(model = "theta=0.3, interaction_resource=1.3")
)
save_table(theta_probe_df_s2, "day39_theta_resboost_effort_probe.csv")

probe_sub_s2 <- "All three models share the same gear (knife_edge=10) -- only theta/interaction_resource differ"

save_plot(effort_plot(
  theta_probe_df_s2, "mean_total", "Mean total biomass",
  "Does the resource-interaction bump change the fishing-resistance pattern?",
  probe_sub_s2, x = probe_x, xlab = probe_xlab, colour = "model"
), "day39_theta_resboost_effort_probe.png")

save_plot(effort_plot(
  theta_probe_df_s2, "mean_yield", "Mean yield",
  "Yield vs. effort, all three variants -- does the resource bump change where yield peaks?",
  probe_sub_s2, x = probe_x, xlab = probe_xlab, colour = "model"
), "day39_theta_resboost_effort_yield.png")

save_plot(effort_plot(
  theta_probe_df_s2, "rel_amplitude", "Relative amplitude of biomass",
  "Relative amplitude vs. effort, all three variants -- alongside the yield plot above",
  probe_sub_s2, x = probe_x, xlab = probe_xlab, colour = "model"
), "day39_theta_resboost_effort_amplitude.png")

cat("Section 2 (theta=0.3+interaction_resource=1.3 effort probe, knife_edge=10): see day39_theta_resboost_effort_probe.png / day39_theta_resboost_effort_yield.png / day39_theta_resboost_effort_amplitude.png / .csv.\n")
print(theta_probe_df_s2)

# Same three-way comparison, knife_edge=5 -- Day 38's own collapse-capable gear,
# vs. the knife_edge=10 isolation above. Unfished trajectory is gear-independent
# (Day 38 Section 7), so the three forks above are reused.
ke5_models <- list(
  "Baseline (theta=1)"                  = list(1,         1,                          sim_fork),
  "theta=0.3"                           = list(theta_low, 1,                          sim_fork_theta_low),
  "theta=0.3, interaction_resource=1.3" = list(theta_low, interaction_resource_boost, sim_fork_resboost)
)

theta_probe_df_ke5 <- bind_rows(lapply(names(ke5_models), function(nm) {
  spec <- ke5_models[[nm]]
  params <- anchovy_params(spec[[1]], spec[[2]], knife_edge_size = 5)
  cannibalism_fishing_probe(params, spec[[3]], fish_level_seq_theta) %>% mutate(model = nm)
}))
save_table(theta_probe_df_ke5, "day39_theta_resboost_effort_yield_ke5.csv")

probe_xlab_ke5 <- "Fishing effort (Constant schedule, knife_edge=5)"

save_plot(effort_plot(
  theta_probe_df_ke5, "mean_yield", "Mean yield",
  "Yield vs. effort at knife_edge=5, all three variants -- does gear that actually reaches the collapse zone change the picture?",
  "Day 38's own collapse-capable gear setting, vs. Sections 1/2's knife_edge=10 isolation above",
  x = probe_x, xlab = probe_xlab_ke5, colour = "model"
), "day39_theta_resboost_effort_yield_ke5.png")

save_plot(effort_plot(
  theta_probe_df_ke5, "rel_amplitude", "Relative amplitude of biomass",
  "Relative amplitude vs. effort at knife_edge=5, all three variants -- alongside the yield plot above",
  "Day 38's own collapse-capable gear setting",
  x = probe_x, xlab = probe_xlab_ke5, colour = "model"
), "day39_theta_resboost_effort_amplitude_ke5.png")

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
idx_grid_check       <- vapply(w_grid_check, function(x) which.min(abs(w(p_scan) - x)), integer(1))
growth_baseline_grid <- fork_rate(getEGrowth, p_scan, sim_fork, idx_grid_check)

run_interaction_grid_case <- function(theta_val, ir_val) {
  params        <- anchovy_params(theta_val, ir_val, knife_edge_size = p2$w_mat)
  sim_fork_grid <- make_anchovy_fork_sim(params, t_fork)
  params        <- seed_from(params, sim_fork_grid)

  bp_grid        <- selected_biomass(sim_fork_grid, t_cut = t_fork - scan_summary_window)
  unfished_ratio <- max(bp_grid) / max(min(bp_grid), 1e-12)

  # 1 = fully recovered vs. the ORIGINAL baseline
  growth_recovery <- getEGrowth(params)[1, idx_grid_check] / growth_baseline_grid

  total_unfished <- mean(after_cut(getBiomass(sim_fork_grid), sim_fork_grid,
                                   t_fork - scan_summary_window))

  sim_fished <- project(params, t_max = scan_summary_window, dt = p2$dt, t_save = 0.2,
                        t_start = t_fork, progress_bar = FALSE, effort = grid_probe_effort)
  fished <- sim_metrics(sim_fished, t_fork + scan_summary_window / 2)

  data.frame(theta = theta_val, interaction_resource = ir_val,
            unfished_cycle_ratio = unfished_ratio,
            growth_recovery_w0.01 = growth_recovery[1], growth_recovery_w1 = growth_recovery[2],
            growth_recovery_w10 = growth_recovery[3], growth_recovery_w30 = growth_recovery[4],
            total_unfished = total_unfished, total_at_effort10 = fished$mean_total,
            mean_yield_at_effort10 = fished$mean_yield,
            # fraction of unfished biomass lost
            fishing_sensitivity = 1 - fished$mean_total / total_unfished)
}

interaction_grid_df <- bind_rows(lapply(theta_grid_seq, function(th) {
  bind_rows(lapply(resource_grid_seq, function(ir) run_interaction_grid_case(th, ir)))
}))
save_table(interaction_grid_df, "day39_interaction_grid.csv")
cat("Section 3 (theta x interaction_resource grid): see day39_interaction_grid.csv.\n")
print(interaction_grid_df)

grid_heatmap <- function(fill, fill_lab, title, subtitle, option = "viridis") {
  ggplot(interaction_grid_df,
         aes(x = factor(theta), y = factor(interaction_resource), fill = .data[[fill]])) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", .data[[fill]])), size = 3, colour = "white") +
    scale_fill_viridis_c(option = option, limits = c(0, NA)) +
    labs(x = "theta (cannibalism)", y = "interaction_resource", fill = fill_lab,
         title = title, subtitle = subtitle) +
    theme_minimal()
}

save_plot(grid_heatmap(
  "growth_recovery_w10", "recovery", "Growth recovery at w=10.1g",
  "1.0 = fully back to the original theta=1/interaction_resource=1 baseline's own growth there"
), "day39_grid_growth_heatmap.png", width = 8, height = 6)

save_plot(grid_heatmap(
  "fishing_sensitivity", "sensitivity",
  "Fishing sensitivity: fraction of unfished biomass lost at effort=10",
  "Higher = more responsive to fishing -- Section 1's original goal", option = "magma"
), "day39_grid_sensitivity_heatmap.png", width = 8, height = 6)

grid_tradeoff_plot <- ggplot(interaction_grid_df,
                             aes(x = fishing_sensitivity, y = growth_recovery_w10,
                                 colour = factor(theta), shape = factor(interaction_resource))) +
  geom_point(size = 3) +
  labs(x = "Fishing sensitivity (fraction of unfished biomass lost at effort=10)",
       y = "Growth recovery at w=10.1g (1.0 = original baseline)",
       colour = "theta", shape = "interaction_resource",
       title = "The trade-off this section is actually navigating",
       subtitle = "Top-right is the ideal quadrant: growth restored AND still fishing-sensitive") +
  theme_minimal()
save_plot(grid_tradeoff_plot, "day39_grid_tradeoff.png")

baseline_grid_row <- interaction_grid_df %>% filter(theta == 1, interaction_resource == 1)
cat(sprintf(
  "Section 3 baseline cell (theta=1, interaction_resource=1): fishing_sensitivity=%.3f at effort=10 -- the bar Section 1 was trying to clear.\n",
  baseline_grid_row$fishing_sensitivity
))

grid_candidates <- interaction_grid_df %>%
  filter(growth_recovery_w10 > 0.5, fishing_sensitivity > baseline_grid_row$fishing_sensitivity) %>%
  arrange(desc(growth_recovery_w10))
save_table(grid_candidates, "day39_grid_candidates.csv")
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

baseline_effort_s4      <- 10
# Extended to 150/200 -- knife_edge=10/9 plateaued rather than collapsed at the old max of 100.
boost_fish_level_seq_s4 <- c(20, 30, 50, 75, 100, 150, 200)
knife_edge_seq_s4       <- c(10, 9, 5)

# Day 36/38's own threshold-rule machinery. The "on" fraction is a soft switch on
# the gear-selected biomass; it is shared by the rate function and by the
# after-the-fact diagnostics so the two can never drift apart.
on_fraction <- function(selected_biomass, threshold, mode, sharpness,
                        hard_step = FALSE) {
  direction <- if (mode == "above") 1 else -1
  x <- direction * (selected_biomass - threshold)
  if (isTRUE(hard_step)) as.numeric(x > 0) else plogis(x / max(sharpness, 1e-12))
}

thresholdFMort <- function(params, n, n_pp, n_other, t, effort, e_growth, pred_mort, ...) {
  p <- other_params(params)
  on_frac <- on_fraction(sum(p$biomass_weight * n),
                         p$threshold, p$mode, p$sharpness, p$hard_step)
  result <- (on_frac * p$fish_level + (1 - on_frac) * p$background_level) * p$f_ref
  dim(result)      <- dim(n)
  dimnames(result) <- dimnames(n)
  result
}

attach_threshold_rule <- function(params, threshold, fish_level, background_level = 0,
                                  mode = c("above", "below"), sharpness, hard_step = FALSE) {
  mode <- match.arg(mode)
  # mizerFMortGear() is linear in effort, so the unit-effort mortality can be
  # computed once here instead of twice per time step inside the rate function.
  f_ref <- colSums(mizerFMortGear(params, effort = 1))
  other_params(params) <- list(threshold = threshold, fish_level = fish_level,
                               background_level = background_level, mode = mode,
                               sharpness = sharpness, hard_step = hard_step,
                               f_ref = f_ref,
                               biomass_weight = sweep(f_ref, 2, w(params) * dw(params), "*"))
  setRateFunction(params, "FMort", "thresholdFMort")
}

threshold_diagnostics <- function(sim, threshold, fish_level, background_level, mode,
                                  sharpness, t_cut, hard_step = FALSE) {
  tv      <- getTimes(sim)
  on_frac <- on_fraction(selected_biomass(sim, t_cut), threshold, mode, sharpness, hard_step)

  runs   <- rle(on_frac > 0.5)
  n_runs <- length(runs$values)

  # First/last run in the window are censored -- see Day 36's own note.
  interior <- rep(TRUE, n_runs)
  if (n_runs > 0) { interior[1] <- FALSE; interior[n_runs] <- FALSE }
  on_lengths <- runs$lengths[runs$values & interior]
  dt_save    <- mean(diff(tv[tv > t_cut]))

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

p_sanity <- seed_from(anchovy_params(theta_low, interaction_resource_boost,
                                     knife_edge_size = sanity_ke),
                      sim_fork_resboost)
sim_const_sanity <- project(p_sanity, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                            t_start = t_fork, progress_bar = FALSE, effort = baseline_effort_s4)
bp_const_sanity  <- selected_biomass(sim_const_sanity, scan_t_cut)
sharpness_sanity <- 0.02 * (max(bp_const_sanity) - min(bp_const_sanity))
threshold_sanity <- unname(quantile(bp_const_sanity, probs = 0.5))

sanity_modes <- c("Constant" = NA_character_,
                  "Threshold (peaks)" = "above",
                  "Threshold (troughs)" = "below")

sanity_series_df <- bind_rows(lapply(names(sanity_modes), function(nm) {
  mode <- sanity_modes[[nm]]
  sim <- if (is.na(mode)) {
    sim_const_sanity
  } else {
    p_rule <- attach_threshold_rule(p_sanity, threshold = threshold_sanity,
                                    fish_level = sanity_boost,
                                    background_level = baseline_effort_s4, mode = mode,
                                    sharpness = sharpness_sanity)
    project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
            t_start = t_fork, progress_bar = FALSE, effort = 1)
  }
  bp <- selected_biomass(sim, t_fork)
  data.frame(t = getTimes(sim)[getTimes(sim) > t_fork],
             selected_biomass = bp,
             on_frac = if (is.na(mode)) NA_real_ else
               on_fraction(bp, threshold_sanity, mode, sharpness_sanity),
             schedule = nm)
}))
save_table(sanity_series_df, "day39_schedule_sanity_series.csv")

sanity_check_plot <- ggplot(sanity_series_df, aes(x = t, y = selected_biomass)) +
  geom_rect(data = sanity_series_df %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = threshold_sanity, linetype = "dashed", colour = "grey40") +
  facet_wrap(~schedule, ncol = 1, scales = "free_y") +
  labs(x = "Time (years)", y = "Selected biomass",
       title = "Sanity check: does thresholdFMort() actually fire where it should?",
       subtitle = sprintf("knife_edge=%.0f, floor=%.2g, boost=%.2g -- dashed = calibration threshold, red = boosted 'on'",
                          sanity_ke, baseline_effort_s4, sanity_boost)) +
  theme_minimal()
save_plot(sanity_check_plot, "day39_schedule_sanity_check.png", width = 9, height = 8)

cat(sprintf(
  "Schedule sanity check (knife_edge=%.0f, boost=%.2g): see day39_schedule_sanity_check.png / day39_schedule_sanity_series.csv.\n",
  sanity_ke, sanity_boost
))

# Generalised over (theta_val, ir_val, fork_sim) so Section 5 can rerun this
# verbatim on the original model -- only the model differs, not the schedule code.
# baseline_effort/boost_seq default to Sections 4/5's own floor+boost range,
# but are parameters (not globals) so Section 6 can reuse this with a
# different, much lower range instead of duplicating the whole function.
run_layered_ke_case <- function(knife_edge_size, theta_val, ir_val, fork_sim,
                                baseline_effort = baseline_effort_s4,
                                boost_seq = boost_fish_level_seq_s4) {
  params <- seed_from(anchovy_params(theta_val, ir_val,
                                     knife_edge_size = knife_edge_size), fork_sim)

  # Constant = the shared floor every other schedule is layered on top of; its own
  # cycle also calibrates the threshold/sharpness used below.
  sim_const <- project(params, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                       t_start = t_fork, progress_bar = FALSE, effort = baseline_effort)
  bp_const     <- selected_biomass(sim_const, scan_t_cut)
  sharpness_ke <- 0.02 * (max(bp_const) - min(bp_const))
  threshold_ke <- unname(quantile(bp_const, probs = 0.5))

  const_row <- cbind(sim_metrics(sim_const, scan_t_cut),
                     knife_edge_size = knife_edge_size, schedule = "Constant",
                     boost_level = NA_real_, mean_effort = baseline_effort,
                     effective_window = NA_real_, n_bursts = NA_integer_)

  boosted_rows <- bind_rows(lapply(boost_seq, function(boost) {
    bind_rows(lapply(c("above", "below"), function(mode) {
      p_rule <- attach_threshold_rule(params, threshold = threshold_ke, fish_level = boost,
                                      background_level = baseline_effort, mode = mode,
                                      sharpness = sharpness_ke)
      sim <- project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                     t_start = t_fork, progress_bar = FALSE, effort = 1)
      diag <- threshold_diagnostics(sim, threshold_ke, boost, baseline_effort, mode,
                                    sharpness_ke, scan_t_cut)
      cbind(sim_metrics(sim, scan_t_cut),
           knife_edge_size = knife_edge_size,
           schedule = if (mode == "above") "Threshold (peaks)" else "Threshold (troughs)",
           boost_level = boost, mean_effort = diag$mean_effort,
           effective_window = diag$effective_window, n_bursts = diag$n_bursts)
    }))
  }))

  bind_rows(const_row, boosted_rows)
}

# Sections 4 and 5 run the identical schedule sweep on two different models, then
# plot the identical four figures, so both are driven from one description.
layered_sweep <- function(theta_val, ir_val, fork_sim) {
  bind_rows(lapply(knife_edge_seq_s4, function(ke) {
    run_layered_ke_case(ke, theta_val, ir_val, fork_sim)
  }))
}

layered_plots <- function(summary_df, prefix, title_stem, subtitle_stem) {
  floor_note <- sprintf("baseline_effort=%.2g (dashed line) -- every schedule fishes at least this much, always.%s",
                        baseline_effort_s4, subtitle_stem)

  save_plot(effort_plot(
    summary_df, "mean_total", "Mean total biomass (log scale)",
    paste0(title_stem, "Constant floor vs. boosting at peaks vs. boosting at troughs -- same floor throughout"),
    floor_note, vline = baseline_effort_s4, facet = ~knife_edge_size,
    log_y = TRUE, floor_y = 1e-12
  ), paste0(prefix, "_layered_collapse.png"), width = 11, height = 5)

  save_plot(effort_plot(
    summary_df, "mean_yield", "Mean yield",
    paste0(title_stem, "yield vs. realised effort, same floor throughout"),
    sprintf("baseline_effort=%.2g (dashed line).%s", baseline_effort_s4, subtitle_stem),
    vline = baseline_effort_s4, facet = ~knife_edge_size, facet_scales = "free_y"
  ), paste0(prefix, "_layered_yield.png"), width = 11, height = 5)

  save_plot(effort_plot(
    summary_df, "rel_amplitude", "Relative amplitude of biomass",
    paste0(title_stem, "relative amplitude vs. realised effort -- alongside the yield plot above"),
    sprintf("baseline_effort=%.2g (dashed line).%s", baseline_effort_s4, subtitle_stem),
    vline = baseline_effort_s4, facet = ~knife_edge_size
  ), paste0(prefix, "_layered_amplitude.png"), width = 11, height = 5)

  # Does "Threshold (peaks)" still self-limit (Day 38's finding) now there's a nonzero floor?
  boost_plot <- effort_plot(
    summary_df %>% filter(schedule != "Constant"), "mean_effort", "Realised mean effort",
    paste0(title_stem, "how much of the nominal boost actually gets applied, on average?"),
    "Dotted line = boost applied in full, all the time (upper bound); floor=10 is the lower bound",
    x = "boost_level", xlab = "Nominal boost (on-period fish_level)", facet = ~knife_edge_size
  ) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60")
  save_plot(boost_plot, paste0(prefix, "_boost_realised.png"), width = 11, height = 5)
}

# Correlation of realised mean_effort with mean_total: negative = fishing-dependent
# / collapsing, positive or NA = still resistant or too few points.
layered_trend <- function(summary_df) {
  summary_df %>%
    group_by(knife_edge_size, schedule) %>%
    summarise(total_trend = cor(mean_effort, mean_total), .groups = "drop") %>%
    arrange(desc(knife_edge_size), schedule)
}

theta_ke_layered_summary <- layered_sweep(theta_low, interaction_resource_boost,
                                          sim_fork_resboost)
save_table(theta_ke_layered_summary, "day39_theta_ke_layered_summary.csv")
layered_plots(theta_ke_layered_summary, "day39_theta_ke", "",
              " Reduced-cannibalism model (theta=0.3, interaction_resource=1.3).")

cat(sprintf(
  "Section 4 (layered schedules -- Constant floor=%.2g, boosted at peaks/troughs, knife_edge=10/9/5): see day39_theta_ke_layered_collapse.png / day39_theta_ke_layered_yield.png / day39_theta_ke_layered_amplitude.png / day39_theta_ke_boost_realised.png / day39_theta_ke_layered_summary.csv.\n",
  baseline_effort_s4
))
print(theta_ke_layered_summary)

cat("Correlation of realised mean_effort with mean_total, by knife_edge_size and schedule (negative = fishing-dependent/collapsing, positive/NA = still resistant or too few points):\n")
print(layered_trend(theta_ke_layered_summary))

################################################################################
# Section 5: same layered schedules, on the ORIGINAL paper model (theta=1,
# interaction_resource=1) -- Section 4 only tested Section 2's modified model,
# and Sections 1-3 found theta/interaction_resource change fishing response a lot,
# so the original model's own schedule comparison can't just be inferred.
################################################################################

original_ke_layered_summary <- layered_sweep(1, 1, sim_fork)
save_table(original_ke_layered_summary, "day39_original_ke_layered_summary.csv")
layered_plots(original_ke_layered_summary, "day39_original_ke",
              "Original paper model (theta=1, interaction_resource=1): ", "")

cat("Section 5 (layered schedules on the ORIGINAL model, theta=1/interaction_resource=1): see day39_original_ke_layered_collapse.png / day39_original_ke_layered_yield.png / day39_original_ke_layered_amplitude.png / day39_original_ke_boost_realised.png / day39_original_ke_layered_summary.csv.\n")
print(original_ke_layered_summary)

cat("Correlation of realised mean_effort with mean_total, ORIGINAL model, by knife_edge_size and schedule:\n")
print(layered_trend(original_ke_layered_summary))

# Direct side-by-side: does cutting cannibalism change which schedule wins, or just shift the picture?
model_comparison_df <- bind_rows(
  original_ke_layered_summary %>% mutate(model = "Original (theta=1, ir=1)"),
  theta_ke_layered_summary %>% mutate(model = "Reduced-cannibalism (theta=0.3, ir=1.3)")
)
save_table(model_comparison_df, "day39_schedule_model_comparison.csv")

save_plot(effort_plot(
  model_comparison_df, "mean_total", "Mean total biomass (log scale)",
  "Does cutting cannibalism change which fishing schedule wins?",
  "Same schedules, same knife-edge sweep, same effort range -- original paper model (top) vs. Section 2's reduced-cannibalism model (bottom)",
  facet = model ~ knife_edge_size, facet_scales = "free_y",
  log_y = TRUE, floor_y = 1e-12
), "day39_schedule_model_comparison.png", width = 12, height = 7)

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

zoom_summary <- run_layered_ke_case(zoom_ke, theta_low, 1, sim_fork_theta_low,
                                    baseline_effort = zoom_baseline_effort,
                                    boost_seq = zoom_boost_seq)
save_table(zoom_summary, "day39_zoom_ke5_summary.csv")

zoom_sub <- sprintf("floor=%.2g (dashed line)", zoom_baseline_effort)
for (metric in list(
  c("mean_total",    "Mean total biomass",              "biomass",           "collapse"),
  c("mean_yield",    "Mean yield",                      "yield",             "yield"),
  c("rel_amplitude", "Relative amplitude of biomass",   "relative amplitude", "amplitude")
)) {
  save_plot(effort_plot(
    zoom_summary, metric[1], metric[2],
    sprintf("theta=0.3, knife_edge=5, effort in [0,7]: %s vs. realised effort", metric[3]),
    zoom_sub, xlab = "Realised mean effort", vline = zoom_baseline_effort
  ), sprintf("day39_zoom_ke5_%s.png", metric[4]))
}

cat("Section 6 (theta=0.3, knife_edge=5, effort in [0,7]): see day39_zoom_ke5_collapse.png / day39_zoom_ke5_yield.png / day39_zoom_ke5_amplitude.png / day39_zoom_ke5_summary.csv.\n")
print(zoom_summary)

################################################################################
# Section 7: yield/relative-amplitude scan, the most fishing-sensitive model
# built so far -- theta=0.3+interaction_resource=1.3 (Section 2's own finding:
# the only variant whose total biomass declines monotonically across the
# whole [1,100] range rather than plateauing). Section 2 already scanned it
# at knife_edge=10 and 5; this fills the one gap, knife_edge=9.
################################################################################

sensitive_probe_df <- cannibalism_fishing_probe(
  anchovy_params(theta_low, interaction_resource_boost, knife_edge_size = 9),
  sim_fork_resboost, fish_level_seq_theta
) %>% mutate(model = "theta=0.3, interaction_resource=1.3")
save_table(sensitive_probe_df, "day39_sensitive_ke9_scan.csv")

sensitive_xlab <- "Fishing effort (Constant schedule, knife_edge=9)"

save_plot(effort_plot(
  sensitive_probe_df, "mean_yield", "Mean yield",
  "Most fishing-sensitive model (theta=0.3, interaction_resource=1.3): yield vs. effort",
  "knife_edge=9 -- the one gear setting Section 2 didn't already cover for this model",
  x = probe_x, xlab = sensitive_xlab, colour = "model"
), "day39_sensitive_ke9_yield.png")

save_plot(effort_plot(
  sensitive_probe_df, "rel_amplitude", "Relative amplitude of biomass",
  "Most fishing-sensitive model (theta=0.3, interaction_resource=1.3): relative amplitude vs. effort",
  "knife_edge=9", x = probe_x, xlab = sensitive_xlab, colour = "model"
), "day39_sensitive_ke9_amplitude.png")

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
cat("Section 5 (the same layered schedules, on the ORIGINAL paper model theta=1/interaction_resource=1): day39_original_ke_layered_collapse.png / day39_original_ke_layered_yield.png / day39_original_ke_layered_amplitude.png / day39_original_ke_boost_realised.png / day39_original_ke_layered_summary.csv -- plus a direct side-by-side against Section 4's model, day39_schedule_model_comparison.png / .csv.\n")
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
