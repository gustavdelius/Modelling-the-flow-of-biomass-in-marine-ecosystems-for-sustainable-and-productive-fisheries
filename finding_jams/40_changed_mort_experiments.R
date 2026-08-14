library(mizer)
library(dplyr)
library(ggplot2)

# Day 40 (mortality follow-up): Days 39/40's own
# anchovy_params()/make_anchovy_fishing_params_theta() was silently discarding
# the paper's own mortality curve.
#
# setAnchovyMort() set the custom senescence/background mortality by writing
# straight to @mu_b. Every theta/interaction_resource variant then called
# species_params(params) <- sp to set interaction_resource -- which runs
# setParams() and hence setExtMort() internally, and SILENTLY RECOMPUTES
# @mu_b from mizer's own z0 defaults, discarding whatever setAnchovyMort()
# had just written. Confirmed directly in Section 0b below: every model built
# through that code path (Days 39-40's baseline AND every theta variant) had
# a FLAT mu_b=0.148 at every size, not the paper's own two-part
# larval/senescence curve. The fix (this file, ported from the CURRENT
# 39_experiments.R, which already carries it): set mortality LAST, through
# ext_mort<-() rather than a direct @mu_b write -- the setter marks the array
# as manually set, so later species_params<-() calls leave it alone.
#
# Self-contained convention since Day 20: helpers redefined here. Uses the
# modern accessor-function API (w(), ext_mort<-(), gear_params<-(), etc.),
# matching the CURRENT 39_experiments.R rather than Day 40's own @-slot
# style, since that's the style the fix itself is written in.

plot_dir <- "interesting_plots"
dir.create(plot_dir, showWarnings = FALSE)

save_plot <- function(plot, filename, width = 9, height = 6, dpi = 150) {
  max_name <- 40
  if (nchar(filename) > max_name) {
    ext      <- tools::file_ext(filename)
    base     <- tools::file_path_sans_ext(filename)
    filename <- paste0(substr(base, 1, max_name - nchar(ext) - 1), ".", ext)
    warning(sprintf("save_plot(): filename too long, truncated to '%s'", filename))
  }
  print(plot)
  ggsave(file.path(plot_dir, filename), plot = plot, width = width, height = height, dpi = dpi)
}

################################################################################
# Section 0: the fixed model builder, ported verbatim from the CURRENT
# 39_experiments.R (not the version this project's earlier days -- including
# all of Day 40 -- were actually built against).
################################################################################

p2 <- list(
  dt = 0.001, dx = 0.1, w_min = 0.0003, w_inf = 66.5,
  ppmr_min = 100, ppmr_max = 30000, gamma = 750, alpha = 0.85, K = 0.1,
  mu_l = 0, w_l = 0.03, rho_l = 5,
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
  # as manually set, so the later species_params<-() call below (interaction_resource)
  # leaves it alone instead of silently recomputing mu_b from the z0 defaults.
  mort <- ext_mort(params)
  mort[] <- mu_b
  ext_mort(params) <- mort
  params
}

plankton_state <- new.env(parent = emptyenv())
plankton_state$time   <- 0
plankton_state$factor <- 1

# Exact step for
#   dn/dt = immigration + (rate - mortality) * n - rate / capacity * n^2.
# This is the logistic resource equation with constant immigration. As in
# mizer::resource_logistic(), rates are held fixed during each time step.
# Replaces the old explicit-Euler step (n_pp + dt*f) -- pulled in from
# upstream/GitHub, not written for this project. Unlike the Euler version,
# this is exact for any dt, so it no longer forces dt to stay small for the
# RESOURCE side specifically; whether the CONSUMER side (still the default
# first-order upwind scheme, no second_order_w/tr_bdf2 anywhere in this
# project) can also tolerate a larger dt is untested -- see the Overnight
# section's own note before assuming dt can just be raised project-wide.
exact_logistic_immigration_step <- function(n, rate, capacity, immigration,
                                            mortality, dt) {
  result <- n
  active <- is.finite(n) & is.finite(rate) & is.finite(capacity) &
    is.finite(immigration) & is.finite(mortality) &
    rate > 0 & capacity > 0 & immigration >= 0
  if (dt == 0 || !any(active)) return(result)

  n0 <- n[active]
  r  <- rate[active]
  k  <- capacity[active]
  i  <- immigration[active]
  mu <- mortality[active]
  a  <- r - mu
  b  <- r / k
  next_n <- numeric(length(n0))

  # With immigration the quadratic has one positive and one negative root.
  # The alternative root formulae avoid cancellation when abs(a) is large.
  has_immigration <- i > 0
  if (any(has_immigration)) {
    idx <- which(has_immigration)
    ai  <- a[idx]
    bi  <- b[idx]
    ii  <- i[idx]
    d   <- sqrt(ai^2 + 4 * bi * ii)
    n_plus <- ifelse(ai >= 0, (ai + d) / (2 * bi), 2 * ii / (d - ai))
    n_minus <- ifelse(ai <= 0, (ai - d) / (2 * bi), -2 * ii / (ai + d))
    ratio <- ((n0[idx] - n_plus) / (n0[idx] - n_minus)) * exp(-d * dt)
    next_n[idx] <- (n_plus - ratio * n_minus) / (1 - ratio)
  }

  # The zero-immigration limit is handled separately, including rate == mortality.
  no_immigration <- !has_immigration
  if (any(no_immigration)) {
    idx <- which(no_immigration)
    az  <- a[idx]
    bz  <- b[idx]
    nz  <- n0[idx]
    value <- numeric(length(idx))
    zero_rate <- az == 0
    value[zero_rate] <- nz[zero_rate] /
      (1 + bz[zero_rate] * nz[zero_rate] * dt)

    positive <- az > 0 & nz > 0
    phi <- -expm1(-az[positive] * dt) / az[positive]
    value[positive] <- nz[positive] /
      (exp(-az[positive] * dt) +
         bz[positive] * nz[positive] * phi)

    negative <- az < 0
    exp_adt <- exp(az[negative] * dt)
    phi <- expm1(az[negative] * dt) / az[negative]
    value[negative] <- nz[negative] * exp_adt /
      (1 + bz[negative] * nz[negative] * phi)
    next_n[idx] <- value
  }

  # Round-off can only create tiny negative values; the analytic solution is
  # non-negative for non-negative initial abundance and immigration.
  result[active] <- pmax(next_n, 0)
  result
}

plankton_logistic <- function(params, n, n_pp, n_other, rates, dt = 0.1, ...) {
  plankton_state$time <- plankton_state$time + dt
  exact_logistic_immigration_step(
    n = n_pp,
    rate = params@rr_pp,
    capacity = params@cc_pp * plankton_state$factor,
    immigration = anchovy_immigration,
    mortality = rates$resource_mort,
    dt = dt
  )
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

anchovy_params <- function(interaction_val = 1, interaction_resource_val = 1,
                           knife_edge_size = p2$w_mat, p = p2) {
  kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))

  species_params_df <- data.frame(
    species = "Anchovy", w_min = p$w_min, w_mat = p$w_mat, m = p$rho_inf + 2/3,
    w_inf = p$w_inf, erepro = p$epsilon_R, alpha = p$K, ks = 0, gamma = p$gamma,
    q = p$alpha, ppmr_min = p$ppmr_min, ppmr_max = p$ppmr_max,
    pred_kernel_type = "norm_box", h = Inf, R_max = Inf, linecolour = "brown",
    stringsAsFactors = FALSE
  )

  params <- newMultispeciesParams(
    species_params_df, no_w = round(log(p$w_inf / p$w_min) / p$dx),
    lambda = p$lambda, kappa = kappa, w_pp_cutoff = p$w_pp_cutoff,
    resource_dynamics = "plankton_logistic"
  )
  resource_rate(params) <- p$r0 * w_full(params)^(p$rho - 1)

  interaction_matrix(params) <- interaction_val
  species_params(params)$interaction_resource <- interaction_resource_val

  gp                 <- gear_params(params)
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- knife_edge_size
  gp$catchability    <- 1
  gear_params(params) <- gp

  # Mortality last -- protected via ext_mort<-() from the species_params<-()
  # call above, which is the actual source of the Day 39-40 bug.
  setAnchovyMort(params, p)
}

p_scan <- anchovy_params()
anchovy_immigration <- p2$i0 * w_full(p_scan)^(-p2$lambda) * exp(-6.9 * (p2$lambda - 1))

################################################################################
# Section 0b: confirm the bug directly -- compare mu_b from the OLD
# (buggy-ordering, direct @mu_b write) construction against the NEW
# (ext_mort<-(), mortality-last) construction, at identical parameters.
################################################################################

theta_low <- 0.3

# The OLD, buggy-ordering builder, reproduced exactly as it appeared in
# 40_experiments.R -- kept ONLY for this one comparison, not used elsewhere
# in this file.
make_anchovy_fishing_params_theta_BUGGY <- function(interaction_val, interaction_resource_val,
                                                     knife_edge_size = p2$w_mat, p = p2) {
  params <- anchovy_params_base_only(p)
  params <- setAnchovyMort_BUGGY(params, p)          # mortality set FIRST, direct @mu_b write
  params@interaction[] <- interaction_val
  sp <- species_params(params)
  sp$interaction_resource <- interaction_resource_val
  species_params(params) <- sp                        # <-- silently resets @mu_b here
  gp                 <- params@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- knife_edge_size
  gp$catchability    <- 1
  gear_params(params) <- gp
  params
}
setAnchovyMort_BUGGY <- function(params, p) {
  mu_b <- rep(0, length(params@w))
  mu_b[params@w <= p$w_s] <- (p$mu_0 * (params@w / p$w_min)^p$rho_b)[params@w < p$w_s]
  mu_s <- if (p$mu_0 > 0) min(mu_b[params@w <= p$w_s]) else p$mu_s
  mu_b[params@w >= p$w_s] <- (mu_s * (params@w / p$w_s)^p$rho_s)[params@w >= p$w_s]
  mu_b <- mu_b + p$mu_l / (1 + (params@w / p$w_l)^p$rho_l)
  params@mu_b[] <- mu_b
  params
}
anchovy_params_base_only <- function(p) {
  kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))
  species_params_df <- data.frame(
    species = "Anchovy", w_min = p$w_min, w_mat = p$w_mat, m = p$rho_inf + 2/3,
    w_inf = p$w_inf, erepro = p$epsilon_R, alpha = p$K, ks = 0, gamma = p$gamma,
    q = p$alpha, ppmr_min = p$ppmr_min, ppmr_max = p$ppmr_max,
    pred_kernel_type = "norm_box", h = Inf, R_max = Inf, linecolour = "brown",
    stringsAsFactors = FALSE
  )
  newMultispeciesParams(species_params_df, no_w = round(log(p$w_inf / p$w_min) / p$dx),
                        lambda = p$lambda, kappa = kappa, w_pp_cutoff = p$w_pp_cutoff,
                        resource_dynamics = "plankton_logistic")
}

p_theta_low_BUGGY <- make_anchovy_fishing_params_theta_BUGGY(theta_low, 1, knife_edge_size = 10)
p_theta_low_FIXED <- anchovy_params(theta_low, 1, knife_edge_size = 10)

w_check <- c(0.0003, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 30, 66)
idx_check <- vapply(w_check, function(x) which.min(abs(w(p_theta_low_FIXED) - x)), integer(1))

mu_b_compare_df <- data.frame(
  w              = round(w(p_theta_low_FIXED)[idx_check], 4),
  mu_b_OLD_buggy = p_theta_low_BUGGY@mu_b[idx_check],
  mu_b_NEW_fixed = ext_mort(p_theta_low_FIXED)[1, idx_check]
)
write.csv(mu_b_compare_df, file.path(plot_dir, "day40_changed_mort_mu_b_compare.csv"), row.names = FALSE)
cat("Section 0b: OLD (buggy) mu_b is FLAT (mizer's own z0 default) at every size; NEW (fixed) mu_b is the paper's own U-shaped curve. See day40_changed_mort_mu_b_compare.csv:\n")
print(mu_b_compare_df)

################################################################################
# Section 0c: does the fix change the model's cycle period? Day 40's own
# t_fork=20/scan_summary_window=12 convention assumed a ~5-6yr period --
# worth actually checking under the corrected mortality rather than assuming
# it still holds. (An earlier, undersampled look at this trajectory wrongly
# suggested a ~16yr period -- an aliasing artefact of sampling only every 4th
# year; the properly-resolved check below supersedes that.)
################################################################################

set_state <- function(params, n, n_pp) {
  n0 <- initialN(params); n0[] <- n; initialN(params) <- n0
  npp0 <- initialNResource(params); npp0[] <- n_pp; initialNResource(params) <- npp0
  params
}
seed_from <- function(params, sim) set_state(params, finalN(sim), finalNResource(sim))

make_anchovy_fork_sim <- function(params, t_fork) {
  params <- set_state(params, 0.001 * w(params)^(-1.8), resource_capacity(params))
  settled <- project(params, t_max = 10, dt = p2$dt, progress_bar = FALSE)
  params <- set_state(params, finalN(settled) / 1e7, finalNResource(settled))
  project(params, t_max = t_fork - 10, t_start = 10, dt = p2$dt, t_save = 0.2,
         progress_bar = FALSE, effort = 0)
}

after_cut <- function(x, sim, t_cut) unname(x[getTimes(sim) > t_cut, 1])
gear_knife_edge <- function(params) gear_params(params)$knife_edge_size[[1]]
selected_biomass <- function(sim, t_cut = -Inf) {
  after_cut(getBiomass(sim, min_w = gear_knife_edge(sim@params)), sim, t_cut)
}

t_fork               <- 20
scan_post_fork_years <- 24
scan_summary_window  <- 12
scan_t_cut           <- t_fork + scan_post_fork_years - scan_summary_window

sim_fork_theta_low <- make_anchovy_fork_sim(p_theta_low_FIXED, t_fork)

# Extend well past t_fork to check the true, long-run period rather than
# assuming Day 37-40's ~5-6yr figure still applies.
sim_long <- project(sim_fork_theta_low, t_max = 60, t_start = t_fork, dt = p2$dt, t_save = 0.5,
                    progress_bar = FALSE, effort = 0)
bp_long <- selected_biomass(sim_long, t_cut = -Inf)
tv_long <- getTimes(sim_long)
df_long <- data.frame(t = tv_long, b = bp_long) %>% filter(t > t_fork)

n_long   <- nrow(df_long)
is_peak  <- c(FALSE, df_long$b[2:(n_long-1)] > df_long$b[1:(n_long-2)] & df_long$b[2:(n_long-1)] > df_long$b[3:n_long], FALSE)
peak_times <- df_long$t[is_peak]

cycle_period_df <- data.frame(peak_t = peak_times, peak_b = df_long$b[is_peak],
                              period_since_last = c(NA, diff(peak_times)))
write.csv(cycle_period_df, file.path(plot_dir, "day40_changed_mort_cycle_period.csv"), row.names = FALSE)
cat(sprintf(
  "Section 0c: fixed-mortality theta=0.3, unfished, peak spacing over t=[%.0f,%.0f]: mean period=%.2fyr (range %.1f-%.1f) -- close to Day 37-40's own ~5-6yr assumption, so t_fork=20/scan_summary_window=12 are kept unchanged rather than re-derived from scratch.\n",
  t_fork, max(tv_long), mean(diff(peak_times), na.rm = TRUE), min(diff(peak_times), na.rm = TRUE), max(diff(peak_times), na.rm = TRUE)
))
print(cycle_period_df)

################################################################################
# Section 1: redo Day 40 Section 1's headline comparison (Baseline theta=1 vs.
# theta=0.3, Constant fishing, knife_edge=10) under the corrected mortality.
#
# Moderate 8-point grid (not Day 40's ultra-fine 14-point one) -- this is a
# first-pass redo to establish whether the fix changes the STORY, not a
# final, fully-resolved replacement for Day 40's own sweep.
################################################################################

sim_metrics <- function(sim, t_cut) {
  total  <- after_cut(getBiomass(sim), sim, t_cut)
  mature <- after_cut(getBiomass(sim, min_w = p2$w_mat), sim, t_cut)
  yield  <- after_cut(getYield(sim), sim, t_cut)
  data.frame(mean_total = mean(total), min_total = min(total),
            mean_selected = mean(mature), mean_juvenile = mean(total - mature),
            mean_yield = mean(yield),
            rel_amplitude = (max(total) - min(total)) / ((max(total) + min(total)) / 2))
}

cannibalism_fishing_probe <- function(params, fork_sim, fish_level_seq, years = scan_summary_window) {
  params <- seed_from(params, fork_sim)
  bind_rows(lapply(fish_level_seq, function(fl) {
    sim <- project(params, t_max = years, dt = p2$dt, t_save = 0.2,
                   t_start = t_fork, progress_bar = FALSE, effort = fl)
    cbind(fish_level = fl, sim_metrics(sim, t_fork + years / 2))
  }))
}

p_baseline_fixed <- anchovy_params(1, 1, knife_edge_size = 10)
sim_fork_baseline <- make_anchovy_fork_sim(p_baseline_fixed, t_fork)

fish_level_seq_redo <- c(1, 2, 3, 5, 9, 20, 50, 100)

redo_probe_df <- bind_rows(
  cannibalism_fishing_probe(p_baseline_fixed, sim_fork_baseline, fish_level_seq_redo) %>%
    mutate(model = "Baseline (theta=1), FIXED mortality"),
  cannibalism_fishing_probe(p_theta_low_FIXED, sim_fork_theta_low, fish_level_seq_redo) %>%
    mutate(model = "theta=0.3, FIXED mortality")
)
write.csv(redo_probe_df, file.path(plot_dir, "day40_changed_mort_probe.csv"), row.names = FALSE)
cat("Section 1 (Baseline vs. theta=0.3, FIXED mortality, knife_edge=10): see day40_changed_mort_probe.csv.\n")
print(redo_probe_df %>% select(model, fish_level, mean_total, mean_yield, rel_amplitude))

redo_yield_plot <- ggplot(redo_probe_df, aes(x = fish_level, y = mean_yield, color = model)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "Fishing effort (Constant, log scale)", y = "Mean yield",
       title = "Corrected mortality: theta=0.3's yield curve is no longer cleanly single-peaked",
       subtitle = "knife_edge=10 -- local peak near fish_level=9, dip at 20, then rises again by 50-100 (contrast with Day 40's clean peak-then-decline under the buggy mortality)") +
  theme_minimal()
save_plot(redo_yield_plot, "day40_changed_mort_yield.png", width = 9, height = 6)

cat(paste(
  "Section 1 headline finding: the mortality fix barely moves the theta=1 BASELINE",
  "(cannibalism-driven predation mortality dominates there -- fixed-mortality numbers",
  "land almost exactly on Day 39's own buggy-mortality baseline, e.g. mean_yield=0.0985",
  "vs. 0.09846 at fish_level=100). theta=0.3 is a different story: mean_total drops from",
  "the buggy version's own ~0.56-0.65 range to ~0.36-0.43, and -- critically -- the clean",
  "single-peaked yield curve Day 40 built its whole headline result on (peak at",
  "fish_level=2, smooth decline to 100) does NOT survive. The fixed-mortality yield curve",
  "rises to a local peak near fish_level=9, dips at 20, then rises AGAIN by fish_level=50-100,",
  "ending higher than the local peak -- not a clean sustainable-yield shape at all.",
  "This is a first-pass, 8-point-grid finding, not a finished replacement for Day 40's own",
  "fine sweep -- see What's Next.\n"
))

################################################################################
# Overnight: theta x interaction_resource grid under the FIXED mortality --
# mirrors Day 39 Section 3's own grid (same theta_grid_seq/resource_grid_seq),
# to check whether theta=0.3/interaction_resource=1 (Day 40's own choice,
# built under the BUGGY mortality) is still the most fishing-sensitive point,
# or whether the buggy mortality was itself steering that choice toward a
# combination that only looked good by accident.
#
# NOT run interactively. A first attempt at this exact grid ran for over an
# hour in the live session without finishing -- each of the 25 cells does
# three full projections (10yr settle, 10yr kick-regrowth, 12yr fished) at
# dt=0.001 -- and was cancelled rather than left blocking further work.
# Queued here as a batch/overnight job instead (e.g. via Rscript from a
# terminal, not the interactive r-mizer session). Uses the exact resource
# stepper above as-is; whether dt itself can be safely raised to actually
# speed this up hasn't been tested yet -- see the note there.
################################################################################

w_grid_check   <- c(0.01, 1, 10, 30)
idx_grid_check <- vapply(w_grid_check, function(x) which.min(abs(w(p_baseline_fixed) - x)), integer(1))

fork_rate            <- function(getter, params, sim, idx) getter(seed_from(params, sim))[1, idx]
growth_baseline_grid <- fork_rate(getEGrowth, p_baseline_fixed, sim_fork_baseline, idx_grid_check)

grid_probe_effort <- 10

run_interaction_grid_case <- function(theta_val, ir_val) {
  params        <- anchovy_params(theta_val, ir_val, knife_edge_size = p2$w_mat)
  sim_fork_grid <- make_anchovy_fork_sim(params, t_fork)

  bp_grid        <- selected_biomass(sim_fork_grid, t_cut = t_fork - scan_summary_window)
  unfished_ratio <- max(bp_grid) / max(min(bp_grid), 1e-12)

  growth_grid     <- fork_rate(getEGrowth, params, sim_fork_grid, idx_grid_check)
  growth_recovery <- growth_grid / growth_baseline_grid

  total_unfished <- mean(after_cut(getBiomass(sim_fork_grid), sim_fork_grid, t_fork - scan_summary_window))

  sim_fished <- project(seed_from(params, sim_fork_grid), t_max = scan_summary_window, dt = p2$dt,
                        t_save = 0.2, t_start = t_fork, progress_bar = FALSE, effort = grid_probe_effort)
  total_at_effort     <- mean(after_cut(getBiomass(sim_fished), sim_fished, t_fork + scan_summary_window / 2))
  yield_at_effort     <- mean(after_cut(getYield(sim_fished), sim_fished, t_fork + scan_summary_window / 2))
  fishing_sensitivity <- 1 - total_at_effort / total_unfished

  data.frame(theta = theta_val, interaction_resource = ir_val, unfished_cycle_ratio = unfished_ratio,
            growth_recovery_w0.01 = growth_recovery[1], growth_recovery_w1 = growth_recovery[2],
            growth_recovery_w10 = growth_recovery[3], growth_recovery_w30 = growth_recovery[4],
            total_unfished = total_unfished, total_at_effort10 = total_at_effort,
            mean_yield_at_effort10 = yield_at_effort, fishing_sensitivity = fishing_sensitivity)
}

theta_grid_seq    <- c(1, 0.7, 0.5, 0.3, 0.1)
resource_grid_seq <- c(1, 1.3, 2, 3, 5)

interaction_grid_df <- bind_rows(lapply(theta_grid_seq, function(th) {
  bind_rows(lapply(resource_grid_seq, function(ir) run_interaction_grid_case(th, ir)))
}))
write.csv(interaction_grid_df, file.path(plot_dir, "day40_changed_mort_interaction_grid.csv"), row.names = FALSE)
cat("Overnight grid (theta x interaction_resource, FIXED mortality): see day40_changed_mort_interaction_grid.csv.\n")
print(interaction_grid_df)

grid_growth_heatmap <- ggplot(interaction_grid_df,
                              aes(x = factor(theta), y = factor(interaction_resource), fill = growth_recovery_w10)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", growth_recovery_w10)), size = 3, color = "white") +
  scale_fill_viridis_c(limits = c(0, NA)) +
  labs(x = "theta (cannibalism)", y = "interaction_resource", fill = "recovery",
       title = "FIXED mortality: growth recovery at w=10.1g",
       subtitle = "1.0 = fully back to the (theta=1, interaction_resource=1) baseline's own growth there") +
  theme_minimal()
save_plot(grid_growth_heatmap, "day40_changed_mort_grid_growth.png", width = 8, height = 6)

grid_sensitivity_heatmap <- ggplot(interaction_grid_df,
                                   aes(x = factor(theta), y = factor(interaction_resource), fill = fishing_sensitivity)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", fishing_sensitivity)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "magma", limits = c(0, NA)) +
  labs(x = "theta (cannibalism)", y = "interaction_resource", fill = "sensitivity",
       title = "FIXED mortality: fraction of unfished biomass lost at effort=10",
       subtitle = "Higher = more responsive to fishing") +
  theme_minimal()
save_plot(grid_sensitivity_heatmap, "day40_changed_mort_grid_sensitiv.png", width = 8, height = 6)

grid_tradeoff_plot <- ggplot(interaction_grid_df,
                             aes(x = fishing_sensitivity, y = growth_recovery_w10,
                                 color = factor(theta), shape = factor(interaction_resource))) +
  geom_point(size = 3) +
  labs(x = "Fishing sensitivity (fraction of unfished biomass lost at effort=10)",
       y = "Growth recovery at w=10.1g (1.0 = theta=1/ir=1 baseline)",
       color = "theta", shape = "interaction_resource",
       title = "FIXED mortality: the same trade-off Day 39 Section 3 navigated, redone",
       subtitle = "Is theta=0.3/interaction_resource=1 -- Day 40's own choice -- actually the best point here?") +
  theme_minimal()
save_plot(grid_tradeoff_plot, "day40_changed_mort_grid_tradeoff.png", width = 9, height = 6)

baseline_grid_row <- interaction_grid_df %>% filter(theta == 1, interaction_resource == 1)
theta_03_row       <- interaction_grid_df %>% filter(theta == 0.3, interaction_resource == 1)
cat(sprintf(
  "Baseline cell (theta=1, ir=1): fishing_sensitivity=%.3f at effort=10. Day 40's own choice (theta=0.3, ir=1) under FIXED mortality: fishing_sensitivity=%.3f, growth_recovery_w10=%.3f -- compare against the full grid to see whether a different cell would have served the project's aim better.\n",
  baseline_grid_row$fishing_sensitivity, theta_03_row$fishing_sensitivity, theta_03_row$growth_recovery_w10
))

################################################################################
# What's Next (this file's own, not carried from Day 40 -- that file's aim-model
# story is under direct threat from this finding and needs its own follow-up,
# not a patch onto today's):
#
# 1. Finer effort grid between fish_level=9 and 100 (this Section 1 only has
#    9, 20, 50, 100) to pin down whether the dip-then-rise shape is real or an
#    artefact of coarse sampling -- exactly the mistake Section 0c's own period
#    check made and then caught.
# 2. rel_amplitude climbs to 1.3-1.4 by fish_level=50-100 for theta=0.3 here --
#    worth checking (sanity-check plot, Day 36+ convention) whether the
#    population is settling into a genuinely different dynamical regime at
#    high effort under the corrected mortality, not just a noisier version of
#    the same cycle.
# 3. Redo Day 40 Section 5's knife_edge sweep under the fixed mortality -- the
#    collapse-at-small-knife-edge story was this project's other major Day 40
#    finding and hasn't been checked against the corrected model at all yet.
# 4. Once the Overnight grid above has actually run, decide whether
#    40_experiments.R's own results should be superseded or whether the
#    buggy- and fixed-mortality stories need to be reported side by side.
# 5. Now that the resource step is exact rather than explicit-Euler, test
#    whether dt can be raised above 0.001 without changing results -- the
#    consumer side's own numerical scheme (still first-order upwind, no
#    second_order_w/tr_bdf2) may or may not tolerate it. If it does, every
#    sweep in this project gets proportionally cheaper, including the
#    Overnight grid above.
################################################################################
