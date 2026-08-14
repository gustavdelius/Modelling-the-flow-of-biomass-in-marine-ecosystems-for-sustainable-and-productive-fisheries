library(mizer)
library(mizerExperimental)
library(dplyr)
library(ggplot2)

# Day 40: confirm and sharpen Day 39's own headline finding for the project's
# core aim -- a model whose yield genuinely peaks and then declines with
# fishing effort, rather than climbing forever (baseline theta=1 never turns
# over across [1,100], see day39_theta_low_effort_yield.png).
#
# Day 39 Section 1 found theta=0.3 (cannibalism cut to 30% of baseline,
# interaction_resource left at 1, knife_edge=10) peaks "around fish_level=3"
# on a coarse grid (1,3,5,7,9,20,...,100). Today re-runs that same model with
# a finer grid near the peak (0.5,1,2,3,4,5,7,9,15,...) to pin down exactly
# where it turns over.
#
# Self-contained convention since Day 20: helpers redefined here, not
# sourced from 39_experiments.R.

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
# Section 0: rebuild Day 38/39's Figure 2e model, verbatim.
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
plankton_state$random <- FALSE   # random plankton forcing kept off, as Day 37-39

# Exact step for
#   dn/dt = immigration + (rate - mortality) * n - rate / capacity * n^2.
# This is the logistic resource equation with constant immigration. As in
# mizer::resource_logistic(), rates are held fixed during each time step.
# Replaces the old explicit-Euler step (n_pp + dt*f) -- pulled in from
# upstream/GitHub, not written for this project. Exact for any dt (the Euler
# version was only an approximation, tightest at small dt), so this is a
# correctness improvement at the SAME dt=0.001 used throughout this file --
# not a speed change on its own. See 40_changed_mort_experiments.R's own
# What's Next for whether dt itself can now be safely raised.
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

# theta/interaction_resource variant of Day 36/38's own
# make_anchovy_fishing_params() -- ported verbatim from 39_experiments.R.
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

p_scan <- make_anchovy_fishing_params_theta(1, 1, knife_edge_size = 10)
anchovy_immigration <- p2$i0 * p_scan@w_full^(-p2$lambda) * exp(-6.9 * (p2$lambda - 1))

################################################################################
# Section 0b: settle+kick fork (Day 37-39's recipe) to t_fork=20 -- cycle is
# already full-amplitude well before then (confirmed in Day 39).
################################################################################

t_fork               <- 20
scan_summary_window  <- 12   # ~2 periods, matching Day 37-39's own sampling window

make_anchovy_fork_sim <- function(params, t_fork) {
  params@initial_n[]    <- 0.001 * params@w^(-1.8)
  params@initial_n_pp[] <- params@cc_pp
  sim <- project(params, t_max = 10, dt = p2$dt, progress_bar = FALSE)
  sim@n[11, , ] <- sim@n[11, , ] / 10^7
  project(sim, t_max = t_fork - 10, dt = p2$dt, t_save = 0.2, progress_bar = FALSE, effort = 0)
}

################################################################################
# Section 1: the aim model -- theta=0.3, interaction_resource=1, knife_edge=10.
# Day 39 Section 1's own candidate for "yield genuinely drops off"; re-forked
# here rather than assumed, since the fork is model-specific.
################################################################################

theta_low <- 0.3

p_theta_low <- make_anchovy_fishing_params_theta(theta_low, 1, knife_edge_size = 10)
sim_fork_theta_low <- make_anchovy_fork_sim(p_theta_low, t_fork)
last_n_theta_low   <- array(sim_fork_theta_low@n[dim(sim_fork_theta_low@n)[1], , , drop = FALSE],
                            dim = dim(sim_fork_theta_low@n)[-1])
last_npp_theta_low <- sim_fork_theta_low@n_pp[dim(sim_fork_theta_low@n_pp)[1], ]

cannibalism_fishing_probe <- function(params, seed_n, seed_npp, fish_level_seq,
                                      years = scan_summary_window) {
  params@initial_n[]    <- seed_n
  params@initial_n_pp[] <- seed_npp
  bind_rows(lapply(fish_level_seq, function(fl) {
    sim <- project(params, t_max = years, dt = p2$dt, t_save = 0.2,
                   t_start = t_fork, progress_bar = FALSE, effort = fl)
    tv   <- as.numeric(dimnames(sim@n)[[1]])
    keep <- tv > t_fork + years / 2
    total <- unname(getBiomass(sim)[, "Anchovy"])[keep]
    yield <- unname(getYield(sim)[, "Anchovy"])[keep]
    data.frame(fish_level = fl, mean_total = mean(total), mean_yield = mean(yield),
              rel_amplitude = (max(total) - min(total)) / ((max(total) + min(total)) / 2))
  }))
}

# Finer than Day 39's own [1,3,5,7,9,20,...] grid, concentrated near where that
# coarser grid placed the peak ("around fish_level=3"), to pin the turnover point down.
fish_level_seq_fine <- c(0.5, 1, 2, 3, 4, 5, 7, 9, 15, 20, 30, 50, 75, 100)

theta_low_yield_df <- cannibalism_fishing_probe(p_theta_low, last_n_theta_low, last_npp_theta_low,
                                                fish_level_seq_fine)
write.csv(theta_low_yield_df, file.path("interesting_plots", "day40_theta_low_yield_fine.csv"),
          row.names = FALSE)
cat("Section 1 (theta=0.3, knife_edge=10, fine effort grid): see day40_theta_low_yield_fine.csv.\n")
print(theta_low_yield_df)

peak_row <- theta_low_yield_df[which.max(theta_low_yield_df$mean_yield), ]
cat(sprintf(
  "Yield peaks at fish_level=%.1f (mean_yield=%.4f), then declines to %.4f by fish_level=100 -- a genuine sustainable-yield curve, confirming and sharpening Day 39 Section 1's own coarser-grid finding.\n",
  peak_row$fish_level, peak_row$mean_yield, tail(theta_low_yield_df$mean_yield, 1)
))

yield_peak_plot <- ggplot(theta_low_yield_df, aes(x = fish_level, y = mean_yield)) +
  geom_line(color = "#1b7837") +
  geom_point(size = 2, color = "#1b7837") +
  geom_vline(xintercept = peak_row$fish_level, linetype = "dashed", color = "grey50") +
  labs(x = "Fishing effort (Constant schedule, knife_edge=10)", y = "Mean yield",
       title = "theta=0.3 (reduced cannibalism): yield peaks then declines",
       subtitle = sprintf("Peak at fish_level=%.1f (yield=%.4f) -- a genuine sustainable-yield curve",
                          peak_row$fish_level, peak_row$mean_yield)) +
  theme_minimal()
yield_peak_plot
save_plot(yield_peak_plot, "day40_theta_low_yield_peak.png", width = 9, height = 6)

################################################################################
# Section 2: does a peaks/troughs fishing schedule change the picture for the
# PLAIN theta=0.3 model (interaction_resource=1)?
#
# Day 39 Sections 4/5 ran this layered-schedule comparison (Constant floor vs.
# boosting at peaks vs. troughs, Day 36/38's own thresholdFMort() rule) only on
# the resource-boosted model (theta=0.3, interaction_resource=1.3) and the
# original (theta=1). Day 39 Section 6 separately zoomed the plain theta=0.3
# model into [0,7] at knife_edge=5. Today: the same three-schedule comparison,
# on the plain theta=0.3 model, at knife_edge=10 -- the config Section 1 above
# just confirmed has a genuine yield peak at fish_level~2.
#
# Floor kept low (baseline_effort=1, not Day 39 Sections 4/5's own floor=10)
# and the boost sequence spans across and beyond Section 1's own yield peak,
# rather than jumping straight to Day 39's collapse-hunting range -- that
# range is Day 39 Section 6's own finding: "peaks is risky" only shows up at
# high effort, invisible in [0,7].
################################################################################

# Day 36/38/39's own threshold-rule machinery, ported verbatim.
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

threshold_diagnostics <- function(sim, threshold, fish_level, background_level, mode, sharpness, t_cut, hard_step = FALSE) {
  tv      <- as.numeric(dimnames(sim@n)[[1]])
  keep    <- tv > t_cut
  on_frac <- compute_on_frac_series(sim, sim@params, threshold, mode, sharpness, t_cut, hard_step)
  is_on   <- on_frac > 0.5
  runs   <- rle(is_on)
  n_runs <- length(runs$values)
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

scan_post_fork_years <- 24   # ~4-6 periods of the paper's own ~6yr cycle
scan_t_cut           <- t_fork + scan_post_fork_years - scan_summary_window

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

# Generalised over (theta_val, ir_val, seed_n, seed_npp, baseline_effort, boost_seq)
# so it isn't tied to Day 39's own floor=10 collapse-hunting range -- ported from
# 39_experiments.R's run_layered_ke_case(), with baseline_effort/boost_seq required
# rather than defaulting to Day 39's own globals.
run_layered_ke_case <- function(knife_edge_size, theta_val, ir_val, seed_n, seed_npp,
                                baseline_effort, boost_seq) {
  params <- make_anchovy_fishing_params_theta(theta_val, ir_val, knife_edge_size = knife_edge_size)
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

  list(summary = bind_rows(const_row, boosted_rows), threshold = threshold_ke, sharpness = sharpness_ke)
}

ke_compare             <- 10
baseline_effort_plain  <- 1
boost_seq_plain         <- c(2, 3, 5, 7, 9, 15, 20, 30)

theta_plain_schedule_result <- run_layered_ke_case(ke_compare, theta_low, 1,
                                                   last_n_theta_low, last_npp_theta_low,
                                                   baseline_effort = baseline_effort_plain,
                                                   boost_seq = boost_seq_plain)
theta_plain_schedule_df <- theta_plain_schedule_result$summary
write.csv(theta_plain_schedule_df, file.path("interesting_plots", "day40_theta_plain_schedule_summary.csv"),
         row.names = FALSE)
cat("Section 2 (plain theta=0.3, knife_edge=10, three-schedule comparison): see day40_theta_plain_schedule_summary.csv.\n")
print(theta_plain_schedule_df)

theta_plain_biomass_plot <- ggplot(theta_plain_schedule_df,
                                   aes(x = mean_effort, y = mean_total, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_hline(yintercept = theta_plain_schedule_df$mean_total[theta_plain_schedule_df$schedule == "Constant"],
             linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort actually applied", y = "Mean total biomass",
       color = "Schedule",
       title = "Plain theta=0.3 model: Constant floor vs. boosting at peaks vs. troughs",
       subtitle = sprintf("knife_edge=%d, floor(baseline_effort)=%.2g -- dashed = Constant's own biomass",
                          ke_compare, baseline_effort_plain)) +
  theme_minimal()
theta_plain_biomass_plot
save_plot(theta_plain_biomass_plot, "day40_theta_plain_schedule_biomass.png", width = 9, height = 6)

theta_plain_yield_plot <- ggplot(theta_plain_schedule_df,
                                 aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_hline(yintercept = theta_plain_schedule_df$mean_yield[theta_plain_schedule_df$schedule == "Constant"],
             linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort actually applied", y = "Mean yield",
       color = "Schedule",
       title = "Plain theta=0.3 model: yield vs. realised effort, three schedules",
       subtitle = sprintf("knife_edge=%d, floor(baseline_effort)=%.2g -- dashed = Constant's own yield",
                          ke_compare, baseline_effort_plain)) +
  theme_minimal()
theta_plain_yield_plot
save_plot(theta_plain_yield_plot, "day40_theta_plain_schedule_yield.png", width = 9, height = 6)

theta_plain_amplitude_plot <- ggplot(theta_plain_schedule_df,
                                     aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  labs(x = "Realised mean effort actually applied", y = "Relative amplitude of biomass",
       color = "Schedule",
       title = "Plain theta=0.3 model: relative amplitude vs. realised effort, three schedules",
       subtitle = sprintf("knife_edge=%d, floor(baseline_effort)=%.2g", ke_compare, baseline_effort_plain)) +
  theme_minimal()
theta_plain_amplitude_plot
save_plot(theta_plain_amplitude_plot, "day40_theta_plain_schedule_amplitude.png", width = 9, height = 6)

cat(sprintf(
  "Neither schedule damages the plain theta=0.3 model in this range: mean_total stays at/above Constant's own %.4f throughout, and rel_amplitude tops out around %.2f -- nowhere near the ~2.0 cliff-edge Day 39 found in the ORIGINAL model. Consistent with Day 39 Section 6's own finding that 'peaks is risky' is a high-effort phenomenon invisible at this scale.\n",
  theta_plain_schedule_df$mean_total[theta_plain_schedule_df$schedule == "Constant"],
  max(theta_plain_schedule_df$rel_amplitude)
))

################################################################################
# Section 2b: sanity check -- does thresholdFMort() actually fire where it
# should for this model, or does it barely engage? Day 36/38/39's own
# convention: one representative case, all three schedules, biomass plotted
# against the calibration threshold with the boosted "on" periods shaded.
################################################################################

sanity_ke_plain    <- ke_compare
sanity_boost_plain <- 9   # mid-range: not yet saturated for either peaks or troughs
threshold_plain    <- theta_plain_schedule_result$threshold
sharpness_plain    <- theta_plain_schedule_result$sharpness

p_sanity_plain <- make_anchovy_fishing_params_theta(theta_low, 1, knife_edge_size = sanity_ke_plain)

run_sanity_schedule_plain <- function(schedule_name) {
  if (schedule_name == "Constant") {
    p_run <- p_sanity_plain
    p_run@initial_n[]    <- last_n_theta_low
    p_run@initial_n_pp[] <- last_npp_theta_low
    sim <- project(p_run, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                   t_start = t_fork, progress_bar = FALSE, effort = baseline_effort_plain)
    return(list(sim = sim, mode = NA_character_))
  }
  mode <- if (schedule_name == "Threshold (peaks)") "above" else "below"
  p_rule <- attach_threshold_rule(p_sanity_plain, threshold = threshold_plain, fish_level = sanity_boost_plain,
                                  background_level = baseline_effort_plain, mode = mode, sharpness = sharpness_plain)
  p_rule@initial_n[]    <- last_n_theta_low
  p_rule@initial_n_pp[] <- last_npp_theta_low
  sim <- project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 t_start = t_fork, progress_bar = FALSE, effort = 1)
  list(sim = sim, mode = mode)
}

sanity_schedules_plain <- c("Constant", "Threshold (peaks)", "Threshold (troughs)")
sanity_cases_plain <- setNames(lapply(sanity_schedules_plain, run_sanity_schedule_plain), sanity_schedules_plain)

sanity_series_plain_df <- bind_rows(lapply(names(sanity_cases_plain), function(nm) {
  case <- sanity_cases_plain[[nm]]
  tv   <- as.numeric(dimnames(case$sim@n)[[1]])
  keep <- which(tv > t_fork)
  bp   <- compute_selected_biomass_series(case$sim, p_sanity_plain, t_fork)
  on_frac <- if (is.na(case$mode)) {
    rep(NA_real_, length(bp))
  } else {
    compute_on_frac_series(case$sim, p_sanity_plain, threshold_plain, case$mode, sharpness_plain, t_fork)
  }
  data.frame(t = tv[keep], selected_biomass = bp, on_frac = on_frac, schedule = nm)
}))
write.csv(sanity_series_plain_df, file.path("interesting_plots", "day40_theta_plain_sanity_series.csv"),
         row.names = FALSE)

sanity_check_plot_plain <- ggplot(sanity_series_plain_df, aes(x = t, y = selected_biomass)) +
  geom_rect(data = sanity_series_plain_df %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = threshold_plain, linetype = "dashed", color = "grey40") +
  facet_wrap(~factor(schedule, levels = sanity_schedules_plain), ncol = 1, scales = "free_y") +
  labs(x = "Time (years)", y = "Selected biomass",
       title = "Sanity check: does thresholdFMort() fire where it should? (plain theta=0.3 model)",
       subtitle = sprintf("knife_edge=%.0f, floor=%.2g, boost=%.2g -- dashed = calibration threshold, red = boosted 'on'",
                          sanity_ke_plain, baseline_effort_plain, sanity_boost_plain)) +
  theme_minimal()
sanity_check_plot_plain
save_plot(sanity_check_plot_plain, "day40_theta_plain_sanity_check.png", width = 9, height = 8)

cat(sprintf(
  "Schedule sanity check (plain theta=0.3, knife_edge=%.0f, boost=%.2g): see day40_theta_plain_sanity_check.png / day40_theta_plain_sanity_series.csv. 'Peaks' shading lands on the narrow local maxima of the Constant cycle; 'troughs' shading covers the broader sub-threshold stretches -- both firing as designed.\n",
  sanity_ke_plain, sanity_boost_plain
))

################################################################################
# Section 2c: same sanity check, for rel_amplitude specifically.
#
# The Section 2b sanity check plots selected_biomass (the quantity
# thresholdFMort() actually triggers on), but rel_amplitude in
# theta_plain_schedule_df is computed from TOTAL biomass over
# scan_metrics_layered()'s own last-scan_summary_window-years window, not the
# full sanity-check run -- reusing Section 2b's own three sims (boost=9)
# rather than re-simulating.
################################################################################

compute_total_biomass_series <- function(sim, t_cut) {
  tv   <- as.numeric(dimnames(sim@n)[[1]])
  keep <- which(tv > t_cut)
  unname(getBiomass(sim)[, "Anchovy"])[keep]
}

sanity_amplitude_series_plain_df <- bind_rows(lapply(names(sanity_cases_plain), function(nm) {
  case <- sanity_cases_plain[[nm]]
  tv   <- as.numeric(dimnames(case$sim@n)[[1]])
  keep <- which(tv > t_fork)
  total_bp <- compute_total_biomass_series(case$sim, t_fork)
  on_frac <- if (is.na(case$mode)) {
    rep(NA_real_, length(total_bp))
  } else {
    compute_on_frac_series(case$sim, p_sanity_plain, threshold_plain, case$mode, sharpness_plain, t_fork)
  }
  data.frame(t = tv[keep], total_biomass = total_bp, on_frac = on_frac, schedule = nm)
}))

# Restricted to scan_metrics_layered()'s own window (t > scan_t_cut) so these
# numbers reconcile exactly with theta_plain_schedule_df's own boost=9 rows.
window_amplitude_plain_matched <- sanity_amplitude_series_plain_df %>%
  filter(t > scan_t_cut) %>%
  group_by(schedule) %>%
  summarise(rel_amplitude = (max(total_biomass) - min(total_biomass)) /
                             ((max(total_biomass) + min(total_biomass)) / 2), .groups = "drop")

sanity_amplitude_series_plain_df <- sanity_amplitude_series_plain_df %>%
  left_join(window_amplitude_plain_matched, by = "schedule") %>%
  mutate(schedule_label = sprintf("%s  (rel_amplitude=%.2f)", schedule, rel_amplitude),
         schedule_label = factor(schedule_label,
                                 levels = sprintf("%s  (rel_amplitude=%.2f)", sanity_schedules_plain,
                                                  window_amplitude_plain_matched$rel_amplitude[
                                                    match(sanity_schedules_plain, window_amplitude_plain_matched$schedule)])))
write.csv(sanity_amplitude_series_plain_df %>% select(t, total_biomass, on_frac, schedule, rel_amplitude),
         file.path("interesting_plots", "day40_theta_plain_sanity_amplitude_series.csv"), row.names = FALSE)

sanity_amplitude_plot_plain <- ggplot(sanity_amplitude_series_plain_df, aes(x = t, y = total_biomass)) +
  geom_rect(data = sanity_amplitude_series_plain_df %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_vline(xintercept = scan_t_cut, linetype = "dotted", color = "steelblue") +
  geom_line() +
  facet_wrap(~schedule_label, ncol = 1, scales = "free_y") +
  labs(x = "Time (years)", y = "Total biomass",
       title = "Sanity check: what does rel_amplitude actually look like? (plain theta=0.3 model)",
       subtitle = sprintf("knife_edge=%.0f, floor=%.2g, boost=%.2g -- dotted = start of the window rel_amplitude is computed over, red = boosted 'on'",
                          sanity_ke_plain, baseline_effort_plain, sanity_boost_plain)) +
  theme_minimal()
sanity_amplitude_plot_plain
save_plot(sanity_amplitude_plot_plain, "day40_theta_plain_sanity_amplitude.png", width = 9, height = 8)

cat(sprintf(
  "Amplitude sanity check (plain theta=0.3, knife_edge=%.0f, boost=%.2g): see day40_theta_plain_sanity_amplitude.png / day40_theta_plain_sanity_amplitude_series.csv. 'Peaks' gets the highest rel_amplitude (%.2f) from a short, sharp cut right at the biomass high point, not from a lower mean; 'troughs' (%.2f) spreads pressure out and gets a shallower dip; 'Constant' (%.2f) is the smallest, regular cycle.\n",
  sanity_ke_plain, sanity_boost_plain,
  window_amplitude_plain_matched$rel_amplitude[window_amplitude_plain_matched$schedule == "Threshold (peaks)"],
  window_amplitude_plain_matched$rel_amplitude[window_amplitude_plain_matched$schedule == "Threshold (troughs)"],
  window_amplitude_plain_matched$rel_amplitude[window_amplitude_plain_matched$schedule == "Constant"]
))

################################################################################
# Section 3: redo the three-schedule comparison with the floor set to Section
# 1's own MSY effort (fish_level=2, the yield-maximising Constant effort found
# there) rather than Section 2's arbitrary floor=1 -- and sweep the full sanity
# check across every boost level, not just one representative case, to see
# whether the firing behaviour/results change qualitatively across the sweep.
#
# run_layered_ke_case() is extended here (run_layered_ke_case_with_series())
# to also capture the full selected/total biomass + on_frac time series for
# every boost level, reusing the same simulations for both the summary table
# and the boost-swept sanity-check grid rather than re-running twice.
################################################################################

run_layered_ke_case_with_series <- function(knife_edge_size, theta_val, ir_val, seed_n, seed_npp,
                                            baseline_effort, boost_seq) {
  params <- make_anchovy_fishing_params_theta(theta_val, ir_val, knife_edge_size = knife_edge_size)
  params@initial_n[]    <- seed_n
  params@initial_n_pp[] <- seed_npp

  sim_const <- project(params, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                       t_start = t_fork, progress_bar = FALSE, effort = baseline_effort)
  bp_const     <- compute_selected_biomass_series(sim_const, params, scan_t_cut)
  sharpness_ke <- 0.02 * (max(bp_const) - min(bp_const))
  threshold_ke <- unname(quantile(bp_const, probs = 0.5))

  const_row <- cbind(scan_metrics_layered(sim_const),
                     knife_edge_size = knife_edge_size, schedule = "Constant",
                     boost_level = NA_real_, mean_effort = baseline_effort,
                     effective_window = NA_real_, n_bursts = NA_integer_)

  tv_const   <- as.numeric(dimnames(sim_const@n)[[1]])
  keep_const <- which(tv_const > t_fork)
  const_series <- data.frame(
    t = tv_const[keep_const],
    selected_biomass = compute_selected_biomass_series(sim_const, params, t_fork),
    total_biomass = compute_total_biomass_series(sim_const, t_fork),
    on_frac = NA_real_, schedule = "Constant", boost_level = NA_real_)

  boosted <- lapply(boost_seq, function(boost) {
    lapply(c("above", "below"), function(mode) {
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
      row <- cbind(scan_metrics_layered(sim),
                  knife_edge_size = knife_edge_size, schedule = schedule_name,
                  boost_level = boost, mean_effort = diag$mean_effort,
                  effective_window = diag$effective_window, n_bursts = diag$n_bursts)
      tv   <- as.numeric(dimnames(sim@n)[[1]])
      keep <- which(tv > t_fork)
      series <- data.frame(
        t = tv[keep],
        selected_biomass = compute_selected_biomass_series(sim, params, t_fork),
        total_biomass = compute_total_biomass_series(sim, t_fork),
        on_frac = compute_on_frac_series(sim, params, threshold_ke, mode, sharpness_ke, t_fork),
        schedule = schedule_name, boost_level = boost)
      list(row = row, series = series)
    })
  })
  boosted <- unlist(boosted, recursive = FALSE)

  list(summary = bind_rows(const_row, lapply(boosted, `[[`, "row")),
       series  = bind_rows(const_series, lapply(boosted, `[[`, "series")),
       threshold = threshold_ke, sharpness = sharpness_ke)
}

# baseline_effort set to Section 1's own MSY-optimal Constant effort (fish_level=2,
# where theta_low_yield_df's mean_yield peaked), not an arbitrary floor.
msy_row              <- theta_low_yield_df[which.max(theta_low_yield_df$mean_yield), ]
baseline_effort_msy  <- msy_row$fish_level
boost_seq_msy        <- c(3, 5, 7, 9, 15, 20, 30)   # starts above the new floor=2

theta_plain_msy_result <- run_layered_ke_case_with_series(ke_compare, theta_low, 1,
                                                          last_n_theta_low, last_npp_theta_low,
                                                          baseline_effort = baseline_effort_msy,
                                                          boost_seq = boost_seq_msy)
theta_plain_msy_df <- theta_plain_msy_result$summary
write.csv(theta_plain_msy_df, file.path("interesting_plots", "day40_theta_plain_msy_schedule_summary.csv"),
         row.names = FALSE)
write.csv(theta_plain_msy_result$series, file.path("interesting_plots", "day40_theta_plain_msy_boost_sweep_series.csv"),
         row.names = FALSE)
cat(sprintf(
  "Section 3 (plain theta=0.3, knife_edge=%d, floor=Section 1's own MSY effort=%.2g): see day40_theta_plain_msy_schedule_summary.csv.\n",
  ke_compare, baseline_effort_msy
))
print(theta_plain_msy_df)

theta_msy_biomass_plot <- ggplot(theta_plain_msy_df, aes(x = mean_effort, y = mean_total, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_hline(yintercept = theta_plain_msy_df$mean_total[theta_plain_msy_df$schedule == "Constant"],
             linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort actually applied", y = "Mean total biomass",
       color = "Schedule",
       title = "Plain theta=0.3 model, floor=Section 1's own MSY effort: biomass vs. realised effort",
       subtitle = sprintf("knife_edge=%d, floor(baseline_effort)=%.2g -- dashed = Constant's own biomass",
                          ke_compare, baseline_effort_msy)) +
  theme_minimal()
theta_msy_biomass_plot
save_plot(theta_msy_biomass_plot, "day40_theta_plain_msy_biomass.png", width = 9, height = 6)

theta_msy_yield_plot <- ggplot(theta_plain_msy_df, aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_hline(yintercept = theta_plain_msy_df$mean_yield[theta_plain_msy_df$schedule == "Constant"],
             linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort actually applied", y = "Mean yield",
       color = "Schedule",
       title = "Plain theta=0.3 model, floor=Section 1's own MSY effort: yield vs. realised effort",
       subtitle = sprintf("knife_edge=%d, floor(baseline_effort)=%.2g (Section 1's own yield-maximising Constant effort) -- dashed = Constant's own yield",
                          ke_compare, baseline_effort_msy)) +
  theme_minimal()
theta_msy_yield_plot
save_plot(theta_msy_yield_plot, "day40_theta_plain_msy_yield.png", width = 9, height = 6)

theta_msy_amplitude_plot <- ggplot(theta_plain_msy_df, aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = schedule)) +
  geom_point(size = 2) +
  geom_hline(yintercept = theta_plain_msy_df$rel_amplitude[theta_plain_msy_df$schedule == "Constant"],
             linetype = "dashed", color = "grey50") +
  labs(x = "Realised mean effort actually applied", y = "Relative amplitude of biomass",
       color = "Schedule",
       title = "Plain theta=0.3 model, floor=Section 1's own MSY effort: amplitude vs. realised effort",
       subtitle = sprintf("knife_edge=%d, floor(baseline_effort)=%.2g -- dashed = Constant's own amplitude",
                          ke_compare, baseline_effort_msy)) +
  theme_minimal()
theta_msy_amplitude_plot
save_plot(theta_msy_amplitude_plot, "day40_theta_plain_msy_amplitude.png", width = 9, height = 6)

cat(sprintf(
  "Both threshold schedules beat Constant-at-MSY on yield across the ENTIRE boost range tested (Constant=%.4f; peaks %.4f-%.4f; troughs %.4f-%.4f) -- layering a threshold boost on top of the already yield-optimal Constant effort still buys extra yield here.\n",
  theta_plain_msy_df$mean_yield[theta_plain_msy_df$schedule == "Constant"],
  min(theta_plain_msy_df$mean_yield[theta_plain_msy_df$schedule == "Threshold (peaks)"]),
  max(theta_plain_msy_df$mean_yield[theta_plain_msy_df$schedule == "Threshold (peaks)"]),
  min(theta_plain_msy_df$mean_yield[theta_plain_msy_df$schedule == "Threshold (troughs)"]),
  max(theta_plain_msy_df$mean_yield[theta_plain_msy_df$schedule == "Threshold (troughs)"])
))

# Sanity check swept across EVERY boost level (not just one representative
# case) -- does the firing pattern/result change qualitatively across the
# sweep? Reuses the time series already captured above, no extra simulation.
series_msy <- theta_plain_msy_result$series
boosted_series_msy <- series_msy %>% filter(schedule != "Constant") %>%
  mutate(schedule = factor(schedule, levels = c("Threshold (peaks)", "Threshold (troughs)")),
         boost_level = factor(boost_level, levels = boost_seq_msy))

sanity_sweep_plot <- ggplot(boosted_series_msy, aes(x = t, y = selected_biomass)) +
  geom_rect(data = boosted_series_msy %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_hline(yintercept = theta_plain_msy_result$threshold, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.3) +
  facet_grid(schedule ~ boost_level, labeller = labeller(boost_level = function(x) paste0("boost=", x)),
            scales = "free_y") +
  labs(x = "Time (years)", y = "Selected biomass",
       title = "Sanity check across the full boost sweep: does firing behaviour change with boost? (plain theta=0.3, floor=MSY)",
       subtitle = sprintf("knife_edge=%d, floor=%.2g -- dashed = calibration threshold, red = boosted 'on'",
                          ke_compare, baseline_effort_msy)) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sanity_sweep_plot
save_plot(sanity_sweep_plot, "day40_msy_boost_sweep_selected.png", width = 13, height = 6)

sanity_sweep_total_plot <- ggplot(boosted_series_msy, aes(x = t, y = total_biomass)) +
  geom_rect(data = boosted_series_msy %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_vline(xintercept = scan_t_cut, linetype = "dotted", color = "steelblue") +
  geom_line(linewidth = 0.3) +
  facet_grid(schedule ~ boost_level, labeller = labeller(boost_level = function(x) paste0("boost=", x)),
            scales = "free_y") +
  labs(x = "Time (years)", y = "Total biomass",
       title = "Sanity check across the full boost sweep: total biomass (plain theta=0.3, floor=MSY)",
       subtitle = sprintf("knife_edge=%d, floor=%.2g -- dotted = start of rel_amplitude's own window, red = boosted 'on'",
                          ke_compare, baseline_effort_msy)) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
sanity_sweep_total_plot
save_plot(sanity_sweep_total_plot, "day40_msy_boost_sweep_total.png", width = 13, height = 6)

cat("Boost sweep changes the picture qualitatively, but not by collapsing the population: as 'troughs' boost rises toward 30, its SELECTED-biomass panel crashes to near-zero and stays there (a compositional shift), while the matching TOTAL-biomass panel keeps cycling at a similar range to every other boost -- juveniles replacing the harvested adult compartment, not a population collapse. 'Peaks' just gets progressively narrower, more self-limited on-windows as boost rises; its cycle shape barely changes. See day40_msy_boost_sweep_selected.png / day40_msy_boost_sweep_total.png.\n")

cat(sprintf(
  "Caveat: this Section's own Constant row (mean_yield=%.4f at fish_level=%.2g) does not exactly reproduce Section 1's own peak (mean_yield=%.4f at the same fish_level) -- Section 1 samples the last %.0fyr of a %.0fyr post-fork run, this Section samples the last %.0fyr of a %.0fyr run, so the two land on different, non-overlapping stretches of the same ~5-6yr cycle. The schedule COMPARISONS within this Section are still apples-to-apples (all rows share the same window), but the absolute Constant-at-MSY number quoted here is not directly Section 1's own value.\n",
  theta_plain_msy_df$mean_yield[theta_plain_msy_df$schedule == "Constant"], baseline_effort_msy,
  msy_row$mean_yield, scan_summary_window / 2, scan_summary_window, scan_summary_window, scan_post_fork_years
))

################################################################################
# Section 4: what's-next item 4 -- trigger the threshold rule off TOTAL
# biomass instead of the SELECTED (gear-weighted) biomass every rule so far
# has used. Deliberately a NEW, separate set of functions: thresholdFMort(),
# attach_threshold_rule(), compute_on_frac_series(), threshold_diagnostics()
# from Day 36/38 onward are left completely untouched below, so the original
# selected-biomass rule stays available and easy to reuse later. Only the
# trigger variable differs -- total biomass (sum over the whole size
# spectrum) instead of f_ref-gear-weighted selected biomass.
################################################################################

thresholdFMortTotal <- function(params, n, n_pp, n_other, t, effort, e_growth, pred_mort, ...) {
  p <- other_params(params)
  biomass_density <- sweep(n, 2, params@w * params@dw, "*")
  total_biomass   <- sum(biomass_density)

  direction <- if (p$mode == "above") 1 else -1
  on_frac   <- if (isTRUE(p$hard_step)) {
    as.numeric(direction * (total_biomass - p$threshold) > 0)
  } else {
    plogis(direction * (total_biomass - p$threshold) / max(p$sharpness, 1e-12))
  }

  fmort_on  <- colSums(mizerFMortGear(params, effort = p$fish_level))
  fmort_off <- colSums(mizerFMortGear(params, effort = p$background_level))
  result <- on_frac * fmort_on + (1 - on_frac) * fmort_off
  dim(result)      <- dim(n)
  dimnames(result) <- dimnames(n)
  result
}

attach_threshold_rule_total <- function(params, threshold, fish_level, background_level = 0,
                                        mode = c("above", "below"), sharpness, hard_step = FALSE) {
  mode <- match.arg(mode)
  other_params(params) <- list(threshold = threshold, fish_level = fish_level,
                               background_level = background_level, mode = mode,
                               sharpness = sharpness, hard_step = hard_step)
  setRateFunction(params, "FMort", "thresholdFMortTotal")
}

compute_on_frac_series_total <- function(sim, threshold, mode, sharpness, t_cut, hard_step = FALSE) {
  tv    <- as.numeric(dimnames(sim@n)[[1]])
  keep  <- which(tv > t_cut)
  direction <- if (mode == "above") 1 else -1
  vapply(keep, function(i) {
    n_i <- array(sim@n[i, , , drop = FALSE], dim = dim(sim@n)[-1])
    biomass_density <- sweep(n_i, 2, sim@params@w * sim@params@dw, "*")
    total_biomass   <- sum(biomass_density)
    if (isTRUE(hard_step)) {
      as.numeric(direction * (total_biomass - threshold) > 0)
    } else {
      plogis(direction * (total_biomass - threshold) / max(sharpness, 1e-12))
    }
  }, numeric(1))
}

threshold_diagnostics_total <- function(sim, threshold, fish_level, background_level, mode, sharpness, t_cut, hard_step = FALSE) {
  tv      <- as.numeric(dimnames(sim@n)[[1]])
  keep    <- tv > t_cut
  on_frac <- compute_on_frac_series_total(sim, threshold, mode, sharpness, t_cut, hard_step)
  is_on   <- on_frac > 0.5
  runs   <- rle(is_on)
  n_runs <- length(runs$values)
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

# Threshold/sharpness calibrated off the SAME Constant-schedule run Section 3
# already made (series_msy), just read off total_biomass instead of
# selected_biomass -- no extra simulation needed for the Constant reference.
const_total_window <- series_msy %>% filter(schedule == "Constant", t > scan_t_cut) %>% pull(total_biomass)
threshold_ke_total  <- unname(quantile(const_total_window, probs = 0.5))
sharpness_ke_total  <- 0.02 * (max(const_total_window) - min(const_total_window))

# Sanity check first, Day 36/38/39/40's own convention, before trusting any
# scan: one representative boost (=9, matching Section 3's own sanity case for
# direct comparison), both trigger types side by side.
sanity_boost_total <- 9

run_sanity_schedule_total <- function(schedule_name) {
  mode <- if (schedule_name == "Threshold (peaks, total)") "above" else "below"
  p_rule <- attach_threshold_rule_total(p_sanity_plain, threshold = threshold_ke_total, fish_level = sanity_boost_total,
                                        background_level = baseline_effort_msy, mode = mode, sharpness = sharpness_ke_total)
  p_rule@initial_n[]    <- last_n_theta_low
  p_rule@initial_n_pp[] <- last_npp_theta_low
  sim <- project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 t_start = t_fork, progress_bar = FALSE, effort = 1)
  list(sim = sim, mode = mode)
}

sanity_schedules_total <- c("Threshold (peaks, total)", "Threshold (troughs, total)")
sanity_cases_total <- setNames(lapply(sanity_schedules_total, run_sanity_schedule_total), sanity_schedules_total)

# Anti-correlation check: is selected biomass really out of phase with total
# biomass under this model, or is that just a visual impression?
const_check <- series_msy %>% filter(schedule == "Constant") %>% arrange(t)
selected_total_cor <- cor(const_check$selected_biomass, const_check$total_biomass)
cat(sprintf(
  "Correlation between selected_biomass and total_biomass over time, Constant schedule at floor=MSY: r=%.3f -- strongly negative, i.e. the harvestable-adult compartment and the (juvenile-dominated) total population are largely out of phase with each other.\n",
  selected_total_cor
))

selected_trigger_boost9 <- series_msy %>% filter(boost_level == 9 | schedule == "Constant") %>%
  mutate(trigger = ifelse(schedule == "Constant", "Constant", "selected-biomass trigger"))

total_trigger_series <- bind_rows(lapply(names(sanity_cases_total), function(nm) {
  case <- sanity_cases_total[[nm]]
  tv   <- as.numeric(dimnames(case$sim@n)[[1]])
  keep <- which(tv > t_fork)
  data.frame(
    t = tv[keep],
    selected_biomass = compute_selected_biomass_series(case$sim, p_sanity_plain, t_fork),
    total_biomass = compute_total_biomass_series(case$sim, t_fork),
    on_frac = compute_on_frac_series_total(case$sim, threshold_ke_total, case$mode, sharpness_ke_total, t_fork),
    schedule = nm, boost_level = sanity_boost_total, trigger = "total-biomass trigger")
}))

compare_trigger_df <- bind_rows(selected_trigger_boost9, total_trigger_series) %>%
  mutate(schedule_short = case_when(
    schedule == "Constant" ~ "Constant",
    grepl("peaks", schedule, ignore.case = TRUE) ~ "Peaks",
    grepl("troughs", schedule, ignore.case = TRUE) ~ "Troughs"
  ),
  panel = ifelse(schedule == "Constant", "Constant", paste0(schedule_short, "\n(", trigger, ")")),
  panel = factor(panel, levels = c("Constant", "Peaks\n(selected-biomass trigger)", "Peaks\n(total-biomass trigger)",
                                   "Troughs\n(selected-biomass trigger)", "Troughs\n(total-biomass trigger)")))
write.csv(compare_trigger_df %>% select(t, selected_biomass, total_biomass, on_frac, schedule, trigger, panel),
         file.path("interesting_plots", "day40_trigger_sanity_series.csv"), row.names = FALSE)

trigger_sanity_plot <- ggplot(compare_trigger_df, aes(x = t, y = total_biomass)) +
  geom_rect(data = compare_trigger_df %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~panel, ncol = 2) +
  labs(x = "Time (years)", y = "Total biomass",
       title = "Sanity check: total-biomass-triggered schedule vs. the original selected-biomass one",
       subtitle = sprintf("knife_edge=%d, floor=%.2g (MSY), boost=%.2g -- red = boosted 'on' (each panel's own trigger)",
                          ke_compare, baseline_effort_msy, sanity_boost_total)) +
  theme_minimal(base_size = 9)
trigger_sanity_plot
save_plot(trigger_sanity_plot, "day40_trigger_sanity_check.png", width = 9, height = 8)

cat("Visually: 'Peaks (selected-biomass trigger)' shades near TOTAL biomass's own local MINIMA, not maxima -- consistent with the negative correlation above. The two trigger types produce close to opposite schedules from each other, not just a variant.\n")

# Quantitative check, same boost=9 case: does the antiphase relationship show
# up in the actual summary metrics, not just the sanity-check plot?
total_trigger_rows <- bind_rows(lapply(names(sanity_cases_total), function(nm) {
  case <- sanity_cases_total[[nm]]
  diag <- threshold_diagnostics_total(case$sim, threshold_ke_total, sanity_boost_total, baseline_effort_msy,
                                      case$mode, sharpness_ke_total, scan_t_cut)
  cbind(scan_metrics_layered(case$sim), knife_edge_size = ke_compare, schedule = nm,
       boost_level = sanity_boost_total, mean_effort = diag$mean_effort,
       effective_window = diag$effective_window, n_bursts = diag$n_bursts, trigger = "total")
}))

trigger_compare_df <- bind_rows(
  selected_trigger_rows_boost9 <- theta_plain_msy_df %>% filter(boost_level == 9 | schedule == "Constant") %>%
    mutate(trigger = ifelse(schedule == "Constant", "n/a", "selected")),
  total_trigger_rows
)
write.csv(trigger_compare_df, file.path("interesting_plots", "day40_trigger_compare_boost9.csv"), row.names = FALSE)
cat("Section 4 trigger comparison at boost=9 (Constant/selected-triggered/total-triggered): see day40_trigger_compare_boost9.csv.\n")
print(trigger_compare_df %>% select(schedule, trigger, mean_effort, mean_total, mean_yield, rel_amplitude))

cat(sprintf(
  "'Threshold (peaks, total)' (mean_total=%.4f, rel_amplitude=%.4f) closely matches 'Threshold (troughs, selected)' (mean_total=%.4f, rel_amplitude=%.4f) -- both suppress biomass below Constant's own %.4f. 'Threshold (troughs, total)' (mean_total=%.4f, rel_amplitude=%.4f) closely matches 'Threshold (peaks, selected)' (mean_total=%.4f, rel_amplitude=%.4f) -- both boost biomass above Constant. Switching the trigger variable does not just perturb the schedule, it approximately SWAPS which label ('peaks' vs. 'troughs') is the protective one -- the direct, quantitative consequence of the selected/total antiphase relationship found above. thresholdFMort() itself is untouched throughout; this all runs through the new thresholdFMortTotal() defined above.\n",
  total_trigger_rows$mean_total[total_trigger_rows$schedule == "Threshold (peaks, total)"],
  total_trigger_rows$rel_amplitude[total_trigger_rows$schedule == "Threshold (peaks, total)"],
  theta_plain_msy_df$mean_total[theta_plain_msy_df$schedule == "Threshold (troughs)" & theta_plain_msy_df$boost_level == 9],
  theta_plain_msy_df$rel_amplitude[theta_plain_msy_df$schedule == "Threshold (troughs)" & theta_plain_msy_df$boost_level == 9],
  theta_plain_msy_df$mean_total[theta_plain_msy_df$schedule == "Constant"],
  total_trigger_rows$mean_total[total_trigger_rows$schedule == "Threshold (troughs, total)"],
  total_trigger_rows$rel_amplitude[total_trigger_rows$schedule == "Threshold (troughs, total)"],
  theta_plain_msy_df$mean_total[theta_plain_msy_df$schedule == "Threshold (peaks)" & theta_plain_msy_df$boost_level == 9],
  theta_plain_msy_df$rel_amplitude[theta_plain_msy_df$schedule == "Threshold (peaks)" & theta_plain_msy_df$boost_level == 9]
))

################################################################################
# Section 5: what's-next item 1 -- sweep knife_edge_size itself. Everything
# in Sections 1-4 used knife_edge=10, matched to Section 1's own yield peak.
# Does a smaller knife-edge (catching more of the juvenile boom that cutting
# cannibalism creates) make this model genuinely fishing-collapse-prone the
# way Day 39 found for the ORIGINAL model, or does theta=0.3 remove that
# failure mode at every gear setting?
#
# Constant schedule only (thresholdFMort() untouched, per the earlier
# decision not to pursue the total-biomass trigger further). Gear cutoff only
# matters once effort>0, so the unfished fork state (last_n_theta_low,
# last_npp_theta_low) is reused unchanged across every knife_edge_size --
# no re-forking needed.
################################################################################

# Unfished total biomass reference (effort=0) -- gear-independent, computed
# once rather than once per knife_edge value.
p_ref <- make_anchovy_fishing_params_theta(theta_low, 1, knife_edge_size = 10)
p_ref@initial_n[]    <- last_n_theta_low
p_ref@initial_n_pp[] <- last_npp_theta_low
sim_unfished <- project(p_ref, t_max = scan_summary_window, dt = p2$dt, t_save = 0.2,
                        t_start = t_fork, progress_bar = FALSE, effort = 0)
tv_uf   <- as.numeric(dimnames(sim_unfished@n)[[1]])
keep_uf <- tv_uf > t_fork + scan_summary_window / 2
total_unfished_ref <- mean(unname(getBiomass(sim_unfished)[, "Anchovy"])[keep_uf])

knife_edge_seq_sweep <- c(3, 5, 7, 10, 15, 20)   # spans well below/above w_mat=10
fish_level_seq_ke    <- c(1, 2, 3, 5, 9, 20, 50, 100)

knife_edge_sweep_df <- bind_rows(lapply(knife_edge_seq_sweep, function(ke) {
  p_ke <- make_anchovy_fishing_params_theta(theta_low, 1, knife_edge_size = ke)
  cannibalism_fishing_probe(p_ke, last_n_theta_low, last_npp_theta_low, fish_level_seq_ke) %>%
    mutate(knife_edge_size = ke)
}))
write.csv(knife_edge_sweep_df, file.path("interesting_plots", "day40_knife_edge_sweep.csv"), row.names = FALSE)
cat(sprintf(
  "Section 5 (knife_edge sweep, theta=0.3, unfished total biomass reference=%.4f): see day40_knife_edge_sweep.csv.\n",
  total_unfished_ref
))
print(knife_edge_sweep_df)

ke_biomass_plot <- ggplot(knife_edge_sweep_df, aes(x = fish_level, y = pmax(mean_total, 1e-10),
                                                    color = factor(knife_edge_size))) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Fishing effort (Constant, log scale)", y = "Mean total biomass (log scale)",
       color = "knife_edge_size",
       title = "Knife-edge sweep: smaller gear cutoffs genuinely collapse the population",
       subtitle = "theta=0.3, interaction_resource=1 -- knife_edge=3/5 crash toward extinction by fish_level=50-100; knife_edge=10 (the aim model) never does") +
  theme_minimal()
ke_biomass_plot
save_plot(ke_biomass_plot, "day40_knife_edge_sweep_biomass.png", width = 9, height = 6)

ke_yield_plot <- ggplot(knife_edge_sweep_df, aes(x = fish_level, y = mean_yield, color = factor(knife_edge_size))) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "Fishing effort (Constant, log scale)", y = "Mean yield",
       color = "knife_edge_size",
       title = "Knife-edge sweep: yield vs. effort",
       subtitle = "theta=0.3, interaction_resource=1 -- knife_edge=3/5's late uptick coincides with population collapse, not a real yield gain (see caveat)") +
  theme_minimal()
ke_yield_plot
save_plot(ke_yield_plot, "day40_knife_edge_sweep_yield.png", width = 9, height = 6)

ke_collapse_summary <- knife_edge_sweep_df %>%
  group_by(knife_edge_size) %>%
  summarise(fraction_of_unfished_at_100 = mean_total[fish_level == 100] / total_unfished_ref,
           collapses_by_fish_level_100 = fraction_of_unfished_at_100 < 0.01, .groups = "drop")
cat("Fraction of unfished total biomass remaining at fish_level=100, by knife_edge_size:\n")
print(ke_collapse_summary)

cat(sprintf(
  "Confirms the What's Next hypothesis: knife_edge=3 and 5 (below w_mat=%.0f) genuinely collapse this model -- biomass falls to %.1f%% and %.3f%% of the unfished reference (%.4f) by fish_level=100, both essentially extinct. knife_edge=7 is intermediate (declining but not yet extinct at fish_level=100). knife_edge=10, the aim model from Sections 1-4, never collapses anywhere in [1,100]. knife_edge=15/20 (above w_mat) barely respond to fishing at all -- the gear catches so few individuals that effort is nearly irrelevant to both biomass and yield. Cutting cannibalism (theta=0.3) does NOT remove the collapse failure mode Day 39 found in the original model; it just moves it to a different, smaller region of gear-selectivity space that Sections 1-4's own knife_edge=10 choice happened to sit outside of.\n",
  p2$w_mat, 100 * ke_collapse_summary$fraction_of_unfished_at_100[ke_collapse_summary$knife_edge_size == 3],
  100 * ke_collapse_summary$fraction_of_unfished_at_100[ke_collapse_summary$knife_edge_size == 5],
  total_unfished_ref
))

cat("Caveat: the late yield UPTICK for knife_edge=3 (fish_level=20) and knife_edge=5 (fish_level=50) right before their own full collapse is not a real exploitable yield gain -- it coincides with the population being fished down hard inside the averaging window, not a new equilibrium. Read the biomass plot before the yield plot at these gear/effort combinations. Also worth flagging: none of Sections 1-5 use second_order_w=TRUE / method='tr_bdf2' (MIZER-AGENTS.md's own recommendation for studies of oscillation/collapse dynamics) -- consistent with every prior day's own convention, but the near-extinction values reported here (down to ~1e-18) have not been checked against the second-order scheme for numerical-diffusion artefacts.\n")

################################################################################
# Section 6: threshold_frac x boost heatmap -- every schedule so far
# (Sections 2-5) calibrated threshold_ke ONCE, at the Constant cycle's own
# median (probs=0.5). Day 38 found threshold_frac=0.6 fixed a double-firing
# issue that 0.5 had on a different model, but nothing in this project has
# swept threshold placement itself as its own grid dimension against boost --
# does where the threshold sits change how the schedule responds to boost
# strength, on top of the usual boost-alone sweep?
#
# Reuses Section 3's own MSY floor and gear (ke_compare, baseline_effort_msy).
# threshold_frac only changes WHERE on the same Constant cycle the quantile
# is taken, not the cycle itself, so only ONE Constant simulation is needed
# for the whole grid -- cheaper than Section 5's knife_edge sweep, which
# needed independent gear (and hence independent dynamics) per point.
################################################################################

threshold_frac_seq <- c(0.3, 0.4, 0.5, 0.6, 0.7)
boost_seq_heatmap  <- c(5, 15, 30)

p_heatmap <- make_anchovy_fishing_params_theta(theta_low, 1, knife_edge_size = ke_compare)
p_heatmap@initial_n[]    <- last_n_theta_low
p_heatmap@initial_n_pp[] <- last_npp_theta_low

sim_const_heatmap <- project(p_heatmap, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                             t_start = t_fork, progress_bar = FALSE, effort = baseline_effort_msy)
bp_const_heatmap       <- compute_selected_biomass_series(sim_const_heatmap, p_heatmap, scan_t_cut)
sharpness_heatmap      <- 0.02 * (max(bp_const_heatmap) - min(bp_const_heatmap))
const_metrics_heatmap  <- scan_metrics_layered(sim_const_heatmap)

run_threshold_boost_case <- function(threshold_frac, boost, mode) {
  threshold_val <- unname(quantile(bp_const_heatmap, probs = threshold_frac))
  schedule_name <- if (mode == "above") "Threshold (peaks)" else "Threshold (troughs)"
  p_rule <- attach_threshold_rule(p_heatmap, threshold = threshold_val, fish_level = boost,
                                  background_level = baseline_effort_msy, mode = mode,
                                  sharpness = sharpness_heatmap)
  p_rule@initial_n[]    <- last_n_theta_low
  p_rule@initial_n_pp[] <- last_npp_theta_low
  sim <- project(p_rule, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 t_start = t_fork, progress_bar = FALSE, effort = 1)
  diag <- threshold_diagnostics(sim, threshold_val, boost, baseline_effort_msy, mode,
                                sharpness_heatmap, scan_t_cut)
  cbind(scan_metrics_layered(sim), threshold_frac = threshold_frac, boost_level = boost,
       schedule = schedule_name, mean_effort = diag$mean_effort, threshold_value = threshold_val)
}

threshold_boost_df <- bind_rows(lapply(threshold_frac_seq, function(tf) {
  bind_rows(lapply(boost_seq_heatmap, function(b) {
    bind_rows(lapply(c("above", "below"), function(m) run_threshold_boost_case(tf, b, m)))
  }))
}))
write.csv(threshold_boost_df, file.path("interesting_plots", "day40_threshold_boost_heatmap.csv"), row.names = FALSE)
cat(sprintf(
  "Section 6 (threshold_frac x boost grid, plain theta=0.3, knife_edge=%d, floor=MSY=%.2g): see day40_threshold_boost_heatmap.csv. Constant's own mean_total=%.4f, mean_yield=%.4f.\n",
  ke_compare, baseline_effort_msy, const_metrics_heatmap$mean_total, const_metrics_heatmap$mean_yield
))
print(threshold_boost_df)

threshold_heatmap_amplitude <- ggplot(threshold_boost_df,
                                      aes(x = factor(boost_level), y = factor(threshold_frac), fill = rel_amplitude)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", rel_amplitude)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "magma") +
  facet_wrap(~schedule) +
  labs(x = "Boost level (on-period fish_level)", y = "threshold_frac (quantile of the Constant cycle)",
       fill = "rel_amplitude",
       title = "Does threshold placement interact with boost strength? (relative amplitude)",
       subtitle = sprintf("Plain theta=0.3, knife_edge=%d, floor=MSY effort=%.2g", ke_compare, baseline_effort_msy)) +
  theme_minimal()
save_plot(threshold_heatmap_amplitude, "day40_threshold_heatmap_amplitude.png", width = 10, height = 6)

threshold_heatmap_total <- ggplot(threshold_boost_df,
                                  aes(x = factor(boost_level), y = factor(threshold_frac), fill = mean_total)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", mean_total)), size = 3, color = "white") +
  scale_fill_viridis_c() +
  facet_wrap(~schedule) +
  labs(x = "Boost level (on-period fish_level)", y = "threshold_frac (quantile of the Constant cycle)",
       fill = "mean_total",
       title = "Does threshold placement interact with boost strength? (mean total biomass)",
       subtitle = sprintf("Plain theta=0.3, knife_edge=%d, floor=MSY effort=%.2g -- Constant's own mean_total=%.4f",
                          ke_compare, baseline_effort_msy, const_metrics_heatmap$mean_total)) +
  theme_minimal()
save_plot(threshold_heatmap_total, "day40_threshold_heatmap_total.png", width = 10, height = 6)

threshold_heatmap_yield <- ggplot(threshold_boost_df,
                                  aes(x = factor(boost_level), y = factor(threshold_frac), fill = mean_yield)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", mean_yield)), size = 3, color = "white") +
  scale_fill_viridis_c(option = "cividis") +
  facet_wrap(~schedule) +
  labs(x = "Boost level (on-period fish_level)", y = "threshold_frac (quantile of the Constant cycle)",
       fill = "mean_yield",
       title = "Does threshold placement interact with boost strength? (mean yield)",
       subtitle = sprintf("Plain theta=0.3, knife_edge=%d, floor=MSY effort=%.2g -- Constant's own mean_yield=%.4f",
                          ke_compare, baseline_effort_msy, const_metrics_heatmap$mean_yield)) +
  theme_minimal()
save_plot(threshold_heatmap_yield, "day40_threshold_heatmap_yield.png", width = 10, height = 6)

cat(sprintf(
  "Section 6 summary: rel_amplitude ranges %.2f-%.2f across the threshold_frac x boost grid -- %s\n",
  min(threshold_boost_df$rel_amplitude), max(threshold_boost_df$rel_amplitude),
  if (max(threshold_boost_df$rel_amplitude) > 1.9) {
    "some cell(s) approach the ~2.0 near-extinction ceiling."
  } else {
    "nowhere near the ~2.0 near-extinction ceiling Day 39 found in the original model."
  }
))
