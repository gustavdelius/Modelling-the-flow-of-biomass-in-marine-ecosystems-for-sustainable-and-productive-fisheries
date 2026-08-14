library(mizer)
library(mizerExperimental)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)

# Day 38 ports Day 36's threshold-triggered-fishing machinery onto the one
# Day 37 configuration that actually sustains a limit cycle: the
# plankton-anchovy model (Canales, Delius & Law 2020), cannibalism only,
# mu_l=0 -- Day 37's own "Figure 2e" build. Day 37 checked three variants at
# t>=25 (well clear of the post-kick transient) and only this one keeps
# oscillating: max/min ratio 5.82, against ~1.00-1.03 for cannibalism+larval
# mortality together or larval mortality alone. Those two are dropped here,
# not carried forward "just in case" -- per Day 37's own What's Next, why
# larval mortality kills the cycle is an open question but not this day's.
#
# Day 36's own model (single-species Anchovy, erepro~80,000, broken) and Day
# 4's replacement (checked in Day 37, no oscillation at all) are both gone.
# The threshold-rule machinery itself (thresholdFMort(), attach_threshold_
# rule(), the selected-biomass/on-frac helpers) is generic in params/sim and
# survives the swap unchanged -- it doesn't know or care what model it's
# fishing.
#
# Self-contained convention since Day 20: helpers redefined here, not sourced.

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
# Section 0: rebuilding the Figure 2e plankton-anchovy model, plus fishing
#
# Model-building code (p2 list, setAnchovyMort(), plankton_logistic(),
# norm_box_pred_kernel(), setAnchovyModel()) is Day 37's, verbatim -- Day 37
# already checked the kernel and mortality curve actually built match the
# paper (day37_anchovy_kernel_check.png, day37_anchovy_mort_check.png), so
# that diagnostic isn't repeated here. What Day 37 never needed and this day
# does is a fishing gear: knife-edge at w_mat, catchability=1, the same
# convention Day 36 used (make_fishing_params()) -- w_mat=10 in both models,
# so it's the same cut point.
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
plankton_state$random <- FALSE   # random plankton forcing kept off, as Day 37

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

# Figure 2e config + knife-edge gear at w_mat, catchability=1 -- Day 36's own
# make_fishing_params() convention, carried over.
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
sim <- project(p_scan,t_max = 100,dt=0.01,method="tr_bdf2")
anyNA(initialN(sim))
anyNA(finalN(sim))
plotSpectra(sim)
animateSpectra(sim)

species_params(p_scan)$w_mat25
gear_params(p_scan)$catchability
plotGrowthCurves(p_scan)

anchovy_immigration <- p2$i0 * p_scan@w_full^(-p2$lambda) * exp(-6.9 * (p2$lambda - 1))

cat(sprintf("Figure 2e plankton-anchovy erepro: %.4g (should sit in [0,1])\n",
           species_params(p_scan)$erepro))

################################################################################
# Section 0b: a numerics check, before trusting any of the above
#
# MIZER-AGENTS.md is explicit that the default first-order-upwind flux
# scheme can silently damp real oscillations, and that oscillation studies
# should use second_order_w=TRUE + method="tr_bdf2". Day 37 didn't do this --
# it used the paper's own dt=0.001 with mizer's default scheme throughout.
# Tried switching before building anything else on top of it: rebuilding
# with second_order_w=TRUE and projecting with tr_bdf2 at dt=0.05 produces
# non-finite values within the first 10-year settle; dropping to dt=0.01
# removes the immediate error but the returned array still contains NAs
# (122 of them in the last time slice alone) -- silently, same failure mode
# as thresholdFMort() would inherit if this were used further downstream.
# Given the second-order scheme is unstable here rather than just untried,
# this project keeps Day 37's own dt=0.001, default (first-order) scheme --
# already shown (Day 37) to reproduce the paper's own Figure 2e oscillation
# quantitatively, which is the standard that actually matters here, not
# which scheme is used to get there.
################################################################################

# Same build as setAnchovyModel() above, second_order_w=TRUE added.
kappa_2nd_order <- p2$a0 * exp(-6.9 * (p2$lambda - 1))
no_w_2nd_order  <- round(log(p2$w_inf / p2$w_min) / p2$dx)
p_scan_2nd_order <- newMultispeciesParams(
  species_params(p_scan), no_w = no_w_2nd_order, lambda = p2$lambda,
  kappa = kappa_2nd_order, w_pp_cutoff = p2$w_pp_cutoff,
  resource_dynamics = "plankton_logistic", second_order_w = TRUE
)
p_scan_2nd_order@rr_pp[] <- p2$r0 * p_scan_2nd_order@w_full^(p2$rho - 1)
p_scan_2nd_order <- setAnchovyMort(p_scan_2nd_order, p2)
p_scan_2nd_order@interaction[] <- 1
p_scan_2nd_order@initial_n[]    <- 0.001 * p_scan_2nd_order@w^(-1.8)
p_scan_2nd_order@initial_n_pp[] <- p_scan_2nd_order@cc_pp

numerics_check <- tryCatch(
  project(p_scan_2nd_order, t_max = 10, dt = 0.01, method = "tr_bdf2", progress_bar = FALSE),
  error = function(e) e
)
if (inherits(numerics_check, "error")) {
  cat(sprintf("Numerics check: second_order_w+tr_bdf2 errored during the 10yr settle (%s) -- keeping Day 37's default scheme.\n",
             conditionMessage(numerics_check)))
} else {
  n_bad <- sum(!is.finite(numerics_check@n))
  cat(sprintf("Numerics check: second_order_w+tr_bdf2 ran but returned %d non-finite value(s) -- keeping Day 37's default scheme.\n",
             n_bad))
}

################################################################################
# Section 1: the threshold rule machinery -- Day 36's, unchanged
#
# Generic in params/sim: doesn't know or assume anything about the
# single-species Anchovy-only model it was written against. Ports as-is.
################################################################################

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

################################################################################
# Section 1b: building and calibrating the fork
#
# Day 37's own settle+kick recipe (10yr settle from a power-law abundance,
# 10^7 knockdown across the whole grid at t=10), continued unfished up to
# t_fork. t_fork=20 (not Day 36's t_fork=500) because the cycle is already
# at full amplitude well before then: Day 37's own fig2e trace shows the
# post-kick transient dying out by about t=14-18 (biomass already swinging
# 0.09-0.54 by t=18), not growing further after that. Confirmed below rather
# than assumed.
################################################################################

t_fork               <- 20
scan_post_fork_years <- 24   # ~4-6 periods of the paper's own ~6yr cycle
scan_summary_window  <- 12   # ~2 periods, matching Day 37's own sampling window
scan_t_cut           <- t_fork + scan_post_fork_years - scan_summary_window

make_anchovy_fork_sim <- function(params, t_fork) {
  params@initial_n[]    <- 0.001 * params@w^(-1.8)
  params@initial_n_pp[] <- params@cc_pp
  sim <- project(params, t_max = 10, dt = p2$dt, progress_bar = FALSE)
  sim@n[11, , ] <- sim@n[11, , ] / 10^7
  project(sim, t_max = t_fork - 10, dt = p2$dt, t_save = 0.2, progress_bar = FALSE, effort = 0)
}

sim_fork <- make_anchovy_fork_sim(p_scan, t_fork)
fork_bp  <- compute_selected_biomass_series(sim_fork, p_scan, t_cut = t_fork - scan_summary_window)
cat(sprintf("Fork check (t in [%.0f,%.0f], last %.0fyr of the %.0fyr fork): selected biomass min=%.4g max=%.4g -- should already span close to the full cycle, not still climbing.\n",
           t_fork - scan_summary_window, t_fork, scan_summary_window, t_fork,
           min(fork_bp), max(fork_bp)))

last_n   <- array(sim_fork@n[dim(sim_fork@n)[1], , , drop = FALSE], dim = dim(sim_fork@n)[-1])
last_npp <- sim_fork@n_pp[dim(sim_fork@n_pp)[1], ]

################################################################################
# Section 2: the scan
#
# Day 36's run_window_effort_scan(), adapted: t_start=t_fork (not restarted
# at 0) since branches are fresh MizerParams objects, not sim continuations;
# dt=p2$dt throughout, no method= override (Section 0b's numerics check).
# threshold_frac calibrated as a quantile of the fork's own unfished
# continuation, same reasoning as Day 36.
################################################################################

run_anchovy_threshold_scan <- function(fish_level_seq, threshold_frac_seq = NULL,
                                       background_level_seq = 0,
                                       schedules = c("Threshold (peaks)", "Threshold (troughs)", "Constant")) {
  scan_metrics <- function(sim) {
    bmv <- getBiomass(sim)[, "Anchovy"]
    ylv <- getYield(sim)[, "Anchovy"]
    tv  <- as.numeric(names(bmv))
    bm_late <- bmv[tv > scan_t_cut]
    yl_late <- ylv[tv > scan_t_cut]
    data.frame(
      rel_amplitude = (max(bm_late) - min(bm_late)) / ((max(bm_late) + min(bm_late)) / 2),
      mean_yield    = mean(yl_late)
    )
  }

  out <- list()

  if ("Constant" %in% schedules) {
    out$constant <- bind_rows(lapply(fish_level_seq, function(fl) {
      sim <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                    progress_bar = FALSE, effort = fl)
      cbind(scan_metrics(sim), threshold_frac = NA_real_, threshold_bp = NA_real_,
           fish_level = fl, schedule = "Constant",
           mean_effort = fl, background_level = NA_real_,
           effective_window = NA_real_, n_bursts = NA_integer_)
    }))
  }

  needs_threshold <- any(c("Threshold (peaks)", "Threshold (troughs)") %in% schedules) &&
    !is.null(threshold_frac_seq)

  if (needs_threshold) {
    sim_ref_unfished <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                                progress_bar = FALSE, effort = 0)
    bp_ref <- compute_selected_biomass_series(sim_ref_unfished, p_scan, scan_t_cut)
    sharpness_bp <- 0.02 * (max(bp_ref) - min(bp_ref))

    out$threshold <- bind_rows(lapply(fish_level_seq, function(fl) {
      bind_rows(lapply(threshold_frac_seq, function(frac) {
        threshold_bp <- unname(quantile(bp_ref, probs = frac))
        bind_rows(lapply(background_level_seq, function(bg) {
          rows <- list()
          if ("Threshold (peaks)" %in% schedules) {
            p_peaks <- attach_threshold_rule(p_scan, threshold = threshold_bp, fish_level = fl,
                                             background_level = bg, mode = "above",
                                             sharpness = sharpness_bp)
            p_peaks@initial_n[]    <- last_n
            p_peaks@initial_n_pp[] <- last_npp
            sim_p <- project(p_peaks, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                             t_start = t_fork, progress_bar = FALSE, effort = 1)
            diag_p <- threshold_diagnostics(sim_p, threshold_bp, fl, bg, "above",
                                            sharpness_bp, scan_t_cut)
            rows$peaks <- cbind(scan_metrics(sim_p), threshold_frac = frac, threshold_bp = threshold_bp,
                                fish_level = fl, schedule = "Threshold (peaks)", mean_effort = diag_p$mean_effort,
                                background_level = bg, effective_window = diag_p$effective_window,
                                n_bursts = diag_p$n_bursts)
          }
          if ("Threshold (troughs)" %in% schedules) {
            p_troughs <- attach_threshold_rule(p_scan, threshold = threshold_bp, fish_level = fl,
                                               background_level = bg, mode = "below",
                                               sharpness = sharpness_bp)
            p_troughs@initial_n[]    <- last_n
            p_troughs@initial_n_pp[] <- last_npp
            sim_t <- project(p_troughs, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                             t_start = t_fork, progress_bar = FALSE, effort = 1)
            diag_t <- threshold_diagnostics(sim_t, threshold_bp, fl, bg, "below",
                                            sharpness_bp, scan_t_cut)
            rows$troughs <- cbind(scan_metrics(sim_t), threshold_frac = frac, threshold_bp = threshold_bp,
                                  fish_level = fl, schedule = "Threshold (troughs)", mean_effort = diag_t$mean_effort,
                                  background_level = bg, effective_window = diag_t$effective_window,
                                  n_bursts = diag_t$n_bursts)
          }
          bind_rows(rows)
        }))
      }))
    }))
  }

  bind_rows(out)
}

################################################################################
# Section 2b: sanity check -- does the switch fire where it should?
#
# Day 36's own Section 1b, before spending time on the grids below: one
# representative fish_level/threshold_frac, all three schedules, selected
# biomass plotted against the calibration threshold with fishing "on"
# periods shaded.
################################################################################

sanity_fish_level     <- 3
sanity_threshold_frac <- 0.6   # raised from 0.5 -- median fired on secondary, spurious
                               # crossings as well as the real peak (two firing windows
                               # per cycle in Threshold (peaks) instead of one). At 0.6
                               # each cycle gets exactly one clean window right at the
                               # peak, without materially touching the peak's own height.

sim_ref_sanity   <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                            progress_bar = FALSE, effort = 0)
bp_ref_sanity    <- compute_selected_biomass_series(sim_ref_sanity, p_scan, scan_t_cut)
sharpness_sanity <- 0.02 * (max(bp_ref_sanity) - min(bp_ref_sanity))
threshold_sanity <- unname(quantile(bp_ref_sanity, probs = sanity_threshold_frac))

run_sanity_case <- function(schedule) {
  if (schedule == "Constant") {
    sim  <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                    progress_bar = FALSE, effort = sanity_fish_level)
    mode <- NA_character_
  } else {
    mode <- if (schedule == "Threshold (peaks)") "above" else "below"
    p <- attach_threshold_rule(p_scan, threshold = threshold_sanity, fish_level = sanity_fish_level,
                               background_level = 0, mode = mode, sharpness = sharpness_sanity)
    p@initial_n[]    <- last_n
    p@initial_n_pp[] <- last_npp
    sim <- project(p, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                   t_start = t_fork, progress_bar = FALSE, effort = 1)
  }
  list(sim = sim, schedule = schedule, mode = mode)
}

sanity_schedules <- c("Constant", "Threshold (peaks)", "Threshold (troughs)")
sanity_cases     <- lapply(sanity_schedules, run_sanity_case)
names(sanity_cases) <- sanity_schedules

sanity_series_df <- bind_rows(lapply(sanity_cases, function(case) {
  tv   <- as.numeric(dimnames(case$sim@n)[[1]])
  keep <- which(tv > t_fork)
  bp   <- compute_selected_biomass_series(case$sim, p_scan, t_fork)
  on_frac <- if (is.na(case$mode)) {
    rep(NA_real_, length(bp))
  } else {
    compute_on_frac_series(case$sim, p_scan, threshold_sanity, case$mode, sharpness_sanity, t_fork)
  }
  data.frame(t = tv[keep], selected_biomass = bp, on_frac = on_frac, schedule = case$schedule)
}))

write.csv(sanity_series_df, file.path("interesting_plots", "day38_sanity_series.csv"), row.names = FALSE)

sanity_check_plot <- ggplot(sanity_series_df, aes(x = t, y = selected_biomass)) +
  geom_rect(data = sanity_series_df %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = threshold_sanity, linetype = "dashed", color = "grey40") +
  facet_wrap(~schedule, ncol = 1) +
  labs(x = "Time (years)", y = "Selected biomass",
       title = "Sanity check: does thresholdFMort() actually switch on where it should?",
       subtitle = sprintf("Figure 2e plankton-anchovy model, fish_level=%.2g, threshold_frac=%.2g -- dashed line = calibration threshold, red shading = fishing 'on'",
                          sanity_fish_level, sanity_threshold_frac)) +
  theme_minimal()
sanity_check_plot
save_plot(sanity_check_plot, "day38_sanity_check.png", width = 9, height = 8)

cat(sprintf(
  "Section 2b (sanity check, fish_level=%.2g, threshold_frac=%.2g): day38_sanity_check.png shows selected_biomass(t) against the calibration threshold, fishing 'on' periods shaded.\n",
  sanity_fish_level, sanity_threshold_frac
))

################################################################################
# Section 3 (Improvement A): effort scan, threshold_frac fixed at the median
#
# Mirrors Day 36 Section 2. fish_level_seq is much smaller than Day 36's
# (5 points vs. 20) and reaches higher (up to 9, not 8) -- each point here
# costs a ~24yr run at dt=0.001 (~1.6s/simulated year), not Day 36's cheap
# dt=0.1/tr_bdf2 point, and Day 37's own effort sweep on this exact model
# found the cycle still going at effort=3, so fish_level needs to clear that
# before there's anything to see.
################################################################################

threshold_frac_a <- 0.5
fish_level_seq_a <- c(1, 3, 5, 7, 9)

scan_a_df <- run_anchovy_threshold_scan(fish_level_seq = fish_level_seq_a,
                                        threshold_frac_seq = threshold_frac_a)

write.csv(scan_a_df, file.path("interesting_plots", "day38_effort_scan.csv"), row.names = FALSE)

effort_amplitude_plot <- ggplot(scan_a_df, aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Relative amplitude of biomass",
       title = "Does threshold-triggered fishing damp the cannibalism-driven cycle?",
       subtitle = sprintf("Figure 2e plankton-anchovy model, threshold_frac=%.2g", threshold_frac_a)) +
  theme_minimal()
effort_amplitude_plot
save_plot(effort_amplitude_plot, "day38_effort_amplitude.png", width = 9, height = 6)

effort_yield_plot <- ggplot(scan_a_df, aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Mean yield",
       title = "Yield vs. mean effort, Figure 2e plankton-anchovy model",
       subtitle = sprintf("threshold_frac=%.2g", threshold_frac_a)) +
  theme_minimal()
effort_yield_plot
save_plot(effort_yield_plot, "day38_effort_yield.png", width = 9, height = 6)

cat(sprintf(
  "Section 3 (effort scan, fish_level=[%.2g,%.2g], threshold_frac=%.2g): see day38_effort_amplitude.png / day38_effort_yield.png / day38_effort_scan.csv.\n",
  min(fish_level_seq_a), max(fish_level_seq_a), threshold_frac_a
))

################################################################################
# Section 3b: Constant-only wide scan -- where does yield start decreasing?
#
# Day 37's own effort sweep on this model only went to effort=3, and yield
# was still rising at that boundary -- item 1 on Day 37's own What's Next.
# Constant only, not the threshold schedules -- no threshold calibration run
# needed, so cost is just fish_level_seq's own length x scan_post_fork_years,
# which is what buys room for a genuinely wide, fine scan here rather than
# the 5-point grids above.
#
# Two passes rather than one large fixed grid: a coarse sweep across a wide
# range to bracket where mean_yield stops rising, then a fine sweep zoomed
# on that bracket. Guessing the right resolution for the whole range up
# front would mean either wasting runs far from the peak or missing it
# between two coarse points.
################################################################################

constant_scan_metrics <- function(fl) {
  sim <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 progress_bar = FALSE, effort = fl)
  ylv <- getYield(sim)[, "Anchovy"]
  bmv <- getBiomass(sim)[, "Anchovy"]
  tv  <- as.numeric(names(ylv))
  yl_late <- ylv[tv > scan_t_cut]
  bm_late <- bmv[tv > scan_t_cut]
  data.frame(fish_level = fl, mean_yield = mean(yl_late),
            yield_min = min(yl_late), yield_max = max(yl_late),
            mean_biomass = mean(bm_late))
}

constant_coarse_seq <- seq(0, 30, by = 2)
constant_coarse_df  <- bind_rows(lapply(constant_coarse_seq, constant_scan_metrics))
write.csv(constant_coarse_df, file.path("interesting_plots", "day38_constant_coarse_scan.csv"), row.names = FALSE)

coarse_peak_idx <- which.max(constant_coarse_df$mean_yield)
coarse_peak_fl  <- constant_coarse_df$fish_level[coarse_peak_idx]
cat(sprintf(
  "Constant-only coarse pass (fish_level=[%.2g,%.2g] by %.2g): mean_yield peaks at fish_level=%.2g (%.4g). Refining around it.\n",
  min(constant_coarse_seq), max(constant_coarse_seq), diff(constant_coarse_seq)[1],
  coarse_peak_fl, constant_coarse_df$mean_yield[coarse_peak_idx]
))

# Fine pass: one coarse step either side of the coarse peak, fine resolution.
# If the coarse peak sits at either end of constant_coarse_seq, that end
# isn't a real interior maximum -- flagged below rather than silently
# reported as if it were.
coarse_step        <- diff(constant_coarse_seq)[1]
constant_fine_seq  <- seq(max(0, coarse_peak_fl - coarse_step), coarse_peak_fl + coarse_step, by = 0.25)
constant_fine_df   <- bind_rows(lapply(constant_fine_seq, constant_scan_metrics))
write.csv(constant_fine_df, file.path("interesting_plots", "day38_constant_fine_scan.csv"), row.names = FALSE)

constant_scan_df <- bind_rows(constant_coarse_df, constant_fine_df) %>%
  distinct(fish_level, .keep_all = TRUE) %>%
  arrange(fish_level)
write.csv(constant_scan_df, file.path("interesting_plots", "day38_constant_scan_combined.csv"), row.names = FALSE)

fine_peak_idx <- which.max(constant_scan_df$mean_yield)
fine_peak_fl  <- constant_scan_df$fish_level[fine_peak_idx]
at_edge       <- coarse_peak_fl %in% range(constant_coarse_seq)

constant_collapse_plot <- ggplot(constant_scan_df, aes(x = fish_level, y = mean_yield)) +
  geom_ribbon(aes(ymin = yield_min, ymax = yield_max), fill = "steelblue", alpha = 0.25) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_vline(xintercept = fine_peak_fl, linetype = "dashed", color = "grey40") +
  labs(x = "Fishing effort (Constant schedule)", y = "Mean yield",
       title = "Where does Constant-fishing yield start decreasing?",
       subtitle = sprintf("Figure 2e plankton-anchovy model -- coarse [%.2g,%.2g] by %.2g, fine step 0.25 around the coarse peak; band = min/max over the last %.0fyr of each run; dashed = peak at fish_level=%.2g",
                          min(constant_coarse_seq), max(constant_coarse_seq), coarse_step,
                          scan_summary_window, fine_peak_fl)) +
  theme_minimal()
constant_collapse_plot
save_plot(constant_collapse_plot, "day38_constant_collapse.png", width = 9, height = 6)

if (at_edge) {
  cat(sprintf(
    "Section 3b (Constant-only scan): coarse peak sits at the edge of the range tested (fish_level=%.2g) -- not a real interior maximum, extend constant_coarse_seq further before trusting this number.\n",
    coarse_peak_fl
  ))
} else {
  cat(sprintf(
    "Section 3b (Constant-only scan): mean_yield peaks at fish_level=%.2g (mean_yield=%.4g) -- see day38_constant_collapse.png / day38_constant_scan_combined.csv.\n",
    fine_peak_fl, constant_scan_df$mean_yield[fine_peak_idx]
  ))
}

################################################################################
# Section 3c: diagnostics -- does the population actually collapse at high
# effort, or is Section 3b's near-flat yield curve hiding something?
#
# Section 3b's own plot is not what naive intuition expects: mean_yield
# barely moves between fish_level=10 and fish_level=30 rather than turning
# over and heading toward zero the way heavy sustained fishing "should".
# Worth checking directly rather than either trusting the plot blind or
# assuming it's wrong. Two real, different explanations are on the table --
# a genuine ecological effect, or a scan artefact -- and only looking at the
# population itself (not just the yield summary) can tell them apart:
#
#   1. A "cultivation effect". The gear is knife-edge at w_mat=10 and never
#      touches anything smaller, and Day 37's whole finding was that
#      cannibalism (large Anchovy eating small Anchovy) is what drives this
#      model's cycle. Fishing out exactly the size class that does the
#      cannibalising could relieve predation on juveniles rather than just
#      remove biomass -- a real, published effect in cannibalistic
#      size-spectrum models, not a mizer quirk.
#   2. A scan artefact -- gear/effort not actually scaling as intended,
#      or scan_post_fork_years=24 not being long enough to reach a new
#      equilibrium at high effort, so the "flat" yield is a transient that
#      hasn't finished declining yet, not the real long-run outcome.
#
# CONFIRMED (2), then (1): the gear check below is clean (F(w) knife-edges
# exactly at w_mat, scales exactly linearly with effort -- no bug), and the
# long-run check confirms 24yr is enough time (0.379 at 24yr vs. 0.415 at
# 100yr, well within noise). The mechanism is (1). mean_juvenile tracks
# mean_total almost exactly at every fish_level tested (juveniles are >95%
# of biomass throughout) and is NOT monotonically falling with effort --
# 0.374 unfished, 0.485 at fish_level=10 (fishing the cannibalistic adults
# measurably *helps* the juveniles), still 0.291 at a fish_level=1000 stress
# test that crushes mean_selected by >200x (0.0178 -> 0.0000836). Total
# biomass keeps cycling with healthy amplitude at every fish_level tested,
# including 1000 -- see day38_extinction_check_series.png. The juvenile
# compartment is structurally insulated from fishing twice over: the gear
# can't reach it, and removing its main predator (adult conspecifics)
# helps it rather than hurting it. That's why Section 3b's yield never
# turns over within the range tested -- not a bug, a real consequence of
# harvesting a cannibal.
################################################################################

# 1. Gear check: does F(w) actually knife-edge at w_mat and scale linearly
# with effort, as make_anchovy_fishing_params() intended?
gear_check_df <- bind_rows(lapply(c(1, 10, 30), function(fl) {
  # as.numeric(unname(...)) -- mizerFMortGear() returns a gear x size
  # matrix, possibly a Matrix::Matrix rather than a base matrix; colSums()
  # on it can come back with a names attribute and/or a class that survives
  # unname() alone, which trips data.frame()'s row-name auto-detection and
  # silently drops the column. as.numeric() forces a plain atomic vector.
  data.frame(w = p_scan@w, effort = fl, Fw = as.numeric(unname(colSums(mizerFMortGear(p_scan, effort = fl)))))
}))

gear_check_plot <- ggplot(gear_check_df, aes(x = w, y = Fw, color = factor(effort))) +
  geom_line() +
  geom_vline(xintercept = p_scan@species_params$w_mat, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(x = "Body mass (g)", y = "Fishing mortality F(w) [1/year]", color = "effort",
       title = "Gear check: does F(w) actually knife-edge at w_mat and scale with effort?",
       subtitle = sprintf("Dashed = w_mat=%.3g -- should be exactly 0 below it, exactly effort x catchability above it",
                          p_scan@species_params$w_mat)) +
  theme_minimal()
gear_check_plot
save_plot(gear_check_plot, "day38_gear_check.png", width = 8, height = 5)

f_above <- vapply(c(1, 10, 30), function(fl) {
  gear_check_df$Fw[gear_check_df$effort == fl & gear_check_df$w > p_scan@species_params$w_mat][1]
}, numeric(1))
cat(sprintf(
  "Section 3c gear check: F(w) just above w_mat at effort=1/10/30 -> %.4g/%.4g/%.4g (should be exactly 1/10/30, catchability=1).\n",
  f_above[1], f_above[2], f_above[3]
))

# 2. Population structure at high effort -- total, selected (w>=w_mat,
# fishable), and juvenile (w<w_mat, never fished) biomass over time, plus
# the final size spectrum. fish_level=1000 is a deliberate stress test: if
# even that doesn't crash the population, persistence at effort=30 is
# structural, not "just needed a bit more effort".
diag_fish_levels <- c(0, 10, 20, 30, 1000)

run_diag_case <- function(fl) {
  sim <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 progress_bar = FALSE, effort = fl)
  tv          <- as.numeric(dimnames(sim@n)[[1]])
  # unname() throughout -- same row-name auto-detection trap as the gear
  # check above; these are all named (by time or by size) vectors coming
  # straight out of mizer accessors.
  total_bm    <- unname(getBiomass(sim)[, "Anchovy"])
  selected_bm <- unname(getBiomass(sim, min_w = p_scan@species_params$w_mat)[, "Anchovy"])
  list(
    series = data.frame(t = tv, total = total_bm, selected = selected_bm,
                        juvenile = total_bm - selected_bm, fish_level = fl),
    final_spectrum = data.frame(w = p_scan@w, abundance = unname(sim@n[dim(sim@n)[1], 1, ]),
                                fish_level = fl)
  )
}

diag_cases <- lapply(diag_fish_levels, run_diag_case)

diag_series_df   <- bind_rows(lapply(diag_cases, `[[`, "series"))
diag_spectrum_df <- bind_rows(lapply(diag_cases, `[[`, "final_spectrum"))

write.csv(diag_series_df, file.path("interesting_plots", "day38_extinction_check_series.csv"), row.names = FALSE)
write.csv(diag_spectrum_df, file.path("interesting_plots", "day38_extinction_check_spectrum.csv"), row.names = FALSE)

diag_series_long <- diag_series_df %>%
  tidyr::pivot_longer(c(total, selected, juvenile), names_to = "compartment", values_to = "biomass")

extinction_series_plot <- ggplot(diag_series_long %>% filter(t >= t_fork),
                                 aes(x = t, y = pmax(biomass, 1e-12), color = compartment)) +
  geom_line() +
  facet_wrap(~fish_level, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  labs(x = "Time (years)", y = "Biomass (log scale, floored at 1e-12 for plotting)",
       title = "Does the population actually collapse at high fishing effort?",
       subtitle = "Total vs. selected (w>=w_mat, fishable) vs. juvenile (w<w_mat, never fished) biomass") +
  theme_minimal()
extinction_series_plot
save_plot(extinction_series_plot, "day38_extinction_check_series.png", width = 11, height = 7)

extinction_spectrum_plot <- ggplot(diag_spectrum_df %>% filter(abundance > 0),
                                   aes(x = w, y = abundance, color = factor(fish_level))) +
  geom_line() +
  geom_vline(xintercept = p_scan@species_params$w_mat, linetype = "dashed", color = "grey40") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Body mass (g)", y = "Abundance density (log scale)", color = "fish_level",
       title = "Final size spectrum at each fish_level",
       subtitle = sprintf("Dashed = w_mat=%.3g, the fishing cutoff -- does abundance actually crash above it?",
                          p_scan@species_params$w_mat)) +
  theme_minimal()
extinction_spectrum_plot
save_plot(extinction_spectrum_plot, "day38_extinction_check_spectrum.png", width = 9, height = 6)

extinction_summary <- diag_series_df %>%
  filter(t > scan_t_cut) %>%
  group_by(fish_level) %>%
  summarise(min_total = min(total), mean_total = mean(total),
           min_selected = min(selected), mean_selected = mean(selected),
           mean_juvenile = mean(juvenile), .groups = "drop")
write.csv(extinction_summary, file.path("interesting_plots", "day38_extinction_check_summary.csv"), row.names = FALSE)

cat(sprintf(
  "Section 3c population check (last %.0fyr, by fish_level): see day38_extinction_check_summary.csv / day38_extinction_check_series.png / day38_extinction_check_spectrum.png.\n",
  scan_summary_window
))
print(extinction_summary)

extreme_row <- extinction_summary[extinction_summary$fish_level == max(diag_fish_levels), ]
cat(sprintf(
  "At the stress-test effort (fish_level=%.0f): mean total biomass=%.4g, mean selected (fishable) biomass=%.4g, mean juvenile biomass=%.4g. %s\n",
  extreme_row$fish_level, extreme_row$mean_total, extreme_row$mean_selected, extreme_row$mean_juvenile,
  if (extreme_row$mean_total < 1e-6) {
    "Population does crash to near-zero eventually -- effort=30's own persistence is a matter of degree, not a floor that can't be crossed."
  } else if (extreme_row$mean_juvenile > extreme_row$mean_selected) {
    "Selected/fishable biomass is suppressed but juvenile biomass persists (or grows) -- consistent with the cultivation-effect hypothesis: fishing removes cannibalistic adults faster than it removes the population's ability to reproduce."
  } else {
    "Neither a clean crash nor a clean juvenile refuge -- worth a closer look at the series/spectrum plots directly before concluding either way."
  }
))

# 3. Is scan_post_fork_years=24 actually long enough to reach a new
# equilibrium at high effort, or is the population still declining when the
# window closes (a transient, not the real long-run outcome)? One long run
# at a representative high effort, comparing Section 3b's own window against
# a much longer one.
long_check_years <- 100
sim_long_check <- project(sim_fork, t_max = long_check_years, dt = p2$dt, t_save = 0.5,
                          progress_bar = FALSE, effort = 30)
tv_long <- as.numeric(dimnames(sim_long_check@n)[[1]])
total_bm_long <- getBiomass(sim_long_check)[, "Anchovy"]
early_mean <- mean(total_bm_long[tv_long > scan_t_cut & tv_long <= t_fork + scan_post_fork_years])
late_mean  <- mean(total_bm_long[tv_long > t_fork + long_check_years - scan_summary_window])

cat(sprintf(
  "Section 3c long-run check (fish_level=30): mean total biomass over Section 3b's own %.0fyr window=%.4g; mean over the last %.0fyr of a %.0fyr run=%.4g. %s\n",
  scan_post_fork_years, early_mean, scan_summary_window, long_check_years, late_mean,
  if (abs(late_mean - early_mean) / early_mean < 0.1) {
    "Consistent -- scan_post_fork_years=24 looks like enough time to see the real long-run behaviour, not a transient."
  } else {
    "NOT consistent -- 24yr may be too short at this effort; Section 3b's numbers there could still be moving."
  }
))

################################################################################
# Section 4 (Improvement B): threshold-range scan
#
# Mirrors Day 36 Section 3, scaled down the same way: 5 threshold_frac x
# 3 fish_level x 2 schedules = 30 branch runs, not Day 36's 5 x 10 = 50.
################################################################################

threshold_frac_seq_b <- seq(0.1, 0.9, by = 0.2)
fish_level_seq_b     <- c(1, 3, 5)

scan_b_df <- run_anchovy_threshold_scan(fish_level_seq = fish_level_seq_b,
                                        threshold_frac_seq = threshold_frac_seq_b,
                                        schedules = c("Threshold (peaks)", "Threshold (troughs)"))

write.csv(scan_b_df, file.path("interesting_plots", "day38_threshold_scan.csv"), row.names = FALSE)

threshold_amplitude_plot <- ggplot(scan_b_df, aes(x = effective_window, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = interaction(schedule, fish_level), alpha = fish_level), linewidth = 0.8) +
  geom_point(aes(alpha = fish_level), size = 1.5) +
  labs(x = "Effective window (years per on-burst, reconstructed)",
       y = "Relative amplitude of biomass",
       title = "Oscillation amplitude vs. the implicit window length",
       subtitle = "Figure 2e plankton-anchovy model -- window is not set directly, this is what threshold_frac actually produced") +
  theme_minimal()
threshold_amplitude_plot
save_plot(threshold_amplitude_plot, "day38_threshold_amplitude.png", width = 10, height = 6)

threshold_yield_plot <- ggplot(scan_b_df, aes(x = effective_window, y = mean_yield, color = schedule)) +
  geom_line(aes(group = interaction(schedule, fish_level), alpha = fish_level), linewidth = 0.8) +
  geom_point(aes(alpha = fish_level), size = 1.5) +
  labs(x = "Effective window (years per on-burst, reconstructed)",
       y = "Mean yield",
       title = "Yield vs. the implicit window length",
       subtitle = "Figure 2e plankton-anchovy model") +
  theme_minimal()
threshold_yield_plot
save_plot(threshold_yield_plot, "day38_threshold_yield.png", width = 10, height = 6)

cat(sprintf(
  "Section 4 (threshold scan, threshold_frac=[%.2g,%.2g], fish_level=[%.2g,%.2g]): see day38_threshold_amplitude.png / day38_threshold_yield.png / day38_threshold_scan.csv.\n",
  min(threshold_frac_seq_b), max(threshold_frac_seq_b), min(fish_level_seq_b), max(fish_level_seq_b)
))

################################################################################
# Section 5: summary
################################################################################

cat("\n===== Day 38 summary =====\n")
cat("Threshold-triggered fishing (Day 36) rebuilt on the Figure 2e plankton-anchovy model (Day 37) -- the one configuration of the three checked that actually sustains a limit cycle.\n")
cat(sprintf(
  "Section 3 (effort scan, fish_level=[%.2g,%.2g], threshold_frac=%.2g): day38_effort_amplitude.png / day38_effort_yield.png / day38_effort_scan.csv.\n",
  min(fish_level_seq_a), max(fish_level_seq_a), threshold_frac_a
))
cat(sprintf(
  "Section 3b (Constant-only scan, coarse [%.2g,%.2g] + fine around the peak): day38_constant_collapse.png / day38_constant_scan_combined.csv.\n",
  min(constant_coarse_seq), max(constant_coarse_seq)
))
cat("Section 3c (diagnostics -- gear check, population structure at high effort, transient-length check): day38_gear_check.png / day38_extinction_check_series.png / day38_extinction_check_spectrum.png / day38_extinction_check_summary.csv.\n")
cat(sprintf(
  "Section 4 (threshold scan, threshold_frac=[%.2g,%.2g], fish_level=[%.2g,%.2g]): day38_threshold_amplitude.png / day38_threshold_yield.png / day38_threshold_scan.csv.\n",
  min(threshold_frac_seq_b), max(threshold_frac_seq_b), min(fish_level_seq_b), max(fish_level_seq_b)
))
cat("Deferred from Day 36's fuller scan set, given this model's cost per run (~1.6s/simulated year vs. Day 36's cheap dt=0.1/tr_bdf2): hard-step-vs-smooth, background-effort-while-off, and the three-improvements-stacked heatmap. Worth doing once the effort/threshold scans above say where the interesting region actually is, rather than covering the whole space blind.\n")

################################################################################
# Section 6: toward a collapse-capable single-species model
#
# Section 3c found the Figure 2e model structurally can't be fished toward
# extinction: the plankton resource is topped up every year by a constant
# external immigration term (anchovy_immigration, from p2$i0=100) regardless
# of grazing pressure, so once cannibalism eases the juvenile compartment is
# never actually resource-limited -- a textbook "cultivation effect"
# (de Roos & Persson 2002, PNAS; de Roos, Persson & Thieme 2003, Proc. Roy.
# Soc.). The same literature shows size-structured cannibalistic
# populations can sit in the opposite regime too -- an emergent Allee
# effect, where harvesting destabilises the population instead of helping
# it -- but that regime needs juveniles and adults genuinely competing for
# a shared, FINITE resource, not one topped up externally every year
# regardless of what the population does to it.
#
# Staying single-species, staying with cannibalism -- just removing the one
# thing in this model's resource dynamics that isn't actually finite. Most
# surgical change available: i0=0, nothing else touched. Exploratory: built
# and sanity-checked here (does it even still oscillate unfished? does the
# Section 3c juvenile-refuge pattern still hold?), not yet run through the
# full Section 3b/3c scan apparatus -- that's the natural next step if this
# probe shows the picture actually changing.
################################################################################

p2_no_immigration <- p2
p2_no_immigration$i0 <- 0

make_anchovy_fishing_params_no_immigration <- function(p = p2_no_immigration) {
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

p_scan_ni <- make_anchovy_fishing_params_no_immigration()

# anchovy_immigration is looked up BY NAME inside plankton_logistic() (see
# Section 0) -- every project() call reads whatever the global currently
# holds, not a value baked into p_scan_ni. Swapping it here changes what any
# subsequent project() call sees, on ANY params object, until it's swapped
# back -- restored at the end of this section. Same formula Section 0 uses,
# not a hand-built zero vector, so this stays correct if p2's other fields
# ever change.
anchovy_immigration_with_subsidy <- anchovy_immigration
anchovy_immigration <- p2_no_immigration$i0 * p_scan_ni@w_full^(-p2_no_immigration$lambda) *
  exp(-6.9 * (p2_no_immigration$lambda - 1))

cat(sprintf("No-immigration variant erepro: %.4g (should sit in [0,1])\n",
           species_params(p_scan_ni)$erepro))

# 1. Does it still oscillate at all without fishing, or did removing the
# subsidy kill the cycle outright? Same settle-and-kick recipe as the rest
# of the day, unfished.
sim_fork_ni <- make_anchovy_fork_sim(p_scan_ni, t_fork)
fork_bp_ni  <- compute_selected_biomass_series(sim_fork_ni, p_scan_ni, t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "No-immigration variant, unfished, last %.0fyr of the %.0fyr fork: selected biomass min=%.4g max=%.4g.\n",
  scan_summary_window, t_fork, min(fork_bp_ni), max(fork_bp_ni)
))

total_bm_ni_fork <- unname(getBiomass(sim_fork_ni)[, "Anchovy"])
tv_ni_fork       <- as.numeric(dimnames(sim_fork_ni@n)[[1]])
keep_ni_fork     <- tv_ni_fork > t_fork - scan_summary_window
cat(sprintf(
  "No-immigration variant, unfished, total biomass over the same window: min=%.4g max=%.4g (ratio=%.3g).\n",
  min(total_bm_ni_fork[keep_ni_fork]), max(total_bm_ni_fork[keep_ni_fork]),
  max(total_bm_ni_fork[keep_ni_fork]) / min(total_bm_ni_fork[keep_ni_fork])
))

# 2. Small Constant-only effort probe -- not Section 3b's full coarse+fine
# scan, just enough points to see whether the juvenile-refuge pattern from
# Section 3c changes shape before committing to a bigger scan.
ni_probe_levels <- c(0, 5, 10, 20, 50)

run_ni_diag_case <- function(fl) {
  sim <- project(sim_fork_ni, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 progress_bar = FALSE, effort = fl)
  tv          <- as.numeric(dimnames(sim@n)[[1]])
  total_bm    <- unname(getBiomass(sim)[, "Anchovy"])
  selected_bm <- unname(getBiomass(sim, min_w = p_scan_ni@species_params$w_mat)[, "Anchovy"])
  data.frame(t = tv, total = total_bm, selected = selected_bm,
            juvenile = total_bm - selected_bm, fish_level = fl)
}

ni_series_df <- bind_rows(lapply(ni_probe_levels, run_ni_diag_case))
write.csv(ni_series_df, file.path("interesting_plots", "day38_no_immigration_series.csv"), row.names = FALSE)

ni_summary <- ni_series_df %>%
  filter(t > scan_t_cut) %>%
  group_by(fish_level) %>%
  summarise(min_total = min(total), mean_total = mean(total),
           mean_selected = mean(selected), mean_juvenile = mean(juvenile), .groups = "drop")
write.csv(ni_summary, file.path("interesting_plots", "day38_no_immigration_summary.csv"), row.names = FALSE)

ni_series_long <- ni_series_df %>%
  tidyr::pivot_longer(c(total, selected, juvenile), names_to = "compartment", values_to = "biomass")

ni_series_plot <- ggplot(ni_series_long %>% filter(t >= t_fork),
                         aes(x = t, y = pmax(biomass, 1e-12), color = compartment)) +
  geom_line() +
  facet_wrap(~fish_level, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  labs(x = "Time (years)", y = "Biomass (log scale, floored at 1e-12 for plotting)",
       title = "No-immigration variant (i0=0): does removing the plankton subsidy change the picture?",
       subtitle = "Total vs. selected (w>=w_mat, fishable) vs. juvenile (w<w_mat, never fished) biomass") +
  theme_minimal()
ni_series_plot
save_plot(ni_series_plot, "day38_no_immigration_series.png", width = 11, height = 7)

cat("Section 6 (no-immigration probe, fish_level=[0,5,10,20,50]): see day38_no_immigration_series.png / day38_no_immigration_summary.csv.\n")
print(ni_summary)

# mean_juvenile still tracking mean_total (and still rising or flat with
# effort, the way Section 3c's did) means the subsidy wasn't the load-
# bearing piece after all -- next lever to try would be resource
# productivity itself (kappa/a0), not just the immigration term. Falling
# mean_juvenile with effort -- especially if min_total starts approaching
# zero at the higher probe levels -- means this is the direction that
# actually breaks the cultivation effect, and it's worth running the full
# Section 3b/3c apparatus on p_scan_ni properly.
ni_juvenile_trend <- cor(ni_summary$fish_level, ni_summary$mean_juvenile)
cat(sprintf(
  "No-immigration variant: correlation of fish_level with mean_juvenile = %.3f. %s\n",
  ni_juvenile_trend,
  if (ni_juvenile_trend > 0.3) {
    "Still positive -- juveniles still being helped, not hurt, by fishing. The immigration term wasn't what was masking the collapse regime; try lowering resource productivity (kappa/a0) next, not just removing the subsidy."
  } else if (ni_juvenile_trend < -0.3) {
    "Negative -- juvenile biomass now falling as effort rises, unlike Section 3c. Worth running the full coarse+fine Constant scan (Section 3b's own machinery) on p_scan_ni to see whether yield actually turns over this time."
  } else {
    "Weak/no trend either way over just 5 probe points -- inconclusive at this resolution, not evidence the effect is gone. Widen ni_probe_levels before concluding anything."
  }
))

anchovy_immigration <- anchovy_immigration_with_subsidy   # restore Section 0's original global

################################################################################
# Section 6b: resource RATE limitation, not just removing the subsidy
#
# Section 6's i0=0 probe removed the external top-up but left the
# plankton's own regrowth RATE untouched (rr_pp = r0 * w_full^(rho-1),
# r0=10) -- and the correlation there was ambiguous (weakly negative but
# non-monotonic: mean_juvenile dipped at fish_level=20 then partly
# recovered at 50), consistent with a resource that still recovers fast
# enough to look effectively unlimited to the fish regardless of the
# subsidy. This is exactly the lever Day 4/18's own resource_semichemostat
# setups used -- resource_rate scaled by resource_decrease=0.001 -- to make
# a resource genuinely limiting rather than merely present. Same idea,
# applied to r0 here since this model's resource dynamics are the paper's
# own plankton_logistic(), not resource_semichemostat. i0=0 carried over
# from Section 6 rather than re-adding the subsidy on top.
################################################################################

resource_decrease_6b <- 0.001   # Day 4/18's own convention, reused for comparability

p2_slow_resource <- p2_no_immigration
p2_slow_resource$r0 <- p2$r0 * resource_decrease_6b

make_anchovy_fishing_params_slow_resource <- function(p = p2_slow_resource) {
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

p_scan_sr <- make_anchovy_fishing_params_slow_resource()

anchovy_immigration_with_subsidy <- anchovy_immigration   # re-save in case Section 6 didn't run first
anchovy_immigration <- p2_slow_resource$i0 * p_scan_sr@w_full^(-p2_slow_resource$lambda) *
  exp(-6.9 * (p2_slow_resource$lambda - 1))   # still all-zero: i0=0 carried over

cat(sprintf("Slow-resource variant erepro: %.4g (should sit in [0,1])\n",
           species_params(p_scan_sr)$erepro))

# 1. Does it even survive the settle? A 1000x slower-regrowing resource,
# starting from the same power-law abundance and a FULL-cc_pp initial
# plankton state, is a much harsher environment than anything built so far
# today -- worth checking the population doesn't just starve out during the
# unfished settle itself, before asking anything about fishing.
sim_fork_sr <- make_anchovy_fork_sim(p_scan_sr, t_fork)
fork_bp_sr  <- compute_selected_biomass_series(sim_fork_sr, p_scan_sr, t_cut = t_fork - scan_summary_window)
cat(sprintf(
  "Slow-resource variant, unfished, last %.0fyr of the %.0fyr fork: selected biomass min=%.4g max=%.4g.\n",
  scan_summary_window, t_fork, min(fork_bp_sr), max(fork_bp_sr)
))

total_bm_sr_fork <- unname(getBiomass(sim_fork_sr)[, "Anchovy"])
tv_sr_fork       <- as.numeric(dimnames(sim_fork_sr@n)[[1]])
keep_sr_fork     <- tv_sr_fork > t_fork - scan_summary_window
cat(sprintf(
  "Slow-resource variant, unfished, total biomass over the same window: min=%.4g max=%.4g (ratio=%.3g).\n",
  min(total_bm_sr_fork[keep_sr_fork]), max(total_bm_sr_fork[keep_sr_fork]),
  max(total_bm_sr_fork[keep_sr_fork]) / min(total_bm_sr_fork[keep_sr_fork])
))

# 2. Same small Constant-only probe as Section 6, same fish_level grid --
# directly comparable side by side.
sr_probe_levels <- c(0, 5, 10, 20, 50)

run_sr_diag_case <- function(fl) {
  sim <- project(sim_fork_sr, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 progress_bar = FALSE, effort = fl)
  tv          <- as.numeric(dimnames(sim@n)[[1]])
  total_bm    <- unname(getBiomass(sim)[, "Anchovy"])
  selected_bm <- unname(getBiomass(sim, min_w = p_scan_sr@species_params$w_mat)[, "Anchovy"])
  data.frame(t = tv, total = total_bm, selected = selected_bm,
            juvenile = total_bm - selected_bm, fish_level = fl)
}

sr_series_df <- bind_rows(lapply(sr_probe_levels, run_sr_diag_case))
write.csv(sr_series_df, file.path("interesting_plots", "day38_slow_resource_series.csv"), row.names = FALSE)

sr_summary <- sr_series_df %>%
  filter(t > scan_t_cut) %>%
  group_by(fish_level) %>%
  summarise(min_total = min(total), mean_total = mean(total),
           mean_selected = mean(selected), mean_juvenile = mean(juvenile), .groups = "drop")
write.csv(sr_summary, file.path("interesting_plots", "day38_slow_resource_summary.csv"), row.names = FALSE)

sr_series_long <- sr_series_df %>%
  tidyr::pivot_longer(c(total, selected, juvenile), names_to = "compartment", values_to = "biomass")

sr_series_plot <- ggplot(sr_series_long %>% filter(t >= t_fork),
                         aes(x = t, y = pmax(biomass, 1e-12), color = compartment)) +
  geom_line() +
  facet_wrap(~fish_level, scales = "free_y", labeller = label_both) +
  scale_y_log10() +
  labs(x = "Time (years)", y = "Biomass (log scale, floored at 1e-12 for plotting)",
       title = sprintf("Slow-resource variant (i0=0, r0 x %.3g): does a genuinely rate-limited resource change the picture?",
                       resource_decrease_6b),
       subtitle = "Total vs. selected (w>=w_mat, fishable) vs. juvenile (w<w_mat, never fished) biomass") +
  theme_minimal()
sr_series_plot
save_plot(sr_series_plot, "day38_slow_resource_series.png", width = 11, height = 7)

cat("Section 6b (slow-resource probe, fish_level=[0,5,10,20,50]): see day38_slow_resource_series.png / day38_slow_resource_summary.csv.\n")
print(sr_summary)

sr_juvenile_trend <- cor(sr_summary$fish_level, sr_summary$mean_juvenile)
cat(sprintf(
  "Slow-resource variant: correlation of fish_level with mean_juvenile = %.3f. %s\n",
  sr_juvenile_trend,
  if (sr_juvenile_trend > 0.3) {
    "Still positive -- even a genuinely rate-limited resource doesn't break the cultivation effect here."
  } else if (sr_juvenile_trend < -0.3) {
    "Negative, and hopefully more cleanly monotonic than Section 6's i0=0-only probe -- slowing the resource's own regrowth rate looks like the more load-bearing lever. Worth the full coarse+fine Constant scan on p_scan_sr next."
  } else {
    "Weak/no trend at this resolution -- inconclusive, not evidence the effect is gone."
  }
))

anchovy_immigration <- anchovy_immigration_with_subsidy   # restore Section 0's original global

################################################################################
# Section 6c: calibrating resource_decrease -- how far can it go before the
# unfished population dies from the settle+kick itself, before fishing is
# even in the picture?
#
# Section 6b's resource_decrease=0.001 (Day 4/18's own convention) was too
# severe combined with this recipe's own 10^7 post-settle kick: min_total
# was 1e-12 at EVERY fish_level tested, including unfished. That's the
# population dying from the artificial perturbation, not from fishing --
# meaningless for what Section 6b was actually trying to test. Before
# probing fishing effects again, sweep resource_decrease and check ONLY
# whether the UNFISHED population survives the standard settle+kick recipe
# at each value. No fishing in this section at all -- purely calibration,
# so the next fishing probe starts from a value known to be viable.
################################################################################

resource_decrease_sweep <- c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)

check_resource_decrease <- function(rd) {
  p2_rd <- p2_no_immigration
  p2_rd$r0 <- p2$r0 * rd

  params <- setAnchovyModel(p2_rd)
  params <- setAnchovyMort(params, p2_rd)
  params@interaction[] <- 1
  gp                  <- params@gear_params
  gp$sel_func         <- "knife_edge"
  gp$knife_edge_size  <- params@species_params$w_mat
  gp$catchability     <- 1
  gear_params(params) <- gp

  # <<- : plankton_logistic() looks anchovy_immigration up BY NAME in the
  # global environment on every project() call below, same as Sections 0/6/
  # 6b. i0=0 throughout this sweep (p2_no_immigration carries that over),
  # so this always evaluates to all-zero, but computed properly rather than
  # hand-built, same reasoning as Section 6's own comment.
  anchovy_immigration <<- p2_rd$i0 * params@w_full^(-p2_rd$lambda) * exp(-6.9 * (p2_rd$lambda - 1))

  sim         <- make_anchovy_fork_sim(params, t_fork)
  tv          <- as.numeric(dimnames(sim@n)[[1]])
  keep        <- tv > t_fork - scan_summary_window
  total_bm    <- unname(getBiomass(sim)[, "Anchovy"])
  selected_bm <- unname(getBiomass(sim, min_w = params@species_params$w_mat)[, "Anchovy"])

  data.frame(
    resource_decrease = rd,
    min_total = min(total_bm[keep]), mean_total = mean(total_bm[keep]),
    min_selected = min(selected_bm[keep]), mean_selected = mean(selected_bm[keep])
  )
}

anchovy_immigration_with_subsidy <- anchovy_immigration   # save before the sweep overwrites it repeatedly

resource_decrease_df <- bind_rows(lapply(resource_decrease_sweep, check_resource_decrease))
write.csv(resource_decrease_df, file.path("interesting_plots", "day38_resource_decrease_calibration.csv"), row.names = FALSE)

cat(sprintf(
  "Section 6c (resource_decrease calibration, UNFISHED only, last %.0fyr of the %.0fyr fork):\n",
  scan_summary_window, t_fork
))
print(resource_decrease_df)

resource_decrease_plot <- ggplot(resource_decrease_df, aes(x = resource_decrease, y = pmax(mean_total, 1e-15))) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "resource_decrease (log scale)", y = "Mean unfished total biomass (log scale, floored at 1e-15)",
       title = "Where does the unfished population stop surviving the settle+kick?",
       subtitle = "i0=0 throughout -- looking for the range where the resource is genuinely limiting but the population still survives") +
  theme_minimal()
resource_decrease_plot
save_plot(resource_decrease_plot, "day38_resource_decrease_calibration.png", width = 8, height = 5)

# 1e-6 threshold: comfortably above Section 6b's 1e-12 collapse floor, well
# below every genuinely-alive value seen so far today (all O(0.01-0.5)).
viable_rd <- resource_decrease_df$resource_decrease[resource_decrease_df$mean_total > 1e-6]
cat(sprintf(
  "Viable resource_decrease values (mean unfished total biomass > 1e-6): %s. %s\n",
  paste(viable_rd, collapse = ", "),
  if (length(viable_rd) == 0) {
    "None survived at any tested value -- even resource_decrease=1 (no rate change at all) may need checking, or the 1e-6 survival threshold itself needs revisiting."
  } else {
    sprintf("Smallest viable value found: %.4g -- use that (or a little above it, not right at the edge) as the next fishing-effort probe, replacing Section 6b's resource_decrease=0.001.",
           min(viable_rd))
  }
))

anchovy_immigration <- anchovy_immigration_with_subsidy   # restore Section 0's original global

################################################################################
# Section 7: changing the knife-edge size on the ORIGINAL model
#
# Sections 6/6b/6c tried to make the population collapse-capable by making
# the RESOURCE scarce -- awkward here, since the resource is topped up by
# immigration regardless of grazing, and dialing that down risked killing
# the population outright before fishing even entered the picture (Section
# 6b's own resource_decrease=0.001 result: min_total=1e-12 even unfished).
#
# A more direct lever, on the ORIGINAL model (i0=100, r0=10 -- Section 0's
# own p_scan, not the no-immigration/slow-resource variants above): Section
# 3c found the refuge is structural, not resource-driven. The gear is
# knife-edge at w_mat=10 and NEVER touches anything smaller, so no matter
# how hard adults are fished, the juvenile compartment (>95% of biomass
# throughout Section 3c) is untouchable by construction. Lowering
# knife_edge_size below w_mat removes that blind spot directly instead of
# trying to starve the whole system.
#
# "juvenile"/"selected" below are still split at the BIOLOGICAL w_mat=10
# throughout, regardless of where a given knife_edge_size cuts -- comparing
# the same two compartments across variants, not redefining them each time.
################################################################################

knife_edge_sweep <- c(10, 5, 1, 0.1, p2$w_l)   # w_mat (current default), then
                                                # down toward the larval-mortality
                                                # peak size (p2$w_l = 0.03)
ke_probe_levels  <- c(0, 10, 30, 100)

make_anchovy_fishing_params_at <- function(knife_edge_size, p = p2) {
  params <- setAnchovyModel(p)
  params <- setAnchovyMort(params, p)
  params@interaction[] <- 1

  gp                 <- params@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- knife_edge_size
  gp$catchability    <- 1
  gear_params(params) <- gp

  params
}

# Explicit reset, not an assumption that Sections 6/6b/6c ran first (each
# restores its own subsidy at its own end, but this section doesn't depend
# on that order holding).
anchovy_immigration <- p2$i0 * p_scan@w_full^(-p2$lambda) * exp(-6.9 * (p2$lambda - 1))

run_ke_diag_case <- function(knife_edge_size, fl) {
  params <- make_anchovy_fishing_params_at(knife_edge_size)
  params@initial_n[]    <- last_n     # Section 1b's own fork state -- same
  params@initial_n_pp[] <- last_npp   # starting point for every variant here
  sim <- project(params, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 t_start = t_fork, progress_bar = FALSE, effort = fl)
  tv          <- as.numeric(dimnames(sim@n)[[1]])
  total_bm    <- unname(getBiomass(sim)[, "Anchovy"])
  selected_bm <- unname(getBiomass(sim, min_w = p2$w_mat)[, "Anchovy"])
  data.frame(t = tv, total = total_bm, selected = selected_bm,
            juvenile = total_bm - selected_bm,
            knife_edge_size = knife_edge_size, fish_level = fl)
}

ke_series_df <- bind_rows(lapply(knife_edge_sweep, function(ke) {
  bind_rows(lapply(ke_probe_levels, function(fl) run_ke_diag_case(ke, fl)))
}))
write.csv(ke_series_df, file.path("interesting_plots", "day38_knife_edge_series.csv"), row.names = FALSE)

ke_summary <- ke_series_df %>%
  filter(t > scan_t_cut) %>%
  group_by(knife_edge_size, fish_level) %>%
  summarise(min_total = min(total), mean_total = mean(total),
           mean_selected = mean(selected), mean_juvenile = mean(juvenile), .groups = "drop")
write.csv(ke_summary, file.path("interesting_plots", "day38_knife_edge_summary.csv"), row.names = FALSE)

ke_collapse_plot <- ggplot(ke_summary, aes(x = fish_level, y = pmax(mean_total, 1e-12),
                                           color = factor(knife_edge_size))) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  labs(x = "Fishing effort (Constant schedule)", y = "Mean total biomass (log scale)",
       color = "knife_edge_size",
       title = "Does lowering the knife-edge size make the population collapse-capable?",
       subtitle = "Original model (i0=100, r0=10) -- w_mat=10 is Section 3c's own default; smaller values cut progressively into the juvenile refuge") +
  theme_minimal()
ke_collapse_plot
save_plot(ke_collapse_plot, "day38_knife_edge_collapse.png", width = 9, height = 6)

cat("Section 7 (knife-edge sweep, original model): see day38_knife_edge_collapse.png / day38_knife_edge_summary.csv.\n")
print(ke_summary)

ke_trend <- ke_summary %>%
  group_by(knife_edge_size) %>%
  summarise(juvenile_trend = cor(fish_level, mean_juvenile), .groups = "drop") %>%
  arrange(desc(knife_edge_size))
cat("Correlation of fish_level with mean_juvenile, by knife_edge_size (more negative = more collapse-capable, less negative/positive = still a cultivation effect):\n")
print(ke_trend)

################################################################################
# Section 8: rebuilding the threshold rule on the collapse-capable model
#
# Section 7 found knife_edge_size=5 gives the original model a genuine,
# graduated collapse risk that knife_edge_size=10 (w_mat, the default)
# never had: healthy at fish_level=10 (mean_total=0.381 vs. 0.391 unfished),
# visibly stressed at 30 (0.083, ~5x down), extinct by 100 (3.5e-10, 9+
# orders of magnitude down). That's exactly the condition Day 36's
# threshold-rule apparatus needs to say anything meaningful: a Constant
# schedule that CAN push the population over an edge, so a Threshold
# schedule's lower REALISED average effort (it only fishes near cycle
# peaks) has something real to be measured against.
#
# Reuses Section 1-4's own machinery unchanged (thresholdFMort(),
# attach_threshold_rule(), run_anchovy_threshold_scan(), the sanity-check
# pattern) -- only the gear differs. p_scan/sim_fork/last_n/last_npp are
# all swapped to the knife_edge_size=5 build for the duration of this
# section and restored at the end.
#
# sim_fork specifically has to be REBUILT, not just have p_scan swapped:
# run_anchovy_threshold_scan()'s Constant branch calls
# project(sim_fork, ..., effort=fl), which continues using sim_fork@params'
# OWN embedded gear, not whatever the global p_scan currently points to.
# Swapping p_scan alone would leave Constant silently still fishing with
# knife_edge_size=10 while the Threshold branches correctly used 5 --
# an inconsistent comparison. The unfished trajectory itself is identical
# either way (gear only matters once effort>0, confirmed by Section 7's own
# fish_level=0 rows), so this is a cheap rebuild, not new information --
# just needed so sim_fork's embedded params are correct for continuation.
################################################################################

p_scan_original   <- p_scan
sim_fork_original <- sim_fork
last_n_original    <- last_n
last_npp_original  <- last_npp

p_scan   <- make_anchovy_fishing_params_at(5)
sim_fork <- make_anchovy_fork_sim(p_scan, t_fork)
last_n   <- array(sim_fork@n[dim(sim_fork@n)[1], , , drop = FALSE], dim = dim(sim_fork@n)[-1])
last_npp <- sim_fork@n_pp[dim(sim_fork@n_pp)[1], ]

cat(sprintf("knife_edge_size=5 variant erepro: %.4g (should sit in [0,1])\n",
           species_params(p_scan)$erepro))

# 1. Sanity check -- fish_level=50, deliberately in the "at risk" zone
# Section 7 found (healthy at 10, stressed at 30, extinct at 100). Does
# Constant fishing at 50 crash it, while Threshold (peaks) -- which only
# fires near the cycle's own peaks -- survives on a much lower realised
# average effort?
sanity_fish_level_ke5     <- 50
sanity_threshold_frac_ke5 <- 0.5

sim_ref_sanity_ke5   <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                                progress_bar = FALSE, effort = 0)
bp_ref_sanity_ke5    <- compute_selected_biomass_series(sim_ref_sanity_ke5, p_scan, scan_t_cut)
sharpness_sanity_ke5 <- 0.02 * (max(bp_ref_sanity_ke5) - min(bp_ref_sanity_ke5))
threshold_sanity_ke5 <- unname(quantile(bp_ref_sanity_ke5, probs = sanity_threshold_frac_ke5))

run_sanity_case_ke5 <- function(schedule) {
  if (schedule == "Constant") {
    sim  <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                    progress_bar = FALSE, effort = sanity_fish_level_ke5)
    mode <- NA_character_
  } else {
    mode <- if (schedule == "Threshold (peaks)") "above" else "below"
    p <- attach_threshold_rule(p_scan, threshold = threshold_sanity_ke5, fish_level = sanity_fish_level_ke5,
                               background_level = 0, mode = mode, sharpness = sharpness_sanity_ke5)
    p@initial_n[]    <- last_n
    p@initial_n_pp[] <- last_npp
    sim <- project(p, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                   t_start = t_fork, progress_bar = FALSE, effort = 1)
  }
  list(sim = sim, schedule = schedule, mode = mode)
}

sanity_schedules_ke5    <- c("Constant", "Threshold (peaks)", "Threshold (troughs)")
sanity_cases_ke5        <- lapply(sanity_schedules_ke5, run_sanity_case_ke5)
names(sanity_cases_ke5) <- sanity_schedules_ke5

cat(sprintf("Section 8 sanity check (knife_edge_size=5, fish_level=%.2g): total biomass by schedule --\n",
           sanity_fish_level_ke5))
for (nm in names(sanity_cases_ke5)) {
  bm <- getBiomass(sanity_cases_ke5[[nm]]$sim)[, "Anchovy"]
  cat(sprintf("  %-22s min=%.4g max=%.4g mean=%.4g\n", nm, min(bm), max(bm), mean(bm)))
}

sanity_series_df_ke5 <- bind_rows(lapply(sanity_cases_ke5, function(case) {
  tv   <- as.numeric(dimnames(case$sim@n)[[1]])
  keep <- which(tv > t_fork)
  bp   <- compute_selected_biomass_series(case$sim, p_scan, t_fork)
  on_frac <- if (is.na(case$mode)) {
    rep(NA_real_, length(bp))
  } else {
    compute_on_frac_series(case$sim, p_scan, threshold_sanity_ke5, case$mode, sharpness_sanity_ke5, t_fork)
  }
  data.frame(t = tv[keep], selected_biomass = bp, on_frac = on_frac, schedule = case$schedule)
}))
write.csv(sanity_series_df_ke5, file.path("interesting_plots", "day38_ke5_sanity_series.csv"), row.names = FALSE)

sanity_check_plot_ke5 <- ggplot(sanity_series_df_ke5, aes(x = t, y = pmax(selected_biomass, 1e-15))) +
  geom_rect(data = sanity_series_df_ke5 %>% filter(!is.na(on_frac), on_frac > 0.5),
           aes(xmin = t, xmax = t + 0.2, ymin = -Inf, ymax = Inf),
           inherit.aes = FALSE, fill = "tomato", alpha = 0.25) +
  geom_line() +
  geom_hline(yintercept = threshold_sanity_ke5, linetype = "dashed", color = "grey40") +
  facet_wrap(~schedule, ncol = 1, scales = "free_y") +
  scale_y_log10() +
  labs(x = "Time (years)", y = "Selected biomass (log scale, floored at 1e-15)",
       title = "knife_edge_size=5: does threshold-triggered fishing survive where Constant doesn't?",
       subtitle = sprintf("fish_level=%.2g, threshold_frac=%.2g -- dashed = calibration threshold, red shading = fishing 'on'",
                          sanity_fish_level_ke5, sanity_threshold_frac_ke5)) +
  theme_minimal()
sanity_check_plot_ke5
save_plot(sanity_check_plot_ke5, "day38_ke5_sanity_check.png", width = 9, height = 8)

cat(sprintf(
  "Section 8 sanity check: see day38_ke5_sanity_check.png / day38_ke5_sanity_series.csv.\n"
))

# 2. Effort scan -- fish_level rescaled for THIS gear's own collapse zone.
# Section 3's old 1-9 range meant nothing here (Section 7 found this
# model's transition sits around 10-100, not 1-9).
threshold_frac_ke5 <- 0.5
fish_level_seq_ke5 <- c(10, 30, 50, 70, 100)

scan_ke5_df <- run_anchovy_threshold_scan(fish_level_seq = fish_level_seq_ke5,
                                          threshold_frac_seq = threshold_frac_ke5)
write.csv(scan_ke5_df, file.path("interesting_plots", "day38_ke5_effort_scan.csv"), row.names = FALSE)

# rel_amplitude can be NaN/degenerate once a run has actually gone extinct
# (max and min both ~0) -- expected, not a bug, and mean_yield near zero is
# itself the meaningful signal in those rows.
ke5_yield_plot <- ggplot(scan_ke5_df, aes(x = fish_level, y = pmax(mean_yield, 1e-15), color = schedule)) +
  geom_line() +
  geom_point(size = 2) +
  scale_y_log10() +
  labs(x = "Nominal fish_level", y = "Mean yield (log scale, floored at 1e-15)",
       title = "knife_edge_size=5: yield by schedule across the collapse zone",
       subtitle = sprintf("threshold_frac=%.2g -- a schedule crashing to ~0 mean_yield has gone extinct, not just quiet", threshold_frac_ke5)) +
  theme_minimal()
ke5_yield_plot
save_plot(ke5_yield_plot, "day38_ke5_effort_yield.png", width = 9, height = 6)

cat(sprintf(
  "Section 8 effort scan (knife_edge_size=5, fish_level=[%.2g,%.2g], threshold_frac=%.2g): see day38_ke5_effort_yield.png / day38_ke5_effort_scan.csv.\n",
  min(fish_level_seq_ke5), max(fish_level_seq_ke5), threshold_frac_ke5
))
print(scan_ke5_df %>% select(schedule, fish_level, mean_effort, mean_yield, rel_amplitude, effective_window))

p_scan   <- p_scan_original     # restore Section 0's original gear/model
sim_fork <- sim_fork_original
last_n   <- last_n_original
last_npp <- last_npp_original

################################################################################
# Section 9: finding a knife-edge size that collapses within fish_level=[1,9]
#
# Section 7/8's knife_edge_size=5 needs effort up to 30-100 before actually
# collapsing -- well outside the fish_level range every earlier scan this
# project has used (Day 34 Section 3, Day 36 Sections 2-5, Day 38 Section 3
# all worked in roughly [0.2,9]). A model that only breaks at effort=100
# isn't testing anything against those established efforts. Narrower sweep,
# between knife_edge_size=5 (survives to 30+) and 0.03 (dead by 10):
# knife_edge_size x fish_level=[1,3,5,7,9] specifically, to find where the
# collapse actually lands inside the range this project's own fishing
# scans already use. Uses last_n/last_npp (Section 1b's own fork state,
# not Section 8's ke5 one, which was already restored above) since gear
# doesn't affect the unfished trajectory either way.
################################################################################
knife_edge_sweep_fine <- c(1, 0.5, 0.3, 0.1)
fish_level_sweep_fine <- c(1, 3, 5, 7, 9)

run_ke_fine_case <- function(knife_edge_size, fl) {
  params <- make_anchovy_fishing_params_at(knife_edge_size)
  params@initial_n[]    <- last_n
  params@initial_n_pp[] <- last_npp
  sim <- project(params, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                 t_start = t_fork, progress_bar = FALSE, effort = fl)
  tv          <- as.numeric(dimnames(sim@n)[[1]])
  total_bm    <- unname(getBiomass(sim)[, "Anchovy"])
  selected_bm <- unname(getBiomass(sim, min_w = p2$w_mat)[, "Anchovy"])
  data.frame(t = tv, total = total_bm, selected = selected_bm,
            juvenile = total_bm - selected_bm,
            knife_edge_size = knife_edge_size, fish_level = fl)
}

ke_fine_series_df <- bind_rows(lapply(knife_edge_sweep_fine, function(ke) {
  bind_rows(lapply(fish_level_sweep_fine, function(fl) run_ke_fine_case(ke, fl)))
}))
write.csv(ke_fine_series_df, file.path("interesting_plots", "day38_knife_edge_fine_series.csv"), row.names = FALSE)

ke_fine_summary <- ke_fine_series_df %>%
  filter(t > scan_t_cut) %>%
  group_by(knife_edge_size, fish_level) %>%
  summarise(min_total = min(total), mean_total = mean(total),
           mean_selected = mean(selected), mean_juvenile = mean(juvenile), .groups = "drop")
write.csv(ke_fine_summary, file.path("interesting_plots", "day38_knife_edge_fine_summary.csv"), row.names = FALSE)

ke_fine_plot <- ggplot(ke_fine_summary, aes(x = fish_level, y = pmax(mean_total, 1e-15),
                                            color = factor(knife_edge_size))) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  labs(x = "Fishing effort (Constant schedule)", y = "Mean total biomass (log scale)",
       color = "knife_edge_size",
       title = "Does the population collapse within fish_level=[1,9], this project's own established range?",
       subtitle = "Original model (i0=100, r0=10)") +
  theme_minimal()
ke_fine_plot
save_plot(ke_fine_plot, "day38_knife_edge_fine.png", width = 9, height = 6)

cat("Section 9 (fine knife-edge sweep, fish_level=[1,9]): see day38_knife_edge_fine.png / day38_knife_edge_fine_summary.csv.\n")
print(ke_fine_summary)

# 0.01 threshold: ~2.5% of the ~0.39 unfished baseline seen throughout this
# day -- a clear, substantial collapse, not just noise.
ke_fine_collapse_point <- ke_fine_summary %>%
  filter(mean_total < 0.01) %>%
  group_by(knife_edge_size) %>%
  summarise(first_collapse_fish_level = min(fish_level), .groups = "drop") %>%
  arrange(first_collapse_fish_level)
cat("First fish_level at which mean_total drops below 0.01 (~2.5% of the unfished baseline), by knife_edge_size -- smallest is the one to build Section 8's rule comparison on:\n")
print(ke_fine_collapse_point)

################################################################################
# Section 10: size spectra sanity plot
#
# plotSpectra(p_scan) plots p_scan@initial_n directly -- and p_scan's own
# initial_n has never actually been used anywhere in this day's work. It's
# whatever mizer's default community-spectrum initialiser (get_initial_n())
# produced when newMultispeciesParams() first built the model. That default
# guess calls getEReproAndGrowth(), which computes feeding_level as a
# saturating function of h -- as h -> Inf (the paper's own choice for this
# species), feeding_level evaluates to exactly 0 in floating point, and
# intake is then computed as feeding_level * h in a SEPARATE step, so
# 0 * Inf = NaN (an IEEE754 floating-point trap, not a mathematical
# necessity -- the actual limit is finite). That NaN corrupts get_initial_
# n()'s guess. It does NOT affect project()'s own dynamics, which is why
# every actual simulation this day has produced sane, finite numbers --
# every one of them seeds initial_n explicitly via make_anchovy_fork_sim(),
# bypassing get_initial_n() entirely. p_scan@initial_n is the one place
# that never happened.
#
# Plotting sim_fork instead -- the settled, cycling population every other
# section of this day's work is actually built on.
################################################################################

p_scan_spectra_plot <- plotSpectra(sim_fork)
p_scan_spectra_plot
save_plot(p_scan_spectra_plot, "day38_p_scan_spectra.png", width = 8, height = 6)

# Same, at effort=10 -- Section 3c's own finding was that this barely
# touches the original (knife_edge_size=10) model at all (mean_total=0.489
# vs. 0.391 unfished), so this should look close to sim_fork's own spectrum
# above, not degenerate or crashed.
sim_effort10 <- project(sim_fork, t_max = scan_post_fork_years, dt = p2$dt, t_save = 0.2,
                        progress_bar = FALSE, effort = 10)
p_scan_effort10_spectra_plot <- plotSpectra(sim_effort10)
plotSpectra2(sim_effort10, sim_fork)
save_plot(p_scan_effort10_spectra_plot, "day38_p_scan_effort10_spectra.png", width = 8, height = 6)

species_params(sim_effort10)$erepro
species_params(sim_fork)$erepro
(getRDI(sim_fork,time_range = c(10,20)))
getRDI(getParams(sim_fork))
summary(validSim(sim_fork))

