library(mizer)
library(mizerExperimental)
library(dplyr)
library(tidyverse)
library(ggplot2)

# Day 36 replaces Day 34 Section 4's open-loop window x fish_level scan with
# a closed-loop one: a custom fishing-mortality (FMort) rate function,
# registered via setRateFunction(), decides every solver step whether to
# fish from the CURRENT state -- on once selected biomass rises above a
# threshold, off once it drops below. No window parameter any more; how
# long fishing stays on each trigger is emergent (effective_window below).
#
# Custom rate functions must be top-level named functions (mizer looks them
# up by name), so thresholdFMort() reads its per-point config off
# other_params(params) rather than from closure state.
#
# Sections 2-5 reuse Day 34's three-improvements structure -- greater effort
# for the non-constant schedules, changing the threshold (was: the window),
# lighter background effort -- then stack all three.
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

# Day 18/34's own anchovy setup, unchanged.
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

# Knife-edge gear pinned at w_mat -- selects everything at/above maturity.
make_fishing_params <- function(lambda = 2.05, resource_decrease = 0.001) {
  p  <- make_params(lambda = lambda, resource_decrease = resource_decrease)
  gp <- p@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- p@species_params$w_mat
  gp$catchability    <- 1
  gear_params(p)     <- gp
  p
}

# Day 18/20/34's settle+kick recipe: 10yr unfished settle, 1000x knock-down
# of w in [10,100] at t=10, then run forward under the given effort.
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

# Same regime/fork point as Day 34 throughout, for direct comparison.
rd_focus             <- 0.005
t_fork               <- 500
scan_post_fork_years <- 100
scan_summary_window  <- 30

################################################################################
# Section 0: the custom FMort rate function
################################################################################

thresholdFMort <- function(params, n, n_pp, n_other, t, effort, e_growth, pred_mort, ...) {
  p <- other_params(params)

  # Selected biomass at reference effort=1 -- fish_level-independent since
  # mizerFMortGear() is linear in effort.
  f_ref            <- mizerFMortGear(params, effort = 1)
  biomass_density  <- sweep(n, 2, params@w * params@dw, "*")
  selected_biomass <- sum(colSums(f_ref) * biomass_density)

  # hard_step=TRUE reverts to a genuine 0/1 switch -- an explicit toggle for
  # comparison against the smooth version, rather than a silent replacement.
  # The smoothing was a hypothesis about why results looked wrong (Newton
  # chatter on a discontinuity in n), never actually confirmed against real
  # data -- the t_start bug (fixed separately) reproduces that exact "NA
  # everywhere" symptom regardless of hard vs. smooth, so it's worth
  # comparing directly now that that's out of the way.
  direction <- if (p$mode == "above") 1 else -1
  on_frac   <- if (isTRUE(p$hard_step)) {
    as.numeric(direction * (selected_biomass - p$threshold) > 0)
  } else {
    plogis(direction * (selected_biomass - p$threshold) / max(p$sharpness, 1e-12))
  }

  # effort is ignored -- decided fresh from state, not set externally.
  # dim/dimnames forced explicitly so the return shape always matches what
  # setRateFunction() validates against, regardless of what the arithmetic
  # above does to it.
  fmort_on  <- colSums(mizerFMortGear(params, effort = p$fish_level))
  fmort_off <- colSums(mizerFMortGear(params, effort = p$background_level))
  result <- on_frac * fmort_on + (1 - on_frac) * fmort_off
  dim(result)      <- dim(n)
  dimnames(result) <- dimnames(n)
  result
}

# Selected-biomass series from a completed sim's saved n, past t_cut --
# used to calibrate threshold_frac and to reconstruct what a fished run
# actually did (threshold_diagnostics()).
compute_selected_biomass_series <- function(sim, params, t_cut) {
  tv    <- as.numeric(dimnames(sim@n)[[1]])
  keep  <- which(tv > t_cut)
  f_ref <- mizerFMortGear(params, effort = 1)

  vapply(keep, function(i) {
    # drop=FALSE + reshape: sim@n[i,,] alone would also drop the (length-1)
    # species dimension for this single-species model.
    n_i <- array(sim@n[i, , , drop = FALSE], dim = dim(sim@n)[-1])
    biomass_density <- sweep(n_i, 2, params@w * params@dw, "*")
    sum(colSums(f_ref) * biomass_density)
  }, numeric(1))
}

# Writes one grid point's config onto params and registers thresholdFMort().
# other_params() must be set before setRateFunction(), which test-calls the
# function immediately. sharpness has no default: it must be scaled off the
# cycle's own biomass range (in run_window_effort_scan()), not off threshold
# itself, which can sit near zero at low threshold_frac.
attach_threshold_rule <- function(params, threshold, fish_level, background_level = 0,
                                  mode = c("above", "below"), sharpness, hard_step = FALSE) {
  mode <- match.arg(mode)
  other_params(params) <- list(threshold = threshold, fish_level = fish_level,
                               background_level = background_level, mode = mode,
                               sharpness = sharpness, hard_step = hard_step)
  setRateFunction(params, "FMort", "thresholdFMort")
}

# on_frac series, replayed with the exact threshold/mode/sharpness/hard_step
# a given sim was actually run with.
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

# mean_effort/effective_window/n_bursts reconstructed after the fact from
# on_frac -- there's no precomputed schedule to read back any more.
threshold_diagnostics <- function(sim, threshold, fish_level, background_level, mode, sharpness, t_cut, hard_step = FALSE) {
  tv      <- as.numeric(dimnames(sim@n)[[1]])
  keep    <- tv > t_cut
  on_frac <- compute_on_frac_series(sim, sim@params, threshold, mode, sharpness, t_cut, hard_step)
  is_on   <- on_frac > 0.5

  runs   <- rle(is_on)
  n_runs <- length(runs$values)

  # The first and last run in the observed window are CENSORED -- we don't
  # know how long an "on" state had already run before the window opened,
  # or how much longer it would have run after it closed. Counting a
  # censored run as if it were a complete burst reports the observation
  # window's own length (scan_summary_window) as the "period" whenever a
  # population gets stuck on one side of the threshold for the whole
  # window -- which looks like an implausible multi-decade cycle that was
  # never actually there. Only interior runs (bounded by a transition on
  # both sides) are real, timeable bursts.
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
# Section 1: the scan
#
# mode="above" brackets peaks (on while biomass is high), mode="below"
# brackets troughs. threshold_frac_seq replaces window_seq, calibrated as a
# QUANTILE of the cycle's own biomass distribution -- not linear min/max
# interpolation, which for a skewed cycle can put threshold_frac=0.5
# somewhere the population barely visits. background_level_seq is the "off"
# level, instead of a hard 0.
################################################################################

run_window_effort_scan <- function(resource_decrease, fish_level_seq, t_fork,
                                   scan_post_fork_years, scan_summary_window,
                                   threshold_frac_seq = NULL, background_level_seq = 0,
                                   schedules = c("Threshold (peaks)", "Threshold (troughs)", "Constant"),
                                   hard_step = FALSE) {
  scan_t_cut <- t_fork + scan_post_fork_years - scan_summary_window

  p_scan <- make_fishing_params(resource_decrease = resource_decrease)

  sim_fork <- make_limit_cycle_sim(p_scan, t_total = t_fork, effort = 0)
  last     <- dim(sim_fork@n)[1]
  last_n   <- array(sim_fork@n[last, , , drop = FALSE], dim = dim(sim_fork@n)[-1])
  last_npp <- sim_fork@n_pp[last, ]

  needs_threshold <- any(c("Threshold (peaks)", "Threshold (troughs)") %in% schedules) &&
    !is.null(threshold_frac_seq)

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
      sim <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                    progress_bar = FALSE, effort = fl, method = "tr_bdf2")
      cbind(scan_metrics(sim), threshold_frac = NA_real_, threshold_bp = NA_real_,
           fish_level = fl, schedule = "Constant",
           mean_effort = fl, background_level = NA_real_,
           effective_window = NA_real_, n_bursts = NA_integer_)
    }))
  }

  if (needs_threshold) {
    # Unfished continuation -- calibrates threshold_frac; no detection step
    # needed since the rule reacts to state directly.
    sim_ref_unfished <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                                progress_bar = FALSE, effort = 0, method = "tr_bdf2")
    bp_ref <- compute_selected_biomass_series(sim_ref_unfished, p_scan, scan_t_cut)

    # Scaled off the cycle's own swing, not off threshold itself -- see
    # attach_threshold_rule()'s comment.
    sharpness_bp <- 0.02 * (max(bp_ref) - min(bp_ref))

    out$threshold <- bind_rows(lapply(fish_level_seq, function(fl) {
      bind_rows(lapply(threshold_frac_seq, function(frac) {
        # Quantile, not linear interpolation -- frac=0.5 is the level the
        # cycle is above exactly half the time.
        threshold_bp <- unname(quantile(bp_ref, probs = frac))
        bind_rows(lapply(background_level_seq, function(bg) {
          rows <- list()
          if ("Threshold (peaks)" %in% schedules) {
            p_peaks <- attach_threshold_rule(p_scan, threshold = threshold_bp, fish_level = fl,
                                             background_level = bg, mode = "above",
                                             sharpness = sharpness_bp, hard_step = hard_step)
            p_peaks@initial_n[]    <- last_n
            p_peaks@initial_n_pp[] <- last_npp
            # effort=1 is a placeholder -- thresholdFMort() ignores it.
            # t_start=t_fork: project() on a MizerParams object (as opposed
            # to continuing a MizerSim) restarts its time axis at 0 by
            # default, which would silently break every t > scan_t_cut
            # filter downstream (scan_t_cut is an absolute time, assuming
            # continuity from t_fork).
            sim_p <- project(p_peaks, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                             t_start = t_fork, progress_bar = FALSE, effort = 1, method = "tr_bdf2")
            diag_p <- threshold_diagnostics(sim_p, threshold_bp, fl, bg, "above",
                                            other_params(p_peaks)$sharpness, scan_t_cut, hard_step)
            rows$peaks <- cbind(scan_metrics(sim_p), threshold_frac = frac, threshold_bp = threshold_bp,
                                fish_level = fl, schedule = "Threshold (peaks)", mean_effort = diag_p$mean_effort,
                                background_level = bg, effective_window = diag_p$effective_window,
                                n_bursts = diag_p$n_bursts)
          }
          if ("Threshold (troughs)" %in% schedules) {
            p_troughs <- attach_threshold_rule(p_scan, threshold = threshold_bp, fish_level = fl,
                                               background_level = bg, mode = "below",
                                               sharpness = sharpness_bp, hard_step = hard_step)
            p_troughs@initial_n[]    <- last_n
            p_troughs@initial_n_pp[] <- last_npp
            sim_t <- project(p_troughs, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                             t_start = t_fork, progress_bar = FALSE, effort = 1, method = "tr_bdf2")
            diag_t <- threshold_diagnostics(sim_t, threshold_bp, fl, bg, "below",
                                            other_params(p_troughs)$sharpness, scan_t_cut, hard_step)
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

  bind_rows(out) %>% mutate(resource_decrease = resource_decrease)
}

################################################################################
# Section 2 (Improvement A): greater effort for the non-constant schedules
#
# threshold_frac held at the cycle's own median (isolates the effort lever);
# fish_level_seq extended to 0.2-8 for every schedule, not just Constant.
################################################################################

threshold_frac_a <- 0.5                    # median of the cycle's own biomass distribution
fish_level_seq_a <- seq(0.2, 8, by = 0.4)  # 0.2-8, 8x Day 34's own upper bound

scan_a_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_a, threshold_frac_seq = threshold_frac_a,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window
)

write.csv(scan_a_df, file.path("interesting_plots", "day36_effort_ext_scan.csv"), row.names = FALSE)

# Does Threshold (peaks)/(troughs) catch up to Constant at higher effort?
effort_ext_amplitude_plot <- ggplot(scan_a_df, aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_point(alpha = 0.6) +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Relative amplitude of biomass",
       title = "Does threshold-triggered fishing catch up to Constant at higher effort?",
       subtitle = sprintf("resource_decrease=%.4g, threshold_frac=%.2g -- fish_level extended to 8 for every schedule (Day 34 only extended Constant past 1)", rd_focus, threshold_frac_a)) +
  theme_minimal()
effort_ext_amplitude_plot
save_plot(effort_ext_amplitude_plot, "day36_effort_ext_amplitude.png", width = 9, height = 6)

effort_ext_yield_plot <- ggplot(scan_a_df, aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_point(alpha = 0.6) +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Mean yield",
       title = "Yield vs. mean effort, fish_level extended to 8 for every schedule",
       subtitle = sprintf("resource_decrease=%.4g, threshold_frac=%.2g", rd_focus, threshold_frac_a)) +
  theme_minimal()
effort_ext_yield_plot
save_plot(effort_ext_yield_plot, "day36_effort_ext_yield.png", width = 9, height = 6)

cat(sprintf(
  "Section 2 (Improvement A, fish_level=[%.2g,%.2g] for every schedule): Constant's own best rel_amplitude at fish_level=%.2g = %.4g; best threshold-triggered rel_amplitude at the same fish_level = %.4g.\n",
  min(fish_level_seq_a), max(fish_level_seq_a), max(fish_level_seq_a),
  min(scan_a_df$rel_amplitude[scan_a_df$schedule == "Constant" & scan_a_df$fish_level == max(fish_level_seq_a)]),
  min(scan_a_df$rel_amplitude[scan_a_df$schedule != "Constant" & scan_a_df$fish_level == max(fish_level_seq_a)])
))

################################################################################
# Section 2b: hard step vs. smooth, same grid as Section 2
#
# The smoothing (thresholdFMort()'s on_frac) was a hypothesis about why
# results looked wrong, never confirmed against real data -- the actual bug
# (t_start) would have produced the same "NA everywhere" symptom regardless
# of hard vs. smooth. Now that it's fixed, rerun Section 2's exact grid with
# hard_step=TRUE and compare directly rather than guessing which one is
# right.
################################################################################

scan_a_hard_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_a, threshold_frac_seq = threshold_frac_a,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window,
  hard_step = TRUE
)

scan_a_compare_df <- bind_rows(
  scan_a_df %>% filter(schedule != "Constant") %>% mutate(step_type = "Smooth"),
  scan_a_hard_df %>% filter(schedule != "Constant") %>% mutate(step_type = "Hard")
)

write.csv(scan_a_compare_df, file.path("interesting_plots", "day36_hard_vs_smooth_scan.csv"), row.names = FALSE)

hard_vs_smooth_amplitude_plot <- ggplot(scan_a_compare_df, aes(x = mean_effort, y = rel_amplitude, color = schedule, shape = step_type)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "Mean effort actually applied over the post-fork run", y = "Relative amplitude of biomass",
       title = "Hard step vs. smooth logistic switch, same grid as Section 2",
       subtitle = sprintf("resource_decrease=%.4g, threshold_frac=%.2g -- does the hard step actually behave differently now that t_start is fixed?", rd_focus, threshold_frac_a)) +
  theme_minimal()
hard_vs_smooth_amplitude_plot
save_plot(hard_vs_smooth_amplitude_plot, "day36_hard_vs_smooth_amplitude.png", width = 9, height = 6)

cat("Section 2b (hard step vs. smooth): rel_amplitude/mean_effort by schedule and step_type, at the highest fish_level tested:\n")
print(scan_a_compare_df %>% filter(fish_level == max(fish_level_seq_a)) %>%
       select(schedule, step_type, fish_level, mean_effort, rel_amplitude, mean_yield, n_bursts))

################################################################################
# Section 3 (Improvement B): changing the threshold (was: the window)
#
# fish_level_seq held at Day 34's own 0.2-1 range; threshold_frac_seq spans
# nearly the full 0-1 range. effective_window (reconstructed) reports what
# duration actually emerged at each threshold_frac.
################################################################################

threshold_frac_seq_b <- seq(0.05, 0.95, by = 0.1)
fish_level_seq_b     <- seq(0.2, 1, by = 0.2)   # Day 34's own range, unchanged

scan_b_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_b, threshold_frac_seq = threshold_frac_seq_b,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window
)

write.csv(scan_b_df, file.path("interesting_plots", "day36_threshold_ext_scan.csv"), row.names = FALSE)

# x = effective_window (years), not threshold_frac -- a quantile is an
# abstract statistical position with no intuitive scale; effective_window is
# the same "years the schedule stays on" quantity Day 34's own window_seq
# used, just reconstructed as an output instead of set as an input. Ties at
# effective_window=0 (threshold_frac values that never triggered at all)
# stack at the left edge, which is itself informative -- it's exactly the
# saturation boundary.
threshold_ext_amplitude_plot <- ggplot(scan_b_df %>% filter(schedule != "Constant"),
                                       aes(x = effective_window, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = interaction(schedule, fish_level), alpha = fish_level), linewidth = 0.8) +
  geom_point(aes(alpha = fish_level), size = 1.5) +
  labs(x = "Effective window (years per on-burst, reconstructed)",
       y = "Relative amplitude of biomass",
       title = "Oscillation amplitude vs. the implicit window length",
       subtitle = sprintf("resource_decrease=%.4g -- window is no longer set directly, this is what threshold_frac actually produced; see day36_threshold_ext_window.png for the frac-to-window mapping itself", rd_focus)) +
  theme_minimal()
threshold_ext_amplitude_plot
save_plot(threshold_ext_amplitude_plot, "day36_threshold_ext_amplitude.png", width = 10, height = 6)

threshold_ext_yield_plot <- ggplot(scan_b_df %>% filter(schedule != "Constant"),
                                   aes(x = effective_window, y = mean_yield, color = schedule)) +
  geom_line(aes(group = interaction(schedule, fish_level), alpha = fish_level), linewidth = 0.8) +
  geom_point(aes(alpha = fish_level), size = 1.5) +
  labs(x = "Effective window (years per on-burst, reconstructed)",
       y = "Mean yield",
       title = "Yield vs. the implicit window length",
       subtitle = sprintf("resource_decrease=%.4g", rd_focus)) +
  theme_minimal()
threshold_ext_yield_plot
save_plot(threshold_ext_yield_plot, "day36_threshold_ext_yield.png", width = 10, height = 6)

# The frac-to-window mapping itself -- threshold_frac stays on the x-axis
# here on purpose, since this is the one plot whose job is to show what
# window length a given threshold_frac setting actually produces.
threshold_ext_window_plot <- ggplot(scan_b_df %>% filter(schedule != "Constant"),
                                    aes(x = threshold_frac, y = effective_window, color = schedule)) +
  geom_line(aes(group = interaction(schedule, fish_level), alpha = fish_level), linewidth = 0.8) +
  geom_point(aes(alpha = fish_level), size = 1.5) +
  labs(x = "threshold_frac", y = "Effective window (years per on-burst, reconstructed)",
       title = "The implicit window: how long fishing stays on, as a function of threshold_frac",
       subtitle = sprintf("resource_decrease=%.4g -- no window_seq exists any more; this is what emerges from the threshold rule itself", rd_focus)) +
  theme_minimal()
threshold_ext_window_plot
save_plot(threshold_ext_window_plot, "day36_threshold_ext_window.png", width = 10, height = 6)

cat("Section 3 (Improvement B, threshold_frac=[0.05,0.95]): rel_amplitude/mean_yield/effective_window at the most extreme threshold_frac tested, by schedule and fish_level:\n")
print(scan_b_df %>% filter(schedule != "Constant", threshold_frac %in% range(threshold_frac_seq_b)) %>%
       select(schedule, threshold_frac, fish_level, rel_amplitude, mean_yield, effective_window, n_bursts))

################################################################################
# Section 4 (Improvement C): lighter background effort away from the trigger
#
# threshold_frac/fish_level held near Day 34 Section 3's own point; sweeps
# background_level instead of a hard-0 "off" state.
################################################################################

threshold_frac_c <- 0.5     # median of the cycle's own biomass distribution
fish_level_c     <- 0.5     # Day 34 Section 3's own value
background_seq_c <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)   # 0 = hard-off baseline

scan_c_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_c, threshold_frac_seq = threshold_frac_c,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window,
  background_level_seq = background_seq_c
)

write.csv(scan_c_df, file.path("interesting_plots", "day36_background_scan.csv"), row.names = FALSE)

constant_ref_c <- scan_c_df %>% filter(schedule == "Constant") %>% slice(1)

background_amplitude_plot <- ggplot(scan_c_df %>% filter(schedule != "Constant"),
                                    aes(x = background_level, y = rel_amplitude, color = schedule)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = constant_ref_c$rel_amplitude, linetype = "dashed", color = "grey40") +
  labs(x = "Background effort (applied while the threshold rule is 'off')", y = "Relative amplitude of biomass",
       title = "Does a light background effort ease the oscillation further?",
       subtitle = sprintf("resource_decrease=%.4g, threshold_frac=%.2g, fish_level=%.2g -- dashed line = Constant's own amplitude at fish_level=%.2g", rd_focus, threshold_frac_c, fish_level_c, fish_level_c)) +
  theme_minimal()
background_amplitude_plot
save_plot(background_amplitude_plot, "day36_background_amplitude.png", width = 9, height = 6)

background_yield_plot <- ggplot(scan_c_df %>% filter(schedule != "Constant"),
                                aes(x = background_level, y = mean_yield, color = schedule)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = constant_ref_c$mean_yield, linetype = "dashed", color = "grey40") +
  labs(x = "Background effort (applied while the threshold rule is 'off')", y = "Mean yield",
       title = "Does a light background effort raise yield toward Constant's own?",
       subtitle = sprintf("resource_decrease=%.4g, threshold_frac=%.2g, fish_level=%.2g -- dashed line = Constant's own yield at fish_level=%.2g", rd_focus, threshold_frac_c, fish_level_c, fish_level_c)) +
  theme_minimal()
background_yield_plot
save_plot(background_yield_plot, "day36_background_yield.png", width = 9, height = 6)

cat("Section 4 (Improvement C, background_level=[0,0.4] at threshold_frac=0.5, fish_level=0.5):\n")
print(scan_c_df %>% filter(schedule != "Constant") %>%
       select(schedule, background_level, rel_amplitude, mean_yield, effective_window))

################################################################################
# Section 5: all three improvements stacked
#
# Extended fish_level_seq + wide threshold_frac_seq + light background_level,
# combined in one grid. Coarser steps than Sections 2-3 -- 130 project() calls.
################################################################################

fish_level_seq_stacked      <- seq(0.2, 8, length.out = 10)
threshold_frac_seq_stacked  <- seq(0.05, 0.95, length.out = 6)
background_level_stacked    <- 0.1

scan_stacked_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_stacked, threshold_frac_seq = threshold_frac_seq_stacked,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window,
  background_level_seq = background_level_stacked
)

write.csv(scan_stacked_df, file.path("interesting_plots", "day36_stacked_scan.csv"), row.names = FALSE)

stacked_tradeoff_plot <- ggplot(scan_stacked_df, aes(x = rel_amplitude, y = mean_yield, color = schedule, size = fish_level)) +
  geom_point(alpha = 0.7) +
  labs(x = "Relative amplitude of biomass (lower = more eased)", y = "Mean yield (higher = better)",
       title = "Amplitude vs. yield trade-off with all three improvements stacked",
       subtitle = sprintf("resource_decrease=%.4g -- fish_level up to 8, threshold_frac spans [0.05,0.95], background_level=%.2g while 'off'; ideal points sit toward the top-left", rd_focus, background_level_stacked)) +
  theme_minimal()
stacked_tradeoff_plot
save_plot(stacked_tradeoff_plot, "day36_stacked_tradeoff.png", width = 9, height = 6)

# x stays threshold_frac here (not effective_window like Section 3's
# amplitude/yield plots) -- geom_tile() needs the actual regular grid
# variable to render clean tiles; effective_window is emergent/continuous
# and would produce uneven, overlapping tiles.
stacked_amplitude_heatmap <- ggplot(scan_stacked_df %>% filter(schedule != "Constant"),
                                    aes(x = threshold_frac, y = fish_level, fill = rel_amplitude)) +
  geom_tile() +
  facet_wrap(~schedule) +
  scale_fill_viridis_c() +
  labs(x = "threshold_frac (quantile of the cycle's own biomass distribution)", y = "fish_level (effort)", fill = "Relative\namplitude",
       title = "Oscillation amplitude, all three improvements stacked",
       subtitle = sprintf("resource_decrease=%.4g, background_level=%.2g -- darker = lower amplitude = more eased", rd_focus, background_level_stacked)) +
  theme_minimal()
stacked_amplitude_heatmap
save_plot(stacked_amplitude_heatmap, "day36_stacked_amplitude_heatmap.png", width = 10, height = 5)

# Does any stacked-improvement point beat Constant's own best, on either axis?
constant_best_amplitude_stacked <- min(scan_stacked_df$rel_amplitude[scan_stacked_df$schedule == "Constant"])
constant_best_yield_stacked     <- max(scan_stacked_df$mean_yield[scan_stacked_df$schedule == "Constant"])
event_stacked_df                <- scan_stacked_df %>% filter(schedule != "Constant")

best_stacked_points <- event_stacked_df %>%
  filter(rel_amplitude <= constant_best_amplitude_stacked) %>%
  arrange(desc(mean_yield)) %>%
  head(5)

write.csv(best_stacked_points, file.path("interesting_plots", "day36_stacked_best_points.csv"), row.names = FALSE)

cat(sprintf(
  "Section 5 (all three improvements stacked): Constant's own best rel_amplitude=%.4g, best mean_yield=%.4g. Any stacked-improvement threshold-triggered point beats it on amplitude: %s; on yield: %s.\n",
  constant_best_amplitude_stacked, constant_best_yield_stacked,
  any(event_stacked_df$rel_amplitude < constant_best_amplitude_stacked),
  any(event_stacked_df$mean_yield > constant_best_yield_stacked)
))
cat("Top stacked-improvement (schedule, threshold_frac, fish_level, effective_window) combinations at or below Constant's own best amplitude, ranked by yield:\n")
print(best_stacked_points %>% select(schedule, threshold_frac, fish_level, effective_window, rel_amplitude, mean_yield))

################################################################################
# Section 6: summary
################################################################################

cat("\n===== Day 36 summary =====\n")
cat(sprintf(
  "Section 2 (Improvement A -- greater effort for threshold-triggered schedules, fish_level=[%.2g,%.2g], threshold_frac=%.2g): see day36_effort_ext_amplitude.png / day36_effort_ext_yield.png / day36_effort_ext_scan.csv.\n",
  min(fish_level_seq_a), max(fish_level_seq_a), threshold_frac_a
))
cat(sprintf(
  "Section 3 (Improvement B -- wider threshold range, threshold_frac=[%.2g,%.2g]): see day36_threshold_ext_amplitude.png / day36_threshold_ext_yield.png / day36_threshold_ext_window.png / day36_threshold_ext_scan.csv.\n",
  min(threshold_frac_seq_b), max(threshold_frac_seq_b)
))
cat(sprintf(
  "Section 4 (Improvement C -- light background effort, background_level=[%.2g,%.2g] at threshold_frac=%.2g, fish_level=%.2g): see day36_background_amplitude.png / day36_background_yield.png / day36_background_scan.csv.\n",
  min(background_seq_c), max(background_seq_c), threshold_frac_c, fish_level_c
))
cat(sprintf(
  "Section 5 (all three stacked, fish_level=[%.2g,%.2g], threshold_frac=[%.2g,%.2g], background_level=%.2g): see day36_stacked_tradeoff.png / day36_stacked_amplitude_heatmap.png / day36_stacked_best_points.csv.\n",
  min(fish_level_seq_stacked), max(fish_level_seq_stacked),
  min(threshold_frac_seq_stacked), max(threshold_frac_seq_stacked), background_level_stacked
))
