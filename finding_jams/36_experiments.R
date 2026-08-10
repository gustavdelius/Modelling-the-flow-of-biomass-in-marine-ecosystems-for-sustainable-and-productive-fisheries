library(mizer)
library(mizerExperimental)
library(dplyr)
library(tidyverse)
library(ggplot2)

# Day 36 revisits Day 34 Section 4's window x fish_level scan (the grid that
# found Constant fishing dominates Peaks-only/Troughs-only on both amplitude
# and yield at resource_decrease=0.005) and tests three improvements to that
# scan individually, then stacks all three together in one final grid:
#
#   A. Fishing at greater efforts for the schedules that AREN'T fished
#      constantly -- Day 34 Section 4b extended fish_level past 1 for
#      Constant only (up to 4), on the reasoning that it was still improving
#      on both axes at the top of the range originally tested. Peaks-only/
#      Troughs-only never got the same extension -- the window x fish_level
#      grid stopped at fish_level=1 for every schedule. Section 2 below
#      pushes fish_level out to 8 for all three schedules (matching Day 34's
#      own "What's Next" item 3, which named 8 as the next target).
#   B. Changing the window -- Day 34's own window_seq only ever ran 0.5-3 by
#      0.5, "from a genuinely narrow burst up to nearly the whole cycle" by
#      that section's own comment (period ~2.9yr there). Section 3 below
#      pushes both ends: narrower down to 0.1 and wider out to 9 (three full
#      cycles' worth).
#   C. Fishing at lighter intensities away from the peak/trough window --
#      every event-triggered schedule since Day 9 uses a hard 0 outside its
#      window. Section 4 below replaces that hard floor with a light,
#      constant background effort applied everywhere, with the full
#      fish_level only added on top of it near each detected event.
#
# Section 1 generalises Day 34's own run_window_effort_scan() with one new
# argument (background_level_seq, defaulting to 0) rather than writing three
# separate scan functions -- passing the default reproduces Day 34's own
# results exactly. Section 5 stacks all three levers in one combined grid.
#
# Same self-contained convention as every script since Day 20: helper
# functions redefined here rather than sourced from `34_experiments.R`.

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

# Knife-edge gear pinned at w_mat -- every size class at or above maturity
# fully selected, everything below it fully protected. Unchanged from Day
# 18/34.
make_fishing_params <- function(lambda = 2.05, resource_decrease = 0.001) {
  p  <- make_params(lambda = lambda, resource_decrease = resource_decrease)
  gp <- p@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- p@species_params$w_mat
  gp$catchability    <- 1
  gear_params(p)     <- gp
  p
}

# Day 18/20/34's settle+kick recipe, unchanged: 10yr unfished settle, then a
# 1000x knock-down of the w in [10,100] mature range at t=10, then run
# forward under the given effort for the rest of t_total.
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

# Same regime and fork point as Day 34 throughout, so every section below
# compares directly against that day's own numbers. scan_post_fork_years/
# scan_summary_window match Day 34 Section 4's own broad-scan settings
# (shorter than Section 3's careful 300/50 there -- a broad scan needs many
# runs, and each grid point only needs to reach its own settled cycle).
rd_focus             <- 0.005
t_fork               <- 500
scan_post_fork_years <- 100
scan_summary_window  <- 30

################################################################################
# Section 1: generalised scan machinery
#
# Identical to Day 34's own run_window_effort_scan() except for one new
# argument, background_level_seq (default 0): the effort level applied
# OUTSIDE the peak/trough window, instead of Day 34's hard 0 floor. Capped
# per-point at min(background_level, fish_level) inside build_schedule() so
# a background_level passed in above the fish_level being tested at that
# point can never invert the schedule into fishing harder away from the
# window than during it. Swept internally alongside fish_level_seq/
# window_seq (not via repeated top-level calls) for the same efficiency
# reason Day 34 already swept fish_level_seq/window_seq internally: the
# settle+detection step (sim_fork, sim_ref, peak/trough detection) is
# expensive and doesn't depend on any of the three, so it happens ONCE per
# resource_decrease regardless of how many (fish_level, window,
# background_level) combinations get tested against it.
################################################################################

run_window_effort_scan <- function(resource_decrease, fish_level_seq, t_fork,
                                   scan_post_fork_years, scan_summary_window,
                                   window_seq = NULL, background_level_seq = 0,
                                   schedules = c("Peaks only", "Troughs only", "Constant")) {
  scan_t_cut <- t_fork + scan_post_fork_years - scan_summary_window

  p_scan    <- make_fishing_params(resource_decrease = resource_decrease)
  gp_scan   <- gear_params(p_scan)
  gear_name <- gp_scan$gear[gp_scan$species == "Anchovy"][1]

  sim_fork <- make_limit_cycle_sim(p_scan, t_total = t_fork, effort = 0)

  needs_events <- any(c("Peaks only", "Troughs only") %in% schedules) && !is.null(window_seq)
  peak_times   <- numeric(0)
  trough_times <- numeric(0)

  if (needs_events) {
    sim_ref <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                       progress_bar = FALSE, effort = 0, method = "tr_bdf2")

    all_times <- as.numeric(dimnames(sim_ref@n)[[1]])
    keep      <- all_times >= t_fork
    t_p       <- all_times[keep]
    time_lbls <- dimnames(sim_ref@n)[[1]][keep]

    bm <- getBiomass(sim_ref)[, "Anchovy"][time_lbls]
    n  <- length(bm)
    is_peak   <- logical(n)
    is_trough <- logical(n)
    is_peak[2:(n - 1)]   <- bm[2:(n - 1)] > bm[1:(n - 2)] & bm[2:(n - 1)] > bm[3:n]
    is_trough[2:(n - 1)] <- bm[2:(n - 1)] < bm[1:(n - 2)] & bm[2:(n - 1)] < bm[3:n]
    peak_times   <- t_p[is_peak]
    trough_times <- t_p[is_trough]

    # floor_level replaces Day 34's hard 0 outside the window with
    # background_level (capped at fish_level -- see Section 1's own header
    # comment above).
    build_schedule <- function(event_times, window, fish_level, floor_level) {
      ev <- rep(floor_level, length(t_p))
      for (et in event_times) {
        ev[t_p >= et - window / 2 & t_p <= et + window / 2] <- fish_level
      }
      matrix(ev, ncol = 1, dimnames = list(time_lbls, gear_name))
    }
  }

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

  # Constant is unaffected by background_level -- it's already "on" at
  # fish_level everywhere -- so it's swept over fish_level_seq only, same as
  # Day 34.
  if ("Constant" %in% schedules) {
    out$constant <- bind_rows(lapply(fish_level_seq, function(fl) {
      sim <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                    progress_bar = FALSE, effort = fl, method = "tr_bdf2")
      cbind(scan_metrics(sim), window = NA_real_, fish_level = fl, schedule = "Constant",
           mean_effort = fl, background_level = NA_real_)
    }))
  }

  if (needs_events) {
    out$events <- bind_rows(lapply(fish_level_seq, function(fl) {
      bind_rows(lapply(window_seq, function(w) {
        bind_rows(lapply(background_level_seq, function(bg) {
          floor_level <- min(bg, fl)
          rows <- list()
          if ("Peaks only" %in% schedules) {
            effort_p <- build_schedule(peak_times, w, fl, floor_level)
            sim_p <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                             progress_bar = FALSE, effort = effort_p, method = "tr_bdf2")
            rows$peaks <- cbind(scan_metrics(sim_p), window = w, fish_level = fl, schedule = "Peaks only",
                                mean_effort = mean(effort_p), background_level = bg)
          }
          if ("Troughs only" %in% schedules) {
            effort_t <- build_schedule(trough_times, w, fl, floor_level)
            sim_t <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                             progress_bar = FALSE, effort = effort_t, method = "tr_bdf2")
            rows$troughs <- cbind(scan_metrics(sim_t), window = w, fish_level = fl, schedule = "Troughs only",
                                  mean_effort = mean(effort_t), background_level = bg)
          }
          bind_rows(rows)
        }))
      }))
    }))
  }

  bind_rows(out) %>%
    mutate(resource_decrease = resource_decrease, n_peaks = length(peak_times), n_troughs = length(trough_times))
}

################################################################################
# Section 2 (Improvement A): fishing at greater efforts for the schedules
# that AREN'T fished constantly
#
# window_seq held at Day 34's own 0.5-3 by 0.5 so this section isolates the
# effort lever alone; fish_level_seq extended from Day 34's 0.2-1 to 0.2-8
# for every schedule, not just Constant.
################################################################################

window_seq_a     <- seq(0.5, 3, by = 0.5)   # Day 34's own range, unchanged
fish_level_seq_a <- seq(0.2, 8, by = 0.4)   # 0.2-8, 8x Day 34's own upper bound

scan_a_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_a, window_seq = window_seq_a,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window
)

write.csv(scan_a_df, file.path("interesting_plots", "day36_effort_ext_scan.csv"), row.names = FALSE)

# Same amplitude/yield-vs-mean_effort read as Day 34 Section 4d, now over
# the full 0.2-8 range for every schedule -- does Peaks-only/Troughs-only
# ever catch up to Constant at high enough effort, or does the gap just
# widen?
effort_ext_amplitude_plot <- ggplot(scan_a_df, aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_point(alpha = 0.6) +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Relative amplitude of biomass",
       title = "Does Peaks-only/Troughs-only catch up to Constant at higher effort?",
       subtitle = sprintf("resource_decrease=%.4g -- fish_level extended to 8 for every schedule (Day 34 only extended Constant past 1)", rd_focus)) +
  theme_minimal()
effort_ext_amplitude_plot
save_plot(effort_ext_amplitude_plot, "day36_effort_ext_amplitude.png", width = 9, height = 6)

effort_ext_yield_plot <- ggplot(scan_a_df, aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_point(alpha = 0.6) +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Mean yield",
       title = "Yield vs. mean effort, fish_level extended to 8 for every schedule",
       subtitle = sprintf("resource_decrease=%.4g", rd_focus)) +
  theme_minimal()
effort_ext_yield_plot
save_plot(effort_ext_yield_plot, "day36_effort_ext_yield.png", width = 9, height = 6)

cat(sprintf(
  "Section 2 (Improvement A, fish_level=[%.2g,%.2g] for every schedule): Constant's own best rel_amplitude at fish_level=%.2g = %.4g; best Peaks/Troughs rel_amplitude at the same fish_level = %.4g.\n",
  min(fish_level_seq_a), max(fish_level_seq_a), max(fish_level_seq_a),
  min(scan_a_df$rel_amplitude[scan_a_df$schedule == "Constant" & scan_a_df$fish_level == max(fish_level_seq_a)]),
  min(scan_a_df$rel_amplitude[scan_a_df$schedule != "Constant" & scan_a_df$fish_level == max(fish_level_seq_a)])
))

################################################################################
# Section 3 (Improvement B): changing the window
#
# fish_level_seq held at Day 34's own 0.2-1 by 0.2 so this section isolates
# the window lever alone; window_seq widened from Day 34's 0.5-3 to 0.1-9 --
# a much tighter burst at the low end, and wide enough at the high end (~3x
# the detected ~2.9yr period) that the schedule should start behaving like
# Constant, since it's "on" for most or all of the run by then.
################################################################################

window_seq_b     <- c(0.1, 0.25, seq(0.5, 3, by = 0.5), 4.5, 6, 9)
fish_level_seq_b <- seq(0.2, 1, by = 0.2)   # Day 34's own range, unchanged

scan_b_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_b, window_seq = window_seq_b,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window
)

write.csv(scan_b_df, file.path("interesting_plots", "day36_window_ext_scan.csv"), row.names = FALSE)

window_ext_amplitude_plot <- ggplot(scan_b_df %>% filter(schedule != "Constant"),
                                    aes(x = window, y = rel_amplitude, color = schedule)) +
  geom_line(aes(group = interaction(schedule, fish_level), alpha = fish_level), linewidth = 0.8) +
  geom_point(aes(alpha = fish_level), size = 1.5) +
  labs(x = "Window (years either side of event)", y = "Relative amplitude of biomass",
       title = "Oscillation amplitude across a much wider window range",
       subtitle = sprintf("resource_decrease=%.4g -- window now runs 0.1-9 (Day 34 tested 0.5-3); at wide windows the event schedule should start converging toward Constant's own amplitude", rd_focus)) +
  theme_minimal()
window_ext_amplitude_plot
save_plot(window_ext_amplitude_plot, "day36_window_ext_amplitude.png", width = 10, height = 6)

window_ext_yield_plot <- ggplot(scan_b_df %>% filter(schedule != "Constant"),
                                aes(x = window, y = mean_yield, color = schedule)) +
  geom_line(aes(group = interaction(schedule, fish_level), alpha = fish_level), linewidth = 0.8) +
  geom_point(aes(alpha = fish_level), size = 1.5) +
  labs(x = "Window (years either side of event)", y = "Mean yield",
       title = "Yield across a much wider window range",
       subtitle = sprintf("resource_decrease=%.4g", rd_focus)) +
  theme_minimal()
window_ext_yield_plot
save_plot(window_ext_yield_plot, "day36_window_ext_yield.png", width = 10, height = 6)

cat("Section 3 (Improvement B, window=[0.1,9]): rel_amplitude/mean_yield at the widest window tested, by schedule and fish_level:\n")
print(scan_b_df %>% filter(schedule != "Constant", window == max(window_seq_b)) %>%
       select(schedule, fish_level, rel_amplitude, mean_yield))

################################################################################
# Section 4 (Improvement C): fishing at lighter intensities away from the
# peak/trough window
#
# window and fish_level held at Day 34 Section 3's own point (window=3,
# fish_level=0.5) so this section isolates the background lever alone.
# background_level_seq=0 reproduces Day 34's own hard-0-outside-the-window
# schedule exactly, included here as the baseline the rest of the sweep is
# measured against.
################################################################################

window_c          <- 3      # Day 34 Section 3's own value
fish_level_c      <- 0.5    # Day 34 Section 3's own value
background_seq_c <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)   # 0 = Day 34's own schedule

scan_c_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_c, window_seq = window_c,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window,
  background_level_seq = background_seq_c
)

write.csv(scan_c_df, file.path("interesting_plots", "day36_background_scan.csv"), row.names = FALSE)

# Constant is unaffected by background_level (it's already "on" everywhere
# at fish_level_c) -- one flat reference line rather than one point per
# background_level.
constant_ref_c <- scan_c_df %>% filter(schedule == "Constant") %>% slice(1)

background_amplitude_plot <- ggplot(scan_c_df %>% filter(schedule != "Constant"),
                                    aes(x = background_level, y = rel_amplitude, color = schedule)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = constant_ref_c$rel_amplitude, linetype = "dashed", color = "grey40") +
  labs(x = "Background effort (applied outside the window)", y = "Relative amplitude of biomass",
       title = "Does a light background effort ease the oscillation further?",
       subtitle = sprintf("resource_decrease=%.4g, window=%g, fish_level=%.2g -- dashed line = Constant's own amplitude at fish_level=%.2g", rd_focus, window_c, fish_level_c, fish_level_c)) +
  theme_minimal()
background_amplitude_plot
save_plot(background_amplitude_plot, "day36_background_amplitude.png", width = 9, height = 6)

background_yield_plot <- ggplot(scan_c_df %>% filter(schedule != "Constant"),
                                aes(x = background_level, y = mean_yield, color = schedule)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = constant_ref_c$mean_yield, linetype = "dashed", color = "grey40") +
  labs(x = "Background effort (applied outside the window)", y = "Mean yield",
       title = "Does a light background effort raise yield toward Constant's own?",
       subtitle = sprintf("resource_decrease=%.4g, window=%g, fish_level=%.2g -- dashed line = Constant's own yield at fish_level=%.2g", rd_focus, window_c, fish_level_c, fish_level_c)) +
  theme_minimal()
background_yield_plot
save_plot(background_yield_plot, "day36_background_yield.png", width = 9, height = 6)

cat("Section 4 (Improvement C, background_level=[0,0.4] at window=3, fish_level=0.5):\n")
print(scan_c_df %>% filter(schedule != "Constant") %>% select(schedule, background_level, rel_amplitude, mean_yield))

################################################################################
# Section 5: all three improvements stacked
#
# Sections 2-4 tested extended effort, extended window, and a light
# background floor one at a time, each holding the other two at Day 34's
# original settings. This section runs them together -- extended
# fish_level_seq (Section 2's range), extended window_seq (Section 3's
# range), and a light background_level (0.1, the low end of Section 4's own
# range, picked so it stays a genuine "light background" relative to the
# fish_level range being tested here rather than swamping it) -- in one
# combined grid, at the same rd_focus used throughout this file.
#
# fish_level_seq_stacked x window_seq_stacked is 10 x 6 = 60 grid points x 2
# schedules (peaks + troughs) + 10 Constant points = 130 project() calls --
# the most expensive section in this file, hence the coarser step sizes
# below relative to Sections 2-3's own finer sweeps.
################################################################################

fish_level_seq_stacked   <- seq(0.2, 8, length.out = 10)
window_seq_stacked       <- c(0.1, 0.5, 1.5, 3, 6, 9)
background_level_stacked <- 0.1

scan_stacked_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_stacked, window_seq = window_seq_stacked,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window,
  background_level_seq = background_level_stacked
)

write.csv(scan_stacked_df, file.path("interesting_plots", "day36_stacked_scan.csv"), row.names = FALSE)

stacked_tradeoff_plot <- ggplot(scan_stacked_df, aes(x = rel_amplitude, y = mean_yield, color = schedule, size = fish_level)) +
  geom_point(alpha = 0.7) +
  labs(x = "Relative amplitude of biomass (lower = more eased)", y = "Mean yield (higher = better)",
       title = "Amplitude vs. yield trade-off with all three improvements stacked",
       subtitle = sprintf("resource_decrease=%.4g -- fish_level up to 8, window up to 9, background_level=%.2g outside the window; ideal points sit toward the top-left", rd_focus, background_level_stacked)) +
  theme_minimal()
stacked_tradeoff_plot
save_plot(stacked_tradeoff_plot, "day36_stacked_tradeoff.png", width = 9, height = 6)

stacked_amplitude_heatmap <- ggplot(scan_stacked_df %>% filter(schedule != "Constant"),
                                    aes(x = window, y = fish_level, fill = rel_amplitude)) +
  geom_tile() +
  facet_wrap(~schedule) +
  scale_fill_viridis_c() +
  labs(x = "Window (years either side of event)", y = "fish_level (effort)", fill = "Relative\namplitude",
       title = "Oscillation amplitude, all three improvements stacked",
       subtitle = sprintf("resource_decrease=%.4g, background_level=%.2g -- darker = lower amplitude = more eased", rd_focus, background_level_stacked)) +
  theme_minimal()
stacked_amplitude_heatmap
save_plot(stacked_amplitude_heatmap, "day36_stacked_amplitude_heatmap.png", width = 10, height = 5)

# Same direct generalisation check Day 34 Section 4c used: does ANY
# stacked-improvement point beat Constant's own best on amplitude, and
# separately on yield, now that Peaks/Troughs have the same effort range
# Constant already had (Section 2), a much wider window range (Section 3),
# and a light background floor (Section 4) all at once?
constant_best_amplitude_stacked <- min(scan_stacked_df$rel_amplitude[scan_stacked_df$schedule == "Constant"])
constant_best_yield_stacked     <- max(scan_stacked_df$mean_yield[scan_stacked_df$schedule == "Constant"])
event_stacked_df                <- scan_stacked_df %>% filter(schedule != "Constant")

best_stacked_points <- event_stacked_df %>%
  filter(rel_amplitude <= constant_best_amplitude_stacked) %>%
  arrange(desc(mean_yield)) %>%
  head(5)

write.csv(best_stacked_points, file.path("interesting_plots", "day36_stacked_best_points.csv"), row.names = FALSE)

cat(sprintf(
  "Section 5 (all three improvements stacked): Constant's own best rel_amplitude=%.4g, best mean_yield=%.4g. Any stacked-improvement Peaks/Troughs point beats it on amplitude: %s; on yield: %s.\n",
  constant_best_amplitude_stacked, constant_best_yield_stacked,
  any(event_stacked_df$rel_amplitude < constant_best_amplitude_stacked),
  any(event_stacked_df$mean_yield > constant_best_yield_stacked)
))
cat("Top stacked-improvement (schedule, window, fish_level) combinations at or below Constant's own best amplitude, ranked by yield:\n")
print(best_stacked_points)

################################################################################
# Section 6: summary
################################################################################

cat("\n===== Day 36 summary =====\n")
cat(sprintf(
  "Section 2 (Improvement A -- greater effort for Peaks/Troughs, fish_level=[%.2g,%.2g]): see day36_effort_ext_amplitude.png / day36_effort_ext_yield.png / day36_effort_ext_scan.csv.\n",
  min(fish_level_seq_a), max(fish_level_seq_a)
))
cat(sprintf(
  "Section 3 (Improvement B -- wider window range, window=[%.2g,%.2g]): see day36_window_ext_amplitude.png / day36_window_ext_yield.png / day36_window_ext_scan.csv.\n",
  min(window_seq_b), max(window_seq_b)
))
cat(sprintf(
  "Section 4 (Improvement C -- light background effort, background_level=[%.2g,%.2g] at window=%g, fish_level=%.2g): see day36_background_amplitude.png / day36_background_yield.png / day36_background_scan.csv.\n",
  min(background_seq_c), max(background_seq_c), window_c, fish_level_c
))
cat(sprintf(
  "Section 5 (all three stacked, fish_level=[%.2g,%.2g], window=[%.2g,%.2g], background_level=%.2g): see day36_stacked_tradeoff.png / day36_stacked_amplitude_heatmap.png / day36_stacked_best_points.csv.\n",
  min(fish_level_seq_stacked), max(fish_level_seq_stacked),
  min(window_seq_stacked), max(window_seq_stacked), background_level_stacked
))
