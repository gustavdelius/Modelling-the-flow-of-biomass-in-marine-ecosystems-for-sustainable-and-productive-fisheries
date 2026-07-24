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

################################################################################
# Day 34: resource_decrease bifurcation sweep -- max/min biomass (Section 1)
# and max/min yield (Section 2) plotted separately, same convention as Day
# 20/32/33's own bifurcation diagrams: branches collapsing onto one curve
# means a fixed point, fanning apart means a real limit cycle, read directly
# off the plot rather than off a single collapsed amplitude number.
#
# Section 1 re-runs Day 20's own resource_decrease sweep (0.0001-0.5,
# log-spaced, 40 points, forward + backward, warm-started, tr_bdf2,
# effort=0) unchanged. Section 2 reuses Day 18's own knife-edge-at-w_mat
# fishing setup (make_fishing_params(), effort=0.5) -- knife_edge selects
# every size AT OR ABOVE knife_edge_size with no upper cutoff, so pinning
# knife_edge_size at w_mat selects every mature size class and nothing
# below it, unlike Day 33's cod-specific double_sigmoid_length() dome
# (narrow window) or knife_edge(w_inf) (largest sizes only). Both sections
# share one generic sweep function, parametrised by the params-builder and
# the metric extracted from each completed sim.
################################################################################

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

# Knife-edge gear pinned at w_mat -- Day 18's make_fishing_params(),
# unchanged: knife_edge is a hard 0/1 step at knife_edge_size (not a dome
# like Day 33's double_sigmoid_length()), so setting knife_edge_size to
# w_mat means every size class at or above maturity is fully selected and
# everything below it is fully protected.
make_fishing_params <- function(lambda = 2.05, resource_decrease = 0.001) {
  p  <- make_params(lambda = lambda, resource_decrease = resource_decrease)
  gp <- p@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- p@species_params$w_mat
  gp$catchability    <- 1
  gear_params(p)     <- gp
  p
}

# Sanity check: selectivity should be a hard 0/1 step at w_mat -- 0 just
# below it, 1 just above it, and still 1 at the largest size in the model
# (no upper cutoff, unlike a dome) -- confirmed explicitly rather than
# assumed, the same defensive habit Day 32/33 applied to their own gear
# changes.
p_fishing_probe    <- make_fishing_params()
gp_probe           <- gear_params(p_fishing_probe)
gear_name_anchovy  <- gp_probe$gear[gp_probe$species == "Anchovy"][1]
sel_probe          <- getSelectivity(p_fishing_probe)[gear_name_anchovy, "Anchovy", ]
w_mat_probe        <- p_fishing_probe@species_params$w_mat[1]
idx_below          <- which.min(abs(p_fishing_probe@w - w_mat_probe * 0.9))
idx_above          <- which.min(abs(p_fishing_probe@w - w_mat_probe * 1.1))
cat(sprintf(
  "Gear sanity check: selectivity at 90%% of w_mat (w=%.4g g) = %.4g (should be 0); at 110%% of w_mat (w=%.4g g) = %.4g (should be 1); at the largest size in the model = %.4g (should be 1, unlike a dome).\n",
  p_fishing_probe@w[idx_below], sel_probe[idx_below],
  p_fishing_probe@w[idx_above], sel_probe[idx_above],
  tail(sel_probe, 1)
))

# Generic forward/backward resource_decrease sweep -- params_fn builds the
# params object for a given resource_decrease (make_params for the biomass
# side, make_fishing_params for the yield side), metric_fn extracts the
# per-timestep metric from a completed sim. Warm-started throughout: every
# point after the first inherits the previous point's own final (n, n_pp)
# state rather than re-cold-starting, same convention as every other
# bifurcation sweep in this project. tryCatch records a crash as an NA
# point rather than killing the whole sweep.
run_bifurcation_sweep_rd <- function(rd_seq, params_fn, effort, metric_fn,
                                     t_run = 300, lambda = 2.05) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      p <- params_fn(lambda = lambda, resource_decrease = seq_vals[i])
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }

      result <- tryCatch({
        sim  <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                        progress_bar = FALSE, effort = effort, method = "tr_bdf2")
        mv   <- metric_fn(sim)
        tv   <- as.numeric(names(mv))
        # Late window only, guards against the leading transient.
        late <- mv[tv > t_run * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late), error = NA_character_,
            n = sim@n[last, , ], npp = sim@n_pp[last, ])
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_, error = conditionMessage(e),
            n = state_n, npp = state_npp)
      })

      if (!is.na(result$error)) {
        warning(sprintf("resource_decrease=%.5g: %s", seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(rd_seq)
  bwd    <- run_one_direction(rev(rd_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

# Log-spaced x-axis -- resource_decrease spans 0.0001-0.5, 3.7 orders of
# magnitude, same range and resolution as Day 20's own sweep so the two
# are directly comparable.
plot_bifurcation_rd <- function(df, y_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_x_log10() +
    labs(x = "resource_decrease", y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}

# Log-spaced, 0.0001 to 0.5 -- same range Day 20 used (40 points / 3.7
# decades =~ 11 points per decade, rather than thinning out badly at the
# low end the way a linear seq() would over this much range).
rd_seq_bif <- exp(seq(log(0.0001), log(0.5), length.out = 40))

################################################################################
# Section 1: biomass bifurcation -- unfished (effort=0), so this is pure
# population dynamics under a shrinking resource, not a fished-yield
# question.
################################################################################

bif_biomass_df <- run_bifurcation_sweep_rd(
  rd_seq_bif, params_fn = make_params, effort = 0, t_run = 300,
  metric_fn = function(sim) getBiomass(sim)[, "Anchovy"]
)

write.csv(bif_biomass_df, file.path("interesting_plots", "day34_rd_biomass_bifurcation.csv"),
          row.names = FALSE)

bif_biomass_plot <- plot_bifurcation_rd(
  bif_biomass_df, "Biomass",
  "Anchovy bifurcation diagram: biomass vs. resource_decrease",
  "effort=0 (unfished), swept forward and backward, log-spaced 0.0001-0.5 -- branches collapsing onto one curve = fixed point, fanning apart = limit cycle"
)
bif_biomass_plot

save_plot(bif_biomass_plot, "day34_rd_biomass_bif.png")

################################################################################
# Section 2: yield bifurcation -- fished (effort=0.5, Day 18's own value),
# knife-edge gear pinned at w_mat so every mature size class is fully
# selected and nothing below maturity is caught.
################################################################################

bif_yield_df <- run_bifurcation_sweep_rd(
  rd_seq_bif, params_fn = make_fishing_params, effort = 0.5, t_run = 300,
  metric_fn = function(sim) getYield(sim)[, "Anchovy"]
)

write.csv(bif_yield_df, file.path("interesting_plots", "day34_rd_yield_bifurcation.csv"),
          row.names = FALSE)

bif_yield_plot <- plot_bifurcation_rd(
  bif_yield_df, "Yield",
  "Anchovy bifurcation diagram: yield vs. resource_decrease, knife-edge gear at w_mat",
  "Knife-edge selects every size at or above maturity (Day 18's make_fishing_params()), effort=0.5 fixed, swept forward and backward, log-spaced 0.0001-0.5 -- branches collapsing onto one curve = fixed point, fanning apart = limit cycle"
)
bif_yield_plot

save_plot(bif_yield_plot, "day34_rd_yield_bif.png")

################################################################################
# Section 3: fishing only at peaks vs. only at troughs, one fixed
# resource_decrease
#
# Sections 1-2 sweep resource_decrease broadly; this section instead fixes
# resource_decrease at a single value inside the oscillating regime and
# asks a different question: does *when* you fish within a cycle matter,
# not just *how much*. attempt_9/10/11 (Days 9-10) already built a
# peaks-only effort schedule against a hand-rolled sub-adult gear and found
# it worsened the juvenile bottleneck rather than relieving it -- but never
# built the mirror-image troughs-only schedule to compare against, and used
# an older, less-accurate predictor-corrector/box-selectivity setup rather
# than this file's own knife-edge-at-w_mat gear. This section redoes that
# comparison with both schedules, on this file's own conventions.
################################################################################

# Day 18/20's settle+kick recipe, unchanged: 10yr unfished settle, then a
# 1000x knock-down of the w in [10,100] mature range at t=10, then run
# forward under the given effort for the rest of t_total. Not needed until
# now because Sections 1-2's sweeps reach the oscillating regime gradually,
# warm-starting each resource_decrease point off the previous one, rather
# than starting from a single hard kick.
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

# effort=0.5 matches Section 2's own fixed yield-sweep effort, so this
# section's results sit on the same scale; window=3 (years either side of a
# detected event) matches attempt_9-11's original peaks-only branch.
window     <- 3
fish_level <- 0.5

# t_fork raised from 300 to 500: plotting the unfished reference trajectory
# showed the oscillation doesn't actually become regular until roughly
# t=400 (a longer settling transient than Day 18-20's original 300yr fork
# point assumed) -- forking at t=300 meant both the detected peak/trough
# times AND the state every branch below continues from were contaminated
# by that still-settling transient, not a genuine limit cycle. 500 gives a
# 100yr margin past the observed regularisation point. post_fork_years and
# summary_window keep the same *shape* as before (each branch still runs
# 300 more years past the fork, and the summary window is still the last
# 50 years of that run -- previously t>550 out of a 300-600 run, now
# t>750 out of a 500-800 run) -- only the absolute fork point moved.
t_fork          <- 500
post_fork_years <- 300
summary_window  <- 50

# Runs one resource_decrease value through the requested forks -- an
# unfished reference (used only for peak/trough detection, and skipped
# entirely if neither event schedule is requested), plus whichever of
# peaks-only/troughs-only/constant `schedules` asks for -- all sharing the
# same t=t_fork pre-fork state, so any difference between schedules is
# attributable only to *when* fishing happens, not to a different starting
# point. `schedules` defaults to all three, but a confirmed non-oscillating
# resource_decrease has no real peaks or troughs to key off, only the
# solver's own residual noise floor registering as spurious tiny local
# max/min -- fishing "at" those isn't a meaningful strategy, so the
# non-oscillating call below passes schedules="Constant" to skip building
# them at all rather than reporting a number that isn't measuring anything.
run_focus_yield_comparison <- function(resource_decrease, window, fish_level, regime_label,
                                       t_fork, post_fork_years, summary_window,
                                       schedules = c("Peaks only", "Troughs only", "Constant")) {
  p_focus   <- make_fishing_params(resource_decrease = resource_decrease)
  gp_focus  <- gear_params(p_focus)
  gear_name <- gp_focus$gear[gp_focus$species == "Anchovy"][1]

  # Settle + perturb using Day 18-20's own make_limit_cycle_sim() recipe --
  # the recipe this project already confirmed produces a genuine limit
  # cycle under tr_bdf2, not the numerical ringing Day 6 ruled out, at
  # resource_decrease values where a cycle actually exists. Forked at
  # t=t_fork (past the observed regularisation point): sim_fork is the
  # pre-fork state every branch below continues from; sim_ref continues
  # unfished for a further post_fork_years, used only to detect peak/trough
  # times off an already-regular trajectory, not the settling transient.
  sim_fork <- make_limit_cycle_sim(p_focus, t_total = t_fork, effort = 0)
  sim_ref  <- project(sim_fork, t_max = post_fork_years, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0, method = "tr_bdf2")

  all_times <- as.numeric(dimnames(sim_ref@n)[[1]])
  keep      <- all_times >= t_fork
  t_p       <- all_times[keep]
  time_lbls <- dimnames(sim_ref@n)[[1]][keep]

  needs_events <- any(c("Peaks only", "Troughs only") %in% schedules)
  peak_times   <- numeric(0)
  trough_times <- numeric(0)

  if (needs_events) {
    # Local maxima/minima of the regular (t >= t_fork) unfished biomass
    # trajectory -- same detection method attempt_9-11's original
    # peaks-only branch used, mirrored here to also find troughs. Time
    # labels are reused directly from sim_ref's own dimnames (not
    # re-parsed/re-formatted) to avoid floating-point drift when building
    # the effort schedules below.
    bm <- getBiomass(sim_ref)[, "Anchovy"][time_lbls]
    n  <- length(bm)

    is_peak   <- logical(n)
    is_trough <- logical(n)
    is_peak[2:(n - 1)]   <- bm[2:(n - 1)] > bm[1:(n - 2)] & bm[2:(n - 1)] > bm[3:n]
    is_trough[2:(n - 1)] <- bm[2:(n - 1)] < bm[1:(n - 2)] & bm[2:(n - 1)] < bm[3:n]

    peak_times   <- t_p[is_peak]
    trough_times <- t_p[is_trough]
  }

  # 0/fish_level effort schedule -- "on" only within `window` years either
  # side of each detected event, "off" everywhere else.
  build_schedule <- function(event_times) {
    ev <- rep(0, length(t_p))
    for (et in event_times) {
      ev[t_p >= et - window / 2 & t_p <= et + window / 2] <- fish_level
    }
    matrix(ev, ncol = 1, dimnames = list(time_lbls, gear_name))
  }

  sims <- list()
  if ("Peaks only" %in% schedules) {
    sims[["Peaks only"]] <- project(sim_fork, t_max = post_fork_years, dt = 0.1, t_save = 0.2,
                                    progress_bar = FALSE, effort = build_schedule(peak_times),
                                    method = "tr_bdf2")
  }
  if ("Troughs only" %in% schedules) {
    sims[["Troughs only"]] <- project(sim_fork, t_max = post_fork_years, dt = 0.1, t_save = 0.2,
                                      progress_bar = FALSE, effort = build_schedule(trough_times),
                                      method = "tr_bdf2")
  }
  if ("Constant" %in% schedules) {
    sims[["Constant"]] <- project(sim_fork, t_max = post_fork_years, dt = 0.1, t_save = 0.2,
                                  progress_bar = FALSE, effort = fish_level, method = "tr_bdf2")
  }

  sim_to_df <- function(metric_fn, value_name) {
    bind_rows(lapply(names(sims), function(nm) {
      v <- metric_fn(sims[[nm]])
      setNames(data.frame(as.numeric(dimnames(sims[[nm]]@n)[[1]]), v, nm),
              c("time", value_name, "schedule"))
    }))
  }

  biomass_df <- bind_rows(
    data.frame(time = as.numeric(dimnames(sim_ref@n)[[1]]),
              biomass = getBiomass(sim_ref)[, "Anchovy"], schedule = "Unfished reference"),
    sim_to_df(function(s) getBiomass(s)[, "Anchovy"], "biomass")
  ) %>% filter(time >= t_fork) %>% mutate(regime = regime_label)

  yield_df <- sim_to_df(function(s) getYield(s)[, "Anchovy"], "yield") %>%
    filter(time >= t_fork) %>% mutate(regime = regime_label)

  # The actual number this section was asked for: average yield over the
  # last summary_window years of each run, by schedule.
  summary_t_cut <- t_fork + post_fork_years - summary_window
  summary_df <- yield_df %>%
    filter(time > summary_t_cut) %>%
    group_by(schedule) %>%
    summarise(mean_yield = mean(yield), .groups = "drop") %>%
    mutate(regime = regime_label, n_peaks = length(peak_times), n_troughs = length(trough_times))

  list(biomass_df = biomass_df, yield_df = yield_df, summary_df = summary_df,
      peak_times = peak_times, trough_times = trough_times)
}

# Two regimes: rd_focus sits inside the oscillating range Section 1's
# biomass bifurcation maps out; rd_undepleted=0.1 was confirmed (running
# this section) to settle to a fixed point rather than a limit cycle, so
# it's used as the non-oscillating comparison case Sections 1-2's own
# fixed-point-vs-limit-cycle bifurcation diagrams predict should exist --
# but only its constant-effort schedule is run: with no real cycle,
# "peaks-only"/"troughs-only" would just be fishing wherever the solver's
# own residual noise happens to look like a local max/min, not a real
# strategy.
rd_focus      <- 0.005
rd_undepleted <- 0.1

comparison_oscillating <- run_focus_yield_comparison(
  rd_focus, window, fish_level, regime_label = "Oscillating (resource_decrease=0.005)",
  t_fork = t_fork, post_fork_years = post_fork_years, summary_window = summary_window
)
comparison_nonoscillating <- run_focus_yield_comparison(
  rd_undepleted, window, fish_level, regime_label = "Non-oscillating (resource_decrease=0.1)",
  schedules = "Constant",
  t_fork = t_fork, post_fork_years = post_fork_years, summary_window = summary_window
)

summary_t_cut <- t_fork + post_fork_years - summary_window

cat(sprintf(
  "Section 3: oscillating regime -- %d peaks / %d troughs detected in t=[%g,%g]; non-oscillating regime -- %d peaks / %d troughs detected.\n",
  comparison_oscillating$summary_df$n_peaks[1], comparison_oscillating$summary_df$n_troughs[1],
  t_fork, t_fork + post_fork_years,
  comparison_nonoscillating$summary_df$n_peaks[1], comparison_nonoscillating$summary_df$n_troughs[1]
))

focus_biomass_df <- bind_rows(comparison_oscillating$biomass_df, comparison_nonoscillating$biomass_df)
focus_yield_df    <- bind_rows(comparison_oscillating$yield_df,    comparison_nonoscillating$yield_df)
focus_summary_df  <- bind_rows(comparison_oscillating$summary_df,  comparison_nonoscillating$summary_df)

write.csv(focus_biomass_df, file.path("interesting_plots", "day34_focus_biomass.csv"), row.names = FALSE)
write.csv(focus_yield_df,   file.path("interesting_plots", "day34_focus_yield.csv"),   row.names = FALSE)
write.csv(focus_summary_df, file.path("interesting_plots", "day34_focus_summary.csv"), row.names = FALSE)

# Biomass and yield over time, faceted by regime -- zoomed to the same
# settled t>summary_t_cut window the summary stats themselves use (last
# summary_window years of each run), not the full post-fork range. The
# full range packs ~35 cycles into one plot (post_fork_years=300 /
# ~8.5yr-ish spacing once fishing slows the cycle down), which is dense
# enough to look like solid colour rather than a readable line -- the
# zoomed window shows a handful of cycles clearly instead. CSVs above are
# still written unfiltered; only the plots themselves are restricted.
zoom_biomass_df <- focus_biomass_df %>% filter(time > summary_t_cut)
zoom_yield_df    <- focus_yield_df   %>% filter(time > summary_t_cut)

# Peak/trough markers, restricted to the same zoom window, so it's visible
# *why* the peaks-only/troughs-only schedules turn on where they do,
# rather than just seeing the resulting biomass/yield curves in isolation.
# Only the oscillating regime has these (Section 3's non-oscillating call
# used schedules="Constant", so its peak_times/trough_times are empty).
event_marker_df <- bind_rows(
  data.frame(time = comparison_oscillating$peak_times,   event = "Peak",
            regime = "Oscillating (resource_decrease=0.005)"),
  data.frame(time = comparison_oscillating$trough_times, event = "Trough",
            regime = "Oscillating (resource_decrease=0.005)")
) %>% filter(time > summary_t_cut)

focus_biomass_plot <- ggplot(zoom_biomass_df, aes(x = time, y = biomass, color = schedule)) +
  geom_vline(data = event_marker_df, aes(xintercept = time, linetype = event),
            color = "grey50", linewidth = 0.3, alpha = 0.6) +
  geom_line(linewidth = 1) +
  facet_wrap(~regime, scales = "free_y") +
  scale_linetype_manual(values = c(Peak = "solid", Trough = "dotted")) +
  labs(x = "Time (years)", y = "Biomass", linetype = "Detected event",
       title = "Anchovy biomass: constant vs. peaks-only vs. troughs-only fishing",
       subtitle = sprintf("Zoomed to t=[%g,%g] (the same settled window the summary stats use); effort=%.2g, window=%g years either side of each event; grey lines mark the detected peaks/troughs driving the event-triggered schedules", summary_t_cut, t_fork + post_fork_years, fish_level, window)) +
  theme_minimal()
focus_biomass_plot

save_plot(focus_biomass_plot, "day34_focus_biomass.png", width = 11, height = 6)

focus_yield_plot <- ggplot(zoom_yield_df, aes(x = time, y = yield, color = schedule)) +
  geom_vline(data = event_marker_df, aes(xintercept = time, linetype = event),
            color = "grey50", linewidth = 0.3, alpha = 0.6) +
  geom_line(linewidth = 1) +
  facet_wrap(~regime, scales = "free_y") +
  scale_linetype_manual(values = c(Peak = "solid", Trough = "dotted")) +
  labs(x = "Time (years)", y = "Yield", linetype = "Detected event",
       title = "Anchovy yield: constant vs. peaks-only vs. troughs-only fishing",
       subtitle = sprintf("Zoomed to t=[%g,%g] (the same settled window the summary stats use); effort=%.2g, window=%g years either side of each event", summary_t_cut, t_fork + post_fork_years, fish_level, window)) +
  theme_minimal()
focus_yield_plot

save_plot(focus_yield_plot, "day34_focus_yield.png", width = 11, height = 6)

# The actual comparison this section was built to make: average yield over
# the last summary_window years of each run, by schedule, faceted by
# regime, as a bar chart rather than read off the time-series lines above.
focus_summary_plot <- ggplot(focus_summary_df, aes(x = schedule, y = mean_yield, fill = schedule)) +
  geom_col() +
  facet_wrap(~regime) +
  labs(x = NULL, y = sprintf("Mean yield (t > %g)", summary_t_cut),
       title = sprintf("Average yield after t=%g: constant vs. peaks-only vs. troughs-only fishing", summary_t_cut),
       subtitle = "Compared across an oscillating and a non-oscillating (resource_decrease=0.1) regime") +
  theme_minimal() +
  theme(legend.position = "none")
focus_summary_plot

save_plot(focus_summary_plot, "day34_focus_yield_summary.png", width = 9, height = 6)

print(focus_summary_df)

# Does event-triggered fishing actually EASE the oscillation, or just shift
# its average level? Two complementary answers, both restricted to the
# oscillating regime -- the non-oscillating regime never ran peaks-only/
# troughs-only (schedules="Constant" in that call), so there's nothing to
# compare there.
osc_biomass_wide <- focus_biomass_df %>%
  filter(regime == "Oscillating (resource_decrease=0.005)",
        schedule %in% c("Constant", "Peaks only", "Troughs only")) %>%
  select(time, schedule, biomass) %>%
  pivot_wider(names_from = schedule, values_from = biomass) %>%
  mutate(
    rel_diff_peaks   = (`Peaks only`   - Constant) / Constant,
    rel_diff_troughs = (`Troughs only` - Constant) / Constant
  )

# 1. Relative difference in biomass vs. constant fishing, at every point in
# time -- not the oscillation itself, but how far each event-triggered
# schedule pulls away from the constant-effort baseline as the cycle
# progresses. A relative-difference line that itself keeps oscillating with
# undiminished amplitude means that schedule is shifting the cycle's phase
# or level, not damping it; one that flattens toward a constant offset (or
# toward 0) would mean it's converging onto constant fishing's own
# trajectory.
rel_diff_df <- bind_rows(
  data.frame(time = osc_biomass_wide$time, rel_diff = osc_biomass_wide$rel_diff_peaks,
            comparison = "Peaks only vs. Constant"),
  data.frame(time = osc_biomass_wide$time, rel_diff = osc_biomass_wide$rel_diff_troughs,
            comparison = "Troughs only vs. Constant")
)

# CSV keeps the full t_fork-to-end range (useful for spotting the slow
# beat/envelope pattern between schedules, which only shows up over many
# cycles); the plot itself is zoomed to the same settled window as the
# biomass/yield plots above, for the same readability reason.
write.csv(rel_diff_df, file.path("interesting_plots", "day34_focus_rel_diff.csv"), row.names = FALSE)

zoom_rel_diff_df <- rel_diff_df %>% filter(time > summary_t_cut)

rel_diff_plot <- ggplot(zoom_rel_diff_df, aes(x = time, y = rel_diff, color = comparison)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Time (years)", y = "Relative difference in biomass vs. constant fishing",
       title = "Does event-triggered fishing ease the oscillation, or just shift it?",
       subtitle = sprintf("(schedule - constant) / constant, resource_decrease=%.4g, zoomed to t=[%g,%g] -- the full t=[%g,%g] range is in day34_focus_rel_diff.csv if you want the slower beat pattern across more cycles", rd_focus, summary_t_cut, t_fork + post_fork_years, t_fork, t_fork + post_fork_years)) +
  theme_minimal()
rel_diff_plot

save_plot(rel_diff_plot, "day34_focus_rel_diff.png", width = 10, height = 6)

# 2. The more direct number: relative amplitude of biomass itself,
# (max-min)/((max+min)/2), same convention Sections 1-2's bifurcation
# sweeps use throughout, computed per schedule over the same settled
# t>summary_t_cut window used everywhere else in this section. Lower than
# Constant's own value means that schedule eases the oscillation; higher
# means it makes it worse.
osc_amplitude_df <- focus_biomass_df %>%
  filter(regime == "Oscillating (resource_decrease=0.005)",
        schedule %in% c("Constant", "Peaks only", "Troughs only"),
        time > summary_t_cut) %>%
  group_by(schedule) %>%
  summarise(
    max_biomass   = max(biomass), min_biomass = min(biomass),
    rel_amplitude = (max_biomass - min_biomass) / ((max_biomass + min_biomass) / 2),
    .groups = "drop"
  )

write.csv(osc_amplitude_df, file.path("interesting_plots", "day34_focus_rel_amplitude.csv"), row.names = FALSE)
print(osc_amplitude_df)

osc_amplitude_plot <- ggplot(osc_amplitude_df, aes(x = schedule, y = rel_amplitude, fill = schedule)) +
  geom_col() +
  labs(x = NULL, y = sprintf("Relative amplitude of biomass (t > %g)", summary_t_cut),
       title = "Oscillation amplitude by fishing schedule",
       subtitle = sprintf("(max-min)/((max+min)/2), resource_decrease=%.4g -- lower than Constant's own bar means that schedule eases the oscillation", rd_focus)) +
  theme_minimal() +
  theme(legend.position = "none")
osc_amplitude_plot

save_plot(osc_amplitude_plot, "day34_focus_rel_amplitude.png", width = 8, height = 6)

cat(sprintf("Section 3 (oscillation amplitude, t>%g, oscillating regime only):\n", summary_t_cut))
for (i in seq_len(nrow(osc_amplitude_df))) {
  cat(sprintf("  %s: rel_amplitude=%.4g\n", osc_amplitude_df$schedule[i], osc_amplitude_df$rel_amplitude[i]))
}

################################################################################
# Section 4: scanning window and fish_level for what minimises oscillation
# amplitude while maximising yield
#
# Section 3's peaks-only/troughs-only comparison was one arbitrary
# (window=3, fish_level=0.5) point -- this section scans a grid of both
# instead, at the same rd_focus, to see whether some combination actually
# trades amplitude against yield better than others.
#
# The unfished settle+detection step (make_limit_cycle_sim() out to t_fork,
# then peak/trough detection) doesn't depend on window or fish_level at
# all, so it's done ONCE here rather than once per grid point. Only the
# post-fork project() calls are repeated per combination. Even so, this is
# the most expensive section in the file -- length(window_seq) *
# length(fish_level_seq) * 2 (peaks + troughs) + length(fish_level_seq)
# (constant) separate project() calls, each scan_post_fork_years years
# long. scan_post_fork_years/scan_summary_window are deliberately shorter
# than Section 3's own post_fork_years/summary_window (100/30 vs 300/50):
# a broad scan needs many runs, and each grid point only needs to reach
# ITS OWN settled cycle, not match Section 3's long-run precision.
################################################################################

scan_post_fork_years <- 100
scan_summary_window  <- 30
scan_t_cut           <- t_fork + scan_post_fork_years - scan_summary_window

# window up to 3 (the ~2.9yr period Section 3 detected, so this range spans
# from a genuinely narrow burst up to "nearly the whole cycle"); fish_level
# up to 1 (double Section 3's own 0.5).
window_seq     <- seq(0.5, 3, by = 0.5)
fish_level_seq <- seq(0.2, 1, by = 0.2)

# Generalises the single-resource_decrease scan into a reusable function --
# used first below for the original rd_focus-only full grid, then reused
# twice more further down: once for an extended-effort Constant-only scan
# (window_seq irrelevant there, schedules="Constant" skips detection
# entirely), and once for a coarser grid repeated across several
# resource_decrease values. Settle+detection (the expensive, ~t_fork-year
# part) happens once per call, same efficiency reasoning as before -- it's
# the repeated calls across resource_decrease that make the multi-rd scan
# below expensive, not this function itself.
run_window_effort_scan <- function(resource_decrease, fish_level_seq, t_fork,
                                   scan_post_fork_years, scan_summary_window,
                                   window_seq = NULL,
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

    build_schedule <- function(event_times, window, fish_level) {
      ev <- rep(0, length(t_p))
      for (et in event_times) {
        ev[t_p >= et - window / 2 & t_p <= et + window / 2] <- fish_level
      }
      matrix(ev, ncol = 1, dimnames = list(time_lbls, gear_name))
    }
  }

  # Same two metrics as Sections 1-3: relative amplitude of biomass and
  # mean yield, both over the last scan_summary_window years of the run.
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

  # mean_effort: the time-averaged effort actually applied over the
  # post-fork run -- for Constant this is just fish_level itself (always
  # on); for Peaks/Troughs only, it depends on how much of the run
  # window/fish_level's schedule was actually "on" for. Recorded so
  # amplitude and yield can each be plotted against this one proxy
  # variable directly, rather than only against each other.
  if ("Constant" %in% schedules) {
    out$constant <- bind_rows(lapply(fish_level_seq, function(fl) {
      sim <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                    progress_bar = FALSE, effort = fl, method = "tr_bdf2")
      cbind(scan_metrics(sim), window = NA_real_, fish_level = fl, schedule = "Constant",
           mean_effort = fl)
    }))
  }

  if (needs_events) {
    out$events <- bind_rows(lapply(fish_level_seq, function(fl) {
      bind_rows(lapply(window_seq, function(w) {
        rows <- list()
        if ("Peaks only" %in% schedules) {
          effort_p <- build_schedule(peak_times, w, fl)
          sim_p <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                           progress_bar = FALSE, effort = effort_p, method = "tr_bdf2")
          rows$peaks <- cbind(scan_metrics(sim_p), window = w, fish_level = fl, schedule = "Peaks only",
                              mean_effort = mean(effort_p))
        }
        if ("Troughs only" %in% schedules) {
          effort_t <- build_schedule(trough_times, w, fl)
          sim_t <- project(sim_fork, t_max = scan_post_fork_years, dt = 0.1, t_save = 0.2,
                           progress_bar = FALSE, effort = effort_t, method = "tr_bdf2")
          rows$troughs <- cbind(scan_metrics(sim_t), window = w, fish_level = fl, schedule = "Troughs only",
                                mean_effort = mean(effort_t))
        }
        bind_rows(rows)
      }))
    }))
  }

  bind_rows(out) %>%
    mutate(resource_decrease = resource_decrease, n_peaks = length(peak_times), n_troughs = length(trough_times))
}

scan_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq, window_seq = window_seq,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window
)
event_scan_df    <- scan_df %>% filter(schedule != "Constant")
constant_scan_df <- scan_df %>% filter(schedule == "Constant")

write.csv(scan_df, file.path("interesting_plots", "day34_window_effort_scan.csv"), row.names = FALSE)

# Heatmaps across the (window, fish_level) grid -- Constant is dropped here
# (no window dimension to put on an axis) and shown instead as points in
# the trade-off scatter below.
amplitude_heatmap <- ggplot(event_scan_df, aes(x = window, y = fish_level, fill = rel_amplitude)) +
  geom_tile() +
  facet_wrap(~schedule) +
  scale_fill_viridis_c() +
  labs(x = "Window (years either side of event)", y = "fish_level (effort)",
       fill = "Relative\namplitude",
       title = "Oscillation amplitude across window x fish_level",
       subtitle = sprintf("resource_decrease=%.4g -- darker = lower amplitude = more eased", rd_focus)) +
  theme_minimal()
amplitude_heatmap

save_plot(amplitude_heatmap, "day34_scan_amplitude_heatmap.png", width = 10, height = 5)

yield_heatmap <- ggplot(event_scan_df, aes(x = window, y = fish_level, fill = mean_yield)) +
  geom_tile() +
  facet_wrap(~schedule) +
  scale_fill_viridis_c(option = "magma") +
  labs(x = "Window (years either side of event)", y = "fish_level (effort)",
       fill = "Mean yield",
       title = "Yield across window x fish_level",
       subtitle = sprintf("resource_decrease=%.4g", rd_focus)) +
  theme_minimal()
yield_heatmap

save_plot(yield_heatmap, "day34_scan_yield_heatmap.png", width = 10, height = 5)

# The actual trade-off: amplitude (x, lower=better) vs. yield (y,
# higher=better) -- ideal points sit toward the top-left. Constant's own
# points (one per fish_level, no window) included as the reference to beat.
tradeoff_plot <- ggplot(scan_df, aes(x = rel_amplitude, y = mean_yield, color = schedule, size = fish_level)) +
  geom_point(alpha = 0.7) +
  labs(x = "Relative amplitude of biomass (lower = more eased)",
       y = "Mean yield (higher = better)",
       title = "Amplitude vs. yield trade-off across the window x fish_level scan",
       subtitle = sprintf("resource_decrease=%.4g -- ideal points sit toward the top-left; point size = fish_level", rd_focus)) +
  theme_minimal()
tradeoff_plot

save_plot(tradeoff_plot, "day34_scan_tradeoff.png", width = 9, height = 6)

# Best available combinations: among event-triggered points at or below
# Constant's own median amplitude (i.e. "at least as calm as constant"),
# ranked by yield.
constant_median_amplitude <- median(constant_scan_df$rel_amplitude)
best_points <- scan_df %>%
  filter(schedule != "Constant", rel_amplitude <= constant_median_amplitude) %>%
  arrange(desc(mean_yield)) %>%
  head(5)

write.csv(best_points, file.path("interesting_plots", "day34_scan_best_points.csv"), row.names = FALSE)

cat(sprintf(
  "Section 4: constant fishing's own median rel_amplitude across the fish_level scan = %.4g.\n",
  constant_median_amplitude
))
cat("Top event-triggered (schedule, window, fish_level) combinations at or below that amplitude, ranked by yield:\n")
print(best_points)

################################################################################
# Section 4b: does Constant's amplitude-down/yield-up trend keep going past
# fish_level=1, plateau, or reverse?
#
# The grid above only tested fish_level up to 1, where Constant was still
# improving on both axes with no sign of turning over. schedules="Constant"
# skips peak/trough detection entirely (needs_events short-circuits on it),
# so this is cheap -- one project() call per fish_level, not up to three.
################################################################################

fish_level_seq_constant_ext <- seq(0.2, 4, by = 0.2)

constant_ext_df <- run_window_effort_scan(
  rd_focus, fish_level_seq = fish_level_seq_constant_ext,
  t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window,
  schedules = "Constant"
)

write.csv(constant_ext_df, file.path("interesting_plots", "day34_constant_extended_scan.csv"), row.names = FALSE)

constant_ext_long <- bind_rows(
  data.frame(fish_level = constant_ext_df$fish_level, value = constant_ext_df$rel_amplitude, metric = "Relative amplitude"),
  data.frame(fish_level = constant_ext_df$fish_level, value = constant_ext_df$mean_yield,   metric = "Mean yield")
)

constant_ext_plot <- ggplot(constant_ext_long, aes(x = fish_level, y = value)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 1.5, color = "steelblue") +
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "fish_level (constant effort)", y = NULL,
       title = "Constant fishing: amplitude and yield across a wider effort range",
       subtitle = sprintf("resource_decrease=%.4g -- does the amplitude-down/yield-up trend from fish_level<=1 continue, plateau, or reverse?", rd_focus)) +
  theme_minimal()
constant_ext_plot

save_plot(constant_ext_plot, "day34_constant_extended_scan.png", width = 10, height = 5)

peak_yield_row <- constant_ext_df[which.max(constant_ext_df$mean_yield), ]
min_amp_row    <- constant_ext_df[which.min(constant_ext_df$rel_amplitude), ]
cat(sprintf(
  "Section 4b: across fish_level=[%.2g,%.2g], constant yield peaks at fish_level=%.2g (mean_yield=%.4g); amplitude is lowest at fish_level=%.2g (rel_amplitude=%.4g).\n",
  min(fish_level_seq_constant_ext), max(fish_level_seq_constant_ext),
  peak_yield_row$fish_level, peak_yield_row$mean_yield,
  min_amp_row$fish_level, min_amp_row$rel_amplitude
))

################################################################################
# Section 4c: does the rd_focus=0.005 pattern (Constant dominates,
# Peaks-only can resonate/worsen the oscillation, Troughs-only approaches
# but doesn't beat Constant) hold at other resource_decrease values, or is
# it specific to this one point?
#
# A full 6x5 window/fish_level grid at every resource_decrease would
# multiply Section 4's own already-substantial cost by the number of rd
# values tested -- this uses a coarser 3x3 grid instead (enough to check
# whether the qualitative pattern recurs, not to remeasure it precisely).
# Every resource_decrease value is comfortably inside Section 1's own
# oscillating region (below its ~0.012-0.015 fixed-point boundary).
#
# 12 points, log-spaced from 0.0005 to 0.011 -- densified from the
# original 5-point version once that scan showed a real crossover between
# 0.0025 (event-triggered fishing beats Constant on amplitude) and 0.005
# (Constant wins outright); this is a finer read of where and how sharply
# that crossover actually happens, not just its two endpoints. This is now
# the single most expensive block in the file: 12 x (9 grid points x 2
# schedules + 3 constant points) = 252 project() calls.
################################################################################

rd_scan_seq <- exp(seq(log(0.0005), log(0.011), length.out = 12))

window_seq_coarse     <- c(0.5, 1.5, 3)
fish_level_seq_coarse <- c(0.2, 0.6, 1)

rd_scan_df <- bind_rows(lapply(rd_scan_seq, function(rd) {
  run_window_effort_scan(
    rd, fish_level_seq = fish_level_seq_coarse, window_seq = window_seq_coarse,
    t_fork = t_fork, scan_post_fork_years = scan_post_fork_years, scan_summary_window = scan_summary_window
  )
}))

write.csv(rd_scan_df, file.path("interesting_plots", "day34_rd_window_effort_scan.csv"), row.names = FALSE)

# Same trade-off reading as the single-point scatter above, now faceted by
# resource_decrease -- whether Constant stays top-left in every panel is
# the actual question this section asks.
rd_scan_tradeoff_plot <- ggplot(rd_scan_df, aes(x = rel_amplitude, y = mean_yield, color = schedule, size = fish_level)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~sprintf("resource_decrease=%.4g", resource_decrease), scales = "free", ncol = 4) +
  labs(x = "Relative amplitude of biomass (lower = more eased)",
       y = "Mean yield (higher = better)",
       title = "Amplitude vs. yield trade-off across resource_decrease",
       subtitle = "Coarse 3x3 window x fish_level grid at each resource_decrease -- does Constant stay dominant everywhere?") +
  theme_minimal()
rd_scan_tradeoff_plot

save_plot(rd_scan_tradeoff_plot, "day34_rd_scan_tradeoff.png", width = 16, height = 12)

# Checking specifically whether the peaks-only "resonance hot spot" found
# at rd_focus=0.005 (moderate window, high effort) recurs at other
# resource_decrease values or was specific to that one point.
rd_scan_amplitude_heatmap <- ggplot(rd_scan_df %>% filter(schedule != "Constant"),
                                    aes(x = window, y = fish_level, fill = rel_amplitude)) +
  geom_tile() +
  facet_grid(schedule ~ sprintf("rd=%.4g", resource_decrease)) +
  scale_fill_viridis_c() +
  labs(x = "Window (years either side of event)", y = "fish_level",
       fill = "Relative\namplitude",
       title = "Does the peaks-only resonance hot spot recur across resource_decrease?",
       subtitle = "Darker = lower amplitude = more eased; coarse 3x3 grid per resource_decrease") +
  theme_minimal()
rd_scan_amplitude_heatmap

save_plot(rd_scan_amplitude_heatmap, "day34_rd_scan_amplitude_heatmap.png", width = 24, height = 6)

# The direct generalisation check: for each resource_decrease, does ANY
# event-triggered point beat Constant's own best on amplitude, and
# separately on yield?
rd_scan_verdict <- rd_scan_df %>%
  group_by(resource_decrease) %>%
  summarise(
    constant_best_amplitude = min(rel_amplitude[schedule == "Constant"]),
    constant_best_yield     = max(mean_yield[schedule == "Constant"]),
    any_event_beats_constant_amplitude = any(rel_amplitude[schedule != "Constant"] < constant_best_amplitude),
    any_event_beats_constant_yield     = any(mean_yield[schedule != "Constant"] > constant_best_yield),
    .groups = "drop"
  )

write.csv(rd_scan_verdict, file.path("interesting_plots", "day34_rd_scan_verdict.csv"), row.names = FALSE)
print(rd_scan_verdict)

cat("Section 4c (does Constant dominate at every resource_decrease tested?):\n")
for (i in seq_len(nrow(rd_scan_verdict))) {
  cat(sprintf(
    "  resource_decrease=%.4g: constant's own best amplitude=%.4g, best yield=%.4g -- any event-triggered point beats it on amplitude: %s; on yield: %s\n",
    rd_scan_verdict$resource_decrease[i], rd_scan_verdict$constant_best_amplitude[i], rd_scan_verdict$constant_best_yield[i],
    rd_scan_verdict$any_event_beats_constant_amplitude[i], rd_scan_verdict$any_event_beats_constant_yield[i]
  ))
}

################################################################################
# Section 4d: does the trade-off scatter's own arc/circle-like shape mean
# amplitude and yield are both driven by one underlying "effective fishing
# pressure" variable, rather than window and fish_level acting as two
# independent knobs? mean_effort (the time-averaged effort actually
# applied over the post-fork run -- just fish_level for Constant, but
# depending on how much of the run each window/fish_level combination's
# schedule was actually "on" for Peaks/Troughs only) is that single proxy.
# If both metrics collapse onto one smooth curve against mean_effort, with
# Peaks-only and Troughs-only points landing on the same curve rather than
# separating out, that confirms one hidden variable is doing the work
# rather than window and fish_level mattering independently of it.
################################################################################

mean_effort_amplitude_plot <- ggplot(rd_scan_df, aes(x = mean_effort, y = rel_amplitude, color = schedule)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sprintf("resource_decrease=%.4g", resource_decrease), ncol = 4, scales = "free") +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Relative amplitude of biomass",
       title = "Amplitude vs. mean effort actually applied -- one curve, or separate by schedule?",
       subtitle = "Peaks-only/Troughs-only points overlapping Constant's own curve here would mean window/fish_level only matter through how much they raise mean effort") +
  theme_minimal()
mean_effort_amplitude_plot

save_plot(mean_effort_amplitude_plot, "day34_mean_effort_amplitude.png", width = 16, height = 12)

mean_effort_yield_plot <- ggplot(rd_scan_df, aes(x = mean_effort, y = mean_yield, color = schedule)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~sprintf("resource_decrease=%.4g", resource_decrease), ncol = 4, scales = "free") +
  labs(x = "Mean effort actually applied over the post-fork run",
       y = "Mean yield",
       title = "Yield vs. mean effort actually applied -- one curve, or separate by schedule?",
       subtitle = "Same test as the amplitude version above, applied to yield") +
  theme_minimal()
mean_effort_yield_plot

save_plot(mean_effort_yield_plot, "day34_mean_effort_yield.png", width = 16, height = 12)

################################################################################
# Section 4e: does the "event-triggered fishing beats Constant on amplitude"
# result at low resource_decrease survive Section 3's own longer, more
# careful settled-window methodology, or was it an artefact of Section 4's
# shorter, coarser scan?
#
# Section 4c's scan uses scan_post_fork_years=100 / scan_summary_window=30
# -- deliberately short, for a broad scan. Section 3 itself used
# post_fork_years=300 / summary_window=50 (and the same t_fork=500) for
# exactly this reason: a coarse scan can only ever be a first pass. This
# reruns the two specific (resource_decrease, window, fish_level) points
# that beat Constant on amplitude in the original 5-point version of
# Section 4c, through run_focus_yield_comparison() -- the same function,
# same settle length, Section 3's own headline numbers came from -- rather
# than trusting the coarse scan's numbers at face value.
################################################################################

verify_points <- data.frame(
  resource_decrease = c(0.001, 0.0025),
  window             = c(0.5, 3),
  fish_level         = c(0.6, 0.6)
)

verify_df <- bind_rows(lapply(seq_len(nrow(verify_points)), function(i) {
  rd <- verify_points$resource_decrease[i]
  w  <- verify_points$window[i]
  fl <- verify_points$fish_level[i]

  comparison <- run_focus_yield_comparison(
    rd, window = w, fish_level = fl, regime_label = sprintf("resource_decrease=%.4g", rd),
    t_fork = t_fork, post_fork_years = post_fork_years, summary_window = summary_window
  )

  comparison$biomass_df %>%
    filter(time > summary_t_cut) %>%
    group_by(schedule) %>%
    summarise(
      max_biomass   = max(biomass), min_biomass = min(biomass),
      rel_amplitude = (max_biomass - min_biomass) / ((max_biomass + min_biomass) / 2),
      .groups = "drop"
    ) %>%
    mutate(resource_decrease = rd, window = w, fish_level = fl)
}))

write.csv(verify_df, file.path("interesting_plots", "day34_verify_amplitude.csv"), row.names = FALSE)
print(verify_df)

cat("Section 4e (verifying low-resource_decrease amplitude wins with the careful t_fork=500/post_fork_years=300/summary_window=50 methodology):\n")
for (rd in unique(verify_df$resource_decrease)) {
  sub       <- verify_df[verify_df$resource_decrease == rd, ]
  const_amp <- sub$rel_amplitude[sub$schedule == "Constant"]
  peak_amp  <- sub$rel_amplitude[sub$schedule == "Peaks only"]
  cat(sprintf(
    "  resource_decrease=%.4g: Constant rel_amplitude=%.4g, Peaks only rel_amplitude=%.4g -- still beats Constant: %s\n",
    rd, const_amp, peak_amp, peak_amp < const_amp
  ))
}

################################################################################
# Section 5: summary
################################################################################

check_bifurcation_rd <- function(df) {
  wide <- df %>%
    tidyr::pivot_wider(names_from = c(direction, branch), values_from = metric) %>%
    mutate(
      rel_amplitude_fwd = (Forward_max - Forward_min) / ((Forward_max + Forward_min) / 2),
      rel_amplitude_bwd = (Backward_max - Backward_min) / ((Backward_max + Backward_min) / 2),
      hysteresis_gap    = abs(Forward_max - Backward_max) + abs(Forward_min - Backward_min)
    )
  n_osc     <- sum(wide$rel_amplitude_fwd > 1e-6 | wide$rel_amplitude_bwd > 1e-6, na.rm = TRUE)
  max_gap_v <- wide$value[which.max(wide$hysteresis_gap)]
  list(check = wide, n_oscillating = n_osc, max_hysteresis_value = max_gap_v)
}

biomass_result <- check_bifurcation_rd(bif_biomass_df)
yield_result   <- check_bifurcation_rd(bif_yield_df)

cat("\n===== Day 34 summary =====\n")
cat(sprintf(
  "Section 1 (biomass bifurcation, effort=0, %d-point forward/backward sweep): %d/%d resource_decrease values cross the 1e-6 relative-amplitude threshold (largest hysteresis gap at resource_decrease=%.5g) -- see day34_rd_biomass_bif.png.\n",
  length(rd_seq_bif), biomass_result$n_oscillating, nrow(biomass_result$check), biomass_result$max_hysteresis_value
))
cat(sprintf(
  "Section 2 (yield bifurcation, knife-edge gear at w_mat, effort=0.5, %d-point forward/backward sweep): %d/%d resource_decrease values cross the 1e-6 relative-amplitude threshold (largest hysteresis gap at resource_decrease=%.5g) -- see day34_rd_yield_bif.png.\n",
  length(rd_seq_bif), yield_result$n_oscillating, nrow(yield_result$check), yield_result$max_hysteresis_value
))
cat(sprintf(
  "Section 3 (constant vs. peaks-only vs. troughs-only fishing, effort=%.2g in a %g-year window, oscillating vs. non-oscillating resource_decrease): average yield after t=%g --\n",
  fish_level, window, summary_t_cut
))
for (i in seq_len(nrow(focus_summary_df))) {
  cat(sprintf("  %s / %s: mean_yield=%.4g\n",
             focus_summary_df$regime[i], focus_summary_df$schedule[i], focus_summary_df$mean_yield[i]))
}
cat("  see day34_focus_yield_summary.png / day34_focus_summary.csv.\n")
cat(sprintf(
  "Section 4 (window x fish_level scan, %d windows x %d fish_levels, resource_decrease=%.4g): constant fishing's own median rel_amplitude = %.4g; best event-triggered combinations at or below that --\n",
  length(window_seq), length(fish_level_seq), rd_focus, constant_median_amplitude
))
for (i in seq_len(nrow(best_points))) {
  cat(sprintf("  %s, window=%.2g, fish_level=%.2g: rel_amplitude=%.4g, mean_yield=%.4g\n",
             best_points$schedule[i], best_points$window[i], best_points$fish_level[i],
             best_points$rel_amplitude[i], best_points$mean_yield[i]))
}
cat("  see day34_scan_tradeoff.png / day34_scan_amplitude_heatmap.png / day34_scan_yield_heatmap.png.\n")
cat(sprintf(
  "Section 4b (constant fishing, fish_level=[%.2g,%.2g]): yield peaks at fish_level=%.2g (mean_yield=%.4g); amplitude lowest at fish_level=%.2g (rel_amplitude=%.4g) -- see day34_constant_extended_scan.png.\n",
  min(fish_level_seq_constant_ext), max(fish_level_seq_constant_ext),
  peak_yield_row$fish_level, peak_yield_row$mean_yield, min_amp_row$fish_level, min_amp_row$rel_amplitude
))
cat(sprintf(
  "Section 4c (does Constant dominate across %d resource_decrease values from %.4g to %.4g?): %d/%d had an event-triggered point beat Constant on amplitude; %d/%d beat it on yield -- see day34_rd_scan_tradeoff.png / day34_rd_scan_amplitude_heatmap.png / day34_rd_scan_verdict.csv.\n",
  length(rd_scan_seq), min(rd_scan_seq), max(rd_scan_seq),
  sum(rd_scan_verdict$any_event_beats_constant_amplitude), nrow(rd_scan_verdict),
  sum(rd_scan_verdict$any_event_beats_constant_yield), nrow(rd_scan_verdict)
))
cat("Section 4e (verifying low-resource_decrease amplitude wins with Section 3's own longer settle):\n")
for (rd in unique(verify_df$resource_decrease)) {
  sub       <- verify_df[verify_df$resource_decrease == rd, ]
  const_amp <- sub$rel_amplitude[sub$schedule == "Constant"]
  peak_amp  <- sub$rel_amplitude[sub$schedule == "Peaks only"]
  cat(sprintf(
    "  resource_decrease=%.4g: Constant=%.4g, Peaks only=%.4g -- still beats Constant: %s\n",
    rd, const_amp, peak_amp, peak_amp < const_amp
  ))
}
cat("  see day34_verify_amplitude.csv.\n")

