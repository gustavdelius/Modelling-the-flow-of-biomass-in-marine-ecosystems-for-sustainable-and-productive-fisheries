library(mizer)
library(mizerExperimental)
library(dplyr)
library(ggplot2)
library(future)
library(future.apply)
library(scales)

# Independent grid points/simulations below don't share state -- same
# multisession convention every script since 28_experiments.R.
plan(multisession)

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH (260 chars) truncation guard -- unchanged from
# 29_experiments.R.
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 150) {
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

# NULL-safe scalar extraction -- resource_params()-style lookups have
# silently returned NULL rather than erroring before (24_experiments.R), so
# tryCatch alone isn't enough.
safe_scalar <- function(expr) {
  v <- tryCatch(expr, error = function(e) NULL)
  if (is.null(v) || length(v) == 0) NA_real_ else v[1]
}

################################################################################
# Day 29's "What's Next" left two specific, unfinished threads on cod:
#
#   1. Sweep cod's lambda. Every anchovy-template build in this project
#      fixes lambda=2.05 by convention; cod_params.rds carries whatever
#      value it was actually fitted with, never checked against a swept
#      range. Section 1 reruns the exact perturbed-amplitude heatmap from
#      Day 29 Section 9 (resource_decrease x capacity_mult grid,
#      rel_amplitude as a continuous log-scaled fill) once per lambda in
#      {1.95, 2.00, 2.05, 2.10, 2.15}, to see whether cod's flatness holds
#      across that range or the community slope turns out to matter.
#
#   2. Investigate cod's "bottom left" corner properly. Day 29's own grid
#      found cod settled at 120/120 points but not flat -- a hot spot at
#      low resource_decrease, low capacity_mult (rel_amplitude up to
#      6.3e-7, still 3+ orders of magnitude above the ~1e-11 to 1e-14
#      deep-grid baseline). _day29_cod_corner_check.R spot-checked exactly
#      ONE point there (the grid's own hottest, resource_decrease=
#      0.00187382, capacity_mult=1) at two t_max values and stopped:
#      6.3e-7 at t_max=600, 4.0e-7 at t_max=1200 -- decaying, consistent
#      with a slow transient, but never pushed further, never checked
#      whether nearby points are hotter still, and never saved or plotted.
#      Section 2 widens that properly: a finer scan of the corner region
#      itself, then a longer t_max ladder (up to 19200) run on several of
#      its hottest points, not just one.
################################################################################

################################################################################
# Section 0: rebuild cod + shared helpers (self-contained convention, same
# as every script since Day 20 -- redefined here rather than sourced).
################################################################################

# Read-only: cod_params.rds is never written back to.
cod_params <- readRDS("cod_params.rds")

# Identical to Day 28/29's version. balance=FALSE: setResource() errors
# outright if both resource_rate and resource_capacity are given while
# balance resolves to TRUE, since balance=TRUE always computes the other
# one FOR you to preserve the current steady state -- exactly the
# perturbation these sweeps introduce.
resource_limitation_2d_cod <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}
build_cod <- function(resource_decrease, capacity_mult) {
  resource_limitation_2d_cod(cod_params, resource_decrease, capacity_mult)
}

# Identical to Day 29's make_limit_cycle_sim(): resource dropped to 10% of
# capacity, run 10 time units, then a slice of the consumer spectrum cut by
# 1000x, then run forward to t_total.
make_limit_cycle_sim <- function(params, t_total = 600, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1

  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 1e3

  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "predictor-corrector")
}

# Day 27's own noise floor: ~1e-12 for a settled run, O(0.01)-O(1) for a
# real oscillation -- a 1e-6 cutoff sits well inside that 12-order-of-
# magnitude gap.
classify_perturbed_run <- function(sim, window_frac = 0.4) {
  bm    <- rowSums(getBiomass(sim))
  times <- as.numeric(names(bm))
  late  <- bm[times >= max(times) * (1 - window_frac)]
  rel_amplitude <- (max(late) - min(late)) / mean(late)
  data.frame(rel_amplitude = rel_amplitude,
             perturbed_verdict = if (rel_amplitude > 1e-6) "oscillating" else "settled",
             error = NA_character_)
}

perturbed_stability_point <- function(params, t_total = 600, effort = 0) {
  result <- tryCatch({
    sim <- make_limit_cycle_sim(params, t_total = t_total, effort = effort)
    classify_perturbed_run(sim)
  }, error = function(e) {
    data.frame(rel_amplitude = NA_real_, perturbed_verdict = NA_character_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("perturbed_stability_point(): %s", result$error))
  }
  result
}

perturbed_grid_point <- function(build_fn, resource_decrease, capacity_mult, effort = 0, t_total = 600) {
  result <- tryCatch({
    p <- build_fn(resource_decrease, capacity_mult)
    perturbed_stability_point(p, t_total = t_total, effort = effort)
  }, error = function(e) {
    data.frame(rel_amplitude = NA_real_, perturbed_verdict = NA_character_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("resource_decrease=%.4g, capacity_mult=%.4g: %s",
                    resource_decrease, capacity_mult, result$error))
  }
  data.frame(resource_decrease = resource_decrease, capacity_mult = capacity_mult, result)
}

################################################################################
# Section 1: Sweep cod's lambda through the exact Day 29 Section 9 heatmap
#
# cod_params.rds is opaque -- reconstructing cod at a swept lambda means
# rebuilding via newSingleSpeciesParams() from cod's own pulled traits
# (w_min/w_max/w_mat/gamma/alpha/ks), the same approach Day 29 Section 3's
# "body=cod, rates=cod" hybrid used. Same caveat applies here as there:
# this reproduces those traits at each lambda, not cod_params' full
# calibration (mortality, erepro, R_max, etc. untouched). unname() is
# load-bearing, not defensive styling -- species_params columns carry the
# species name as a names() attribute, and a NAMED no_w breaks gear setup
# downstream with an opaque "number of fishing gears" error (confirmed by
# Day 29 before writing that section).
#
# Cost note: 5 lambdas x the full 12x10 grid = 600 perturbed simulations,
# five times Day 29 Section 9's single-species cost -- split across workers
# via future_lapply, same as every grid in this project.
################################################################################

cod_sp    <- cod_params@species_params
cod_gamma <- unname(cod_sp$gamma[1])
cod_alpha <- unname(cod_sp$alpha[1])
cod_ks    <- unname(cod_sp$ks[1])
cod_w_max <- unname(cod_sp$w_max[1])
cod_w_mat <- unname(cod_sp$w_mat[1])
cod_w_min <- unname(min(cod_params@w))

# Reference context only -- not used to build anything below, so there's no
# risk of the sweep accidentally reproducing whatever this comes back as.
native_cod_lambda <- safe_scalar(resource_params(cod_params)$lambda)

lambda_seq <- seq(1.95, 2.15, by = 0.01)

cat(sprintf(
  "cod_params.rds's own native lambda (resource_params()$lambda): %s -- swept range is %.2f-%.2f (%s).\n",
  ifelse(is.na(native_cod_lambda), "NOT FOUND", sprintf("%.4g", native_cod_lambda)),
  min(lambda_seq), max(lambda_seq),
  if (!is.na(native_cod_lambda) && native_cod_lambda >= min(lambda_seq) && native_cod_lambda <= max(lambda_seq)) {
    "brackets it"
  } else {
    "may not bracket it"
  }
))

make_cod_with_lambda <- function(lambda, resource_decrease, capacity_mult,
                                 second_order = TRUE, ext_diff = 0.00) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))  # same formula every anchovy build uses
  no_w  <- round(log(cod_w_max / cod_w_min) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = sprintf("cod_lambda_%.2f", lambda),
    w_min = cod_w_min, w_max = cod_w_max, w_mat = cod_w_mat,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = cod_alpha, gamma = cod_gamma, ks = cod_ks
  )

  default_capacity <- getResourceCapacity(params)
  r  <- getResourceRate(params) * resource_decrease
  cc <- default_capacity * capacity_mult

  params <- setResource(params, resource_rate = r, resource_capacity = cc,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE)

  params
}

# Factory, not a single closure over a mutable `lambda` -- each
# future_lapply() worker below needs a build_fn bound to a FIXED lambda
# value at creation time, not a shared variable that could change between
# the closure's creation and its eventual call.
make_build_cod_lambda <- function(lambda) {
  function(resource_decrease, capacity_mult) {
    make_cod_with_lambda(lambda, resource_decrease, capacity_mult)
  }
}

# Same grid as Day 29 Section 9 (resource_decrease_seq x capacity_mult_seq,
# 12x10=120 points), so this heatmap is directly comparable point-for-point
# against Day 29's cod panel (which used lambda=cod_params' own native
# value throughout).
resource_decrease_seq <- exp(seq(log(0.001), log(1), length.out = 12))
capacity_mult_seq     <- exp(seq(log(1), log(20), length.out = 10))
lambda_grid_params    <- expand.grid(resource_decrease = resource_decrease_seq,
                                     capacity_mult = capacity_mult_seq)

perturbed_grid_generic <- function(build_fn, label) {
  df <- bind_rows(future_lapply(seq_len(nrow(lambda_grid_params)), function(i) {
    perturbed_grid_point(build_fn, lambda_grid_params$resource_decrease[i],
                         lambda_grid_params$capacity_mult[i])
  }, future.seed = TRUE))
  df$lambda <- label
  df
}

cod_lambda_grid_df <- bind_rows(lapply(lambda_seq, function(lam) {
  perturbed_grid_generic(make_build_cod_lambda(lam), lam)
}))

print(cod_lambda_grid_df %>% filter(is.na(error)) %>% group_by(lambda) %>%
        summarise(min_amplitude = min(rel_amplitude), max_amplitude = max(rel_amplitude),
                  n_oscillating = sum(perturbed_verdict == "oscillating"), .groups = "drop"))
write.csv(cod_lambda_grid_df, file.path("interesting_plots", "day30_cod_lambda_grid.csv"),
          row.names = FALSE)

# Floored, not the raw column, so a perfectly-settled point doesn't break
# the log10 fill scale -- the CSV above keeps the unfloored value. Same
# convention as Day 29 Section 9.
cod_lambda_plot_df <- cod_lambda_grid_df %>%
  filter(is.na(error)) %>%
  mutate(rel_amplitude_floored = pmax(rel_amplitude, 1e-15))

cod_lambda_plot <- ggplot(cod_lambda_plot_df,
                          aes(x = resource_decrease, y = capacity_mult, fill = rel_amplitude_floored)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(trans = "log10", labels = scales::label_scientific(),
                       name = "rel_amplitude\n(perturbed sim)") +
  facet_wrap(~lambda, labeller = labeller(lambda = function(x) paste0("lambda = ", x))) +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       title = "Cod's perturbed-amplitude heatmap across a swept lambda",
       subtitle = "Same grid/method as Day 29 Section 9 -- cod's own gamma/alpha/ks/w_max/w_mat held fixed, lambda only varied") +
  theme_minimal()
cod_lambda_plot

save_plot(cod_lambda_plot, "day30_cod_lambda_heatmap.png", width = 13, height = 9)

cat(sprintf(
  "\nSection 1 summary: %d/%d total points oscillating across all %d lambda values tested (%s).\n",
  sum(cod_lambda_grid_df$perturbed_verdict == "oscillating", na.rm = TRUE),
  nrow(cod_lambda_grid_df), length(lambda_seq), paste(lambda_seq, collapse = ", ")
))

################################################################################
# Section 2: Investigate cod's "bottom left" corner properly
#
# Three things _day29_cod_corner_check.R didn't do:
#   (a) A finer scan of the corner region itself, at Day 29's own t_max=600,
#       to confirm the coarse grid's "hottest point" actually is the local
#       peak rather than just the nearest grid line to one.
#   (b) The same t_max escalation as the spot-check, but for several corner
#       points (not just one) and out to a longer ladder (up to t_max=
#       19200, one doubling past the spot-check's furthest point).
#   (c) One long simulation per point (to the ladder's longest t_max),
#       sliced at each shorter cutoff, rather than independently
#       re-simulating from t=0 at every t_max -- the same trick Day 26
#       Section 3b used for the delay-DE tau sweep: a deterministic,
#       fixed-step integration's trajectory up to any given time doesn't
#       depend on how much further it's eventually run, so slicing one
#       t_max=19200 run reproduces what independent runs at 600, 1200,
#       2400, 4800, and 9600 would each give, at roughly a fifth of the
#       simulation cost. (Day 26 established this reasoning for a
#       fixed-step DDE solver; extending it to mizer's fixed-dt
#       predictor-corrector project() is analogous but not independently
#       reverified here.)
################################################################################

# Zoomed into the region Day 29's coarse 12x10 grid flagged as hot --
# resource_decrease from the grid's own minimum up to roughly its 5th
# point, capacity_mult from its minimum up to roughly its 3rd -- at finer
# resolution than the original grid had there (8x5=40 points vs. the
# original's 5x3=15 covering the same span).
corner_resource_decrease_seq <- exp(seq(log(0.001), log(0.02), length.out = 8))
corner_capacity_mult_seq     <- exp(seq(log(1), log(3), length.out = 5))
corner_grid_params <- expand.grid(resource_decrease = corner_resource_decrease_seq,
                                  capacity_mult = corner_capacity_mult_seq)

corner_scan_df <- bind_rows(future_lapply(seq_len(nrow(corner_grid_params)), function(i) {
  perturbed_grid_point(build_cod, corner_grid_params$resource_decrease[i],
                       corner_grid_params$capacity_mult[i])
}, future.seed = TRUE))

cat("\nTop 10 points from the finer corner scan (t_max=600):\n")
print(corner_scan_df %>% filter(is.na(error)) %>% arrange(desc(rel_amplitude)) %>% head(10))
write.csv(corner_scan_df, file.path("interesting_plots", "day30_corner_scan.csv"),
          row.names = FALSE)

corner_scan_plot_df <- corner_scan_df %>%
  filter(is.na(error)) %>%
  mutate(rel_amplitude_floored = pmax(rel_amplitude, 1e-15))

corner_scan_plot <- ggplot(corner_scan_plot_df,
                           aes(x = resource_decrease, y = capacity_mult, fill = rel_amplitude_floored)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(trans = "log10", labels = scales::label_scientific(),
                       name = "rel_amplitude\n(t_max=600)") +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       title = "Cod's bottom-left corner, at finer resolution than Day 29's grid",
       subtitle = "t_max=600 quick scan -- locates the corner's actual hot region before the t_max escalation below") +
  theme_minimal()
corner_scan_plot
save_plot(corner_scan_plot, "day30_corner_scan_heatmap.png")

# Top 5 hottest points from the finer scan -- these, not just the coarse
# grid's single hottest point, get the full t_max escalation below.
top_corner_points <- corner_scan_df %>%
  filter(is.na(error)) %>%
  arrange(desc(rel_amplitude)) %>%
  slice_head(n = 5) %>%
  mutate(point_id = sprintf("rd=%.4g, cm=%.3g", resource_decrease, capacity_mult))

cat("\nTop 5 corner points carried into the t_max escalation:\n")
print(top_corner_points %>% select(point_id, resource_decrease, capacity_mult, rel_amplitude))

# Late-window amplitude read off a slice of an already-computed long
# trajectory, rather than a fresh sim per t_max -- see (c) above. Same
# "late window = last 40% of the range up to this cutoff" convention as
# Day 26's amplitude_at_cutoff().
amplitude_at_cutoff_mizer <- function(t_max_effective, sim, window_frac = 0.4) {
  bm    <- rowSums(getBiomass(sim))
  times <- as.numeric(names(bm))
  late  <- bm[times > t_max_effective * (1 - window_frac) & times <= t_max_effective]
  data.frame(t_max = t_max_effective,
             rel_amplitude = (max(late) - min(late)) / mean(late))
}

t_max_ladder <- c(300,600,900,1200)

run_corner_point_ladder <- function(resource_decrease, capacity_mult, point_id) {
  p   <- build_cod(resource_decrease, capacity_mult)
  sim <- make_limit_cycle_sim(p, t_total = max(t_max_ladder))
  ladder_df <- bind_rows(lapply(t_max_ladder, amplitude_at_cutoff_mizer, sim = sim))
  ladder_df$point_id           <- point_id
  ladder_df$resource_decrease  <- resource_decrease
  ladder_df$capacity_mult      <- capacity_mult

  # Cross-check against the linear verdict at the same point, same
  # convention as Day 29 Section 8's dual_stability_check().
  stab <- tryCatch({
    p_steady <- steadyNewton(p, stability = TRUE)
    attr(p_steady, "stability")
  }, error = function(e) NULL)
  ladder_df$spectral_radius <- if (is.null(stab)) NA_real_ else stab$spectral_radius
  ladder_df
}

corner_ladder_df <- bind_rows(future_lapply(seq_len(nrow(top_corner_points)), function(i) {
  run_corner_point_ladder(top_corner_points$resource_decrease[i],
                          top_corner_points$capacity_mult[i],
                          top_corner_points$point_id[i])
}, future.seed = TRUE))

print(corner_ladder_df %>% select(point_id, t_max, rel_amplitude, spectral_radius))
write.csv(corner_ladder_df, file.path("interesting_plots", "day30_corner_ladder.csv"),
          row.names = FALSE)

corner_ladder_plot <- ggplot(corner_ladder_df, aes(x = t_max, y = rel_amplitude, color = point_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 1e-6, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_scientific()) +
  labs(x = "t_max (how long the run is allowed to settle)",
       y = "Relative amplitude, late window",
       color = "corner point",
       title = "Does cod's bottom-left corner decay to the settled floor, or plateau?",
       subtitle = "One long simulation per point, sliced at each t_max -- dashed line = the oscillating/settled threshold used throughout this project") +
  theme_minimal()
corner_ladder_plot
save_plot(corner_ladder_plot, "day30_corner_decay_ladder.png")

# Verdict per point: compare the two longest t_max steps. Still shrinking
# by more than 15% at the ladder's far end -> genuinely still decaying,
# consistent with a transient. Flat within that, but still under the 1e-6
# oscillating/settled threshold -> plateaued, worth flagging rather than
# assuming away. Dropped near the deep-grid noise floor -> confirmed
# transient.
corner_verdict_df <- corner_ladder_df %>%
  group_by(point_id, resource_decrease, capacity_mult) %>%
  arrange(t_max, .by_group = TRUE) %>%
  summarise(
    rel_amplitude_first = first(rel_amplitude),
    rel_amplitude_last  = last(rel_amplitude),
    last_step_ratio      = last(rel_amplitude) / nth(rel_amplitude, -2),
    verdict = case_when(
      rel_amplitude_last < 1e-9 ~ "decayed to settled floor -- transient confirmed",
      last_step_ratio    < 0.85 ~ "still shrinking at t_max=19200 -- transient, not yet fully settled",
      TRUE                       ~ "plateaued -- below the oscillating threshold, but not decaying further; genuine small-amplitude structure can't be ruled out"
    ),
    .groups = "drop"
  )

print(corner_verdict_df %>% select(point_id, rel_amplitude_first, rel_amplitude_last, last_step_ratio, verdict))
write.csv(corner_verdict_df, file.path("interesting_plots", "day30_corner_verdict.csv"),
          row.names = FALSE)

cat("\nSection 2 verdict, per corner point:\n")
for (i in seq_len(nrow(corner_verdict_df))) {
  cat(sprintf("  %s: %.4g -> %.4g (t_max 600 -> 19200). %s\n",
             corner_verdict_df$point_id[i], corner_verdict_df$rel_amplitude_first[i],
             corner_verdict_df$rel_amplitude_last[i], corner_verdict_df$verdict[i]))
}

################################################################################
# Section 3: Broad coarse scan -- much wider resource_decrease/capacity_mult
# range than any grid so far in this project
#
# Section 1 found that cod's own traits (gamma/alpha/ks/w_max/w_mat), once
# rebuilt at a low enough lambda, DO cross into genuine oscillation
# (rel_amplitude > 1e-6) -- but only within the same resource_decrease
# (1e-3 to 1) x capacity_mult (1 to 20) window every grid in this project
# has used since Day 24. Before assuming that window is the right one to
# keep refining, this runs one much wider, much coarser reconnaissance pass
# on NATIVE cod_params (build_cod, no lambda substitution) -- resource_
# decrease down to 1e-7 (four orders of magnitude further into resource
# starvation than anything tested before) and capacity_mult from 1e-3 (an
# actively SHRUNK carrying capacity, never tested -- every grid so far only
# ever multiplied capacity up) to 1e2 (5x past the previous 20x ceiling).
#
# Coarse on purpose -- 8x6=48 points, a reconnaissance pass meant to flag
# which corner of this much bigger space deserves a properly resolved
# follow-up grid, not a finished answer on its own.
#
# At resource_decrease=1e-7 the resource renews at ~0 -- a real risk of
# population collapse (biomass -> 0) rather than a sustained oscillation.
# rel_amplitude = (max-min)/mean blows up to NaN/Inf under that condition
# (0/0), which would silently break the log10 fill scale below, so each
# point also checks final biomass against a collapse floor and reports it
# as its own category rather than letting it masquerade as "oscillating".
################################################################################

broad_resource_decrease_seq <- exp(seq(log(1e-7), log(1), length.out = 8))
broad_capacity_mult_seq     <- exp(seq(log(1e-3), log(1e2), length.out = 6))
broad_grid_params <- expand.grid(resource_decrease = broad_resource_decrease_seq,
                                 capacity_mult = broad_capacity_mult_seq)

# Same shape as perturbed_grid_point(), but also reads the settled-window's
# mean biomass so a collapsed run doesn't get classified as "oscillating"
# off a NaN/Inf rel_amplitude.
broad_scan_point <- function(build_fn, resource_decrease, capacity_mult, t_total = 600) {
  result <- tryCatch({
    p     <- build_fn(resource_decrease, capacity_mult)
    sim   <- make_limit_cycle_sim(p, t_total = t_total)
    bm    <- rowSums(getBiomass(sim))
    times <- as.numeric(names(bm))
    late  <- bm[times >= max(times) * 0.6]
    mean_late_biomass <- mean(late)
    # Collapse floor: four orders of magnitude below the run's own starting
    # biomass is well past "still oscillating" and into "effectively extinct".
    initial_biomass <- bm[1]
    collapsed <- mean_late_biomass < initial_biomass * 1e-6
    rel_amplitude <- if (collapsed) NA_real_ else (max(late) - min(late)) / mean_late_biomass
    data.frame(rel_amplitude = rel_amplitude,
               perturbed_verdict = if (collapsed) {
                 "collapsed"
               } else if (rel_amplitude > 1e-6) {
                 "oscillating"
               } else {
                 "settled"
               },
               mean_late_biomass = mean_late_biomass, error = NA_character_)
  }, error = function(e) {
    data.frame(rel_amplitude = NA_real_, perturbed_verdict = NA_character_,
               mean_late_biomass = NA_real_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("resource_decrease=%.4g, capacity_mult=%.4g: %s",
                    resource_decrease, capacity_mult, result$error))
  }
  data.frame(resource_decrease = resource_decrease, capacity_mult = capacity_mult, result)
}

broad_scan_df <- bind_rows(future_lapply(seq_len(nrow(broad_grid_params)), function(i) {
  broad_scan_point(build_cod, broad_grid_params$resource_decrease[i],
                   broad_grid_params$capacity_mult[i])
}, future.seed = TRUE))

print(broad_scan_df %>% filter(is.na(error)) %>% count(perturbed_verdict))
write.csv(broad_scan_df, file.path("interesting_plots", "day30_broad_scan.csv"),
          row.names = FALSE)

# Continuous amplitude fill, collapsed points excluded (their rel_amplitude
# is NA by construction) -- see the categorical plot below for where they
# actually sit.
broad_scan_plot_df <- broad_scan_df %>%
  filter(is.na(error), perturbed_verdict != "collapsed") %>%
  mutate(rel_amplitude_floored = pmax(rel_amplitude, 1e-15))

broad_scan_plot <- ggplot(broad_scan_plot_df,
                          aes(x = resource_decrease, y = capacity_mult, fill = rel_amplitude_floored)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(trans = "log10", labels = scales::label_scientific(),
                       name = "rel_amplitude\n(perturbed sim)") +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       title = "Native cod_params, coarse reconnaissance across a much wider range",
       subtitle = sprintf("resource_decrease %.0e-%.0e, capacity_mult %.0e-%.0e -- %d/%d points collapsed (excluded from this fill)",
                          min(broad_resource_decrease_seq), max(broad_resource_decrease_seq),
                          min(broad_capacity_mult_seq), max(broad_capacity_mult_seq),
                          sum(broad_scan_df$perturbed_verdict == "collapsed", na.rm = TRUE),
                          nrow(broad_scan_df))) +
  theme_minimal()
broad_scan_plot
save_plot(broad_scan_plot, "day30_broad_scan_heatmap.png", width = 10)

# Collapsed points shown as their own category here -- extinction is a
# different phenomenon from oscillation, worth seeing where the boundary
# between the two falls rather than just excluding it silently.
broad_scan_verdict_plot <- ggplot(broad_scan_df %>% filter(is.na(error)),
                                  aes(x = resource_decrease, y = capacity_mult, fill = perturbed_verdict)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c(settled = "#4C72B0", oscillating = "#C44E52", collapsed = "grey30"),
                    na.value = "grey80") +
  labs(x = "resource_decrease", y = "capacity_mult", fill = NULL,
       title = "Same coarse scan, categorical: settled / oscillating / collapsed") +
  theme_minimal()
broad_scan_verdict_plot
save_plot(broad_scan_verdict_plot, "day30_broad_scan_verdict.png", width = 10)

cat(sprintf(
  "\nSection 3 summary: %d/%d settled, %d/%d oscillating, %d/%d collapsed (native cod_params, coarse %dx%d scan).\n",
  sum(broad_scan_df$perturbed_verdict == "settled", na.rm = TRUE), nrow(broad_scan_df),
  sum(broad_scan_df$perturbed_verdict == "oscillating", na.rm = TRUE), nrow(broad_scan_df),
  sum(broad_scan_df$perturbed_verdict == "collapsed", na.rm = TRUE), nrow(broad_scan_df),
  length(broad_resource_decrease_seq), length(broad_capacity_mult_seq)
))

################################################################################
# Section 4: Sweep cod's native q (mizer's search-volume exponent) through
# the same perturbed-amplitude heatmap
#
# NOT de Roos & Persson's q (the juvenile/adult competitiveness ratio Day
# 29's q_equivalent computed as a diagnostic) -- that's an unrelated
# concept that happens to share a letter. This is mizer's own species_params
# column: SearchVolume = gamma * w^q, consumed directly by getEncounter().
#
# Unlike Section 1's lambda sweep, q is a genuine column on the REAL,
# fully-calibrated cod_params object, so this is applied straight to it via
# given_species_params() -- no newSingleSpeciesParams() reconstruction, no
# Section 1-style caveat about only reproducing some traits. Mortality,
# erepro, R_max, gear, and cod's own native lambda all stay exactly as
# fitted; only q changes.
################################################################################

native_cod_q <- unname(cod_params@species_params$q[1])
cat(sprintf("cod_params.rds's own native q (search-volume exponent): %.4g\n", native_cod_q))

# given_species_params() has been a silent no-op before in this exact
# project: D_ext gets attached but is never actually consumed without a
# follow-up setExtDiffusion() call (Day 20's own notes flag this as
# unresolved), and Day 29 separately found capacity_mult sitting unused in
# an earlier make_anchovy_params() draft. q feeds into SearchVolume, which
# mizer may cache as params@search_vol at construction time rather than
# recomputing it live from species_params -- so before trusting an entire
# grid built on this, check directly whether search_vol/getEncounter()
# actually move when q changes.
q_probe_base <- cod_params
q_probe_low  <- cod_params
given_species_params(q_probe_low)$q <- native_cod_q - 0.3

search_vol_changed <- !isTRUE(all.equal(q_probe_base@search_vol, q_probe_low@search_vol))
encounter_changed  <- !isTRUE(all.equal(getEncounter(q_probe_base), getEncounter(q_probe_low)))

cat(sprintf("q propagation check: search_vol changed = %s, getEncounter() changed = %s.\n",
           search_vol_changed, encounter_changed))
if (!search_vol_changed && !encounter_changed) {
  warning(paste(
    "given_species_params(params)$q <- value does NOT appear to change",
    "search_vol or getEncounter() -- q may be cached at construction time,",
    "and the sweep below could be silently doing nothing. Consider",
    "species_params(params)$q <- value (full replacement, triggers",
    "recalculation) instead before trusting the grid below."
  ))
}

make_build_cod_q <- function(q_value) {
  function(resource_decrease, capacity_mult) {
    p <- cod_params
    given_species_params(p)$q <- q_value
    resource_limitation_2d_cod(p, resource_decrease, capacity_mult)
  }
}

# Centred on cod's own fitted value, +/- 0.3 in steps of 0.1 -- same
# "bracket the native value" convention as Section 1's lambda_seq. Reuses
# the exact resource_decrease x capacity_mult grid from Sections 1/2
# (lambda_grid_params) so this is directly comparable point-for-point.
q_seq <- round(native_cod_q + seq(-0.3, 0.3, by = 0.1), 4)

perturbed_grid_generic_q <- function(build_fn, label) {
  df <- bind_rows(future_lapply(seq_len(nrow(lambda_grid_params)), function(i) {
    perturbed_grid_point(build_fn, lambda_grid_params$resource_decrease[i],
                         lambda_grid_params$capacity_mult[i])
  }, future.seed = TRUE))
  df$q <- label
  df
}

cod_q_grid_df <- bind_rows(lapply(q_seq, function(qv) {
  perturbed_grid_generic_q(make_build_cod_q(qv), qv)
}))

print(cod_q_grid_df %>% filter(is.na(error)) %>% group_by(q) %>%
        summarise(min_amplitude = min(rel_amplitude), max_amplitude = max(rel_amplitude),
                  n_oscillating = sum(perturbed_verdict == "oscillating"), .groups = "drop"))
write.csv(cod_q_grid_df, file.path("interesting_plots", "day30_cod_q_grid.csv"),
          row.names = FALSE)

cod_q_plot_df <- cod_q_grid_df %>%
  filter(is.na(error)) %>%
  mutate(rel_amplitude_floored = pmax(rel_amplitude, 1e-15))

cod_q_plot <- ggplot(cod_q_plot_df,
                     aes(x = resource_decrease, y = capacity_mult, fill = rel_amplitude_floored)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(trans = "log10", labels = scales::label_scientific(),
                       name = "rel_amplitude\n(perturbed sim)") +
  facet_wrap(~q, labeller = labeller(q = function(x) paste0("q = ", x))) +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       title = "Cod's perturbed-amplitude heatmap across a swept search-volume exponent q",
       subtitle = sprintf(
         "Native cod_params.rds throughout (mortality/erepro/R_max/lambda untouched) -- only q overridden via given_species_params(); native q=%.3g",
         native_cod_q
       )) +
  theme_minimal()
cod_q_plot

save_plot(cod_q_plot, "day30_cod_q_heatmap.png", width = 13, height = 9)

cat(sprintf(
  "\nSection 4 summary: %d/%d total points oscillating across all %d q values tested (%s).\n",
  sum(cod_q_grid_df$perturbed_verdict == "oscillating", na.rm = TRUE),
  nrow(cod_q_grid_df), length(q_seq), paste(q_seq, collapse = ", ")
))

################################################################################
# Section 5: Day 29's juvenile/adult metric was ITSELF unweighted across
# unequal-width bins -- a second bug stacked on top of the difference-vs-
# ratio fix Day 29 already made
#
# mean_mass_specific_rate() (Day 29) computed a plain mean() over
# mass_specific_rate[idx], where idx selects the juvenile or adult weight
# bins. mizer's weight grid is logarithmically spaced -- dw grows with w --
# so a plain mean() silently treats every BIN as equally important
# regardless of how much weight range it actually represents. Small-w bins
# are far more numerous than large-w bins on a log grid, so an unweighted
# mean over either the juvenile or the adult subset systematically
# overweights the smallest individuals in that stage (the newest
# juveniles, the just-matured adults) relative to the largest -- the
# opposite of a population-representative average. Exactly the kind of bug
# this project has hit before with a plain mean()/sum() over a
# non-uniform mizer grid.
#
# Fix: weight each bin by dw (the actual width of weight it represents)
# instead of counting every bin equally -- a Riemann-sum-consistent
# average rather than an average-over-indices. Self-contained: redefines
# make_anchovy_params()/anchovy_params and Day 29's ORIGINAL (unweighted)
# mean_mass_specific_rate() here too, so the before/after comparison is
# generated fresh rather than assumed.
################################################################################

make_anchovy_params <- function(lambda = 2.05, resource_decrease = 0.01,
                                capacity_mult = 10, second_order = TRUE,
                                ext_diff = 0.00, alpha = 0.1, gamma = 750,
                                ks = 0, kappa_override = NULL) {
  a0    <- 100
  kappa <- if (is.null(kappa_override)) a0 * exp(-6.9 * (lambda - 1)) else kappa_override
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = alpha, gamma = gamma, ks = ks
  )

  given_species_params(params)$D_ext <- ext_diff
  params <- setBevertonHolt(params)

  default_capacity <- getResourceCapacity(params)
  r  <- getResourceRate(params) * resource_decrease
  cc <- default_capacity * capacity_mult

  params <- setResource(params, resource_rate = r, resource_capacity = cc,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE)

  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }

  params
}
anchovy_params <- make_anchovy_params()

# Day 29's ORIGINAL metric, unchanged -- kept here so the corrected version
# below can be compared against it directly rather than against a
# remembered number.
mean_mass_specific_rate <- function(params, idx = TRUE) {
  E <- getEncounter(params)
  f <- getFeedingLevel(params)
  per_capita_rate    <- (E * (1 - f))[1, ]
  mass_specific_rate <- per_capita_rate / params@w
  mean(mass_specific_rate[idx])
}

# The fix: weighted.mean(..., w = dw) instead of mean(). A biomass-weighted
# version (weighting by initialN()*w*dw, i.e. by how much of the
# population's actual mass sits in each bin) would be an even more
# defensible choice -- dw-weighting alone corrects the grid-density bias
# but still treats an empty bin the same as a crowded one -- but dw is the
# direct fix for the specific bug described (bins of different sizes going
# into an unweighted mean), so that's what's implemented here.
mean_mass_specific_rate_weighted <- function(params, idx = TRUE) {
  E <- getEncounter(params)
  f <- getFeedingLevel(params)
  per_capita_rate    <- (E * (1 - f))[1, ]
  mass_specific_rate <- per_capita_rate / params@w
  weighted.mean(mass_specific_rate[idx], w = params@dw[idx])
}

corrected_timescale_summary <- function(params, label) {
  w_mat <- params@species_params$w_mat[1]
  w     <- params@w

  msr_juv_unweighted <- mean_mass_specific_rate(params, w < w_mat)
  msr_ad_unweighted  <- mean_mass_specific_rate(params, w >= w_mat)
  msr_juv_weighted   <- mean_mass_specific_rate_weighted(params, w < w_mat)
  msr_ad_weighted    <- mean_mass_specific_rate_weighted(params, w >= w_mat)

  ratio_unweighted <- msr_juv_unweighted / msr_ad_unweighted
  ratio_weighted   <- msr_juv_weighted / msr_ad_weighted

  data.frame(
    species                    = label,
    juvenile_adult_ratio_unweighted = ratio_unweighted,
    q_equivalent_unweighted          = 2 / (ratio_unweighted + 1),
    juvenile_adult_ratio_dw_weighted = ratio_weighted,
    q_equivalent_dw_weighted          = 2 / (ratio_weighted + 1)
  )
}

corrected_timescale_df <- bind_rows(
  corrected_timescale_summary(cod_params, "cod"),
  corrected_timescale_summary(anchovy_params, "anchovy")
)
print(corrected_timescale_df)
write.csv(corrected_timescale_df, file.path("interesting_plots", "day30_corrected_timescale.csv"),
          row.names = FALSE)

cat(sprintf(
  "\nSection 5 summary: dw-weighting flips juvenile_adult_ratio for cod from %.4g to %.4g, and for anchovy from %.4g to %.4g (>1 = adult-driven, <1 = juvenile-driven).\n",
  corrected_timescale_df$juvenile_adult_ratio_unweighted[corrected_timescale_df$species == "cod"],
  corrected_timescale_df$juvenile_adult_ratio_dw_weighted[corrected_timescale_df$species == "cod"],
  corrected_timescale_df$juvenile_adult_ratio_unweighted[corrected_timescale_df$species == "anchovy"],
  corrected_timescale_df$juvenile_adult_ratio_dw_weighted[corrected_timescale_df$species == "anchovy"]
))

################################################################################
# Section 6: Summary -- programmatic readout, not asserted conclusions.
################################################################################

cat("\n===== Day 30 summary =====\n")
cat(sprintf(
  "Section 1 (cod lambda sweep, %d values x %d-point grid): %d/%d points oscillating overall.\n",
  length(lambda_seq), nrow(lambda_grid_params),
  sum(cod_lambda_grid_df$perturbed_verdict == "oscillating", na.rm = TRUE), nrow(cod_lambda_grid_df)
))
cat(sprintf(
  "Section 2 (bottom-left corner, %d points x t_max up to %d): %d transient/still-decaying, %d plateaued.\n",
  nrow(top_corner_points), max(t_max_ladder),
  sum(grepl("transient|decayed", corner_verdict_df$verdict)),
  sum(grepl("plateaued", corner_verdict_df$verdict))
))
cat(sprintf(
  "Section 3 (broad coarse scan, native cod_params, %d points): %d oscillating, %d collapsed.\n",
  nrow(broad_scan_df),
  sum(broad_scan_df$perturbed_verdict == "oscillating", na.rm = TRUE),
  sum(broad_scan_df$perturbed_verdict == "collapsed", na.rm = TRUE)
))
cat(sprintf(
  "Section 4 (cod q sweep, %d values x %d-point grid): %d/%d points oscillating overall.\n",
  length(q_seq), nrow(lambda_grid_params),
  sum(cod_q_grid_df$perturbed_verdict == "oscillating", na.rm = TRUE), nrow(cod_q_grid_df)
))
cat(sprintf(
  "Section 5 (corrected competitiveness metric): cod ratio %.4g -> %.4g, anchovy ratio %.4g -> %.4g (unweighted -> dw-weighted).\n",
  corrected_timescale_df$juvenile_adult_ratio_unweighted[corrected_timescale_df$species == "cod"],
  corrected_timescale_df$juvenile_adult_ratio_dw_weighted[corrected_timescale_df$species == "cod"],
  corrected_timescale_df$juvenile_adult_ratio_unweighted[corrected_timescale_df$species == "anchovy"],
  corrected_timescale_df$juvenile_adult_ratio_dw_weighted[corrected_timescale_df$species == "anchovy"]
))
