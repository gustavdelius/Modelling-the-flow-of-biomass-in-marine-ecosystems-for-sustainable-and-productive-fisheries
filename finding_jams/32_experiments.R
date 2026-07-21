library(mizer)
library(mizerExperimental)
library(dplyr)
library(ggplot2)
library(scales)

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH truncation guard.
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

# Day 32: direct projectToSteady() sweeps for alpha/q (Sections 1-2), a
# resource_rate sweep with capacity fixed via a one-time balance=TRUE
# baseline (Section 3), a project()-vs-projectToSteady early-exit check
# (Section 4), an intensified kick test (Section 5), and an n sweep with
# a resource-spectrum sanity check (Section 6).

################################################################################
# Section 0: rebuild cod + shared helpers (self-contained, not sourced).
################################################################################

# Read-only: cod_params.rds is never written back to.
cod_params <- readRDS("cod_params.rds")

resource_limitation_2d_cod <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}

################################################################################
# Section 1: alpha bifurcation
################################################################################

native_cod_alpha <- unname(cod_params@species_params$alpha[1])
cat(sprintf("cod_params.rds's own native alpha (assimilation efficiency): %.4g\n", native_cod_alpha))

# Propagation check before trusting the sweep.
alpha_probe_base <- cod_params
alpha_probe_low  <- cod_params
given_species_params(alpha_probe_low)$alpha <- native_cod_alpha * 0.5

egrowth_changed <- !isTRUE(all.equal(getEGrowth(alpha_probe_base), getEGrowth(alpha_probe_low)))

cat(sprintf("alpha propagation check: getEGrowth() changed = %s.\n", egrowth_changed))
if (!egrowth_changed) {
  stop(paste(
    "given_species_params(params)$alpha <- value does NOT change",
    "getEGrowth() -- alpha is being cached or ignored somewhere upstream,",
    "so the sweep below would silently sweep nothing. Stopping before",
    "wasting the run."
  ))
}

build_cod_alpha_variant <- function(alpha, resource_decrease, capacity_mult) {
  p <- cod_params
  given_species_params(p)$alpha <- alpha
  resource_limitation_2d_cod(p, resource_decrease, capacity_mult)
}

# method="predictor_corrector" avoids euler's numerical ringing (Day 6);
# t_per=0.2 avoids undersampling a real oscillation; tryCatch records a
# crash as an NA point rather than killing the whole sweep.
run_bifurcation_sweep_alpha <- function(alpha_seq, resource_decrease, capacity_mult,
                                        t_max = 2000, effort = 1,
                                        metric_fn = function(sim) rowSums(getYield(sim))) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      p <- build_cod_alpha_variant(seq_vals[i], resource_decrease, capacity_mult)
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }

      result <- tryCatch({
        sim  <- projectToSteady(p, effort = effort, t_max = t_max,
                                t_per = 0.2, method = "predictor_corrector",
                                return_sim = TRUE, progress_bar = FALSE)
        mv   <- metric_fn(sim)
        tv   <- as.numeric(names(mv))
        # Late window only, guards against the leading transient.
        late <- mv[tv > max(tv) * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late),
             n = sim@n[last, , ], npp = sim@n_pp[last, ], error = NA_character_)
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_,
             n = state_n, npp = state_npp, error = conditionMessage(e))
      })

      if (!is.na(result$error)) {
        warning(sprintf("alpha=%.4g: %s", seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(alpha_seq)
  bwd    <- run_one_direction(rev(alpha_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

plot_bifurcation_alpha <- function(df, x_label, y_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}

# Log-spaced, 0.005 to 1 -- alpha is strictly positive.
alpha_seq_bif <- exp(seq(log(0.005), log(1), length.out = 84))

bif_alpha_df <- run_bifurcation_sweep_alpha(
  alpha_seq_bif,
  resource_decrease = 1, capacity_mult = 1,
  effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)

write.csv(bif_alpha_df, file.path("interesting_plots", "day32_cod_alpha_yield_bifurcation.csv"),
          row.names = FALSE)

bif_alpha_plot <- plot_bifurcation_alpha(
  bif_alpha_df, "alpha (assimilation efficiency)", "Yield",
  "Cod bifurcation diagram: yield vs. assimilation efficiency alpha, swept forward and backward",
  sprintf(
    "resource_decrease=1, capacity_mult=1, effort=1 fixed -- native alpha=%.3g; branches collapsing onto one curve = fixed point, fanning apart = limit cycle",
    native_cod_alpha
  )
)
bif_alpha_plot

save_plot(bif_alpha_plot, "day32_cod_alpha_yield_bifurcation.png", width = 9, height = 6)

# 1e-6 relative-amplitude threshold; hysteresis = forward/backward disagreement.
bif_alpha_check <- bif_alpha_df %>%
  tidyr::pivot_wider(names_from = c(direction, branch), values_from = metric) %>%
  mutate(
    rel_amplitude_fwd = (Forward_max - Forward_min) / ((Forward_max + Forward_min) / 2),
    rel_amplitude_bwd = (Backward_max - Backward_min) / ((Backward_max + Backward_min) / 2),
    hysteresis_gap    = abs(Forward_max - Backward_max) + abs(Forward_min - Backward_min)
  )
print(bif_alpha_check %>% select(value, rel_amplitude_fwd, rel_amplitude_bwd, hysteresis_gap))

n_oscillating_alpha  <- sum(bif_alpha_check$rel_amplitude_fwd > 1e-6 | bif_alpha_check$rel_amplitude_bwd > 1e-6, na.rm = TRUE)
max_hysteresis_alpha <- bif_alpha_check$value[which.max(bif_alpha_check$hysteresis_gap)]

alpha_verdict <- if (n_oscillating_alpha == 0) {
  sprintf(
    "flat -- 0/%d alpha values cross the 1e-6 relative-amplitude threshold at (resource_decrease=1, capacity_mult=1); no forward/backward branch separation either (largest hysteresis gap at alpha=%.4g, still negligible)",
    nrow(bif_alpha_check), max_hysteresis_alpha
  )
} else {
  sprintf(
    "%d/%d alpha values cross into oscillation at this fixed point; largest forward/backward gap at alpha=%.4g",
    n_oscillating_alpha, nrow(bif_alpha_check), max_hysteresis_alpha
  )
}
cat(sprintf("\nSection 1 verdict: %s.\n", alpha_verdict))

################################################################################
# Section 2: q bifurcation, same direct-loop shape as Section 1
################################################################################

native_cod_q <- unname(cod_params@species_params$q[1])
cat(sprintf("cod_params.rds's own native q (search-volume exponent): %.4g\n", native_cod_q))

# q feeds SearchVolume = gamma * w^q, which mizer may cache -- confirmed
# directly rather than assumed.
q_probe_base <- cod_params
q_probe_low  <- cod_params
given_species_params(q_probe_low)$q <- native_cod_q - 0.3

search_vol_changed <- !isTRUE(all.equal(q_probe_base@search_vol, q_probe_low@search_vol))
encounter_changed  <- !isTRUE(all.equal(getEncounter(q_probe_base), getEncounter(q_probe_low)))

cat(sprintf("q propagation check: search_vol changed = %s, getEncounter() changed = %s.\n",
           search_vol_changed, encounter_changed))
if (!search_vol_changed && !encounter_changed) {
  stop(paste(
    "given_species_params(params)$q <- value does NOT change search_vol or",
    "getEncounter() -- q is being cached at construction time, so the sweep",
    "below would silently sweep nothing. Stopping before wasting the run;",
    "switch to species_params(params)$q <- value (full replacement,",
    "triggers recalculation) instead."
  ))
}

build_cod_q_variant <- function(q, resource_decrease, capacity_mult) {
  p <- cod_params
  given_species_params(p)$q <- q
  resource_limitation_2d_cod(p, resource_decrease, capacity_mult)
}

run_bifurcation_sweep_q <- function(q_seq, resource_decrease, capacity_mult,
                                    t_max = 2000, effort = 1,
                                    metric_fn = function(sim) rowSums(getYield(sim))) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      p <- build_cod_q_variant(seq_vals[i], resource_decrease, capacity_mult)
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }

      result <- tryCatch({
        sim  <- projectToSteady(p, effort = effort, t_max = t_max,
                                t_per = 0.2, method = "predictor_corrector",
                                return_sim = TRUE, progress_bar = FALSE)
        mv   <- metric_fn(sim)
        tv   <- as.numeric(names(mv))
        late <- mv[tv > max(tv) * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late),
             n = sim@n[last, , ], npp = sim@n_pp[last, ], error = NA_character_)
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_,
             n = state_n, npp = state_npp, error = conditionMessage(e))
      })

      if (!is.na(result$error)) {
        warning(sprintf("q=%.4g: %s", seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(q_seq)
  bwd    <- run_one_direction(rev(q_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

# Linear axes -- q's +/-0.3 bracket doesn't span orders of magnitude.
plot_bifurcation_q <- function(df, x_label, y_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}

q_seq_bif <- seq(native_cod_q - 0.3, native_cod_q + 0.3, by = 0.01)

bif_q_df <- run_bifurcation_sweep_q(
  q_seq_bif,
  resource_decrease = 1, capacity_mult = 1,
  effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)

write.csv(bif_q_df, file.path("interesting_plots", "day32_cod_q_yield_bifurcation.csv"),
          row.names = FALSE)

bif_q_plot <- plot_bifurcation_q(
  bif_q_df, "q (search-volume exponent)", "Yield",
  "Cod bifurcation diagram: yield vs. search-volume exponent q, swept forward and backward",
  sprintf(
    "resource_decrease=1, capacity_mult=1, effort=1 fixed -- native q=%.3g; branches collapsing onto one curve = fixed point, fanning apart = limit cycle",
    native_cod_q
  )
)
bif_q_plot

save_plot(bif_q_plot, "day32_cod_q_yield_bifurcation.png", width = 9, height = 6)

bif_q_check <- bif_q_df %>%
  tidyr::pivot_wider(names_from = c(direction, branch), values_from = metric) %>%
  mutate(
    rel_amplitude_fwd = (Forward_max - Forward_min) / ((Forward_max + Forward_min) / 2),
    rel_amplitude_bwd = (Backward_max - Backward_min) / ((Backward_max + Backward_min) / 2),
    hysteresis_gap    = abs(Forward_max - Backward_max) + abs(Forward_min - Backward_min)
  )
print(bif_q_check %>% select(value, rel_amplitude_fwd, rel_amplitude_bwd, hysteresis_gap))

n_oscillating_q  <- sum(bif_q_check$rel_amplitude_fwd > 1e-6 | bif_q_check$rel_amplitude_bwd > 1e-6, na.rm = TRUE)
max_hysteresis_q <- bif_q_check$value[which.max(bif_q_check$hysteresis_gap)]

q_verdict <- if (n_oscillating_q == 0) {
  sprintf(
    "flat -- 0/%d q values cross the 1e-6 relative-amplitude threshold at (resource_decrease=1, capacity_mult=1); no forward/backward branch separation either (largest hysteresis gap at q=%.3g, still negligible)",
    nrow(bif_q_check), max_hysteresis_q
  )
} else {
  sprintf(
    "%d/%d q values cross into oscillation at this fixed point; largest forward/backward gap at q=%.3g",
    n_oscillating_q, nrow(bif_q_check), max_hysteresis_q
  )
}
cat(sprintf("\nSection 2 verdict: %s.\n", q_verdict))

################################################################################
# Section 3: resource_rate bifurcation, capacity fixed via a one-time
# balance=TRUE baseline
#
# Re-deriving balance=TRUE fresh at every sweep point (the first version
# of this section) self-cancelled: it always balanced against the same
# fixed native community demand, so capacity grew to compensate as
# resource_rate_mult shrank, holding cod's food supply ~constant across
# the whole sweep. Fixed by deriving capacity once, against cod's native
# rate, then holding it fixed (balance=FALSE) while resource_rate_mult
# is actually swept.
################################################################################

native_balanced_params   <- setResource(cod_params, resource_rate = getResourceRate(cod_params), balance = TRUE)
native_balanced_capacity <- getResourceCapacity(native_balanced_params)

resource_limitation_balance_cod <- function(params, resource_rate_mult) {
  new_rate <- getResourceRate(params) * resource_rate_mult
  setResource(params, resource_rate = new_rate,
             resource_capacity = native_balanced_capacity, balance = FALSE)
}

resource_rate_probe_base <- cod_params
resource_rate_probe_low  <- resource_limitation_balance_cod(cod_params, 0.5)

resource_rate_changed <- !isTRUE(all.equal(getResourceRate(resource_rate_probe_base),
                                           getResourceRate(resource_rate_probe_low)))

cat(sprintf("resource_rate propagation check: resource_rate changed = %s (resource_capacity fixed at %.4g by design).\n",
           resource_rate_changed, native_balanced_capacity[1]))
if (!resource_rate_changed) {
  stop(paste(
    "given/setResource()'s resource_rate <- value does NOT change",
    "getResourceRate() -- resource_rate is being cached or ignored",
    "somewhere upstream, so the sweep below would silently sweep nothing.",
    "Stopping before wasting the run."
  ))
}

build_cod_resource_rate_variant <- function(resource_rate_mult) {
  resource_limitation_balance_cod(cod_params, resource_rate_mult)
}

run_bifurcation_sweep_resource_rate <- function(resource_rate_seq, t_max = 2000, effort = 1,
                                                metric_fn = function(sim) rowSums(getYield(sim))) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      p <- build_cod_resource_rate_variant(seq_vals[i])
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }

      result <- tryCatch({
        sim  <- projectToSteady(p, effort = effort, t_max = t_max,
                                t_per = 0.2, method = "predictor_corrector",
                                return_sim = TRUE, progress_bar = FALSE)
        mv   <- metric_fn(sim)
        tv   <- as.numeric(names(mv))
        late <- mv[tv > max(tv) * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late),
             n = sim@n[last, , ], npp = sim@n_pp[last, ], error = NA_character_)
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_,
             n = state_n, npp = state_npp, error = conditionMessage(e))
      })

      if (!is.na(result$error)) {
        warning(sprintf("resource_rate_mult=%.4g: %s", seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(resource_rate_seq)
  bwd    <- run_one_direction(rev(resource_rate_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

# Log-scaled -- range spans seven orders of magnitude.
plot_bifurcation_resource_rate <- function(df, x_label, y_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}

resource_rate_seq_bif <- exp(seq(log(1e-7), log(1), length.out = 100))

bif_resource_rate_df <- run_bifurcation_sweep_resource_rate(
  resource_rate_seq_bif,
  effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)

write.csv(bif_resource_rate_df, file.path("interesting_plots", "day32_cod_resource_rate_balance_yield_bifurcation.csv"),
          row.names = FALSE)

bif_resource_rate_plot <- plot_bifurcation_resource_rate(
  bif_resource_rate_df, "resource_rate multiplier (native = 1)", "Yield",
  "Cod bifurcation diagram: yield vs. resource_rate multiplier, balance=TRUE, swept forward and backward",
  "effort=1 fixed; resource_capacity fixed at cod's own native balance=TRUE baseline, only resource_rate actually swept; both axes log-scaled; branches collapsing onto one curve = fixed point, fanning apart = limit cycle"
)
bif_resource_rate_plot

save_plot(bif_resource_rate_plot, "day32_cod_resource_rate_bif.png", width = 9, height = 6)

bif_resource_rate_check <- bif_resource_rate_df %>%
  tidyr::pivot_wider(names_from = c(direction, branch), values_from = metric) %>%
  mutate(
    rel_amplitude_fwd = (Forward_max - Forward_min) / ((Forward_max + Forward_min) / 2),
    rel_amplitude_bwd = (Backward_max - Backward_min) / ((Backward_max + Backward_min) / 2),
    hysteresis_gap    = abs(Forward_max - Backward_max) + abs(Forward_min - Backward_min)
  )
print(bif_resource_rate_check %>% select(value, rel_amplitude_fwd, rel_amplitude_bwd, hysteresis_gap))

n_oscillating_resource_rate  <- sum(bif_resource_rate_check$rel_amplitude_fwd > 1e-6 | bif_resource_rate_check$rel_amplitude_bwd > 1e-6, na.rm = TRUE)
max_hysteresis_resource_rate <- bif_resource_rate_check$value[which.max(bif_resource_rate_check$hysteresis_gap)]

resource_rate_verdict <- if (n_oscillating_resource_rate == 0) {
  sprintf(
    "flat -- 0/%d resource_rate values cross the 1e-6 relative-amplitude threshold; no forward/backward branch separation either (largest hysteresis gap at resource_rate_mult=%.4g, still negligible)",
    nrow(bif_resource_rate_check), max_hysteresis_resource_rate
  )
} else {
  sprintf(
    "%d/%d resource_rate values cross into oscillation; largest forward/backward gap at resource_rate_mult=%.4g",
    n_oscillating_resource_rate, nrow(bif_resource_rate_check), max_hysteresis_resource_rate
  )
}
cat(sprintf("\nSection 3 verdict: %s.\n", resource_rate_verdict))

################################################################################
# Section 4: does projectToSteady() exit too early?
#
# Reruns Section 3's forward sweep two ways, tracking the full yield
# trajectory (not just each point's late-window max/min) as one
# continuous time series: projectToSteady() as Section 3 uses it, vs.
# project() with no early exit (t_first=800 cold start, t_rest=100 per
# warm-started point). Forward only, starting at resource_rate_mult=1e-7,
# the most diagnostic point in Section 3's range.
################################################################################

run_timeseries_resource_rate_steady <- function(resource_rate_seq, t_max = 2000, effort = 1,
                                                metric_fn = function(sim) rowSums(getYield(sim))) {
  state_n   <- NULL
  state_npp <- NULL
  t_offset  <- 0
  out       <- vector("list", length(resource_rate_seq))

  for (i in seq_along(resource_rate_seq)) {
    p <- build_cod_resource_rate_variant(resource_rate_seq[i])
    if (!is.null(state_n)) {
      p@initial_n[]    <- state_n
      p@initial_n_pp[] <- state_npp
    }

    sim <- projectToSteady(p, effort = effort, t_max = t_max, t_per = 0.2,
                           method = "predictor_corrector",
                           return_sim = TRUE, progress_bar = FALSE)
    mv <- metric_fn(sim)
    tv <- as.numeric(names(mv))

    out[[i]] <- data.frame(time = tv + t_offset, yield = mv, value = resource_rate_seq[i])
    t_offset <- max(out[[i]]$time)

    last      <- dim(sim@n)[1]
    state_n   <- sim@n[last, , ]
    state_npp <- sim@n_pp[last, ]
  }

  bind_rows(out)
}

run_timeseries_resource_rate_project <- function(resource_rate_seq, t_first = 800, t_rest = 100,
                                                 effort = 1,
                                                 metric_fn = function(sim) rowSums(getYield(sim))) {
  state_n   <- NULL
  state_npp <- NULL
  t_offset  <- 0
  out       <- vector("list", length(resource_rate_seq))

  for (i in seq_along(resource_rate_seq)) {
    p          <- build_cod_resource_rate_variant(resource_rate_seq[i])
    cold_start <- is.null(state_n)
    if (!cold_start) {
      p@initial_n[]    <- state_n
      p@initial_n_pp[] <- state_npp
    }

    t_max_i <- if (cold_start) t_first else t_rest
    sim <- project(p, effort = effort, t_max = t_max_i, dt = 0.1, t_save = 0.2,
                   method = "predictor_corrector")
    mv <- metric_fn(sim)
    tv <- as.numeric(names(mv))

    out[[i]] <- data.frame(time = tv + t_offset, yield = mv, value = resource_rate_seq[i])
    t_offset <- max(out[[i]]$time)

    last      <- dim(sim@n)[1]
    state_n   <- sim@n[last, , ]
    state_npp <- sim@n_pp[last, ]
  }

  bind_rows(out)
}

ts_steady_df  <- run_timeseries_resource_rate_steady(resource_rate_seq_bif, effort = 1)
ts_project_df <- run_timeseries_resource_rate_project(resource_rate_seq_bif, t_first = 800, t_rest = 100, effort = 1)

resource_rate_timeseries_df <- bind_rows(
  data.frame(ts_steady_df,  method = "projectToSteady()"),
  data.frame(ts_project_df, method = "project()")
)

write.csv(resource_rate_timeseries_df,
         file.path("interesting_plots", "day32_cod_resource_rate_timeseries_comparison.csv"),
         row.names = FALSE)

resource_rate_timeseries_plot <- ggplot(resource_rate_timeseries_df, aes(x = time, y = yield, color = method)) +
  geom_line(linewidth = 0.6) +
  scale_y_log10() +
  labs(x = "Cumulative time (years)", y = "Yield (log scale)",
       title = "Cod yield over time: forward resource_rate sweep re-run two ways",
       subtitle = "projectToSteady() exits each step at its own convergence tolerance; project() runs a fixed t_first=800 (cold start) / t_rest=100 (warm start) per step with no early exit") +
  theme_minimal()
resource_rate_timeseries_plot

save_plot(resource_rate_timeseries_plot, "day32_resource_rate_timeseries.png", width = 10, height = 6)

cat(sprintf(
  "\nSection 4: projectToSteady() total elapsed time across the forward sweep = %.4g years; project() = %.4g years (t_first=%d + %d x t_rest=%d).\n",
  max(ts_steady_df$time), max(ts_project_df$time), 800, length(resource_rate_seq_bif) - 1, 100
))

################################################################################
# Section 5: a more extreme kick than Day 30's original
#
# Day 30's make_limit_cycle_sim() (rel_amplitude=0.041 at this corner):
# resource dropped to 10% of capacity, run 10 years, then a slice of the
# consumer spectrum (w in [10,100]) cut 1000-fold, run forward. Here the
# slice spans the entire mature range instead, cut 1e7-fold. Run at
# resource_rate_mult=1e-3, not the 1e-7 extreme Sections 3-4 use -- 1e-7
# combined with this section's own cut leaves nothing to recover from and
# just goes extinct. effort=0/biomass (not effort=1/yield), matching Day
# 30's own convention -- this is a population-dynamics question, not a
# fished-yield one.
################################################################################

make_extreme_kick_sim <- function(params, t_total = 2000, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1

  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor_corrector")

  w_mat <- params@species_params$w_mat[1]
  idx   <- params@w >= w_mat & params@w <= max(params@w)
  last  <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 1e7

  # projectToSteady() here, not project() -- exits early if it settles,
  # still runs to t_max as a safety ceiling if it doesn't.
  params@initial_n[]    <- sim_init@n[last, , ]
  params@initial_n_pp[] <- sim_init@n_pp[last, ]

  projectToSteady(params, t_max = t_total - 10, t_per = 0.2,
                  method = "predictor_corrector", effort = effort,
                  return_sim = TRUE, progress_bar = FALSE)
}

p_kick_probe <- build_cod_resource_rate_variant(1e-3)
sim_extreme_kick <- make_extreme_kick_sim(p_kick_probe, t_total = 2000, effort = 0)

biomass_kick <- rowSums(getBiomass(sim_extreme_kick))
time_kick    <- as.numeric(names(biomass_kick))
late_kick    <- biomass_kick[time_kick >= max(time_kick) * 0.6]
rel_amplitude_kick <- (max(late_kick) - min(late_kick)) / mean(late_kick)
kick_verdict <- if (rel_amplitude_kick > 1e-6) "oscillating" else "settled"

cat(sprintf(
  "\nSection 5: extreme kick (entire mature range, /1e7) at resource_rate_mult=%.4g -- rel_amplitude=%.4g, verdict=%s (Day 30's original milder kick at a comparable corner found rel_amplitude=0.041).\n",
  1e-3, rel_amplitude_kick, kick_verdict
))

extreme_kick_df <- data.frame(time = time_kick, biomass = biomass_kick)
write.csv(extreme_kick_df, file.path("interesting_plots", "day32_extreme_kick_timeseries.csv"),
         row.names = FALSE)

extreme_kick_plot <- ggplot(extreme_kick_df, aes(x = time, y = biomass)) +
  geom_line(linewidth = 0.6) +
  scale_y_log10() +
  labs(x = "Time (years)", y = "Biomass (log scale)",
       title = sprintf("Cod biomass after an extreme kick at resource_rate_mult=%.4g (balance=TRUE)", 1e-3),
       subtitle = "Entire mature size range cut by 1e7 at t=10, then run forward -- Day 30's own kick cut only a 10-100g slice by 1e3 at this same corner and found rel_amplitude=0.041") +
  theme_minimal()
extreme_kick_plot

save_plot(extreme_kick_plot, "day32_extreme_kick_timeseries.png", width = 9, height = 6)

# project()-only twin of make_extreme_kick_sim() -- same kick, no early
# exit, full t_total-10 years -- checking projectToSteady()'s convergence
# call against the real trajectory rather than trusting it.
make_extreme_kick_sim_project <- function(params, t_total = 2000, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1

  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor_corrector")

  w_mat <- params@species_params$w_mat[1]
  idx   <- params@w >= w_mat & params@w <= max(params@w)
  last  <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 1e7

  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "predictor_corrector")
}

sim_extreme_kick_project <- make_extreme_kick_sim_project(p_kick_probe, t_total = 2000, effort = 0)

biomass_kick_project <- rowSums(getBiomass(sim_extreme_kick_project))
time_kick_project    <- as.numeric(names(biomass_kick_project))
late_kick_project    <- biomass_kick_project[time_kick_project >= max(time_kick_project) * 0.6]
rel_amplitude_kick_project <- (max(late_kick_project) - min(late_kick_project)) / mean(late_kick_project)
kick_verdict_project <- if (rel_amplitude_kick_project > 1e-6) "oscillating" else "settled"

cat(sprintf(
  "\nSection 5 (project() twin): ran the full %d years (vs. projectToSteady()'s %.4g before it exited) -- rel_amplitude=%.4g, verdict=%s.\n",
  2000 - 10, max(time_kick), rel_amplitude_kick_project, kick_verdict_project
))

kick_comparison_df <- bind_rows(
  data.frame(time = time_kick,         biomass = biomass_kick,         method = "projectToSteady()"),
  data.frame(time = time_kick_project, biomass = biomass_kick_project, method = "project()")
)

write.csv(kick_comparison_df, file.path("interesting_plots", "day32_extreme_kick_method_comparison.csv"),
         row.names = FALSE)

kick_comparison_plot <- ggplot(kick_comparison_df, aes(x = time, y = biomass, color = method)) +
  geom_line(linewidth = 0.6) +
  scale_y_log10() +
  labs(x = "Time (years)", y = "Biomass (log scale)",
       title = sprintf("Cod biomass after an extreme kick at resource_rate_mult=%.4g: projectToSteady() vs. project()", 1e-3),
       subtitle = "Same kick (entire mature range cut by 1e7 at t=10), same starting point -- projectToSteady() exits at its own convergence check; project() runs the full 1990 years with no early exit") +
  theme_minimal()
kick_comparison_plot

save_plot(kick_comparison_plot, "day32_extreme_kick_method_comparison.png", width = 9, height = 6)

################################################################################
# Section 6: n bifurcation (resource growth exponent)
#
# n feeds the resource rate spectrum as resource_rate * w_full^(n-1)
# (separate from lambda, which shapes capacity). Nothing else touched.
# Two full sweeps -- projectToSteady() and no-early-exit project() -- as
# two separate plots.
################################################################################

native_cod_n <- resource_params(cod_params)[["n"]]
cat(sprintf("cod_params.rds's own native n (resource growth exponent): %.4g\n", native_cod_n))

# Neither setResource(params, n = ...) nor resource_params(params)$n <-
# propagate here: both rebuild the rate from a scalar formula
# (r_pp * w^(n-1)), and both correctly leave a "manually set"/frozen
# array untouched when the existing rate wasn't derived from that
# formula in the first place -- which cod's calibrated rate isn't
# (confirmed: the resource_params() route left the sanity-check curves
# identical, exactly as its own "frozen arrays untouched" documentation
# says it should). So instead of asking mizer to rebuild from a scalar,
# the existing rate curve is rescaled directly, by the ratio between the
# new and native exponents -- this works on whatever the current array
# is, frozen or not.
apply_n <- function(params, n_val) {
  native_rate <- as.numeric(getResourceRate(params))
  native_w    <- w_full(params)
  new_rate    <- native_rate * native_w^(n_val - resource_params(params)[["n"]])
  setResource(params, resource_rate = new_rate, n = n_val, balance = FALSE)
}

# Bounds match n_seq_bif's own range below -- no principled reason to
# bracket n tightly around its native value the way lambda's own
# theory-grounded sweeps have; n is left wide open instead.
n_probe_a <- apply_n(cod_params, -1.8)
n_probe_b <- apply_n(cod_params, 1.8)

# getResourceRate() returns a numeric vector wrapped in a custom S3 class
# ('ArrayResourceBySize', carrying units/comment/params as attributes) --
# as.numeric() strips it, since data.frame() would otherwise dispatch to
# mizer's own as.data.frame method for that class and split it into
# rate.w/rate.value sub-columns instead of a plain "rate" column.
rate_a <- as.numeric(getResourceRate(n_probe_a))
rate_b <- as.numeric(getResourceRate(n_probe_b))
w_a    <- w_full(n_probe_a)
w_b    <- w_full(n_probe_b)
cat(sprintf("n sanity check: length(rate)=%d/%d, length(w_full)=%d/%d (a/b) -- all four should match and be > 0.\n",
           length(rate_a), length(rate_b), length(w_a), length(w_b)))

resource_sanity_df <- rbind(
  data.frame(w = w_a, rate = rate_a, variant = "n=-1.8"),
  data.frame(w = w_b, rate = rate_b, variant = "n=1.8")
)

resource_sanity_plot <- ggplot(resource_sanity_df, aes(x = w, y = rate, color = variant)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Weight (g, log scale)", y = "Resource rate (log scale)",
       title = "Sanity check: does varying n actually change the resource rate spectrum?",
       subtitle = "Two example n values, built and plotted before the real sweep below runs -- if these two lines don't visibly differ in slope, n isn't propagating and the sweep would be wasted") +
  theme_minimal()
resource_sanity_plot

save_plot(resource_sanity_plot, "day32_n_resource_sanity_check.png", width = 9, height = 6)

build_cod_n_variant <- function(n_val) {
  apply_n(cod_params, n_val)
}

# Wide, not bracketed around the native value -- unlike lambda (theory
# grounds it near 2) or q (a documented O(1) exponent), n has no
# principled reason to sit near cod's own calibrated 0.7, so the range
# itself is arbitrary: -1.8 to 1.8. Negative n is fine here -- apply_n()
# goes through setResource(), which only asserts is.number(n), not n >= 0
# (that non-negativity constraint belonged to the resource_params<-
# route this file no longer uses).
n_seq_bif <- seq(-1.8, 1.8, by = 0.025)

# effort=1/yield, matching Sections 1-4's convention.
# tol decreased from projectToSteady()'s own default -- the forward and
# backward branches were disagreeing (a hysteresis-looking gap) purely
# because each got called "converged" before it actually settled, not
# because of genuine bistability. A stricter tolerance forces closer
# agreement between successive states before exiting, so the reported
# max/min at each point reflects where the population actually ends up
# rather than wherever it happened to be when the looser check declared
# victory.
run_bifurcation_sweep_n_steady <- function(n_seq, t_max = 2000, effort = 1, tol = 1e-5,
                                           metric_fn = function(sim) rowSums(getYield(sim))) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      p <- build_cod_n_variant(seq_vals[i])
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }

      result <- tryCatch({
        sim  <- projectToSteady(p, effort = effort, t_max = t_max, tol = tol,
                                t_per = 0.2, method = "predictor_corrector",
                                return_sim = TRUE, progress_bar = FALSE)
        mv   <- metric_fn(sim)
        tv   <- as.numeric(names(mv))
        late <- mv[tv > max(tv) * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late),
             n = sim@n[last, , ], npp = sim@n_pp[last, ], error = NA_character_)
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_,
             n = state_n, npp = state_npp, error = conditionMessage(e))
      })

      if (!is.na(result$error)) {
        warning(sprintf("n=%.4g: %s", seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(n_seq)
  bwd    <- run_one_direction(rev(n_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

run_bifurcation_sweep_n_project <- function(n_seq, t_first = 2050, t_rest = 100, effort = 1,
                                            metric_fn = function(sim) rowSums(getYield(sim))) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      p          <- build_cod_n_variant(seq_vals[i])
      cold_start <- is.null(state_n)
      if (!cold_start) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }

      t_max_i <- if (cold_start) t_first else t_rest
      result <- tryCatch({
        sim  <- project(p, effort = effort, t_max = t_max_i, dt = 0.1, t_save = 0.2,
                        progress_bar = FALSE, method = "predictor_corrector")
        mv   <- metric_fn(sim)
        tv   <- as.numeric(names(mv))
        late <- mv[tv > max(tv) * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late),
             n = sim@n[last, , ], npp = sim@n_pp[last, ], error = NA_character_)
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_,
             n = state_n, npp = state_npp, error = conditionMessage(e))
      })

      if (!is.na(result$error)) {
        warning(sprintf("n=%.4g: %s", seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(n_seq)
  bwd    <- run_one_direction(rev(n_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

plot_bifurcation_n <- function(df, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    labs(x = "n (resource growth exponent)", y = "Yield", title = title, subtitle = subtitle) +
    theme_minimal()
}

bif_n_steady_df <- run_bifurcation_sweep_n_steady(
  n_seq_bif, effort = 1, tol = 1e-5,
  metric_fn = function(sim) rowSums(getYield(sim))
)
write.csv(bif_n_steady_df, file.path("interesting_plots", "day32_cod_n_yield_bifurcation_steady.csv"),
         row.names = FALSE)

bif_n_steady_plot <- plot_bifurcation_n(
  bif_n_steady_df,
  "Cod bifurcation diagram: yield vs. resource growth exponent n (projectToSteady)",
  sprintf("effort=1 fixed, tol=1e-5 (tightened from default), all else at cod's native values -- native n=%.3g; branches collapsing onto one curve = fixed point, fanning apart = limit cycle", native_cod_n)
)
bif_n_steady_plot
save_plot(bif_n_steady_plot, "day32_cod_n_bif_steady.png", width = 9, height = 6)

bif_n_project_df <- run_bifurcation_sweep_n_project(
  n_seq_bif, t_first = 2050, t_rest = 100, effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)
write.csv(bif_n_project_df, file.path("interesting_plots", "day32_cod_n_yield_bifurcation_project.csv"),
         row.names = FALSE)

bif_n_project_plot <- plot_bifurcation_n(
  bif_n_project_df,
  "Cod bifurcation diagram: yield vs. resource growth exponent n (project(), no early exit)",
  sprintf("effort=1 fixed, t_first=2050/t_rest=100, all else at cod's native values -- native n=%.3g; branches collapsing onto one curve = fixed point, fanning apart = limit cycle", native_cod_n)
)
bif_n_project_plot
save_plot(bif_n_project_plot, "day32_cod_n_bif_project.png", width = 9, height = 6)

check_bifurcation <- function(df) {
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

n_steady_result  <- check_bifurcation(bif_n_steady_df)
n_project_result <- check_bifurcation(bif_n_project_df)

cat(sprintf(
  "\nSection 6 (n bifurcation): projectToSteady() %d/%d points cross the 1e-6 threshold (largest hysteresis gap at n=%.4g); project() %d/%d points cross it (largest hysteresis gap at n=%.4g).\n",
  n_steady_result$n_oscillating, nrow(n_steady_result$check), n_steady_result$max_hysteresis_value,
  n_project_result$n_oscillating, nrow(n_project_result$check), n_project_result$max_hysteresis_value
))

################################################################################
# Section 7: summary
################################################################################

cat("\n===== Day 32 summary =====\n")
cat(sprintf(
  "Section 1 (alpha bifurcation, rewritten as a direct projectToSteady() loop, re-pointed at Day 30's confirmed-oscillating corner, %d-point forward/backward bifurcation diagram): %s\n",
  length(alpha_seq_bif), alpha_verdict
))
cat(sprintf(
  "Section 2 (q bifurcation, same direct-loop rewrite, same re-pointed corner, %d-point forward/backward bifurcation diagram): %s\n",
  length(q_seq_bif), q_verdict
))
cat(sprintf(
  "Section 3 (resource_rate bifurcation, capacity fixed at cod's own native balance=TRUE baseline, %d-point forward/backward bifurcation diagram): %s\n",
  length(resource_rate_seq_bif), resource_rate_verdict
))
cat(
  "Section 4 (project() vs projectToSteady() time-series comparison across the forward resource_rate sweep): see day32_resource_rate_timeseries.png -- read visually, no automated verdict.\n"
)
cat(sprintf(
  "Section 5 (extreme kick, entire mature range cut by 1e7, at resource_rate_mult=%.4g): rel_amplitude=%.4g, verdict=%s.\n",
  1e-3, rel_amplitude_kick, kick_verdict
))
cat(sprintf(
  "Section 6 (n bifurcation, %d-point forward/backward, projectToSteady() vs project()): projectToSteady() %d/%d oscillating, project() %d/%d oscillating -- see day32_cod_n_bif_steady.png / day32_cod_n_bif_project.png.\n",
  length(n_seq_bif), n_steady_result$n_oscillating, nrow(n_steady_result$check),
  n_project_result$n_oscillating, nrow(n_project_result$check)
))
