library(mizer)
library(mizerExperimental)
library(dplyr)
library(ggplot2)
library(scales)

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH (260 chars) truncation guard -- unchanged since
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

################################################################################
# Day 32: Day 31 Section 1b's alpha bifurcation ran through the generic
# run_bifurcation_sweep_cod(param_seq, param_name, params_fn = ...)
# dispatcher -- the same function the q sweep in Section 1 used, just called
# with params_fn = build_cod_alpha_variant instead of the default
# build_cod_q_variant. That indirection is correct but one step removed from
# how the q sweep actually reads on the page: q's own bifurcation loop is
# what run_bifurcation_sweep_cod() was written FROM (Day 23's Follow-up 5),
# so it's a direct, un-indirected projectToSteady() loop hardcoded to the
# parameter being swept. This rewrites alpha's sweep the same direct way --
# no param_name string, no params_fn argument, no fixed_params list --
# build_cod_alpha_variant() called explicitly inside the loop, same as the
# q sweep would look with q instead of alpha substituted in by hand. The
# projectToSteady() call itself (method, t_per, tryCatch convention) is
# unchanged from Day 31, since that part was already right.
################################################################################

################################################################################
# Section 0: rebuild cod + shared helpers (self-contained convention, same
# as every script since Day 20 -- redefined here rather than sourced).
################################################################################

# Read-only: cod_params.rds is never written back to.
cod_params <- readRDS("cod_params.rds")

resource_limitation_2d_cod <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}

################################################################################
# Section 1: alpha bifurcation, rewritten as a direct projectToSteady() loop
#
# alpha is mizer's assimilation efficiency (the fraction of consumed energy
# that's actually assimilated, before growth/reproduction allocation) --
# see Day 31 Section 1b for the full background on why alpha specifically.
################################################################################

native_cod_alpha <- unname(cod_params@species_params$alpha[1])
cat(sprintf("cod_params.rds's own native alpha (assimilation efficiency): %.4g\n", native_cod_alpha))

# Same propagation caution as every given_species_params() write in this
# project (unresolved silent no-ops before, e.g. D_ext/capacity_mult) --
# confirmed directly rather than assumed. getEGrowth() is alpha's own
# natural analogue of getEncounter(): energy available for growth is
# alpha * encounter, net of metabolic costs, so it has to move if alpha
# actually propagated.
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

# alpha applied straight to the real cod_params object -- same pattern as
# Day 31's build_cod_q_variant()/build_cod_alpha_variant().
build_cod_alpha_variant <- function(alpha, resource_decrease, capacity_mult) {
  p <- cod_params
  given_species_params(p)$alpha <- alpha
  resource_limitation_2d_cod(p, resource_decrease, capacity_mult)
}

# Direct port of Day 31 Section 1's sweep loop (itself Day 23's Follow-up 5
# generalised, un-generalised back down) -- hardcoded to alpha instead of
# routed through run_bifurcation_sweep_cod()'s param_name/params_fn
# indirection. Every setting below is unchanged from Day 31: method is
# explicit ("predictor_corrector", not projectToSteady()'s own "euler"
# default -- Day 6 showed euler produces numerical ringing that can
# masquerade as a real oscillation), t_per is explicit (0.2, matching this
# project's t_save elsewhere, not the 1.5-year default that would
# undersample a real oscillation), and tryCatch records a crash as an
# NA/error point rather than taking down the whole sweep (Day 17's
# steady()/projectToSteady() crash under an earlier mizer build, not
# independently re-confirmed fixed since PR #452).
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
        # Late window, not the whole returned trajectory -- guards against
        # a leading transient in the very first chunk; projectToSteady()
        # has already made the "has it settled" call itself (via `tol`)
        # before this window is even taken.
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
    scale_x_log10() +
    scale_y_log10() +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}

# Fixed absolute bracket, log-spaced (0.05 to 0.9) rather than Day 31's
# proportional 0.5x-1.5x-native bracket -- still log-spaced since alpha is
# strictly positive (an additive bracket, q's own +/-0.3, risks pushing it
# negative or barely moving it), but the bounds themselves are now explicit
# rather than centred on cod's native alpha (0.6), giving a wider low-end
# range down to near-zero assimilation efficiency.
alpha_seq_bif <- exp(seq(log(0.005), log(1), length.out = 84))
# Same fixed point and fishing intensity as Day 31's q sweep, so the two
# bifurcation diagrams stay directly comparable.
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
    # Stated directly rather than "same fixed point as Day 31's q sweep" --
    # Day 31's own q-sweep subtitle claims resource_decrease=1e-7,
    # capacity_mult=1000, even though Day 31's actual run_bifurcation_sweep_cod()
    # call used fixed_params = list(resource_decrease = 1, capacity_mult = 1);
    # deferring to "same as Day 31" is ambiguous given that pre-existing
    # contradiction, so the values actually used here are given explicitly.
    "resource_decrease=1, capacity_mult=1, effort=1 fixed; native alpha=%.3g; both axes log-scaled; branches collapsing onto one curve = fixed point, fanning apart = limit cycle",
    native_cod_alpha
  )
)
bif_alpha_plot

save_plot(bif_alpha_plot, "day32_cod_alpha_yield_bifurcation.png", width = 9, height = 6)

# Same verdict logic as Day 31's bif_q_check/bif_alpha_check: is the max/min
# spread big enough to call a limit cycle (1e-6 relative-amplitude
# threshold, same as every sweep in this project), and do the
# forward/backward branches disagree (hysteresis) rather than retrace each
# other.
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
# Section 2: summary -- programmatic readout, not asserted conclusions.
################################################################################

cat("\n===== Day 32 summary =====\n")
cat(sprintf(
  "Section 1 (alpha bifurcation, rewritten as a direct projectToSteady() loop matching Day 31's q sweep structure, %d-point forward/backward bifurcation diagram): %s\n",
  length(alpha_seq_bif), alpha_verdict
))
