library(mizer)
library(mizerExperimental)
library(dplyr)
library(ggplot2)
library(scales)

# No future/future.apply here, unlike every grid-based script since Day 28
# -- the bifurcation sweep in Section 1 carries state from each step into
# the next, so its points are inherently sequential and can't be farmed
# out to parallel workers the way an independent grid can.

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
# Day 30 closed with two threads explicitly left as methodology-only:
#
#   1. The q sweep (Section 4 there) had its grid, its build function, and
#      its propagation safety-check written, but the 7x120-point grid
#      itself hadn't finished running -- no heatmap, no verdict on whether
#      mizer's own search-volume exponent produces a threshold the way
#      lambda did. Section 1 below finishes that thread, but as a classic
#      forward/backward bifurcation diagram (the format every
#      resource_decrease/capacity_mult sweep has used since Day 18) rather
#      than a perturbed-amplitude heatmap -- q swept directly with min/max
#      late-window yield plotted, not a grid of oscillating/settled
#      verdicts.
#
#   2. The corrected competitiveness metric -- dw-weighting
#      mean_mass_specific_rate() instead of an unweighted mean() over
#      mizer's unequal-width bins -- had the fix written and a single
#      before/after ratio printed, but no explicit final call on whether
#      cod and anchovy actually flip from adult- to juvenile-driven (or
#      stay put) once corrected. Section 2 below carries that through to
#      an explicit verdict for both species.
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
# Section 1: finish the q sweep -- as a classic forward/backward bifurcation
# diagram, not a perturbed-amplitude heatmap
#
# NOT de Roos & Persson's q (the juvenile/adult competitiveness ratio
# Section 2 below deals with) -- this is mizer's own species_params column,
# SearchVolume = gamma * w^q, consumed directly by getEncounter(). q is a
# genuine column on the real, fully-calibrated cod_params object, so it's
# applied straight to it via given_species_params() -- no
# newSingleSpeciesParams() reconstruction, no partial-calibration caveat.
# Mortality, erepro, R_max, gear, and cod's own native lambda all stay
# exactly as fitted; only q changes.
#
# Every earlier q pass asked "does this (resource_decrease, capacity_mult,
# q) point oscillate under an artificial kick" -- a grid of binary
# verdicts. This instead reuses the exact bifurcation-diagram method every
# resource_decrease/capacity_mult sweep has used since Day 18 (Day 23's
# Follow-up 5 generalised it into run_bifurcation_sweep()/
# plot_bifurcation()): sweep q itself forward then backward, carrying each
# step's own converged state into the next rather than re-perturbing from
# scratch each time, and plot max/min late-window yield per step.
# Branches collapsing onto a single curve = fixed point; branches fanning
# apart into separate max/min lines = a genuine limit cycle; forward and
# backward branches disagreeing at the same q = hysteresis, not just a
# q-dependent threshold.
#
# resource_decrease and capacity_mult held fixed at (1e-7, 1000) -- Day 30's
# broad scan found native cod_params' own largest oscillation
# (rel_amplitude=0.041) at resource_decrease=1e-7, capacity_mult=1
# (extreme resource starvation, at cod's native carrying capacity); this
# pushes capacity_mult an order of magnitude past the 100x that scan's own
# widest range covered, giving q the best chance of revealing similar
# structure in a 1D slice rather than picking an arbitrary point that's
# null for every mechanism tried so far.
################################################################################

native_cod_q <- unname(cod_params@species_params$q[1])
cat(sprintf("cod_params.rds's own native q (search-volume exponent): %.4g\n", native_cod_q))

# given_species_params() has been a silent no-op before in this project --
# D_ext gets attached but never consumed without a follow-up
# setExtDiffusion() call (unresolved since Day 20), and Day 29 separately
# found capacity_mult sitting unused in an early make_anchovy_params()
# draft. q feeds into SearchVolume, which mizer may cache at construction
# time rather than recompute live -- so before trusting the sweep below,
# confirm directly that search_vol/getEncounter() actually move when q
# changes.
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

# q applied straight to the real cod_params object -- see the header note.
build_cod_q_variant <- function(q, resource_decrease, capacity_mult) {
  p <- cod_params
  given_species_params(p)$q <- q
  resource_limitation_2d_cod(p, resource_decrease, capacity_mult)
}

cod_params_low <- build_cod_q_variant(0.5,1,1)
cod_steady_params_low <- projectToSteady(cod_q_params_low,effort=1,t_max=800)
params <- cod_steady_params_low
given_species_params(params)$q <- 0.51

cod_steady_params_olow <- projectToSteady(params,effort=1,t_max=800)


cod_steady_params_high <- projectToSteady(cod_q_params_high,effort=1,t_max=800)
plotSpectra(cod_steady_params_high,log="xy")
plotSpectra(cod_steady_params_low,log="xy")
plotHover(getEGrowth(cod_steady_params_low),log="xy")
plotHover(getEGrowth(cod_steady_params_high),log="xy")

plot2(search_vol(cod_steady_params_high),search_vol(cod_steady_params_low),log="xy")

given_species_params(cod_steady_params_high)$q
given_species_params(cod_steady_params_low)$q
getYield(cod_steady_params_high)
getYield(cod_steady_params_low)
attr(cod_steady_params_low,"convergence")
# Generic forward/backward bifurcation sweep -- direct port of Day 23's
# run_bifurcation_sweep() (Follow-up 5 there) onto cod's own convention,
# now built on projectToSteady() rather than project(). Every version of
# this sweep before today ran project() for a fixed t_max and just hoped
# that was long enough for the system to settle -- t_run/t_run_first were
# exactly that guess, made explicit rather than removed. projectToSteady()
# removes the guess itself: it keeps projecting forward in chunks and
# checking its own convergence tolerance (consecutive states agreeing to
# within `tol`, mizer's own distanceSSLogN by default) between them,
# stopping the moment they do rather than at some fixed clock time picked
# in advance. t_max below is now just a generous safety ceiling for a run
# that never converges at all -- which is exactly what a genuine limit
# cycle looks like, since a truly oscillating point has no fixed point to
# settle onto -- not a per-point estimate of settling time.
#
# If a point simply fails to converge by t_max, projectToSteady() surfaces
# that as a WARNING, not an error, and still returns whatever trajectory
# it ran -- which is exactly the case this sweep wants to see, since a
# non-converged, still-oscillating trajectory is what makes the late-
# window max/min genuinely split apart below. tryCatch here is for a
# different, narrower purpose: Day 17's own notes
# (day_17_experiments.R) flag a real crash in
# steady()/projectToSteady() under an earlier mizer build (a missing
# r$rdd entry breaking project_n_no_diffusion), which forced project()
# to be used directly back then. That was under mizer 3.1.0.9000, before
# PR #452 (getStability()/steadyNewton(), installed Day 28) updated the
# package, and hasn't been independently re-confirmed fixed since -- so a
# crash here is caught and recorded as an NA/error point (the same
# convention every fallible sweep point in this project has used since
# Day 25) rather than taking down the whole sweep.
#
# metric_fn/effort generalise this beyond biomass: getYield() is
# identically zero under effort=0 (Day 23's own convention, since those
# sweeps were about stability, not fishing), so plotting yield's min/max
# requires actually fishing the population -- effort is a parameter here,
# not hardcoded, so the caller has to make that choice explicitly rather
# than the sweep silently reusing a setting tuned for a different metric.
run_bifurcation_sweep_cod <- function(param_seq, param_name, fixed_params = list(),
                                      params_fn = build_cod_q_variant, t_max = 2000,
                                      effort = 0,
                                      metric_fn = function(sim) rowSums(getBiomass(sim))) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      args <- fixed_params
      args[[param_name]] <- seq_vals[i]
      p <- do.call(params_fn, args)
      if (!is.null(state_n)) {
        p@initial_n[]    <- state_n
        p@initial_n_pp[] <- state_npp
      }

      result <- tryCatch({
        # method is explicit and NOT the default: projectToSteady()'s own
        # default is "euler", the exact integrator Day 6 showed produces
        # numerical ringing that can masquerade as a real oscillation --
        # every project() call in this project has used
        # "predictor-corrector" since that finding, and projectToSteady()
        # needs the same override (spelled with an underscore here,
        # unlike project()'s own hyphenated method names) or it would
        # silently undo that six-day-old lesson. t_per is also explicit
        # (0.2, not the 1.5-year default): t_per sets both how often
        # convergence is checked AND the spacing of points in the
        # returned sim, and this project's own project() calls elsewhere
        # in this file all save at t_save=0.2 -- leaving t_per at 1.5
        # would undersample any real oscillation with a period shorter
        # than that, and make the late-window max/min below unreliable.
        sim  <- projectToSteady(p, effort = effort, t_max = t_max,
                                t_per = 0.2, method = "predictor_corrector",
                                return_sim = TRUE, progress_bar = FALSE)
        mv   <- metric_fn(sim)
        tv   <- as.numeric(names(mv))
        # Late window, not the whole returned trajectory -- guards against
        # a leading transient in the very first chunk, but no longer a
        # guess at whether the window is long ENOUGH to have settled:
        # projectToSteady() has already made that call (via `tol`) before
        # this window is even taken.
        late <- mv[tv > max(tv) * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late),
             n = sim@n[last, , ], npp = sim@n_pp[last, ], error = NA_character_)
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_,
             n = state_n, npp = state_npp, error = conditionMessage(e))
      })

      if (!is.na(result$error)) {
        warning(sprintf("%s=%.4g: %s", param_name, seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(param_seq)
  bwd    <- run_one_direction(rev(param_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

plot_bifurcation_cod <- function(df, x_label, y_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}

# Centred on cod's own fitted value, +/- 0.3 -- same "bracket the native
# value" convention as Day 30's lambda_seq. Stepped by 0.01 (61 points
# across the bracket) rather than the coarser 21-point/~0.03 step used
# before projectToSteady() replaced project() above: each point now needs
# projectToSteady() to converge from the previous point's own steady
# state, and a small step keeps that starting guess close to the new
# answer, the same reason numerical continuation methods take small steps
# in the swept parameter rather than large ones.
q_seq_bif <- seq(native_cod_q - 0.3, native_cod_q + 0.3, by = 0.01)

# effort=1 -- cod_params' own calibrated fishing gear at its native
# intensity, not the effort=0 used for every stability-only sweep so far.
# Yield is identically zero at effort=0 (nothing is being caught), so
# there's no "unfished yield bifurcation diagram" to fall back on the way
# there was for biomass; this is also the first cod sweep in this project
# to actually fish the population, which Day 29's "What's Next" flagged as
# still outstanding ("introduce fishing mortality against yield").
#
bif_q_df <- run_bifurcation_sweep_cod(
  q_seq_bif, "q",
  fixed_params = list(resource_decrease = 1, capacity_mult = 1),
  effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)

write.csv(bif_q_df, file.path("interesting_plots", "day31_cod_q_yield_bifurcation.csv"),
          row.names = FALSE)

bif_q_plot <- plot_bifurcation_cod(
  bif_q_df, "q (search-volume exponent)", "Yield",
  "Cod bifurcation diagram: yield vs. search-volume exponent q, swept forward and backward",
  sprintf(
    "resource_decrease=1e-7, capacity_mult=1000, effort=1 fixed (Day 30's broad-scan starvation corner, pushed past its own capacity_mult range, at cod's native fishing intensity) -- native q=%.3g; branches collapsing onto one curve = fixed point, fanning apart = limit cycle",
    native_cod_q
  )
)
bif_q_plot

save_plot(bif_q_plot, "day31_cod_q_yield_bifurcation.png", width = 9, height = 6)

# Verdict: per q value, is the max/min spread big enough to call a limit
# cycle (same 1e-6 relative-amplitude threshold used throughout this
# project), and do the forward/backward branches actually disagree
# (hysteresis) rather than retracing each other. Same logic as the
# biomass version, just applied to yield's own max/min instead.
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
    "flat -- 0/%d q values cross the 1e-6 relative-amplitude threshold at (resource_decrease=1e-7, capacity_mult=1000); no forward/backward branch separation either (largest hysteresis gap at q=%.3g, still negligible)",
    nrow(bif_q_check), max_hysteresis_q
  )
} else {
  sprintf(
    "%d/%d q values cross into oscillation at this fixed point; largest forward/backward gap at q=%.3g",
    n_oscillating_q, nrow(bif_q_check), max_hysteresis_q
  )
}
cat(sprintf("\nSection 1 verdict: %s.\n", q_verdict))

################################################################################
# Section 1b: the same bifurcation-diagram sweep, but varying alpha instead
# of q
#
# alpha is mizer's assimilation efficiency (the fraction of consumed energy
# that's actually assimilated, before growth/reproduction allocation) --
# Day 25 found it was the single dominant lever for anchovy's juvenile-
# pileup metric, a genuine viability threshold spanning 27.5 orders of
# magnitude in catchable_fraction, with anchovy's own native alpha sitting
# right at the foot of that climb. That sweep was never carried over to
# cod, and never asked the oscillation question specifically -- this does
# both, reusing Section 1's exact bifurcation-diagram machinery
# (run_bifurcation_sweep_cod()/plot_bifurcation_cod()) with alpha as the
# swept parameter instead of q.
#
# Bracket shape differs from q's deliberately: q's +/- 0.3 is an additive
# bracket appropriate to q's own O(1) scale, but alpha sits on a much
# smaller, strictly-positive scale, so an additive +/-0.3 bracket risks
# pushing alpha negative or barely moving it at all depending on where the
# native value sits. A proportional, log-spaced bracket around cod's own
# native alpha (half to 1.5x) is used instead -- still centred on the
# native value, same "bracket it" logic as every sweep since Day 25, just
# scaled to alpha's own range rather than copying q's absolute numbers.
################################################################################

native_cod_alpha <- unname(cod_params@species_params$alpha[1])
cat(sprintf("cod_params.rds's own native alpha (assimilation efficiency): %.4g\n", native_cod_alpha))

# Same caution as the q propagation check above, adapted to alpha: alpha
# feeds directly into assimilated energy (and therefore growth) at every
# timestep, not into a cached construction-time array the way q feeds
# search_vol -- but this project has been burned before by assuming a
# given_species_params() write actually took effect, so it's confirmed
# directly rather than assumed. getEGrowth() is alpha's own natural
# analogue of getEncounter(): energy available for growth is
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

# alpha applied straight to the real cod_params object, same pattern as
# build_cod_q_variant().
build_cod_alpha_variant <- function(alpha, resource_decrease, capacity_mult) {
  p <- cod_params
  given_species_params(p)$alpha <- alpha
  resource_limitation_2d_cod(p, resource_decrease, capacity_mult)
}

# Proportional bracket, log-spaced, half to 1.5x cod's own native alpha --
# see the header note for why this isn't q's additive +/-0.3.
alpha_seq_bif <- native_cod_alpha * exp(seq(log(0.5), log(1.5), length.out = 21))

# Same fixed point and fishing intensity as Section 1's q sweep, so the two
# bifurcation diagrams are directly comparable against each other, not just
# each internally consistent.
bif_alpha_df <- run_bifurcation_sweep_cod(
  alpha_seq_bif, "alpha",
  fixed_params = list(resource_decrease = 1, capacity_mult = 1),
  params_fn = build_cod_alpha_variant,
  effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)

write.csv(bif_alpha_df, file.path("interesting_plots", "day31_cod_alpha_yield_bifurcation.csv"),
          row.names = FALSE)

bif_alpha_plot <- plot_bifurcation_cod(
  bif_alpha_df, "alpha (assimilation efficiency)", "Yield",
  "Cod bifurcation diagram: yield vs. assimilation efficiency alpha, swept forward and backward",
  sprintf(
    "resource_decrease=1, capacity_mult=1, effort=1 fixed (same fixed point as the q sweep) -- native alpha=%.3g; branches collapsing onto one curve = fixed point, fanning apart = limit cycle",
    native_cod_alpha
  )
)
bif_alpha_plot

save_plot(bif_alpha_plot, "day31_cod_alpha_yield_bifurcation.png", width = 9, height = 6)

# Same verdict logic as Section 1's bif_q_check, applied to alpha instead.
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
cat(sprintf("\nSection 1b verdict: %s.\n", alpha_verdict))

################################################################################
# Section 2: finish the corrected competitiveness metric
#
# Day 29's mean_mass_specific_rate() took a plain mean() over
# mass_specific_rate[idx], where idx selects the juvenile or adult weight
# bins. mizer's weight grid is logarithmically spaced -- dw grows with w --
# so an unweighted mean() treats every bin as equally important regardless
# of how much weight range it actually represents, systematically
# overweighting the smallest individuals in each stage. Day 30 wrote the
# dw-weighted fix (weighted.mean(..., w = dw)) and printed one ratio
# comparison; this section is that same fix, carried through to an
# explicit adult- vs juvenile-driven call for both species -- the thing
# Day 29's original metric was actually being used for.
#
# Self-contained: redefines make_anchovy_params()/anchovy_params here too,
# same convention as every script since Day 20.
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
# below is compared against a freshly-generated number, not a remembered
# one.
mean_mass_specific_rate <- function(params, idx = TRUE) {
  E <- getEncounter(params)
  f <- getFeedingLevel(params)
  per_capita_rate    <- (E * (1 - f))[1, ]
  mass_specific_rate <- per_capita_rate / params@w
  mean(mass_specific_rate[idx])
}

# The fix: weight each bin by dw (the actual width of weight it
# represents) instead of counting every bin equally -- a Riemann-sum-
# consistent average rather than an average-over-indices.
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

  # >1 = juveniles out-eat adults on a mass-specific basis (juvenile-
  # driven, de Roos & Persson's q<1 regime); <1 = adult-driven (q>1).
  # Read directly off the weighted ratio, not the unweighted one -- that's
  # the entire point of the fix.
  stage_driven <- if (ratio_weighted > 1) "juvenile-driven" else "adult-driven"
  flipped      <- (ratio_unweighted > 1) != (ratio_weighted > 1)

  data.frame(
    species                          = label,
    juvenile_adult_ratio_unweighted  = ratio_unweighted,
    q_equivalent_unweighted          = 2 / (ratio_unweighted + 1),
    juvenile_adult_ratio_dw_weighted = ratio_weighted,
    q_equivalent_dw_weighted         = 2 / (ratio_weighted + 1),
    stage_driven_corrected           = stage_driven,
    verdict_flipped                  = flipped
  )
}

corrected_timescale_df <- bind_rows(
  corrected_timescale_summary(cod_params, "cod"),
  corrected_timescale_summary(anchovy_params, "anchovy")
)
print(corrected_timescale_df)
write.csv(corrected_timescale_df, file.path("interesting_plots", "day31_corrected_timescale.csv"),
          row.names = FALSE)

cat("\nSection 2 verdict, per species (dw-weighted, the corrected metric):\n")
for (i in seq_len(nrow(corrected_timescale_df))) {
  row <- corrected_timescale_df[i, ]
  cat(sprintf(
    "  %s: ratio %.4g -> %.4g (unweighted -> dw-weighted), q_equivalent=%.4g -- %s%s\n",
    row$species, row$juvenile_adult_ratio_unweighted, row$juvenile_adult_ratio_dw_weighted,
    row$q_equivalent_dw_weighted, row$stage_driven_corrected,
    if (row$verdict_flipped) " (FLIPPED from the unweighted call)" else " (unchanged from the unweighted call)"
  ))
}

corrected_timescale_plot_df <- corrected_timescale_df %>%
  select(species, juvenile_adult_ratio_unweighted, juvenile_adult_ratio_dw_weighted) %>%
  tidyr::pivot_longer(cols = starts_with("juvenile_adult_ratio"),
                      names_to = "method", values_to = "ratio") %>%
  mutate(method = ifelse(method == "juvenile_adult_ratio_unweighted", "unweighted mean()", "dw-weighted mean"))

corrected_timescale_plot <- ggplot(corrected_timescale_plot_df,
                                   aes(x = species, y = ratio, fill = method)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  labs(x = NULL, y = "juvenile / adult mass-specific intake ratio",
       fill = NULL,
       title = "Correcting Day 29's competitiveness metric for mizer's unequal-width bins",
       subtitle = "Above the dashed line = juvenile-driven, below = adult-driven") +
  theme_minimal()
corrected_timescale_plot

save_plot(corrected_timescale_plot, "day31_corrected_timescale.png", width = 8)

################################################################################
# Section 3: summary -- programmatic readout, not asserted conclusions.
################################################################################

cat("\n===== Day 31 summary =====\n")
cat(sprintf(
  "Section 1 (finished cod q sweep, %d-point forward/backward bifurcation diagram): %s\n",
  length(q_seq_bif), q_verdict
))
cat(sprintf(
  "Section 1b (cod alpha sweep, %d-point forward/backward bifurcation diagram): %s\n",
  length(alpha_seq_bif), alpha_verdict
))
cat(sprintf(
  "Section 2 (finished corrected competitiveness metric): cod %s (ratio %.4g -> %.4g)%s; anchovy %s (ratio %.4g -> %.4g)%s.\n",
  corrected_timescale_df$stage_driven_corrected[corrected_timescale_df$species == "cod"],
  corrected_timescale_df$juvenile_adult_ratio_unweighted[corrected_timescale_df$species == "cod"],
  corrected_timescale_df$juvenile_adult_ratio_dw_weighted[corrected_timescale_df$species == "cod"],
  if (corrected_timescale_df$verdict_flipped[corrected_timescale_df$species == "cod"]) ", FLIPPED" else ", unchanged",
  corrected_timescale_df$stage_driven_corrected[corrected_timescale_df$species == "anchovy"],
  corrected_timescale_df$juvenile_adult_ratio_unweighted[corrected_timescale_df$species == "anchovy"],
  corrected_timescale_df$juvenile_adult_ratio_dw_weighted[corrected_timescale_df$species == "anchovy"],
  if (corrected_timescale_df$verdict_flipped[corrected_timescale_df$species == "anchovy"]) ", FLIPPED" else ", unchanged"
))
