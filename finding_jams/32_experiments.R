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
# Day 32: two fixes to Day 31's q and alpha bifurcation sweeps.
#
# 1. Both sweeps went through the generic run_bifurcation_sweep_cod(
#    param_seq, param_name, params_fn = ...) dispatcher -- one step removed
#    from how the sweep loop actually reads on the page, since that
#    dispatcher was generalised FROM a q-specific loop in the first place
#    (Day 23's Follow-up 5). Both are rewritten here as direct,
#    un-indirected projectToSteady() loops hardcoded to the parameter being
#    swept -- no param_name string, no params_fn argument, no fixed_params
#    list. The projectToSteady() call itself (method, t_per, tryCatch
#    convention) is unchanged from Day 31, since that part was already
#    right.
#
# 2. Neither sweep had actually been run at the corner Day 30's own broad
#    scan found genuinely oscillating (resource_decrease~=1e-7,
#    capacity_mult=1, rel_amplitude=0.041) -- both sat instead at
#    resource_decrease=1, capacity_mult=1, the corner every sweep so far
#    has found flat. Both are re-pointed at the oscillating corner here.
#
# 3. A third sweep (Section 3) tries a different lever entirely:
#    resource_rate on its own, with balance=TRUE so resource_capacity is
#    derived implicitly from it rather than set independently the way
#    resource_limitation_2d_cod() sets both. Sections 1 and 2 both held
#    the resource state fixed and varied a species trait instead; this
#    asks whether letting resource_rate and resource_capacity move
#    together, rather than independently, is itself what surfaces an
#    oscillation the other two didn't find.
#
# 4. Section 3's own sweep never actually perturbs the system -- it's
#    either a cold start or a warm start carried from the adjacent point,
#    and projectToSteady() exits the moment its own convergence tolerance
#    is satisfied. Section 4 tests whether that tolerance is too generous
#    -- too quick to call a still-oscillating point "converged" -- by
#    rerunning the forward sweep two ways, projectToSteady() vs a
#    no-early-exit project() (t_first=800 cold start, t_rest=100 per
#    warm-started step), and comparing the full yield time series rather
#    than each point's collapsed late-window max/min.
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
# Re-pointed at the corner that actually oscillates, not the flat one --
# resource_decrease=1, capacity_mult=1 (used by Day 31's q sweep, and by
# this file's own first alpha run) was always the boring corner; Day 30's
# own broad scan found a confirmed real oscillation (rel_amplitude=0.041)
# at resource_decrease~=1e-7, capacity_mult=1 -- extreme resource
# starvation at cod's native carrying capacity -- and neither the q nor
# the alpha sweep had actually been run there until now.
bif_alpha_df <- run_bifurcation_sweep_alpha(
  alpha_seq_bif,
  resource_decrease = 1e-7, capacity_mult = 1,
  effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)

write.csv(bif_alpha_df, file.path("interesting_plots", "day32_cod_alpha_yield_bifurcation.csv"),
          row.names = FALSE)

bif_alpha_plot <- plot_bifurcation_alpha(
  bif_alpha_df, "alpha (assimilation efficiency)", "Yield",
  "Cod bifurcation diagram: yield vs. assimilation efficiency alpha, swept forward and backward",
  sprintf(
    # Stated directly, not "same fixed point as Day 31's q sweep" -- it
    # deliberately isn't anymore. Day 31's q sweep (and this file's own
    # first run) sat at resource_decrease=1, capacity_mult=1, the corner
    # every sweep so far has found flat; this run is at Day 30's own
    # confirmed-oscillating corner instead.
    "resource_decrease=1e-7, capacity_mult=1, effort=1 fixed (Day 30's confirmed-oscillating corner, rel_amplitude=0.041) -- native alpha=%.3g; both axes log-scaled; branches collapsing onto one curve = fixed point, fanning apart = limit cycle",
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
    "flat -- 0/%d alpha values cross the 1e-6 relative-amplitude threshold at (resource_decrease=1e-7, capacity_mult=1, Day 30's confirmed-oscillating corner); no forward/backward branch separation either (largest hysteresis gap at alpha=%.4g, still negligible)",
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
# Section 2: q bifurcation, same direct-loop rewrite as Section 1's alpha
# sweep, re-pointed at the same confirmed-oscillating corner
#
# Day 31's q sweep never actually tested the corner Day 30's own broad scan
# found genuinely oscillating -- it ran at resource_decrease=1,
# capacity_mult=1 (the corner every sweep so far, q and alpha both, has
# found flat) instead of resource_decrease~=1e-7, capacity_mult=1
# (rel_amplitude=0.041). This section is q's own half of that fix,
# written the same direct, un-indirected way Section 1 rewrote alpha's --
# no run_bifurcation_sweep_cod() dispatcher, no param_name/params_fn
# indirection, build_cod_q_variant() called explicitly inside the loop.
################################################################################

native_cod_q <- unname(cod_params@species_params$q[1])
cat(sprintf("cod_params.rds's own native q (search-volume exponent): %.4g\n", native_cod_q))

# Same propagation caution as Section 1's alpha check, adapted to q: q
# feeds into SearchVolume = gamma * w^q, which mizer may cache at
# construction time rather than recompute live, so this is confirmed
# directly (via search_vol/getEncounter()) rather than assumed -- the same
# check Day 31 ran before trusting its own q sweep.
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

# q applied straight to the real cod_params object -- same pattern as
# Section 1's build_cod_alpha_variant().
build_cod_q_variant <- function(q, resource_decrease, capacity_mult) {
  p <- cod_params
  given_species_params(p)$q <- q
  resource_limitation_2d_cod(p, resource_decrease, capacity_mult)
}

# Direct port of Section 1's run_bifurcation_sweep_alpha(), hardcoded to q
# instead -- same method/t_per/tryCatch conventions, unchanged from Day 31.
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

# Linear axes, not log -- unlike alpha's bracket, q's own +/- 0.3 additive
# bracket around its native O(1) value never spans more than about a unit,
# so there's no order-of-magnitude range for a log scale to usefully
# compress the way there was for alpha's 0.005-1 sweep.
plot_bifurcation_q <- function(df, x_label, y_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}

# Same bracket as Day 31: centred on cod's own fitted value, +/- 0.3,
# stepped by 0.01 -- unchanged here since today's fix is about the fixed
# point (resource_decrease/capacity_mult), not the swept range itself.
q_seq_bif <- seq(native_cod_q - 0.3, native_cod_q + 0.3, by = 0.01)

# Re-pointed at the corner that actually oscillates -- see the header note
# above and Section 1's identical fix for alpha.
bif_q_df <- run_bifurcation_sweep_q(
  q_seq_bif,
  resource_decrease = 1e-7, capacity_mult = 1,
  effort = 1,
  metric_fn = function(sim) rowSums(getYield(sim))
)

write.csv(bif_q_df, file.path("interesting_plots", "day32_cod_q_yield_bifurcation.csv"),
          row.names = FALSE)

bif_q_plot <- plot_bifurcation_q(
  bif_q_df, "q (search-volume exponent)", "Yield",
  "Cod bifurcation diagram: yield vs. search-volume exponent q, swept forward and backward",
  sprintf(
    "resource_decrease=1e-7, capacity_mult=1, effort=1 fixed (Day 30's confirmed-oscillating corner, rel_amplitude=0.041) -- native q=%.3g; branches collapsing onto one curve = fixed point, fanning apart = limit cycle",
    native_cod_q
  )
)
bif_q_plot

save_plot(bif_q_plot, "day32_cod_q_yield_bifurcation.png", width = 9, height = 6)

# Same verdict logic as Section 1's bif_alpha_check, applied to q instead.
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
    "flat -- 0/%d q values cross the 1e-6 relative-amplitude threshold at (resource_decrease=1e-7, capacity_mult=1, Day 30's confirmed-oscillating corner); no forward/backward branch separation either (largest hysteresis gap at q=%.3g, still negligible)",
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
# Section 3: resource_rate bifurcation with balance=TRUE
#
# Sections 1 and 2 both varied a species-level trait (alpha, q) at a fixed
# resource state, reached by scaling resource_rate and resource_capacity
# independently with balance=FALSE (resource_limitation_2d_cod() -- the two
# move separately, capacity_mult never touching rate and vice versa). This
# section asks a different question: what if resource_capacity is never
# set directly at all, but left for mizer to derive from resource_rate via
# setResource(..., balance = TRUE) -- does letting the two move together,
# rather than independently, surface an oscillation neither Section 1 nor
# Section 2 (nor Day 31's identically-shaped sweeps) found?
#
# Same forward/backward bifurcation-diagram machinery as Sections 1 and 2,
# same projectToSteady() settings (method, t_per, tryCatch convention) --
# only the parameter build function changes.
################################################################################

# Only resource_rate is reassigned -- no resource_capacity argument passed
# to setResource() at all, balance = TRUE is what lets mizer derive it
# implicitly rather than this function choosing it.
resource_limitation_balance_cod <- function(params, resource_rate_mult) {
  new_rate <- getResourceRate(params) * resource_rate_mult
  setResource(params, resource_rate = new_rate, balance = TRUE)
}

# Same propagation caution as Sections 1 and 2's checks, adapted here: the
# entire premise of this section is that resource_capacity moves as a
# side-effect of resource_rate under balance = TRUE, so both getResourceRate()
# and getResourceCapacity() have to change, or there's nothing to sweep.
resource_rate_probe_base <- cod_params
resource_rate_probe_low  <- resource_limitation_balance_cod(cod_params, 0.5)

resource_rate_changed     <- !isTRUE(all.equal(getResourceRate(resource_rate_probe_base),
                                               getResourceRate(resource_rate_probe_low)))
resource_capacity_changed <- !isTRUE(all.equal(getResourceCapacity(resource_rate_probe_base),
                                               getResourceCapacity(resource_rate_probe_low)))

cat(sprintf("resource_rate (balance=TRUE) propagation check: resource_rate changed = %s, resource_capacity changed = %s.\n",
           resource_rate_changed, resource_capacity_changed))
if (!resource_rate_changed || !resource_capacity_changed) {
  stop(paste(
    "setResource(params, resource_rate = ..., balance = TRUE) is NOT moving",
    "both resource_rate and resource_capacity -- either balance = TRUE isn't",
    "deriving capacity the way this section assumes, or resource_rate itself",
    "isn't propagating. Stopping before wasting the run."
  ))
}

# resource_rate_mult applied straight to the real cod_params object -- same
# pattern as Sections 1 and 2's build_cod_alpha_variant()/build_cod_q_variant(),
# just a single argument here since there's no separate capacity_mult to pass.
build_cod_resource_rate_variant <- function(resource_rate_mult) {
  resource_limitation_balance_cod(cod_params, resource_rate_mult)
}

# Direct port of Sections 1 and 2's sweep loops, hardcoded to resource_rate_mult.
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

# Log-scaled, both axes -- the swept range spans seven orders of magnitude
# in resource_rate_mult (1e-7 to 1), the same reason Section 1's widened
# alpha bracket got log axes rather than linear ones.
plot_bifurcation_resource_rate <- function(df, x_label, y_label, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = x_label, y = y_label, title = title, subtitle = subtitle) +
    theme_minimal()
}
# 1e-7 to 1, log-spaced, ~100 points -- 1e-7 is the same extreme-starvation
# multiplier Sections 1 and 2 tested, 1 is cod's own native resource_rate
# (no change at all), so this brackets the entire range from native down
# to the most extreme starvation tested so far, rather than bracketing a
# native value the way alpha/q's own brackets did.
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
  "effort=1 fixed; resource_capacity derived implicitly from resource_rate via setResource(balance=TRUE), not set independently; both axes log-scaled; branches collapsing onto one curve = fixed point, fanning apart = limit cycle"
)
bif_resource_rate_plot

save_plot(bif_resource_rate_plot, "day32_cod_resource_rate_bif.png", width = 9, height = 6)

# Same verdict logic as Sections 1 and 2's bif_alpha_check/bif_q_check,
# applied to resource_rate_mult instead. Same caveat those sections' own
# results surfaced today applies here too: the 1e-6 relative-amplitude
# threshold is only trustworthy at yield scales comparable to where it was
# first calibrated -- read the printed verdict alongside the plot itself,
# not instead of it.
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
# Section 4: does projectToSteady() exit too early? project() vs
# projectToSteady(), as a continuous time series across the whole sweep
#
# Section 3's bifurcation diagram never actually perturbs the system --
# every point after the first is a warm start carried over from its
# adjacent neighbour, and projectToSteady() exits the moment ITS OWN
# convergence tolerance is satisfied. If that tolerance is too generous --
# declaring "converged" while the population is still slowly oscillating
# -- every point in Section 3 would read as a fixed point even where a
# genuine, slower limit cycle exists underneath, regardless of which
# corner the sweep is re-pointed at.
#
# This reruns the forward direction of Section 3's own sweep two ways,
# tracking the full yield trajectory -- not just each point's late-window
# max/min -- as one continuous time series across the whole run:
#
#   1. projectToSteady() -- exactly as Section 3 uses it (method =
#      "predictor_corrector", t_per = 0.2, t_max = 2000 as a safety
#      ceiling) -- exits each point the moment its own convergence check
#      is satisfied.
#   2. project() directly, with no early exit at all -- t_first = 800
#      years for the first, cold-started point (the same long cold-start
#      budget earlier sweeps in this project gave a fresh start), t_rest
#      = 100 years for every subsequent, warm-started point.
#
# Only the forward direction is run here, not forward+backward -- this is
# a targeted test of the early-exit hypothesis, not a second full
# bifurcation diagram, and forward alone already starts at
# resource_rate_mult=1e-7, the single most extreme, most diagnostic point
# in Section 3's own range.
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
# Section 5: summary -- programmatic readout, not asserted conclusions.
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
  "Section 3 (resource_rate bifurcation, balance=TRUE deriving capacity implicitly, %d-point forward/backward bifurcation diagram): %s\n",
  length(resource_rate_seq_bif), resource_rate_verdict
))
cat(
  "Section 4 (project() vs projectToSteady() time-series comparison across the forward resource_rate sweep): see day32_resource_rate_timeseries.png -- read visually, no automated verdict.\n"
)
