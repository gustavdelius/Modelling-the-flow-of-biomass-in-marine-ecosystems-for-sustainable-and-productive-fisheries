library(mizer)
#library(mizerEcopath)
library(dplyr)
library(ggplot2)
library(scales)
library(mizerExperimental)
p <- readParams("cod_params.rds")

# Avoid automatic recalculation of defaults
p@given_species_params <- p@species_params

# Set power-law resource rate while keeping steady state
rp <- resource_params(p)
rr <- rp$r_pp * w_full(p) ^ rp$n
p <- setResource(p, resource_rate = rr, balance = TRUE)
# We have set `balance = TRUE` to stay at steady state
# This means that the carrying capacity is no longer a
# perfect power law now

# Check that model is still at steady state
sim <- project(p, t_max = 300)
plotBiomass(sim)
# Note the range of the y axis.
# Nevertheless, let's move it closer to steady state
p <- finalParams(sim)

# See how sensitive the model is to fishing
#plotYieldVsF(p, species = "Cod")
# Looks reasonable

plotSpectra(p, power = 2, log = "", size_axis = "l")


# Now you can change the resource rate
rr <- rp$r_pp * w_full(p) ^ 0.5
pp <- setResource(p, resource_rate = rr, balance = FALSE)
pps <- projectToSteady(pp)
attr(pps, "convergence")$type
plotSpectra2(p, pps, power = 2, log = "")
# Unforutnately this does not seem to trigger an instability

# Sweep the exponent instead of testing a single value, to see whether any
# choice -- not just 0.5 -- is enough to trigger the instability that this
# one comparison didn't. Widened to -5..5, and warm-started the same way
# yesterday's n bifurcation sweep was (run_bifurcation_sweep_n_project() in
# 32_experiments.R): only the first point is cold-started from p's own
# steady state, every point after starts from the previous point's final
# state instead of re-cold-starting each time.
n_seq <- seq(-5, 5, by = 0.1)
n_sweep <- data.frame(n = n_seq, convergence_type = NA_character_,
                      error = NA_character_)

state_n   <- NULL
state_npp <- NULL

for (i in seq_along(n_seq)) {
  # tryCatch so one exponent erroring out (e.g. extinction) doesn't kill the
  # whole sweep -- recorded as an NA/error point instead, same convention
  # used throughout this project's other sweeps
  result <- tryCatch({
    rr <- rp$r_pp * w_full(p) ^ n_seq[i]
    pp <- setResource(p, resource_rate = rr, balance = FALSE)
    if (!is.null(state_n)) {
      pp@initial_n[]    <- state_n
      pp@initial_n_pp[] <- state_npp
    }
    pps  <- projectToSteady(pp)
    conv <- attr(pps, "convergence")
    # conv$type has been seen to come back NULL for some exponents even on
    # a successful return, so guard against that rather than let a length-0
    # RHS crash the assignment below
    list(type = if (is.null(conv$type)) NA_character_ else conv$type,
        error = NA_character_,
        n = pps@initial_n, npp = pps@initial_n_pp)
  }, error = function(e) {
    list(type = NA_character_, error = conditionMessage(e),
        n = state_n, npp = state_npp)
  })
  n_sweep$convergence_type[i] <- result$type
  n_sweep$error[i]            <- result$error
  state_n   <- result$n
  state_npp <- result$npp
}

n_sweep
table(n_sweep$convergence_type, useNA = "ifany")

# Extreme fishing constrained to a narrow window around maturity size
#
# A paper predicts that fishing selectivity concentrated narrowly around
# w_mat (rather than mizer's knife_edge default, which selects everything
# above a size with no upper cutoff at all) can destabilise a size-
# structured population into a real limit cycle: narrow-window fishing
# removes exactly the size class that would otherwise regulate
# recruitment, decoupling juveniles from adults. mizer's
# double_sigmoid_length() gives a dome-shaped selectivity (a rising
# sigmoid times a falling sigmoid) that can be centred tightly on l_mat,
# unlike knife_edge's open-ended cutoff.
sp_cod    <- p@species_params
w_mat_cod <- sp_cod$w_mat[sp_cod$species == "Cod"]
a_cod     <- sp_cod$a[sp_cod$species == "Cod"]
b_cod     <- sp_cod$b[sp_cod$species == "Cod"]
l_mat_cod <- (w_mat_cod / a_cod) ^ (1 / b_cod)
cat(sprintf("cod's native w_mat = %.4g g -> l_mat = %.4g cm (via a=%.4g, b=%.4g).\n",
           w_mat_cod, l_mat_cod, a_cod, b_cod))

# Narrow dome: rising sigmoid ramps from 15% below l_mat to 5% below it;
# falling sigmoid ramps from 5% above l_mat back down to 15% above it --
# a tight window centred on maturity length, not knife_edge's "everything
# above a size" cutoff.
gp <- gear_params(p)
gear_name_cod <- gp$gear[gp$species == "Cod"][1]
gp$sel_func[gp$species == "Cod"] <- "double_sigmoid_length"
for (col in c("l25", "l50", "l50_right", "l25_right")) {
  if (!col %in% names(gp)) gp[[col]] <- NA_real_
}
gp$l25[gp$species == "Cod"]       <- l_mat_cod * 0.85
gp$l50[gp$species == "Cod"]       <- l_mat_cod * 0.95
gp$l50_right[gp$species == "Cod"] <- l_mat_cod * 1.05
gp$l25_right[gp$species == "Cod"] <- l_mat_cod * 1.15

p_mat_window <- p
gear_params(p_mat_window) <- gp

# Sanity check: does the new gear actually only select around maturity,
# rather than mizer's default knife_edge (everything above a size)?
sel_at_w <- getSelectivity(p_mat_window)[gear_name_cod, "Cod", ]
cat(sprintf(
  "Selectivity peaks at w=%.4g g (native w_mat=%.4g g); selectivity at the largest size in the model = %.4g (should be ~0, unlike knife_edge's 1).\n",
  p_mat_window@w[which.max(sel_at_w)], w_mat_cod, tail(sel_at_w, 1)
))

# Extreme effort, swept from none up to far beyond anything used
# elsewhere in this project (effort=1 was the fixed baseline throughout
# Day 30-32's sweeps) -- warm-started the same way the n sweep above is,
# only the first point cold-started from p's own steady state. Forward
# only (no backward pass): the question here is whether an oscillation
# appears at all as effort rises, not forward/backward hysteresis.
#
# t_max raised from 2000 to 8000: the first pass through this sweep
# showed convergence time climbing steeply as effort rose (1.8 years at
# effort=5 up to 280 years at effort=16 -- classic critical slowing down
# approaching a bifurcation), so a point just below where the population
# goes extinct could plausibly need much longer than 2000 years to settle
# for real, rather than being cut short mid-transient.
#
# attr(pps, "convergence") came back NULL on every call in the first
# pass -- not just missing $type, the whole attribute -- almost
# certainly because mizerExperimental is loaded after mizer above and
# shadows mizer's own projectToSteady() with an older version that
# doesn't attach it. Calling mizer::projectToSteady() explicitly below
# sidesteps that regardless of load order, and the output method no
# longer depends on the attribute at all: same late-window max/min
# convention as the alpha/q bifurcation sweeps (Day 32), plotted the same
# way -- max and min collapsing onto one line means a fixed point,
# visibly splitting apart means a real oscillation, read directly off
# the plot rather than off a single collapsed rel_amplitude number.
plot_bifurcation_maturity_fishing <- function(df, title, subtitle = NULL) {
  ggplot(df, aes(x = value, y = metric, color = direction, linetype = branch)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    labs(x = "Fishing effort", y = "Yield", title = title, subtitle = subtitle) +
    theme_minimal()
}

# p_variant is passed in explicitly (rather than closed over) so this same
# loop can be reused against any narrow-window selectivity variant --
# maturity size below, largest sizes further down -- without duplicating it.
#
# Forward AND backward now, same convention as the alpha/q bifurcation
# sweeps (Day 32): the backward pass reverses effort_seq and warm-starts
# from the forward pass's own final state, so the two branches collapsing
# onto one curve means a fixed point and visibly splitting apart means
# real hysteresis, not just a one-directional read.
#
# tol tightened from projectToSteady()'s own default (0.1*t_per = 0.02)
# to 1e-5 -- same fix Day 32 applied to the n bifurcation sweep, since a
# loose tolerance can call a point "converged" before it's actually
# settled, which shows up as spurious forward/backward disagreement
# rather than real hysteresis.
run_maturity_fishing_sweep <- function(effort_seq, p_variant, t_max = 8000, tol = 1e-5) {
  run_one_direction <- function(seq_vals, init_n = NULL, init_n_pp = NULL) {
    out       <- data.frame(value = seq_vals, max_metric = NA_real_, min_metric = NA_real_)
    state_n   <- init_n
    state_npp <- init_n_pp

    for (i in seq_along(seq_vals)) {
      result <- tryCatch({
        pm <- p_variant
        if (!is.null(state_n)) {
          pm@initial_n[]    <- state_n
          pm@initial_n_pp[] <- state_npp
        }
        sim  <- mizer::projectToSteady(pm, effort = seq_vals[i], t_max = t_max,
                                       t_per = 0.2, tol = tol, method = "predictor_corrector",
                                       return_sim = TRUE, progress_bar = FALSE)
        mv   <- rowSums(getYield(sim))
        tv   <- as.numeric(names(mv))
        late <- mv[tv > max(tv) * 0.6]
        last <- dim(sim@n)[1]
        list(max_metric = max(late), min_metric = min(late), error = NA_character_,
            n = sim@n[last, , ], npp = sim@n_pp[last, ])
      }, error = function(e) {
        list(max_metric = NA_real_, min_metric = NA_real_, error = conditionMessage(e),
            n = state_n, npp = state_npp)
      })

      if (!is.na(result$error)) {
        warning(sprintf("effort=%.4g: %s", seq_vals[i], result$error))
      }

      state_n   <- result$n
      state_npp <- result$npp
      out$max_metric[i] <- result$max_metric
      out$min_metric[i] <- result$min_metric
    }
    list(df = out, n_final = state_n, npp_final = state_npp)
  }

  fwd    <- run_one_direction(effort_seq)
  bwd    <- run_one_direction(rev(effort_seq), init_n = fwd$n_final, init_n_pp = fwd$npp_final)
  bwd_df <- bwd$df[order(bwd$df$value), ]

  bind_rows(
    data.frame(value = fwd$df$value, metric = fwd$df$max_metric, direction = "Forward",  branch = "max"),
    data.frame(value = fwd$df$value, metric = fwd$df$min_metric, direction = "Forward",  branch = "min"),
    data.frame(value = bwd_df$value, metric = bwd_df$max_metric, direction = "Backward", branch = "max"),
    data.frame(value = bwd_df$value, metric = bwd_df$min_metric, direction = "Backward", branch = "min")
  )
}

# Capped at 14 -- the previous, wider pass showed effort=16 already deep
# in critical slowing down (280 years to converge) and effort>=18 outright
# extinct, so 0-14 keeps this sweep inside the region where the population
# reliably survives, rather than the collapse itself dominating the plot.
effort_seq <- c(seq(0, 5, by = 0.5), seq(6, 14, by = 2))
maturity_fishing_df <- run_maturity_fishing_sweep(effort_seq, p_mat_window)

maturity_fishing_plot <- plot_bifurcation_maturity_fishing(
  maturity_fishing_df,
  "Cod bifurcation diagram: yield vs. fishing effort, dome selectivity centred on maturity size",
  "Narrow double_sigmoid_length window on l_mat (not knife_edge), swept forward and backward -- max/min collapsing onto one curve = fixed point, branches splitting apart = a real oscillation/hysteresis"
)
maturity_fishing_plot

# Biomass trajectories at a handful of the higher-effort points, cold-
# started fresh each time (not warm-started off each other) so each plot
# shows what that specific effort level actually does on its own, rather
# than inheriting a head start from a neighbouring point.
maturity_high_effort <- c(8, 10, 12, 14)
for (e in maturity_high_effort) {
  sim_e <- mizer::project(p_mat_window, effort = e, t_max = 600,
                          method = "predictor_corrector", progress_bar = FALSE)
  print(plotBiomass(sim_e,tlim=c(550,600)) +
         labs(subtitle = sprintf("Maturity-window selectivity, effort=%.4g", e)))
}

# Same narrow-window construction, but anchored on the largest sizes in
# the model (w_inf, cod's asymptotic weight) instead of maturity size --
# but knife_edge instead of double_sigmoid_length here, not the same dome
# shape as the maturity case. At w_mat a knife-edge is wrong because it
# sweeps up everything above it (most of the adult population); at w_inf
# there's essentially nothing left above it to inadvertently include --
# growth is asymptotic, so density beyond w_inf is already thin -- so
# knife_edge(w_inf) isolates the largest sizes at least as surgically as
# a dome, and more so, since it doesn't ramp up starting below w_inf the
# way double_sigmoid_length's l25 (0.85*l_inf) did. knife_edge is weight-
# based directly, so no length-weight (a, b) conversion is needed here.
#
# Rochet & Benoit (2012, Proc. R. Soc. B 279:284-292) found oscillations
# appear at lower fishing intensity and larger amplitude specifically
# when fishing is BOTH selective (a narrow size range) AND targets large
# fish, rather than either alone -- this tests that combination directly,
# reusing the same sweep machinery as the maturity-size case above.
w_inf_cod <- sp_cod$w_inf[sp_cod$species == "Cod"]
cat(sprintf("cod's native w_inf = %.4g g.\n", w_inf_cod))

gp_large <- gear_params(p)
gp_large$sel_func[gp_large$species == "Cod"] <- "knife_edge"
if (!"knife_edge_size" %in% names(gp_large)) gp_large$knife_edge_size <- NA_real_
gp_large$knife_edge_size[gp_large$species == "Cod"] <- w_inf_cod

p_large_window <- p
gear_params(p_large_window) <- gp_large

# Sanity check: selectivity should be a hard 0/1 step at w_inf now, not a
# dome -- 0 everywhere below it, 1 everywhere at or above it (including
# at the largest size in the model, unlike the dome's ~0 there).
sel_at_w_large <- getSelectivity(p_large_window)[gear_name_cod, "Cod", ]
w_just_below_inf <- which.min(abs(p_large_window@w - 0.9 * w_inf_cod))
cat(sprintf(
  "Selectivity at 90%% of w_inf (w=%.4g g) = %.4g (should be 0); selectivity at the largest size in the model = %.4g (should be 1).\n",
  p_large_window@w[w_just_below_inf], sel_at_w_large[w_just_below_inf],
  tail(sel_at_w_large, 1)
))

# Extended to 30 (not capped at 14 like the maturity-window sweep above):
# large fish are a much smaller share of the spectrum's abundance than
# fish at maturity size, so the same effort level bites less hard here --
# worth checking the fuller range rather than assuming it collapses at
# the same point the maturity window did.
effort_seq_large <- c(seq(0, 5, by = 0.5), seq(6, 30, by = 2))
large_fishing_df <- run_maturity_fishing_sweep(effort_seq_large, p_large_window)

large_fishing_plot <- plot_bifurcation_maturity_fishing(
  large_fishing_df,
  "Cod bifurcation diagram: yield vs. fishing effort, knife_edge selectivity at w_inf",
  "knife_edge(w_inf) -- fully selects the largest sizes only, unlike w_mat where a knife-edge would sweep up the whole adult population -- swept forward and backward -- max/min collapsing onto one curve = fixed point, branches splitting apart = a real oscillation/hysteresis"
)
large_fishing_plot

# Same handful-of-high-effort biomass trajectories as the maturity-window
# case, cold-started fresh each time, spread across the wider 0-30 range.
large_high_effort <- c(14, 20, 25, 30)
for (e in large_high_effort) {
  sim_e <- mizer::project(p_large_window, effort = e, t_max = 600,
                          method = "predictor_corrector", progress_bar = FALSE)
  print(plotBiomass(sim_e,tlim=c(500,600)) +
         labs(subtitle = sprintf("knife_edge(w_inf) selectivity, effort=%.4g", e)))
}



