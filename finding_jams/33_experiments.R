library(mizer)
library(mizerEcopath)
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
plotYieldVsF(p, species = "Cod")
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
# only the first point cold-started from p's own steady state.
effort_seq <- c(seq(0, 5, by = 0.5), seq(6, 30, by = 2))
maturity_fishing_sweep <- data.frame(effort = effort_seq,
                                     convergence_type = NA_character_,
                                     error = NA_character_)

state_n   <- NULL
state_npp <- NULL

for (i in seq_along(effort_seq)) {
  result <- tryCatch({
    pm <- p_mat_window
    if (!is.null(state_n)) {
      pm@initial_n[]    <- state_n
      pm@initial_n_pp[] <- state_npp
    }
    pps  <- projectToSteady(pm, effort = effort_seq[i], t_max = 2000,
                            t_per = 0.2, method = "predictor_corrector",
                            progress_bar = FALSE)
    conv <- attr(pps, "convergence")
    list(type = if (is.null(conv$type)) NA_character_ else conv$type,
        error = NA_character_,
        n = pps@initial_n, npp = pps@initial_n_pp)
  }, error = function(e) {
    list(type = NA_character_, error = conditionMessage(e),
        n = state_n, npp = state_npp)
  })
  maturity_fishing_sweep$convergence_type[i] <- result$type
  maturity_fishing_sweep$error[i]            <- result$error
  state_n   <- result$n
  state_npp <- result$npp
}

maturity_fishing_sweep
table(maturity_fishing_sweep$convergence_type, useNA = "ifany")

# mizer's own projectToSteady() cycle detector (autocorrelation-based,
# requires 3 persistent periods -- see Day 32) firing type=="cycle" is a
# much stronger signal than the ad hoc 1e-6 relative-amplitude threshold
# Day 32 found unreliable at small yields.
cycle_efforts <- maturity_fishing_sweep$effort[maturity_fishing_sweep$convergence_type == "cycle"]
if (length(cycle_efforts) > 0) {
  cat(sprintf("Oscillation reproduced: convergence type == 'cycle' at effort = %s.\n",
             paste(cycle_efforts, collapse = ", ")))
} else {
  cat("No effort level in this sweep converged to a genuine limit cycle (type=='cycle').\n")
}

# Visualise the most extreme point directly, regardless of the verdict above
pps_extreme <- projectToSteady(p_mat_window, effort = max(effort_seq), t_max = 2000,
                               t_per = 0.2, method = "predictor_corrector",
                               return_sim = TRUE, progress_bar = FALSE)
plotBiomass(pps_extreme)
attr(pps_extreme, "convergence")



