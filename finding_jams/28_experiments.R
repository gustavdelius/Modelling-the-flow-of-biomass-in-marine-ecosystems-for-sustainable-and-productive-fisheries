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
library(future)
library(future.apply)
library(scales)
library(deSolve)

# Each check_hopf_grid_point() call is independent (its own steadyNewton()
# solve), so the grid sweeps below run across workers rather than
# sequentially. multisession works cross-platform (unlike multicore, which
# is a no-op on Windows).
plan(multisession)

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH (260 chars) has silently truncated filenames -- even the
# .png extension itself -- once combined with this repo's long, deeply
# nested folder path. Keep names short at the call site; this is a
# defensive last resort so a long one truncates safely instead of silently.
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

# Read-only: never saveRDS() back to this path. resource_limitation() always
# returns a new object, so this binding itself is never reassigned.
cod_params <- readRDS("cod_params.rds")
plotSpectra(cod_params)

# balance=FALSE deliberately: mizer's default would partly auto-compensate
# for the perturbation this sweep is meant to introduce.
resource_limitation <- function(params, resource_decrease) {
  new_rate <- getResourceRate(params) * resource_decrease
  setResource(params, resource_rate = new_rate, balance = FALSE)
}

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

limited_cod_params <- resource_limitation(cod_params,0.05)
sim <- make_limit_cycle_sim(limited_cod_params)
animateSpectra(sim,tlim=c(550,600))


################################################################################
# Bifurcation diagram: does slowing the resource's renewal destabilise cod?
# Every local max/min of biomass in the late window becomes its own point --
# one dot = settled to a fixed point; two branches = a genuine limit cycle.
################################################################################

# Single-species, so getBiomass() is a one-column matrix -- rowSums() just
# extracts it, kept general in case cod_params ever becomes multi-species.
cod_biomass_series <- function(sim) {
  bm <- getBiomass(sim)
  data.frame(time = as.numeric(rownames(bm)), biomass = rowSums(bm))
}

find_local_extrema <- function(times, y) {
  is_max <- c(FALSE, y[-1] > y[-length(y)]) & c(y[-length(y)] > y[-1], FALSE)
  is_min <- c(FALSE, y[-1] < y[-length(y)]) & c(y[-length(y)] < y[-1], FALSE)
  data.frame(
    time  = c(times[is_max], times[is_min]),
    value = c(y[is_max], y[is_min]),
    type  = c(rep("max", sum(is_max)), rep("min", sum(is_min)))
  )
}

resource_bifurcation_point <- function(resource_decrease, params = cod_params, t_total = 600,
                                       window_frac = 0.4) {
  p    <- resource_limitation(params, resource_decrease)
  s    <- make_limit_cycle_sim(p, t_total = t_total)
  bm   <- cod_biomass_series(s)
  late <- bm[bm$time >= t_total * (1 - window_frac), ]

  extrema <- find_local_extrema(late$time, late$biomass)
  if (nrow(extrema) == 0) {
    # Settled to a fixed point -- report the final value as one dot.
    extrema <- data.frame(time = tail(late$time, 1), value = tail(late$biomass, 1),
                          type = "settled")
  }
  extrema$resource_decrease <- resource_decrease
  extrema
}

# Log-spaced 0.001 (usual default) to 1 (unchanged); kept to 12 points since
# each one is a full mizer simulation, not a cheap closed-form scan.
resource_decrease_seq <- exp(seq(log(0.001), log(1), length.out = 12))

resource_bifurcation_df <- bind_rows(lapply(resource_decrease_seq, resource_bifurcation_point))

print(resource_bifurcation_df)
write.csv(resource_bifurcation_df, file.path("interesting_plots", "cod_resource_bifurcation.csv"),
          row.names = FALSE)

resource_bifurcation_plot <- ggplot(resource_bifurcation_df, aes(x = resource_decrease, y = value)) +
  geom_point(aes(color = type), alpha = 0.6, size = 1.5) +
  scale_x_log10() +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "Cod biomass -- local extrema in the late window",
       color = NULL,
       title = "Bifurcation diagram: does slowing the resource's renewal destabilise cod?",
       subtitle = "One settled dot = fixed point; a max/min branch split = a genuine limit cycle") +
  theme_minimal()
resource_bifurcation_plot

save_plot(resource_bifurcation_plot, "cod_bif_diagram.png")

################################################################################
# Testing the new (experimental) getStability()/steadyNewton(stability=)
# from mizer PR #452 -- eigenvalues of dG/dN at the steady state, a
# mizer-native counterpart to the closed-form Hopf conditions from Days
# 25-27. Checked against NS_params (PR's own test) and against the
# simulated bifurcation diagram above, before trusting it on cod_params.
################################################################################

# 1. Sanity check against NS_params, matching the PR's own test.
ns_steady    <- steadyNewton(NS_params, stability = TRUE)
ns_stability <- attr(ns_steady, "stability")
cat(sprintf(
  "NS_params sanity check: stable = %s, spectral_radius = %.4f (expect TRUE / < 1).\n",
  ns_stability$stable, ns_stability$spectral_radius
))

# 2. Reduced vs full agreement on cod_params (PR's tests use 5% tolerance).
cod_steady   <- steadyNewton(cod_params, stability = TRUE)
stab_reduced <- attr(cod_steady, "stability")
stab_full    <- getStability(cod_steady, include_resource = TRUE)

cat(sprintf(
  paste0(
    "cod_params reduced vs full: spectral_radius reduced = %.4f, full = %.4f ",
    "(relative difference %.2f%%, n_active %d vs %d).\n"
  ),
  stab_reduced$spectral_radius, stab_full$spectral_radius,
  100 * abs(stab_reduced$spectral_radius - stab_full$spectral_radius) / stab_reduced$spectral_radius,
  stab_reduced$n_active, stab_full$n_active
))

# 3. Predicted stability across resource_decrease_seq, wrapped since
# steadyNewton() can fail to converge at extreme resource_decrease.
stability_sweep_point <- function(resource_decrease, params = cod_params) {
  result <- tryCatch({
    p        <- resource_limitation(params, resource_decrease)
    p_steady <- steadyNewton(p, stability = TRUE)
    stab     <- attr(p_steady, "stability")
    data.frame(spectral_radius = stab$spectral_radius,
               stable = stab$stable,
               hopf_period = if (is.null(stab$hopf_period)) NA_real_ else stab$hopf_period,
               dominant_period = stab$dominant_period,
               n_active = stab$n_active,
               error = NA_character_)
  }, error = function(e) {
    data.frame(spectral_radius = NA_real_, stable = NA, hopf_period = NA_real_,
               dominant_period = NA_real_, n_active = NA_integer_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("resource_decrease = %.4g: %s", resource_decrease, result$error))
  }
  data.frame(resource_decrease = resource_decrease, result)
}

stability_sweep_df <- bind_rows(lapply(resource_decrease_seq, stability_sweep_point))

print(stability_sweep_df)
write.csv(stability_sweep_df, file.path("interesting_plots", "cod_stability_sweep.csv"),
          row.names = FALSE)

spectral_radius_plot <- ggplot(stability_sweep_df %>% filter(is.na(error)),
                               aes(x = resource_decrease, y = spectral_radius)) +
  geom_line(color = "#8172B2", linewidth = 1) +
  geom_point(color = "#8172B2", size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "Spectral radius of dG/dN at the steady state",
       title = "getStability()'s predicted stability boundary",
       subtitle = "Dashed = spectral_radius = 1 -- stable below, unstable above") +
  theme_minimal()
spectral_radius_plot

# Stacked on the simulated bifurcation diagram, same x-axis, to compare the
# predicted crossing directly against the simulated one.
bifurcation_vs_stability_plot <- resource_bifurcation_plot / spectral_radius_plot +
  plot_layout(guides = "collect")
bifurcation_vs_stability_plot

save_plot(bifurcation_vs_stability_plot, "cod_bif_vs_stab.png",
          height = 10)

# 4. Predicted vs measured period -- reuses the "max" extrema already in
# resource_bifurcation_df, no rerun needed.
measured_period_df <- resource_bifurcation_df %>%
  filter(type == "max") %>%
  group_by(resource_decrease) %>%
  summarise(measured_period = if (n() >= 2) mean(diff(sort(time))) else NA_real_, .groups = "drop")

period_comparison_df <- stability_sweep_df %>%
  filter(is.na(error)) %>%
  left_join(measured_period_df, by = "resource_decrease")

print(period_comparison_df)
write.csv(period_comparison_df, file.path("interesting_plots", "cod_period_comparison.csv"),
          row.names = FALSE)

period_comparison_plot <- ggplot(period_comparison_df, aes(x = resource_decrease)) +
  geom_line(aes(y = hopf_period, color = "predicted (getStability hopf_period)"), linewidth = 1) +
  geom_point(aes(y = measured_period, color = "measured (peak spacing)"), size = 2.5) +
  scale_x_log10() +
  labs(x = "resource_decrease (multiplier on resource renewal rate)", y = "Oscillation period (time steps)",
       title = "getStability()'s predicted Hopf period vs the simulated oscillation period",
       subtitle = "Only meaningful where the bifurcation diagram actually shows a limit cycle", color = NULL) +
  theme_minimal()
period_comparison_plot

save_plot(period_comparison_plot, "cod_period_pred_meas.png")

################################################################################
# Does a Hopf bifurcation appear when varying BOTH resource carrying
# capacity and resource rate together? "carrying capacity" here is
# resource_capacity (mizer's capacity_mult from Days 20-24 -- never
# confirmed there as specifically Hopf, just "bifurcation"/"paradox of
# enrichment"). spectral_radius alone can't distinguish a Hopf bifurcation
# from a real-eigenvalue instability (e.g. plain collapse/extinction) --
# only a complex dominant eigenvalue crossing the unit circle is Hopf, so
# that's checked directly via the sign of its imaginary part.
################################################################################

# balance=FALSE isn't just stylistic here: setResource() errors outright if
# both resource_rate and resource_capacity are given while balance is TRUE
# ("You should only provide either..."), since with balance=TRUE the other
# one is always computed FOR you. Setting both independently requires
# balance=FALSE.
resource_limitation_2d <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}

# effort defaults to 0 (unfished), matching the original grid below exactly;
# passed straight through to steadyNewton() the same way project()/steady()
# take it elsewhere in mizer, so the steady state (and its stability) is
# solved for under that constant fishing mortality rather than unfished.
check_hopf_grid_point <- function(resource_decrease, capacity_mult, params = cod_params, effort = 0) {
  result <- tryCatch({
    p        <- resource_limitation_2d(params, resource_decrease, capacity_mult)
    p_steady <- steadyNewton(p, effort = effort, stability = TRUE)
    stab     <- attr(p_steady, "stability")
    dominant <- stab$eigenvalues[1]
    is_complex <- abs(Im(dominant)) > 1e-8

    data.frame(
      spectral_radius  = stab$spectral_radius,
      stable           = stab$stable,
      hopf_period      = if (is.null(stab$hopf_period)) NA_real_ else stab$hopf_period,
      bifurcation_type = if (stab$stable) "stable" else if (is_complex) "hopf" else "non-oscillatory",
      error            = NA_character_
    )
  }, error = function(e) {
    data.frame(spectral_radius = NA_real_, stable = NA, hopf_period = NA_real_,
               bifurcation_type = NA_character_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("resource_decrease=%.4g, capacity_mult=%.4g, effort=%.4g: %s",
                    resource_decrease, capacity_mult, effort, result$error))
  }
  data.frame(resource_decrease = resource_decrease, capacity_mult = capacity_mult,
             effort = effort, result)
}

# capacity_mult spans the same direction as Days 22-24's oscillation find
# (destabilisation came from RAISING capacity, not lowering it).
capacity_mult_seq <- exp(seq(log(1), log(20), length.out = 10))

# stability = TRUE makes each point far from cheap: getStability() builds
# its Jacobian by numerical differentiation, one project_n_loop() call per
# active weight-class cell (x2, centred differences) ON TOP OF whatever
# steadyNewton() itself needs to converge. With n_active potentially in the
# hundreds, that's hundreds of full solver evaluations per grid point --
# flattened to one grid of 120 independent points and run across workers
# (plan(multisession), set above) rather than sequentially.
hopf_grid_params <- expand.grid(resource_decrease = resource_decrease_seq,
                                capacity_mult = capacity_mult_seq)

hopf_grid_df <- bind_rows(future_lapply(seq_len(nrow(hopf_grid_params)), function(i) {
  check_hopf_grid_point(hopf_grid_params$resource_decrease[i], hopf_grid_params$capacity_mult[i])
}, future.seed = TRUE))

print(hopf_grid_df)
write.csv(hopf_grid_df, file.path("interesting_plots", "cod_hopf_grid.csv"), row.names = FALSE)

cat(sprintf(
  "Any grid point with a genuine Hopf-type instability (complex dominant eigenvalue)? %s\n",
  any(hopf_grid_df$bifurcation_type == "hopf", na.rm = TRUE)
))

hopf_grid_plot <- ggplot(hopf_grid_df %>% filter(is.na(error)),
                         aes(x = resource_decrease, y = capacity_mult, fill = bifurcation_type)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c(stable = "#4C72B0", hopf = "#C44E52", `non-oscillatory` = "#DD8452"),
                    na.value = "grey80") +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       fill = NULL,
       title = "Where does varying resource rate and capacity together cause a Hopf bifurcation?",
       subtitle = "hopf = spectral_radius >= 1 with a complex dominant eigenvalue") +
  theme_minimal()
hopf_grid_plot

save_plot(hopf_grid_plot, "cod_hopf_grid.png")



################################################################################
# Testing out fishing pressure's effect on whether it can cause a Hopf
# bifurcation: same (resource_decrease, capacity_mult) grid, same
# classification, just re-solved at a constant fishing effort of 0.5
# instead of the unfished baseline (effort = 0) used above.
################################################################################

FISHING_EFFORT_TEST <- 0.5

hopf_grid_fished_df <- bind_rows(future_lapply(seq_len(nrow(hopf_grid_params)), function(i) {
  check_hopf_grid_point(hopf_grid_params$resource_decrease[i], hopf_grid_params$capacity_mult[i],
                        effort = FISHING_EFFORT_TEST)
}, future.seed = TRUE))

print(hopf_grid_fished_df)
write.csv(hopf_grid_fished_df, file.path("interesting_plots", "cod_hopf_grid_fished.csv"),
          row.names = FALSE)

cat(sprintf(
  "Unfished grid: %d/%d points classified 'hopf'. Fished (effort=%.1f): %d/%d.\n",
  sum(hopf_grid_df$bifurcation_type == "hopf", na.rm = TRUE), nrow(hopf_grid_df),
  FISHING_EFFORT_TEST,
  sum(hopf_grid_fished_df$bifurcation_type == "hopf", na.rm = TRUE), nrow(hopf_grid_fished_df)
))

hopf_grid_compare_df <- bind_rows(
  hopf_grid_df        %>% mutate(effort_label = "effort = 0 (unfished)"),
  hopf_grid_fished_df %>% mutate(effort_label = sprintf("effort = %.1f", FISHING_EFFORT_TEST))
)

hopf_grid_fished_plot <- ggplot(hopf_grid_compare_df %>% filter(is.na(error)),
                                aes(x = resource_decrease, y = capacity_mult, fill = bifurcation_type)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c(stable = "#4C72B0", hopf = "#C44E52", `non-oscillatory` = "#DD8452"),
                    na.value = "grey80") +
  facet_wrap(~effort_label) +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       fill = NULL,
       title = "Does fishing pressure help trigger a Hopf bifurcation?",
       subtitle = "hopf = spectral_radius >= 1 with a complex dominant eigenvalue") +
  theme_minimal()
hopf_grid_fished_plot

save_plot(hopf_grid_fished_plot, "cod_hopf_grid_fished_cmp.png", width = 12)



E <- getEncounter(cod_params)      # encounter rate, g/year -- E_i(w)
f <- getFeedingLevel(cod_params)   # feeding level (satiation, 0-1) -- f_i(w) = E/(E+h)

per_capita_rate <- E * (1 - f)     # realized per-capita consumption, g/year

w     <- cod_params@w
w_mat <- cod_params@species_params$w_mat

juvenile_rate <- per_capita_rate[1, w < w_mat]
adult_rate    <- per_capita_rate[1, w >= w_mat]

juvenile_rate - adult_rate
mean(juvenile_rate) - mean(adult_rate)


