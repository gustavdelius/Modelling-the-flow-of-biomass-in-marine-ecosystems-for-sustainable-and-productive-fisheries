suppressPackageStartupMessages({library(mizer); library(mizerExperimental); library(ggplot2)})

dir.create("interesting_plots", showWarnings = FALSE)

cod_params <- readRDS("cod_params.rds")

resource_limitation_2d_cod <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}
build_cod <- function(resource_decrease, capacity_mult) {
  resource_limitation_2d_cod(cod_params, resource_decrease, capacity_mult)
}

make_limit_cycle_sim <- function(params, t_total = 600, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1
  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0, method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 1e3
  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort, method = "predictor-corrector")
}

# The broad scan's day30_broad_scan.csv (Section 3 of 30_experiments.R) is
# an exact log-decade grid (resource_decrease = 1e-7..1, capacity_mult =
# 1e-3..1e2, one grid point per decade), so "the 1e-6 and 1e2 corner" lands
# exactly on a tested point rather than needing interpolation: rd=1e-06,
# cm=100 came back rel_amplitude=1.62e-06, perturbed_verdict="oscillating"
# -- just barely over the project's 1e-6 threshold, and
# mean_late_biomass=13044 (nowhere near the collapse floor). Worth actually
# looking at the trajectory rather than trusting a single borderline number.
resource_decrease <- 1e-6
capacity_mult      <- 100

cat(sprintf("Rebuilding cod at resource_decrease=%.4g, capacity_mult=%.4g...\n",
           resource_decrease, capacity_mult))
p   <- build_cod(resource_decrease, capacity_mult)
sim <- make_limit_cycle_sim(p, t_total = 600)

bm    <- rowSums(getBiomass(sim))
times <- as.numeric(names(bm))
late  <- bm[times >= max(times) * 0.6]
rel_amplitude <- (max(late) - min(late)) / mean(late)
cat(sprintf("Rechecked directly: rel_amplitude = %.6g (broad scan reported 1.62e-06)\n", rel_amplitude))

full_biomass_plot <- plotBiomass(sim)
full_biomass_plot
ggsave(file.path("interesting_plots", "day30_broad_corner_biomass_full.png"),
      plot = full_biomass_plot, width = 9, height = 5, dpi = 150)

# Late window only (matches the window classify_perturbed_run() actually
# scores) -- the full-trajectory plot above can hide a small late
# oscillation behind the much larger early transient from the perturbation
# kick at t=10.
late_biomass_plot <- plotBiomass(sim, tlim = c(360, 600))
late_biomass_plot
ggsave(file.path("interesting_plots", "day30_broad_corner_biomass_late.png"),
      plot = late_biomass_plot, width = 9, height = 5, dpi = 150)

# Size spectrum at the end of the run -- shows whether this corner's
# oscillation (if real) is a whole-spectrum cycle or concentrated in one
# body-size band, same diagnostic Day 22/28 used via plotSpectra().
spectrum_plot <- plotSpectra(sim, time_range = c(560, 600), power = 2)
spectrum_plot
ggsave(file.path("interesting_plots", "day30_broad_corner_spectrum.png"),
      plot = spectrum_plot, width = 7, height = 5, dpi = 150)

cat("Saved: day30_broad_corner_biomass_full.png, day30_broad_corner_biomass_late.png, day30_broad_corner_spectrum.png\n")
