library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)

params <- newSingleSpeciesParams(lambda = 2.05)

# Baseline
sim_base <- project(params, t_max = 100)

# Resource rate experiments
params_low_r <- setResource(params, resource_rate = 0.001,resource_dynamics = "resource_semichemostat")

sim_low_r <- project(params_low_r, t_max = 100)

# Growth curve comparison
p1 <- plotGrowthCurves(sim_base) + ggtitle("Baseline")
p2 <- plotGrowthCurves(sim_low_r) + ggtitle("Low resource rate")
wrap_plots(p1, p2)

# Spectra comparison
p3 <- plotSpectra(sim_base, time_range = 90:100) + ggtitle("Baseline")
p4 <- plotSpectra(sim_low_r, time_range = 90:100) + ggtitle("Low resource rate")
wrap_plots(p3,p4)
plotSpectraRelative(sim_base,sim_low_r)

# Deviation from power law
deviation_plot <- function(sim, title) {
  spec <- finalNResource(sim)
  expected <- params@w_full^(-2.05)
  deviation <- spec / expected
  plot(params@w_full, deviation, type = "l", log = "x",
          main = title, xlab = "Size (g)", ylab = "Deviation from power law")
  abline(h = 1, lty = 2)
}

par(mfrow = c(1, 2))
deviation_plot(sim_base, "Baseline")
deviation_plot(sim_low_r, "Low resource rate")

#Animate sim
animate(sim_low_r)
