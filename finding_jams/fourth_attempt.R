library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse) 
library(glue)
################# Testing what Gustav asked me to ###########################
deviation_plot <- function(sim, title) {
  spec <- finalNResource(sim)
  expected <- params@w_full^(-2.05)
  deviation <- spec / expected
  plot(params@w_full, deviation, type = "l", log = "x",
       main = title, xlab = "Size (g)", ylab = "Deviation from power law")
  abline(h = 1, lty = 2)
}

test_sizes <- list(
  "small"  = c(0.1, 1),
  "medium" = c(1, 10),
  "large"  = c(10, 100)
)#size ranges in which we are changing the initial population
par(mfrow = c(5, 2), mar = c(4, 4, 2, 1))  # 3 rows (one per size range), 1 columns (deviation)

varying_beta_sigma_test <- function(beta_new = 100, sigma_new = 1.3) {
  test_sizes <- list(
    "small"  = c(0.1, 1),
    "medium" = c(1, 10),
    "large"  = c(10, 100)
  )
  par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
  p <- newSingleSpeciesParams(lambda = 2.05, beta = beta_new, sigma = sigma_new)
  
  r <- resource_rate(p)
  r <- r * 0.001
  p <- setResource(p, resource_rate = r, resource_dynamics = "resource_semichemostat")
  rng <- test_sizes[["large"]]
  idx <- p@w >= rng[1] & p@w <= rng[2]
  initialN(p)[, idx] <- initialN(p)[, idx] * 5
  
  project(p, t_max = 100, t_save = 0.5)
}



sim <- varying_beta_sigma_test(30,0.5)
#plotHover(getBiomass(sim))
plotHover(getEGrowth(sim))
flux <- getFlux(sim)
biomass_flux <- flux
biomass_flux[] <- flux*w(getParams(sim))
N(sim)[201,1,98:100]
plotHover(biomass_flux)#,wlim=c(10,100))
plot2(biomass_flux,flux,log="xy")
#plot(getFlux(sim))
plotHover(getEncounter(sim))
plot

          