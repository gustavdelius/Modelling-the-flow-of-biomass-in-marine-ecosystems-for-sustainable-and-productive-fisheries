library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse) 
library(glue)
####################### First Improvement - Just decrease the resource rate by a lot more ########################
# params <- newSingleSpeciesParams(lambda = 2.05)
# # Knock the initial population off steady state
# 
# # Resource rate experiments
# params_low_r <- setResource(params, resource_rate = 0.0000001,resource_dynamics = "resource_semichemostat")
# 
# initialN(params) <- initialN(params) * 0.01
# initialN(params_low_r) <- initialN(params_low_r) * 0.01
# 
# # Baseline
# sim_base <- project(params, t_max = 100)
# 
# #Low Resource rate
# sim_low_r <- project(params_low_r, t_max = 100)
# 
# #SpectraComparison
# plotSpectraRelative(sim_base,sim_low_r)
# 
# # Deviation from power law
# deviation_plot <- function(sim, title) {
#   spec <- finalNResource(sim)
#   expected <- params@w_full^(-2.05)
#   deviation <- spec / expected
#   plot(params@w_full, deviation, type = "l", log = "x",
#        main = title, xlab = "Size (g)", ylab = "Deviation from power law")
#   abline(h = 1, lty = 2)
# }
# 
par(mfrow = c(1, 3))
# deviation_plot(sim_base, "Baseline")
# deviation_plot(sim_low_r, "Low resource rate")
# 
# #Animate sim
# animateSpectra(sim_low_r)

########### Results - all it does is make the peaks much more extreme - not very informative ###



############# Fix 2 - change the abundance more locally ##########################

# Deviation from power law
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
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))  # 3 rows (one per size range), 1 columns (deviation)



#no longer need to loop
#Got Osiciliatory motion!!!!!!!!!!!!!
p <- newSingleSpeciesParams(lambda = 2.05)

r <- resource_rate(p)
r <- r*0.0001
p <- setResource(p,resource_rate = r,resource_dynamics = "resource_semichemostat")
rng <- test_sizes[["large"]]
idx <- p@w >= rng[1] & p@w <= rng[2]
initialN(p)[, idx] <- initialN(p)[, idx] * 5

sim <- project(p, t_max = 100,t_save = 0.5)


# Deviation from power law
deviation_plot(sim, paste("Deviation — bump at", label, "range"))
animateSpectra(sim,resource=FALSE,log_x=FALSE,log_y=FALSE,power=2)

for (label in names(test_sizes)) {
  p <- newSingleSpeciesParams(lambda = 2.05)
  p <- setResource(p,resource_dynamics = "resource_semichemostat")
  
  rng <- test_sizes[[label]]
  idx <- p@w >= rng[1] & p@w <= rng[2]
  initialN(p)[, idx] <- initialN(p)[, idx] * 5
  
  sim <- project(p, t_max = 10,t_save = 0.5)
  
  
  # Deviation from power law
  deviation_plot(sim, paste("Deviation — bump at", label, "range"))
  animateSpectra(sim,resource=FALSE,log_x=FALSE,log_y=FALSE,power=2)
  n_pp <- initialNResource(params)
  n_pp[] <- NResource(sim)[15,]

  growth_rate <- getEGrowth(params,n_pp=n_pp)
  #plot(growth_rate, log = "")
  
  n <- initialN(params)
  
  growth_rate_2 <- getEGrowth(params)
  #plot(growth_rate, log = "")
  growth_rate-growth_rate_2
  plot2(growth_rate,growth_rate_2,log="")
}


