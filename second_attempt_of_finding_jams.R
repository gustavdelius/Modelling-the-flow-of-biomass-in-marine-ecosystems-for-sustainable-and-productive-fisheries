library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse) 
library(glue)

############# Fix that worked - change abundance locally and reduce resource rate ##########################

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
r <- r*0.001
p <- setResource(p,resource_rate = r,resource_dynamics = "resource_semichemostat")
rng <- test_sizes[["large"]]
idx <- p@w >= rng[1] & p@w <= rng[2]
initialN(p)[, idx] <- initialN(p)[, idx] * 5

sim <- project(p, t_max = 100,t_save = 0.5)

animateSpectra(sim,resource=FALSE,log_x=TRUE,log_y=FALSE,power=1.5)



# Deviation from power law
deviation_plot(sim, paste("Deviation — bump at", label, "range"))




