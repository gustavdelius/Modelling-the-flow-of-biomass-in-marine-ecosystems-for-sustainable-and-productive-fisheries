library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse) 
library(glue)

############# Iterating on previous work - gonna vary sigma and beta and (maybe) kappa ##########################
########################### Varying Beta #########################################
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
  
  sim <- project(p, t_max = 100, t_save = 0.5)
  animateSpectra(sim, resource = FALSE, log_x = FALSE, log_y = FALSE, power = 2)
}

varying_beta_sigma_test(sigma_new = 0.65)
a <- varying_beta_sigma_test(beta_new = 30,sigma_new=0.5)

test_sizes <- list(
  "small"  = c(0.1, 1),
  "medium" = c(1, 10),
  "large"  = c(10, 100)
)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
p <- newSingleSpeciesParams(lambda = 2.05, beta = 30, sigma = 0.5)

r <- resource_rate(p)
r <- r * 0.001
p <- setResource(p, resource_rate = r, resource_dynamics = "resource_semichemostat")
rng <- test_sizes[["large"]]
idx <- p@w >= rng[1] & p@w <= rng[2]
initialN(p)[, idx] <- initialN(p)[, idx] * 5

sim <- project(p, t_max = 100, t_save = 0.5)

plotBiomass(sim,start_time=80)
species_params(p)
# Deviation from power law
deviation_plot(sim, paste("Deviation — bump at", label, "range"))

######################## Iterating up on the above ##############
varying_sigma <- function(beta_new=30,sigma_new) {
  p <- newSingleSpeciesParams(lambda = 2.05, beta = beta_new, sigma = sigma_new)
  r <- resource_rate(p)
  r <- r * 0.001
  p <- setResource(p, resource_rate = r, resource_dynamics = "resource_semichemostat")
  idx <- p@w >= 10 & p@w <= 100
  initialN(p)[, idx] <- initialN(p)[, idx] * 5
  project(p, t_max = 100, t_save = 0.5)
}

amplitude_finder <- function(t_max, sims) {
  times    <- as.numeric(dimnames(sims@n)[[1]])
  time_idx <- times >= t_max-10 & times <= t_max
  n_sub    <- sims@n[time_idx, , , drop = FALSE]
  list(
    n_max = apply(n_sub, 3, max) * sims@params@w^2,
    n_min = apply(n_sub, 3, min) * sims@params@w^2
  )
}

sigma_values <- c(0.8,0.75,0.7,0.65,0.6,0.55,0.5)

# Compute amplitude (max - min biomass density) for each sigma
amplitudes <- lapply(sigma_values, function(s) {
  sim <- varying_sigma(s)
  res <- amplitude_finder(t_max = 100, sims = sim)
  list(
    w         = sim@params@w,
    amplitude = res$n_max - res$n_min
  )
})

# Plot
par(mfrow = c(1, 1))

amplitude_against_size_various_sigmas <- function(amplitudes, sigma_values) {
  par(mfrow = c(1, 1))
  cols <- rainbow(length(sigma_values))
  
  # plot largest sigma first so smaller sigma (larger amplitude) draws on top
  y_max <- max(sapply(amplitudes, function(a) max(a$amplitude)))
  plot(amplitudes[[1]]$w, amplitudes[[1]]$amplitude, type = "l",
       log = "x", col = cols[1], lwd = 2,
       ylim = c(0, y_max),
       xlab = "Size (g)", ylab = "Amplitude (max - min biomass density)",
       main = "Effect of sigma on oscillation amplitude")
  for (i in 2:length(sigma_values)) {
    lines(amplitudes[[i]]$w, amplitudes[[i]]$amplitude, col = cols[i], lwd = 2)
  }
  legend("topleft", legend = paste("sigma =", sigma_values),
         col = cols, lty = 1, lwd = 2, cex = 0.7)
}

amplitude_against_size_various_sigmas(amplitudes, sigma_values)



########################Plotting sigma against beta in a phase diagram plot #######

beta_values  <- c(10, 15, 30, 50, 100)
sigma_values <- c(0.3, 0.5,0.6,0.7,0.8, 1.0, 1.3, 1.5)

varying_beta_sigma <- function(beta_new, sigma_new) {
  p <- newSingleSpeciesParams(lambda = 2.05, beta = beta_new, sigma = sigma_new)
  r <- resource_rate(p)
  r <- r * 0.001
  p <- setResource(p, resource_rate = r, resource_dynamics = "resource_semichemostat")
  idx <- p@w >= 10 & p@w <= 100
  initialN(p)[, idx] <- initialN(p)[, idx] * 2.5
  project(p, t_max = 100, t_save = 0.5)
}

max_amplitude <- function(sim) {
  res <- amplitude_finder(t_max = 100, sims = sim)
  n_mean <- (res$n_max + res$n_min) / 2
  max((res$n_max - res$n_min) / n_mean, na.rm = TRUE)
}

# Run all combinations
phase_grid <- expand.grid(beta = beta_values, sigma = sigma_values)
phase_grid$amplitude <- mapply(function(b, s) {
  sim <- varying_beta_sigma(beta_new = b, sigma_new = s)
  max_amplitude(sim)
}, phase_grid$beta, phase_grid$sigma)

ggplot(phase_grid, aes(x = beta, y = sigma, z = amplitude)) +
  geom_contour_filled() +
  scale_x_log10() +
  labs(x = expression(beta), y = expression(sigma),
       fill = "Relative amplitude",
       title = "Phase diagram: oscillation amplitude") +
  theme_minimal()