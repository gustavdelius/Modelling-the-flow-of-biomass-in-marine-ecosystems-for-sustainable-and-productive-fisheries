library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)
library(glue)
library(ggplot2)
################## Way to calclulate the age #################

p <- list(
  dt     = 0.001,
  dx     = 0.1,
  w_min  = 0.0003,
  w_max  = 66.5,     # was w_inf
  w_mat  = 10,
  alpha  = 0.1,      # assimilation efficiency (was p$K)
  gamma  = 750,
  lambda = 2,
  a0     = 100
)

kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))
no_w  <- round(log(p$w_max / p$w_min) / p$dx)

params <- newSingleSpeciesParams(
  species_name = "Anchovy",
  w_min  = p$w_min,
  w_max  = p$w_max,
  w_mat  = p$w_mat,
  no_w   = no_w,
  lambda = p$lambda,
  kappa  = kappa,
  alpha  = p$alpha,
  gamma  = p$gamma,
  ks     = 0
)

r <- getResourceRate(params)
r <- r * 0.001
params <- setResource(params, resource_rate = r, resource_dynamics = "resource_semichemostat")

test_sizes <- list(
  "small"  = c(0.1, 1),
  "medium" = c(1, 10),
  "large"  = c(10, 100)
)

params@initial_n_pp[] <- params@cc_pp * 0.1

# Run up to T1
sim <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

# Select the size range to perturb
rng <- test_sizes[["large"]]
idx <- params@w >= rng[1] & params@w <= rng[2]

# Divide that size range, at the current end of sim, by 10^3
last <- dim(sim@n)[1]
sim@n[last, , idx] <- sim@n[last, , idx] / 10^3

# Continue for the rest of the run
sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

# Growth rate at each size, from a params object (steady-state)
g  <- getEGrowth(params)["Anchovy", ]
w  <- params@w
dw <- params@dw

# age(w) = integral_{w_min}^{w} 1 / g(w') dw'  (cumulative trapezoidal rule)
integrand <- 1 / g
integrand
integrand[-1]
age <- c(0, cumsum((integrand[-1] + integrand[-length(integrand)]) / 2 * dw[-length(dw)]))
age_at_mat <- approx(x = w, y = age, xout = params@species_params$w_mat)$y
age_at_mat
times <- as.numeric(dimnames(sim_600@n)[[1]])

bm   <- getBiomass(sim_600)[, "Anchovy"]

period_dom <- find_period(bm, times)
omega      <- 2 * pi / period_dom
1/(age_at_mat/period_dom)
df <- data.frame(w = head(w, -1), age = head(age, -1))

ggplot(df, aes(x = w, y = age)) +
  geom_line() +
  geom_vline(xintercept = params@species_params$w_mat, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = age_at_mat, linetype = "dashed", colour = "grey50") +
  geom_point(aes(x = params@species_params$w_mat, y = age_at_mat), colour = "red", size = 2) +
  scale_x_log10() +
  labs(x = "body mass (g)", y = "age (years)",
       title = "Age-at-size, Anchovy")

p <- list(
  dt     = 0.001,
  dx     = 0.1,
  w_min  = 0.0003,
  w_max  = 66.5,     # was w_inf
  w_mat  = 10,
  alpha  = 0.1,      # assimilation efficiency (was p$K)
  gamma  = 750,
  lambda = 2.05,
  a0     = 100
)

kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))
no_w  <- round(log(p$w_max / p$w_min) / p$dx)

params <- newSingleSpeciesParams(
  species_name = "Anchovy",
  w_min  = p$w_min,
  w_max  = p$w_max,
  w_mat  = p$w_mat,
  no_w   = no_w,
  lambda = p$lambda,
  kappa  = kappa,
  alpha  = p$alpha,
  gamma  = p$gamma,
  ks     = 0
)

r <- getResourceRate(params)
r <- r * 0.001
params <- setResource(params, resource_rate = r, resource_dynamics = "resource_semichemostat")

test_sizes <- list(
  "small"  = c(0.1, 1),
  "medium" = c(1, 10),
  "large"  = c(10, 100)
)

params@initial_n_pp[] <- params@cc_pp * 0.1

# Run up to T1
sim <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

# Select the size range to perturb
rng <- test_sizes[["large"]]
idx <- params@w >= rng[1] & params@w <= rng[2]

# Divide that size range, at the current end of sim, by 10^3
last <- dim(sim@n)[1]
sim@n[last, , idx] <- sim@n[last, , idx] / 10^3

# Continue for the rest of the run
sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

#Fourier transform work

# Reusable helper: takes a signal + sample interval, returns period/power df
# Removes DC before transforming, skips the zero-frequency bin
power_spectrum <- function(signal, dt) {
  signal <- signal - mean(signal)
  n      <- length(signal)
  freqs  <- (1:(n - 1)) / (n * dt)          # skip DC (freq = 0)
  power  <- Mod(fft(signal))^2
  half   <- seq_len(floor((n - 1) / 2))     # positive frequencies only
  data.frame(freq = freqs[half], period = 1 / freqs[half], power = power[half + 1])
}

######################## Iteration 2 ##############################
#Seeing how much size of maturity sets the period size

# Reuses power_spectrum() from above
find_period <- function(bm, times, t_min = 300, dt = 0.2) {
  post <- times > t_min
  ps   <- power_spectrum(bm[post], dt)
  ps$period[which.max(ps$power)]
}


###################### Using a smaller period of time, so that the heatmap is more informative ##########

w_mat_vals <- 10^seq(0.5, 1.7, length.out = 12)   # ~3 g to ~50 g

times <- as.numeric(dimnames(sim_600@n)[[1]])

## 3a. Space-time heatmap
# Diagonal bands = cohort wave; horizontal bands = all sizes pulsing together
t_win  <- times > 560 & times < 580   # covers ~4 cycles

bm_win <- sweep(sim_600@n[t_win, "Anchovy", ], 2, params@w^2, "*")

df_heat <- melt(bm_win)
names(df_heat) <- c("time", "w", "bm")
df_heat$time <- as.numeric(as.character(df_heat$time))
df_heat$w    <- as.numeric(as.character(df_heat$w))

ggplot(df_heat, aes(x = time, y = log10(w), fill = log10(pmax(bm, 1e-12)))) +
  geom_raster() +
  scale_fill_viridis_c(name = "log10\nbiomass") +
  labs(x = "Time (years)", y = "log10(size [g])",
       title = "Size–time heatmap")

############## Illustrates the weirdness at 2.05 ############################
bm   <- getBiomass(sim_600)[, "Anchovy"]

period_dom <- find_period(bm, times)
omega      <- 2 * pi / period_dom
1/(age_at_mat/period_dom)
post   <- times > 300
t_post <- times[post]
t_vec  <- t_post - t_post[1]
n_post <- sim_600@n[post, "Anchovy", ]

# DFT coefficient at the dominant frequency for each size bin
coeffs <- apply(n_post, 2, function(x) {
  signal <- x - mean(x)
  sum(signal * exp(-1i * omega * t_vec))
})


df_phase <- data.frame(
  w         = params@w,
  amplitude = Mod(coeffs),
  phase     = Arg(coeffs)
)


p_amp <- ggplot(df_phase, aes(x = w, y = amplitude)) +
  geom_line() + scale_x_log10() + scale_y_log10() +
  labs(x = "Size (g)", y = "Amplitude",
       title = "Oscillation amplitude by size class")

p_phs <- ggplot(df_phase, aes(x = w, y = phase)) +
  geom_line() + scale_x_log10() +
  labs(x = "Size (g)", y = "Phase (radians)",
       title = "Oscillation phase by size class")

p_amp + p_phs


#################### Trying to find parameters where the oscillations are most "sine" like
nice_biomass_plot <- function(sim,t){
  abm <- melt(getBiomass(sim))
  abmr <- melt(getBiomass(sim, min_w = 0.01, max_w = 0.4))
  abmr$sp <- "small Anchovy"
  pbm <- sim@n_pp %*% (params@w_full * params@dw_full)
  pbm <- melt(pbm)
  names(pbm)[names(pbm) == "Var1"] <- "time"
  pbm$Var2 <- NULL
  pbm$sp <- "Plankton"
  bm <- rbind(pbm, abm, abmr)
  plot_ly(bm) %>%
    filter(time >= t) %>%
    add_lines(x = ~time, y = ~value, color = ~sp) %>%
    # Use logarithmic axes
    layout(yaxis = list(type = "log", exponentformat = "power",
                        title_text = "biomass (g/m^3)",
                        range = c(-7, 2)),
           xaxis = list(title_text = "time (year)"))
}


p <- list(
  dt     = 0.001,
  dx     = 0.1,
  w_min  = 0.0003,
  w_max  = 66.5,     # was w_inf
  w_mat  = 10,
  alpha  = 0.1,      # assimilation efficiency (was p$K)
  gamma  = 750,
  lambda = 2.05,
  a0     = 100
)

kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))
no_w  <- round(log(p$w_max / p$w_min) / p$dx)

params <- newSingleSpeciesParams(
  species_name = "Anchovy",
  w_min  = p$w_min,
  w_max  = p$w_max,
  w_mat  = p$w_mat,
  no_w   = no_w,
  lambda = p$lambda,
  kappa  = kappa,
  alpha  = p$alpha,
  gamma  = p$gamma,
  ks     = 0
)

r <- getResourceRate(params)
r <- r * 0.01
params <- setResource(params, resource_rate = r, resource_dynamics = "resource_semichemostat")

test_sizes <- list(
  "small"  = c(0.1, 1),
  "medium" = c(1, 10),
  "large"  = c(10, 100)
)

params@initial_n_pp[] <- params@cc_pp * 0.1

# Run up to T1
sim <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

# Select the size range to perturb
rng <- test_sizes[["large"]]
idx <- params@w >= rng[1] & params@w <= rng[2]

# Divide that size range, at the current end of sim, by 10^3
last <- dim(sim@n)[1]
sim@n[last, , idx] <- sim@n[last, , idx] / 10^3

# Continue for the rest of the run
sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")

sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")

# nice_biomass_plot(sim_600,550)
# plotHover(getBiomass(sim_600),tlim=c(550,600))
times <- as.numeric(dimnames(sim_600@n)[[1]])

period_dom <- find_period(bm, times)
omega      <- 2 * pi / period_dom

post   <- times > 300
t_post <- times[post]
t_vec  <- t_post - t_post[1]
n_post <- sim_600@n[post, "Anchovy", ]

# DFT coefficient at the dominant frequency for each size bin
coeffs <- apply(n_post, 2, function(x) {
  signal <- x - mean(x)
  sum(signal * exp(-1i * omega * t_vec))
})


df_phase <- data.frame(
  w         = params@w,
  amplitude = Mod(coeffs),
  phase     = Arg(coeffs)
)


p_amp <- ggplot(df_phase, aes(x = w, y = amplitude)) +
  geom_line() + scale_x_log10() + scale_y_log10() +
  labs(x = "Size (g)", y = "Amplitude",
       title = "Oscillation amplitude by size class")

p_phs <- ggplot(df_phase, aes(x = w, y = phase)) +
  geom_line() + scale_x_log10() +
  labs(x = "Size (g)", y = "Phase (radians)",
       title = "Oscillation phase by size class")

p_amp + p_phs
