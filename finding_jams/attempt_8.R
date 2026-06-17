library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse) 
library(glue)
library(ggplot2)
################## Code from yesterday that I'll experiment on #################

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

r <- resource_rate(params)
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
sim@n[last, , idx] <- sim@n[last, , idx] /10^3

# Continue for the rest of the run
sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "predictor-corrector")#can use tr_bdf2,but needs much finer dt

#Will make it go on for along time to better see where the attractor  is

sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                   progress_bar = FALSE, method = "predictor-corrector")#can use tr_bdf2,but needs much finer dt

sim_600 <-  project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                    progress_bar = FALSE, method = "predictor-corrector")#can use tr_bdf2,but needs much finer dt

plotHover(getBiomass(sim_600),tlim=c(500,550))
animateSpectra(sim_600,power=2,tlim=c(500,550))

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

dt    <- 0.2   # t_save interval
bm    <- getBiomass(sim_600)[, "Anchovy"]
times <- as.numeric(names(bm))

ps_bm <- power_spectrum(bm[times > 300], dt)

plot_ly(ps_bm |> filter(period > 1, period < 50)) |>
  add_lines(x = ~period, y = ~power) |>
  layout(
    xaxis = list(title = "Period (years)"),
    yaxis = list(title = "Power", type = "log"),
    title = "Power spectrum — Anchovy biomass"
  )

# Feeding level at the smallest size class
fl       <- getFeedingLevel(sim_600)
fl_times <- as.numeric(dimnames(fl)[[1]])
fl_small <- fl[, "Anchovy", 1]

ps_fl <- power_spectrum(fl_small[fl_times > 300], dt)

plot_ly(ps_fl |> filter(period > 1, period < 50)) |>
  add_lines(x = ~period, y = ~power) |>
  layout(
    xaxis = list(title = "Period (years)"),
    yaxis = list(title = "Power", type = "log"),
    title = "Power spectrum — feeding level (smallest size class)"
  )


######################## Iteration 2 ##############################
#Seeing how much size of maturity sets the period size

# Reuses power_spectrum() from above
find_period <- function(bm, times, t_min = 300, dt = 0.2) {
  post <- times > t_min
  ps   <- power_spectrum(bm[post], dt)
  ps$period[which.max(ps$power)]
}

w_mat_vals <- 10^seq(0.5, 1.7, length.out = 12)   # ~3 g to ~50 g

# Note: 12 runs × ~600 years each - takes several minutes
results_wmat <- lapply(w_mat_vals, function(wm) {
  params_wm <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = p$w_min, w_max = p$w_max, w_mat = wm,
    no_w  = no_w,    lambda = p$lambda, kappa = kappa,
    alpha = p$alpha, gamma  = p$gamma,  ks    = 0
  )
  r_wm      <- resource_rate(params_wm) * 0.001
  params_wm <- setResource(params_wm, resource_rate = r_wm,
                           resource_dynamics = "resource_semichemostat")
  params_wm@initial_n_pp[] <- params_wm@cc_pp * 0.1
  
  sim_wm  <- project(params_wm, t_max = 10, dt = 0.1, t_save = 0.2,
                     progress_bar = FALSE, method = "predictor-corrector")
  
  # Crash adults (everything above w_mat)
  idx_wm  <- params_wm@w >= wm & params_wm@w <= p$w_max
  last_wm <- dim(sim_wm@n)[1]
  sim_wm@n[last_wm, , idx_wm] <- sim_wm@n[last_wm, , idx_wm] / 10^3
  
  sim_wm  <- project(sim_wm, t_max = 590, dt = 0.1, t_save = 0.2,
                     progress_bar = FALSE, method = "predictor-corrector")
  
  bm_wm    <- getBiomass(sim_wm)[, "Anchovy"]
  times_wm <- as.numeric(names(bm_wm))
  settled  <- bm_wm[times_wm > 300]
  cv       <- sd(settled) / mean(settled)
  
  data.frame(
    w_mat  = wm,
    period = if (cv > 0.01) find_period(bm_wm, times_wm) else NA,
    cv     = cv
  )
})

df_wmat <- do.call(rbind, results_wmat)

plot_ly(df_wmat[!is.na(df_wmat$period), ]) |>
  add_markers(x = ~w_mat, y = ~period) |>
  layout(
    xaxis = list(type = "log", title = "w_mat (g)"),
    yaxis = list(type = "log", title = "Period (years)"),
    title = "Period vs maturation size"
  )

# Power-law fit: if period ~ w_mat^b, the slope on log-log is b
fit <- lm(log10(period) ~ log10(w_mat), data = df_wmat[!is.na(df_wmat$period), ])
summary(fit)   # slope gives the scaling exponent

###################### Iteration 3 ##########################################
#Cohort tracking

times <- as.numeric(dimnames(sim_600@n)[[1]])

## 3a. Space-time heatmap
# Diagonal bands = cohort wave; horizontal bands = all sizes pulsing together
t_win  <- times > 500 & times < 580   # covers ~4 cycles

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

## 3b. Phase of oscillation across size classes
# Uses period_dom from section 1 — update if needed
omega  <- 2 * pi / period_dom

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


params_cohort <- params
params_cohort@species_params$RDD <- "noRDD"   # kills reproduction

params_cohort@initial_n[] <- 0
params_cohort@initial_n["Anchovy", 1] <- 1e10  # pulse at w_min only

sim_cohort <- project(params_cohort, t_max = 600, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, method = "predictor-corrector")

#Cohort tracking

times <- as.numeric(dimnames(sim_cohort@n)[[1]])

## 3a. Space-time heatmap
# Diagonal bands = cohort wave; horizontal bands = all sizes pulsing together
times  <- as.numeric(dimnames(sim_cohort@n)[[1]])
t_win  <- times >= 0 & times <= 20

bm_win <- sweep(sim_cohort@n[t_win, "Anchovy", ], 2, params@w^2, "*")
df_heat <- melt(bm_win)
names(df_heat) <- c("time", "w", "bm")
df_heat$time <- as.numeric(as.character(df_heat$time))
df_heat$w    <- as.numeric(as.character(df_heat$w))

ggplot(df_heat, aes(x = time, y = log10(w), fill = log10(pmax(bm, 1e-12)))) +
  geom_raster() +
  scale_fill_viridis_c(name = "log10\nbiomass") +
  labs(x = "Time (years)", y = "log10(size [g])",
       title = "Size–time heatmap")
animateSpectra(sim_cohort,tlim=c(0,20))

#Testing to see whether seasonal reproduction can help with this

