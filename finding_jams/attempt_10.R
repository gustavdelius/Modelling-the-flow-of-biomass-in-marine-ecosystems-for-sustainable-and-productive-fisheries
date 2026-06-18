library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)
library(glue)
library(ggplot2)
#Useful Plots 
nice_animation <- function(sim,t){
  nf <- melt(sim@n)
  n_ppf <- melt(sim@n_pp)
  n_ppf$sp <- "Plankton"
  nf <- rbind(nf, n_ppf)
  
  plot_ly(nf) %>%
    # show only part of plankton spectrum
    filter(w > 10^-5) %>%
    # start at time 50
    filter(time >= t) %>%
    # calculate biomass density with respect to log size
    mutate(b = value * w^2) %>%
    # Plot lines
    add_lines(
      x = ~w, y = ~b,
      color = ~sp,
      frame = ~time,
      line = list(simplify = FALSE)
    ) %>%
    # Use logarithmic axes
    layout(xaxis = list(type = "log", exponentformat = "power",
                        title_text = "body mass (g)"),
           yaxis = list(type = "log", exponentformat = "power",
                        title_text = "biomass (g/m^3)",
                        range = c(-8, 0)))
}

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

large_idx    <- params@w_full >= 1e-1 & params@w_full <1e1  # w_full >= 10
r[large_idx] <- r[large_idx] * 0.01

params <- setResource(params, resource_rate = r, resource_dynamics = "resource_semichemostat")
params@initial_n[] <- 0.001 * params@w^(-1.8)

test_sizes <- list(
  "small"  = c(0.1, 1),
  "medium" = c(1, 10),
  "large"  = c(10, 100)
)

# --- Shared phase: perturbation + settling, no fishing ---
sim <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
               effort = 0, progress_bar = FALSE, method = "predictor-corrector")

# rng <- test_sizes[["large"]]
# idx <- params@w >= rng[1] & params@w <= rng[2]
# last <- dim(sim@n)[1]
# sim@n[last, , idx] <- sim@n[last, , idx] / 10^3

sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               effort = 0, progress_bar = FALSE, method = "predictor-corrector")

sim_300 <- project(sim, t_max = 200, dt = 0.1, t_save = 0.2,
                   effort = 0, progress_bar = FALSE, method = "predictor-corrector")

# --- Branch 1: no fishing (this is your existing pipeline) ---
sim_600 <- project(sim_300, t_max = 300, dt = 0.1, t_save = 0.2,
                   effort = 0, progress_bar = FALSE, method = "predictor-corrector")


plotHover(getBiomass(sim_600),tlim=c(550,600))
plotHover(getFlux(sim_600,power=2),tlim=c(550,600),log="xy")
animateSpectra(sim_600,power=2,tlim=c(550,600),resource = FALSE,wlim=c(2e-3,2e2))


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
omega  <- 2 * pi / 7.89

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



# --- Branch 2: fishing targeted at the bottleneck size range ---
# Bottleneck = juveniles below maturation (w_mat = 10g), where the
# Day 8 phase analysis put the food-limitation queue (amplitude peak ~0.01-0.05g).

bottleneck_sel <- function(w, w_low, w_high) {
  as.numeric(w >= 10 & w <= 100)
}

params_fished <- sim_300@params
params_fished@gear_params$sel_func     <- "bottleneck_sel"
params_fished@gear_params$w_low        <- 0.01
params_fished@gear_params$w_high       <- p$w_mat
params_fished@gear_params$catchability <- 1

sim_300_fished        <- sim_300
sim_300_fished@params <- params_fished

sim_fished <- project(sim_300_fished, t_max = 300, dt = 0.1, t_save = 0.2,
                      effort = 0.3, progress_bar = FALSE,
                      method = "predictor-corrector")

# --- Compare yield ---
plotYield(sim_fished, sim_600, tlim = c(550, 600))

y_fished <- getYield(sim_fished)
y_base   <- getYield(sim_600)
t_f      <- as.numeric(rownames(y_fished))
t_b      <- as.numeric(rownames(y_base))

mean(y_fished[t_f > 550, ])   # average yield, fishing the bottleneck
mean(y_base[t_b > 550, ])     # average yield, no fishing (will be 0)


# --- Branch 3: fish only at peaks, never at troughs ---

# 1. Find peak times from the settled, unfished trajectory (sim_600),
#    reusing sim_600's own dimnames directly to avoid floating-point
#    drift from re-parsing/re-formatting times, and including the
#    fork point t = 300 itself (not just times strictly greater than it)
all_times <- as.numeric(dimnames(sim_600@n)[[1]])
keep      <- all_times >= 300
t_p       <- all_times[keep]
time_lbls <- dimnames(sim_600@n)[[1]][keep]

bm   <- getBiomass(sim_600)[, "Anchovy"]
bm_p <- bm[time_lbls]

is_peak    <- c(FALSE, diff(sign(diff(bm_p))) == -2, FALSE)
peak_times <- t_p[is_peak]

# 2. Build a 0/fish_level effort schedule: "on" only within `window`
#    years either side of each peak, "off" everywhere else (troughs included)
window     <- 3.4    # years either side of peak to fish
fish_level <- 0.3

effort_vec <- rep(0, length(t_p))
for (pt in peak_times) {
  effort_vec[t_p >= pt - window / 2 & t_p <= pt + window / 2] <- fish_level
}

gear_name  <- params_fished@gear_params$gear[1]
effort_arr <- matrix(effort_vec, ncol = 1,
                     dimnames = list(time_lbls, gear_name))

# 3. Run the bottleneck-gear branch with this schedule instead of constant effort
sim_peaks_only <- project(sim_300_fished, t_max = 300, dt = 0.1, t_save = 0.2,
                          effort = effort_arr, progress_bar = FALSE,
                          method = "predictor-corrector")

plotYield(sim_peaks_only, sim_fished, tlim = c(550, 600))
y_fished <- getYield(sim_peaks_only)
y_base   <- getYield(sim_fished)
t_f      <- as.numeric(rownames(y_fished))
t_b      <- as.numeric(rownames(y_base))

mean(y_fished[t_f > 550, ])   # average yield, fishing the bottleneck
mean(y_base[t_b > 550, ])     # average yield, no fishing (will be 0)
plotHover(getBiomass(sim_600),tlim=c(550,600))
plotHover(getBiomass(sim_peaks_only),tlim=c(550,600))
plotRelative(getBiomass(sim_peaks_only),getBiomass(sim_600),tlim=c(550,600))
plotRelative(getFlux(sim_600,power=2),getFlux(sim_peaks_only,power=2),tlim=c(550,600))
plotHover(getFlux(sim_600,power=2),tlim=c(550,600))
plotHover(getFlux(sim_peaks_only,power=2),tlim=c(550,600))
animateSpectra(sim_peaks_only,power=2,tlim=c(550,600),resource = FALSE,wlim=c(2e-3,2e2))

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
bm    <- getBiomass(sim_peaks_only)[, "Anchovy"]
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
fl       <- getFeedingLevel(sim_peaks_only)
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

