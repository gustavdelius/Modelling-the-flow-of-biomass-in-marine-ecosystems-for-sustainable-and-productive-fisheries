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





for (label in names(test_sizes)) {
  p <- newSingleSpeciesParams(lambda = 2.05)
  p <- setResource(p,resource_dynamics = "resource_semichemostat")
  
  rng <- test_sizes[[label]]
  idx <- p@w >= rng[1] & p@w <= rng[2]
  initialN(p)[, idx] <- initialN(p)[, idx] * 5
  
  sim <- project(p, t_max = 6,t_save = 0.2)
  
  
  # Deviation from power law
  deviation_plot(sim, paste("Deviation — bump at", label, "range"))
  animateSpectra(sim,resource=FALSE,log_x=FALSE,log_y=FALSE,power=2)
  n_pp <- initialNResource(params)
  n_pp[] <- NResource(sim)[30,]

  growth_rate <- getEGrowth(params,n_pp=n_pp)
  #plot(growth_rate, log = "")
  
  n <- initialN(params)
  
  growth_rate_2 <- getEGrowth(params)
  #plot(growth_rate, log = "")
  growth_rate-growth_rate_2
  plot2(growth_rate,growth_rate_2,log="")
}

# resource_rates <- c(0.001, 0.0001, 0.00001, 0.000001, 0.0000001)
# 
# plot_list <- list()
# 
# for (rate in resource_rates) {
#   p <- newSingleSpeciesParams(lambda = 2.05)
#   p <- setResource(p, resource_rate = rate,
#                    resource_dynamics = "resource_semichemostat")
#   
#   idx <- p@w >= 10 & p@w <= 100
#   initialN(p)[, idx] <- initialN(p)[, idx] * 5
#   
#   t_max <- min(ceiling(5 / rate), 2000)
#   sim <- project(p, t_max = t_max, dt = 1)
#   
#   w_full <- sim@params@w_full
#   deviation <- finalNResource(sim) / w_full^(-2.05)
#   
#   plot_list[[as.character(rate)]] <- data.frame(
#     w         = w_full,
#     deviation = deviation,
#     rate      = as.character(rate)
#   )
# }
# 
# bind_rows(plot_list) |>
#   ggplot(aes(x = w, y = deviation, colour = rate)) +
#   geom_line() +
#   scale_x_log10(limits = c(1e-3, 1e3)) +
#   labs(x = "Size (g)", y = "Deviation from power law",
#        colour = "Resource rate",
#        title = "Resource depletion vs resource rate (large bump)") +
#   theme_minimal()


#resource_rates <- c(0.001, 0.0001, 0.00001, 0.000001, 0.0000001)
resource_rates <- c(0.001)

for (rate in resource_rates) {
  p <- newSingleSpeciesParams(lambda = 2.05)
  p <- setResource(p, resource_rate = rate,
                   resource_dynamics = "resource_semichemostat")
  sim_no_bump <- project(p, t_max = 2000, dt = 1)
  idx <- p@w >= 10 & p@w <= 100
  initialN(p)[, idx] <- initialN(p)[, idx] * 5
  
  sim_bump    <- project(p,    t_max = 2000, dt = 1)
  
  
  print(plotSpectraRelative(sim_no_bump, sim_bump))
  
  
}

# times <- as.numeric(dimnames(sim_no_bump@n)[[1]])
# 
#  
# p_base <- newSingleSpeciesParams(lambda = 2.05)
# p_base <- setResource(p_base, resource_rate = 0.001,
#                       resource_dynamics = "resource_semichemostat")
# 
# p_bump <- p_base
# idx <- p_bump@w >= 10 & p_bump@w <= 100
# initialN(p_bump)[, idx] <- initialN(p_bump)[, idx] * 5
# 
# sim_no_bump <- project(p_base, t_max = 2000, dt = 1)
# sim_bump    <- project(p_bump, t_max = 2000, dt = 1)
# 
# n_base <- sim_no_bump@n
# n_bump <- sim_bump@n
# times  <- as.numeric(dimnames(sim_no_bump@n)[[1]])
# w      <- p_base@w
# 
# bump_idx <- w >= 10 & w <= 100
# 
# rel_diff_over_time <- sapply(seq_along(times), function(t) {
#   base_t <- n_base[t, 1, bump_idx]
#   bump_t <- n_bump[t, 1, bump_idx]
#   mean((bump_t - base_t) / (base_t + 1e-10))
# })
# 
# p1 <- data.frame(time = times, rel_diff = rel_diff_over_time) |>
#   ggplot(aes(x = time, y = rel_diff)) +
#   geom_line() +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "Time", y = "Mean relative difference",
#        title = "Bump at 10–100g: does it persist or dissipate?") +
#   theme_minimal()
# 
# time_idx <- seq(1, length(times), by = 20)
# 
# heatmap_data <- do.call(rbind, lapply(time_idx, function(ti) {
#   base_t <- n_base[ti, 1, ]
#   bump_t <- n_bump[ti, 1, ]
#   data.frame(
#     time     = times[ti],
#     w        = w,
#     rel_diff = (bump_t - base_t) / (base_t + 1e-10)
#   )
# }))
# 
# p2 <- heatmap_data |>
#   filter(w >= 0.1, w <= 1000) |>
#   ggplot(aes(x = w, y = time, fill = rel_diff)) +
#   geom_tile() +
#   scale_x_log10() +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   labs(x = "Size (g)", y = "Time", fill = "Relative diff",
#        title = "Relative difference across size and time") +
#   theme_minimal()
# 
# p1 / p2