library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse)
library(glue)
library(ggplot2)
library(colorRamps)   

# Helper Functions
# Builds the standard single-species Anchovy params.
make_params <- function(lambda = 2.05, resource_decrease = 0.001) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)
  
  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  r      <- getResourceRate(params) * resource_decrease
  params <- setResource(params, resource_rate = r,
                        resource_dynamics = "resource_semichemostat",balance=FALSE)
  params
}

# Runs the two-stage perturbation needed to kick the system onto the limit cycle.
# Stage 1 (t=0–10): run with depleted resource to destabilise.
# Stage 2 (t=10): knock down mature fish by 1000×, then run to t_total.
make_limit_cycle_sim <- function(params, t_total = 600, effort = 0,perturbation=1e3) {
  params@initial_n_pp[] <- params@cc_pp * 0.1
  
  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / perturbation
  
  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "tr_bdf2")
}

#Experiment: Testing Perturbation to see whether or not it's stable

test_perturbation <- function(perturbation, lambda = 2.05, t_total = 600) {
  p   <- make_params(lambda = lambda)   # resource_decrease stays at 0.001
  sim <- make_limit_cycle_sim(p, t_total = t_total, perturbation = perturbation)
  bm  <- getBiomass(sim)[, "Anchovy"]
  t   <- as.numeric(names(bm))
  list(sim = sim, bm = bm, t = t)
}
perturbation <- c(1e3,1e2,1e1,5,1)
results_fishing <- lapply(perturbation, test_perturbation)

# Trajectory panels
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (i in seq_along(results_fishing)) {  # fixed variable name
  r <- results_fishing[[i]]
  plot(r$t, r$bm, type = "l",
       xlab = "Time (years)", ylab = "Biomass",
       main = paste0("Perturbation = ", perturbation[i]))
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Filter plot for t > 550
plot(NULL, xlim = c(550, 600), ylim = range(sapply(results_fishing, function(r) range(r$bm[r$t > 550]))),
     xlab = "Time (years)", ylab = "Biomass", main = "Biomass trajectories (t > 550)")
cols <- colorRampPalette(c("blue", "red"))(length(perturbation))
for (i in seq_along(results_fishing)) {
  r    <- results_fishing[[i]]
  mask <- r$t > 550
  lines(r$t[mask], r$bm[mask], col = cols[i])
}
legend("topright", legend = paste0("Perturb = ", perturbation), col = cols, lty = 1, cex = 0.8)

#Looking at this, but yield against time
test_perturbation <- function(perturbation, lambda = 2.05, t_total = 600, effort = 0.5) {
  p   <- make_params(lambda = lambda,resource_decrease = perturbation)
  sim <- make_limit_cycle_sim(p, t_total = t_total, perturbation = perturbation, effort = effort)
  bm  <- getBiomass(sim)[, "Anchovy"]
  y   <- getYield(sim)[, "Anchovy"]
  t   <- as.numeric(names(bm))
  list(sim = sim, bm = bm, y = y, t = t)
}

perturbation <- c(1e3,1e2,1e1,5,1)
results_fishing <- lapply(perturbation, test_perturbation)

# Trajectory panels
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))

# Trajectory panels
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (i in seq_along(results_fishing)) {
  r <- results_fishing[[i]]
  plot(r$t, r$y, type = "l",
       xlab = "Time (years)", ylab = "Yield",
       main = paste0("Perturbation = ", perturbation[i]))
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Filter plot for t > 550
cols <- colorRampPalette(c("blue", "red"))(length(perturbation))
plot(NULL,
     xlim = c(550, 600),
     ylim = range(sapply(results_fishing, function(r) range(r$y[r$t > 550]))),
     xlab = "Time (years)", ylab = "Yield",
     main = "Yield trajectories (t > 550)")
for (i in seq_along(results_fishing)) {
  r    <- results_fishing[[i]]
  mask <- r$t > 550
  lines(r$t[mask], r$y[mask], col = cols[i])
}
legend("topright", legend = paste0("Perturb = ", perturbation), col = cols, lty = 1, cex = 0.8)

#-------------------------------------------------------------------------------
############ Switch off diffusion and look for hysteresis ######################
#-------------------------------------------------------------------------------

make_params_nodiff <- function(lambda = 2.05, resource_decrease = 0.001) {
  a0    <- 100
  kappa <- a0 * exp(-6.9 * (lambda - 1))
  no_w  <- round(log(66.5 / 0.0003) / 0.1)
  
  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = 0.1, gamma = 750, ks = 0
  )
  r      <- getResourceRate(params) * resource_decrease
  params <- setResource(params, resource_rate = r,
                        resource_dynamics = "resource_semichemostat")
  params <- setExtDiffusion(params,
                                      array(0, dim      = dim(params@ext_diffusion),
                                            dimnames = dimnames(params@ext_diffusion)))
  params
}

# ── Sweep helper (unchanged except uses make_params_nodiff) ───────────────────
run_rd_sweep <- function(rd_seq, init_n = NULL, init_n_pp = NULL,
                         t_run = 300, lambda = 2.05, label = "") {
  out       <- data.frame(rd = rd_seq, amp = NA_real_, bm_mean = NA_real_)
  state_n   <- init_n
  state_npp <- init_n_pp
  
  for (i in seq_along(rd_seq)) {
    p <- make_params_nodiff(lambda = lambda, resource_decrease = rd_seq[i])
    if (!is.null(state_n)) {
      p@initial_n[]    <- state_n
      p@initial_n_pp[] <- state_npp
    }
    sim <- project(p, t_max = t_run, dt = 0.1, t_save = 0.5,
                   progress_bar = FALSE, effort = 0, method = "tr_bdf2")
    bm   <- getBiomass(sim)[, "Anchovy"]
    tv   <- as.numeric(names(bm))
    late <- bm[tv > t_run * 0.6]
    last <- dim(sim@n)[1]
    
    state_n   <- sim@n[last, , ]
    state_npp <- sim@n_pp[last, ]
    
    out$amp[i]     <- max(late) - min(late)
    out$bm_mean[i] <- mean(late)
    cat(sprintf("[%s] rd = %.5f  amp = %.5f\n", label, rd_seq[i], out$amp[i]))
  }
  list(df = out, n_final = state_n, npp_final = state_npp)
}

rd_seq <- seq(0.001, 0.03, by = 0.001)

fwd <- run_rd_sweep(rd_seq, label = "FWD")
bwd <- run_rd_sweep(rev(rd_seq),
                    init_n    = fwd$n_final,
                    init_n_pp = fwd$npp_final,
                    label     = "BWD")

# ── Plot ──────────────────────────────────────────────────────────────────────
fwd_df <- fwd$df
bwd_df <- bwd$df[order(bwd$df$rd), ]

thresh     <- 0.005
fwd_onset  <- fwd_df$rd[which(fwd_df$amp  > thresh)[1]]
bwd_offset <- bwd_df$rd[tail(which(bwd_df$amp > thresh), 1)]

plot(fwd_df$rd, fwd_df$amp,
     type = "l", col = "steelblue", lwd = 2,
     xlab = "resource_decrease",
     ylab = "Oscillation amplitude (max − min biomass)",
     main = "Hysteresis test (no diffusion)",
     ylim = c(0, max(fwd_df$amp, bwd_df$amp) * 1.1))
lines(bwd_df$rd, bwd_df$amp, col = "firebrick", lwd = 2)
abline(v = fwd_onset,  lty = 2, col = "steelblue")
abline(v = bwd_offset, lty = 2, col = "firebrick")
legend("topleft",
       legend = c(
         sprintf("Forward  — onset  rd ≈ %.4f", fwd_onset),
         sprintf("Backward — offset rd ≈ %.4f", bwd_offset)
       ),
       col = c("steelblue", "firebrick"), lwd = 2)


#Seeing how yield varies with respect to resource

make_fishing_params <- function(lambda = 2.05, resource_decrease = 0.001) {
  p  <- make_params(lambda = lambda, resource_decrease = resource_decrease)
  gp <- p@gear_params
  gp$sel_func        <- "knife_edge"
  gp$knife_edge_size <- p@species_params$w_mat
  gp$catchability    <- 1
  gear_params(p)     <- gp
  p
}

get_lc_yield <- function(p, effort, t_total = 600) {
  sim <- make_limit_cycle_sim(p, t_total = t_total, effort = effort)
  y   <- getYield(sim)[, "Anchovy"]
  t   <- as.numeric(rownames(getYield(sim)))
  list(y = y, t = t)
}

get_fp_yield <- function(p, effort) {
  p_fp   <- steadyNewton(p,reproduction="dynamic")
  sim_fp <- project(p_fp, t_max = 1, dt = 0.1, t_save = 1,
                    effort = effort, progress_bar = FALSE)
  getYield(sim_fp)[nrow(getYield(sim_fp)), "Anchovy"]
}

run_yield_sweep <- function(rd_seq, effort, t_run = 600, lambda = 2.05) {
  out <- data.frame(
    rd       = rd_seq,
    lc_mean  = NA_real_,
    lc_max   = NA_real_,
    lc_min   = NA_real_,
    fp_yield = NA_real_
  )
  
  for (i in seq_along(rd_seq)) {
    p       <- make_fishing_params(lambda = lambda, resource_decrease = rd_seq[i])
    lc      <- get_lc_yield(p, effort = effort, t_total = t_run)
    settled <- lc$y[lc$t > 580]
    
    out$lc_mean[i]  <- mean(settled)
    out$lc_max[i]   <- max(settled)
    out$lc_min[i]   <- min(settled)
    out$fp_yield[i] <- get_fp_yield(p, effort = effort)
    
    cat(sprintf("rd = %.5f  lc_mean = %.5f  fp = %.5f\n",
                rd_seq[i], out$lc_mean[i], out$fp_yield[i]))
  }
  out
}

plot_yield_sweep <- function(df) {
  ylim <- range(df$lc_max, df$lc_min, df$fp_yield, na.rm = TRUE)
  
  plot(df$rd, df$lc_mean,
       type = "l", col = "steelblue", lwd = 2,
       xlab = "resource_decrease", ylab = "Yield",
       main = "Yield vs resource_decrease: limit cycle vs fixed point",
       ylim = ylim)
  polygon(c(df$rd, rev(df$rd)),
          c(df$lc_max, rev(df$lc_min)),
          col = adjustcolor("steelblue", alpha.f = 0.2), border = NA)
  lines(df$rd, df$fp_yield, col = "firebrick", lwd = 2)
  legend("topleft",
         legend = c("Limit cycle mean", "Limit cycle range", "Fixed point"),
         col    = c("steelblue", adjustcolor("steelblue", alpha.f = 0.4), "firebrick"),
         lty    = c(1, 1, 1), lwd = c(2, 8, 2))
}

# ── Usage ─────────────────────────────────────────────────────────────────────
rd_seq     <- exp(seq(log(0.001), log(0.03), length.out = 20))
effort_val <- 0.5

sweep <- run_yield_sweep(rd_seq, effort = effort_val)
plot_yield_sweep(sweep)
