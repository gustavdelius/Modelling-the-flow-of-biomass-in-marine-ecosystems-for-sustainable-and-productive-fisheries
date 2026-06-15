
library(mizer)
library(patchwork)
library(reshape2)
library(plotly)
library(dplyr)
library(mizerExperimental)
library(tidyverse) 
library(glue)

############# Using a different method in project to see whether phenomenon effect ###############
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
  
  sim <- project(p, t_max = 30, t_save = 0.01)#,method = "predictor-corrector")
  animateSpectra(sim, resource = FALSE, log_x = FALSE, log_y = FALSE, power = 2)
}

varying_beta_sigma_test(beta_new = 30,sigma_new=0.5)

####################### No Longer does oscillate ####################################

########################## Now checking Delius's code ##################

p <- list(
  dt = 0.001,
  dx = 0.1,
  w_min = 0.0003,
  w_inf = 66.5,
  ppmr_min = 100,
  ppmr_max = 30000,
  gamma = 750,
  alpha = 0.85, # q
  K = 0.1, # alpha
  # Larval mortality
  mu_l = 0,
  w_l = 0.03,
  rho_l = 5,
  # background mortality
  mu_0 = 1,
  rho_b = -0.25,
  # Senescent mortality
  w_s = 0.5,
  rho_s = 1,
  # reproduction
  w_mat = 10,
  rho_m = 15,
  rho_inf = 0.2,
  epsilon_R = 0.1,
  # plankton
  w_pp_cutoff = 0.1,
  r0 = 10,
  a0 = 100,
  i0 = 100,
  rho = 0.85,
  lambda = 2
)

setAnchovyMort <- 
  function(params, p) {
    mu_b <- rep(0, length(params@w))
    mu_b[params@w <= p$w_s] <- 
      (p$mu_0 * (params@w / p$w_min)^p$rho_b)[params@w < p$w_s]
    if (p$mu_0 > 0) {
      mu_s <- min(mu_b[params@w <= p$w_s])
    } else {
      mu_s <- p$mu_s
    }
    mu_b[params@w >= p$w_s] <- 
      (mu_s * (params@w / p$w_s)^p$rho_s)[params@w >= p$w_s]
    # Add larval mortality
    mu_b <- mu_b + p$mu_l / (1 + (params@w / p$w_l)^p$rho_l)
    
    params@mu_b[] <- mu_b
    return(params)
  }

plankton_state <- new.env(parent = emptyenv())
plankton_state$time <- 0
plankton_state$factor <- 1
plankton_state$random <- FALSE
plankton_state$phi <- 0
plankton_state$sigma <- 0.5

plankton_logistic <- function(params, n, n_pp, n_other, rates, dt = 0.1, ...) {
  plankton_state$time <- plankton_state$time + dt
  if (plankton_state$random == "paper" && plankton_state$time >= 0.5) {
    # This is the random factor by which we multiply the carrying capacity
    # in the paper, which changes once every six months to a new
    # independent random value
    plankton_state$factor <- exp(runif(1, log(1/2), log(2)))
    plankton_state$time <- 0
  } else if (plankton_state$random == "red") {
    # Here the random factor multiplying the carrying capacity changes
    # at every time step and is given as the exponential of an AR(1)
    # process, i.e., red noise.
    plankton_state$factor <- plankton_state$factor ^ plankton_state$phi * 
      exp(rnorm(1, 0, plankton_state$sigma))
  }
  f <- params@rr_pp * n_pp * (1 - n_pp / params@cc_pp / plankton_state$factor) + 
    i - rates$resource_mort * n_pp 
  f[is.na(f)] <- 0
  return(n_pp + dt * f)
}

norm_box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
  phi <- rep(1, length(ppmr))
  phi[ppmr > ppmr_max] <- 0
  phi[ppmr < ppmr_min] <- 0
  # Do not allow feeding at own size
  phi[1] <- 0
  # normalise in log space
  logppmr <- log(ppmr)
  dl <- logppmr[2] - logppmr[1]
  N <- sum(phi) * dl
  phi <- phi / N
  return(phi)
}

setModel <- function(p) {
  kappa = p$a0 * exp(-6.9*(p$lambda - 1))
  n = 2/3 # irrelevant value
  
  species_params <- data.frame(
    species = "Anchovy",
    w_min = p$w_min,
    w_mat = p$w_mat,
    m = p$rho_inf + n,
    w_inf = p$w_inf,
    erepro = p$epsilon_R,
    alpha = p$K,
    ks = 0,
    gamma = p$gamma,
    ppmr_min = p$ppmr_min,
    ppmr_max = p$ppmr_max,
    pred_kernel_type = "norm_box",
    h = Inf,
    R_max = Inf,
    linecolour = "brown",
    stringsAsFactors = FALSE)
  
  no_w <- round(log(p$w_inf / p$w_min) / p$dx)
  
  params <- set_multispecies_model(
    species_params,
    no_w = no_w,
    lambda = p$lambda,
    kappa = kappa,
    w_pp_cutoff = p$w_pp_cutoff,
    q = p$alpha,
    resource_dynamics = "plankton_logistic")
  
  params@rr_pp[] <- p$r0 * params@w_full^(p$rho - 1)
  return(params)
}

params <- setModel(p)

i <- p$i0 * params@w_full^(-p$lambda) * exp(-6.9*(p$lambda - 1))

p$mu_l <- 0
params <- setAnchovyMort(params, p)
params@interaction[] <- 0

params@initial_n[] <- 0.001 * params@w^(-1.8)
params@initial_n_pp[] <- params@cc_pp
sim <- project(params, t_max = 10, dt = p$dt, progress_bar = FALSE,method = "predictor-corrector")

sim@n[11, , ] <- sim@n[11, , ] / 10^7
sim <- project(sim, t_max = 30, dt = p$dt, t_save = 0.2, progress_bar = FALSE,method = "predictor-corrector")

params@interaction[] <- 1
params@initial_n[] <- 0.001 * params@w^(-1.8)
params@initial_n_pp[] <- params@cc_pp
simc <- project(params, t_max = 10, dt = p$dt, progress_bar = FALSE,method = "predictor-corrector")
simc@n[11, , ] <- simc@n[11, , ] / 10^7
simc <- project(simc, t_max = 30, dt = p$dt, t_save = 0.2, progress_bar = FALSE,method = "predictor-corrector")

nf <- melt(simc@n)
n_ppf <- melt(simc@n_pp)
n_ppf$sp <- "Plankton"
nf <- rbind(nf, n_ppf)

plot_ly(nf) %>%
  # show only part of plankton spectrum
  filter(w > 10^-5) %>% 
  # start at time 20
  filter(time >= 26) %>% 
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
  layout(p, xaxis = list(type = "log", exponentformat = "power",
                         title_text = "body mass (g)"),
         yaxis = list(type = "log", exponentformat = "power",
                      title_text = "biomass (g/m^3)",
                      range = c(-8, 0)))


deviation_plot(sim,"")
animateSpectra(simc,resource=FALSE, log_x = TRUE, log_y = FALSE, power = 2)


############### Hopefully can make the code less weird here ####################
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
sim@n[last, , idx] <- sim@n[last, , idx] / 10^3

# Continue for the rest of the run
sim <- project(sim, t_max = 90, dt = 0.1, t_save = 0.2,
               progress_bar = FALSE, method = "tr_bdf2")#can use tr_bdf2


nf <- melt(sim@n)
n_ppf <- melt(sim@n_pp)
n_ppf$sp <- "Plankton"
nf <- rbind(nf, n_ppf)

plot_ly(nf) %>%
  # show only part of plankton spectrum
  filter(w > 10^-5) %>%
  # start at time 50
  filter(time >= 50) %>%
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

################ 

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
  filter(time >= 10) %>%
  add_lines(x = ~time, y = ~value, color = ~sp) %>%
  # Use logarithmic axes
  layout(yaxis = list(type = "log", exponentformat = "power",
                      title_text = "biomass (g/m^3)",
                      range = c(-7, 2)),
         xaxis = list(title_text = "time (year)"))

animateSpectra(sim,log="xy",power=2)
w_full(params)
w(params)
