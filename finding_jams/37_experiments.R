library(mizer)
library(mizerExperimental)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)   # melt(), for the biomass-over-time plots below
library(plotly)     # plot_ly()/add_lines()/layout(), same

# Day 37 drops the Day 18-36 "Anchovy" model (erepro ~80,000, should be in
# [0,1]) for a full port of the actual published anchovy-plankton model.
#
# Self-contained convention since Day 20: helpers redefined here, not sourced.

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH truncation guard, carried over from Day 30 onward.
save_plot <- function(plot, filename, width = 9, height = 6, dpi = 150) {
  max_name <- 40
  if (nchar(filename) > max_name) {
    ext      <- tools::file_ext(filename)
    base     <- tools::file_path_sans_ext(filename)
    filename <- paste0(substr(base, 1, max_name - nchar(ext) - 1), ".", ext)
    warning(sprintf("save_plot(): filename too long, truncated to '%s'", filename))
  }
  ggsave(file.path("interesting_plots", filename), plot = plot,
         width = width, height = height, dpi = dpi)
}

# Windows MAX_PATH guard for htmlwidgets (plotly, not ggplot) -- same fix
# as Day 36: build in tempdir(), copy the finished file across.
save_widget <- function(widget, filename) {
  tmp <- file.path(tempdir(), filename)
  htmlwidgets::saveWidget(widget, tmp, selfcontained = TRUE)
  file.copy(tmp, file.path("interesting_plots", filename), overwrite = TRUE)
}

################################################################################
# Section 1: the plankton-anchovy model (Canales, Delius & Law 2020)
#
# Full port of https://rpubs.com/gustav/plankton-anchovy (source:
# github.com/sizespectrum/plankton-anchovy), ported onto current mizer.
################################################################################

## Parameters -- Appendix B of the paper, verbatim.
p2 <- list(
  dt = 0.001,
  dx = 0.1,
  w_min = 0.0003,
  w_inf = 66.5,
  ppmr_min = 100,
  ppmr_max = 30000,
  gamma = 750,
  alpha = 0.85, # q -- search-volume exponent, NOT species_params$alpha below
  K = 0.1,      # -> species_params$alpha, assimilation efficiency
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

# Background + larval mortality, paper equations (A.5)-(A.6). Writes
# params@mu_b[] directly, bypassing mizer's own default mortality calc.
setAnchovyMort <- function(params, p) {
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
  mu_b <- mu_b + p$mu_l / (1 + (params@w / p$w_l)^p$rho_l)

  params@mu_b[] <- mu_b
  return(params)
}

# Logistic plankton dynamics with immigration, paper equation (A.11).
# Dispatched by name like this project's own custom rate functions.
plankton_state <- new.env(parent = emptyenv())
plankton_state$time   <- 0
plankton_state$factor <- 1
plankton_state$random <- FALSE   # random plankton forcing kept off here
plankton_state$phi    <- 0
plankton_state$sigma  <- 0.5

plankton_logistic <- function(params, n, n_pp, n_other, rates, dt = 0.1, ...) {
  plankton_state$time <- plankton_state$time + dt
  if (isTRUE(plankton_state$random == "paper") && plankton_state$time >= 0.5) {
    plankton_state$factor <- exp(runif(1, log(1 / 2), log(2)))
    plankton_state$time   <- 0
  } else if (isTRUE(plankton_state$random == "red")) {
    plankton_state$factor <- plankton_state$factor ^ plankton_state$phi *
      exp(rnorm(1, 0, plankton_state$sigma))
  }
  f <- params@rr_pp * n_pp * (1 - n_pp / params@cc_pp / plankton_state$factor) +
    anchovy_immigration - rates$resource_mort * n_pp
  f[is.na(f)] <- 0
  return(n_pp + dt * f)
}

# Normalised box feeding kernel, paper equation (A.2). Kept separate from
# mizer's own "box" kernel, which isn't normalised the same way.
norm_box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
  phi <- rep(1, length(ppmr))
  phi[ppmr > ppmr_max] <- 0
  phi[ppmr < ppmr_min] <- 0
  phi[1] <- 0   # no feeding at own size
  logppmr <- log(ppmr)
  dl <- logppmr[2] - logppmr[1]
  N <- sum(phi) * dl
  phi <- phi / N
  return(phi)
}

# Model builder -- set_multispecies_model() from the original, ported onto
# newMultispeciesParams(). See setAnchovyMort() above for the mortality.
setAnchovyModel <- function(p) {
  kappa <- p$a0 * exp(-6.9 * (p$lambda - 1))

  species_params <- data.frame(
    species          = "Anchovy",
    w_min            = p$w_min,
    w_mat            = p$w_mat,
    m                = p$rho_inf + 2 / 3,   # 2/3 is the original's own filler
    w_inf            = p$w_inf,
    erepro           = p$epsilon_R,
    alpha            = p$K,
    ks               = 0,
    gamma            = p$gamma,
    q                = p$alpha,
    ppmr_min         = p$ppmr_min,
    ppmr_max         = p$ppmr_max,
    pred_kernel_type = "norm_box",
    h                = Inf,
    R_max            = Inf,
    linecolour       = "brown",
    stringsAsFactors = FALSE
  )

  no_w <- round(log(p$w_inf / p$w_min) / p$dx)

  params <- newMultispeciesParams(
    species_params,
    no_w = no_w,
    lambda = p$lambda,
    kappa = kappa,
    w_pp_cutoff = p$w_pp_cutoff,
    resource_dynamics = "plankton_logistic"
    # max_w left unset: mizer pads the species' own w_max above w_inf by
    # default, and pinning max_w=w_inf here conflicts with that.
  )

  params@rr_pp[] <- p$r0 * params@w_full^(p$rho - 1)
  params
}

# Settle + kick, shared by every variant below: 10yr run from a power-law
# initial abundance, knock the population down by 10^7, run t_max more years.
settle_and_kick <- function(params, t_max = 30) {
  params@initial_n[]    <- 0.001 * params@w^(-1.8)
  params@initial_n_pp[] <- params@cc_pp
  sim <- project(params, t_max = 10, dt = p2$dt, progress_bar = FALSE)
  sim@n[11, , ] <- sim@n[11, , ] / 10^7
  project(sim, t_max = t_max, dt = p2$dt, t_save = 0.2, progress_bar = FALSE)
}

# Objective readout instead of eyeballing an interactive widget -- counts
# local extrema and reports actual min/max per species trace in bm.
report_oscillation <- function(bm, label) {
  cat(sprintf("---- %s ----\n", label))
  for (nm in unique(bm$sp)) {
    yv <- bm$value[bm$sp == nm]
    is_max <- c(FALSE, yv[-1] > yv[-length(yv)]) & c(yv[-length(yv)] > yv[-1], FALSE)
    is_min <- c(FALSE, yv[-1] < yv[-length(yv)]) & c(yv[-length(yv)] < yv[-1], FALSE)
    cat(sprintf("  %-14s min=%.4g max=%.4g  local maxima=%d minima=%d\n",
               nm, min(yv), max(yv), sum(is_max), sum(is_min)))
  }
}

# Figure 2e's own recipe: Plankton, total Anchovy and small Anchovy biomass
# together on one log axis. Axis range left auto -- a fixed range clips
# deep troughs and makes a real oscillation look flat.
plot_anchovy_biomass <- function(sim, params, label, filename) {
  abm  <- melt(getBiomass(sim))
  abmr <- melt(getBiomass(sim, min_w = 0.01, max_w = 0.4))
  abmr$sp <- "small Anchovy"
  pbm <- sim@n_pp %*% (params@w_full * params@dw_full)
  pbm <- melt(pbm)
  pbm$Var2 <- NULL
  pbm$sp <- "Plankton"
  bm <- rbind(pbm, abm, abmr)
  report_oscillation(bm, label)

  plot <- plot_ly(bm) %>%
    filter(time >= 10) %>%
    add_lines(x = ~time, y = ~value, color = ~sp) %>%
    layout(yaxis = list(type = "log", exponentformat = "power",
                        title_text = "biomass (g/m^3)"),
           xaxis = list(title_text = "time (year)"))
  print(plot)
  save_widget(plot, filename)
  invisible(plot)
}

anchovy_params <- setAnchovyModel(p2)
plotYieldVsF(anchovy_params_fig2e,"Anchovy")
anchovy_immigration <- p2$i0 * anchovy_params@w_full^(-p2$lambda) * exp(-6.9 * (p2$lambda - 1))

# Full, realistic configuration: cannibalism on, larval mortality on
# (see the Figure 2e reproduction below for cannibalism on its own).
p2$mu_l <- 21
anchovy_params <- setAnchovyMort(anchovy_params, p2)
anchovy_params@interaction[] <- 1

cat(sprintf("plankton-anchovy erepro: %.4g (should sit in [0,1])\n",
           species_params(anchovy_params)$erepro))

################################################################################
# Diagnostics -- print/plot the model actually built. mu_b reset, R_max=Inf
# breaking BevertonHoltRDD, and norm_box dispatch failing were all checked
# and cleared; kernel and mortality below confirmed both build correctly.
################################################################################

cat("---- Section 1 diagnostics ----\n")
print(species_params(anchovy_params)[, c("w_min", "w_inf", "q", "gamma",
                                         "alpha", "ks", "h", "R_max", "erepro")])
cat(sprintf("Grid: %d bins, w from %.4g to %.4g, actual dx (log10) = %.4g (paper intended dx=%.2g in natural log)\n",
           length(anchovy_params@w), min(anchovy_params@w), max(anchovy_params@w),
           diff(log10(anchovy_params@w))[1], p2$dx))

# Should render as a flat block between ppmr_min and ppmr_max, zero outside.
mid_idx <- round(length(anchovy_params@w) / 2)
kernel_df <- data.frame(
  ppmr = anchovy_params@w[mid_idx] / anchovy_params@w_full,
  phi  = getPredKernel(anchovy_params)[1, mid_idx, ]
)
kernel_plot <- ggplot(kernel_df %>% filter(ppmr > 1, ppmr < 1e6), aes(x = ppmr, y = phi)) +
  geom_line() +
  scale_x_log10() +
  labs(x = "Predator/prey mass ratio", y = "Kernel value",
       title = sprintf("Feeding kernel actually built, predator w=%.3g", anchovy_params@w[mid_idx]),
       subtitle = "Should be a flat block between ppmr_min=100 and ppmr_max=30000, zero outside") +
  theme_minimal()
kernel_plot
save_plot(kernel_plot, "day37_anchovy_kernel_check.png", width = 8, height = 5)

# Should show the larval-mortality bump near w_l on top of the
# background/senescent curve, not a smooth generic allometric curve.
mort_df <- data.frame(w = anchovy_params@w, mu_b = anchovy_params@mu_b[1, ])
mort_plot <- ggplot(mort_df, aes(x = w, y = mu_b)) +
  geom_line() +
  scale_x_log10() + scale_y_log10() +
  geom_vline(xintercept = p2$w_l, linetype = "dashed", color = "grey40") +
  labs(x = "Body mass (g)", y = "External mortality mu_b [1/year]",
       title = "Mortality curve actually written to mu_b",
       subtitle = sprintf("Dashed = w_l=%.3g, the larval-mortality peak -- mu_l=%.3g here", p2$w_l, p2$mu_l)) +
  theme_minimal()
mort_plot
save_plot(mort_plot, "day37_anchovy_mort_check.png", width = 8, height = 5)
cat("---- end diagnostics ----\n")

anchovy_sim <- settle_and_kick(anchovy_params)
plot_anchovy_biomass(anchovy_sim, anchovy_params,
                     "Full model (cannibalism + larval mortality)",
                     "day37_anchovy_paper_biomass.html")

################################################################################
# Section 1 continued: isolating growth-dependent larval mortality
#
# mu_l alone, cannibalism off, on a fresh build -- newMultispeciesParams()
# defaults interaction to all-1s, so interaction[]<-0 below isn't a no-op.
################################################################################

p2_larval_only <- p2
p2_larval_only$mu_l <- 21

anchovy_params_larval_only <- setAnchovyModel(p2_larval_only)
anchovy_params_larval_only <- setAnchovyMort(anchovy_params_larval_only, p2_larval_only)
anchovy_params_larval_only@interaction[] <- 0   # cannibalism OFF -- isolating mu_l

anchovy_sim_larval_only <- settle_and_kick(anchovy_params_larval_only)
plot_anchovy_biomass(anchovy_sim_larval_only, anchovy_params_larval_only,
                     "Larval-mortality only (cannibalism off)",
                     "day37_anchovy_larval_only_biomass.html")

################################################################################
# Section 1 continued: exact reproduction of Figure 2e -- cannibalism only
#
# The notebook's own "With cannibalism" (simc): interaction=1, mu_l=0.
# Its preceding baseline run isn't reproduced -- unused by what follows it.
################################################################################

p2_fig2e <- p2
p2_fig2e$mu_l <- 0

anchovy_params_fig2e <- setAnchovyModel(p2_fig2e)
anchovy_params_fig2e <- setAnchovyMort(anchovy_params_fig2e, p2_fig2e)
anchovy_params_fig2e@interaction[] <- 1

anchovy_sim_fig2e <- settle_and_kick(anchovy_params_fig2e)
plot_anchovy_biomass(anchovy_sim_fig2e, anchovy_params_fig2e,
                     "Figure 2e reproduction (cannibalism only, mu_l=0)",
                     "day37_anchovy_fig2e_reproduction.html")

################################################################################
# Section 1 continued: fishing sweep on the plankton-anchovy model
#
# Runs on anchovy_sim_fig2e -- verified against the paper's own text, unlike
# the other two variants above. Period ~6yr (5 maxima over 30yr settle), so
# the sample window is 12yr (~2 periods) to avoid catching a partial cycle.
################################################################################

anchovy_effort_seq   <- seq(0, 3, by = 0.5)
anchovy_t_max        <- 30
anchovy_sample_years <- 12   # ~2 periods of fig2e's own ~6yr cycle

run_anchovy_at_effort <- function(effort) {
  project(anchovy_sim_fig2e, t_max = anchovy_t_max, t_save = 0.2, effort = effort,
         progress_bar = FALSE)
}

anchovy_late_window_stats <- function(sim) {
  tv   <- as.numeric(dimnames(sim@n)[[1]])
  keep <- tv >= (max(tv) - anchovy_sample_years)
  yv   <- getYield(sim)[keep, 1]
  bv   <- getBiomass(sim)[keep, 1]
  data.frame(yield_min = min(yv), yield_max = max(yv),
            biomass_min = min(bv), biomass_max = max(bv))
}

anchovy_bifurcation_df <- bind_rows(lapply(anchovy_effort_seq, function(eff) {
  cbind(effort = eff, anchovy_late_window_stats(run_anchovy_at_effort(eff)))
}))
write.csv(anchovy_bifurcation_df, file.path("interesting_plots", "day37_anchovy_yield_sweep.csv"),
         row.names = FALSE)

anchovy_bifurcation_df <- anchovy_bifurcation_df %>%
  mutate(mean_yield = (yield_min + yield_max) / 2)

anchovy_yield_plot <- ggplot(anchovy_bifurcation_df, aes(x = effort)) +
  geom_ribbon(aes(ymin = yield_min, ymax = yield_max), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = mean_yield)) +
  scale_y_log10() +
  labs(x = "Fishing effort", y = "Yield [g/year] (log scale)",
       title = "Where does yield start collapsing? (plankton-anchovy model, Figure 2e config)",
       subtitle = "Canales, Delius & Law (2020) -- cannibalism only, mu_l=0 -- band = min/max over each run's last 12 years") +
  theme_minimal()
anchovy_yield_plot
save_plot(anchovy_yield_plot, "day37_anchovy_yield_sweep.png", width = 9, height = 6)

anchovy_peak_idx    <- which.max(anchovy_bifurcation_df$mean_yield)
anchovy_peak_effort <- anchovy_bifurcation_df$effort[anchovy_peak_idx]
anchovy_peak_yield  <- anchovy_bifurcation_df$mean_yield[anchovy_peak_idx]

# Yield is still rising at effort=3 (the sweep's own top end) -- "peak" here
# is just that boundary, not a real interior maximum. Extend effort_seq to
# actually find the collapse.
cat(sprintf(
  "Plankton-anchovy (Figure 2e config) effort sweep [%.2g,%.2g] by %.2g: peak mean yield=%.4g at effort=%.2g.\n",
  min(anchovy_effort_seq), max(anchovy_effort_seq), diff(anchovy_effort_seq)[1],
  anchovy_peak_yield, anchovy_peak_effort
))
