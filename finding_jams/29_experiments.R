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
library(future)
library(future.apply)
library(scales)
library(deSolve)

# check_hybrid_point()/check_hopf_grid_point_generic() calls are independent
# getStability() solves, same as Day 28's grid -- run across workers rather
# than sequentially.
plan(multisession)

dir.create("interesting_plots", showWarnings = FALSE)

# Windows MAX_PATH (260 chars) has silently truncated filenames -- even the
# .png extension itself -- once combined with this repo's long, deeply
# nested folder path. Keep names short at the call site; this is a
# defensive last resort so a long one truncates safely instead of silently.
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 150) {
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

################################################################################
# Where Day 28 left us
#
# Day 28 found cod_params stable everywhere tested (spectral_radius pinned
# near 0.51-0.76 across a 120-point resource_decrease x capacity_mult grid),
# and confirmed both cod and Day 22's anchovy parameters are "adult-driven"
# by de Roos & Persson's (2003) mass-specific-ingestion-rate criterion --
# adults out-eat juveniles in both. That's a poor match for the earlier
# read of the paper's asymmetry-alone story: if raw distance from q=1 (or
# from de Roos & Persson's "adults vs juveniles" sign) predicted cycling,
# cod -- whose adult/juvenile imbalance is ~44x starker than anchovy's --
# should be the MORE cycle-prone of the two, not the one that stayed flat.
#
# A follow-up read of the de Roos & Persson lineage (specifically Soudijn,
# de Roos et al. 2017, Theor Ecol 10:73-90, which reuses the 2003 mechanism
# in a stage-structured biomass model) points at two things the raw
# asymmetry sign alone misses:
#   1. Adult-driven cycling additionally requires a *food-dependent
#      maturation delay* -- it vanishes if juvenile development time is
#      held fixed. Both mizer models already have food-dependent growth,
#      so this alone can't distinguish cod from anchovy.
#   2. Adult-driven cycling is fragile to *body size*: in their sweep,
#      juvenile-driven cycles (q<1) persisted across the whole tested body
#      size range, but adult-driven cycles (q>1) disappeared at small body
#      size. Their parameterisation makes the reason legible: consumer
#      rates scale allometrically (~W_A^-0.25) while the resource's own
#      turnover rate doesn't scale with consumer body size at all -- so
#      body size is really a proxy for the RATIO of resource-turnover rate
#      to consumer rate, i.e. a timescale-separation condition.
#
# Today tests that timescale-separation reading directly, in three steps:
#   Section 1 -- compute a timescale-separation number for both cod_params
#     and the anchovy calibration and just compare them.
#   Section 2 -- run the exact getStability() grid Day 28 used on cod, on
#     anchovy too, so the two are compared with identical machinery instead
#     of Day 22-27's simulation/catchable_fraction proxies (which Day 27
#     already found gave an ambiguous answer for anchovy's own reference
#     point: rel_amplitude ~1e-12, i.e. flat, at capacity_mult=10,
#     resource_decrease=0.001).
#   Section 3 -- a body-size x intake-rate-parameter 2x2 factorial, built
#     inside the anchovy single-species template, mirroring the "cross
#     species params with gamma alone" 2x2 design Day 19/20 already used
#     for isolating a different lever. This is the section that actually
#     answers "is it body size, is it gamma/alpha/ks, or both" rather than
#     just re-confirming that cod and anchovy differ somehow.
################################################################################

################################################################################
# Section 0: Rebuild both species objects, self-contained (same convention
# every script since Day 20 has used -- redefined here rather than sourced
# from an earlier file).
################################################################################

# Read-only: cod_params.rds is never written back to.
cod_params <- readRDS("cod_params.rds")

# Identical to Day 27's make_second_order_params_kr(), renamed for this
# file. Defaults reproduce the capacity_mult=10, resource_decrease=0.001
# reference operating point every juvenile-pileup sweep has used since
# Day 24.
make_anchovy_params <- function(lambda = 2.05, resource_decrease = 0.01,
                                capacity_mult = 10, second_order = TRUE,
                                ext_diff = 0.00, alpha = 0.1, gamma = 750,
                                ks = 0, kappa_override = NULL) {
  a0    <- 100
  kappa <- if (is.null(kappa_override)) a0 * exp(-6.9 * (lambda - 1)) else kappa_override
  no_w  <- round(log(66.5 / 0.0003) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = "Anchovy",
    w_min = 0.0003, w_max = 66.5, w_mat = 10,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = alpha, gamma = gamma, ks = ks
  )

  given_species_params(params)$D_ext <- ext_diff
  params <- setBevertonHolt(params)

  default_capacity <- getResourceCapacity(params)
  r  <- getResourceRate(params) * resource_decrease
  cc <- default_capacity * capacity_mult

  params <- setResource(params, resource_rate = r, resource_capacity = cc,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE)

  params
}
make_limit_cycle_sim <- function(params, t_total = 600, effort = 0) {
  params@initial_n_pp[] <- params@cc_pp * 0.1
  
  sim_init <- project(params, t_max = 10, dt = 0.1, t_save = 0.2,
                      progress_bar = FALSE, effort = 0,
                      method = "predictor-corrector")
  idx  <- params@w >= 10 & params@w <= 100
  last <- dim(sim_init@n)[1]
  sim_init@n[last, , idx] <- sim_init@n[last, , idx] / 1e3
  
  project(sim_init, t_max = t_total - 10, dt = 0.1, t_save = 0.2,
          progress_bar = FALSE, effort = effort,
          method = "predictor-corrector")
}

anchovy_params <- make_anchovy_params()
anch_sim <- make_limit_cycle_sim(anchovy_params)
plotBiomass(anch_sim,tlim = c(550,600))

resource_limitation_2d_cod <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}
edited_cod_prams <- resource_limitation_2d_cod(cod_params,0.00013,10)

cod_sim <- make_limit_cycle_sim(edited_cod_prams)
plotBiomass(cod_sim,tlim=c(550,600))
animate(cod_sim)
################################################################################
# Section 1: Competitiveness ratio + timescale-separation metric
#
# juvenile_adult_rate_ratio = mean juvenile mass-specific intake rate /
# mean adult mass-specific intake rate -- a RATIO, not a difference, so it's
# directly comparable across species regardless of each one's own rate
# units/scale. Close to 1 means juveniles and adults are near-equally
# competitive; far from 1 (either direction) means one stage dominates.
# This is exactly de Roos & Persson's own q comparison, and q_equivalent
# converts it onto their literal q scale for a direct reading against
# Soudijn et al. (2017)'s q=0.5/1/1.5 sweep.
#
# timescale_ratio = resource turnover rate / mean mass-specific consumer
# intake rate. This is the same ratio Soudijn et al.'s body-size sweep
# varies implicitly (consumer rates ~W_A^-0.25, resource rate fixed) --
# computing it directly sidesteps needing to run an actual body-size sweep
# just to see whether cod and anchovy already differ on this axis.
################################################################################

mean_mass_specific_rate <- function(params, idx = TRUE) {
  E <- getEncounter(params)      # encounter rate, g/year
  f <- getFeedingLevel(params)   # feeding level (satiation), 0-1
  per_capita_rate    <- (E * (1 - f))[1, ]
  mass_specific_rate <- per_capita_rate / params@w
  mean(mass_specific_rate[idx])
}

timescale_summary <- function(params, label) {
  w     <- params@w
  w_mat <- params@species_params$w_mat[1]

  msr_all <- mean_mass_specific_rate(params)
  msr_juv <- mean_mass_specific_rate(params, w < w_mat)
  msr_ad  <- mean_mass_specific_rate(params, w >= w_mat)

  # getResourceRate() returns a per-weight-bin VECTOR (the resource is
  # itself a spectrum, not a single pool), not a scalar -- confirmed by
  # actually running this against cod_params, where it's a 361-element
  # vector, ~0 outside the active foraging range. mean() over the nonzero
  # bins only gives one representative rate per species; including the
  # structural zeros would understate cod's rate for no biological reason.
  resource_rate_vec <- getResourceRate(params)
  resource_rate      <- mean(resource_rate_vec[resource_rate_vec > 0])

  # RATIO, not difference: a difference isn't scale-free (its magnitude
  # depends on whatever units the mass-specific rate happens to come out
  # in for that species' calibration -- not comparable across cod and
  # anchovy on its own). The ratio of the two MEANS (not a mean of
  # per-bin ratios -- juvenile and adult have different bin counts, 156
  # vs 44 for cod, so there's no natural pairing to average over) is
  # exactly de Roos & Persson's own q comparison: ratio = 1 means
  # juveniles and adults are equally competitive, which is the regime
  # their paper (and Soudijn et al. 2017's q sweep) puts at the
  # stable/no-cycling end -- ratio far from 1 in either direction is
  # where cycling shows up. q_equivalent inverts Soudijn et al.'s own
  # ratio-to-q relationship (nu_J/nu_A = (2-q)/q, so q = 2/(ratio+1)) to
  # put cod and anchovy directly on their q=0.5/1/1.5 scale.
  juvenile_adult_rate_ratio <- msr_juv / msr_ad

  data.frame(
    species                    = label,
    resource_rate               = resource_rate,
    mean_mass_specific_rate     = msr_all,
    juvenile_rate                = msr_juv,
    adult_rate                   = msr_ad,
    juvenile_adult_rate_ratio    = juvenile_adult_rate_ratio,
    q_equivalent                 = 2 / (juvenile_adult_rate_ratio + 1),
    # >1: resource turns over fast relative to the consumer (large
    # timescale separation, the regime Soudijn et al. tie to large body
    # size and persistent adult-driven cycling). <1: resource is the slow,
    # limiting side.
    timescale_ratio              = resource_rate / msr_all
  )
}

timescale_df <- bind_rows(
  timescale_summary(cod_params, "cod"),
  timescale_summary(anchovy_params, "anchovy")
)

print(timescale_df)
write.csv(timescale_df, file.path("interesting_plots", "day29_timescale_summary.csv"),
          row.names = FALSE)

################################################################################
# Section 2: Method-matched getStability() grid for cod AND anchovy
#
# Day 28 ran this exact grid (resource_decrease_seq x capacity_mult_seq,
# unfished, complex-eigenvalue check for genuine Hopf vs plain collapse) on
# cod_params only. Running the identical grid on anchovy settles Day 27's
# "is the reference point actually oscillating" ambiguity with the same
# analytic tool used for cod, instead of the raw time-series/
# catchable_fraction proxies Days 22-27 relied on for anchovy.
################################################################################

# balance=FALSE: setResource() errors outright if both resource_rate and
# resource_capacity are given while balance resolves to TRUE, since
# balance=TRUE always computes the other one FOR you to preserve the
# current steady state -- exactly the perturbation this sweep introduces.
resource_limitation_2d_cod <- function(params, resource_decrease, capacity_mult) {
  new_rate     <- getResourceRate(params) * resource_decrease
  new_capacity <- getResourceCapacity(params) * capacity_mult
  setResource(params, resource_rate = new_rate, resource_capacity = new_capacity, balance = FALSE)
}

# Two build_fn's, same (resource_decrease, capacity_mult) signature: cod's
# multiplies its OWN existing resource setup (Day 28's convention);
# anchovy's constructs the resource setup directly (Days 22-27's
# convention). check_hopf_grid_point_generic() below doesn't care which.
build_cod     <- function(resource_decrease, capacity_mult) {
  resource_limitation_2d_cod(cod_params, resource_decrease, capacity_mult)
}
build_anchovy <- function(resource_decrease, capacity_mult) {
  make_anchovy_params(resource_decrease = resource_decrease, capacity_mult = capacity_mult)
}

# Same classification as Day 28: spectral_radius alone can't distinguish a
# genuine Hopf bifurcation from a real-eigenvalue (collapse-type)
# instability -- only a complex dominant eigenvalue crossing the unit
# circle is Hopf.
check_hopf_grid_point_generic <- function(build_fn, resource_decrease, capacity_mult, effort = 0) {
  result <- tryCatch({
    p        <- build_fn(resource_decrease, capacity_mult)
    p_steady <- steadyNewton(p, effort = effort, stability = TRUE)
    stab     <- attr(p_steady, "stability")
    dominant <- stab$eigenvalues[1]
    is_complex <- abs(Im(dominant)) > 1e-8

    data.frame(
      spectral_radius  = stab$spectral_radius,
      stable           = stab$stable,
      hopf_period      = if (is.null(stab$hopf_period)) NA_real_ else stab$hopf_period,
      bifurcation_type = if (stab$stable) "stable" else if (is_complex) "hopf" else "non-oscillatory",
      error            = NA_character_
    )
  }, error = function(e) {
    data.frame(spectral_radius = NA_real_, stable = NA, hopf_period = NA_real_,
               bifurcation_type = NA_character_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("resource_decrease=%.4g, capacity_mult=%.4g, effort=%.4g: %s",
                    resource_decrease, capacity_mult, effort, result$error))
  }
  data.frame(resource_decrease = resource_decrease, capacity_mult = capacity_mult,
             effort = effort, result)
}

# Identical to Day 28's grid, so the two species are compared at matched
# resolution.
resource_decrease_seq <- exp(seq(log(0.001), log(1), length.out = 12))
capacity_mult_seq     <- exp(seq(log(1), log(20), length.out = 10))
hopf_grid_params       <- expand.grid(resource_decrease = resource_decrease_seq,
                                      capacity_mult = capacity_mult_seq)

# Each point is expensive (getStability() differentiates numerically, one
# project_n_loop() per active weight-class cell, x2 for centred
# differences, ON TOP of the steadyNewton() solve itself -- see Day 28's
# notes on cost) -- flattened across workers rather than run sequentially,
# same as Day 28.
hopf_grid_generic <- function(build_fn, species_label) {
  df <- bind_rows(future_lapply(seq_len(nrow(hopf_grid_params)), function(i) {
    check_hopf_grid_point_generic(build_fn, hopf_grid_params$resource_decrease[i],
                                  hopf_grid_params$capacity_mult[i])
  }, future.seed = TRUE))
  df$species <- species_label
  df
}

cod_hopf_grid_29     <- hopf_grid_generic(build_cod, "cod")
anchovy_hopf_grid_29 <- hopf_grid_generic(build_anchovy, "anchovy")

species_hopf_grid_df <- bind_rows(cod_hopf_grid_29, anchovy_hopf_grid_29)

print(species_hopf_grid_df %>% filter(is.na(error)) %>% count(species, bifurcation_type))
write.csv(species_hopf_grid_df, file.path("interesting_plots", "day29_species_hopf_grid.csv"),
          row.names = FALSE)

cat(sprintf(
  paste0(
    "Method-matched grid (%d points each, unfished): cod hopf points = %d, ",
    "anchovy hopf points = %d.\n",
    "If anchovy shows hopf points and cod shows none on this IDENTICAL grid, that's direct\n",
    "evidence something about each species' own calibration -- not just the resource_decrease/\n",
    "capacity_mult axis itself -- is the lever. Section 3 tries to isolate what.\n"
  ),
  nrow(hopf_grid_params),
  sum(cod_hopf_grid_29$bifurcation_type == "hopf", na.rm = TRUE),
  sum(anchovy_hopf_grid_29$bifurcation_type == "hopf", na.rm = TRUE)
))

species_hopf_grid_plot <- ggplot(species_hopf_grid_df %>% filter(is.na(error)),
                                 aes(x = resource_decrease, y = capacity_mult, fill = bifurcation_type)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c(stable = "#4C72B0", hopf = "#C44E52", `non-oscillatory` = "#DD8452"),
                    na.value = "grey80") +
  facet_wrap(~species) +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       fill = NULL,
       title = "Same grid, same method (getStability()), cod vs anchovy",
       subtitle = "hopf = spectral_radius >= 1 with a complex dominant eigenvalue") +
  theme_minimal()
species_hopf_grid_plot

save_plot(species_hopf_grid_plot, "day29_species_hopf_grid.png", width = 12)

################################################################################
# Section 3: Isolating the lever -- a body-size x intake-rate-parameter 2x2
# factorial, built inside the anchovy single-species template. Mirrors the
# "species params (anchovy vs default) crossed with gamma alone" 2x2 design
# Day 19/20 already used to isolate a different lever.
#
# Caveat worth keeping in mind reading the results: "cod body"/"cod rates"
# reproduce cod_params' own w_max/w_mat/gamma/alpha/ks values INSIDE the
# anchovy template (same functional-response type, same semichemostat
# resource, same second-order scheme) -- this is not a reconstruction of
# cod_params itself, which may differ in other respects (mortality
# calibration, erepro, R_max, etc. not touched here). It isolates these
# four named traits specifically, nothing more.
################################################################################

cod_sp <- cod_params@species_params
cat(sprintf(
  "cod_params traits pulled into the factorial: gamma=%.4g, alpha=%.4g, ks=%.4g, w_max=%.4g, w_mat=%.4g\n",
  cod_sp$gamma[1], cod_sp$alpha[1], cod_sp$ks[1], cod_sp$w_max[1], cod_sp$w_mat[1]
))

# unname() is load-bearing, not defensive styling: species_params columns
# carry the species name ("Cod") as a names() attribute, which arithmetic
# (w_min <- w_max * ratio, no_w <- round(log(...))) silently propagates.
# Passing a NAMED no_w/w_max/etc. into newSingleSpeciesParams() breaks gear
# setup downstream with "The number of fishing gears must be consistent
# across the catchability and selectivity (dim 1) slots" -- confirmed by
# running this exact construction with and without unname() on real
# cod_params values before writing the rest of this section.
cod_gamma <- unname(cod_sp$gamma[1])
cod_alpha <- unname(cod_sp$alpha[1])
cod_ks    <- unname(cod_sp$ks[1])
cod_w_max <- unname(cod_sp$w_max[1])
cod_w_mat <- unname(cod_sp$w_mat[1])

# Used in Section 5 (selective fishing); pulled here since it's the same
# unname()-the-species_params-column pattern as above, needed for both.
anchovy_w_mat <- unname(anchovy_params@species_params$w_mat[1])

make_hybrid_params <- function(body = c("anchovy", "cod"), rates = c("anchovy", "cod"),
                               resource_decrease = 0.001, capacity_mult = 10,
                               second_order = TRUE, ext_diff = 0.00) {
  body  <- match.arg(body)
  rates <- match.arg(rates)

  w_max <- if (body == "anchovy") 66.5 else cod_w_max
  w_mat <- if (body == "anchovy") 10   else cod_w_mat
  # Same w_min:w_max ratio as the anchovy template (0.0003/66.5), so grid
  # resolution (no_w) scales consistently with body size rather than the
  # discretisation itself changing shape between combos.
  w_min <- w_max * (0.0003 / 66.5)

  gamma <- if (rates == "anchovy") 750 else cod_gamma
  alpha <- if (rates == "anchovy") 0.1 else cod_alpha
  ks    <- if (rates == "anchovy") 0   else cod_ks

  lambda <- 2.05
  a0     <- 100
  kappa  <- a0 * exp(-6.9 * (lambda - 1))
  no_w   <- round(log(w_max / w_min) / 0.1)

  params <- newSingleSpeciesParams(
    species_name = sprintf("body_%s_rates_%s", body, rates),
    w_min = w_min, w_max = w_max, w_mat = w_mat,
    no_w = no_w, lambda = lambda, kappa = kappa,
    alpha = alpha, gamma = gamma, ks = ks
  )

  given_species_params(params)$D_ext <- ext_diff
  params <- setBevertonHolt(params)

  default_capacity <- getResourceCapacity(params)
  r  <- getResourceRate(params) * resource_decrease
  cc <- default_capacity * capacity_mult

  params <- setResource(params, resource_rate = r, resource_capacity = cc,
                        resource_dynamics = "resource_semichemostat",
                        balance = FALSE)

  if (second_order) {
    second_order_w(params) <- c(flux = "centred", bin_average = TRUE)
  }

  params
}

check_hybrid_point <- function(body, rates, resource_decrease, capacity_mult, effort = 0) {
  result <- tryCatch({
    p        <- make_hybrid_params(body = body, rates = rates,
                                   resource_decrease = resource_decrease,
                                   capacity_mult = capacity_mult)
    p_steady <- steadyNewton(p, effort = effort, stability = TRUE)
    stab     <- attr(p_steady, "stability")
    dominant <- stab$eigenvalues[1]
    is_complex <- abs(Im(dominant)) > 1e-8

    data.frame(
      spectral_radius  = stab$spectral_radius,
      stable           = stab$stable,
      hopf_period      = if (is.null(stab$hopf_period)) NA_real_ else stab$hopf_period,
      bifurcation_type = if (stab$stable) "stable" else if (is_complex) "hopf" else "non-oscillatory",
      error            = NA_character_
    )
  }, error = function(e) {
    data.frame(spectral_radius = NA_real_, stable = NA, hopf_period = NA_real_,
               bifurcation_type = NA_character_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("body=%s, rates=%s, resource_decrease=%.4g, capacity_mult=%.4g: %s",
                    body, rates, resource_decrease, capacity_mult, result$error))
  }
  data.frame(body = body, rates = rates, resource_decrease = resource_decrease,
             capacity_mult = capacity_mult, effort = effort, result)
}

factorial_combos <- expand.grid(body = c("anchovy", "cod"), rates = c("anchovy", "cod"),
                                stringsAsFactors = FALSE)

# Reference operating point first (cheap, 4 points) -- Day 24-27's own
# baseline (resource_decrease = 0.001, capacity_mult = 10).
factorial_reference_df <- bind_rows(lapply(seq_len(nrow(factorial_combos)), function(i) {
  check_hybrid_point(factorial_combos$body[i], factorial_combos$rates[i],
                     resource_decrease = 0.001, capacity_mult = 10)
}))
print(factorial_reference_df)
write.csv(factorial_reference_df, file.path("interesting_plots", "day29_factorial_reference.csv"),
          row.names = FALSE)

# Widen each of the 4 combos across a grid, so a combo that looks stable at
# the single reference point isn't mistaken for stable everywhere. Coarser
# than Section 2's grid (5x4=20 rather than 12x10=120 points per combo,
# 80 points total) purely for cost -- this section already runs 4 combos'
# worth of getStability() solves on top of Section 2's 240.
factorial_resource_decrease_seq <- exp(seq(log(0.001), log(1), length.out = 5))
factorial_capacity_mult_seq     <- exp(seq(log(1), log(20), length.out = 4))

factorial_grid_params <- expand.grid(
  body = c("anchovy", "cod"), rates = c("anchovy", "cod"),
  resource_decrease = factorial_resource_decrease_seq,
  capacity_mult = factorial_capacity_mult_seq,
  stringsAsFactors = FALSE
)

factorial_grid_df <- bind_rows(future_lapply(seq_len(nrow(factorial_grid_params)), function(i) {
  row <- factorial_grid_params[i, ]
  check_hybrid_point(row$body, row$rates, row$resource_decrease, row$capacity_mult)
}, future.seed = TRUE))

write.csv(factorial_grid_df, file.path("interesting_plots", "day29_factorial_grid.csv"),
          row.names = FALSE)

factorial_summary <- factorial_grid_df %>%
  filter(is.na(error)) %>%
  count(body, rates, bifurcation_type)
print(factorial_summary)
write.csv(factorial_summary, file.path("interesting_plots", "day29_factorial_summary.csv"),
          row.names = FALSE)

factorial_grid_plot <- ggplot(factorial_grid_df %>% filter(is.na(error)),
                              aes(x = resource_decrease, y = capacity_mult, fill = bifurcation_type)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_manual(values = c(stable = "#4C72B0", hopf = "#C44E52", `non-oscillatory` = "#DD8452"),
                    na.value = "grey80") +
  facet_grid(rows = vars(rates), cols = vars(body),
            labeller = labeller(rates = function(x) paste0("rates: ", x),
                                body  = function(x) paste0("body: ", x))) +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       fill = NULL,
       title = "Isolating the lever: body size x intake-rate parameters",
       subtitle = "Anchovy template throughout; 'cod' body/rates swap in cod_params' own w_max/w_mat or gamma/alpha/ks") +
  theme_minimal()
factorial_grid_plot

save_plot(factorial_grid_plot, "day29_factorial_grid.png", width = 10, height = 9)

################################################################################
# Section 4: Idea #1 -- relax the reproduction cap
#
# Before sweeping anything: R_max for BOTH cod_params and the anchovy
# template already comes back Inf (species_params$R_max), i.e.
# reproduction_level = R_dd/R_max is already ~0 for both -- there is no
# Beverton-Holt cap currently active to relax. Recruitment is already
# fully density-independent (R_dd = R_di). So "relax the cap" can't be
# tested in the direction it was originally proposed in -- there's
# nowhere further to relax to.
#
# What CAN be tested is the only direction actually available: does
# IMPOSING compensation (raising reproduction_level toward 1, the
# opposite of relaxing) ever change the verdict? Beverton-Holt
# compensation is a saturating, not overshooting, nonlinearity (unlike
# Ricker), so classical stock-recruit theory predicts it should only ever
# add damping, never induce cycling on its own -- this checks that
# expectation rather than assuming it.
################################################################################

set_reproduction_level <- function(params, reproduction_level) {
  setBevertonHolt(params, reproduction_level = reproduction_level)
}

cat(sprintf(
  "Current R_max: cod = %s, anchovy = %s (Inf means no Beverton-Holt cap is active -- 'relax the cap' has nowhere further to go).\n",
  format(cod_params@species_params$R_max[1]), format(anchovy_params@species_params$R_max[1])
))

check_reproduction_level_point <- function(build_fn, reproduction_level, resource_decrease = 0.001,
                                           capacity_mult = 10, effort = 0) {
  result <- tryCatch({
    p        <- set_reproduction_level(build_fn(resource_decrease, capacity_mult), reproduction_level)
    p_steady <- steadyNewton(p, effort = effort, stability = TRUE)
    stab     <- attr(p_steady, "stability")
    dominant <- stab$eigenvalues[1]
    is_complex <- abs(Im(dominant)) > 1e-8
    data.frame(
      spectral_radius  = stab$spectral_radius,
      stable           = stab$stable,
      bifurcation_type = if (stab$stable) "stable" else if (is_complex) "hopf" else "non-oscillatory",
      error            = NA_character_
    )
  }, error = function(e) {
    data.frame(spectral_radius = NA_real_, stable = NA, bifurcation_type = NA_character_,
               error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("reproduction_level=%.4g: %s", reproduction_level, result$error))
  }
  data.frame(reproduction_level = reproduction_level, result)
}

# 0.99 rather than closer to 1: setBevertonHolt() itself caps
# reproduction_level requests below what R_dd needs at 0.99, so this
# already reaches its own practical ceiling.
reproduction_level_seq <- c(1e-4, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)

reproduction_level_df <- bind_rows(
  bind_rows(future_lapply(reproduction_level_seq, check_reproduction_level_point,
                          build_fn = build_cod, future.seed = TRUE)) %>%
    mutate(species = "cod"),
  bind_rows(future_lapply(reproduction_level_seq, check_reproduction_level_point,
                          build_fn = build_anchovy, future.seed = TRUE)) %>%
    mutate(species = "anchovy")
)

print(reproduction_level_df)
write.csv(reproduction_level_df, file.path("interesting_plots", "day29_reproduction_level.csv"),
          row.names = FALSE)

reproduction_level_plot <- ggplot(reproduction_level_df %>% filter(is.na(error)),
                                  aes(x = reproduction_level, y = spectral_radius, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(x = "reproduction_level (R_dd / R_max -- 0 = uncapped, 1 = fully saturated)",
       y = "Spectral radius", color = NULL,
       title = "Does adding Beverton-Holt recruitment compensation ever destabilise?",
       subtitle = "Reference point (resource_decrease=0.001, capacity_mult=10); dashed = stability boundary") +
  theme_minimal()
reproduction_level_plot
save_plot(reproduction_level_plot, "day29_repro_level.png")

################################################################################
# Section 5: Idea #2 -- selective/size-targeted fishing
#
# Day 28 only tested a flat effort=0.5 uniform across all sizes. Andersen
# & Pedersen (2010) -- already in Further_Reading.md #16 -- found that
# fishing NARROWLY targeting a size range is specifically the most
# destabilising fishing pattern, more so than raising uniform effort.
# Tested directly: three selectivity patterns (uniform/full-catch,
# knife-edge at w_mat, and a narrow Gaussian dome centred on w_mat) swept
# across fishing effort, at the reference resource point, for both
# species.
################################################################################

# Weight-native custom selectivity -- mizer resolves sel_func by NAME via
# do.call() inside calc_selectivity(), confirmed directly (built a
# one-off dome-selectivity gear on cod_params and checked getStability()
# ran on it) before writing the rest of this section around it.
narrow_dome_weight <- function(w, w_center, w_width, ...) {
  exp(-((w - w_center) / w_width)^2)
}

set_selectivity <- function(params, pattern, w_center) {
  species_name <- params@species_params$species[1]
  gear_df <- switch(pattern,
    uniform = data.frame(species = species_name, gear = "uniform_gear", sel_func = "knife_edge",
                         knife_edge_size = min(params@w), catchability = 1),
    knife_edge_adult = data.frame(species = species_name, gear = "knife_gear", sel_func = "knife_edge",
                                  knife_edge_size = w_center, catchability = 1),
    dome_narrow = data.frame(species = species_name, gear = "dome_gear", sel_func = "narrow_dome_weight",
                             w_center = w_center, w_width = w_center * 0.1, catchability = 1)
  )
  gear_params(params) <- gear_df
  params
}

check_fishing_point <- function(build_fn, pattern, w_center, effort, resource_decrease = 0.001,
                                capacity_mult = 10) {
  result <- tryCatch({
    p        <- set_selectivity(build_fn(resource_decrease, capacity_mult), pattern, w_center)
    p_steady <- steadyNewton(p, effort = effort, stability = TRUE)
    stab     <- attr(p_steady, "stability")
    dominant <- stab$eigenvalues[1]
    is_complex <- abs(Im(dominant)) > 1e-8
    data.frame(
      spectral_radius  = stab$spectral_radius,
      stable           = stab$stable,
      bifurcation_type = if (stab$stable) "stable" else if (is_complex) "hopf" else "non-oscillatory",
      error            = NA_character_
    )
  }, error = function(e) {
    data.frame(spectral_radius = NA_real_, stable = NA, bifurcation_type = NA_character_,
               error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("pattern=%s, effort=%.4g: %s", pattern, effort, result$error))
  }
  data.frame(pattern = pattern, effort = effort, result)
}

fishing_effort_seq <- c(0, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
fishing_patterns    <- c("uniform", "knife_edge_adult", "dome_narrow")
fishing_grid_params <- expand.grid(pattern = fishing_patterns, effort = fishing_effort_seq,
                                   stringsAsFactors = FALSE)

# Plain lapply(), not future_lapply(): narrow_dome_weight is resolved by
# STRING name inside mizer's calc_selectivity(), not a direct symbol
# reference in check_fishing_point()'s body, so future's static
# globals-detection (which scans for symbol references) has no way to
# know a worker process needs it exported. Running this sequentially
# sidesteps that risk entirely rather than guessing whether the export
# would actually happen.
fishing_sweep_generic <- function(build_fn, w_center, species_label) {
  df <- bind_rows(lapply(seq_len(nrow(fishing_grid_params)), function(i) {
    check_fishing_point(build_fn, fishing_grid_params$pattern[i], w_center, fishing_grid_params$effort[i])
  }))
  df$species <- species_label
  df
}

cod_fishing_df     <- fishing_sweep_generic(build_cod, cod_w_mat, "cod")
anchovy_fishing_df <- fishing_sweep_generic(build_anchovy, anchovy_w_mat, "anchovy")

fishing_sweep_df <- bind_rows(cod_fishing_df, anchovy_fishing_df)

print(fishing_sweep_df %>% filter(is.na(error)) %>% count(species, pattern, bifurcation_type))
write.csv(fishing_sweep_df, file.path("interesting_plots", "day29_fishing_sweep.csv"), row.names = FALSE)

fishing_sweep_plot <- ggplot(fishing_sweep_df %>% filter(is.na(error)),
                             aes(x = effort, y = spectral_radius, color = pattern)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  facet_wrap(~species) +
  labs(x = "Fishing effort", y = "Spectral radius", color = "Selectivity pattern",
       title = "Does narrow size-selective fishing destabilise where uniform fishing doesn't?",
       subtitle = "Andersen & Pedersen (2010): narrow-band selectivity is the most destabilising fishing pattern") +
  theme_minimal()
fishing_sweep_plot
save_plot(fishing_sweep_plot, "day29_fishing_sweep.png", width = 10)

################################################################################
# Section 6: Idea #4 -- push capacity_mult further than Section 2's grid
#
# Section 2 only went up to capacity_mult=20. Days 22-24's original,
# never-confirmed-as-Hopf anchovy destabilisation was reported somewhere
# past ~capacity_mult=7-10, so 20x should already have been comfortable
# headroom if that story were right -- but it's cheap to push further and
# rule out "the grid just stopped short" directly, at the reference
# resource_decrease=0.001 where any paradox-of-enrichment effect would be
# expected to show up first.
################################################################################

capacity_mult_extended_seq <- exp(seq(log(1), log(500), length.out = 20))

extended_capacity_generic <- function(build_fn, species_label) {
  df <- bind_rows(future_lapply(capacity_mult_extended_seq, function(cm) {
    check_hopf_grid_point_generic(build_fn, resource_decrease = 0.001, capacity_mult = cm)
  }, future.seed = TRUE))
  df$species <- species_label
  df
}

cod_extended_df     <- extended_capacity_generic(build_cod, "cod")
anchovy_extended_df <- extended_capacity_generic(build_anchovy, "anchovy")
extended_capacity_df <- bind_rows(cod_extended_df, anchovy_extended_df)

print(extended_capacity_df %>% filter(is.na(error)) %>% count(species, bifurcation_type))
write.csv(extended_capacity_df, file.path("interesting_plots", "day29_extended_capacity.csv"),
          row.names = FALSE)

extended_capacity_plot <- ggplot(extended_capacity_df %>% filter(is.na(error)),
                                 aes(x = capacity_mult, y = spectral_radius, color = species)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(x = "capacity_mult (multiplier on resource carrying capacity, up to 500x)",
       y = "Spectral radius", color = NULL,
       title = "Does pushing capacity_mult well past Section 2's 20x cap ever destabilise?",
       subtitle = "resource_decrease fixed at 0.001 (the reference point)") +
  theme_minimal()
extended_capacity_plot
save_plot(extended_capacity_plot, "day29_extended_capacity.png")

################################################################################
# Section 8: Cross-validating getStability()'s linear verdict against
# actual perturbed simulation
#
# getStability() only proves LOCAL stability: it linearises at the steady
# state and checks whether infinitesimal perturbations grow or decay. That
# can't rule out a fixed point that's linearly stable but still reachable
# into a sustained limit cycle by a large-enough perturbation (a
# subcritical Hopf, or any other form of bistability) -- every "stable"
# verdict in Sections 2-6 rests on the linear check alone, never on
# actually running the dynamics forward from a disturbed state. This
# section does that instead, reusing make_limit_cycle_sim() (resource
# dropped to 10% of capacity, then a slice of the consumer spectrum
# further cut by 1000x, then run forward) at points where getStability()
# has already given a verdict, so the two methods are directly comparable.
################################################################################

# Day 27 established this noise floor directly: a genuinely settled run
# comes back rel_amplitude ~1e-12, a real oscillation O(0.01)-O(1) -- a
# 1e-6 cutoff sits well inside that 12-order-of-magnitude gap.
classify_perturbed_run <- function(sim, window_frac = 0.4) {
  bm    <- rowSums(getBiomass(sim))
  times <- as.numeric(names(bm))
  late  <- bm[times >= max(times) * (1 - window_frac)]
  rel_amplitude <- (max(late) - min(late)) / mean(late)
  data.frame(rel_amplitude = rel_amplitude,
             perturbed_verdict = if (rel_amplitude > 1e-6) "oscillating" else "settled",
             error = NA_character_)
}

perturbed_stability_point <- function(params, t_total = 600, effort = 0) {
  result <- tryCatch({
    sim <- make_limit_cycle_sim(params, t_total = t_total, effort = effort)
    classify_perturbed_run(sim)
  }, error = function(e) {
    data.frame(rel_amplitude = NA_real_, perturbed_verdict = NA_character_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("perturbed_stability_point(): %s", result$error))
  }
  result
}

linear_stability_point <- function(params, effort = 0) {
  result <- tryCatch({
    p_steady <- steadyNewton(params, effort = effort, stability = TRUE)
    stab     <- attr(p_steady, "stability")
    dominant <- stab$eigenvalues[1]
    is_complex <- abs(Im(dominant)) > 1e-8
    data.frame(spectral_radius = stab$spectral_radius, stable = stab$stable,
               linear_verdict = if (stab$stable) "stable" else if (is_complex) "hopf" else "non-oscillatory",
               error = NA_character_)
  }, error = function(e) {
    data.frame(spectral_radius = NA_real_, stable = NA, linear_verdict = NA_character_,
               error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("linear_stability_point(): %s", result$error))
  }
  result
}

dual_stability_check <- function(params, label, effort = 0, t_total = 600) {
  lin  <- linear_stability_point(params, effort = effort)
  pert <- perturbed_stability_point(params, t_total = t_total, effort = effort)
  data.frame(species            = label,
             linear_verdict     = lin$linear_verdict, spectral_radius = lin$spectral_radius,
             perturbed_verdict  = pert$perturbed_verdict, rel_amplitude = pert$rel_amplitude,
             agree              = identical(lin$linear_verdict == "stable", pert$perturbed_verdict == "settled"))
}

# Each species' own current reference calibration, plus the most
# enrichment-heavy point from Section 6's extended sweep -- the regime
# most likely to hide a perturbation-only instability if one exists
# anywhere in this project's search so far.
dual_check_df <- bind_rows(
  dual_stability_check(cod_params, "cod (reference)"),
  dual_stability_check(anchovy_params, "anchovy (reference)"),
  dual_stability_check(build_cod(resource_decrease = 0.001, capacity_mult = max(capacity_mult_extended_seq)),
                       "cod (capacity_mult=500x)"),
  dual_stability_check(build_anchovy(resource_decrease = 0.001, capacity_mult = max(capacity_mult_extended_seq)),
                       "anchovy (capacity_mult=500x)")
)

print(dual_check_df)
write.csv(dual_check_df, file.path("interesting_plots", "day29_dual_stability_check.csv"),
          row.names = FALSE)

cat(sprintf(
  "\nLinear (getStability()) vs perturbed-simulation verdicts agree at %d/%d checked points.\n",
  sum(dual_check_df$agree, na.rm = TRUE), nrow(dual_check_df)
))
if (any(!dual_check_df$agree, na.rm = TRUE)) {
  cat("Disagreement found -- exactly the case getStability() alone can't catch:\n")
  print(dual_check_df %>% filter(!agree))
}

################################################################################
# Section 9: The continuous picture -- perturbed-amplitude heatmap
#
# Section 8 found a real disagreement at anchovy's reference point: stable
# by getStability() (spectral_radius=0.35), but rel_amplitude=0.84 --
# clearly oscillating -- under an actual perturbed simulation. Sections 2
# and 6 only ever reported a three-way categorical call (stable/hopf/
# non-oscillatory) from the LINEAR check across the resource_decrease x
# capacity_mult grid, which is exactly the check that just got caught
# missing something. This reruns that same grid with the perturbed
# simulation instead, keeping rel_amplitude itself as a continuous fill
# color (log-scaled -- Day 27's own numbers span ~1e-12 for settled up to
# O(1) for a real limit cycle, 12+ orders of magnitude, so a binary
# stable/not call would hide most of the structure) rather than collapsing
# it back down to a three-way category.
################################################################################

perturbed_grid_point <- function(build_fn, resource_decrease, capacity_mult, effort = 0, t_total = 600) {
  result <- tryCatch({
    p <- build_fn(resource_decrease, capacity_mult)
    perturbed_stability_point(p, t_total = t_total, effort = effort)
  }, error = function(e) {
    data.frame(rel_amplitude = NA_real_, perturbed_verdict = NA_character_, error = conditionMessage(e))
  })
  if (!is.na(result$error)) {
    warning(sprintf("resource_decrease=%.4g, capacity_mult=%.4g: %s",
                    resource_decrease, capacity_mult, result$error))
  }
  data.frame(resource_decrease = resource_decrease, capacity_mult = capacity_mult, result)
}

# Same grid as Section 2 (hopf_grid_params, resource_decrease_seq x
# capacity_mult_seq), so this heatmap and Section 2's categorical one are
# directly comparable point-for-point. Each point is a full t_max=600
# nonlinear simulation rather than a single Jacobian solve -- costlier
# than Section 2's grid, but Section 8's 3-point timing check ran in
# ~10s total, so the full 240-point grid (120 x 2 species) is tractable
# in the same run, especially split across workers.
perturbed_grid_generic <- function(build_fn, species_label) {
  df <- bind_rows(future_lapply(seq_len(nrow(hopf_grid_params)), function(i) {
    perturbed_grid_point(build_fn, hopf_grid_params$resource_decrease[i], hopf_grid_params$capacity_mult[i])
  }, future.seed = TRUE))
  df$species <- species_label
  df
}

cod_perturbed_grid_df     <- perturbed_grid_generic(build_cod, "cod")
anchovy_perturbed_grid_df <- perturbed_grid_generic(build_anchovy, "anchovy")
perturbed_grid_df <- bind_rows(cod_perturbed_grid_df, anchovy_perturbed_grid_df)

print(perturbed_grid_df %>% filter(is.na(error)) %>% group_by(species) %>%
        summarise(min_amplitude = min(rel_amplitude), max_amplitude = max(rel_amplitude),
                  n_oscillating = sum(perturbed_verdict == "oscillating"), .groups = "drop"))
write.csv(perturbed_grid_df, file.path("interesting_plots", "day29_perturbed_amplitude_grid.csv"),
          row.names = FALSE)

# Floored, not the raw column, so a perfectly-settled point (rel_amplitude
# could in principle land on exactly 0) doesn't break the log10 fill
# scale -- the CSV above keeps the unfloored value.
perturbed_grid_plot_df <- perturbed_grid_df %>%
  filter(is.na(error)) %>%
  mutate(rel_amplitude_floored = pmax(rel_amplitude, 1e-15))

# Sequential, single-hue, perceptually uniform, colorblind-safe (viridis)
# for this continuous magnitude field -- not the categorical
# stable/hopf/non-oscillatory palette Sections 2/6 use, since this is a
# magnitude now, not an identity.
perturbed_amplitude_plot <- ggplot(perturbed_grid_plot_df,
                                   aes(x = resource_decrease, y = capacity_mult, fill = rel_amplitude_floored)) +
  geom_tile() +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(trans = "log10", labels = scales::label_scientific(),
                       name = "rel_amplitude\n(perturbed sim)") +
  facet_wrap(~species) +
  labs(x = "resource_decrease (multiplier on resource renewal rate)",
       y = "capacity_mult (multiplier on resource carrying capacity)",
       title = "Perturbed-simulation oscillation amplitude across resource rate and carrying capacity",
       subtitle = "Continuous rel_amplitude (log scale), not a binary stable/hopf call -- ~1e-12 = settled, O(0.01)-O(1) = a real limit cycle") +
  theme_minimal()
perturbed_amplitude_plot
save_plot(perturbed_amplitude_plot, "day29_perturbed_amplitude.png", width = 12)

################################################################################
# Section 10: Summary -- programmatic readout, not asserted conclusions.
# Read factorial_summary above once it's actually run: if only the
# body=cod rows show hopf points regardless of rates, body size is the
# lever; if only rates=cod rows do regardless of body, gamma/alpha/ks are
# the lever; if only body=cod & rates=cod (i.e. only when both match cod)
# shows hopf points, it's an interaction, not either alone; if NEITHER
# swap reproduces cod's flatness/anchovy's cycling, this q/body-size
# mechanism isn't the lever at all.
################################################################################

cat("\n===== Day 29 summary =====\n")
cat(sprintf(
  "Timescale ratio (resource_rate / mean mass-specific consumer rate): cod = %.4g, anchovy = %.4g\n",
  timescale_df$timescale_ratio[timescale_df$species == "cod"],
  timescale_df$timescale_ratio[timescale_df$species == "anchovy"]
))
cat(sprintf(
  "juvenile_adult_rate_ratio (q_equivalent): cod = %.4g (%.4g), anchovy = %.4g (%.4g)\n",
  timescale_df$juvenile_adult_rate_ratio[timescale_df$species == "cod"],
  timescale_df$q_equivalent[timescale_df$species == "cod"],
  timescale_df$juvenile_adult_rate_ratio[timescale_df$species == "anchovy"],
  timescale_df$q_equivalent[timescale_df$species == "anchovy"]
))
cat(sprintf(
  "Section 2 (method-matched, native cod_params/anchovy calibration): cod hopf points = %d/%d, anchovy hopf points = %d/%d.\n",
  sum(cod_hopf_grid_29$bifurcation_type == "hopf", na.rm = TRUE), nrow(hopf_grid_params),
  sum(anchovy_hopf_grid_29$bifurcation_type == "hopf", na.rm = TRUE), nrow(hopf_grid_params)
))
cat("Section 3 factorial (body x rates) hopf-point counts, out of the sub-grid tested:\n")
print(factorial_summary %>% filter(bifurcation_type == "hopf"))
cat(sprintf(
  "Section 4 (reproduction_level sweep, idea #1): hopf points found = %d/%d.\n",
  sum(reproduction_level_df$bifurcation_type == "hopf", na.rm = TRUE), nrow(reproduction_level_df)
))
cat(sprintf(
  "Section 5 (selective fishing, idea #2): hopf points found = %d/%d.\n",
  sum(fishing_sweep_df$bifurcation_type == "hopf", na.rm = TRUE), nrow(fishing_sweep_df)
))
cat(sprintf(
  "Section 6 (extended capacity_mult to 500x, idea #4): hopf points found = %d/%d.\n",
  sum(extended_capacity_df$bifurcation_type == "hopf", na.rm = TRUE), nrow(extended_capacity_df)
))
cat(sprintf(
  "Section 8 (linear vs perturbed-simulation cross-check): agreement at %d/%d points checked.\n",
  sum(dual_check_df$agree, na.rm = TRUE), nrow(dual_check_df)
))
cat(sprintf(
  "Section 9 (perturbed-amplitude heatmap): cod oscillating points = %d/%d, anchovy oscillating points = %d/%d.\n",
  sum(cod_perturbed_grid_df$perturbed_verdict == "oscillating", na.rm = TRUE), nrow(cod_perturbed_grid_df),
  sum(anchovy_perturbed_grid_df$perturbed_verdict == "oscillating", na.rm = TRUE), nrow(anchovy_perturbed_grid_df)
))
