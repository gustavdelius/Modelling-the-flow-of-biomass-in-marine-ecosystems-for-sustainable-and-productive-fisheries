#| message: false
#| warning: false
library(mizer)

make_params <- function(lambda = 2.05, resource_decrease = 0.001) {
  kappa  <- 100 * exp(-6.9 * (lambda - 1))
  params <- newSingleSpeciesParams(species_name = "Anchovy",
                                   w_min = 0.0003, w_max = 66.5, w_mat = 10,
                                   lambda = lambda, kappa = kappa, alpha = 0.1, gamma = 750)
  r <- getResourceRate(params) * resource_decrease
  setResource(params, resource_rate = r,
              resource_dynamics = "resource_semichemostat", balance = TRUE)
}

anchovy_params <- make_params(resource_decrease = 1)  # native, undepleted resource

anchovy_osc <- make_params(resource_decrease = 0.001)
anchovy_osc@initial_n[]    <- 0.001 * anchovy_osc@w^(-1.8)
anchovy_osc@initial_n_pp[] <- anchovy_osc@cc_pp

sim     <- project(anchovy_osc, t_max = 300, dt = 0.1, t_save = 0.2,
                   effort = 0, progress_bar = FALSE, method = "predictor-corrector")
sim_600 <- project(sim, t_max = 300, dt = 0.1, t_save = 0.2,
                   effort = 0, progress_bar = FALSE, method = "predictor-corrector")

library(reshape2); library(plotly); library(dplyr)

nice_biomass_plot <- function(sim, t) {
  abm  <- melt(getBiomass(sim)); abmr <- melt(getBiomass(sim, min_w = 0.01, max_w = 0.4))
  abmr$sp <- "small Anchovy"
  pbm <- sim@n_pp %*% (sim@params@w_full * sim@params@dw_full)
  pbm <- melt(pbm); names(pbm)[names(pbm) == "Var1"] <- "time"; pbm$Var2 <- NULL; pbm$sp <- "Plankton"
  bm  <- rbind(pbm, abm, abmr)
  plot_ly(bm) %>% filter(time >= t) %>%
    add_lines(x = ~time, y = ~value, color = ~sp) %>%
    layout(yaxis = list(type = "log", title = "biomass (g/m^3)", range = c(-7, 2)),
           xaxis = list(title = "time (year)"))
}

nice_animation <- function(sim, t) {
  nf <- melt(sim@n); n_ppf <- melt(sim@n_pp); n_ppf$sp <- "Plankton"; nf <- rbind(nf, n_ppf)
  plot_ly(nf) %>% filter(w > 10^-5, time >= t) %>%
    mutate(b = value * w^2) %>%
    add_lines(x = ~w, y = ~b, color = ~sp, frame = ~time, line = list(simplify = FALSE)) %>%
    layout(xaxis = list(type = "log", title = "body mass (g)"),
           yaxis = list(type = "log", title = "biomass (g/m^3)", range = c(-8, 0)))
}
nice_animation(sim_600, 550)
nice_biomass_plot(sim_600, 550)