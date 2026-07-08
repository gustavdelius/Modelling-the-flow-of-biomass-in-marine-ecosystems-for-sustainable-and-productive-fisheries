suppressPackageStartupMessages({
  library(mizer); library(patchwork); library(reshape2); library(plotly)
  library(dplyr); library(mizerExperimental); library(tidyverse); library(glue)
  library(ggplot2); library(colorRamps); library(future); library(scales)
})

# pull in just the function defs (lines 14-49 and 149-332) without running the
# full top-level sweep code, by sourcing only those line ranges
lines <- readLines("23_experiments.R")
defs  <- c(lines[14:49], lines[149:332])
eval(parse(text = defs), envir = .GlobalEnv)

t0 <- Sys.time()
small_snake <- run_snake_grid(rd_grid = c(0.001, 0.01), cc_grid = c(3, 10), t_run = 20)
cat("run_snake_grid OK, took", as.numeric(Sys.time() - t0), "s\n")
print(small_snake)

# now hit the MIN_VIABLE_BIOMASS line directly, as the real script does
res <- tryCatch({
  small_snake %>%
    mutate(collapsed = mean_bm < MIN_VIABLE_BIOMASS,
           phase = case_when(collapsed ~ "Collapsed", rel_amplitude > 1e-2 ~ "Oscillating", TRUE ~ "Fixed point"))
}, error = function(e) { cat("ERROR at MIN_VIABLE_BIOMASS mutate:", conditionMessage(e), "\n"); NULL })
