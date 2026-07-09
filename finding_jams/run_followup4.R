lines <- readLines("finding_jams/23_experiments.R")
needed <- c(
  lines[1:12],     # libraries
  lines[14:49],    # make_second_order_params_kr
  lines[153],      # MIN_VIABLE_BIOMASS
  lines[247:252],  # dir.create / save_plot
  lines[278:338],  # build_path_capacity_rows / build_path_rd_rows / run_along_path
  lines[395:401],  # classify_phase
  lines[538],      # rd_grid_fine
  lines[566:607]   # Follow-up 4 block itself
)
eval(parse(text = needed), envir = .GlobalEnv)
cat("FOLLOWUP4 DONE\n")
