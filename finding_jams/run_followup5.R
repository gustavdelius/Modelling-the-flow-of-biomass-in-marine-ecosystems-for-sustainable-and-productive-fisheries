lines <- readLines("finding_jams/23_experiments.R")
needed <- c(
  lines[1:12],     # libraries
  lines[14:49],    # make_second_order_params_kr
  lines[247:252],  # dir.create / save_plot
  lines[609:804]   # Follow-up 5 block itself
)
eval(parse(text = needed), envir = .GlobalEnv)
cat("FOLLOWUP5 DONE\n")
