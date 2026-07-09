lines <- readLines("finding_jams/23_experiments.R")
needed <- c(
  lines[1:12],     # libraries
  lines[247:252],  # dir.create / save_plot
  lines[740],      # rd_seq_bif
  lines[635:677],  # make_second_order_params_balanced (fixed version)
  lines[679:722],  # run_bifurcation_sweep
  lines[767:804]   # Section 2: balance=TRUE diffusion sweep
)
eval(parse(text = needed), envir = .GlobalEnv)
cat("PLOTSET2 DONE\n")
