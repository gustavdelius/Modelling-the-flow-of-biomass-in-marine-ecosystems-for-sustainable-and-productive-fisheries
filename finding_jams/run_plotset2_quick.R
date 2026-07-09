lines <- readLines("finding_jams/23_experiments.R")
eval(parse(text = c(
  lines[1:12],     # libraries
  lines[247:252],  # dir.create / save_plot
  lines[635:677],  # make_second_order_params_balanced (fixed version)
  lines[679:722]   # run_bifurcation_sweep
)), envir = .GlobalEnv)

# Cut-down version, not written back to 23_experiments.R: 12 points (vs 20)
# and t_run = 300 (vs 600) -- t_run = 300 already validated as sufficient to
# settle in today's Follow-ups 3/4, so this trades resolution, not
# reliability, to fit a ~20-minute budget instead of ~90.
rd_seq_quick <- exp(seq(log(0.0001), log(0.5), length.out = 12))
ext_diff_values_bif <- seq(0.001, 0.5, length.out = 4)

bif_balanced_df <- bind_rows(lapply(ext_diff_values_bif, function(ed) {
  run_bifurcation_sweep(rd_seq_quick, "resource_decrease",
                        fixed_params = list(resource_level = 1, ext_diff = ed),
                        params_fn = make_second_order_params_balanced,
                        t_run = 300) %>%
    mutate(ext_diff = ed)
})) %>%
  mutate(ext_diff_label = factor(sprintf("ext_diff = %.3f", ext_diff),
                                 levels = sprintf("ext_diff = %.3f", sort(unique(ext_diff)))))

saveRDS(bif_balanced_df, file.path("interesting_plots", "bifurcation_balanced_df_quick.rds"))

bif_balanced_plot <- ggplot(bif_balanced_df,
                            aes(x = value, y = biomass, color = direction, linetype = branch)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  scale_x_log10() +
  facet_wrap(~ext_diff_label, nrow = 1) +
  labs(x = "resource_decrease", y = "Biomass",
       title = "Bifurcation: resource_decrease swept at 4 diffusion strengths (quick pass)",
       subtitle = "balance = TRUE, resource_level = 1 -- 12 points, t_run = 300") +
  theme_minimal()
bif_balanced_plot

save_plot(bif_balanced_plot, "Bifurcation - balanced setup across diffusion strengths (quick).png",
         width = 16, height = 5)

cat("PLOTSET2_QUICK DONE\n")
