library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstatix)

ratio_fixed_pms <- list()
# Or ratio_pms <- list()
variable_platform <- 0
# Other parameters used are set to default

for (ratio in seq(0.2, 0.8, by = 0.2)) {
  print(paste("Processing ratio:", ratio))
  weight_wall <- ratio # DC weight
  weight_place <- 1 - weight_wall # PC weight

  # Define static platform location
  platform_x0 <- cos(-pi / 4) * pool_diameter / 4
  platform_y0 <- sin(-pi / 4) * pool_diameter / 4
  platform_x_mag <- abs(platform_x0)
  platform_y_mag <- abs(platform_y0)

  # Starting locations of the modeled animal (4 different ones)
  strad <- pool_diameter / 2 * 0.85
  starting_xs <- strad * c(cos(pi / 6), cos(pi / 3), cos(7 * pi / 6), cos(4 * pi / 3))
  starting_ys <- strad * c(sin(pi / 6), sin(pi / 3), sin(7 * pi / 6), sin(4 * pi / 3))

  Nruns <- 30
  PMs_comb <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  track_x_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  track_y_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  name <- paste0("final_clean_pm_fixed_ratio_", ratio)
  ratio_fixed_pms[[name]] <- data.frame()

  for (reps in 1:Nruns) {
    # ==================== Combined Model =======================
    print(reps)
    weights_pc <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult # c
    weights_dc <- matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult # c

    # Generate place and distance cells
    DC <- rep(0, N_dc)
    PC_x <- rep(0, N_pc) # 1xN_pc matrix containing the x = 0 coordinate for each place cell# c
    PC_y <- rep(0, N_pc) # 1xN_pc matrix containing the y = 0 coordinate for each place cell# c

    for (i in 1:N_dc) {
      # For each place cell:
      DC[i] <- runif(1) * (pool_diameter / 2) # Random positions of distance cells
    }
    for (i in 1:N_pc) { # c
      # For each place cell:
      PC_x[i] <- (runif(1) - 0.5) * pool_diameter # Random positions of place cells
      PC_y[i] <- (runif(1) - 0.5) * pool_diameter
      while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter / 2)^2)) {
        # Checks for out of bounds
        PC_x[i] <- (runif(1) - 0.5) * pool_diameter
        PC_y[i] <- (runif(1) - 0.5) * pool_diameter
      }
    }

    # Shared trial schedule for fair model comparison in this run
    start_idx_sched <- matrix(0, nrow = Ndays, ncol = Ntrials)
    platform_x_sign_sched <- matrix(sign(platform_x0), nrow = Ndays, ncol = Ntrials)
    platform_y_sign_sched <- matrix(sign(platform_y0), nrow = Ndays, ncol = Ntrials)
    for (day in 1:Ndays) {
      start_idx_sched[day, ] <- sample(1:4)
      if (variable_platform == 1) {
        platform_x_sign_sched[day, ] <- sample(c(1, -1), Ntrials, replace = TRUE)
        platform_y_sign_sched[day, ] <- sample(c(1, -1), Ntrials, replace = TRUE)
      }
    }

    for (day in 1:Ndays) {
      for (trial in 1:Ntrials) {
        platform_x_trial <- platform_x_sign_sched[day, trial] * platform_x_mag
        platform_y_trial <- platform_y_sign_sched[day, trial] * platform_y_mag
        target_quad <- get_quadrant_index(platform_x_trial, platform_y_trial)
        opposite_quad <- get_opposite_quadrant(target_quad)

        idx <- start_idx_sched[day, trial] # take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]
        # run trial
        modresults <- run_trial_comb(weights_pc, weights_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall, weight_wall, weight_place)
        # update weights and trak trajectories for future steps
        weights_pc <- modresults[[1]] # c
        weights_dc <- modresults[[2]] # c
        track_x <- modresults[[3]]
        track_y <- modresults[[4]]
        vel_x <- modresults[[5]]
        vel_y <- modresults[[6]]
        track_idx <- (reps - 1) * Ndays * Ntrials + (day - 1) * Ntrials + trial
        track_x_sum[[track_idx]] <- track_x
        track_y_sum[[track_idx]] <- track_y

        PMs_comb[1, day, trial] <- modresults[[10]] # latency
        PMs_comb[2, day, trial] <- modresults[[7]] # dist
        PMs_comb[3, day, trial] <- modresults[[9]][target_quad] * 100 # target quadrant
        PMs_comb[4, day, trial] <- modresults[[9]][opposite_quad] * 100 # opposite quadrant
        PMs_comb[5, day, trial] <- modresults[[8]] * 100 # wall zone
      }
    }

    # ==================== Cleaning Data ====================
    dimnames(PMs_comb) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:5), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )
    clean_pm_comb <- as.data.frame.table(PMs_comb, responseName = "value")
    clean_pm_comb$Model <- "Combined"
    clean_pm_comb$Reps <- reps
    # clean_pm$Reps <- reps
    ratio_fixed_pms[[name]] <- rbind(
      ratio_fixed_pms[[name]],
      clean_pm_comb
    )
  }
  ratio_fixed_pms[[name]]$wall_ratio <- ratio
  ratio_fixed_pms[[paste0(name, "_wide")]] <- pivot_wider(
    ratio_fixed_pms[[name]],
    names_from = Parameter,
    values_from = value
  )
}

# Same for variable platform: ratio_pms_wide, ratio_combined_wide ...
ratio_fixed_pms_wide <- ratio_fixed_pms[grep("_wide$", names(ratio_fixed_pms))]
ratio_fixed_combined_wide <- bind_rows(ratio_fixed_pms_wide)
ratio_fixed_combined_wide$wall_ratio <- as.factor(ratio_fixed_combined_wide$wall_ratio)

# The following lines assume that `ratio_combined_wide` is already generated using
# previous steps
ratio_fixed_combined_wide$Platform <- "fixed"
ratio_combined_wide$Platform <- "variable"
ratio_combined_wide_full <- rbind(ratio_fixed_combined_wide, ratio_combined_wide)
ratio_combined_wide_full_day5 <- ratio_combined_wide_full[ratio_combined_wide_full$Day == "Day_5", ] %>%
  group_by(Day, Platform, wall_ratio) %>%
  summarise(
    mean_latency = mean(Latency, na.rm = TRUE),
    se_latency = sd(Latency, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Statistical tests and data visualization
# 1. Boxplots of Target Quadrant percentages across different wall ratios and platforms
ggplot(ratio_combined_wide_full, aes(x = Day, y = Target_quad, fill = wall_ratio)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~Platform) +
  labs(
    x = "Training Day",
    y = "Target Quadrant (%)",
    fill = "Weight"
  ) +
  scale_fill_manual(values = c(
    "0.2" = "#BFDFD2",
    "0.4" = "#51999F",
    "0.6" = "#7BC0CD",
    "0.8" = "#DBCB92"
  ))

# Conduct pairwise t-tests for Day 5, comparing latencies across different wall ratios within each platform
posthoc <- ratio_combined_wide_full[ratio_combined_wide_full$Day == "Day_5", ] %>%
  group_by(Platform) %>%
  pairwise_t_test(
    Latency ~ wall_ratio,
    p.adjust.method = "bonferroni"
  ) %>%
  add_y_position(fun = "max")

# 2. Final Day Latency plot
ggplot(ratio_combined_wide_full_day5, aes(x = wall_ratio, y = mean_latency, fill = wall_ratio)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_latency - se_latency, ymax = mean_latency + se_latency),
    position = position_dodge(width = 0.8),
    width = 0.2,
    alpha = 0.6
  ) +
  facet_wrap(~Platform) +
  stat_pvalue_manual(
    posthoc,
    label = "p.adj.signif",
    y.position = 38,
    step.increase = 0.15,
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  labs(
    y = "Latency (s)",
    x = "weight",
    fill = "Weight"
  ) +
  scale_fill_manual(values = c(
    "0.2" = "#BFDFD2",
    "0.4" = "#51999F",
    "0.6" = "#7BC0CD",
    "0.8" = "#DBCB92"
  )) +
  theme_classic() +
  theme(legend.position = "none")
