library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstatix)

# ========= Alpha tuning (learning rate [0.005..0.02]) =========
# High Alpha = Fast, volatile learning; Low Alpha = slow, stable learning
alpha_pms <- list() # Or alpha_fixed_pms <- list()
variable_platform <- 1

for (a in c(0.005, 0.01, 0.015, 0.02)) {
  print(paste("Processing alpha:", a))
  alpha <- a
  platform_x0 <- cos(-pi / 4) * pool_diameter / 4
  platform_y0 <- sin(-pi / 4) * pool_diameter / 4
  platform_x_mag <- abs(platform_x0)
  platform_y_mag <- abs(platform_y0)

  # Starting locations of the modeled animal (4 different ones)
  strad <- pool_diameter / 2 * 0.85
  starting_xs <- strad * c(cos(pi / 6), cos(pi / 3), cos(7 * pi / 6), cos(4 * pi / 3))
  starting_ys <- strad * c(sin(pi / 6), sin(pi / 3), sin(7 * pi / 6), sin(4 * pi / 3))

  Nruns <- 30
  PMs_pc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  PMs_dc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  PMs_comb <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  track_x_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  track_y_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  name <- paste0("final_clean_pm_alpha_", a)
  alpha_pms[[name]] <- data.frame() # Or alpha_fixed_pms[[name]] <- data.frame()

  for (reps in 1:Nruns) {
    # ==================== Place Cell Model =======================
    # Generate initial weights for each run
    weights <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult
    print(reps)
    # Generate place cells for each run
    PC_x <- rep(0, N_pc) # 1xN_pc matrix containing the x = 0 coordinate for each place cell
    PC_y <- rep(0, N_pc) # 1xN_pc matrix containing the y = 0 coordinate for each place cell
    for (i in 1:N_pc) {
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

        modresults <- run_trial_pc(weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall)
        # run trial
        weights <- modresults[[1]]

        PMs_pc[1, day, trial] <- modresults[[9]] # latency
        PMs_pc[2, day, trial] <- modresults[[6]] # dist
        PMs_pc[3, day, trial] <- modresults[[8]][target_quad] * 100 # target quadrant
        PMs_pc[4, day, trial] <- modresults[[8]][opposite_quad] * 100 # opposite quadrant
        PMs_pc[5, day, trial] <- modresults[[7]] * 100 # wall zone
      }
    }
    # ==================== Distance Cell Model =======================
    weights <- matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult
    # Generate distance cells for each run
    DC <- rep(0, N_dc) # 1xN_dc matrix containing the each place cell
    for (i in 1:N_dc) {
      # For each place cell:
      DC[i] <- runif(1) * (pool_diameter / 2) # Random positions of distance cells
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

        modresults <- run_trial_dc(weights, Wmult, sigma_dc, sigma_ac, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall)
        # run trial
        weights <- modresults[[1]]

        PMs_dc[1, day, trial] <- modresults[[9]] # latency
        PMs_dc[2, day, trial] <- modresults[[6]] # dist
        PMs_dc[3, day, trial] <- modresults[[8]][target_quad] * 100 # target quadrant
        PMs_dc[4, day, trial] <- modresults[[8]][opposite_quad] * 100 # opposite quadrant
        PMs_dc[5, day, trial] <- modresults[[7]] * 100 # wall zone
        # record performance measures
      }
    }

    # ==================== Combined Model =======================
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

    for (day in 1:Ndays) {
      for (trial in 1:Ntrials) {
        platform_x_trial <- platform_x_sign_sched[day, trial] * platform_x_mag
        platform_y_trial <- platform_y_sign_sched[day, trial] * platform_y_mag
        target_quad <- get_quadrant_index(platform_x_trial, platform_y_trial)
        opposite_quad <- get_opposite_quadrant(target_quad)

        idx <- start_idx_sched[day, trial] # take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]

        modresults <- run_trial_comb(weights_pc, weights_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall, weight_wall, weight_place)
        # run trial
        weights_pc <- modresults[[1]] # c
        weights_dc <- modresults[[2]] # c
        track_x <- modresults[[3]]
        track_y <- modresults[[4]]
        vel_x <- modresults[[5]]
        vel_y <- modresults[[6]]
        track_idx <- (reps - 1) * Ndays * Ntrials + (day - 1) * Ntrials + trial
        track_x_sum[[track_idx]] <- track_x
        track_y_sum[[track_idx]] <- track_y
        #        weights <- wres

        PMs_comb[1, day, trial] <- modresults[[10]] # latency
        PMs_comb[2, day, trial] <- modresults[[7]] # dist
        PMs_comb[3, day, trial] <- modresults[[9]][target_quad] * 100 # target quadrant
        PMs_comb[4, day, trial] <- modresults[[9]][opposite_quad] * 100 # opposite quadrant
        PMs_comb[5, day, trial] <- modresults[[8]] * 100 # wall zone
      }
    }

    # ==================== Cleaning Data ====================
    dimnames(PMs_dc) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )
    dimnames(PMs_pc) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )
    dimnames(PMs_comb) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )

    clean_pm_dc <- as.data.frame.table(PMs_dc, responseName = "value")
    clean_pm_dc$Model <- "Distance_cell"
    clean_pm_dc$Reps <- reps
    clean_pm_pc <- as.data.frame.table(PMs_pc, responseName = "value")
    clean_pm_pc$Model <- "Place_cell"
    clean_pm_pc$Reps <- reps
    clean_pm_comb <- as.data.frame.table(PMs_comb, responseName = "value")
    clean_pm_comb$Model <- "Combined"
    clean_pm_comb$Reps <- reps
    # clean_pm$Reps <- reps
    alpha_pms[[name]] <- rbind(
      alpha_pms[[name]], clean_pm_dc, clean_pm_pc,
      clean_pm_comb
    )
  }

  alpha_pms[[name]]$Alpha <- a
  alpha_pms[[paste0(name, "_wide")]] <- pivot_wider(
    alpha_pms[[name]],
    names_from = Parameter,
    values_from = value
  )
}

# Or alpha_combined_fixed_wide
alpha_pms_wide <- alpha_pms[grep("_wide$", names(alpha_pms))]
alpha_combined_wide <- bind_rows(alpha_pms_wide)
alpha_combined_wide$Alpha <- as.factor(alpha_combined_wide$Alpha)

# ========= Gamma tuning (Discount factor [0.75..0.95]) =========
# High Gamma = Long-sighted; Low Gamma = Short-sighted
gamma_pms <- list() # Or gamma_fixed_pms <- list()
variable_platform <- 1

for (g in seq(0.75, 0.95, by = 0.05)) {
  print(paste("Processing gamma:", g))
  gamma <- g
  platform_x0 <- cos(-pi / 4) * pool_diameter / 4
  platform_y0 <- sin(-pi / 4) * pool_diameter / 4
  platform_x_mag <- abs(platform_x0)
  platform_y_mag <- abs(platform_y0)

  # Starting locations of the modeled animal (4 different ones)
  strad <- pool_diameter / 2 * 0.85
  starting_xs <- strad * c(cos(pi / 6), cos(pi / 3), cos(7 * pi / 6), cos(4 * pi / 3))
  starting_ys <- strad * c(sin(pi / 6), sin(pi / 3), sin(7 * pi / 6), sin(4 * pi / 3))

  Nruns <- 30
  PMs_pc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  PMs_dc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  PMs_comb <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  track_x_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  track_y_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  name <- paste0("final_clean_pm_gamma_", g)
  gamma_pms[[name]] <- data.frame() # Or gamma_fixed_pms[[name]] <- data.frame()

  for (reps in 1:Nruns) {
    # ==================== Place Cell Model =======================
    # Generate initial weights for each run
    weights <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult
    print(reps)
    # Generate place cells for each run
    PC_x <- rep(0, N_pc) # 1xN_pc matrix containing the x = 0 coordinate for each place cell
    PC_y <- rep(0, N_pc) # 1xN_pc matrix containing the y = 0 coordinate for each place cell
    for (i in 1:N_pc) {
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

        modresults <- run_trial_pc(weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall)
        # run trial
        weights <- modresults[[1]]

        PMs_pc[1, day, trial] <- modresults[[9]] # latency
        PMs_pc[2, day, trial] <- modresults[[6]] # dist
        PMs_pc[3, day, trial] <- modresults[[8]][target_quad] * 100 # target quadrant
        PMs_pc[4, day, trial] <- modresults[[8]][opposite_quad] * 100 # opposite quadrant
        PMs_pc[5, day, trial] <- modresults[[7]] * 100 # wall zone
      }
    }
    # ==================== Distance Cell Model =======================
    weights <- matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult
    # Generate distance cells for each run
    DC <- rep(0, N_dc) # 1xN_dc matrix containing the each place cell
    for (i in 1:N_dc) {
      # For each place cell:
      DC[i] <- runif(1) * (pool_diameter / 2) # Random positions of distance cells
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

        modresults <- run_trial_dc(weights, Wmult, sigma_dc, sigma_ac, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall)
        # run trial
        weights <- modresults[[1]]

        PMs_dc[1, day, trial] <- modresults[[9]] # latency
        PMs_dc[2, day, trial] <- modresults[[6]] # dist
        PMs_dc[3, day, trial] <- modresults[[8]][target_quad] * 100 # target quadrant
        PMs_dc[4, day, trial] <- modresults[[8]][opposite_quad] * 100 # opposite quadrant
        PMs_dc[5, day, trial] <- modresults[[7]] * 100 # wall zone
        # record performance measures
      }
    }

    # ==================== Combined Model =======================
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

    for (day in 1:Ndays) {
      for (trial in 1:Ntrials) {
        platform_x_trial <- platform_x_sign_sched[day, trial] * platform_x_mag
        platform_y_trial <- platform_y_sign_sched[day, trial] * platform_y_mag
        target_quad <- get_quadrant_index(platform_x_trial, platform_y_trial)
        opposite_quad <- get_opposite_quadrant(target_quad)

        idx <- start_idx_sched[day, trial] # take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]

        modresults <- run_trial_comb(weights_pc, weights_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall, weight_wall, weight_place)
        # run trial
        weights_pc <- modresults[[1]] # c
        weights_dc <- modresults[[2]] # c
        track_x <- modresults[[3]]
        track_y <- modresults[[4]]
        vel_x <- modresults[[5]]
        vel_y <- modresults[[6]]
        track_idx <- (reps - 1) * Ndays * Ntrials + (day - 1) * Ntrials + trial
        track_x_sum[[track_idx]] <- track_x
        track_y_sum[[track_idx]] <- track_y
        #        weights <- wres

        PMs_comb[1, day, trial] <- modresults[[10]] # latency
        PMs_comb[2, day, trial] <- modresults[[7]] # dist
        PMs_comb[3, day, trial] <- modresults[[9]][target_quad] * 100 # target quadrant
        PMs_comb[4, day, trial] <- modresults[[9]][opposite_quad] * 100 # opposite quadrant
        PMs_comb[5, day, trial] <- modresults[[8]] * 100 # wall zone
      }
    }

    # ==================== Cleaning Data ====================
    dimnames(PMs_dc) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )
    dimnames(PMs_pc) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )
    dimnames(PMs_comb) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )

    clean_pm_dc <- as.data.frame.table(PMs_dc, responseName = "value")
    clean_pm_dc$Model <- "Distance_cell"
    clean_pm_dc$Reps <- reps
    clean_pm_pc <- as.data.frame.table(PMs_pc, responseName = "value")
    clean_pm_pc$Model <- "Place_cell"
    clean_pm_pc$Reps <- reps
    clean_pm_comb <- as.data.frame.table(PMs_comb, responseName = "value")
    clean_pm_comb$Model <- "Combined"
    clean_pm_comb$Reps <- reps
    # clean_pm$Reps <- reps
    gamma_pms[[name]] <- rbind(
      gamma_pms[[name]], clean_pm_dc, clean_pm_pc,
      clean_pm_comb
    )
  }

  gamma_pms[[name]]$Gamma <- g
  gamma_pms[[paste0(name, "_wide")]] <- pivot_wider(
    gamma_pms[[name]],
    names_from = Parameter,
    values_from = value
  )
}

# Or gamma_combined_fixed_wide
gamma_pms_wide <- gamma_pms[grep("_wide$", names(gamma_pms))]
gamma_combined_wide <- bind_rows(gamma_pms_wide)
gamma_combined_wide$Gamma <- as.factor(gamma_combined_wide$Gamma)

# ========= Beta tuning (Exploration-exploitation factor [0.5,12]) =========
# High Beta = Exploitation; Low Beta = Exploration
beta_pms <- list() # Or beta_fixed_pms <- list()
variable_platform <- 1

for (b in c(0.5, 1, 6, 9, 12)) {
  print(paste("Processing beta:", b))
  beta <- b
  platform_x0 <- cos(-pi / 4) * pool_diameter / 4
  platform_y0 <- sin(-pi / 4) * pool_diameter / 4
  platform_x_mag <- abs(platform_x0)
  platform_y_mag <- abs(platform_y0)

  # Starting locations of the modeled animal (4 different ones)
  strad <- pool_diameter / 2 * 0.85
  starting_xs <- strad * c(cos(pi / 6), cos(pi / 3), cos(7 * pi / 6), cos(4 * pi / 3))
  starting_ys <- strad * c(sin(pi / 6), sin(pi / 3), sin(7 * pi / 6), sin(4 * pi / 3))

  Nruns <- 30
  PMs_pc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  PMs_dc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  PMs_comb <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
  track_x_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  track_y_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
  name <- paste0("final_clean_pm_beta_", b)
  beta_pms[[name]] <- data.frame() # Or beta_fixed_pms[[name]] <- data.frame()

  for (reps in 1:Nruns) {
    # ==================== Place Cell Model =======================
    # Generate initial weights for each run
    weights <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult
    print(reps)
    # Generate place cells for each run
    PC_x <- rep(0, N_pc) # 1xN_pc matrix containing the x = 0 coordinate for each place cell
    PC_y <- rep(0, N_pc) # 1xN_pc matrix containing the y = 0 coordinate for each place cell
    for (i in 1:N_pc) {
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

        modresults <- run_trial_pc(weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall)
        # run trial
        weights <- modresults[[1]]

        PMs_pc[1, day, trial] <- modresults[[9]] # latency
        PMs_pc[2, day, trial] <- modresults[[6]] # dist
        PMs_pc[3, day, trial] <- modresults[[8]][target_quad] * 100 # target quadrant
        PMs_pc[4, day, trial] <- modresults[[8]][opposite_quad] * 100 # opposite quadrant
        PMs_pc[5, day, trial] <- modresults[[7]] * 100 # wall zone
      }
    }
    # ==================== Distance Cell Model =======================
    weights <- matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult
    # Generate distance cells for each run
    DC <- rep(0, N_dc) # 1xN_dc matrix containing the each place cell
    for (i in 1:N_dc) {
      # For each place cell:
      DC[i] <- runif(1) * (pool_diameter / 2) # Random positions of distance cells
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

        modresults <- run_trial_dc(weights, Wmult, sigma_dc, sigma_ac, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall)
        # run trial
        weights <- modresults[[1]]

        PMs_dc[1, day, trial] <- modresults[[9]] # latency
        PMs_dc[2, day, trial] <- modresults[[6]] # dist
        PMs_dc[3, day, trial] <- modresults[[8]][target_quad] * 100 # target quadrant
        PMs_dc[4, day, trial] <- modresults[[8]][opposite_quad] * 100 # opposite quadrant
        PMs_dc[5, day, trial] <- modresults[[7]] * 100 # wall zone
        # record performance measures
      }
    }

    # ==================== Combined Model =======================
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

    for (day in 1:Ndays) {
      for (trial in 1:Ntrials) {
        platform_x_trial <- platform_x_sign_sched[day, trial] * platform_x_mag
        platform_y_trial <- platform_y_sign_sched[day, trial] * platform_y_mag
        target_quad <- get_quadrant_index(platform_x_trial, platform_y_trial)
        opposite_quad <- get_opposite_quadrant(target_quad)

        idx <- start_idx_sched[day, trial] # take each location
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]

        modresults <- run_trial_comb(weights_pc, weights_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x_trial, platform_y_trial, starting_x, starting_y, speed, hitwall, weight_wall, weight_place)
        # run trial
        weights_pc <- modresults[[1]] # c
        weights_dc <- modresults[[2]] # c
        track_x <- modresults[[3]]
        track_y <- modresults[[4]]
        vel_x <- modresults[[5]]
        vel_y <- modresults[[6]]
        track_idx <- (reps - 1) * Ndays * Ntrials + (day - 1) * Ntrials + trial
        track_x_sum[[track_idx]] <- track_x
        track_y_sum[[track_idx]] <- track_y
        #        weights <- wres

        PMs_comb[1, day, trial] <- modresults[[10]] # latency
        PMs_comb[2, day, trial] <- modresults[[7]] # dist
        PMs_comb[3, day, trial] <- modresults[[9]][target_quad] * 100 # target quadrant
        PMs_comb[4, day, trial] <- modresults[[9]][opposite_quad] * 100 # opposite quadrant
        PMs_comb[5, day, trial] <- modresults[[8]] * 100 # wall zone
      }
    }

    # ==================== Cleaning Data ====================
    dimnames(PMs_dc) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )
    dimnames(PMs_pc) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )
    dimnames(PMs_comb) <- list(
      Parameter = paste0(c("Latency", "Dist", "Target_quad", "Oppo_quad", "Wall_zone")), # Rows
      Day       = paste0("Day_", 1:Ndays), # Columns
      Trial     = paste0("Trial_", 1:4) # Tables (3rd dimension)
    )

    clean_pm_dc <- as.data.frame.table(PMs_dc, responseName = "value")
    clean_pm_dc$Model <- "Distance_cell"
    clean_pm_dc$Reps <- reps
    clean_pm_pc <- as.data.frame.table(PMs_pc, responseName = "value")
    clean_pm_pc$Model <- "Place_cell"
    clean_pm_pc$Reps <- reps
    clean_pm_comb <- as.data.frame.table(PMs_comb, responseName = "value")
    clean_pm_comb$Model <- "Combined"
    clean_pm_comb$Reps <- reps
    # clean_pm$Reps <- reps
    beta_pms[[name]] <- rbind(
      beta_pms[[name]], clean_pm_dc, clean_pm_pc,
      clean_pm_comb
    )
  }

  beta_pms[[name]]$Beta <- b
  beta_pms[[paste0(name, "_wide")]] <- pivot_wider(
    beta_pms[[name]],
    names_from = Parameter,
    values_from = value
  )
}

# or beta_combined_fixed_wide
beta_pms_wide <- beta_pms[grep("_wide$", names(beta_pms))]
beta_combined_wide <- bind_rows(beta_pms_wide)
beta_combined_wide$Beta <- as.factor(beta_combined_wide$Beta)

# ========= Data visualization =========
# ========= 1. Learning curves (quantified by latency) =========
# Same for alpha, gamma tunning
# Or beta_combined_fixed_wide_summary_Latency
beta_combined_wide_summary_Latency <- beta_combined_wide %>%
  group_by(Day, Model, Beta) %>%
  summarise(
    mean_latency = mean(Latency, na.rm = TRUE),
    # Calculate Standard Error (SE) for error bars
    se_latency = sd(Latency, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

beta_combined_wide_summary_Latency$Day <- factor(beta_combined_wide_summary_Latency$Day,
  levels = c("Day_1", "Day_2", "Day_3", "Day_4", "Day_5")
)
beta_combined_wide_summary_Latency$Beta <- as.factor(beta_combined_wide_summary_Latency$Beta)

ggplot(beta_combined_wide_summary_Latency, aes(x = Day, y = mean_latency, group = Model, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  # Add error bars (Mean +/- Standard Error)
  geom_errorbar(aes(ymin = mean_latency - se_latency, ymax = mean_latency + se_latency), width = 0.2) +
  # Create separate panels for each Beta learning rate
  facet_wrap(~Beta) +
  # Labels and aesthetics
  labs(
    x = "Training Day",
    y = "Average Latency (seconds)",
    color = "Model Strategy"
  ) +
  scale_color_manual(values = c(
    "Combined" = "#BFDFD2",
    "Distance_cell" = "#51999F",
    "Place_cell" = "#7BC0CD"
  )) +
  theme_classic()

# ========== 2. Learning gain (Day 1 - Day 5) plot =========
# Same for alpha, gamma tunning
fc_raw <- beta_combined_wide %>%
  group_by(Model, Beta, Reps, Day) %>%
  summarise(
    Latency = mean(Latency),
    .groups = "drop"
  )
fc_raw <- fc_raw %>%
  group_by(Model, Beta, Reps) %>%
  summarise(
    day1 = Latency[Day == "Day_1"],
    day5 = Latency[Day == "Day_5"],
    fold_change = day1 - day5,
    .groups = "drop"
  )
fc_summary <- fc_raw %>%
  group_by(Model, Beta) %>%
  summarise(
    mean_fc = mean(fold_change),
    se_fc = sd(fold_change) / sqrt(n()),
    .groups = "drop"
  )

posthoc <- fc_raw %>%
  group_by(Model) %>%
  pairwise_t_test(
    fold_change ~ Beta,
    p.adjust.method = "bonferroni"
  ) %>%
  add_y_position(fun = "max")

ggplot(fc_summary, aes(x = Beta, y = mean_fc, fill = Model)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_fc - se_fc, ymax = mean_fc + se_fc),
    width = 0.2,
    alpha = 0.6
  ) +
  facet_wrap(~Model) +
  stat_pvalue_manual(
    posthoc,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  ) +
  labs(
    x = "Beta",
    y = "Learning Gain (Day 1 - Day 5)"
  ) +
  scale_fill_manual(values = c(
    "Combined" = "#BFDFD2",
    "Distance_cell" = "#51999F",
    "Place_cell" = "#7BC0CD"
  )) +
  theme_classic() +
  theme(legend.position = "none")
