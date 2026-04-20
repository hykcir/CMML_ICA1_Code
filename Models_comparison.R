library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rstatix)

# Robust action probability calculation, aim at avoiding NaN and Inf issues
safe_action_probs <- function(values, beta) {
  scores <- as.numeric(values)
  scores[!is.finite(scores)] <- 0
  scores <- pmax(scores, 0)
  scores <- scores^beta

  total <- sum(scores)
  if (!is.finite(total) || total <= 0) {
    return(rep(1 / length(scores), length(scores)))
  }

  probs <- scores / total
  probs[!is.finite(probs)] <- 0
  total2 <- sum(probs)
  if (!is.finite(total2) || total2 <= 0) {
    return(rep(1 / length(scores), length(scores)))
  }

  probs / total2
}

# Maps an (x, y) coordinate to a standard quadrant (1-4).
# 1: Top-Right, 2: Top-Left, 3: Bottom-Left, 4: Bottom-Right
get_quadrant_index <- function(x, y) {
  if (x > 0 && y > 0) {
    return(1)
  }
  if (x < 0 && y > 0) {
    return(2)
  }
  if (x < 0 && y < 0) {
    return(3)
  }
  return(4)
}

# Returns the diametrically opposite quadrant
get_opposite_quadrant <- function(q) {
  c(3, 4, 1, 2)[q]
}

# ================== Run Trials Functions ===================
run_trial_dc <- function(weights0, Wmult, sigma_dc, sigma_ac, DC, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun) {
  # FIXED PARAMETERS OF THE EXPERIMENT

  pool_diameter <- 1.4 # Maze diameter in metres (m)
  platform_radius <- 0.06 # Platform radius

  N_dc <- 10 # Population of distance cells
  N_ac <- 36 # Population of action cells
  which <- 0

  dist <- 0
  wall_zone <- 0
  quadrants <- c(0, 0, 0, 0) # Percentage spent on each quadrant

  weights <- weights0 # Initialize modifiable weights

  el_tr <- matrix(rep(0, N_dc * N_ac), nrow = N_dc) # Initialize eligibility traces matrix

  # Initialize trajectories
  track_x <- starting_x # Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0

  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2) {
    weights <- weights * (1 - noise) + matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult * noise

    # current distance to the wall
    dist.to.wall <- pool_diameter / 2 - sqrt(track_x[length(track_x)]^2 + track_y[length(track_y)]^2)

    # Calculate DC activation
    DC_activation <- rep(0, N_dc)
    for (i in 1:N_dc) {
      DC_activation[i] <- exp(-(dist.to.wall - DC[i])^2 / (2 * sigma_dc^2))
    }

    # Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1) {
      prevQ <- AC_activation[which] # Displays the Q value before movement
    }

    AC_activation <- DC_activation %*% weights

    # Make an action
    ACsel <- safe_action_probs(AC_activation, beta)
    ASrand <- runif(1)
    which <- 1
    ASsum <- ACsel[1]
    while (which < N_ac && ASsum < ASrand) {
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }

    # Eligibility traces
    el_tr <- el_tr * etdecay

    for (j in 1:N_ac) {
      itmp <- min(abs(j - which), N_ac - abs(j - which))
      actgaus <- exp(-(itmp * itmp) / (2 * sigma_ac * sigma_ac))
      el_tr[, j] <- el_tr[, j] + actgaus * AC_activation[j] * t(t(DC_activation))
    }

    # moving direction
    # dir_x =
    if (track_x[length(track_x)] >= 0) {
      central.angle <- atan(track_y[length(track_y)] / (track_x[length(track_x)]))
    } else {
      central.angle <- pi + atan(track_y[length(track_y)] / (track_x[length(track_x)]))
    }
    moving.dir <- pi + central.angle + which / N_ac * 2 * pi


    vel_x <- c(vel_x, (vel_x[length(vel_x)] + ac_const * cos(moving.dir)) * Vdecay)
    vel_y <- c(vel_y, (vel_y[length(vel_y)] + ac_const * sin(moving.dir)) * Vdecay)
    # velocity per time step (not second)
    track_x <- c(track_x, track_x[length(track_x)] + vel_x[length(vel_x)])
    track_y <- c(track_y, track_y[length(track_y)] + vel_y[length(vel_y)])

    # Check if not out of bounds, reset location & speed if so
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter / 2)^2) {
      ratio <- (track_x[length(track_x)]^2 + track_y[length(track_y)]^2) / ((pool_diameter / 2)^2)
      track_x[length(track_x)] <- track_x[length(track_x)] / sqrt(ratio)
      track_y[length(track_y)] <- track_y[length(track_y)] / sqrt(ratio)
      vel_x[length(vel_x)] <- track_x[length(track_x)] - track_x[length(track_x) - 1]
      vel_y[length(vel_y)] <- track_y[length(track_y)] - track_y[length(track_y) - 1]
    }


    if (length(track_x) > 2) {
      if ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 < platform_radius^2) {
        rew <- 10
      } # found platform - reward
      else if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (0.99 * pool_diameter / 2)^2) {
        rew <- -wall_pun
      } # hit wall - punishment
      else {
        rew <- 0
      } # didn't find - no reward

      currQ <- AC_activation[which]
      tderr <- rew + discf * currQ - prevQ # temporal difference error
      weights <- pmax(weights + lrate * tderr * el_tr, 0)
    }

    laststep <- sqrt((track_x[length(track_x)] - track_x[length(track_x) - 1])^2 + (track_y[length(track_y)] - track_y[length(track_y) - 1])^2)
    dist <- dist + laststep

    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > 0.8 * (pool_diameter / 2)^2) {
      wall_zone <- wall_zone + 1
    } else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0) {
      quadrants[1] <- quadrants[1] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0) {
      quadrants[2] <- quadrants[2] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0) {
      quadrants[3] <- quadrants[3] + 1
    } else {
      quadrants[4] <- quadrants[4] + 1
    }

    if (length(track_x) > 100) # evaluate latency only after 100+ steps to be accurate
      {
        speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2))) # speed in meters/time step
        latency <- (length(track_x) - 1) * speed_ts / speed # convert to seconds
        if (latency > 60) # if more than a minute, stop
          {
            break
          }
      }
  }

  latency <- length(track_x) - 1 # latency in time steps
  wall_zone <- wall_zone / latency
  quadrants <- quadrants / latency
  speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2))) # speed in meters/time step

  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency))
}

run_trial_pc <- function(weights0, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun) {
  # FIXED PARAMETERS OF THE EXPERIMENT

  pool_diameter <- 1.4 # Maze diameter in metres (m)
  platform_radius <- 0.06 # Platform radius

  N_pc <- 211 # Population of place cells
  N_ac <- 36 # Population of action cells
  which <- 0

  dist <- 0
  wall_zone <- 0
  quadrants <- c(0, 0, 0, 0) # Percentage spent on each quadrant

  weights <- weights0 # Initialize modifiable weights

  el_tr <- matrix(rep(0, N_pc * N_ac), nrow = N_pc) # Initialize eligibility traces matrix

  # Initialize trajectories
  track_x <- starting_x # Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0

  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2) {
    weights <- weights * (1 - noise) + matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult * noise

    # Calculate PC activation
    PC_activation <- rep(0, N_pc)
    for (i in 1:N_pc) {
      PC_activation[i] <- exp(-((track_x[length(track_x)] - PC_x[i])^2 + (track_y[length(track_y)] - PC_y[i])^2) / (2 * sigma_pc^2))
    }

    # Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1) {
      prevQ <- AC_activation[which] # Displays the Q value before movement
    }

    AC_activation <- PC_activation %*% weights

    # Make an action
    ACsel <- safe_action_probs(AC_activation, beta)
    ASrand <- runif(1)
    which <- 1
    ASsum <- ACsel[1]
    while (which < N_ac && ASsum < ASrand) {
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }

    # Eligibility traces
    el_tr <- el_tr * etdecay

    for (j in 1:N_ac) {
      itmp <- min(abs(j - which), N_ac - abs(j - which))
      actgaus <- exp(-(itmp * itmp) / (2 * sigma_ac * sigma_ac))
      el_tr[, j] <- el_tr[, j] + actgaus * AC_activation[j] * t(t(PC_activation))
    }

    vel_x <- c(vel_x, (vel_x[length(vel_x)] + ac_const * cos(which / N_ac * 2 * pi)) * Vdecay)
    vel_y <- c(vel_y, (vel_y[length(vel_y)] + ac_const * sin(which / N_ac * 2 * pi)) * Vdecay)
    # velocity per time step (not second)
    track_x <- c(track_x, track_x[length(track_x)] + vel_x[length(vel_x)])
    track_y <- c(track_y, track_y[length(track_y)] + vel_y[length(vel_y)])

    # Check if not out of bounds, reset location & speed if so
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter / 2)^2) {
      ratio <- (track_x[length(track_x)]^2 + track_y[length(track_y)]^2) / ((pool_diameter / 2)^2)
      track_x[length(track_x)] <- track_x[length(track_x)] / sqrt(ratio)
      track_y[length(track_y)] <- track_y[length(track_y)] / sqrt(ratio)
      vel_x[length(vel_x)] <- track_x[length(track_x)] - track_x[length(track_x) - 1]
      vel_y[length(vel_y)] <- track_y[length(track_y)] - track_y[length(track_y) - 1]
    }


    if (length(track_x) > 2) {
      if ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 < platform_radius^2) {
        rew <- 10
      } # found platform - reward
      else if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (0.99 * pool_diameter / 2)^2) {
        rew <- -wall_pun
      } # hit wall - punishment
      else {
        rew <- 0
      } # didn't find - no reward

      currQ <- AC_activation[which]
      tderr <- rew + discf * currQ - prevQ # temporal difference error
      weights <- pmax(weights + lrate * tderr * el_tr, 0)
    }

    laststep <- sqrt((track_x[length(track_x)] - track_x[length(track_x) - 1])^2 + (track_y[length(track_y)] - track_y[length(track_y) - 1])^2)
    dist <- dist + laststep

    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > 0.8 * (pool_diameter / 2)^2) {
      wall_zone <- wall_zone + 1
    } else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0) {
      quadrants[1] <- quadrants[1] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0) {
      quadrants[2] <- quadrants[2] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0) {
      quadrants[3] <- quadrants[3] + 1
    } else {
      quadrants[4] <- quadrants[4] + 1
    }

    if (length(track_x) > 100) # evaluate latency only after 100+ steps to be accurate
      {
        speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2))) # speed in meters/time step
        latency <- (length(track_x) - 1) * speed_ts / speed # convert to seconds
        if (latency > 60) # if more than a minute, stop
          {
            break
          }
      }
  }

  latency <- length(track_x) - 1 # latency in time steps
  wall_zone <- wall_zone / latency
  quadrants <- quadrants / latency
  speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2))) # speed in meters/time step

  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency))
}

run_trial_comb <- function(weights0_pc, weights0_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun, weight_wall, weight_place) # c
{
  # FIXED PARAMETERS OF THE EXPERIMENT

  pool_diameter <- 1.4 # Maze diameter in metres (m)
  platform_radius <- 0.06 # Platform radius

  N_dc <- 10 # Population of place cells
  N_pc <- 211 # Population of place cells # c
  N_ac <- 36 # Population of action cells
  which_pc <- 0 # c
  which_dc <- 0 # c
  which <- 0 # c

  dist <- 0
  wall_zone <- 0
  quadrants <- c(0, 0, 0, 0) # Percentage spent on each quadrant

  weights_pc <- weights0_pc # Initialize modifiable weights # c
  weights_dc <- weights0_dc # Initialize modifiable weights # c


  el_tr_dc <- matrix(rep(0, N_dc * N_ac), nrow = N_dc) # Initialize eligibility traces matrix # c
  el_tr_pc <- matrix(rep(0, N_pc * N_ac), nrow = N_pc) # Initialize eligibility traces matrix # c

  # Initialize trajectories
  track_x <- starting_x # Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0


  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2) # while out of platform
  {
    weights_dc <- weights_dc * (1 - noise) + matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult * noise # c
    weights_pc <- weights_pc * (1 - noise) + matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult * noise # c

    dist.to.wall <- pool_diameter / 2 - sqrt(track_x[length(track_x)]^2 + track_y[length(track_y)]^2)

    # Calculate DC and PC activation
    DC_activation <- rep(0, N_dc)
    for (i in 1:N_dc) {
      DC_activation[i] <- exp(-(dist.to.wall - DC[i])^2 / (2 * sigma_dc^2))
    }
    PC_activation <- rep(0, N_pc) # c
    for (i in 1:N_pc) { # c
      PC_activation[i] <- exp(-((track_x[length(track_x)] - PC_x[i])^2 + (track_y[length(track_y)] - PC_y[i])^2) / (2 * sigma_pc^2))
    } # c

    # Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1) {
      prevQ_pc <- AC_activation_pc[which_pc] # Displays the Q value before movement # c
      prevQ_dc <- AC_activation_dc[which_dc] # Displays the Q value before movement # c
    }
    # print(DC_activation)
    AC_activation_dc <- DC_activation %*% weights_dc # c
    AC_activation_pc <- PC_activation %*% weights_pc # c


    # Make an action
    ACsel_pc <- safe_action_probs(AC_activation_pc, beta) # c
    ACsel_dc <- safe_action_probs(AC_activation_dc, beta) # c


    shift <- round(
      (180 - atan2(track_y[length(track_y)], track_x[length(track_x)]) /
        pi * 180) / 10
    )
    if (shift == 36) {
      shift <- 0
    }


    if (shift < 0) {
      ACsel_dc <- c(ACsel_dc[(shift + 37):36], ACsel_dc[1:(shift + 36)])
    } else if (shift == 0) {
      ACsel_dc <- ACsel_dc
    } else {
      ACsel_dc <- c(ACsel_dc[(shift + 1):36], ACsel_dc[1:shift])
    }


    ACsel <- ACsel_dc * weight_wall + ACsel_pc * weight_place
    ACsel <- safe_action_probs(ACsel, 1)


    ASrand <- runif(1) # c
    which <- 1 # c
    ASsum <- ACsel[1] # c

    while (which < N_ac && ASsum <= ASrand) { # c
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }

    if (which + shift > 36) {
      which_dc <- which + shift - 36
    } else if (which + shift <= 0) {
      which_dc <- which + shift + 36
    } else {
      which_dc <- which + shift
    }


    which_pc <- which


    # Eligibility traces
    el_tr_pc <- el_tr_pc * etdecay # c
    el_tr_dc <- el_tr_dc * etdecay # c

    for (j in 1:N_ac) { # c
      itmp_pc <- min(abs(j - which_pc), N_ac - abs(j - which_pc))
      itmp_dc <- min(abs(j - which_dc), N_ac - abs(j - which_dc))
      actgaus_pc <- exp(-(itmp_pc * itmp_pc) / (2 * sigma_ac * sigma_ac))
      actgaus_dc <- exp(-(itmp_dc * itmp_dc) / (2 * sigma_ac * sigma_ac))
      el_tr_pc[, j] <- el_tr_pc[, j] + actgaus_pc * AC_activation_pc[j] * t(t(PC_activation))
      el_tr_dc[, j] <- el_tr_dc[, j] + actgaus_dc * AC_activation_dc[j] * t(t(DC_activation))
    }


    vel_x <- c(vel_x, (vel_x[length(vel_x)] + ac_const * cos(which / N_ac * 2 * pi)) * Vdecay)
    vel_y <- c(vel_y, (vel_y[length(vel_y)] + ac_const * sin(which / N_ac * 2 * pi)) * Vdecay)

    # velocity per time step (not second)
    track_x <- c(track_x, track_x[length(track_x)] + vel_x[length(vel_x)])
    track_y <- c(track_y, track_y[length(track_y)] + vel_y[length(vel_y)])

    # Check if not out of bounds, reset location & speed if so
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter / 2)^2) {
      ratio <- (track_x[length(track_x)]^2 + track_y[length(track_y)]^2) / ((pool_diameter / 2)^2)
      track_x[length(track_x)] <- track_x[length(track_x)] / sqrt(ratio)
      track_y[length(track_y)] <- track_y[length(track_y)] / sqrt(ratio)
      vel_x[length(vel_x)] <- track_x[length(track_x)] - track_x[length(track_x) - 1]
      vel_y[length(vel_y)] <- track_y[length(track_y)] - track_y[length(track_y) - 1]
    }


    if (length(track_x) > 2) {
      if ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 < platform_radius^2) {
        rew <- 10
      } # found platform - reward
      else if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (0.99 * pool_diameter / 2)^2) {
        rew <- -wall_pun
      } # hit wall - punishment
      else {
        rew <- 0
      } # didn't find - no reward

      currQ_pc <- AC_activation_pc[which_pc] # c
      currQ_dc <- AC_activation_dc[which_dc] # c
      tderr_pc <- rew + discf * currQ_pc - prevQ_pc # temporal difference error # c
      tderr_dc <- rew + discf * currQ_dc - prevQ_dc # temporal difference error # c
      weights_pc <- pmax(weights_pc + lrate * tderr_pc * el_tr_pc, 0) # c
      weights_dc <- pmax(weights_dc + lrate * tderr_dc * el_tr_dc, 0) # c
    }

    laststep <- sqrt((track_x[length(track_x)] - track_x[length(track_x) - 1])^2 + (track_y[length(track_y)] - track_y[length(track_y) - 1])^2)
    dist <- dist + laststep

    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > 0.8 * (pool_diameter / 2)^2) {
      wall_zone <- wall_zone + 1
    } else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0) {
      quadrants[1] <- quadrants[1] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0) {
      quadrants[2] <- quadrants[2] + 1
    } else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0) {
      quadrants[3] <- quadrants[3] + 1
    } else {
      quadrants[4] <- quadrants[4] + 1
    }

    if (length(track_x) > 100) # evaluate latency only after 100+ steps to be accurate
      {
        speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2))) # speed in meters/time step
        latency <- (length(track_x) - 1) * speed_ts / speed # convert to seconds
        if (latency > 60) # if more than a minute, stop
          {
            break
          }
      }
  }

  latency <- length(track_x) - 1 # latency in time steps
  wall_zone <- wall_zone / latency
  quadrants <- quadrants / latency
  speed_ts <- mean(sqrt((vel_x[-1]^2 + vel_y[-1]^2))) # speed in meters/time step
  # speed per action step from Hanbing
  speed_ps <- (vel_x[-1]^2 + vel_y[-1]^2)^0.5

  # time step
  time_step <- speed_ts / speed

  # mean turning angle
  vel <- as.matrix(data.frame(vel_x, vel_y))
  angle <- c()
  for (steps in 2:(length(vel_x) - 1)) {
    A <- vel[steps, ]
    B <- vel[steps + 1, ]
    angle <- c(angle, acos((A %*% B)[1, 1]) / norm(as.matrix(A) * norm(as.matrix(B))))
  }
  angle <- angle * 180 / pi
  mean_angle <- mean(angle)

  # speed standard deviation
  speed_std <- sd((vel_x[-1]^2 + vel_y[-1]^2)^0.5)
  speed_std <- speed_std / time_step

  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights_pc, weights_dc, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency, speed_std, speed_ps, mean_angle, time_step))
}

# ================== Parameters ===================
# Universal parameters
pool_diameter <- 1.4 # Maze diameter (\phi) in metres (m)
platform_radius <- 0.06 # Platform radius (m)
sigma_pc <- 0.1 # place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac <- 2 # action cell sigma (standard deviation), in action cells [1..3]
sigma_dc <- 0.1 # distance cell sigma (standard deviation), in meters [0.05..0.2]

etdecay <- 0.83 # Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
beta <- 6 # Exploration-exploitation factor (\beta) [0.5..12]
alpha <- 0.01 # Learning rate (\alpha) [0.005..0.02]
gamma <- 0.85 # Discount factor (\gamma) [0.75..0.95]

Vdecay <- 0.82 # velocity decay [0.75..0.95]
ac_const <- 0.02 # acceleration const [0.01..0.03]
Wnoise <- 0.0004 # Weight noise [0.0001..0.0007]
Wmult <- 0.1 # Weight multiplier [0.05..0.15]
hitwall <- 0.5 # punishment for hitting the wall [0..1]
speed <- 0.175 # mouse speed (m/s) [0.1..0.25]

Ntrials <- 4 # number of trials per day
Ndays <- 5 # number of days

N_ac <- 36 # Population of action cells [25..50]

# Place cell model specific parameters
N_pc <- 211 # Population of place cells [100..300]

# Distance cell model specific parameters
N_dc <- 10 # Population of distance cells [100..300]
variable_platform <- 1 # the platform is variable - 1, the platform is fixed - 0

# Combined model specific parameters
weight_wall <- 0.2
weight_place <- 1 - weight_wall

# Platform coordinates (immutable reference)
platform_x0 <- cos(-pi / 4) * pool_diameter / 4
platform_y0 <- sin(-pi / 4) * pool_diameter / 4
platform_x_mag <- abs(platform_x0)
platform_y_mag <- abs(platform_y0)

# Starting locations of the modeled animal
strad <- pool_diameter / 2 * 0.85
starting_xs <- strad * c(cos(pi / 6), cos(pi / 3), cos(7 * pi / 6), cos(4 * pi / 3))
starting_ys <- strad * c(sin(pi / 6), sin(pi / 3), sin(7 * pi / 6), sin(4 * pi / 3))

Nruns <- 40
PMs_pc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
PMs_dc <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
PMs_comb <- array(rep(0, 5 * Ndays * Ntrials), c(5, Ndays, Ntrials))
track_x_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
track_y_sum <- vector(mode = "list", length = Nruns * Ndays * Ntrials)
final_clean_pm_fixed <- data.frame()

# Runing all three models for Nruns runs, Ndays days, and Ntrials trials with environement reset for each run,
# and shared trial schedule for fair comparison between models within each run
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
  final_clean_pm_fixed <- rbind(
    final_clean_pm_fixed, clean_pm_dc, clean_pm_pc,
    clean_pm_comb
  )
}

final_clean_pm_fixed_wide <- pivot_wider(
  final_clean_pm_fixed,
  names_from = Parameter,
  values_from = value
)

# Statistical tests and data visualization
model <- aov(Dist ~ Day + Model, data = final_clean_pm_fixed_wide)
summary(model)

strategy_summary_Latency_fixed <- final_clean_pm_fixed_wide %>%
  group_by(Day, Model) %>%
  summarise(
    mean_latency = mean(Latency, na.rm = TRUE),
    # Calculate Standard Error (SE) for error bars
    se_latency = sd(Latency, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

strategy_summary_Latency_fixed$Day <- factor(strategy_summary_Latency_fixed$Day,
  levels = c("Day_1", "Day_2", "Day_3", "Day_4", "Day_5")
)

# 1. Plotting Latency with error bars for the fixed and variable platform condition (learning curves)
ggplot(
  strategy_summary_Latency_fixed,
  aes(x = Day, y = mean_latency, group = Model, color = Model)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_latency - se_latency, ymax = mean_latency + se_latency),
    width = 0.15,
    alpha = 0.6
  ) +
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

strategy_summary_Latency <- final_clean_pm_fixed_wide %>%
  group_by(Day, Model) %>%
  summarise(
    mean_latency = mean(Latency, na.rm = TRUE),
    # Calculate Standard Error (SE) for error bars
    se_latency = sd(Latency, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

strategy_summary_Latency$Day <- factor(strategy_summary_Latency$Day,
  levels = c("Day_1", "Day_2", "Day_3", "Day_4", "Day_5")
)

ggplot(
  strategy_summary_Latency,
  aes(x = Day, y = mean_latency, group = Model, color = Model)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = mean_latency - se_latency, ymax = mean_latency + se_latency),
    width = 0.15,
    alpha = 0.6
  ) +
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

# 2. target quad v.s. opposite quad preference for the fixed and variable platform condition
target_v_oppo_fixed <- final_clean_pm_fixed_wide %>%
  mutate(target_pref = Target_quad / (Target_quad + Oppo_quad)) %>%
  group_by(Day, Model) %>%
  summarise(
    mean_pref = mean(target_pref, na.rm = TRUE),
    se_pref = sd(target_pref, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(target_v_oppo_fixed, aes(x = Day, y = mean_pref, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_pref - se_pref, ymax = mean_pref + se_pref),
    position = position_dodge(width = 0.8),
    width = 0.2,
    alpha = 0.6
  ) +
  labs(
    y = "Target Preference",
    x = "Day"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c(
    "Combined" = "#BFDFD2",
    "Distance_cell" = "#51999F",
    "Place_cell" = "#7BC0CD"
  )) +
  theme_classic()

target_v_oppo <- final_clean_pm_fixed_wide %>%
  mutate(target_pref = Target_quad / (Target_quad + Oppo_quad)) %>%
  group_by(Day, Model) %>%
  summarise(
    mean_pref = mean(target_pref, na.rm = TRUE),
    se_pref = sd(target_pref, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(target_v_oppo, aes(x = Day, y = mean_pref, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_pref - se_pref, ymax = mean_pref + se_pref),
    position = position_dodge(width = 0.8),
    width = 0.2,
    alpha = 0.6
  ) +
  labs(
    y = "Target Preference",
    x = "Day"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c(
    "Combined" = "#BFDFD2",
    "Distance_cell" = "#51999F",
    "Place_cell" = "#7BC0CD"
  )) +
  theme_classic()
