#' Perform computer simulations of candidate behavior
#'
#' @description `neodowns()` simulates the changes in candidate positions under three electoral strategies.
#'
#' @param data simulated data from `sim_data`
#' @param n_iter Number of iterations.
#' @param mu_c Parameter that controls for the average level of spatial voting.
#' @param mu_b Parameter that controls for the average level of co-ethnic voting.
#' @param sigma_c Parameter that controls for the standard deviation around c
#' @param sigma_b Parameter that controls for the standard deviation around b
#' @param eps_sd Parameter that controls for the level of idiosyncratic factors.
#' @param unit Unit with which candidates move in the two-dimensional space.
#' @param force Logical value for forcing moderate candidates taking more moderate positions than extreme candidates within each group.
#' @param while_max Number of trials to apply `force`
#' @param boost Unit by which a simulation pushes a stuck candidate out of a trap.
#'  @return
#' A list, which has the following components:
#'  \item{d_voters}{a matrix containing the simulated ideological positions of voters.}
#'  \item{d_cands}{a matrix containing the simulated ideological positions of candidates.}
#' @importFrom dplyr case_when `%>%` mutate select left_join rename tibble filter bind_rows pull
#' @importFrom progress progress_bar
#' @importFrom rlang := !!
#' @export


neodowns <- function(data,
                     n_iter = 3000,
                     mu_c = 1, # average spatial voting
                     mu_b = 2, # average co-ethnic voting
                     sigma_c = 0.5,
                     sigma_b = 0.5,
                     eps_sd = 0.5,
                     unit = 0.05,
                     force = TRUE,
                     while_max = 200,
                     boost = 6,
                     seed = 14231) {

# Only for debuggin
  n_iter = 100
  mu_c = 1 # average spatial voting
  mu_b = 2 # average co-ethnic voting
  sigma_c = 0.5
  sigma_b = 0.5
  eps_sd = 0.5
  unit = 0.05
  force = TRUE
  while_max = 200
  boost = 6
  seed = 14231


# Set up -----------------------------------------------------------------------
  set.seed(seed)

  d_voters <- data$gen_voters
  d_cands <- data$gen_cands %>%
    mutate(dist = sqrt((x - 0)^2 + (y - 0)^2))

  n_iter <- n_iter
  J <- dim(d_cands)[1] # number of candidates

  # Fixed parameters, created outside chains
  N_voters <- dim(d_voters)[1]
  m_vec <- rep(1:(J/2), 2) # used in J-loop


# Chain initialization ---------------------------------------------------------

  # Feeding the voter distribution to "initial chain" for both systems
  init <- chain <- d_voters

  init <- init %>%
    mutate(c = rnorm(n = N_voters, mean = mu_c, sd = sigma_c),
           b = rnorm(n = N_voters, mean = mu_b, sd = sigma_b))

  # c <- init$c
  # b <- init$b

  # isseue with c
  # must change the later part simultanesouly

  # c <- rnorm(n = N_voters, mean = mu_c, sd = sigma_c)
  # b <- rnorm(n = N_voters, mean = mu_b, sd = sigma_b)

  # Create distance for all candidate

  # for (j in 1:J) {
  #
  #   init <- init %>%
  #     mutate("D{j}" := sqrt((x - d_cands$x[{j}])^2 + (y - d_cands$y[{j}])^2),
  #            "m{j}" := ifelse(init$ethnic_group != paste0("Group ", m_vec[j]), 1, 0),
  #            "eps{j}" := rnorm(n = N_voters, sd = eps_sd),
  #            "V{j}" := -c * !! as.name(paste0("D",j))
  #            -b * !! as.name(paste0("m",j))
  #            + !! as.name(paste0("eps",j))
  #     )
  # }

  for (j in seq_len(J)) {
    D_j <- sqrt((init$x - d_cands$x[j])^2 + (init$y - d_cands$y[j])^2)
    m_j <- ifelse(init$ethnic_group != paste0("Group ", m_vec[j]), 1, 0)
    eps_j <- rnorm(n = N_voters, sd = eps_sd)
    V_j <- -init$c * D_j - init$b * m_j + eps_j
    init[[paste0("V", j)]] <- V_j
  }




  # # Mismatch Indicator
  # m1 <- ifelse(init$ethnic_group != "Group 1", 1, 0)
  # m2 <- ifelse(init$ethnic_group != "Group 2", 1, 0)
  # m3 <- ifelse(init$ethnic_group != "Group 3", 1, 0)
  # m4 <- ifelse(init$ethnic_group != "Group 1", 1, 0)
  # m5 <- ifelse(init$ethnic_group != "Group 2", 1, 0)
  # m6 <- ifelse(init$ethnic_group != "Group 3", 1, 0)
  #
  # # Random Errors, created outside chains
  # eps1 <- rnorm(n = N_voters, sd = eps_sd)
  # eps2 <- rnorm(n = N_voters, sd = eps_sd)
  # eps3 <- rnorm(n = N_voters, sd = eps_sd)
  # eps4 <- rnorm(n = N_voters, sd = eps_sd)
  # eps5 <- rnorm(n = N_voters, sd = eps_sd)
  # eps6 <- rnorm(n = N_voters, sd = eps_sd)

  # # Computing the observed utility for each party (This is where we specify utility functions)
  # init$V1 <- -1 * c * init$D1 - b * m1 + eps1
  # init$V2 <- -1 * c * init$D2 - b * m2 + eps2
  # init$V3 <- -1 * c * init$D3 - b * m3 + eps3
  # init$V4 <- -1 * c * init$D4 - b * m4 + eps4
  # init$V5 <- -1 * c * init$D5 - b * m5 + eps5
  # init$V6 <- -1 * c * init$D6 - b * m6 + eps6

  # (1): FPTP
  # Computing the first-choice probability that each voter votes for each party
  V_den <- exp(init$V1) + exp(init$V2) + exp(init$V3) + exp(init$V4) + exp(init$V5) + exp(init$V6)
  init$P1 <- exp(init$V1) / V_den
  init$P2 <- exp(init$V2) / V_den
  init$P3 <- exp(init$V3) / V_den
  init$P4 <- exp(init$V4) / V_den
  init$P5 <- exp(init$V5) / V_den
  init$P6 <- exp(init$V6) / V_den

  # Probability of co-ethnic voting for each voter
  init$epsilon1 <- case_when(
    init$ethnic_group == "Group 1" ~ init$P1 + init$P4,
    init$ethnic_group == "Group 2" ~ init$P2 + init$P5,
    init$ethnic_group == "Group 3" ~ init$P3 + init$P6
  )

  # # Check
  # mean(is.na(init$epsilon1))
  # hist(init$epsilon1)

  # Compute the first-choice probability score (sum of all support probabilities)
  P1_score <- P2_score <- P3_score <- P4_score <- P5_score <- P6_score <- NA
  P1_score[1] <- sum(init$P1)
  P2_score[1] <- sum(init$P2)
  P3_score[1] <- sum(init$P3)
  P4_score[1] <- sum(init$P4)
  P5_score[1] <- sum(init$P5)
  P6_score[1] <- sum(init$P6)

  # First-choice probabilities (for FPTP, RCV1, RCV2)
  P_vec <- NA
  P_vec <- rbind(P1_score, P2_score, P3_score, P4_score, P5_score, P6_score)

  # # CHECK: Sum of all first-choice probabilities: they need to sum up to N=1000
  # print(sum(P_vec)) # First-choice FPTP

  # ========================================================================#
  # Updating party positions (initial move_max1)
  move_max1 <- d_cands
  theta <- runif(n = 6, min = 0, max = 360) # Generating random angels

  move_max1$x <- d_cands$x + cos(theta / 180 * pi) * (unit / 2) # Compute x-value given theta SETTING THIS 0.01 (8/3/2022)
  move_max1$y <- d_cands$y + sin(theta / 180 * pi) * (unit / 2) # Compute y-value given theta
  move_max1$dist <- sqrt((move_max1$x - 0)^2 + (move_max1$y - 0)^2)
  move_max1$moderation <- move_max1$dist / d_cands$dist

  # ========================================================================#

  iter <- NA
  cands_max1 <- list()
  voter_max1 <- list()
  new_theta <- NA

  pb <- progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = n_iter,
    complete = "=",   # Completion bar character
    incomplete = "-", # Incomplete bar character
    current = ">",    # Current bar character
    clear = FALSE,    # If TRUE, clears the bar when finish
    width = 100       # Width of the progress bar
  )


  # Loop for Markov Chains -------------------------------------------------------
  # Iterating Updates
  iter <- 1
  for (t in 2:n_iter) {
    # Updates the current state
    pb$tick()
    iter[t] <- t

    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # (1): FPTP
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # Compute the voting probabilities again

    chain$D1 <- sqrt((chain$x - move_max1$x[1])^2 + (chain$y - move_max1$y[1])^2)
    chain$D2 <- sqrt((chain$x - move_max1$x[2])^2 + (chain$y - move_max1$y[2])^2)
    chain$D3 <- sqrt((chain$x - move_max1$x[3])^2 + (chain$y - move_max1$y[3])^2)
    chain$D4 <- sqrt((chain$x - move_max1$x[4])^2 + (chain$y - move_max1$y[4])^2)
    chain$D5 <- sqrt((chain$x - move_max1$x[5])^2 + (chain$y - move_max1$y[5])^2)
    chain$D6 <- sqrt((chain$x - move_max1$x[6])^2 + (chain$y - move_max1$y[6])^2)

    # Computing the observed utility for each party (This is where we specify utility functions)
    chain$V1 <- -1 * c * chain$D1 - b * m1 + eps1
    chain$V2 <- -1 * c * chain$D2 - b * m2 + eps2
    chain$V3 <- -1 * c * chain$D3 - b * m3 + eps3
    chain$V4 <- -1 * c * chain$D4 - b * m4 + eps4
    chain$V5 <- -1 * c * chain$D5 - b * m5 + eps5
    chain$V6 <- -1 * c * chain$D6 - b * m6 + eps6

    # Computing the probability that each voter votes for each party (this is fixed)
    den <- exp(chain$V1) + exp(chain$V2) + exp(chain$V3) + exp(chain$V4) + exp(chain$V5) + exp(chain$V6)

    chain$P1 <- exp(chain$V1) / den
    chain$P2 <- exp(chain$V2) / den
    chain$P3 <- exp(chain$V3) / den
    chain$P4 <- exp(chain$V4) / den
    chain$P5 <- exp(chain$V5) / den
    chain$P6 <- exp(chain$V6) / den

    # Compute the probability score (sum of all support probabilities)
    P1_score[t] <- sum(chain$P1)
    P2_score[t] <- sum(chain$P2)
    P3_score[t] <- sum(chain$P3)
    P4_score[t] <- sum(chain$P4)
    P5_score[t] <- sum(chain$P5)
    P6_score[t] <- sum(chain$P6)

    P_vec <- rbind(P1_score, P2_score, P3_score, P4_score, P5_score, P6_score)
    sum(P_vec[, t]) # This must sum up to N=1000

    # CHECK
    try(if (sum(P_vec[, 1]) != N_voters) stop("FPTP: 1st-choice ranking probabilities do not sum up to one"))

    # ============================================================================#
    #  Updating party locations
    # ============================================================================#
    new_theta <- NA # Initialize
    move_max1$new_dist <- NA # Initialize
    for (i in 1:6) {
      # For Moderate Parties
      if (i < 4) {
        new_theta[i] <- ifelse(P_vec[i, t] >= P_vec[i, t - 1], # If new position has higher prob score
                               theta[i], # Keep going
                               runif(n = 1, min = theta[i] + 90, max = theta[i] + 270)
        ) # New position on the other side

        move_max1$x[i] <- move_max1$x[i] + cos(new_theta[i] / 180 * pi) * unit # Compute x-value given theta
        move_max1$y[i] <- move_max1$y[i] + sin(new_theta[i] / 180 * pi) * unit # Compute y-value given theta
        move_max1$new_dist[i] <- sqrt((move_max1$x[i] - 0)^2 + (move_max1$y[i] - 0)^2) # New Distance
        move_max1$moderation[i] <- move_max1$new_dist[i] / move_max1$dist[i] # New Location / Initial Location
        theta[i] <- new_theta[i] # Update for the next iteration


        # For Extreme Parties
      } else {
        new_theta[i] <- ifelse(P_vec[i, t] >= P_vec[i, t - 1], # If new position has higher prob score
                               theta[i], # Keep going
                               runif(n = 1, min = theta[i] + 90, max = theta[i] + 270)
        ) # New position on the other side

        # Updating party locations
        propose_x <- move_max1$x[i] + cos(new_theta[i] / 180 * pi) * unit # Compute x-value given theta
        propose_y <- move_max1$y[i] + sin(new_theta[i] / 180 * pi) * unit # Compute y-value given theta
        move_max1$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

        # Change the direction until extreme parties become more extreme than moderate parties
        if (force == TRUE) {
          while_iter <- 1
          while (move_max1$new_dist[i] < move_max1$new_dist[i - 3] & while_iter <= while_max) {
            new_theta[i] <- runif(n = 1, min = 0, max = 360) # Anywhere is okay as long as extreme parties can get out of the trap
            propose_x <- move_max1$x[i] + cos(new_theta[i] / 180 * pi) * unit * boost # Compute x-value given theta
            propose_y <- move_max1$y[i] + sin(new_theta[i] / 180 * pi) * unit * boost # Compute y-value given theta
            move_max1$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

            while_iter <- while_iter + 1
          } # Close while () loop

          # If extreme parties cannot find a way out, STAY (NO move_max1)
          if (while_iter == while_max) {
            propose_x <- move_max1$x[i]
            propose_y <- move_max1$y[i]
          }
        } # Closing Assumption 1 condition


        move_max1$x[i] <- propose_x # Saving the accepted location
        move_max1$y[i] <- propose_y # Saving the accepted location
        move_max1$moderation[i] <- move_max1$new_dist[i] / move_max1$dist[i] # New Location / Initial Location
        theta[i] <- new_theta[i] # Update for the next iteration
      } # Close else{}
    } # Close for () loop

    try(if (sum(P_vec[, 1]) != N_voters) stop("FPTP: 1st-choice probabilities do not sum up to one"))



    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # Keep cands_max1 of Party Positions
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #

    move_max1$iter <- t

    # Save candidate information
    cands_max1[[t]] <- move_max1


    # Save voter information
    chain$iter <- t
    voter_max1[[t]] <- chain

    iter[t] <- t

  }




  # Summarizing results ----------------------------------------------------------

  cands_max1[[1]] <- d_cands %>%
    mutate(moderation = 1,
           new_dist = dist,
           iter = 1)

  cands_max1 <- as.data.frame(do.call(rbind, cands_max1)) %>%
    mutate(system = "max1")

  # combine all results
  cands_chains <- rbind(cands_max1) %>%
    tibble()


  voter_max1 <- as.data.frame(do.call(rbind, voter_max1)) %>%
    mutate(system = "max1")

  # combine all results
  voter_chains <- dplyr::bind_rows(voter_max1) %>%
    tibble()

  out <- list(
    voters = voter_chains,
    cands = cands_chains
  )

  return(out)

  # FUTURE WORK: Add the summary table for each experiment
  # ADD some warning based on "check"s
}



# Burn-in
# Random Sampling from the Chains (to address autocorrelation)
# Take every 25th chain
burn_samp <- function(chain, burnin = 1000, int = 1) {
  chain <- pull(chain) # Vectorize the input chain
  burned <- chain[burnin:length(chain)] # Drop burn-in periods ("burnin")
  sampled <- burned[seq(1, length(burned), by = int)] # Sample each "int" statistic
  return(sampled)
}
