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
#' @importFrom dplyr case_when `%>%` mutate select left_join rename tibble filter bind_rows
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


  d_voters <- data$gen_voters
  d_cands <- data$gen_cands %>%
    mutate(dist = sqrt((x - 0)^2 + (y - 0)^2))

  set.seed(seed)

  n_iter <- n_iter
  J <- dim(d_cands)[1] # number of candidates

  # Feeding the voter distribution to "initial chain" for both systems
  init <- chain <- chain_rcv <- chain_rcv_t <- d_voters

  # chain <- d_voters       # Feeding the voter distribution to Markov "chain" for FPTP
  # chain_rcv <- d_voters   # Feeding the voter distribution to Markov "chain" for RCV -- Second choice
  # chain_rcv_t <- d_voters # Feeding the voter distribution to Markov "chain" for RCV -- Third Choice
  # chain and chain_rcv appear in the t loop below (used in the iteration)


  # Fixed parameters, created outside chains
  N_voters <- dim(d_voters)[1]
  m_vec <- rep(1:(J/2), 2) # used in J-loop

  init <- init %>%
    mutate(c = rnorm(n = N_voters, mean = mu_c, sd = sigma_c),
           b = rnorm(n = N_voters, mean = mu_b, sd = sigma_b))

  c <- init$c
  b <- init$b

  # isseue with c
  # must change the later part simultanesouly

  # c <- rnorm(n = N_voters, mean = mu_c, sd = sigma_c)
  # b <- rnorm(n = N_voters, mean = mu_b, sd = sigma_b)

  # Create distance for all candidate

    for (j in 1:J) {

    init <- init %>%
      mutate("D{j}" := sqrt((x - d_cands$x[{j}])^2 + (y - d_cands$y[{j}])^2),
             "m{j}" := ifelse(init$ethnic_group != paste0("Group ", m_vec[j]), 1, 0),
             "eps{j}" := rnorm(n = N_voters, sd = eps_sd),
             "V{j}" := -c * !! as.name(paste0("D",j))
                       -b * !! as.name(paste0("m",j))
                       + !! as.name(paste0("eps",j))
             )
  }

  # Mismatch Indicator
  m1 <- ifelse(init$ethnic_group != "Group 1", 1, 0)
  m2 <- ifelse(init$ethnic_group != "Group 2", 1, 0)
  m3 <- ifelse(init$ethnic_group != "Group 3", 1, 0)
  m4 <- ifelse(init$ethnic_group != "Group 1", 1, 0)
  m5 <- ifelse(init$ethnic_group != "Group 2", 1, 0)
  m6 <- ifelse(init$ethnic_group != "Group 3", 1, 0)

  # Random Errors, created outside chains
  eps1 <- rnorm(n = N_voters, sd = eps_sd)
  eps2 <- rnorm(n = N_voters, sd = eps_sd)
  eps3 <- rnorm(n = N_voters, sd = eps_sd)
  eps4 <- rnorm(n = N_voters, sd = eps_sd)
  eps5 <- rnorm(n = N_voters, sd = eps_sd)
  eps6 <- rnorm(n = N_voters, sd = eps_sd)

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
  P1_score_rcv <- P2_score_rcv <- P3_score_rcv <- P4_score_rcv <- P5_score_rcv <- P6_score_rcv <- NA
  P1_score_rcv_t <- P2_score_rcv_t <- P3_score_rcv_t <- P4_score_rcv_t <- P5_score_rcv_t <- P6_score_rcv_t <- NA
  P1_score[1] <- sum(init$P1)
  P2_score[1] <- sum(init$P2)
  P3_score[1] <- sum(init$P3)
  P4_score[1] <- sum(init$P4)
  P5_score[1] <- sum(init$P5)
  P6_score[1] <- sum(init$P6)
  P1_score_rcv[1] <- sum(init$P1)
  P2_score_rcv[1] <- sum(init$P2)
  P3_score_rcv[1] <- sum(init$P3)
  P4_score_rcv[1] <- sum(init$P4)
  P5_score_rcv[1] <- sum(init$P5)
  P6_score_rcv[1] <- sum(init$P6)
  P1_score_rcv_t[1] <- sum(init$P1)
  P2_score_rcv_t[1] <- sum(init$P2)
  P3_score_rcv_t[1] <- sum(init$P3)
  P4_score_rcv_t[1] <- sum(init$P4)
  P5_score_rcv_t[1] <- sum(init$P5)
  P6_score_rcv_t[1] <- sum(init$P6)

  # First-choice probabilities (for FPTP, RCV1, RCV2)
  P_vec <- NA
  P_vec_rcv <- NA
  P_vec_rcv_t <- NA
  P_vec <- rbind(P1_score, P2_score, P3_score, P4_score, P5_score, P6_score)
  P_vec_rcv <- rbind(P1_score_rcv, P2_score_rcv, P3_score_rcv, P4_score_rcv, P5_score_rcv, P6_score_rcv)
  P_vec_rcv_t <- rbind(P1_score_rcv_t, P2_score_rcv_t, P3_score_rcv_t, P4_score_rcv_t, P5_score_rcv_t, P6_score_rcv_t)

  # # CHECK: Sum of all first-choice probabilities: they need to sum up to N=1000
  # print(sum(P_vec)) # First-choice FPTP
  # print(sum(P_vec_rcv)) # First-choice RCV (only up to second)
  # print(sum(P_vec_rcv_t)) # First-choice RCV  (up to third)

  # (1): RCV - Second Choice Probability
  # Computing the second-choice probability that each voter votes for each party
  V_den_ex1 <- exp(init$V2) + exp(init$V3) + exp(init$V4) + exp(init$V5) + exp(init$V6) # Second-choice prob without P1
  V_den_ex2 <- exp(init$V1) + exp(init$V3) + exp(init$V4) + exp(init$V5) + exp(init$V6) # Second-choice prob without P2
  V_den_ex3 <- exp(init$V1) + exp(init$V2) + exp(init$V4) + exp(init$V5) + exp(init$V6) # Second-choice prob without P3
  V_den_ex4 <- exp(init$V1) + exp(init$V2) + exp(init$V3) + exp(init$V5) + exp(init$V6) # Second-choice prob without P4
  V_den_ex5 <- exp(init$V1) + exp(init$V2) + exp(init$V3) + exp(init$V4) + exp(init$V6) # Second-choice prob without P5
  V_den_ex6 <- exp(init$V1) + exp(init$V2) + exp(init$V3) + exp(init$V4) + exp(init$V5) # Second-choice prob without P6

  # Average second-choice probabilities (the probability of ranking (A,B,...))
  init$P1_s <- (exp(init$V1) / V_den_ex2) * init$P2 + (exp(init$V1) / V_den_ex3) * init$P3 +
    (exp(init$V1) / V_den_ex4) * init$P4 + (exp(init$V1) / V_den_ex5) * init$P5 + (exp(init$V1) / V_den_ex6) * init$P6
  init$P2_s <- (exp(init$V2) / V_den_ex1) * init$P1 + (exp(init$V2) / V_den_ex3) * init$P3 +
    (exp(init$V2) / V_den_ex4) * init$P4 + (exp(init$V2) / V_den_ex5) * init$P5 + (exp(init$V2) / V_den_ex6) * init$P6
  init$P3_s <- (exp(init$V3) / V_den_ex1) * init$P1 + (exp(init$V3) / V_den_ex2) * init$P2 +
    (exp(init$V3) / V_den_ex4) * init$P4 + (exp(init$V3) / V_den_ex5) * init$P5 + (exp(init$V3) / V_den_ex6) * init$P6
  init$P4_s <- (exp(init$V4) / V_den_ex1) * init$P1 + (exp(init$V4) / V_den_ex2) * init$P2 +
    (exp(init$V4) / V_den_ex3) * init$P3 + (exp(init$V4) / V_den_ex5) * init$P5 + (exp(init$V4) / V_den_ex6) * init$P6
  init$P5_s <- (exp(init$V5) / V_den_ex1) * init$P1 + (exp(init$V5) / V_den_ex2) * init$P2 +
    (exp(init$V5) / V_den_ex3) * init$P3 + (exp(init$V5) / V_den_ex4) * init$P4 + (exp(init$V5) / V_den_ex6) * init$P6
  init$P6_s <- (exp(init$V6) / V_den_ex1) * init$P1 + (exp(init$V6) / V_den_ex2) * init$P2 +
    (exp(init$V6) / V_den_ex3) * init$P3 + (exp(init$V6) / V_den_ex4) * init$P4 + (exp(init$V6) / V_den_ex5) * init$P5

  # # CHECK: mean should be 1
  # mean(init$P1_s > 0 & init$P1_s < 1)
  # mean(init$P2_s > 0 & init$P2_s < 1)
  # mean(init$P3_s > 0 & init$P3_s < 1)
  # mean(init$P4_s > 0 & init$P4_s < 1)
  # mean(init$P5_s > 0 & init$P5_s < 1)
  # mean(init$P6_s > 0 & init$P6_s < 1)
  #
  # # Check that maximum probability of choosing candidate 1 1st and choosing candidate 1 second does not exceeed 1
  # max(init$P1 + init$P1_s)
  # max(init$P2 + init$P2_s)
  # max(init$P3 + init$P3_s)
  # max(init$P4 + init$P4_s)
  # max(init$P5 + init$P5_s)
  # max(init$P6 + init$P6_s)

  # Compute Ethnic Polarization in initial second choices
  init$epsilon2 <- case_when(
    init$ethnic_group == "Group 1" ~ init$P1_s + init$P4_s,
    init$ethnic_group == "Group 2" ~ init$P2_s + init$P5_s,
    init$ethnic_group == "Group 3" ~ init$P3_s + init$P6_s
  )


  # Compute the SECOND-choice probability score (sum of all support probabilities)
  P1_s_score <- NA
  P2_s_score <- NA
  P3_s_score <- NA
  P4_s_score <- NA
  P5_s_score <- NA
  P6_s_score <- NA
  P1_st_score <- NA
  P2_st_score <- NA
  P3_st_score <- NA
  P4_st_score <- NA
  P5_st_score <- NA
  P6_st_score <- NA

  P1_s_score[1] <- sum(init$P1_s)
  P2_s_score[1] <- sum(init$P2_s)
  P3_s_score[1] <- sum(init$P3_s)
  P4_s_score[1] <- sum(init$P4_s)
  P5_s_score[1] <- sum(init$P5_s)
  P6_s_score[1] <- sum(init$P6_s)
  P1_st_score[1] <- sum(init$P1_s)
  P2_st_score[1] <- sum(init$P2_s)
  P3_st_score[1] <- sum(init$P3_s)
  P4_st_score[1] <- sum(init$P4_s)
  P5_st_score[1] <- sum(init$P5_s)
  P6_st_score[1] <- sum(init$P6_s)

  P_s_vec <- NA
  P_st_vec <- NA
  P_s_vec <- rbind(P1_s_score, P2_s_score, P3_s_score, P4_s_score, P5_s_score, P6_s_score)
  P_st_vec <- rbind(P1_st_score, P2_st_score, P3_st_score, P4_st_score, P5_st_score, P6_st_score)

  # # They must sum up to N=1000
  # print(sum(P_s_vec)) # GOOD
  # print(sum(P_st_vec)) # GOOD


  # (3): RCV - Third Choice Probability
  # Computing the third-choice probability that each voter votes for each party
  V_den_ex12 <- exp(init$V3) + exp(init$V4) + exp(init$V5) + exp(init$V6) # Third-choice prob without P1 & P2
  V_den_ex13 <- exp(init$V2) + exp(init$V4) + exp(init$V5) + exp(init$V6) # Third-choice prob without P1 & P3
  V_den_ex14 <- exp(init$V2) + exp(init$V3) + exp(init$V5) + exp(init$V6) # Third-choice prob without P1 & P4
  V_den_ex15 <- exp(init$V2) + exp(init$V3) + exp(init$V4) + exp(init$V6) # Third-choice prob without P1 & P5
  V_den_ex16 <- exp(init$V2) + exp(init$V3) + exp(init$V4) + exp(init$V5) # Third-choice prob without P1 & P6
  V_den_ex23 <- exp(init$V1) + exp(init$V4) + exp(init$V5) + exp(init$V6) # Third-choice prob without P2 & P3
  V_den_ex24 <- exp(init$V1) + exp(init$V3) + exp(init$V5) + exp(init$V6) # Third-choice prob without P2 & P4
  V_den_ex25 <- exp(init$V1) + exp(init$V3) + exp(init$V4) + exp(init$V6) # Third-choice prob without P2 & P5
  V_den_ex26 <- exp(init$V1) + exp(init$V3) + exp(init$V4) + exp(init$V5) # Third-choice prob without P2 & P6
  V_den_ex34 <- exp(init$V1) + exp(init$V2) + exp(init$V5) + exp(init$V6) # Third-choice prob without P3 & P4
  V_den_ex35 <- exp(init$V1) + exp(init$V2) + exp(init$V4) + exp(init$V6) # Third-choice prob without P3 & P5
  V_den_ex36 <- exp(init$V1) + exp(init$V2) + exp(init$V4) + exp(init$V5) # Third-choice prob without P3 & P6
  V_den_ex45 <- exp(init$V1) + exp(init$V2) + exp(init$V3) + exp(init$V6) # Third-choice prob without P4 & P5
  V_den_ex46 <- exp(init$V1) + exp(init$V2) + exp(init$V3) + exp(init$V5) # Third-choice prob without P4 & P6
  V_den_ex56 <- exp(init$V1) + exp(init$V2) + exp(init$V3) + exp(init$V4) # Third-choice prob without P5 & P6


  # Weighted Average third-choice probabilities
  # exp(party A)/(sum exp(remaining)) * [prob(second|first)prob(first) + prob(second|first)prob(first)]
  init$P1_t <- (exp(init$V1) / V_den_ex23) * ((exp(init$V3) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex3) * init$P3) + # Prob of choosing 2 then 3 + Prob of choosing 3 then 2
    (exp(init$V1) / V_den_ex24) * ((exp(init$V4) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex4) * init$P4) +
    (exp(init$V1) / V_den_ex25) * ((exp(init$V5) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex5) * init$P5) +
    (exp(init$V1) / V_den_ex26) * ((exp(init$V6) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex6) * init$P6) +
    (exp(init$V1) / V_den_ex34) * ((exp(init$V4) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex4) * init$P4) +
    (exp(init$V1) / V_den_ex35) * ((exp(init$V5) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex5) * init$P5) +
    (exp(init$V1) / V_den_ex36) * ((exp(init$V6) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex6) * init$P6) +
    (exp(init$V1) / V_den_ex45) * ((exp(init$V5) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex5) * init$P5) +
    (exp(init$V1) / V_den_ex46) * ((exp(init$V6) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex6) * init$P6) +
    (exp(init$V1) / V_den_ex56) * ((exp(init$V6) / V_den_ex5) * init$P5 + (exp(init$V5) / V_den_ex6) * init$P6)

  init$P2_t <- (exp(init$V2) / V_den_ex13) * ((exp(init$V3) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex3) * init$P3) + #
    (exp(init$V2) / V_den_ex14) * ((exp(init$V4) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex4) * init$P4) +
    (exp(init$V2) / V_den_ex15) * ((exp(init$V5) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex5) * init$P5) +
    (exp(init$V2) / V_den_ex16) * ((exp(init$V6) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex6) * init$P6) +
    (exp(init$V2) / V_den_ex34) * ((exp(init$V4) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex4) * init$P4) +
    (exp(init$V2) / V_den_ex35) * ((exp(init$V5) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex5) * init$P5) +
    (exp(init$V2) / V_den_ex36) * ((exp(init$V6) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex6) * init$P6) +
    (exp(init$V2) / V_den_ex45) * ((exp(init$V5) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex5) * init$P5) +
    (exp(init$V2) / V_den_ex46) * ((exp(init$V6) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex6) * init$P6) +
    (exp(init$V2) / V_den_ex56) * ((exp(init$V6) / V_den_ex5) * init$P5 + (exp(init$V5) / V_den_ex6) * init$P6)

  init$P3_t <- (exp(init$V3) / V_den_ex12) * ((exp(init$V2) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex2) * init$P2) + #
    (exp(init$V3) / V_den_ex14) * ((exp(init$V4) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex4) * init$P4) +
    (exp(init$V3) / V_den_ex15) * ((exp(init$V5) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex5) * init$P5) +
    (exp(init$V3) / V_den_ex16) * ((exp(init$V6) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex6) * init$P6) +
    (exp(init$V3) / V_den_ex24) * ((exp(init$V4) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex4) * init$P4) +
    (exp(init$V3) / V_den_ex25) * ((exp(init$V5) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex5) * init$P5) +
    (exp(init$V3) / V_den_ex26) * ((exp(init$V6) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex6) * init$P6) +
    (exp(init$V3) / V_den_ex45) * ((exp(init$V5) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex5) * init$P5) +
    (exp(init$V3) / V_den_ex46) * ((exp(init$V6) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex6) * init$P6) +
    (exp(init$V3) / V_den_ex56) * ((exp(init$V6) / V_den_ex5) * init$P5 + (exp(init$V5) / V_den_ex6) * init$P6)

  init$P4_t <- (exp(init$V4) / V_den_ex12) * ((exp(init$V2) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex2) * init$P2) + #
    (exp(init$V4) / V_den_ex13) * ((exp(init$V3) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex3) * init$P3) +
    (exp(init$V4) / V_den_ex15) * ((exp(init$V5) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex5) * init$P5) +
    (exp(init$V4) / V_den_ex16) * ((exp(init$V6) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex6) * init$P6) +
    (exp(init$V4) / V_den_ex23) * ((exp(init$V3) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex3) * init$P3) +
    (exp(init$V4) / V_den_ex25) * ((exp(init$V5) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex5) * init$P5) +
    (exp(init$V4) / V_den_ex26) * ((exp(init$V6) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex6) * init$P6) +
    (exp(init$V4) / V_den_ex35) * ((exp(init$V5) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex5) * init$P5) +
    (exp(init$V4) / V_den_ex36) * ((exp(init$V6) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex6) * init$P6) +
    (exp(init$V4) / V_den_ex56) * ((exp(init$V6) / V_den_ex5) * init$P5 + (exp(init$V5) / V_den_ex6) * init$P6)

  init$P5_t <- (exp(init$V5) / V_den_ex12) * ((exp(init$V2) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex2) * init$P2) + #
    (exp(init$V5) / V_den_ex13) * ((exp(init$V3) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex3) * init$P3) +
    (exp(init$V5) / V_den_ex14) * ((exp(init$V4) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex4) * init$P4) +
    (exp(init$V5) / V_den_ex16) * ((exp(init$V6) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex6) * init$P6) +
    (exp(init$V5) / V_den_ex23) * ((exp(init$V3) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex3) * init$P3) +
    (exp(init$V5) / V_den_ex24) * ((exp(init$V4) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex4) * init$P4) +
    (exp(init$V5) / V_den_ex26) * ((exp(init$V6) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex6) * init$P6) +
    (exp(init$V5) / V_den_ex34) * ((exp(init$V4) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex4) * init$P4) +
    (exp(init$V5) / V_den_ex36) * ((exp(init$V6) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex6) * init$P6) +
    (exp(init$V5) / V_den_ex46) * ((exp(init$V6) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex6) * init$P6)

  init$P6_t <- (exp(init$V6) / V_den_ex12) * ((exp(init$V2) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex2) * init$P2) +
    (exp(init$V6) / V_den_ex13) * ((exp(init$V3) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex3) * init$P3) +
    (exp(init$V6) / V_den_ex14) * ((exp(init$V4) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex4) * init$P4) +
    (exp(init$V6) / V_den_ex15) * ((exp(init$V5) / V_den_ex1) * init$P1 + (exp(init$V1) / V_den_ex5) * init$P5) +
    (exp(init$V6) / V_den_ex23) * ((exp(init$V3) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex3) * init$P3) +
    (exp(init$V6) / V_den_ex24) * ((exp(init$V4) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex4) * init$P4) +
    (exp(init$V6) / V_den_ex25) * ((exp(init$V5) / V_den_ex2) * init$P2 + (exp(init$V2) / V_den_ex5) * init$P5) +
    (exp(init$V6) / V_den_ex34) * ((exp(init$V4) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex4) * init$P4) +
    (exp(init$V6) / V_den_ex35) * ((exp(init$V5) / V_den_ex3) * init$P3 + (exp(init$V3) / V_den_ex5) * init$P5) +
    (exp(init$V6) / V_den_ex45) * ((exp(init$V5) / V_den_ex4) * init$P4 + (exp(init$V4) / V_den_ex5) * init$P5)

  # # CHECK: mean should be 1
  # mean(init$P1_t > 0 & init$P1_t < 1)
  # mean(init$P2_t > 0 & init$P2_t < 1)
  # mean(init$P3_t > 0 & init$P3_t < 1)
  # mean(init$P4_t > 0 & init$P4_t < 1)
  # mean(init$P5_t > 0 & init$P5_t < 1)
  # mean(init$P6_t > 0 & init$P6_t < 1)
  #
  # # Check that maximum probability of choosing candidate 1 1st and choosing candidate 1 second
  # # and choosing candidate 1 third does not exceed 1
  #
  # max(init$P1 + init$P1_s + init$P1_t)
  # max(init$P2 + init$P2_s + init$P2_t)
  # max(init$P3 + init$P3_s + init$P3_t)
  # max(init$P4 + init$P4_s + init$P4_t)
  # max(init$P5 + init$P5_s + init$P5_t)
  # max(init$P6 + init$P6_s + init$P6_t)

  # Compute Ethnic Polarization in initial third choices
  init$epsilon3 <- case_when(
    init$ethnic_group == "Group 1" ~ init$P1_t + init$P4_t,
    init$ethnic_group == "Group 2" ~ init$P2_t + init$P5_t,
    init$ethnic_group == "Group 3" ~ init$P3_t + init$P6_t
  )


  # Compute the THIRD-choice probability score (sum of all support probabilities)
  P1_tt_score <- NA
  P2_tt_score <- NA
  P3_tt_score <- NA
  P4_tt_score <- NA
  P5_tt_score <- NA
  P6_tt_score <- NA
  P1_tt_score[1] <- sum(init$P1_t)
  P2_tt_score[1] <- sum(init$P2_t)
  P3_tt_score[1] <- sum(init$P3_t)
  P4_tt_score[1] <- sum(init$P4_t)
  P5_tt_score[1] <- sum(init$P5_t)
  P6_tt_score[1] <- sum(init$P6_t)

  P_tt_vec <- NA
  P_tt_vec <- rbind(P1_tt_score, P2_tt_score, P3_tt_score, P4_tt_score, P5_tt_score, P6_tt_score)

  # # CHECK
  # print(sum(P_tt_vec)) # GOOD! It needs to be N=1000

  # ========================================================================#
  # Updating party positions (initial move)
  move <- d_cands
  theta <- runif(n = 6, min = 0, max = 360) # Generating random angels

  move$x <- d_cands$x + cos(theta / 180 * pi) * (unit / 2) # Compute x-value given theta SETTING THIS 0.01 (8/3/2022)
  move$y <- d_cands$y + sin(theta / 180 * pi) * (unit / 2) # Compute y-value given theta
  move$dist <- sqrt((move$x - 0)^2 + (move$y - 0)^2)
  move$moderation <- move$dist / d_cands$dist

  move_rcv <- move # Using the same initial conditions
  move_rcv_t <- move # Using the same initial conditions
  theta_rcv <- theta # Using the same initial directions
  theta_rcv_t <- theta # Using the same initial directions
  # ========================================================================#

  iter <- NA
  track <- track_rcv <- track_rcv_t <- list()
  voter_max1 <- voter_max2 <- voter_max3 <- list()
  new_theta <- NA
  new_theta_rcv <- NA
  new_theta_rcv_t <- NA
  delta <- delta_md <- delta_ex <- NA
  delta_rcv <- delta_rcv_md <- delta_rcv_ex <- NA
  delta_rcv_t <- delta_rcv_t_md <- delta_rcv_t_ex <- NA
  dist_Am <- dist_Ae <- dist_Bm <- dist_Be <- dist_Cm <- dist_Ce <- NA
  dist_rcvAm <- dist_rcvAe <- dist_rcvBm <- dist_rcvBe <- dist_rcvCm <- dist_rcvCe <- NA
  dist_rcv_tAm <- dist_rcv_tAe <- dist_rcv_tBm <- dist_rcv_tBe <- dist_rcv_tCm <- dist_rcv_tCe <- NA
  epsilon <- NA
  epsilon_rcv_1 <- NA
  epsilon_rcv_2 <- NA
  epsilon_rcv_t_1 <- NA
  epsilon_rcv_t_2 <- NA
  epsilon_rcv_t_3 <- NA
  check_A <- check_B <- check_C <- NA # To Check if Extreme Parties stay extreme
  check_rcvA <- check_rcvB <- check_rcvC <- NA # To Check if Extreme Parties stay extreme
  check_rcv_tA <- check_rcv_tB <- check_rcv_tC <- NA # To Check if Extreme Parties stay extreme


  pb <- progress_bar$new(
    format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    total = n_iter,
    complete = "=", # Completion bar character
    incomplete = "-", # Incomplete bar character
    current = ">", # Current bar character
    clear = FALSE, # If TRUE, clears the bar when finish
    width = 100
  ) # Width of the progress bar


  # BEGINNING OF t loop ########################################################################
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

    chain$D1 <- sqrt((chain$x - move$x[1])^2 + (chain$y - move$y[1])^2)
    chain$D2 <- sqrt((chain$x - move$x[2])^2 + (chain$y - move$y[2])^2)
    chain$D3 <- sqrt((chain$x - move$x[3])^2 + (chain$y - move$y[3])^2)
    chain$D4 <- sqrt((chain$x - move$x[4])^2 + (chain$y - move$y[4])^2)
    chain$D5 <- sqrt((chain$x - move$x[5])^2 + (chain$y - move$y[5])^2)
    chain$D6 <- sqrt((chain$x - move$x[6])^2 + (chain$y - move$y[6])^2)

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

    # Compute Ethnic Polarization in first choices
    chain$epsilon1 <- case_when(
      chain$ethnic_group == "Group 1" ~ chain$P1 + chain$P4,
      chain$ethnic_group == "Group 2" ~ chain$P2 + chain$P5,
      chain$ethnic_group == "Group 3" ~ chain$P3 + chain$P6
    )


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
    move$new_dist <- NA # Initialize
    for (i in 1:6) {
      # For Moderate Parties
      if (i < 4) {
        new_theta[i] <- ifelse(P_vec[i, t] >= P_vec[i, t - 1], # If new position has higher prob score
          theta[i], # Keep going
          runif(n = 1, min = theta[i] + 90, max = theta[i] + 270)
        ) # New position on the other side

        move$x[i] <- move$x[i] + cos(new_theta[i] / 180 * pi) * unit # Compute x-value given theta
        move$y[i] <- move$y[i] + sin(new_theta[i] / 180 * pi) * unit # Compute y-value given theta
        move$new_dist[i] <- sqrt((move$x[i] - 0)^2 + (move$y[i] - 0)^2) # New Distance
        move$moderation[i] <- move$new_dist[i] / move$dist[i] # New Location / Initial Location
        theta[i] <- new_theta[i] # Update for the next iteration


        # For Extreme Parties
      } else {
        new_theta[i] <- ifelse(P_vec[i, t] >= P_vec[i, t - 1], # If new position has higher prob score
          theta[i], # Keep going
          runif(n = 1, min = theta[i] + 90, max = theta[i] + 270)
        ) # New position on the other side

        # Updating party locations
        propose_x <- move$x[i] + cos(new_theta[i] / 180 * pi) * unit # Compute x-value given theta
        propose_y <- move$y[i] + sin(new_theta[i] / 180 * pi) * unit # Compute y-value given theta
        move$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

        # Change the direction until extreme parties become more extreme than moderate parties
        if (force == TRUE) {
          while_iter <- 1
          while (move$new_dist[i] < move$new_dist[i - 3] & while_iter <= while_max) {
            new_theta[i] <- runif(n = 1, min = 0, max = 360) # Anywhere is okay as long as extreme parties can get out of the trap
            propose_x <- move$x[i] + cos(new_theta[i] / 180 * pi) * unit * boost # Compute x-value given theta
            propose_y <- move$y[i] + sin(new_theta[i] / 180 * pi) * unit * boost # Compute y-value given theta
            move$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

            while_iter <- while_iter + 1
          } # Close while () loop

          # If extreme parties cannot find a way out, STAY (NO MOVE)
          if (while_iter == while_max) {
            propose_x <- move$x[i]
            propose_y <- move$y[i]
          }
        } # Closing Assumption 1 condition


        move$x[i] <- propose_x # Saving the accepted location
        move$y[i] <- propose_y # Saving the accepted location
        move$moderation[i] <- move$new_dist[i] / move$dist[i] # New Location / Initial Location
        theta[i] <- new_theta[i] # Update for the next iteration
      } # Close else{}
    } # Close for () loop

    try(if (sum(P_vec[, 1]) != N_voters) stop("FPTP: 1st-choice probabilities do not sum up to one"))


    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # (2): RCV -- Second Choice
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # Compute the voting probabilities again

    # Computing the Euclidean distance between parties and voters
    chain_rcv$D1 <- sqrt((chain_rcv$x - move_rcv$x[1])^2 + (chain_rcv$y - move_rcv$y[1])^2)
    chain_rcv$D2 <- sqrt((chain_rcv$x - move_rcv$x[2])^2 + (chain_rcv$y - move_rcv$y[2])^2)
    chain_rcv$D3 <- sqrt((chain_rcv$x - move_rcv$x[3])^2 + (chain_rcv$y - move_rcv$y[3])^2)
    chain_rcv$D4 <- sqrt((chain_rcv$x - move_rcv$x[4])^2 + (chain_rcv$y - move_rcv$y[4])^2)
    chain_rcv$D5 <- sqrt((chain_rcv$x - move_rcv$x[5])^2 + (chain_rcv$y - move_rcv$y[5])^2)
    chain_rcv$D6 <- sqrt((chain_rcv$x - move_rcv$x[6])^2 + (chain_rcv$y - move_rcv$y[6])^2)

    # Computing the observed utility for each party (This is where we specify utility functions)
    chain_rcv$V1 <- -1 * c * chain_rcv$D1 - b * m1 + eps1
    chain_rcv$V2 <- -1 * c * chain_rcv$D2 - b * m2 + eps2
    chain_rcv$V3 <- -1 * c * chain_rcv$D3 - b * m3 + eps3
    chain_rcv$V4 <- -1 * c * chain_rcv$D4 - b * m4 + eps4
    chain_rcv$V5 <- -1 * c * chain_rcv$D5 - b * m5 + eps5
    chain_rcv$V6 <- -1 * c * chain_rcv$D6 - b * m6 + eps6

    # Computing the probability that each voter votes for each party (this is fixed)
    den <- exp(chain_rcv$V1) + exp(chain_rcv$V2) + exp(chain_rcv$V3) + exp(chain_rcv$V4) + exp(chain_rcv$V5) + exp(chain_rcv$V6)

    chain_rcv$P1 <- exp(chain_rcv$V1) / den
    chain_rcv$P2 <- exp(chain_rcv$V2) / den
    chain_rcv$P3 <- exp(chain_rcv$V3) / den
    chain_rcv$P4 <- exp(chain_rcv$V4) / den
    chain_rcv$P5 <- exp(chain_rcv$V5) / den
    chain_rcv$P6 <- exp(chain_rcv$V6) / den

    # Compute Ethnic Polarization in first choices
    chain_rcv$epsilon1 <- case_when(
      chain_rcv$ethnic_group == "Group 1" ~ chain_rcv$P1 + chain_rcv$P4,
      chain_rcv$ethnic_group == "Group 2" ~ chain_rcv$P2 + chain_rcv$P5,
      chain_rcv$ethnic_group == "Group 3" ~ chain_rcv$P3 + chain_rcv$P6
    )




    V_den_ex1 <- exp(chain_rcv$V2) + exp(chain_rcv$V3) + exp(chain_rcv$V4) + exp(chain_rcv$V5) + exp(chain_rcv$V6) # Second-choice prob without P1
    V_den_ex2 <- exp(chain_rcv$V1) + exp(chain_rcv$V3) + exp(chain_rcv$V4) + exp(chain_rcv$V5) + exp(chain_rcv$V6) # Second-choice prob without P2
    V_den_ex3 <- exp(chain_rcv$V1) + exp(chain_rcv$V2) + exp(chain_rcv$V4) + exp(chain_rcv$V5) + exp(chain_rcv$V6) # Second-choice prob without P3
    V_den_ex4 <- exp(chain_rcv$V1) + exp(chain_rcv$V2) + exp(chain_rcv$V3) + exp(chain_rcv$V5) + exp(chain_rcv$V6) # Second-choice prob without P4
    V_den_ex5 <- exp(chain_rcv$V1) + exp(chain_rcv$V2) + exp(chain_rcv$V3) + exp(chain_rcv$V4) + exp(chain_rcv$V6) # Second-choice prob without P5
    V_den_ex6 <- exp(chain_rcv$V1) + exp(chain_rcv$V2) + exp(chain_rcv$V3) + exp(chain_rcv$V4) + exp(chain_rcv$V5) # Second-choice prob without P6

    # Average second-choice probabilities (the probability of ranking (A,B,...))
    chain_rcv$P1_s <- (exp(chain_rcv$V1) / V_den_ex2) * chain_rcv$P2 + (exp(chain_rcv$V1) / V_den_ex3) * chain_rcv$P3 +
      (exp(chain_rcv$V1) / V_den_ex4) * chain_rcv$P4 + (exp(chain_rcv$V1) / V_den_ex5) * chain_rcv$P5 + (exp(chain_rcv$V1) / V_den_ex6) * chain_rcv$P6
    chain_rcv$P2_s <- (exp(chain_rcv$V2) / V_den_ex1) * chain_rcv$P1 + (exp(chain_rcv$V2) / V_den_ex3) * chain_rcv$P3 +
      (exp(chain_rcv$V2) / V_den_ex4) * chain_rcv$P4 + (exp(chain_rcv$V2) / V_den_ex5) * chain_rcv$P5 + (exp(chain_rcv$V2) / V_den_ex6) * chain_rcv$P6
    chain_rcv$P3_s <- (exp(chain_rcv$V3) / V_den_ex1) * chain_rcv$P1 + (exp(chain_rcv$V3) / V_den_ex2) * chain_rcv$P2 +
      (exp(chain_rcv$V3) / V_den_ex4) * chain_rcv$P4 + (exp(chain_rcv$V3) / V_den_ex5) * chain_rcv$P5 + (exp(chain_rcv$V3) / V_den_ex6) * chain_rcv$P6
    chain_rcv$P4_s <- (exp(chain_rcv$V4) / V_den_ex1) * chain_rcv$P1 + (exp(chain_rcv$V4) / V_den_ex2) * chain_rcv$P2 +
      (exp(chain_rcv$V4) / V_den_ex3) * chain_rcv$P3 + (exp(chain_rcv$V4) / V_den_ex5) * chain_rcv$P5 + (exp(chain_rcv$V4) / V_den_ex6) * chain_rcv$P6
    chain_rcv$P5_s <- (exp(chain_rcv$V5) / V_den_ex1) * chain_rcv$P1 + (exp(chain_rcv$V5) / V_den_ex2) * chain_rcv$P2 +
      (exp(chain_rcv$V5) / V_den_ex3) * chain_rcv$P3 + (exp(chain_rcv$V5) / V_den_ex4) * chain_rcv$P4 + (exp(chain_rcv$V5) / V_den_ex6) * chain_rcv$P6
    chain_rcv$P6_s <- (exp(chain_rcv$V6) / V_den_ex1) * chain_rcv$P1 + (exp(chain_rcv$V6) / V_den_ex2) * chain_rcv$P2 +
      (exp(chain_rcv$V6) / V_den_ex3) * chain_rcv$P3 + (exp(chain_rcv$V6) / V_den_ex4) * chain_rcv$P4 + (exp(chain_rcv$V6) / V_den_ex5) * chain_rcv$P5

    # Compute Ethnic Polarization in second choices
    chain_rcv$epsilon2 <- case_when(
      chain_rcv$ethnic_group == "Group 1" ~ chain_rcv$P1_s + chain_rcv$P4_s,
      chain_rcv$ethnic_group == "Group 2" ~ chain_rcv$P2_s + chain_rcv$P5_s,
      chain_rcv$ethnic_group == "Group 3" ~ chain_rcv$P3_s + chain_rcv$P6_s
    )



    # Compute the probability score (sum of all support probabilities)
    P1_score_rcv[t] <- sum(chain_rcv$P1)
    P2_score_rcv[t] <- sum(chain_rcv$P2)
    P3_score_rcv[t] <- sum(chain_rcv$P3)
    P4_score_rcv[t] <- sum(chain_rcv$P4)
    P5_score_rcv[t] <- sum(chain_rcv$P5)
    P6_score_rcv[t] <- sum(chain_rcv$P6)

    # Compute the second-choice probability score (sum of all support probabilities)
    P1_s_score[t] <- sum(chain_rcv$P1_s)
    P2_s_score[t] <- sum(chain_rcv$P2_s)
    P3_s_score[t] <- sum(chain_rcv$P3_s)
    P4_s_score[t] <- sum(chain_rcv$P4_s)
    P5_s_score[t] <- sum(chain_rcv$P5_s)
    P6_s_score[t] <- sum(chain_rcv$P6_s)

    P_vec_rcv <- rbind(P1_score_rcv, P2_score_rcv, P3_score_rcv, P4_score_rcv, P5_score_rcv, P6_score_rcv)
    P_s_vec <- rbind(P1_s_score, P2_s_score, P3_s_score, P4_s_score, P5_s_score, P6_s_score)

    # CHECK
    try(if (sum(P_vec_rcv[, 1]) != N_voters) stop("RCV1: 1st-choice ranking probabilities do not sum up to one"))
    try(if (sum(P_s_vec[, 1]) != N_voters) stop("RCV1: 2nd-choice ranking probabilities do not sum up to one"))


    # ============================================================================#
    #  Updating party locations
    # ============================================================================#

    new_theta_rcv <- NA # Initialize
    move_rcv$new_dist <- NA # Initialize
    for (i in 1:6) {
      # For moderate parties
      if (i < 4) {
        new_theta_rcv[i] <- ifelse(P_vec_rcv[i, t] >= P_vec_rcv[i, t - 1] &
          P_s_vec[i, t] >= P_s_vec[i, t - 1], # If new position has higher 1st and 2nd choice probs
        theta_rcv[i], # Keep going
        runif(n = 1, min = theta_rcv[i] + 90, max = theta_rcv[i] + 270)
        ) # New position on the other side

        move_rcv$x[i] <- move_rcv$x[i] + cos(new_theta_rcv[i] / 180 * pi) * unit # Compute x-value given theta
        move_rcv$y[i] <- move_rcv$y[i] + sin(new_theta_rcv[i] / 180 * pi) * unit # Compute y-value given theta
        move_rcv$new_dist[i] <- sqrt((move_rcv$x[i] - 0)^2 + (move_rcv$y[i] - 0)^2) # New Distance
        move_rcv$moderation[i] <- move_rcv$new_dist[i] / move_rcv$dist[i] # New Location / Initial Location
        theta_rcv[i] <- new_theta_rcv[i] # Update for the next iteration

        # For extreme parties
      } else {
        new_theta_rcv[i] <- ifelse(P_vec_rcv[i, t] >= P_vec_rcv[i, t - 1] &
          P_s_vec[i, t] >= P_s_vec[i, t - 1], # Extreme parties stay "extreme"
        theta_rcv[i], # Keep going
        runif(n = 1, min = theta_rcv[i] + 90, max = theta_rcv[i] + 270)
        ) # New position on the other side

        # Updating party locations
        propose_x <- move_rcv$x[i] + cos(new_theta_rcv[i] / 180 * pi) * unit # Compute x-value given theta
        propose_y <- move_rcv$y[i] + sin(new_theta_rcv[i] / 180 * pi) * unit # Compute y-value given theta
        move_rcv$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

        if (force == TRUE) {
          while_iter <- 1
          # Change the direction until extreme parties become more extreme than moderate parties
          while (move_rcv$new_dist[i] < move_rcv$new_dist[i - 3] & while_iter <= while_max) {
            new_theta_rcv[i] <- runif(n = 1, min = 0, max = 360) # Anywhere is okay as long as extreme parties can get out of the trap
            propose_x <- move_rcv$x[i] + cos(new_theta_rcv[i] / 180 * pi) * unit * boost # Compute x-value given theta
            propose_y <- move_rcv$y[i] + sin(new_theta_rcv[i] / 180 * pi) * unit * boost # Compute y-value given theta
            move_rcv$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

            while_iter <- while_iter + 1
          } # Close while () loop
          # if (while_iter >= while_max) {
          #   print(paste0("While Iteration (RCV - 2nd): ", while_iter, " for party ", i, " at iteration = ", t))
          # }

          # If extreme parties cannot find a way out, STAY (NO MOVE)
          if (while_iter == while_max) {
            propose_x <- move$x[i]
            propose_y <- move$y[i]
          }
        } # Closing the Assumption 1 condition


        move_rcv$x[i] <- propose_x # Saving the accepted location
        move_rcv$y[i] <- propose_y # Saving the accepted location

        move_rcv$moderation[i] <- move_rcv$new_dist[i] / move_rcv$dist[i] # New Location / Initial Location
        theta_rcv[i] <- new_theta_rcv[i] # Update for the next iteration
      } # Close else{}
    } # Close for () loop


    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # (3): RCV -- Third Choice
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # Compute the voting probabilities again

    # Computing the Euclidean distance between parties and voters
    chain_rcv_t$D1 <- sqrt((chain_rcv_t$x - move_rcv_t$x[1])^2 + (chain_rcv_t$y - move_rcv_t$y[1])^2)
    chain_rcv_t$D2 <- sqrt((chain_rcv_t$x - move_rcv_t$x[2])^2 + (chain_rcv_t$y - move_rcv_t$y[2])^2)
    chain_rcv_t$D3 <- sqrt((chain_rcv_t$x - move_rcv_t$x[3])^2 + (chain_rcv_t$y - move_rcv_t$y[3])^2)
    chain_rcv_t$D4 <- sqrt((chain_rcv_t$x - move_rcv_t$x[4])^2 + (chain_rcv_t$y - move_rcv_t$y[4])^2)
    chain_rcv_t$D5 <- sqrt((chain_rcv_t$x - move_rcv_t$x[5])^2 + (chain_rcv_t$y - move_rcv_t$y[5])^2)
    chain_rcv_t$D6 <- sqrt((chain_rcv_t$x - move_rcv_t$x[6])^2 + (chain_rcv_t$y - move_rcv_t$y[6])^2)

    # Computing the observed utility for each party (This is where we specify utility functions)
    chain_rcv_t$V1 <- -1 * c * chain_rcv_t$D1 - b * m1 + eps1
    chain_rcv_t$V2 <- -1 * c * chain_rcv_t$D2 - b * m2 + eps2
    chain_rcv_t$V3 <- -1 * c * chain_rcv_t$D3 - b * m3 + eps3
    chain_rcv_t$V4 <- -1 * c * chain_rcv_t$D4 - b * m4 + eps4
    chain_rcv_t$V5 <- -1 * c * chain_rcv_t$D5 - b * m5 + eps5
    chain_rcv_t$V6 <- -1 * c * chain_rcv_t$D6 - b * m6 + eps6

    # Computing the probability that each voter votes for each party (this is fixed)
    den <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6)

    chain_rcv_t$P1 <- exp(chain_rcv_t$V1) / den
    chain_rcv_t$P2 <- exp(chain_rcv_t$V2) / den
    chain_rcv_t$P3 <- exp(chain_rcv_t$V3) / den
    chain_rcv_t$P4 <- exp(chain_rcv_t$V4) / den
    chain_rcv_t$P5 <- exp(chain_rcv_t$V5) / den
    chain_rcv_t$P6 <- exp(chain_rcv_t$V6) / den

    # Compute Ethnic Polarization in first choices (RCV 3 choice model)
    chain_rcv_t$epsilon1 <- case_when(
      chain_rcv_t$ethnic_group == "Group 1" ~ chain_rcv_t$P1 + chain_rcv_t$P4,
      chain_rcv_t$ethnic_group == "Group 2" ~ chain_rcv_t$P2 + chain_rcv_t$P5,
      chain_rcv_t$ethnic_group == "Group 3" ~ chain_rcv_t$P3 + chain_rcv_t$P6
    )


    V_den_ex1 <- exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Second-choice prob without P1
    V_den_ex2 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Second-choice prob without P2
    V_den_ex3 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Second-choice prob without P3
    V_den_ex4 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Second-choice prob without P4
    V_den_ex5 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V6) # Second-choice prob without P5
    V_den_ex6 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) # Second-choice prob without P6

    # Average second-choice probabilities (the probability of ranking (A,B,...))
    chain_rcv_t$P1_s <- (exp(chain_rcv_t$V1) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V1) / V_den_ex3) * chain_rcv_t$P3 +
      (exp(chain_rcv_t$V1) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V1) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V1) / V_den_ex6) * chain_rcv_t$P6
    chain_rcv_t$P2_s <- (exp(chain_rcv_t$V2) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V2) / V_den_ex3) * chain_rcv_t$P3 +
      (exp(chain_rcv_t$V2) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V2) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V2) / V_den_ex6) * chain_rcv_t$P6
    chain_rcv_t$P3_s <- (exp(chain_rcv_t$V3) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V3) / V_den_ex2) * chain_rcv_t$P2 +
      (exp(chain_rcv_t$V3) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V3) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V3) / V_den_ex6) * chain_rcv_t$P6
    chain_rcv_t$P4_s <- (exp(chain_rcv_t$V4) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V4) / V_den_ex2) * chain_rcv_t$P2 +
      (exp(chain_rcv_t$V4) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V4) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V4) / V_den_ex6) * chain_rcv_t$P6
    chain_rcv_t$P5_s <- (exp(chain_rcv_t$V5) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V5) / V_den_ex2) * chain_rcv_t$P2 +
      (exp(chain_rcv_t$V5) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V5) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V5) / V_den_ex6) * chain_rcv_t$P6
    chain_rcv_t$P6_s <- (exp(chain_rcv_t$V6) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V6) / V_den_ex2) * chain_rcv_t$P2 +
      (exp(chain_rcv_t$V6) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V6) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V6) / V_den_ex5) * chain_rcv_t$P5

    # Compute Ethnic Polarization in second choices (RCV 3 choice model)
    chain_rcv_t$epsilon2 <- case_when(
      chain_rcv_t$ethnic_group == "Group 1" ~ chain_rcv_t$P1_s + chain_rcv_t$P4_s,
      chain_rcv_t$ethnic_group == "Group 2" ~ chain_rcv_t$P2_s + chain_rcv_t$P5_s,
      chain_rcv_t$ethnic_group == "Group 3" ~ chain_rcv_t$P3_s + chain_rcv_t$P6_s
    )


    # Computing the third-choice probability that each voter votes for each party
    V_den_ex12 <- exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Third-choice prob without P1 & P2
    V_den_ex13 <- exp(chain_rcv_t$V2) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Third-choice prob without P1 & P3
    V_den_ex14 <- exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Third-choice prob without P1 & P4
    V_den_ex15 <- exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V6) # Third-choice prob without P1 & P5
    V_den_ex16 <- exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) # Third-choice prob without P1 & P6
    V_den_ex23 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Third-choice prob without P2 & P3
    V_den_ex24 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Third-choice prob without P2 & P4
    V_den_ex25 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V6) # Third-choice prob without P2 & P5
    V_den_ex26 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) # Third-choice prob without P2 & P6
    V_den_ex34 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V5) + exp(chain_rcv_t$V6) # Third-choice prob without P3 & P4
    V_den_ex35 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V6) # Third-choice prob without P3 & P5
    V_den_ex36 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V4) + exp(chain_rcv_t$V5) # Third-choice prob without P3 & P6
    V_den_ex45 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V6) # Third-choice prob without P4 & P5
    V_den_ex46 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V5) # Third-choice prob without P4 & P6
    V_den_ex56 <- exp(chain_rcv_t$V1) + exp(chain_rcv_t$V2) + exp(chain_rcv_t$V3) + exp(chain_rcv_t$V4) # Third-choice prob without P5 & P6

    # Weighted Average third-choice probabilities
    # exp(party A)/(sum exp(remaining)) * [prob(second|first)prob(first) + prob(second|first)prob(first)]
    chain_rcv_t$P1_t <- (exp(chain_rcv_t$V1) / V_den_ex23) * ((exp(chain_rcv_t$V3) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex3) * chain_rcv_t$P3) + # Prob of choosing 2 then 3 + Prob of choosing 3 then 2
      (exp(chain_rcv_t$V1) / V_den_ex24) * ((exp(chain_rcv_t$V4) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V1) / V_den_ex25) * ((exp(chain_rcv_t$V5) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V1) / V_den_ex26) * ((exp(chain_rcv_t$V6) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V1) / V_den_ex34) * ((exp(chain_rcv_t$V4) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V1) / V_den_ex35) * ((exp(chain_rcv_t$V5) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V1) / V_den_ex36) * ((exp(chain_rcv_t$V6) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V1) / V_den_ex45) * ((exp(chain_rcv_t$V5) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V1) / V_den_ex46) * ((exp(chain_rcv_t$V6) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V1) / V_den_ex56) * ((exp(chain_rcv_t$V6) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V5) / V_den_ex6) * chain_rcv_t$P6)

    chain_rcv_t$P2_t <- (exp(chain_rcv_t$V2) / V_den_ex13) * ((exp(chain_rcv_t$V3) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex3) * chain_rcv_t$P3) + #
      (exp(chain_rcv_t$V2) / V_den_ex14) * ((exp(chain_rcv_t$V4) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V2) / V_den_ex15) * ((exp(chain_rcv_t$V5) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V2) / V_den_ex16) * ((exp(chain_rcv_t$V6) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V2) / V_den_ex34) * ((exp(chain_rcv_t$V4) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V2) / V_den_ex35) * ((exp(chain_rcv_t$V5) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V2) / V_den_ex36) * ((exp(chain_rcv_t$V6) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V2) / V_den_ex45) * ((exp(chain_rcv_t$V5) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V2) / V_den_ex46) * ((exp(chain_rcv_t$V6) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V2) / V_den_ex56) * ((exp(chain_rcv_t$V6) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V5) / V_den_ex6) * chain_rcv_t$P6)

    chain_rcv_t$P3_t <- (exp(chain_rcv_t$V3) / V_den_ex12) * ((exp(chain_rcv_t$V2) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex2) * chain_rcv_t$P2) + #
      (exp(chain_rcv_t$V3) / V_den_ex14) * ((exp(chain_rcv_t$V4) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V3) / V_den_ex15) * ((exp(chain_rcv_t$V5) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V3) / V_den_ex16) * ((exp(chain_rcv_t$V6) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V3) / V_den_ex24) * ((exp(chain_rcv_t$V4) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V3) / V_den_ex25) * ((exp(chain_rcv_t$V5) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V3) / V_den_ex26) * ((exp(chain_rcv_t$V6) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V3) / V_den_ex45) * ((exp(chain_rcv_t$V5) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V3) / V_den_ex46) * ((exp(chain_rcv_t$V6) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V3) / V_den_ex56) * ((exp(chain_rcv_t$V6) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V5) / V_den_ex6) * chain_rcv_t$P6)

    chain_rcv_t$P4_t <- (exp(chain_rcv_t$V4) / V_den_ex12) * ((exp(chain_rcv_t$V2) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex2) * chain_rcv_t$P2) + #
      (exp(chain_rcv_t$V4) / V_den_ex13) * ((exp(chain_rcv_t$V3) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex3) * chain_rcv_t$P3) +
      (exp(chain_rcv_t$V4) / V_den_ex15) * ((exp(chain_rcv_t$V5) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V4) / V_den_ex16) * ((exp(chain_rcv_t$V6) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V4) / V_den_ex23) * ((exp(chain_rcv_t$V3) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex3) * chain_rcv_t$P3) +
      (exp(chain_rcv_t$V4) / V_den_ex25) * ((exp(chain_rcv_t$V5) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V4) / V_den_ex26) * ((exp(chain_rcv_t$V6) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V4) / V_den_ex35) * ((exp(chain_rcv_t$V5) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V4) / V_den_ex36) * ((exp(chain_rcv_t$V6) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V4) / V_den_ex56) * ((exp(chain_rcv_t$V6) / V_den_ex5) * chain_rcv_t$P5 + (exp(chain_rcv_t$V5) / V_den_ex6) * chain_rcv_t$P6)

    chain_rcv_t$P5_t <- (exp(chain_rcv_t$V5) / V_den_ex12) * ((exp(chain_rcv_t$V2) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex2) * chain_rcv_t$P2) + #
      (exp(chain_rcv_t$V5) / V_den_ex13) * ((exp(chain_rcv_t$V3) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex3) * chain_rcv_t$P3) +
      (exp(chain_rcv_t$V5) / V_den_ex14) * ((exp(chain_rcv_t$V4) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V5) / V_den_ex16) * ((exp(chain_rcv_t$V6) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V5) / V_den_ex23) * ((exp(chain_rcv_t$V3) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex3) * chain_rcv_t$P3) +
      (exp(chain_rcv_t$V5) / V_den_ex24) * ((exp(chain_rcv_t$V4) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V5) / V_den_ex26) * ((exp(chain_rcv_t$V6) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V5) / V_den_ex34) * ((exp(chain_rcv_t$V4) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V5) / V_den_ex36) * ((exp(chain_rcv_t$V6) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex6) * chain_rcv_t$P6) +
      (exp(chain_rcv_t$V5) / V_den_ex46) * ((exp(chain_rcv_t$V6) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex6) * chain_rcv_t$P6)

    chain_rcv_t$P6_t <- (exp(chain_rcv_t$V6) / V_den_ex12) * ((exp(chain_rcv_t$V2) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex2) * chain_rcv_t$P2) +
      (exp(chain_rcv_t$V6) / V_den_ex13) * ((exp(chain_rcv_t$V3) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex3) * chain_rcv_t$P3) +
      (exp(chain_rcv_t$V6) / V_den_ex14) * ((exp(chain_rcv_t$V4) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V6) / V_den_ex15) * ((exp(chain_rcv_t$V5) / V_den_ex1) * chain_rcv_t$P1 + (exp(chain_rcv_t$V1) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V6) / V_den_ex23) * ((exp(chain_rcv_t$V3) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex3) * chain_rcv_t$P3) +
      (exp(chain_rcv_t$V6) / V_den_ex24) * ((exp(chain_rcv_t$V4) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V6) / V_den_ex25) * ((exp(chain_rcv_t$V5) / V_den_ex2) * chain_rcv_t$P2 + (exp(chain_rcv_t$V2) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V6) / V_den_ex34) * ((exp(chain_rcv_t$V4) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex4) * chain_rcv_t$P4) +
      (exp(chain_rcv_t$V6) / V_den_ex35) * ((exp(chain_rcv_t$V5) / V_den_ex3) * chain_rcv_t$P3 + (exp(chain_rcv_t$V3) / V_den_ex5) * chain_rcv_t$P5) +
      (exp(chain_rcv_t$V6) / V_den_ex45) * ((exp(chain_rcv_t$V5) / V_den_ex4) * chain_rcv_t$P4 + (exp(chain_rcv_t$V4) / V_den_ex5) * chain_rcv_t$P5)


    # Compute Ethnic Polarization in first choices (RCV 3 choice model)
    chain_rcv_t$epsilon3 <- case_when(
      chain_rcv_t$ethnic_group == "Group 1" ~ chain_rcv_t$P1_t + chain_rcv_t$P4_t,
      chain_rcv_t$ethnic_group == "Group 2" ~ chain_rcv_t$P2_t + chain_rcv_t$P5_t,
      chain_rcv_t$ethnic_group == "Group 3" ~ chain_rcv_t$P3_t + chain_rcv_t$P6_t
    )


    # Compute the probability score (sum of all support probabilities)
    P1_score_rcv_t[t] <- sum(chain_rcv_t$P1)
    P2_score_rcv_t[t] <- sum(chain_rcv_t$P2)
    P3_score_rcv_t[t] <- sum(chain_rcv_t$P3)
    P4_score_rcv_t[t] <- sum(chain_rcv_t$P4)
    P5_score_rcv_t[t] <- sum(chain_rcv_t$P5)
    P6_score_rcv_t[t] <- sum(chain_rcv_t$P6)

    # Compute the second-choice probability score (sum of all support probabilities)
    P1_st_score[t] <- sum(chain_rcv_t$P1_s)
    P2_st_score[t] <- sum(chain_rcv_t$P2_s)
    P3_st_score[t] <- sum(chain_rcv_t$P3_s)
    P4_st_score[t] <- sum(chain_rcv_t$P4_s)
    P5_st_score[t] <- sum(chain_rcv_t$P5_s)
    P6_st_score[t] <- sum(chain_rcv_t$P6_s)

    # Compute the third-choice probability score (sum of all support probabilities)
    P1_tt_score[t] <- sum(chain_rcv_t$P1_t)
    P2_tt_score[t] <- sum(chain_rcv_t$P2_t)
    P3_tt_score[t] <- sum(chain_rcv_t$P3_t)
    P4_tt_score[t] <- sum(chain_rcv_t$P4_t)
    P6_tt_score[t] <- sum(chain_rcv_t$P6_t)
    P5_tt_score[t] <- sum(chain_rcv_t$P5_t)

    P_vec_rcv_t <- rbind(P1_score_rcv_t, P2_score_rcv_t, P3_score_rcv_t, P4_score_rcv_t, P5_score_rcv_t, P6_score_rcv_t)
    P_st_vec <- rbind(P1_st_score, P2_st_score, P3_st_score, P4_st_score, P5_st_score, P6_st_score)
    P_tt_vec <- rbind(P1_tt_score, P2_tt_score, P3_tt_score, P4_tt_score, P5_tt_score, P6_tt_score)

    # CHECK
    try(if (sum(P_vec_rcv_t[, 1]) != N_voters) stop("RCV3: 1st-choice ranking probabilities do not sum up to one"))
    try(if (sum(P_st_vec[, 1]) != N_voters) stop("RCV3: 2nd-choice ranking probabilities do not sum up to one"))
    try(if (sum(P_tt_vec[, 1]) != N_voters) stop("RCV3: 3rd-choice ranking probabilities do not sum up to one"))

    # ============================================================================#
    #  Updating party locations
    # ============================================================================#

    new_theta_rcv_t <- NA # Initialize
    move_rcv_t$new_dist <- NA # Initialize
    for (i in 1:6) {
      # For moderate parties
      if (i < 4) {
        new_theta_rcv_t[i] <- ifelse(P_vec_rcv_t[i, t] >= P_vec_rcv_t[i, t - 1] &
          P_st_vec[i, t] >= P_st_vec[i, t - 1] &
          P_tt_vec[i, t] >= P_tt_vec[i, t - 1], # If new position has higher 1st and 2nd choice probs
        theta_rcv_t[i], # Keep going
        runif(n = 1, min = theta_rcv_t[i] + 90, max = theta_rcv_t[i] + 270)
        ) # New position on the other side

        move_rcv_t$x[i] <- move_rcv_t$x[i] + cos(new_theta_rcv_t[i] / 180 * pi) * unit # Compute x-value given theta
        move_rcv_t$y[i] <- move_rcv_t$y[i] + sin(new_theta_rcv_t[i] / 180 * pi) * unit # Compute y-value given theta
        move_rcv_t$new_dist[i] <- sqrt((move_rcv_t$x[i] - 0)^2 + (move_rcv_t$y[i] - 0)^2) # New Distance
        move_rcv_t$moderation[i] <- move_rcv_t$new_dist[i] / move_rcv_t$dist[i] # New Location / Initial Location
        theta_rcv_t[i] <- new_theta_rcv_t[i] # Update for the next iteration

        # For extreme parties
      } else {
        new_theta_rcv_t[i] <- ifelse(P_vec_rcv_t[i, t] >= P_vec_rcv_t[i, t - 1] &
          P_st_vec[i, t] >= P_st_vec[i, t - 1] &
          P_tt_vec[i, t] >= P_tt_vec[i, t - 1], # Extreme parties stay "extreme"
        theta_rcv_t[i], # Keep going
        runif(n = 1, min = theta_rcv_t[i] + 90, max = theta_rcv_t[i] + 270)
        ) # New position on the other side

        # Updating party locations
        propose_x <- move_rcv_t$x[i] + cos(new_theta_rcv_t[i] / 180 * pi) * unit # Compute x-value given theta
        propose_y <- move_rcv_t$y[i] + sin(new_theta_rcv_t[i] / 180 * pi) * unit # Compute y-value given theta
        move_rcv_t$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

        if (force == TRUE) {
          while_iter <- 1
          # Change the direction until extreme parties become more extreme than moderate parties
          while (move_rcv_t$new_dist[i] < move_rcv_t$new_dist[i - 3] & while_iter <= while_max) {
            new_theta_rcv_t[i] <- runif(n = 1, min = 0, max = 360) # Anywhere is okay as long as extreme parties can get out of the trap
            propose_x <- move_rcv_t$x[i] + cos(new_theta_rcv[i] / 180 * pi) * unit * boost # Extra Bump (unit x 5) (8/3/2022)
            propose_y <- move_rcv_t$y[i] + sin(new_theta_rcv[i] / 180 * pi) * unit * boost # Extra Bump (unit x 5) (8/3/2022)
            move_rcv_t$new_dist[i] <- sqrt((propose_x - 0)^2 + (propose_y - 0)^2) # New Distance

            while_iter <- while_iter + 1
          } # Close while () loop
          # if (while_iter >= while_max) {
          #   print(paste0("While Iteration (RCV -3rd): ", while_iter, " for party ", i, " at iteration = ", t))
          # }

          # If extreme parties cannot find a way out, STAY
          if (while_iter == while_max) {
            propose_x <- move$x[i]
            propose_y <- move$y[i]
          }
        } # Closing the Assumption 1 condition

        move_rcv_t$x[i] <- propose_x # Saving the accepted location
        move_rcv_t$y[i] <- propose_y # Saving the accepted location

        move_rcv_t$moderation[i] <- move_rcv_t$new_dist[i] / move_rcv_t$dist[i] # New Location / Initial Location
        theta_rcv_t[i] <- new_theta_rcv_t[i] # Update for the next iteration
      } # Close else{}
    } # Close for () loop


    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #
    # Keep Track of Party Positions
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX #

    move$iter <- t
    move_rcv$iter <- t
    move_rcv_t$iter <- t

    # Save candidate information
    track[[t]] <- move
    track_rcv[[t]] <- move_rcv
    track_rcv_t[[t]] <- move_rcv_t

    # Save voter information
    chain$iter <- t
    chain_rcv$iter <- t
    chain_rcv_t$iter <- t
    voter_max1[[t]] <- chain
    voter_max2[[t]] <- chain_rcv
    voter_max3[[t]] <- chain_rcv_t

    # Ideological polarization
    delta[t] <- mean(move$moderation)
    delta_rcv[t] <- mean(move_rcv$moderation)
    delta_rcv_t[t] <- mean(move_rcv_t$moderation)

    ## By parties
    delta_md[t] <- mean(move$moderation[move$type == "mod"])
    delta_rcv_md[t] <- mean(move_rcv$moderation[move$type == "mod"])
    delta_rcv_t_md[t] <- mean(move_rcv_t$moderation[move$type == "mod"])
    delta_ex[t] <- mean(move$moderation[move$type == "ext"])
    delta_rcv_ex[t] <- mean(move_rcv$moderation[move$type == "ext"])
    delta_rcv_t_ex[t] <- mean(move_rcv_t$moderation[move$type == "ext"])


    # Ethnic polarization
    epsilon[t] <- mean(chain$epsilon1) / mean(init$epsilon1)
    epsilon_rcv_1[t] <- mean(chain_rcv$epsilon1) / mean(init$epsilon1)
    epsilon_rcv_2[t] <- mean(chain_rcv$epsilon2) / mean(init$epsilon2)
    epsilon_rcv_t_1[t] <- mean(chain_rcv_t$epsilon1) / mean(init$epsilon1)
    epsilon_rcv_t_2[t] <- mean(chain_rcv_t$epsilon2) / mean(init$epsilon2)
    epsilon_rcv_t_3[t] <- mean(chain_rcv_t$epsilon3) / mean(init$epsilon3)

    iter[t] <- t

    # For keeping track of each patty's location
    dist_Am[t] <- move$new_dist[1]
    dist_Ae[t] <- move$new_dist[1 + 3]
    dist_Bm[t] <- move$new_dist[2]
    dist_Be[t] <- move$new_dist[2 + 3]
    dist_Cm[t] <- move$new_dist[3]
    dist_Ce[t] <- move$new_dist[3 + 3]

    dist_rcvAm[t] <- move_rcv$new_dist[1]
    dist_rcvAe[t] <- move_rcv$new_dist[1 + 3]
    dist_rcvBm[t] <- move_rcv$new_dist[2]
    dist_rcvBe[t] <- move_rcv$new_dist[2 + 3]
    dist_rcvCm[t] <- move_rcv$new_dist[3]
    dist_rcvCe[t] <- move_rcv$new_dist[3 + 3]

    dist_rcv_tAm[t] <- move_rcv_t$new_dist[1]
    dist_rcv_tAe[t] <- move_rcv_t$new_dist[1 + 3]
    dist_rcv_tBm[t] <- move_rcv_t$new_dist[2]
    dist_rcv_tBe[t] <- move_rcv_t$new_dist[2 + 3]
    dist_rcv_tCm[t] <- move_rcv_t$new_dist[3]
    dist_rcv_tCe[t] <- move_rcv_t$new_dist[3 + 3]

    # For checking the violetion of Assumption 1
    check_A[t] <- move$new_dist[4] > move$new_dist[1]
    check_B[t] <- move$new_dist[5] > move$new_dist[2]
    check_C[t] <- move$new_dist[6] > move$new_dist[3]

    check_rcvA[t] <- move_rcv$new_dist[4] > move_rcv$new_dist[1]
    check_rcvB[t] <- move_rcv$new_dist[5] > move_rcv$new_dist[2]
    check_rcvC[t] <- move_rcv$new_dist[6] > move_rcv$new_dist[3]

    check_rcv_tA[t] <- move_rcv_t$new_dist[4] > move_rcv_t$new_dist[1]
    check_rcv_tB[t] <- move_rcv_t$new_dist[5] > move_rcv_t$new_dist[2]
    check_rcv_tC[t] <- move_rcv_t$new_dist[6] > move_rcv_t$new_dist[3]
  }
  # END OF t loop #############################################################################
  # ed <- Sys.time()
  # print(ed - st) # Total time to execute the loop

  # Combining the summary statistics

  # out <- tibble(
  #   iter,
  #   delta, epsilon, delta_md, delta_ex,
  #   delta_rcv, epsilon_rcv_1, epsilon_rcv_2, delta_rcv_md, delta_rcv_ex,
  #   delta_rcv_t, epsilon_rcv_t_1, epsilon_rcv_t_2, epsilon_rcv_t_3, delta_rcv_t_md, delta_rcv_t_ex,
  #   # dist_Am, dist_Ae, dist_Bm, dist_Be, dist_Cm, dist_Ce,
  #   # dist_rcvAm, dist_rcvAe, dist_rcvBm, dist_rcvBe, dist_rcvCm, dist_rcvCe,
  #   # dist_rcv_tAm, dist_rcv_tAe, dist_rcv_tBm, dist_rcv_tBe, dist_rcv_tCm, dist_rcv_tCe
  # ) %>%
  #   dplyr::filter(iter != 1)

  track[[1]] <- d_cands %>%
    mutate(moderation = 1,
           new_dist = dist,
           iter = 1)

  track_rcv[[1]] <- d_cands %>%
    mutate(moderation = 1,
           new_dist = dist,
           iter = 1)

  track_rcv_t[[1]] <- d_cands %>%
    mutate(moderation = 1,
           new_dist = dist,
           iter = 1)

  # voter_max1[[1]] <- init %>%
  #   mutate(moderation = 1,
  #          new_dist = dist,
  #          iter = 1)

  track_max1 <- as.data.frame(do.call(rbind, track)) %>%
    mutate(system = "max1")
  track_max2 <- as.data.frame(do.call(rbind, track_rcv)) %>%
    mutate(system = "max2")
  track_max3 <- as.data.frame(do.call(rbind, track_rcv_t)) %>%
    mutate(system = "max3")

  # combine all results
  track <- rbind(track_max1, track_max2, track_max3) %>%
    tibble()


  voter_max1 <- as.data.frame(do.call(rbind, voter_max1)) %>%
    mutate(system = "max1")
  voter_max2 <- as.data.frame(do.call(rbind, voter_max2)) %>%
    mutate(system = "max2")
  voter_max3 <- as.data.frame(do.call(rbind, voter_max3)) %>%
    mutate(system = "max3")

  # combine all results
  voter_chains <- dplyr::bind_rows(voter_max1, voter_max2, voter_max3) %>%
    tibble()

  out <- list(
    voters = voter_chains,
    cands = track
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


