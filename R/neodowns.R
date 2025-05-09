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
                     strategy = "max1",
                     n_iter = 100,
                     mu_c = 1, # average spatial voting
                     mu_b = 2, # average co-ethnic voting
                     sigma_c = 0.5,
                     sigma_b = 0.5,
                     eps_sd = 0.5,
                     unit = 0.05,
                     # force = TRUE,
                     # while_max = 200,
                     # boost = 6,
                     seed = 14231) {

# Compute voter utility and candidate expected vote share ----------------------

  get_util_max1 <- function(d_voters, d_cands, N_voters, eps_sd, m_vec) {
    J <- nrow(d_cands)
    V_mat <- matrix(NA_real_, nrow = N_voters, ncol = J)
    for (j in seq_len(J)) {
      D_j <- sqrt((d_voters$x - d_cands$x[j])^2 + (d_voters$y - d_cands$y[j])^2)
      m_j <- as.integer(d_voters$ethnic_group != paste0("Group ", m_vec[j]))
      eps_j <- rnorm(N_voters, sd = eps_sd)
      V_mat[, j] <- -d_voters$c * D_j - d_voters$b * m_j + eps_j
    }
    expV <- exp(V_mat)
    P_mat <- expV / rowSums(expV)
    colnames(P_mat) <- paste0("P", seq_len(J))
    d_voters <- dplyr::bind_cols(d_voters, as.data.frame(V_mat), as.data.frame(P_mat))
    P_vec <- colSums(P_mat)
    list(d_voters = d_voters, P_vec = P_vec)
  }


  update_max1 <- function(chain, d_cands, theta, unit = 0.05, p_before, p_now) {
    J <- nrow(d_cands)
    new_theta <- numeric(J)
    improve <- p_now >= p_before
    new_theta[improve] <- theta[improve]
    new_theta[!improve] <- runif(sum(!improve), min = theta[!improve] + 90, max = theta[!improve] + 270)
    rad_theta <- new_theta * pi / 180
    dx <- cos(rad_theta) * unit
    dy <- sin(rad_theta) * unit
    d_cands$x <- d_cands$x + dx
    d_cands$y <- d_cands$y + dy
    list(d_cands = d_cands, theta = new_theta)
  }


  get_util_max2 <- function(d_voters, d_cands, N_voters, eps_sd, m_vec) {
    J <- nrow(d_cands)
    V_mat <- matrix(NA_real_, nrow = N_voters, ncol = J)

    for (j in seq_len(J)) {
      D_j <- sqrt((d_voters$x - d_cands$x[j])^2 + (d_voters$y - d_cands$y[j])^2)
      m_j <- as.integer(d_voters$ethnic_group != paste0("Group ", m_vec[j]))
      eps_j <- rnorm(N_voters, sd = eps_sd)
      V_mat[, j] <- -d_voters$c * D_j - d_voters$b * m_j + eps_j
    }

    expV <- exp(V_mat)
    P_mat <- expV / rowSums(expV)
    colnames(P_mat) <- paste0("P", seq_len(J))

    P_s_mat <- matrix(0, nrow = N_voters, ncol = J)
    for (j in seq_len(J)) {
      for (k in setdiff(seq_len(J), j)) {
        V_den_ex_k <- rowSums(expV[, setdiff(seq_len(J), k), drop = FALSE])
        contrib <- (expV[, j] / V_den_ex_k) * P_mat[, k]
        P_s_mat[, j] <- P_s_mat[, j] + contrib
      }
    }

    colnames(P_s_mat) <- paste0("P", seq_len(J), "_s")
    d_voters <- dplyr::bind_cols(d_voters, as.data.frame(P_mat), as.data.frame(P_s_mat))
    P_vec <- colSums(P_mat)
    P_s_vec <- colSums(P_s_mat)
    P_all <- rbind(first = P_vec, second = P_s_vec)
    list(d_voters = d_voters, P_vec = P_all)
  }


update_max2 <- function(chain, d_cands, theta, unit = 0.05, p_before, p_now) {
  stopifnot(nrow(p_now) == 2)
  J <- ncol(p_now)

  P_before_first <- p_before[1, ]
  P_before_second <- p_before[2, ]
  P_now_first <- p_now[1, ]
  P_now_second <- p_now[2, ]

  new_theta <- numeric(J)
  improve <- (P_now_first >= P_before_first) & (P_now_second >= P_before_second)
  new_theta[improve] <- theta[improve]
  new_theta[!improve] <- runif(sum(!improve), min = theta[!improve] + 90, max = theta[!improve] + 270)

  rad_theta <- new_theta * pi / 180
  dx <- cos(rad_theta) * unit
  dy <- sin(rad_theta) * unit

  d_cands$x <- d_cands$x + dx
  d_cands$y <- d_cands$y + dy

  list(d_cands = d_cands, theta = new_theta)
}

  get_util_max3 <- function(d_voters, d_cands, N_voters, eps_sd, m_vec) {
    J <- nrow(d_cands)
    V_mat <- matrix(NA_real_, nrow = N_voters, ncol = J)

    for (j in seq_len(J)) {
      D_j <- sqrt((d_voters$x - d_cands$x[j])^2 + (d_voters$y - d_cands$y[j])^2)
      m_j <- as.integer(d_voters$ethnic_group != paste0("Group ", m_vec[j]))
      eps_j <- rnorm(N_voters, sd = eps_sd)
      V_mat[, j] <- -d_voters$c * D_j - d_voters$b * m_j + eps_j
    }

    expV <- exp(V_mat)
    P_mat <- expV / rowSums(expV)
    colnames(P_mat) <- paste0("P", seq_len(J))

    # Second-choice
    P_s_mat <- matrix(0, nrow = N_voters, ncol = J)
    for (j in seq_len(J)) {
      for (k in setdiff(seq_len(J), j)) {
        V_den_ex_k <- rowSums(expV[, setdiff(seq_len(J), k), drop = FALSE])
        contrib <- (expV[, j] / V_den_ex_k) * P_mat[, k]
        P_s_mat[, j] <- P_s_mat[, j] + contrib
      }
    }
    colnames(P_s_mat) <- paste0("P", seq_len(J), "_s")

    # Third-choice
    P_t_mat <- matrix(0, nrow = N_voters, ncol = J)
    for (j in seq_len(J)) {
      for (k1 in setdiff(seq_len(J), j)) {
        for (k2 in setdiff(setdiff(seq_len(J), j), k1)) {
          V_den_ex_k1 <- rowSums(expV[, setdiff(seq_len(J), k1), drop = FALSE])
          V_den_ex_k2 <- rowSums(expV[, setdiff(seq_len(J), k2), drop = FALSE])
          prob_k1k2 <- (expV[, k2] / V_den_ex_k1) * P_mat[, k1] +
            (expV[, k1] / V_den_ex_k2) * P_mat[, k2]
          V_den_ex_pair <- rowSums(expV[, setdiff(seq_len(J), c(k1, k2)), drop = FALSE])
          P_t_mat[, j] <- P_t_mat[, j] + (expV[, j] / V_den_ex_pair) * prob_k1k2
        }
      }
    }
    colnames(P_t_mat) <- paste0("P", seq_len(J), "_t")

    d_voters <- dplyr::bind_cols(d_voters, as.data.frame(P_mat), as.data.frame(P_s_mat), as.data.frame(P_t_mat))
    P_vec <- rbind(colSums(P_mat), colSums(P_s_mat), colSums(P_t_mat))
    rownames(P_vec) <- c("first", "second", "third")

    list(d_voters = d_voters, P_vec = P_vec)
  }

update_max3 <- function(chain, d_cands, theta, unit = 0.05, p_before, p_now) {
  stopifnot(nrow(p_now) == 3)
  J <- ncol(p_now)

  P_now_first <- p_now[1, ]
  P_now_second <- p_now[2, ]
  P_now_third <- p_now[3, ]
  P_before_first <- p_before[1, ]
  P_before_second <- p_before[2, ]
  P_before_third <- p_before[3, ]

  improve <- (P_now_first >= P_before_first) &
    (P_now_second >= P_before_second) &
    (P_now_third >= P_before_third)

  new_theta <- numeric(J)
  new_theta[improve] <- theta[improve]
  new_theta[!improve] <- runif(sum(!improve), min = theta[!improve] + 90, max = theta[!improve] + 270)

  rad_theta <- new_theta * pi / 180
  dx <- cos(rad_theta) * unit
  dy <- sin(rad_theta) * unit

  d_cands$x <- d_cands$x + dx
  d_cands$y <- d_cands$y + dy

  list(d_cands = d_cands, theta = new_theta)
}


# Implement the simulation -----------------------------------------------------

pb <- progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = n_iter,
      complete = "=",   # completion bar character
      incomplete = "-", # incomplete bar character
      current = ">",    # current bar character
      clear = FALSE,    # if TRUE, clears the bar when finish
      width = 100       # width of the progress bar
    )

  # # # # # # Only for debuggin
  #   n_iter = 100
  #   mu_c = 1 # average spatial voting
  #   mu_b = 2 # average co-ethnic voting
  #   sigma_c = 0.5
  #   sigma_b = 0.5
  #   eps_sd = 0.5
  #   unit = 0.05
  #   seed = 14231
  #
  # data <- sim_data()

  # Set up -----------------------------------------------------------------------
  set.seed(seed)

  d_voters <- data$gen_voters
  d_cands <- data$gen_cands

  n_iter <- n_iter
  J <- dim(d_cands)[1] # number of candidates

  # Fixed parameters, created outside chains
  N_voters <- dim(d_voters)[1]
  m_vec <- rep(1:(J/2), 2) # used in J-loop

  d_voters$c <- rnorm(N_voters, mu_c, sigma_c)
  d_voters$b <- rnorm(N_voters, mu_b, sigma_b)

  get_util <- switch(strategy,
                     max1 = get_util_max1,
                     max2 = get_util_max2,
                     max3 = get_util_max3,
                     stop("Invalid strategy. Use 'max1', 'max2', or 'max3'."))

  update_positions <- switch(strategy,
                             max1 = update_max1,
                             max2 = update_max2,
                             max3 = update_max3,
                             stop("Invalid strategy. Use 'max1', 'max2', or 'max3'."))


# Initialize before the loop
  p_before <- matrix(NA_real_, nrow = 3, ncol = J)
  rownames(p_before) <- c("first", "second", "third")
  chain_voters <- list()
  chain_cands <- list()

# Loop over iterations
  for (t in seq_len(n_iter)) {

    pb$tick()

    if (t == 1) {
      theta <- runif(J, 0, 360)
    } else {
      theta <- chain_cands[[t - 1]]$theta
      d_cands <- chain_cands[[t - 1]]$d_cands
    }

    chain_voters[[t]] <- get_util(d_voters, d_cands, N_voters, eps_sd, m_vec)
    p_now <- chain_voters[[t]]$P_vec

    chain_cands[[t]] <- update_positions(
      chain_voters[[t]]$d_voters,
      d_cands,
      theta,
      unit,
      p_before = if (t == 1) p_now else p_before,
      p_now = p_now
    )

    chain_voters[[t]]$d_voters$iter <- t
    chain_cands[[t]]$d_cands$iter <- t
    p_before <- p_now
  }



# combine all results
out_voters <- bind_rows(purrr::map(chain_voters, "d_voters")) %>% tibble()
out_cands <- bind_rows(purrr::map(chain_cands, "d_cands")) %>% tibble()


out <- list(
  voters = out_voters,
  cands = out_cands
)

return(out)
}


# # # # Check via visualization
# out1 <- neodowns(data, strategy = "max1", n_iter = 2000)
# out2 <- neodowns(data, strategy = "max2", n_iter = 2000)
# out3 <- neodowns(data, strategy = "max3", n_iter = 2000)
#
# cands1 <- out1$cands %>% mutate(strategy = "max1")
# cands2 <- out2$cands %>% mutate(strategy = "max2")
# cands3 <- out3$cands %>% mutate(strategy = "max3")
#
# cands_all <- bind_rows(cands1, cands2, cands3)
#
#
# cands_all %>%
# #  filter(iter > 1000) %>%
# ggplot(aes(x = x, y = y, color = party, group = party)) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
#   geom_path(linewidth = 0.01, arrow = arrow(type = "open", length = unit(0.15, "cm"))) +
#   geom_point(alpha = 0.7, size = 0.05) +
#   geom_point(data = subset(out1$cands, iter == 1),
#              shape = 21, fill = "black", size = 3, stroke = 1.2) +
#   geom_text(data = subset(out1$cands, iter == 1),
#             aes(label = party), vjust = -1, size = 3, color = "black") +
#   labs(
#     title = "Trajectory of Candidate Positions Over Iterations",
#     x = "", y = "", color = "Candidate ID"
#   ) +
#   coord_equal() +
#   theme_minimal() +
#   facet_wrap(~ strategy) +
#   theme(legend.position = "none")




# Burn-in
# Random Sampling from the Chains (to address autocorrelation)
# Take every 25th chain
burn_samp <- function(chain, burnin = 1000, int = 1) {
  chain <- pull(chain) # Vectorize the input chain
  burned <- chain[burnin:length(chain)] # Drop burn-in periods ("burnin")
  sampled <- burned[seq(1, length(burned), by = int)] # Sample each "int" statistic
  return(sampled)
}
