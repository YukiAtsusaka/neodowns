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
#'  @return
#' A list, which has the following components:
#'  \item{d_voters}{a matrix containing the simulated ideological positions of voters.}
#'  \item{d_cands}{a matrix containing the simulated ideological positions of candidates.}
#' @importFrom dplyr case_when `%>%` mutate select left_join rename tibble filter bind_rows pull
#' @importFrom progress progress_bar
#' @importFrom rlang := !!
#' @export


neodowns <- function(data,
                     strategy = NULL,
                     n_iter = 100,
                     mu_c = 1, # average spatial voting
                     mu_b = 2, # average co-ethnic voting
                     sigma_c = 0.5,
                     sigma_b = 0.5,
                     eps_sd = 0.5,
                     unit = 0.05,
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
    P_vec <- matrix(colSums(P_mat), nrow = 1)
    rownames(P_vec) <- "first"
    list(d_voters = d_voters, P_vec = P_vec)
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
    V_den_ex_list <- lapply(seq_len(J), function(k) rowSums(expV[, -k, drop = FALSE]))
    names(V_den_ex_list) <- seq_len(J)
    P_s_mat <- matrix(0, nrow = N_voters, ncol = J)
    for (j in seq_len(J)) {
      for (k in setdiff(seq_len(J), j)) {
        V_den_ex_k <- V_den_ex_list[[as.character(k)]]
        contrib <- (expV[, j] / V_den_ex_k) * P_mat[, k]
        P_s_mat[, j] <- P_s_mat[, j] + contrib
      }
    }
    colnames(P_s_mat) <- paste0("P", seq_len(J), "_s")
    d_voters <- dplyr::bind_cols(d_voters, as.data.frame(P_mat), as.data.frame(P_s_mat))
    P_vec <- rbind(first = colSums(P_mat), second = colSums(P_s_mat))
    list(d_voters = d_voters, P_vec = P_vec)
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

    V_den_ex_list <- lapply(seq_len(J), function(k) rowSums(expV[, -k, drop = FALSE]))
    names(V_den_ex_list) <- seq_len(J)

    P_s_mat <- matrix(0, nrow = N_voters, ncol = J)
    for (j in seq_len(J)) {
      for (k in setdiff(seq_len(J), j)) {
        V_den_ex_k <- V_den_ex_list[[as.character(k)]]
        contrib <- (expV[, j] / V_den_ex_k) * P_mat[, k]
        P_s_mat[, j] <- P_s_mat[, j] + contrib
      }
    }
    colnames(P_s_mat) <- paste0("P", seq_len(J), "_s")

    combs <- combn(seq_len(J), 2, simplify = FALSE)
    V_den_ex_pair_list <- list()
    for (pair in combs) {
      key <- paste(sort(pair), collapse = ",")
      V_den_ex_pair_list[[key]] <- rowSums(expV[, -pair, drop = FALSE])
    }

    P_t_mat <- matrix(0, nrow = N_voters, ncol = J)
    for (j in seq_len(J)) {
      others <- setdiff(seq_len(J), j)
      combs_j <- combn(others, 2, simplify = FALSE)
      for (pair in combs_j) {
        k1 <- pair[1]
        k2 <- pair[2]
        key <- paste(sort(pair), collapse = ",")
        V_den_ex_k1 <- V_den_ex_list[[as.character(k1)]]
        V_den_ex_k2 <- V_den_ex_list[[as.character(k2)]]
        V_den_ex_pair <- V_den_ex_pair_list[[key]]
        prob_k1k2 <- (expV[, k2] / V_den_ex_k1) * P_mat[, k1] +
          (expV[, k1] / V_den_ex_k2) * P_mat[, k2]
        P_t_mat[, j] <- P_t_mat[, j] + (expV[, j] / V_den_ex_pair) * prob_k1k2
      }
    }
    colnames(P_t_mat) <- paste0("P", seq_len(J), "_t")

    d_voters <- dplyr::bind_cols(d_voters, as.data.frame(P_mat), as.data.frame(P_s_mat), as.data.frame(P_t_mat))
    P_vec <- rbind(colSums(P_mat), colSums(P_s_mat), colSums(P_t_mat))
    rownames(P_vec) <- c("first", "second", "third")
    list(d_voters = d_voters, P_vec = P_vec)
  }

  # Update Functions --------------------------------------------------

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

# Candidate-specific strategy functions

strategy_map <- list(
  max1 = list(util = get_util_max1, update = update_max1),
  max2 = list(util = get_util_max2, update = update_max2),
  max3 = list(util = get_util_max3, update = update_max3)
)


pb <- progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = n_iter,
  complete = "=",   # completion bar character
  incomplete = "-", # incomplete bar character
  current = ">",    # current bar character
  clear = FALSE,    # if TRUE, clears the bar when finish
  width = 100       # width of the progress bar
)


# Implement the simulation -----------------------------------------------------

  # Only for debuggin
  #   strategy = c("max1", "max1", "max1",
  #                "max1", "max1", "max3")
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

  ##  Validate strategy length

  if (is.null(strategy)) strategy <- rep("max1", J)
  if (length(strategy) != J) stop("Length of strategy vector must match number of candidates.")
  if (J < 3 && any(strategy == "max3")) {
    warning("max3 requires at least 3 candidates. Downgrading all 'max3' strategies to 'max2'.")
    strategy[strategy == "max3"] <- "max2"
  } # To avoid errors for two-candidate races


  # Fixed parameters, created outside chains
  N_voters <- dim(d_voters)[1]
  m_vec <- rep(1:(J/2), 2) # used in J-loop

  d_voters$c <- rnorm(N_voters, mu_c, sigma_c)
  d_voters$b <- rnorm(N_voters, mu_b, sigma_b)

  # Initialize chains
  p_before <- matrix(NA_real_, nrow = 3, ncol = J)
  rownames(p_before) <- c("first", "second", "third")
  chain_voters <- vector("list", n_iter)
  chain_cands <- vector("list", n_iter)


  for (t in seq_len(n_iter)) {
    pb$tick()

    if (t == 1) {
      theta <- runif(J, 0, 360)
    } else {
      theta <- chain_cands[[t - 1]]$theta
      d_cands <- chain_cands[[t - 1]]$d_cands
    }


    # Determine which utility function to use
    if (J == 2 || all(strategy == "max1")) {
      full_util <- get_util_max1(d_voters, d_cands, N_voters, eps_sd, m_vec)
    } else if (all(strategy %in% c("max1", "max2"))) {
      full_util <- get_util_max2(d_voters, d_cands, N_voters, eps_sd, m_vec)
    } else {
      full_util <- get_util_max3(d_voters, d_cands, N_voters, eps_sd, m_vec)
    }

    chain_voters[[t]] <- full_util
    P_now_all <- full_util$P_vec


    # <-- MODIFIED: Loop through candidates individually to update positions
    for (j in seq_len(J)) {
      strat_j <- strategy[[j]]
      update_fun <- strategy_map[[strat_j]]$update

      p_now_j <- switch(strat_j,
                        max1 = matrix(P_now_all[1, j], nrow = 1),
                        max2 = matrix(P_now_all[1:2, j], nrow = 2),
                        max3 = matrix(P_now_all[1:3, j], nrow = 3))

      p_before_j <- switch(strat_j,
                           max1 = matrix(if (t == 1) P_now_all[1, j] else p_before[1, j], nrow = 1),
                           max2 = matrix(if (t == 1) P_now_all[1:2, j] else p_before[1:2, j], nrow = 2),
                           max3 = matrix(if (t == 1) P_now_all[1:3, j] else p_before[1:3, j], nrow = 3))

      update_res <- update_fun(
        full_util$d_voters,
        d_cands[j, , drop = FALSE],
        theta[j],
        unit,
        p_before = p_before_j,
        p_now = p_now_j
      )

      d_cands[j, c("x", "y")] <- update_res$d_cands[1, c("x", "y")]
      theta[j] <- update_res$theta[1]
    } # End of j loop


    chain_cands[[t]] <- list(d_cands = d_cands, theta = theta)
    chain_voters[[t]]$d_voters$iter <- t
    chain_cands[[t]]$d_cands$iter <- t
    chain_cands[[t]]$d_cands$strategy <- strategy

    p_before <- P_now_all
  } # End of t loop

  out_voters <- bind_rows(purrr::map(chain_voters, "d_voters")) %>% tibble()
  out_cands <- bind_rows(purrr::map(chain_cands, "d_cands")) %>% tibble()

  out <- list(
    voters = out_voters,
    cands = out_cands
  )
  return(out)
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
