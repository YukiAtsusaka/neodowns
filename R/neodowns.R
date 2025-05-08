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
                     n_iter = 100,
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

# # # # # Only for debuggin
#   n_iter = 100
#   mu_c = 1 # average spatial voting
#   mu_b = 2 # average co-ethnic voting
#   sigma_c = 0.5
#   sigma_b = 0.5
#   eps_sd = 0.5
#   unit = 0.05
#   force = TRUE
#   while_max = 200
#   boost = 6
#   seed = 14231
#
# data <- sim_data()

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


  # Feeding the voter distribution to "initial chain" for both systems
  init <- d_voters

  init <- init %>%
    mutate(c = rnorm(n = N_voters, mean = mu_c, sd = sigma_c),
           b = rnorm(n = N_voters, mean = mu_b, sd = sigma_b))


# Compute voter utility and candidate expected vote share ---------------------


get_util_vs <- function(init, d_cands, N_voters, esp_sd){

  # Utility function
  for (j in seq_len(J)) {
    D_j <- sqrt((init$x - d_cands$x[j])^2 + (init$y - d_cands$y[j])^2)
    m_j <- ifelse(init$ethnic_group != paste0("Group ", m_vec[j]), 1, 0)
    eps_j <- rnorm(n = N_voters, sd = eps_sd)
    V_j <- -init$c * D_j - init$b * m_j + eps_j
    init[[paste0("V", j)]] <- V_j
  }

  # Choice probability
  V_mat <- as.matrix(init[, paste0("V", seq_len(J))])
  expV <- exp(V_mat)
  V_den <- rowSums(expV)

  for (j in seq_len(J)) {
    init[[paste0("P", j)]] <- expV[, j] / V_den
  }

  # Save the sum of probabilities
  P_cols <- paste0("P", seq_len(J))
  P_vec <- sapply(P_cols, function(col) sum(init[[col]]))
  names(P_vec) <- P_cols


# Return
list(
  init = init,
  P_vec = P_vec
)

}


update <- function(chain, d_cands, theta, unit = 0.05,
                   p_before, p_now) {

  J <- length(theta)
  new_theta <- numeric(J)


  for (j in seq_len(J)) {

    # Determine whether to continue or reverse direction
    improve <- p_now[j] >= p_before[j]
    new_theta[j] <- if (improve) {
      theta[j]
    } else {
      runif(1, min = theta[j] + 90, max = theta[j] + 270)
    }

    # Proposed new location
    propose_x <- d_cands$x[j] + cos(new_theta[j] / 180 * pi) * unit
    propose_y <- d_cands$y[j] + sin(new_theta[j] / 180 * pi) * unit

    # Update candidate position
    d_cands$x[j] <- propose_x
    d_cands$y[j] <- propose_y
  }

  # Return
  list(
    d_cands = d_cands,
    theta = new_theta
  )
}



# Repeat the two functions

pb <- progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = n_iter,
      complete = "=",   # Completion bar character
      incomplete = "-", # Incomplete bar character
      current = ">",    # Current bar character
      clear = FALSE,    # If TRUE, clears the bar when finish
      width = 100       # Width of the progress bar
    )


chain_voters <- list()
chain_cands <- list()

for (t in seq_len(n_iter)) {

pb$tick()

# initialization
if(t == 1){

  theta <- runif(n = J, min = 0, max = 360)
  d_cands <- d_cands

}else{
  theta <- chain_cands[[t-1]]$theta
  d_cands <- chain_cands[[t-1]]$d_cands

}



# Get voter utility and candidate vote share
chain_voters[[t]] <- get_util_vs(init, d_cands,
                                 N_voters, esp_sd)

p_now <- chain_voters[[t]]$P_vec
p_before <- if(t == 1){
  chain_voters[[t]]$P_vec}else{
    chain_voters[[t-1]]$P_vec
    }


# Update positions
chain_cands[[t]] <- update(chain_voters[[t]]$init,
                           d_cands,
                           theta,
                           unit = 0.05, p_before, p_now)

# Record iteration
chain_voters[[t]]$init$iter <- t
chain_cands[[t]]$d_cands$iter <- t

}



# combine all results
out_voters <- bind_rows(purrr::map(chain_voters, "init")) %>% tibble()
out_cands <- bind_rows(purrr::map(chain_cands, "d_cands")) %>% tibble()


out <- list(
  voters = out_voters,
  cands = out_cands
)

return(out)
}


# # # Check via visualization
# out$cands %>%
#   filter(iter > 1000) %>%
# ggplot(aes(x = x, y = y, color = party, group = party)) +
#   geom_path(linewidth = 0.01, arrow = arrow(type = "open", length = unit(0.15, "cm"))) +
#   geom_point(alpha = 0.7, size = 0.05) +
#   geom_point(data = subset(out$cands, iter == 1),
#              shape = 21, fill = "black", size = 3, stroke = 1.2) +
#   geom_text(data = subset(out$cands, iter == 1),
#             aes(label = party), vjust = -1, size = 3, color = "black") +
#   labs(
#     title = "Trajectory of Candidate Positions Over Iterations",
#     x = "X Position", y = "Y Position", color = "Candidate ID"
#   ) +
#   coord_equal() +
#   theme_minimal() +
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
