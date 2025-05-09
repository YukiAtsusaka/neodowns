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
                     # force = TRUE,
                     # while_max = 200,
                     # boost = 6,
                     seed = 14231) {

# # # # # # Only for debuggin
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
  d_cands <- data$gen_cands

  n_iter <- n_iter
  J <- dim(d_cands)[1] # number of candidates

  # Fixed parameters, created outside chains
  N_voters <- dim(d_voters)[1]
  m_vec <- rep(1:(J/2), 2) # used in J-loop

  d_voters$c <- rnorm(N_voters, mu_c, sigma_c)
  d_voters$b <- rnorm(N_voters, mu_b, sigma_b)


# Compute voter utility and candidate expected vote share ----------------------


get_util_vs <- function(d_voters, d_cands,
                        N_voters, esp_sd){


  V_mat <- matrix(NA_real_, nrow = N_voters, ncol = J)
  colnames(V_mat) <- paste0("V", seq_len(J))

  # Voter utility
  for (j in seq_len(J)) {
    D_j <- sqrt((d_voters$x - d_cands$x[j])^2 + (d_voters$y - d_cands$y[j])^2)
    m_j <- as.integer(d_voters$ethnic_group != paste0("Group ", m_vec[j]))
    eps_j <- rnorm(N_voters, sd = eps_sd)
    V_mat[, j] <- -d_voters$c * D_j - d_voters$b * m_j + eps_j
  }

  expV <- exp(V_mat) # exponentiated utilities
  P_mat <- expV / rowSums(expV) # random-utility model exp(V)/sum(exp(V))

  d_voters <- dplyr::bind_cols(d_voters, as.data.frame(V_mat))
  colnames(P_mat) <- paste0("P", seq_len(J))
  d_voters <- dplyr::bind_cols(d_voters, as.data.frame(P_mat))

  P_vec <- colSums(P_mat)

# Return
list(
  d_voters = d_voters,
  P_vec = P_vec
)

}


# Update positions -------------------------------------------------------------
update <- function(chain, d_cands,
                   theta, unit = 0.05,
                   p_before, p_now) {

  J <- length(theta)
  new_theta <- numeric(J)

  # Compute updated directions all at once
  improve <- p_now >= p_before
  new_theta[improve] <- theta[improve]
  new_theta[!improve] <- runif(sum(!improve),
                               min = theta[!improve] + 90,
                               max = theta[!improve] + 270)

  # Precompute cos/sin once
  rad_theta <- new_theta * pi / 180
  dx <- cos(rad_theta) * unit
  dy <- sin(rad_theta) * unit

  # Update positions
  d_cands$x <- d_cands$x + dx
  d_cands$y <- d_cands$y + dy

  # Return
  list(
    d_cands = d_cands,
    theta = new_theta
  )
}



# Implement the simulation -----------------------------------------------------

# pb <- progress_bar$new(
#       format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
#       total = n_iter,
#       complete = "=",   # completion bar character
#       incomplete = "-", # incomplete bar character
#       current = ">",    # current bar character
#       clear = FALSE,    # if TRUE, clears the bar when finish
#       width = 100       # width of the progress bar
#     )
#
#   tick_every <- max(1, floor(n_iter / 50))  # aim for ~50 updates max

# Initialize before the loop
p_before <- rep(NA_real_, J)
chain_voters <- list()
chain_cands <- list()

for (t in seq_len(n_iter)) {

  # if (t %% tick_every == 0 || t == n_iter) pb$tick()

  # Use previous candidate positions and theta
  if (t == 1) {
    theta <- runif(J, 0, 360)
    d_cands <- d_cands  # unchanged
  } else {
    theta <- chain_cands[[t - 1]]$theta
    d_cands <- chain_cands[[t - 1]]$d_cands
  }

  # Compute current voter utility and vote shares
  chain_voters[[t]] <- get_util_vs(d_voters, d_cands, N_voters, esp_sd)
  p_now <- chain_voters[[t]]$P_vec

  # Update positions based on previous and current vote shares
  chain_cands[[t]] <- update(chain_voters[[t]]$d_voters,
                             d_cands,
                             theta,
                             unit = 0.05,
                             p_before = if (t == 1) p_now else p_before,
                             p_now = p_now)

  # Store iteration info
  chain_voters[[t]]$d_voters$iter <- t
  chain_cands[[t]]$d_cands$iter <- t

  # Update p_before for next iteration
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


# # Check via visualization
# out$cands %>%
# #  filter(iter > 1000) %>%
# ggplot(aes(x = x, y = y, color = party, group = party)) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
#   geom_path(linewidth = 0.01, arrow = arrow(type = "open", length = unit(0.15, "cm"))) +
#   geom_point(alpha = 0.7, size = 0.05) +
#   geom_point(data = subset(out$cands, iter == 1),
#              shape = 21, fill = "black", size = 3, stroke = 1.2) +
#   geom_text(data = subset(out$cands, iter == 1),
#             aes(label = party), vjust = -1, size = 3, color = "black") +
#   labs(
#     title = "Trajectory of Candidate Positions Over Iterations",
#     x = "", y = "", color = "Candidate ID"
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
