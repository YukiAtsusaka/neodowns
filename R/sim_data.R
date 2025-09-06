#' Simulate voter and candidate data
#'
#' @description `sim_data()` simulates ideological positions of voters and candidates.
#'
#' @param n_voter Number of voters in the simulation.
#' @param n_group Number of groups in the simulation.
#' @param skew Skew---s parameter---of the zipf distribution for generating relative balance of ethnic groups,
#' a high skew implies the biggest group is much bigger then other groups while
#' a low skew implies more balance, skew close to zero such as .1 will result in equal groups
#' .6 ensures no group has an outright majority. 1 results in a majority group twice as large as other groups. 
#' @param dev Parameter that tweaks the variance of the various bits of code that generate voter ideological positions,
#' a high dev (around 1) implies more variance within groups, a low dev implies less variance within groups.
#' @param dist Voter distribution, supporting 'Normal,' 'Flat,' 'Polarized,' or 'Clustered'
#' @param init 'initial distance' between candidates/groups.
#' a high init parameter indicates groups in distributions like the polarized and cluster should form far away from the origin wih no overlap, whereas
#' a low init parameter implies groups clustered around the origin with significant overlap.
#' @param eth Parameter that controls for the separation of ideological space by ethnicity.
#' eth = .5 is the maximum and eth = 1 is the minimum separation.
#' @param seed Random seed for reproducibility.
#' @importFrom dplyr case_when `%>%` mutate select left_join sample_n
#' @importFrom tidyr pivot_longer
#' @importFrom VGAM rzipf rrayleigh
#' @export

sim_data <- function(n_voter = 1000,
                     n_cand = 3,
                     n_group = 3,
                     skew = 2,
                     dev = 0.75,
                     dist = "Normal",
                     init = 1,
                     eth = 1,
                     seed = NULL
                     ) {

  set.seed(seed)

  # Step 1: Generate voters' ethnicity
  ethnic_group <- rzipf(n_voter, N = n_group, s = skew)
  ethnic <- tabulate(ethnic_group)
  freq <- ethnic / sum(ethnic)

  voters <- data.frame(ethnic_group)
  voters$cx <- rnorm(n_voter, 0, dev) / 4
  voters$cy <- rnorm(n_voter, 0, dev) / 4

  shares <- 360 * freq
  sharediv <- shares / 2
#  angles <- c(0, sharediv[1:n_group - 1] + sharediv[2:n_group]) this causes an errro when n.group = 1
  angles <- c(0, head(sharediv, -1) + tail(sharediv, -1))
  starting_angles <- cumsum(angles) + 15

  u <- runif(n = n_voter, min = -1, max = 1)
  rand <- u * n_group * eth * shares[voters$ethnic_group] / 2
  voters$eth <- ifelse(runif(n_voter, 0, 1) < eth / 2, sample(1:3, n_voter, replace = TRUE), voters$ethnic_group)
  voterAngle <- starting_angles[voters$ethnic_group] + rand
  voterAngle_f <- starting_angles[voters$eth]

  myR <- rrayleigh(n_voter, init)

  if (dist != "Fan") {
    v1 <- cos(voterAngle * pi / 180)
    v2 <- sin(voterAngle * pi / 180)
    voters$x <- v1 * myR
    voters$y <- v2 * myR
     if (dist == "Polarized") {
      voters <- voters[order(voters$ethnic_group), ]
      cutoff <- sum(ethnic[1:floor(length(ethnic) / 2)])
      voters$x[1:cutoff] <- voters$x[1:cutoff] / 10 + rnorm(cutoff, init, dev / 2)
      voters$x[(cutoff + 1):n_voter] <- voters$x[(cutoff + 1):n_voter] / 10 + rnorm(n_voter - cutoff, -init, dev / 2)
    }

  } else {
    v1 <- cos(voterAngle_f * pi / 180)
    v2 <- sin(voterAngle_f * pi / 180)
    voters$x <- voters$cx + (v1 * myR)
    voters$y <- voters$cy + (v2 * myR)
  }

  if (dist == "Clustered") {
    clusters <- data.frame(group = c(1:n_group))
    v1 <- cos(starting_angles * pi / 180)
    v2 <- sin(starting_angles * pi / 180)
    clusters$x <- v1 * 1.5 * init
    clusters$y <- v2 * 1.5 * init
    clustx <- rnorm(n_voter, 0, dev / n_group)
    clusty <- rnorm(n_voter, 0, dev / n_group)
    voters$x <- voters$cx + clusters$x[voters$eth] * init * 1.2 + clustx
    voters$y <- voters$cy + clusters$y[voters$eth] * init * 1.2 + clusty
  }


  # Finalize voters
  voters <- voters %>%
    mutate(
      x = 8 * (x - min(x)) / (max(x) - min(x)) - 4,
      y = case_when(
        dist == "Flat" | dist =="Polarized" ~ 2 * (y - min(y)) / (max(y) - min(y)) - 1,
        TRUE ~ 8 * (y - min(y)) / (max(y) - min(y)) - 4),
      ethnic_group = paste0("Group ", ethnic_group),
      group = eth
    ) %>%
    dplyr::select(group, ethnic_group, x, y) %>%
    tibble()

  # Generate candidates
  cands <- voters%>%
    dplyr::sample_n(n_cand) %>%
    dplyr::mutate(candidate = 1:n_cand)


  # Output
  out <- list(
    gen_voters = voters,
    gen_cands = cands
  )

  return(out)
}



