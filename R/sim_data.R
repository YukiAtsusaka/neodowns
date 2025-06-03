#' Simulate voter and candidate data
#'
#' @description `sim_data()` simulates ideological positions of voters and candidates.
#'
#' @param N_voters Number of voters in the simulation.
#' @param N_groups Number of groups in the simulation.
#' @param skew Skew---s parameter---of the zipf distribution for generating relative balance of ethnic groups,
#' a high skew implies the biggest group is much bigger then other groups while
#' a low skew implies more balance,
#' .6 ensures no group has an outright majority.
#' @param dev Parameter that tweaks the variance of the various bits of code that generate voter ideological positions,
#' a high dev implies more variance within groups, a low dev implies less variance within groups.
#' @param dist Voter distribution, supporting 'Normal,' 'Flat,' 'Polarized,' or 'Clustered'
#' @param init 'initial distance' between candidates/groups.
#' a high init parameter indicates should form a circle around the origin with each group far away from it, whereas
#' a low init parameter implies groups clustered around the origin with significant overlap.
#' @param eth Parameter that controls for the separation of ideological space by ethnicity.
#' eth = 0 is the maximum and eth = 1 is the minimum separation.
#' @param seed Random seed for reproducibility.
#' @importFrom dplyr case_when `%>%` mutate select left_join
#' @importFrom tidyr pivot_longer
#' @importFrom VGAM rzipf rrayleigh
#' @export

sim_data <- function(N_voters = 1000,
                     N_groups = 3,
                     skew = 0.3,
                     dev = 0.75,
                     dist = "Normal",
                     init = 1,
                     eth = 1/3,
                     seed = 1524,
                     n_cands_by_group = rep(2, N_groups)  # vector of candidate counts by group
) {
  set.seed(seed)

  # Step 1: Generate voters' ethnicity
  ethnic_group <- rzipf(N_voters, N = N_groups, s = skew)
  ethnic <- tabulate(ethnic_group)
  freq <- ethnic / sum(ethnic)

  voters <- data.frame(ethnic_group)
  voters$cx <- rnorm(N_voters, 0, dev) / 4
  voters$cy <- rnorm(N_voters, 0, dev) / 4

  shares <- 360 * freq
  sharediv <- shares / 2
  angles <- c(0, sharediv[1:N_groups - 1] + sharediv[2:N_groups])
  starting_angles <- cumsum(angles) + 15

  u <- runif(n = N_voters, min = -1, max = 1)
  rand <- u * N_groups * eth * shares[voters$ethnic_group] / 2
  voters$eth <- ifelse(runif(N_voters, 0, 1) < eth / 2, sample(1:3, N_voters, replace = TRUE), voters$ethnic_group)
  voterAngle <- starting_angles[voters$ethnic_group] + rand
  voterAngle_f <- starting_angles[voters$eth]

  myR <- rrayleigh(N_voters, init)

  if (dist != "Fan") {
    v1 <- cos(voterAngle * pi / 180)
    v2 <- sin(voterAngle * pi / 180)
    voters$x <- v1 * myR
    voters$y <- v2 * myR

    if (dist == "Flat") {
      voters$x <- voters$x * 1.25
      voters$y <- voters$y / 5
    } else if (dist == "Polarized") {
      voters <- voters[order(voters$ethnic_group), ]
      cutoff <- sum(ethnic[1:floor(length(ethnic) / 2)])
      voters$x[1:cutoff] <- voters$x[1:cutoff] / 10 + rnorm(cutoff, init, dev / 2)
      voters$x[(cutoff + 1):N_voters] <- voters$x[(cutoff + 1):N_voters] / 10 + rnorm(N_voters - cutoff, -init, dev / 2)
    }

  } else {
    v1 <- cos(voterAngle_f * pi / 180)
    v2 <- sin(voterAngle_f * pi / 180)
    voters$x <- voters$cx + (v1 * myR)
    voters$y <- voters$cy + (v2 * myR)
  }

  if (dist == "Clustered") {
    clusters <- data.frame(group = c(1:N_groups))
    v1 <- cos(starting_angles * pi / 180)
    v2 <- sin(starting_angles * pi / 180)
    clusters$x <- v1 * 1.5 * init
    clusters$y <- v2 * 1.5 * init
    clustx <- rnorm(N_voters, 0, dev / N_groups)
    clusty <- rnorm(N_voters, 0, dev / N_groups)
    voters$x <- voters$cx + clusters$x[voters$eth] * init * 1.2 + clustx
    voters$y <- voters$cy + clusters$y[voters$eth] * init * 1.2 + clusty
  }

  # Validate input length
  if (length(n_cands_by_group) != N_groups) {
    stop("Length of n_cands_by_group must equal N_groups.")
  }

  # Generalized candidate generation (simplified)
  cand_list <- vector("list", N_groups)
  candidate_counter <- 1  # 🆕 Unique candidate ID

  for (g in 1:N_groups) {
    angle_rad <- starting_angles[g] * pi / 180
    direction <- c(cos(angle_rad), sin(angle_rad))
    center <- direction * init
    n_cand <- n_cands_by_group[g]
    pos_seq <- seq(-1, 1, length.out = n_cand)

    cand_df <- data.frame(
      group = g,
      ethnic_group = paste0("Group ", g),
      x = center[1] + pos_seq * init * 0.5,
      y = center[2] + pos_seq * init * 0.5,
      candidate = candidate_counter:(candidate_counter + n_cand - 1)  # 🆕 Unique ID
    )

    candidate_counter <- candidate_counter + n_cand  # 🆕 Increment
    cand_list[[g]] <- cand_df
  }

  cands <- bind_rows(cand_list)

  voters <- voters %>%
    mutate(
      x = case_when(x > 4 ~ 4, x < -4 ~ -4, TRUE ~ x),
      y = case_when(y > 4 ~ 4, y < -4 ~ -4, TRUE ~ y),
      ethnic_group = paste0("Group ", ethnic_group),
      group = eth
    ) %>%
    dplyr::select(group, ethnic_group, x, y) %>%
    tibble()

  cands <- cands %>%
    dplyr::select(group, ethnic_group, x, y, candidate)  # 🆕 Final structure

  out <- list(
    gen_voters = voters,
    gen_cands = cands
  )

  return(out)
}
