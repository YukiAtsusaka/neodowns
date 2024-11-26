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
                     skew = 0.6,
                     dev = 0.75,
                     dist = "Normal",
                     init = 1,
                     eth = 1/3,
                     seed = 1524) {

  set.seed(seed)

# Step 1: Generate voters ethnicity (group affiliation) ========================
  ## based on zipfs law
  ## "s" determines the distribution of groups
  ## Allocate N_voters into N_groups

  ethnic_group <- rzipf(N_voters,
    N = N_groups,
    s = skew
  )

  ## local objects
  ethnic <- tabulate(ethnic_group) # this will used at the end
  freq <- ethnic / sum(ethnic)

# Step 2: Generate voters' ideological positions ==============================
  ## First, generate small, random numbers by drawing from a uniform distribution
  ## that is symmetrically bounded by "dev"
  ## the higher "dev" becomes, the more initial variation

  voters <- data.frame(ethnic_group) # voters in rows

  # Generate cx: noise (random values) in x-axis
  # Generate cy: noise (random values) in y-axis
  # cx and cy are tiny values, which guarantees that each voter has a distinct
  # position with no overlap (idiosyncratic factors)

  voters$cx <- rnorm(N_voters, 0, dev) / 4
  voters$cy <- rnorm(N_voters, 0, dev) / 4


  # create a frame for plotting the midlines of each ethhnic distribution
  cands <- data.frame(
    group = c(1:N_groups)
  )

  # Generate voters' positions by distribution type

  # Create a few arguments
  # For all distributions, except for clustered

  # This determines the starting angle for each group
  # for each ethnic group
  # angles: a vector of length N_group

  # angle from which each group starts
  # First divide the space into shares based on the frequency of groups from the zipf distribution

  shares <- 360 * freq
  sharediv <- shares / 2

  # Find the midpoint of each share and adjust the distribution by 15 degrees

  angles <- c(0, sharediv[1:N_groups - 1] + sharediv[2:N_groups])
  starting_angles <- cumsum(angles) + 15

  # With the starting angle, we do the following

  # Random adjustment to each voter's angle (for technical reason)
  u <- runif(n = N_voters, min = -1, max = 1)

  # group can vary by 60 degree (for N_groups == 3)
  rand <- u * N_groups * eth * shares[voters$ethnic_group] / 2

  voters$eth <- ifelse(runif(N_voters, 0, 1) < eth / 2, sample(1:3, N_voters, replace = TRUE), voters$ethnic_group)

  # This is where we can potentially include ethnic_ideology
  voterAngle <- starting_angles[voters$ethnic_group] + rand
  voterAngle_f <- starting_angles[voters$eth]

  # myR (radius) is the distance from the origin in a polar coordinate system,
  # if we generate it for a given point with a normal distribution we are generating a bivariate distribution that has a 'normal' radius

  myR <-rrayleigh(N_voters,init)
    #abs(rnorm(N_voters,0, .5 * init))

  # Family I (Normal)
  ## 1. Basic Normal distribution
  if (dist != "Fan") {

    # Generate Cartesian position for each voter from circular coordinates
    # v1 and v2 generate x and y adjustments based on angles for polar coordinate to cartesian coordinate adjustment

    v1 <- cos(voterAngle * pi / 180)
    v2 <- sin(voterAngle * pi / 180)

    # Key part --- Determine each voter's position
    ## voters$x: ideology in dimension x
    ## voters$y: ideology in dimension y

    voters$x <- (v1 * myR) # Finish polar to x cartesian coordinate converion with normal radius + angle based conversion factor
    voters$y <- (v2 * myR) # Finish polar to y cartesian coordinate converion with half-normal radius + angle based conversion factor

    ## 2. Flat distribution
    if (dist == "Flat") {
      voters$x <- voters$x * 1.25
      voters$y <- voters$y / 5
      ## 3. Polarized distribution
    } else if (dist == "Polarized") {
      voters <- voters[order(voters$ethnic_group), ]
      cutoff <- sum(ethnic[1:floor(length(ethnic) / 2)])
      voters$x[1:cutoff] <- voters$x[1:cutoff] / 10 + rnorm(cutoff, init, dev / 2)
      voters$x[(cutoff + 1):N_voters] <- voters$x[(cutoff + 1):N_voters] / 10 + rnorm(N_voters - cutoff, -init, dev / 2)
    }

    # Family II (Fan distribution)
  } else {
    # Generate Cartesian position for each voter
    # Generate Cartesian position for each voter from circular coordinates

    # v1 and v2 generate x and y adjustments based on angles for circular coordinate to cartesian coordinate adjustment

    v1 <- cos(voterAngle_f * pi / 180)
    v2 <- sin(voterAngle_f * pi / 180)

    # Key part --- Determine each voter's position
    ## voters$x: ideology in dimension x
    ## voters$y: ideology in dimension

    voters$x <- voters$cx + (v1 * myR) # Finish circular to x cartesian coordinate converion with half-normal radius + angle based conversion factor
    # Because fan distribution has no angle variation, cartesian variation is introduced via cx and cy
    voters$y <- voters$cy + (v2 * myR) # Finish circular to y cartesian coordinate converion with half-normal radius + angle based conversion factor
  }

  # Family III (Clustered distribution)

  if (dist == "Clustered") {
    clusters <- data.frame(
      group = c(1:N_groups)
    )

    v1 <- cos(starting_angles * pi / 180)
    v2 <- sin(starting_angles * pi / 180)

    clusters$x <- v1 * 1.5 * init # circular to cartesian conversion x for N_groups single points to determine the centers of each cluster.
    clusters$y <- v2 * 1.5 * init # circular to cartesian conversion y for N_groups single points to determine the centers of each cluster.
    clustx <- rnorm(N_voters, 0, dev / N_groups)
    clusty <- rnorm(N_voters, 0, dev / N_groups)
    voters$x <- voters$cx + clusters$x[voters$eth] * init * 1.2 + clustx
    voters$y <- voters$cy + clusters$y[voters$eth] * init * 1.2 + clusty
  }
  cands$group <- c(1:N_groups)

# Step 3: Generate candidates' initial positions ===============================
  ## by creating a midline for each group
  ## populate extreme and moderate parties at specific intervals
  ## along this midline and voters symmetrically around it

  cl <- rainbow(N_groups)

  ## Divide the entire space into i (# groups)
  if (dist != "Clustered" & dist != "Polarized") {
    x_vectors <- cos((starting_angles * pi) / 180) # determines how parties are distributed on the angles in X
    y_vectors <- tan((starting_angles * pi) / 180) * ifelse(dist == "Flat", init * .1, 1) # determines how parties are distributed on the angles in Y
    cands$x_mod <- init * .75 * x_vectors # X value for moderate candidates
    cands$x_ext <- init * 1.5 * x_vectors # X value for extreme candidates
    cands$y_mod <- cands$x_mod * y_vectors # Y value for moderate candidate
    cands$y_ext <- cands$x_ext * y_vectors # Y value for extreme candidate
  } else {

    x <- voters %>%
      group_by(ethnic_group) %>%
      summarize(x = mean(x))
    y <- voters %>%
      group_by(ethnic_group) %>%
      summarize(y = mean(y))

    cands$x_mod <- x$x * init * .85
    cands$y_mod <- y$y * init * .85
    cands$x_ext <- x$x * init * 1.3 # X value for extreme candidates
    cands$y_ext <- y$y * init * 1.3 # X value for extreme candidates
  }

  # Limit ideologies between -4 and 4
  voters <- voters %>%
    mutate(
      x = case_when(
        x > 4 ~ 4,
        x < -4 ~ -4,
        TRUE ~ x
      ),
      y = case_when(
        y > 4 ~ 4,
        y < -4 ~ -4,
        TRUE ~ y
      ),
      ethnic_group = paste0("Group ", ethnic_group),
      group = eth
    ) %>%
    dplyr::select(group, ethnic_group, x, y) %>%
    tibble()


  cands <- cands %>%
    pivot_longer(
      x_mod:y_ext,
      cols_vary = "slowest",
      names_to = c(".value", "type"),
      names_sep = "_"
    ) %>%
  mutate(party = paste0("Group_", group, "_", type),
         ethnic_group = paste0("Group ", group)) %>%
    dplyr::select(group, ethnic_group, x, y, type, party)


# Output
  out <- list(
    gen_voters = voters, # voter data
    gen_cands = cands   # candidate data
  )

  return(out)
}
