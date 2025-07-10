#' Simulate voter and candidate data
#'
#' @description `sim_data()` simulates ideological positions of voters and candidates.
#'
#' @param data Simulated data.
#' @importFrom ggplot2 ggplot geom_point geom_hline geom_vline theme_minimal theme
#' @export
#'
#'
#'


see <- function(data){


# Extract voter positions
voters <- data$gen_voters

# Extract candidate positions
cands <- data$gen_cands

# Visualize voters and candidates

ggplot() +
  geom_point(data = voters,
             aes(x, y,
                 shape = ethnic_group,
                 color = ethnic_group),
             size = 1, alpha = 0.5) +
  geom_point(data = cands,
             aes(x, y,
                 shape = ethnic_group,
                 color = ethnic_group),
             size = 4) +
  geom_hline(aes(yintercept = 0), color = "black", alpha = 0.5, linetype = "dashed") +
  geom_vline(aes(xintercept = 0), color = "black", alpha = 0.5, linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "none") +
  ylab("") +
  xlab("")
}
