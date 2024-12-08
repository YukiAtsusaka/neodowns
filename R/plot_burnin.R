#' Plot candidate behavior inside and outside the burn-in period
#'
#' @description `plot_burnin()` visualizes the sequence of candidate positions.
#'
#' @param out Output from `neodowns`
#' @param n_iter Number of iterations to plot
#' @return ggplot output
#' @importFrom dplyr `%>%` filter
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline scale_size scale_color_manual scale_shape_manual facet_wrap ylim xlim theme_minimal
#' @export

# library(gganimate)
# library(gapminder)
# library(ggplot2)

plot_burnin <- function(out,
                        n_iter){


out_cands <- out$cands %>%
  filter(iter <= n_iter)

ggplot(out_cands, aes(x = x, y = y,
                      colour = party,
                      shape = factor(party))) +
  geom_point(aes(alpha = iter/n_iter), show.legend = F) +
  geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed") +
  scale_size(range = c(2, 10)) +
  scale_color_manual(
    values = c(
      "deepskyblue4",
      "deeppink4",
      "slategray",
      "deepskyblue4",
      "deeppink4",
      "slategray")
    ) +
  scale_shape_manual(values = c(1, 2, 0, 16, 17, 15)) +
  facet_wrap(~ system, ncol = 2) +
  ylim(-4, 4) +
  xlim(-4, 4) +
  theme_minimal() -> p

return(p)

}
