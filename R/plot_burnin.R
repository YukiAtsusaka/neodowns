#' Plot candiate behavior inside and outside the burnin period
#'
#' @description `plot_burnin()` visualizes the sequence of candiate positions.
#'
#' @param out Output from `neodowns`
#'  @return
#' @importFrom dplyr case_when `%>%` mutate select left_join rename tibble filter bind_rows
#' @importFrom ggplot2 ggplot
#' @export

# library(gganimate)
# library(gapminder)
# library(ggplot2)

plot_burnin <- function(out){

n_iter <- 100
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
      alpha("deepskyblue4", 0.8),
      alpha("deeppink4", 0.8),
      alpha("slategray", 0.8),
      alpha("deepskyblue4", 0.8),
      alpha("deeppink4", 0.8),
      alpha("slategray", 0.8)
    )) +
  scale_shape_manual(values = c(1, 2, 0, 16, 17, 15)) +
  facet_wrap(~ system, ncol = 2) +
  ylim(-4, 4) +
  xlim(-4, 4) +
  theme_minimal()


q <- extract(out, quantity = "ideology")

ggplot(q, aes(y = `mean(moderation)`, x = iter)) +
  geom_point() +
  facet_wrap(~ system, ncol = 2)

#
#
#   ggplot(out_cands, aes(x = x, y = y,
#                   size = count,
#                   colour = party, shape = factor(party))) +
#     geom_point(alpha = 0.7, show.legend = F) +
#     geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") +
#     geom_vline(aes(xintercept = 0), color = "black", linetype = "dashed") +
#     scale_size(range = c(2, 10)) +
#     scale_color_manual(
#       values = c(
#         alpha("deepskyblue4", 0.8),
#         alpha("deeppink4", 0.8),
#         alpha("slategray", 0.8),
#         alpha("deepskyblue4", 0.8),
#         alpha("deeppink4", 0.8),
#         alpha("slategray", 0.8)
#       )) +
#     scale_shape_manual(values = c(1, 2, 0, 16, 17, 15)) +
#     facet_wrap(~ type, ncol = 2) +
#     ylim(-3, 3) +
#     xlim(-3, 3) +
#     theme_minimal() +
#     # This is where gganimation begins
#     labs(title = 'Iteration: {frame_time}',
#          x = 'First dimension', y = 'Second dimension') +
#     transition_time(as.integer(iter)) +
#     ease_aes('quadratic-in') -> gif
#
#
#   gif
#
#   anim_save("animation.gif", gif)
#


}
