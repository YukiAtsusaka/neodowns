#' Plot candidate behavior inside and outside the burn-in period
#'
#' @description `plot_trace()` makes the trace plot for one or more chains.
#'
#' @param out Output from `neodowns`
#' @return ggplot output
#' @importFrom dplyr `%>%` bind_rows tibble
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap theme_classic xlab ylab theme
#' @export

plot_trace <- function(out){

out_list <- out %>%
  dplyr::bind_rows(.id = 'chain') %>%
  tibble()


ggplot(out_list, aes(x = iter, y = estimate)) +
   geom_line(aes(color = chain)) +
   facet_wrap(~ system) +
   theme_classic() +
   xlab('iteration') +
   ylab('') +
   theme(legend.position="none") -> p

  return(p)

}
