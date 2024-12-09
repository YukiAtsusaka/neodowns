#' Plot candidate behavior inside and outside the burn-in period
#'
#' @description `plot_trace()` makes the trace plot for one or more chains.
#'
#' @param out Output from `neodowns`
#' @return ggplot output
#' @importFrom dplyr `%>%` bind_rows tibble
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap theme_classic xlab ylab theme geom_density
#' @export

plot_trace <- function(out){

out_list <- out %>%
  dplyr::bind_rows(.id = 'chain') %>%
  tibble()

mod <- out_list %>% filter(type == "mod")
ext <- out_list %>% filter(type == "ext")

  ggplot() +
    geom_density(data = ext, aes(x = estimate,
                          color = chain), alpha = 0.7) +
    geom_density(data = mod, aes(x = estimate,
                                 fill = chain), alpha = 0.7, color = F) +
    facet_wrap(~ system) +
    theme_classic() +
    xlab('estimate') +
    ylab('') +
    theme(legend.position="none") -> p

  return(p)

}
