#' Extracting quantities of interest from computer simulations
#'
#' @description `extract()` extract various quantities of interest from the output of `neodowns`.
#'
#' @param out output from `neodowns`
#' @return
#' A list, which has the following components:
#'  \item{d_voters}{a matrix containing the simulated ideological positions of voters.}
#'  \item{d_cands}{a matrix containing the simulated ideological positions of candidates.}
#' @importFrom dplyr mutate case_when `%>%` group_by summarise select ungroup
#' @export


extract_ideology <- function(out,
                             n_burn = 1000,
                             int = 10,
                             subset = FALSE
                    ) {

n_iter <- max(out$cands$iter)

out_cands <- out$cands %>%
    filter(iter > n_burn,
           iter %in% seq(n_burn, n_iter, by = int))

if(subset == FALSE){

out <- out_cands %>%
  dplyr::group_by(system, iter) %>%
  dplyr::summarise(estimate = mean(moderation))

}else{

  out <- out_cands %>%
    dplyr::group_by(system, iter, type) %>%
    dplyr::summarise(estimate = mean(moderation))
}

return(out)
}


