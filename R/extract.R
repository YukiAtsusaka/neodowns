#' Extracting quantities of interest from computer simulations
#'
#' @description `extract()` extract various quantities of interest from the output of `neodowns`.
#'
#' @param out output from `neodowns`
#' @return
#' A list, which has the following components:
#'  \item{d_voters}{a matrix containing the simulated ideological positions of voters.}
#'  \item{d_cands}{a matrix containing the simulated ideological positions of candidates.}
#' @importFrom dplyr mutate case_when `%>%` group_by summarise select
#' @export

extract <- function(out,
                    quantity,
                    n_burn = 1000,
                    int = 10
                    ) {

if(quantity %in% c("ideology", "ethnic")){}
  else{print("Error: quantitiy must be `ideology` or `ethnic` ")}


n_iter <- max(out$cands$iter)

out_cands <- out$cands %>%
    filter(iter > n_burn,
           iter %in% seq(n_burn, n_iter, by = int))


if(quantity == "ideology"){
return(
  out_cands %>%
    dplyr::group_by(system, iter) %>%
    dplyr::summarise(estimate = mean(moderation))
)
}

if(quantity == "ethnic"){
return(
    out$voters %>%
      mutate(
        ethnic_voting_first = case_when(
          ethnic_group == "Group 1" ~ P1 + P4,
          ethnic_group == "Group 2" ~ P2 + P5,
          ethnic_group == "Group 3" ~ P3 + P6,
          TRUE ~ NA),
        ethnic_voting_second = case_when(
          ethnic_group == "Group 1" ~ P1_s + P4_s,
          ethnic_group == "Group 2" ~ P2_s + P5_s,
          ethnic_group == "Group 3" ~ P3_s + P6_s,
          TRUE ~ NA)) %>%
    dplyr::select(iter, ethnic_group,
                  ethnic_voting_first,
                  ethnic_voting_second)
)
  }




}
