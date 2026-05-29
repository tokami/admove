#' Skipjack tuna in the Eastern Pacific Ocean
#'
#' A simulated data set of data-storage and mark-recapture tags for skipjack
#' tuna in the Eastern Pacific Ocean
#'
#' @format
#' A list with:
#' \describe{
#'   \item{sim}{The complete simulated data set of class `admove_sim`}
#'   \item{grid}{The simulation grid of class `admove_grid`}
#'   \item{cov}{The covariance fields of class `admove_cov`}
#'   \item{ctags}{Simulated mark-recapture tags}
#'   \item{dtags}{Simulated data-storage tags}
#'   \item{fit}{Fitted admove object of class `admove`}
#' }
#' @source simulated with admove
"skjepo"



#' Montagu's harrier
#'
#' A simulated data set of mark-resight tags for Montagu's harrier
#'
#' @format
#' The complete simulated data set of class `admove_sim` containing:
#' \describe{
#'   \item{grid}{The simulation grid of class `admove_grid`}
#'   \item{cov}{The covariance fields of class `admove_cov`}
#'   \item{par_sim}{The "true" parameters used for simulation}
#'   \item{tags}{Simulated tags (here: make-resight tags)}
#'   \item{dat}{The `admove` input data of class `admove_data`}
#'   \item{conf}{A list with configuration settings}
#'   \item{par}{A list with initial parameters}
#'   \item{map}{A list with the mapping of parameters}
#' }
#'
#' @source simulated with admove
"montagus_harrier"
