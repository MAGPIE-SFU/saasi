#' Ebola tree example
#'
#' An ebola tree with 1493 tips and 3 states used as example. 
#'
#' @format A `phylo` tree object 
#' @source Downloaded from Nextstrain 
#' @examples
#' data(ebola_tree)
#' ebola_tree$root.edge
"ebola_tree"



#' Simulated tree example
#'
#' A simulated birth-death-sampling tree with 13 tips and 3 states used for demonstration
#' purposes. 
#' The dataset includes the a tree, transition rate matrix, 
#' and birth-death-sampling parameters used for simulation.
#'
#' @format A list containing three objects:
#' \describe{
#'   \item{demo_tree_prepared}{A `phylo` tree object with tip states}
#'   \item{demo_Q}{A 3x3 transition rate matrix}
#'   \item{demo_pars}{A data.frame with birth-death-sampling parameters including
#'     rootprior, lambda, mu, and psi for each state}
#' }
#' @source Simulated tree
#' @examples
#' data(saasi_demo)
#' plot(demo_tree_prepared)
#' demo_Q
#' demo_pars
"saasi_demo"