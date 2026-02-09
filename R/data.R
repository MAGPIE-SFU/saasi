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


#' Simulated demonstration tree example
#'
#' A simulated birth-death-sampling tree with 13 tips and 3 states.
#'
#' @format A `phylo` tree object with the following components
#' @examples
#' data(demo_tree_prepared)
#' demo_tree_prepared$tip.state
"demo_tree_prepared"

#' Transition rate matrix for demonstration tree
#'
#' A 3x3 transition rate matrix (Q matrix) used to simulate the demonstration
#' tree. This matrix defines the instantaneous rates of transition between
#' character states.
#'
#' @format A 3x3 numeric matrix
#' @examples
#' data(demo_Q)
#' demo_Q
#' # Check that rows sum to 0
#' rowSums(demo_Q)
"demo_Q"

#' Birth-death-sampling parameters for demonstration tree
#'
#' A data frame containing the birth-death-sampling parameters used to
#' simulate the demonstration tree. These parameters control diversification
#' and sampling rates for each character state.
#'
#' @format A data frame with 3 rows (one per state) and 5 columns:
#' \describe{
#'   \item{state}{Character state names (1, 2, or 3)}
#'   \item{rootprior}{Prior probability of root being in each state (1/3 each)}
#'   \item{lambda}{Speciation rate (1.5 for all states)}
#'   \item{mu}{Extinction rate (0.3 for all states)}
#'   \item{psi}{Sampling rate (0.1 for state 1, 0.5 for states 2 and 3)}
#' }
#' @examples
#' data(demo_pars)
#' demo_pars
#' # Compare sampling rates across states
#' demo_pars$psi
"demo_pars"



#' Metadata for demonstration tree
#' 
#' A data frame containing the node name and tip states
#' 
#' @format A data frame with 13 rows (one per tip) and 2 columns
#' \describe{
#'    \item{node}{node name for each tip}
#'    \item{state}{Character state for each tip (1, 2, or 3)}
#' }
#' @examples
#' data(demo_metadata)
#' demo_metadata
#' demo_metadata$states
#' 
"demo_metadata"

