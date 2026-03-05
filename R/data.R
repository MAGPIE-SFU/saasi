#' SAASI demonstration: timed phylogenetic tree
#'
#' A simulated birth-death-sampling tree with 13 tips.
#' The tip states are located separately in the [demo_metadata] file.
#'
#' @format An object of class `phylo`
#' @examples
#' demo_tree
"demo_tree"

#' SAASI demonstration: pre-prepared tree
#'
#' A simulated birth-death-sampling tree with 13 tips and 3 states.
#' The tree has pre-prepared to be compatible with `saasi`.
#'
#' @format An object of class `phylo` with `tip.state` attribute
#' @examples
#' data(demo_tree_prepared)
#' demo_tree_prepared$tip.state
"demo_tree_prepared"

#' SAASI demonstration: transition rate matrix
#'
#' A 3x3 transition rate matrix used to simulate discrete states in the demonstration
#' tree. This matrix defines the instantaneous rates of transition between states.
#'
#' @format A 3x3 numeric matrix
#' @examples
#' data(demo_Q)
#' demo_Q
#' # Check that rows sum to 0
#' rowSums(demo_Q)
"demo_Q"

#' SAASI demonstration: birth-death-sampling process parameters
#'
#' A `data.frame` specifying the parameters of the birth-death-sampling process that is used to
#' simulate the demonstration tree. These parameters control diversification
#' and sampling rates for each state.
#'
#' @format A data frame with 3 rows (one per state) and 5 columns:
#' \describe{
#'   \item{`state`}{The character label for each state}
#'   \item{`rootprior`}{Prior probability distribution over states for the root node}
#'   \item{`lambda`}{Numeric speciation rate for each state}
#'   \item{`mu`}{Numeric extinction rate for each state}
#'   \item{`psi`}{Numeric Sampling rate for each state}
#' }
#' @examples
#' data(demo_pars)
#' demo_pars
"demo_pars"



#' SAASI demonstration: tree metadata
#'
#' A `data.frame` linking the tip nodes to their state
#'
#' @format A `data.frame` with 13 rows (one per tip) and 2 columns
#' \describe{
#'    \item{`node`}{Node labels for each tip}
#'    \item{`state`}{Character state for each tip}
#' }
#' @examples
#' data(demo_metadata)
#' demo_metadata
#'
"demo_metadata"

