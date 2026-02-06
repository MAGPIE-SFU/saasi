#' Prepare phylogenetic tree for SAASI
#'
#' @param tree A tree of class `phylo`.
#' @param tip_data 
#'   A named vector where names = tip labels, values = states, or
#'   a dataframe with 2 columns: column 1 = tip labels, column 2 = states.
#'   If tree already contains tip states they should be named `tip.state`.
#' @param resolve_polytomies  Resolve polytomies with [multi2di].
#' @param fix_branches Resolve zero/negative branch lengths.
#' @param min_branch_length Replace value for zero or negative branch lengths (default: 1e-5)
#' @param drop_states States to be removed from the tree. NA is always included.
#'   Add additional values to drop states. Example:c("Not Collected").
#' @return A prepared `phylo` object to be used as input for [saasi] and/or [plot_saasi].
#' @examples
#' data(ebola_tree)
#' # this tree already has tip states 
#' ebola_tree$tip.state
#' ebola_tree_prep <- prepare_tree_for_saasi(ebola_tree)
#' 
#' @export
prepare_tree_for_saasi <- function(tree,
                                   tip_data = NULL,
                                   resolve_polytomies = TRUE,
                                   fix_branches       = TRUE,
                                   min_branch_length  = 1e-5,
                                   drop_states        = character(0)) {
  if(!inherits(tree, "phylo")){
    stop("'tree' must be a phylo object.")
  }

  if(!ape::is.rooted(tree)){
    tree <- ape::root(tree)
  }
  if(is.null(tree$edge.length)){
    stop("Tree must have branch lengths.")
  }
  
  if (is.null(tree$tip.state)) {
  if (!is.null(tip_data)){
    tree <- add_tip_states(tree, tip_data)
    n_total   <- length(tree$tip.label)
    n_na      <- sum(is.na(tree$tip.state))
    if (n_na > 0) {
      message("NA tips:", n_na)
    }
  } 
  else{
    stop("No tip states provided or found on tree.",
            "Attach tip.state.")
  }
  }
  tree <- ape::collapse.singles(tree)
  
  if(resolve_polytomies){
    tree <- ape::multi2di(tree, tol = 1e-8)
  } 

  if(fix_branches){
    n_fix <- sum(tree$edge.length <= 0)
    if(n_fix > 0){
      tree$edge.length[tree$edge.length <= 0] <- min_branch_length
    } 
  } 
  
  if(!is.null(drop_states) && !is.null(tree$tip.state)){
    all_drop <- c(NA, drop_states)
    tree <- drop_tips_by_state(tree, all_drop)
  } 
  return(tree)
}