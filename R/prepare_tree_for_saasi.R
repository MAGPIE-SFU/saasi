#' Reformat a phylogenetic tree to be compatible with `saasi`
#' 
#' `prepare_tree_for_saasi()` will take in a `phylo` object and reformat it so that the output is compabtible
#' with the `saasi` function. This reformatting includes:
#' - Adding tip annotations
#' - Resolving polytomies
#' - Imposing a minimum branch length
#' - Removing undesired states
#' 
#' The reformatting is conducted as follows:
#' - Polytomies are resolved using [ape::multi2di()].
#' - The length of any branches with zero or negative length is set to `min_branch_length`.
#' - All tips with an `NA` value for their state annotation are removed in addition to any states specified in `drop_states`.
#'
#' @param tree An object of class `phylo`.
#' @param tip_data Either a named vector or a `data.frame`. 
#' If a named vector, the elements of the vector must be the state of each tip and the names must be the labels of those tips.
#' If a `data.frame`, there must be 2 columns: the first with the tip labels and the second with the tip states.
#' This can be omitted if `tree` already has a `tip.states` attribute.
#' @param resolve_polytomies Boolean. Indicates whether polytomies should be resolved. Default value is `TRUE`. 
#' Setting to `FALSE` will produce a tree that is incompatible with `saasi`.
#' @param fix_branches Boolean. Indicates whether zero and negative branch lengths should be fixed. Default value is `TRUE`.
#' Setting to `FALSE` will produce a tree that is incompatible with `saasi`.
#' @param min_branch_length Numeric. The minimum branch length for the resultant tree. Any short branch lengths will be increased accordingly. Default value is `1e-5`.
#' @param drop_states An optional character vector of states to remove from the tree.
#' @return An object of class `phylo` that is compatible with the `saasi()` function.
#' @seealso [add_tip_states()] to manually add state annotations to tip states or 
#' [drop_tips_by_state()] to manually remove tips from the tree based on their state.
#' 
#' @examples
#' # Check if the demo tree is compatible
#' check_tree_compatibility(demo_tree)
#' 
#' # Since tip annotations are missing, we need state data for the tips
#' head(demo_metadata)
#' 
#' # Reformat the demo tree with the tip states
#' tree_prepared <- prepare_tree_for_saasi(demo_tree, demo_metadata)
#' 
#' # Check the compability of the reformatted tree
#' check_tree_compatibility(tree_prepared)
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
  
  if(!is.null(tip_data)){
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