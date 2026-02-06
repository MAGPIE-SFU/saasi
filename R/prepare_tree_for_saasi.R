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


#' @keywords internal
drop_tips_by_state <- function(tree, drop_values) {
  
  if(is.null(tree$tip.state)){
    stop("Tree has no tip.state. Attach states first.")
  }
  
  tips_to_drop <- integer(0)
  
  for(dv in drop_values){
    if (length(dv) == 1 && is.na(dv)) {
      # Drop NA
      tips_to_drop <- c(tips_to_drop, which(is.na(tree$tip.state)))
    } else {
      # Drop tips matching this value
      tips_to_drop <- c(tips_to_drop,
                        which(!is.na(tree$tip.state) & tree$tip.state == dv))
    }
  }
  tips_to_drop <- unique(tips_to_drop)
  if(length(tips_to_drop) == 0){
    return(tree)
  }
  tip_labels_to_drop <- tree$tip.label[tips_to_drop]
  new_tree <- ape::drop.tip(tree, tip_labels_to_drop)
  
  tips_to_keep <- setdiff(1:length(tree$tip.label), tips_to_drop)
  new_tree$tip.state <- tree$tip.state[tips_to_keep]
  return(new_tree)
}