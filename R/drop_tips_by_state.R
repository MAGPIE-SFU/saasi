#' Drop tips from a phylogenetic tree by their state values
#'
#' `drop_tips_by_state()` will remove from a `phylo` object all tips in a specified state(s).
#'
#' @param tree An object of class `phylo` with a `tip.state` attribute. 
#' See [add_tip_states()] for adding tip annotations to a `phylo` object.
#' @param drop_values Character vector of states to remove.
#' @return An object of class `phylo` with the specified states removed.
#' @seealso [add_tip_states()] to add tip annotations
#' @examples
#' # This demo tree has 13 tips in 3 possible states
#' demo_tree_prepared
#' table(demo_tree_prepared$tip.state)
#' 
#' # Remove all tips in state 3
#' new_tree <- drop_tips_by_state(demo_tree_prepared, "3")
#' 
#' # The resulting tree has 8 tips, having removed all tips in state 3
#' new_tree
#' table(new_tree$tip.state)
#' @export
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