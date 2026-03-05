#' Check if phylogenetic tree is compatible with `saasi`
#'
#'@description
#' Verifies the compatibility of a `phylo` object with the [saasi()] function without modifying the object.
#' For `tree` to be compatible with `saasi()`, it must:
#' - be rooted,
#' - have the correct number of internal nodes (the number of tips - 1),
#' - have no unary nodes or polytomies,
#' - have positive branch lengths, and
#' - have valid tip state annotations for all tips.
#' 
#' If the tree is not compatible with `saasi()`, the user can run [prepare_tree_for_saasi()] to resolve the issues.
#'
#' @param tree An object of class `phylo`.
#' @return Logical. Returns `TRUE` if the tree is compatible with `saasi` and `FALSE` otherwise.
#' @examples
#' # Check if this demo tree is compatible with saasi:
#' check_tree_compatibility(demo_tree)
#' 
#' # Now check another tree which has been properly prepared:
#' check_tree_compatibility(demo_tree_prepared)
#' @export
check_tree_compatibility <- function(tree) {
  
  if(!inherits(tree, "phylo")){
    stop("'tree' must be a phylo object.")
  }
  
  n_tips <- length(tree$tip.label)
  n_internal <- tree$Nnode
  
  # check binary
  child_counts <- tabulate(tree$edge[, 1], nbins = n_tips + n_internal)
  internal_counts <- child_counts[(n_tips + 1):(n_tips + n_internal)]
  
  is_binary <- all(internal_counts == 2)
  has_unary <- any(internal_counts == 1)
  has_polytomy <- any(internal_counts > 2)
  
  # check number of internal node == number of tips - 1
  expected_internal <- n_tips - 1
  node_count_ok <- (n_internal == expected_internal)
  
  is_rooted <- ape::is.rooted(tree)
  has_branches <- !is.null(tree$edge.length)
  has_tip_states <- !is.null(tree$tip.state) && length(tree$tip.state) == n_tips
  no_zero      <- has_branches && all(tree$edge.length > 0)
  no_negative  <- has_branches && all(tree$edge.length >= 0)
  no_na_states <- has_tip_states && !any(is.na(tree$tip.state))
  
  all_pass <- is_rooted &&
    is_binary &&
    node_count_ok &&
    has_branches &&
    has_tip_states &&
    no_zero &&
    no_negative &&
    no_na_states
  
    if(all_pass){
      message("Tree is compatible with SAASI")
    }
    else{
      message("Tree is not compatible with SAASI")
      if(!is_rooted){
        message("Tree is not rooted")}
      if(!node_count_ok){
        diff <- expected_internal - n_internal
        if(diff > 0){
          message("Missing internal nodes")
        }else{
          message("Extra internal nodes")
        }
      }
      if(has_unary){
        message("Unary nodes present. Suppress with ape::collapse.branches() or re-run multi2di()")
      } 
      if(has_polytomy){
        message("Polytomies present. Resolve with multi2di(tree, tol = 1e-8)")
      } 
      if(!no_zero){
        message("Zero-length branches present. Fix with: tree$edge.length[tree$edge.length <= 0] <- 1e-5")
      } 
      if(!no_negative){
        message("Negative branches present. Fix with: tree$edge.length[tree$edge.length <= 0] <- 1e-5")
      } 
      if(!no_na_states){
        message("NA tip states present. Remove with drop_tips_by_state(tree, NA)")
      } 
      if(!has_tip_states){
        message("No tip states. Attach with attach_tip_states() or prepare_tree_for_saasi()")
      } 
    }
  return(all_pass)
}