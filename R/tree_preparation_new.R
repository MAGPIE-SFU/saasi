#' Add tip states to a phylogenetic tree
#'
#'`add_tip_states()` annotates the tips of a `phylo` object with state values.
#'
#' A `tip.states` attribute is added to `tree` based on the `tip_data` provided by the user. 
#' This function is used within the `prepare_tree_for_saasi()` function.
#' @param tree An object of class `phylo`.
#' @param tip_data Either a named vector or a `data.frame`. 
#' If a named vector, the elements of the vector must be the state of each tip and the names must be the labels of those tips.
#' If a `data.frame`, there must be 2 columns: the first with the tip labels and the second with the tip states.
#'
#' @return An object of class `phylo` with annotated tip states.
#' @seealso [drop_tips_by_state()] to remove specific traits or [prepare_tree_for_saasi()] for reformatting `phylo` objects to be `saasi`-compatible.
#' 
#' @examples
#' demo_tree
#' head(demo_metadata)
#' 
#' add_tip_states(demo_tree, demo_metadata)
#' 
#' @export
add_tip_states <- function(tree, tip_data) {
  
  if(is.atomic(tip_data) && !is.null(names(tip_data))){
    tip_labels <- names(tip_data)
    states <- as.vector(tip_data)
  } 
  else if(is.data.frame(tip_data)){
    if (ncol(tip_data) != 2) {
      stop("Data frame must have exactly 2 columns:\n",
           "column 1 = tip labels\n",
           "column 2 = states")
    }
    tip_labels <- as.character(tip_data[[1]])
    states <- as.character(tip_data[[2]])
  } 
  else{
    stop("'tip_data' must be either:\n",
         "A named vector\n",
         "A data frame with 2 columns (column 1 = tip labels, 
         column 2 = states)")
  }
  if(any(duplicated(tip_labels))){
    stop("Duplicate tip labels found in tip_data \n",
         "Each tip label must appear exactly once.")
  }
  tree_tips <- tree$tip.label
  
  n_matched <- sum(tree_tips %in% tip_labels)
  if (n_matched == 0) {
    stop("Tree$tip.label and tip_data are not matched.")
  }
  tip_states <- states[match(tree_tips, tip_labels)]
  tree$tip.state <- tip_states
  
  return(tree)
}

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