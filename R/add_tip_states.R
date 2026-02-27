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