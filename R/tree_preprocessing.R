#' Read phylogenetic tree from file
#'
#' Automatically detects and reads Newick (.nwk, .tre) or Nexus (.nexus, .nex) 
#' tree files. Includes automatic fixes for Nextstrain formatting.
#'
#' @param filepath Path to tree file
#' @param fix_nextstrain Logical. If TRUE, automatically fixes Nextstrain-specific
#'   formatting issues. Default TRUE.
#' @return A phylo object
#' @export
read_tree_file <- function(filepath, fix_nextstrain = TRUE) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  ext <- tolower(tools::file_ext(filepath))
  
  tree_text <- readLines(filepath, warn = FALSE)
  tree_text <- paste(tree_text, collapse = "")
  
  if (fix_nextstrain) {
    if (grepl("^\\(.*:0\\):0;$", tree_text)) {
      tree_text <- sub("^\\(", "", tree_text)
      tree_text <- sub(":0\\):0;$", ":0;", tree_text)
    }
  }
  
  tmp_file <- tempfile(fileext = paste0(".", ext))
  writeLines(tree_text, tmp_file)
  
  if (ext %in% c("nwk", "tre", "tree", "newick")) {
    tree <- ape::read.tree(tmp_file)
  } else if (ext %in% c("nex", "nexus")) {
    tree <- ape::read.nexus(tmp_file)
  } else {
    stop("Unknown file extension: ", ext, "\n",
         "Supported formats: .nwk, .tre, .tree, .newick, .nex, .nexus")
  }
  unlink(tmp_file)
  return(tree)
}


#' Fix zero or negative branch lengths
#'
#' Replaces zero or negative branch lengths with a small positive value.
#'
#' @param tree A phylo object
#' @param min_length Minimum branch length to use (default: 1e-5)
#' @param replace_zero Logical. Replace zero-length branches (default: TRUE)
#' @param replace_negative Logical. Replace negative branches (default: TRUE)
#' @return A phylo object with fixed branch lengths
#' @export
fix_branch_lengths <- function(tree, 
                               min_length = 1e-5,
                               replace_zero = TRUE,
                               replace_negative = TRUE) {
  
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be a phylo object")
  }
  if (is.null(tree$edge.length)) {
    warning("Tree has no branch lengths. Cannot fix.")
    return(tree)
  }
  if (replace_negative) {
    neg_idx <- which(tree$edge.length < 0)
    if (length(neg_idx) > 0) {
      tree$edge.length[neg_idx] <- min_length
    }
  }
  if (replace_zero) {
    zero_idx <- which(tree$edge.length == 0)
    if (length(zero_idx) > 0) {
      tree$edge.length[zero_idx] <- min_length
    }
  }
  return(tree)
}


#' Attach tip states to phylogenetic tree
#'
#' Adds tip state information to a phylo object.
#'
#' @param tree A phylo object
#' @param tip_data Either:
#'   - A named vector where names = tip labels, values = states
#'   - A data frame with 2 columns: column 1 = tip labels, column 2 = states
#' @param require_all_tips Logical. If TRUE (default), all tree tips must have states.
#' @return A phylo object with tip.state attribute added
#' @export
attach_tip_states <- function(tree, 
                              tip_data,
                              require_all_tips = TRUE) {
  
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be a phylo object")
  }
  
  # Handle named vector
  if (is.atomic(tip_data) && !is.null(names(tip_data))) {
    tip_labels <- names(tip_data)
    states <- as.vector(tip_data)
  } 
  # Handle data frame
  else if (is.data.frame(tip_data)) {
    if (ncol(tip_data) != 2) {
      stop("Data frame must have exactly 2 columns:\n",
           "  Column 1 = tip labels\n",
           "  Column 2 = states")
    }
    tip_labels <- tip_data[[1]]
    states <- tip_data[[2]]
  } 
  else {
    stop("'tip_data' must be either:\n",
         "  - A named vector (names = tip labels, values = states)\n",
         "  - A data frame with 2 columns (column 1 = tip labels, column 2 = states)")
  }
  
  # Check for duplicates in tip_data
  if (any(duplicated(tip_labels))) {
    stop("Duplicate tip labels found in tip_data")
  }
  
  tree_tips <- tree$tip.label
  
  # Check how many tips match
  matched_tips <- tip_labels %in% tree_tips
  n_matched <- sum(matched_tips)
  
  if (n_matched == 0) {
    stop("No matching tip labels found between tree and tip_data.\n",
         "Tree has ", length(tree_tips), " tips.\n",
         "tip_data has ", length(tip_labels), " entries.\n",
         "Example tree tips: ", paste(head(tree_tips, 3), collapse = ", "), "\n",
         "Example tip_data labels: ", paste(head(tip_labels, 3), collapse = ", "))
  }
  
  # Check if all tree tips have states
  missing_in_data <- setdiff(tree_tips, tip_labels)
  if (length(missing_in_data) > 0) {
    if (require_all_tips) {
      stop(length(missing_in_data), " tree tips are missing from tip_data.\n",
           "Missing tips: ", paste(head(missing_in_data, 5), collapse = ", "),
           if (length(missing_in_data) > 5) paste0(" ... and ", length(missing_in_data) - 5, " more"),
           "\n\nSet require_all_tips = FALSE to allow partial matching.")
    } else {
      warning(length(missing_in_data), " tree tips have no states (will be NA)")
    }
  }
  
  # Check if tip_data has labels not in tree
  extra_in_data <- setdiff(tip_labels, tree_tips)
  if (length(extra_in_data) > 0) {
    warning(length(extra_in_data), " labels in tip_data are not in tree (will be ignored)")
  }
  
  # Match states to tree tips
  tip_states <- states[match(tree_tips, tip_labels)]
  tree$tip.state <- tip_states
  
  return(tree)
}

#' Prepare phylogenetic tree for SAASI
#'
#' Validates and prepares a phylogenetic tree for use with SAASI.
#'
#' @param tree A phylo object or path to tree file (newick/nexus)
#' @param tip_data Optional. Either:
#'   - A named vector where names = tip labels, values = states
#'   - A data frame with 2 columns: column 1 = tip labels, column 2 = states
#' @param resolve_polytomies Logical. Automatically resolve polytomies (default: TRUE)
#' @param fix_branch_lengths Logical. Fix zero/negative branch lengths (default: TRUE)
#' @param min_branch_length Minimum branch length for fixing (default: 1e-5)
#' @param require_all_tips Logical. Require all tips to have states (default: TRUE)
#' @return A prepared phylo object ready for SAASI
#' @export
prepare_tree_for_saasi <- function(tree,
                                   tip_data = NULL,
                                   resolve_polytomies = TRUE,
                                   fix_branch_lengths = TRUE,
                                   min_branch_length = 1e-5,
                                   require_all_tips = TRUE) {
  
  if (is.character(tree)) {
    tree <- read_tree_file(tree, fix_nextstrain = TRUE)
  }
  
  if (!inherits(tree, "phylo")) {
    stop("'tree' must be a phylo object or filepath to a tree file")
  }
  
  if (!ape::is.rooted(tree)) {
    stop("Tree must be rooted.\n",
         "Root your tree before using SAASI with ape::root() or phytools::midpoint.root()")
  }
  
  if (is.null(tree$edge.length)) {
    stop("Tree must have branch lengths")
  }
  
  # Resolve polytomies BEFORE fixing branch lengths
  # (multi2di creates zero-length branches)
  if (!ape::is.binary(tree)) {
    if (resolve_polytomies) {
      tree <- ape::multi2di(tree)
    } else {
      stop("Tree must be fully binary (no polytomies).\n",
           "Use resolve_polytomies = TRUE to auto-resolve or use ape::multi2di()")
    }
  }
  
  # Fix branch lengths AFTER polytomy resolution
  if (fix_branch_lengths) {
    n_zero <- sum(tree$edge.length == 0)
    n_neg <- sum(tree$edge.length < 0)
    if (n_zero > 0 || n_neg > 0) {
      tree <- fix_branch_lengths(tree, min_length = min_branch_length)
    }
  }
  
  if (!is.null(tip_data)) {
    tree <- attach_tip_states(tree, tip_data, require_all_tips)
  } else if (is.null(tree$tip.state)) {
    warning("No tip states found. Add tip.state to the tree before using SAASI.")
  }
  
  return(tree)
}

#' Check if tree is compatible with SAASI
#'
#' Quick validation function that checks tree compatibility without modification.
#'
#' @param tree A phylo object
#' @param verbose Logical. Print detailed messages (default: TRUE)
#' @return Logical. TRUE if tree is compatible, FALSE otherwise
#' @export
check_tree_compatibility <- function(tree, verbose = TRUE) {
  
  if (!inherits(tree, "phylo")) {
    if (verbose) message("Not a phylo object")
    return(FALSE)
  }
  
  checks <- list(
    rooted = ape::is.rooted(tree),
    binary = ape::is.binary(tree),
    has_edge_lengths = !is.null(tree$edge.length),
    has_tip_states = !is.null(tree$tip.state),
    no_zero_branches = if(!is.null(tree$edge.length)) all(tree$edge.length > 0) else FALSE,
    no_negative_branches = if(!is.null(tree$edge.length)) all(tree$edge.length >= 0) else FALSE
  )
  
  if (verbose) {
    message("Tree compatibility check:")
    message("  Rooted: ", if(checks$rooted) "True" else "False")
    message("  Binary: ", if(checks$binary) "True" else "False")
    message("  Has branch lengths: ", if(checks$has_edge_lengths) "True" else "False")
    message("  Has tip states: ", if(checks$has_tip_states) "True" else "False")
    if (checks$has_edge_lengths) {
      message("  No zero-length branches: ", if(checks$no_zero_branches) "True" else "False")
      message("  No negative branches: ", if(checks$no_negative_branches) "True" else "False")
    }
  }
  
  all_pass <- checks$rooted && checks$binary && checks$has_edge_lengths && checks$has_tip_states
  
  if (verbose) {
    if (all_pass) {
      message("\n Tree is compatible with SAASI")
    } else {
      message("\n Tree needs preparation. Use prepare_tree_for_saasi()")
    }
  }
  
  return(all_pass)
}