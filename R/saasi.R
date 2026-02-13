#' Sampling-Aware Ancestral State Inference
#'
#' Reconstructs ancestral states for internal nodes of a phylogenetic tree while 
#' accounting for heterogeneous sampling rates across states or locations.
#'
#' @param phy A `phylo` object containing the phylogenetic tree. The tree must be 
#'   rooted, binary, and contain `tip.state`. 
#' @param Q A numeric transition rate matrix (n × n) where n is the number 
#'   of states. Row and column names represent states. Off-diagonal 
#'   elements are transition rates; diagonal elements are set such that rows sum to zero.
#' @param pars A data frame with the other parameters used in the ancestral state reconstruction algorithm. 
#' Must have the following column names: state, prior, lambda, mu, and psi. 
#' The prior values refer to the baseline probabilities of the states (used at the root of the tree).
#'   
#' @return A data frame with state probabilities for each internal node in `phy`. 
#'   
#' @examples
#' \dontrun{
#' # Run saasi with equal sampling rates
#' result <- saasi(phy = tree, 
#'                 q_matrix = Q,
#'                 lambda = rates$lambda,
#'                 mu = rates$mu,
#'                 psi = rates$psi,
#'                 prior = NULL)
#' 
#' # Run saasi with different sampling rates per state
#' result2 <- saasi(phy = tree,
#'                  q_matrix = Q,
#'                  lambda = rates$lambda,
#'                  mu = rates$mu,
#'                  psi = c(rates$psi, rates$psi/2, rates$psi),
#'                  prior = NULL)
#' }
#' 
#' @export
saasi <- function(phy, Q, pars) {
  ## INPUT CHECKS --------------------------------------------------------------
  
  ## phy
  # Check compatibility of phy with SAASI
  if( !check_tree_compatibility(phy) ) {
    stop("`phy` is not compatible with `saasi`. Please use the `prepare_tree_for_saasi` function to reformat your `phylo` object.")
  }
  
  # Check if the tree is likely not a timed tree
  # - If the median branch length is less than 1e-3, flag it as being potentially not a timed tree
  if( quantile(phy$edge.length, probs=0.5) <= 1e-3 )
    warning("The median branch length in your tree is less than 0.001. Please double-check that your tree is a timed tree.")
  
  
  ## Q
  # Check that q_matrix is compatible, it must:
  # - be matrix with no NAs
  # - be a square matrix
  # - have rows summing to 0 (a threshold of 1e-10 is used to account small numerical issues)
  # - both rows and columns must have names
  if( any(is.na(Q)) ) {
    stop("Q contains NAs.")
  }
  if( !is.matrix(Q)) {
    stop("Q must be a matrix. Please consult the documentation for details.")
  }
  if( nrow(Q) != ncol(Q) ) {
    stop(paste0("The transition rate matrix must be square. Current input has ", nrow(Q), " rows and ", ncol(Q), " columns."))
  }
  if( any(abs(rowSums(Q)) > 1e-10) ) {
    stop("The rows of the transition rate matrix must sum to 0.")
  }
  if( is.null(rownames(Q)) || is.null(colnames(Q)) ) {
    stop("Both rows and columns of the transition rate matrix should have names corresponding to the traits.")
  }
  
  #  If Q is entered as NULL, estimate them using the default values
  if( is.null(Q) ){
    warning("No transition rate matrix provided, estimating with ace using an equal rates matrix structure. If the user desires an alternative matrix structure refer to the estimate_transition_rates function.")
    Q <- estimate_transition_rates(phy)
  }
 
  
  ## pars
  # Input should be a data.frame with no NAs
  if( any(is.na(pars)) )
    stop("pars contains NAs.")
  if( !is.data.frame(pars))
    stop("pars must be a data.frame. Please consult the documentation for details.")
  
  # Check that pars has all of the right column names
  par_names <- c("state", "lambda", "mu", "psi")
  missing_pars <- par_names[!(par_names %in% names(pars))] # vector of expected parameters that are missing in the input dataframe
  if( length(missing_pars) > 0 ){
    stop("The following columns are missing from pars: ", paste(missing_pars, sep=","))
  }
  
  # If no prior is provided, add a uniform prior
  if( !("root_prior" %in% names(pars)) ) {
    pars$root_prior <- rep(1, nrow(pars)) / nrow(pars)
  } else if( (abs(sum(pars$root_prior) - 1) > 1e-5) || min(pars$root_prior) < 0 ) { # check that root_prior is a valid probability
    stop("Please ensure that root_prior is a valid probability: all entries must be non-negative and sum to 1.")
  }
  
  # Check birth-death-sampling rates
  # - mu: all must be non-negative
  # - lambda: all must be non-negative and at least one must be positive
  # - psi: all must be non-negative and all observed states must be positive
  with(pars, {
    if( any(mu < 0) ) {
      stop("All death rates (mu) must be non-negative.")
    }
    if( all(!is.null(psi)) ) {
      if( any(psi < 0) ) {
        stop("All sampling rates (psi) must be non-negative.")
      } else if( any(sort(colnames(Q)[sum(psi > 0) > 0]) != sort(unique(phy$tip.state))) ) {
        stop("Sampling rate (psi) must be positive for all observed states.")
      }
    }
    if( all(!is.null(lambda)) ) {
      if( any(lambda < 0) ) {
        stop("All birth rates (lambda) must be non-negative.")
      } else if( sum(lambda > 0) == 0 ) {
        stop("At least one birth rate (lambda) must be positive.")
      }
    }
  })
  
  # Check that pars & Q have the same states
  # - phy can have missing states, but cannot have states not present in Q or pars
  if( any(sort(colnames(Q)) != sort(pars$state)) ){
    stop("The states specified in Q and pars do not match.")
  }
  if( any(!(unique(phy$tip.state) %in% pars$state)) ) {
    stop("There is at least one state in phy that is not in pars or Q.")
  }
  

  ## PRE-PROCESSING ------------------------------------------------------------
  # reformat input and instantiate derived quantities prior to executing main algorithm

  # Convert states to numerical values and preserve state names in separate column/attribute
  if( !is.numeric(phy$tip.state) ) {
    phy$tip.statename <- phy$tip.state
    phy$tip.state <- as.numeric(factor(phy$tip.state))
    pars$statename <- pars$state
    pars$state <- as.numeric(factor(pars$state))
    # Reorder pars and Q so that row/column corresponds to the numeric state value
    pars <- pars[order(pars$state),]
    Q <- Q[order(as.numeric(factor(rownames(Q)))), order(as.numeric(factor(colnames(Q))))]
  }
  
  # if(!is.numeric(phy$tip.state)){
  #   # stop("Convert tip state to numeric values.")
  #     if( !all(sort(rownames(Q))==sort(params_df$state))) {
  #         stop("State names and q names need to agree and be present only once") 
  #     } else {
  #         phy$tip.statename = phy$tip.state
  #         phy$tip.state = as.numeric(factor(phy$tip.state)) # makes it numeric 
  #     params_df$statename = params_df$state
  #     params_df$state = as.numeric(factor(params_df$statename)) # numeric state
  #     params_df = params_df[order(params_df$state), ]
  #     # reorder q if necessary, so its order corresponds to the numeric state 
  #     Q = Q[ order(as.numeric(factor(rownames(Q)))),
  #                          order(as.numeric(factor(colnames(Q))))]
  #     }
  # }
  

  # Check that 'state' is a sequence of 1-based natural numbers
  # if (!is.numeric(params_df$state)) {
  #   stop("'state' column must be a sequence of numerical values (1, 2, ..., n).")
  # }
  
  # Checking the number of unique states in tip.state matches the number of states
  # in params_df$state
  
  # if(!length(unique(phy$tip.state)) == length(params_df$state)){
  #   stop("The number of states does not match.")
  # }
  
  nstate <- nrow(pars)
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- phy[["Nnode"]] * 2 + 1
  # Root node ID == number of leaf nodes + 1 == number of internal nodes + 2.
  root_node <- phy[["Nnode"]] + 2
  
  node_depths <- ape::node.depth.edgelength(phy)
  max_depth <- max(node_depths)
  post_order_edges <- ape::reorder.phylo(phy, "postorder")[["edge"]]
  
  topology_df <- get_topology_df(nnode,
                                 node_depths,
                                 max_depth,
                                 post_order_edges)
  
  
  ## Run SAASI -----------------------------------------------------------------
  backwards_likelihoods_list <- get_backwards_likelihoods_list(phy,
                                                               pars,
                                                               Q,
                                                               nstate,
                                                               nnode,
                                                               post_order_edges,
                                                               topology_df)
  
  state_probabilities_list <- get_state_probabilities_list(
    phy,
    pars,
    Q,
    nstate,
    nnode,
    root_node,
    post_order_edges,
    topology_df,
    backwards_likelihoods_list
  )
  
  # if (plot) {
  #   # https://stackoverflow.com/a/17735894
  #   highest_likelihoods <-
  #     apply(state_probabilities_df, 1, function(e) which(e == max(e)))
  #   
  #   # https://colorbrewer2.org/?type=qualitative&scheme=Set3&n=12
  #   # Will serve up to n states == len of vector; afterwards recycles colors
  #   state_colors <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
  #                     "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
  #                     "#ccebc5", "#ffed6f")
  #   # Map likelihoods to colors
  #   highest_likelihood_colors <-
  #     sapply(highest_likelihoods,
  #            function(e) state_colors[[e %% (length(state_colors) + 1)]])
  #   
  #   ape::plot.phylo(phy, label.offset = 0.05)
  #   ape::tiplabels(highest_likelihoods[1:(phy[["Nnode"]] + 1)],
  #                  frame = "circle",
  #                  cex = cex,
  #                  bg = highest_likelihood_colors)
  #   ape::nodelabels(highest_likelihoods[-(1:(phy[["Nnode"]] + 1))],
  #                   frame = "circle",
  #                   cex = cex,
  #                   bg = highest_likelihood_colors)
  # }
  # 
  return(state_probabilities_list)
}

#' Get data frame representation of tree topology.
#'
#' @return Data frame showing the parents, children, and depth of nodes.
#' @noRd
get_topology_df <- function(nnode, node_depths, max_depth, post_order_edges) {
  topology_df <- data.frame(id = seq_len(nnode),
                            t_root = NA,
                            left = NA,
                            right = NA,
                            parent = NA)
  
  # Populate topology df with time distances from present day
  invisible(lapply(seq_len(length(node_depths)), function(i) {
    topology_df$t_root[topology_df$id == i] <<- max_depth - node_depths[i]
  }))
  
  # Populate topology df with parent/children connections
  invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
    node <- post_order_edges[i, 1]
    left <- post_order_edges[i, 2]
    right <- post_order_edges[i + 1, 2]
    
    topology_df$parent[topology_df$id == left] <<- node
    topology_df$parent[topology_df$id == right] <<- node
    topology_df$left[topology_df$id == node] <<- left
    topology_df$right[topology_df$id == node] <<- right
  }))
  
  return(topology_df)
}

#' Get list of likelihoods calculated during backwards traversal of ancestral
#' state reconstruction algorithm.
#'
#' @return List of state likelihoods used in backwards time equations.
#' `list[[x]][[y]]` is the likelihood for state y in node x.
#' @noRd
get_backwards_likelihoods_list <- function(phy,
                                           params_df,
                                           Q,
                                           nstate,
                                           nnode,
                                           post_order_edges,
                                           topology_df) {
  backwards_likelihoods_list <- rep(list(rep(0, nstate)), nnode)
  
  # Populate leaf node likelihoods for backwards time equations
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_freq <- params_df$psi[params_df$state == state]
    backwards_likelihoods_list[[i]][[state]] <<- state_freq
  }))
  
  # Populate internal node likelihoods for backwards time equations
  invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
    node <- post_order_edges[[i, 1]]
    left <- post_order_edges[[i, 2]]
    right <- post_order_edges[[i + 1, 2]]
    
    tf <- topology_df$t_root[topology_df$id == node]
    left_t0 <- topology_df$t_root[topology_df$id == left]
    right_t0 <- topology_df$t_root[topology_df$id == right]
    
    left_likelihoods <- backwards_likelihoods_list[[left]]
    right_likelihoods <- backwards_likelihoods_list[[right]]
    likelihoods <- get_backwards_likelihoods(abs(left_likelihoods),
                                             abs(right_likelihoods),
                                             left_t0, right_t0, tf,
                                             params_df, Q)
    #abs_likelihoods <- abs(likelihoods)
    #norm_likelihoods <- abs_likelihoods / sum(abs_likelihoods)
    backwards_likelihoods_list[[node]] <<- likelihoods # note changed by CC 
  }))
  
  return(backwards_likelihoods_list)
}

#' Get list of probabilities ultimately generated by ancestral state
#' reconstruction algorithm.
#'
#' @return List of ancestral state probabilities.
#' `list[[x]][[y]]` is the probability of state y in node x.
#' @noRd
get_state_probabilities_list <- function(phy,
                                         params_df,
                                         Q,
                                         nstate,
                                         nnode,
                                         root_node,
                                         post_order_edges,
                                         topology_df,
                                         backwards_likelihoods_list) {
  state_probabilities_list <- rep(list(rep(0, nstate)), nnode)
  
  # Populate leaf node state probabilities
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_probabilities_list[[i]][[state]] <<- 1
  }))
  
  state_probabilities_list[[root_node]] <- (
    backwards_likelihoods_list[[root_node]] * params_df$root_prior
    / (sum(backwards_likelihoods_list[[root_node]] * params_df$root_prior))
  )
  
  # Populate internal node state probabilities
  invisible(lapply(seq(length(post_order_edges[, 1]), 1, -2), function(i) {
    node <- post_order_edges[[i, 1]]
    if (node == root_node) {
      return()
    }
    
    parent <- topology_df$parent[topology_df$id == node]
    parent_state_probabilities <- (
      state_probabilities_list[[parent]] * backwards_likelihoods_list[[node]]
    )
    abs_parent_state_probabilities <- abs(parent_state_probabilities)
    norm_probabilities <- (
      params_df$lambda*abs_parent_state_probabilities / sum(params_df$lambda*abs_parent_state_probabilities)
    )
    
    t0 <- topology_df$t_root[topology_df$id == parent]
    tf <- topology_df$t_root[topology_df$id == node]

    likelihoods <- get_forwards_likelihoods(norm_probabilities,
                                            t0, tf,
                                            params_df, Q)
    #abs_likelihoods <- abs(likelihoods)
    state_probabilities_list[[node]] <<- (
      backwards_likelihoods_list[[node]] * likelihoods
      / sum(backwards_likelihoods_list[[node]] * likelihoods)
    )
  }))
  internal_node_ids <- root_node:nnode
  
  node_prob <- matrix(unlist(state_probabilities_list[internal_node_ids]),ncol=nstate,byrow = TRUE)
  node_lik <- as.data.frame(node_prob,row.names = internal_node_ids)
  colnames(node_lik) = c(1:nstate)
  if(!is.null(phy$tip.statename)){
      colnames(node_lik) = levels(factor(phy$tip.statename))
  }
  return(node_lik)
}

#' Convert state probabilities list to data frame.
#'  List of state probabilities.
#' Data frame of state probabilities in each node in phy.
#' 
# get_state_probabilities_df <- function(phy, nstate, state_probabilities_list) {
#   state_probabilities_df <-
#     as.data.frame(do.call(rbind, state_probabilities_list))
#   row.names(state_probabilities_df) <-
#     c(phy[["tip.label"]], phy[["node.label"]])
#   names(state_probabilities_df) <- seq_len(nstate)
#   return(state_probabilities_df)
# }
