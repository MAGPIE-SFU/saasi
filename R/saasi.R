#' Sampling Aware Ancestral State Inference
#'
#' Get the internal node state probabilities of a tree with defined leaf states.
#'
#' @param phy A `phylo` phylogenetic tree (ape format). Must contain
#' `tip.state`.
#' @param params_df Data frame containing non-q parameters used in ancestral
#' state reconstruction algorithm. Must have the following column names:
#' `state`, `prior`, `lambda`, `mu`, and `psi`. The state values should be a
#' 1-based sequence of natural numbers: `1`, `2`, ..., *`n`*.
#'
#' Example:
#' | **state** | **prior** | **lambda** | **mu** | **psi** |
#' | :--- | :--- | :--- | :--- | :--- |
#' | 1 | 1/3 | 3 | 0.02 | 1 |
#' | 2 | 1/3 | 3 | 0.02 | 1 |
#' | 3 | 1/3 | 3 | 0.02 | 1 |
#' @param q_matrix Numeric q matrix used in ancestral state reconstruction
#' algorithm. Row and column indices represent states.
#' @param plot If `TRUE`, a plot of `phy` with overlaid state labels will be
#' rendered.
#' @param cex See `cex` parameter in \link{nodelabels}.
#' @return A data frame listing the state probabilities of every node in `phy`.
#' @export
saasi <- function(phy, params_df, q_matrix, plot = FALSE, cex = 1) {
  nstate <- nrow(params_df)
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
  
  backwards_likelihoods_list <- get_backwards_likelihoods_list(phy,
                                                               params_df,
                                                               q_matrix,
                                                               nstate,
                                                               nnode,
                                                               post_order_edges,
                                                               topology_df)
  
  state_probabilities_list <- get_state_probabilities_list(
    phy,
    params_df,
    q_matrix,
    nstate,
    nnode,
    root_node,
    post_order_edges,
    topology_df,
    backwards_likelihoods_list
  )
  
  state_probabilities_df <- get_state_probabilities_df(phy,
                                                       nstate,
                                                       state_probabilities_list)
  
  if (plot) {
    # https://stackoverflow.com/a/17735894
    highest_likelihoods <-
      apply(state_probabilities_df, 1, function(e) which(e == max(e)))
    
    # https://colorbrewer2.org/?type=qualitative&scheme=Set3&n=12
    # Will serve up to n states == len of vector; afterwards recycles colors
    state_colors <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
                      "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
                      "#ccebc5", "#ffed6f")
    # Map likelihoods to colors
    highest_likelihood_colors <-
      sapply(highest_likelihoods,
             function(e) state_colors[[e %% (length(state_colors) + 1)]])
    
    ape::plot.phylo(phy, label.offset = 0.05)
    ape::tiplabels(highest_likelihoods[1:(phy[["Nnode"]] + 1)],
                   frame = "circle",
                   cex = cex,
                   bg = highest_likelihood_colors)
    ape::nodelabels(highest_likelihoods[-(1:(phy[["Nnode"]] + 1))],
                    frame = "circle",
                    cex = cex,
                    bg = highest_likelihood_colors)
  }
  
  return(state_probabilities_df)
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
#' list[[x]][[y]] is the likelihood for state y in node x.
#' @noRd
get_backwards_likelihoods_list <- function(phy,
                                           params_df,
                                           q_matrix,
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
                                             params_df, q_matrix)
    abs_likelihoods <- abs(likelihoods)
    norm_likelihoods <- abs_likelihoods / sum(abs_likelihoods)
    backwards_likelihoods_list[[node]] <<- norm_likelihoods
  }))
  
  return(backwards_likelihoods_list)
}

#' Get list of probabilities ultimately generated by ancestral state
#' reconstruction algorithm.
#'
#' @return List of ancestral state probabilities.
#' list[[x]][[y]] is the probability of state y in node x.
#' @noRd
get_state_probabilities_list <- function(phy,
                                         params_df,
                                         q_matrix,
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
    backwards_likelihoods_list[[root_node]] * params_df$prior
    / (sum(backwards_likelihoods_list[[root_node]] * params_df$prior))
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
      abs_parent_state_probabilities / sum(abs_parent_state_probabilities)
    )
    
    t0 <- topology_df$t_root[topology_df$id == parent]
    tf <- topology_df$t_root[topology_df$id == node]
    
    likelihoods <- get_forwards_likelihoods(norm_probabilities,
                                            t0, tf,
                                            params_df, q_matrix)
    abs_likelihoods <- abs(likelihoods)
    state_probabilities_list[[node]] <<- (
      backwards_likelihoods_list[[node]] * abs_likelihoods
      / sum(backwards_likelihoods_list[[node]] * abs_likelihoods)
    )
  }))
  
  return(state_probabilities_list)
}

#' Convert state probabilities list to data frame.
#'
#'  Data frame of state probabilities in each node in phy.
#'
get_state_probabilities_df <- function(phy, nstate, state_probabilities_list) {
  state_probabilities_df <-
    as.data.frame(do.call(rbind, state_probabilities_list))
  row.names(state_probabilities_df) <-
    c(phy[["tip.label"]], phy[["node.label"]])
  names(state_probabilities_df) <- seq_len(nstate)
  return(state_probabilities_df)
}
