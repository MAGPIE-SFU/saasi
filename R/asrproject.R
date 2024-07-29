#' Ancestral state reconstruction.
#'
#' Get the internal node state probabilities of a tree with defined leaf states.
#'
#' @param phy A \code{phylo} phylogenetic tree (ape format). Must contain
#' \code{tip.state}.
#' @param params Data frame containing non-q parameters used in ancestral
#' state reconstruction algorithm. Must have the following column names:
#' \code{state}, \code{freq}, \code{lambda}, \code{mu}, \code{psi}. The state
#' values should be a 1-based sequence of natural numbers:
#' \code{1, 2, ..., \emph{n}}.
#'
#' Example:
#' | **state** | **freq** | **lambda** | **mu** | **psi** |
#' | :--- | :--- | :--- | :--- | :--- |
#' | 1 | 0.2 | 3 | 0.02 | 1 |
#' | 2 | 0.3 | 3 | 0.02 | 1 |
#' | 3 | 0.5 | 3 | 0.02 | 1 |
#' @param q_matrix_file Path to CSV file containing q matrix used in ancestral
#' state reconstruction algorithm. Must have column and row headers
#' listing state values in \code{params_file}.
#'
#' Example:
#'
#' \code{,1,2,3}
#'
#' \code{1,,0.3,0.6}
#'
#' \code{2,0.5,,0.1}
#'
#' \code{3,0.2,0.4,}
#' @return A data frame listing the state probabilities of every node in
#' \code{phy}.
#' @export
asr <- function(phy, params, q_matrix_file) {
  # TODO validate parameters
  q_matrix <- get_q_matrix(q_matrix_file)

  nstate <- nrow(params)
  nleaf_node <- phy[["Nnode"]]
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- nleaf_node * 2 + 1
  # Root node ID == number of leaf nodes + 1 == number of internal nodes + 2.
  root_node <- nleaf_node + 2

  node_depths <- ape::node.depth.edgelength(phy)
  max_depth <- max(node_depths)
  post_order_edges <- ape::reorder.phylo(phy, "postorder")[["edge"]]
  topology_df <- get_topology_df(nnode,
                                 node_depths,
                                 max_depth,
                                 post_order_edges)

  backwards_likelihoods_list <- get_backwards_likelihoods_list(phy,
                                                               params,
                                                               q_matrix,
                                                               nstate,
                                                               nnode,
                                                               post_order_edges,
                                                               topology_df)

  state_probabilities_list <- get_state_probabilities_list(
    phy,
    params,
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

  return(state_probabilities_df)
}

#' Convert q matrix CSV file to matrix object.
#'
#' @param q_matrix_csv Path to CSV file containing q matrix used in ancestral
#' state reconstruction algorithm.
#' @return Matrix representation of \code{q_matrix_csv}.
#' @noRd
get_q_matrix <- function(q_matrix_file) {
  q_matrix <- as.matrix(read.csv(file.path(q_matrix_file)))
  rownames(q_matrix) <- q_matrix[, 1]
  q_matrix <- q_matrix[, -c(1)]
  diag(q_matrix) <- NA
  class(q_matrix) <- "numeric"
  return(q_matrix)
}

#' Get data frame representation of tree topology.
#'
#' @param nnode Number of nodes in tree.
#' @param node_depths Depths of nodes in tree.
#' @param max_depth Depth of tree.
#' @param post_order_edges Post order representation of edges in tree.
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

#' TODO
#'
#' @param phy TODO
#' @param params TODO
#' @param q_matrix TODO
#' @param nstate TODO
#' @param nnode TODO
#' @param post_order_edges TODO
#' @param topology_df TODO
#' @return List of state likelihoods used in backwards time equations.
#' list[[x]][[y]] is the likelihood for state y in node x.
#' @noRd
get_backwards_likelihoods_list <- function(phy,
                                           params,
                                           q_matrix,
                                           nstate,
                                           nnode,
                                           post_order_edges,
                                           topology_df) {
  backwards_likelihoods_list <- rep(list(rep(0, nstate)), nnode)

  # Populate leaf node likelihoods for backwards time equations
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_freq <- params$freq[params$state == state]
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
    likelihoods <- get_backwards_likelihoods(left_likelihoods,
                                             right_likelihoods,
                                             left_t0, right_t0, tf,
                                             params, q_matrix)
    backwards_likelihoods_list[[node]] <<- likelihoods
  }))

  return(backwards_likelihoods_list)
}

#' TODO
#'
#' @param phy TODO
#' @param params TODO
#' @param q_matrix TODO
#' @param nstate TODO
#' @param nnode TODO
#' @param root_node TODO
#' @param post_order_edges TODO
#' @param topology_df TODO
#' @param backwards_likelihoods_list TODO
#' @return List of ancestral state probabilities.
#' list[[x]][[y]] is the probability of state y in node x.
#' @noRd
get_state_probabilities_list <- function(phy,
                                         params,
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
    backwards_likelihoods_list[[root_node]] * params$freq
    / (sum(backwards_likelihoods_list[[root_node]] * params$freq))
  )

  # Populate internal node state probabilities
  invisible(lapply(seq(length(post_order_edges[, 1]), 1, -2), function(i) {
    node <- post_order_edges[[i, 1]]
    if (node == root_node) {
      return()
    }

    parent <- topology_df$parent[topology_df$id == node]
    parent_state_probabilities <- state_probabilities_list[[parent]]

    t0 <- topology_df$t_root[topology_df$id == parent]
    tf <- topology_df$t_root[topology_df$id == node]

    likelihoods <- get_forwards_likelihoods(parent_state_probabilities,
                                            t0, tf,
                                            params, q_matrix)
    state_probabilities_list[[node]] <<- (
      backwards_likelihoods_list[[node]] * likelihoods
      / sum(backwards_likelihoods_list[[node]] * likelihoods)
    )
  }))

  return(state_probabilities_list)
}

#' TODO
#'
#' @param phy TODO
#' @param nstate TODO
#' @param state_probabilities_list TODO
#' @return TODO
#' @noRd
get_state_probabilities_df <- function(phy, nstate, state_probabilities_list) {
  state_probabilities_df <-
    as.data.frame(do.call(rbind, state_probabilities_list))
  row.names(state_probabilities_df) <-
    c(phy[["tip.label"]], phy[["node.label"]])
  names(state_probabilities_df) <- seq_len(nstate)
  return(state_probabilities_df)
}
