#' Test package.
#'
#' @param phy TODO has tip.state
#' @param params_file TODO
#' @param q_matrix_file TODO
#' @return TODO phylo
#' @export
get_tree <- function(phy, params_file, q_matrix_file) {
  # TODO validate parameters
  params_df <- read.csv(file.path(params_file))
  q_matrix <- get_q_matrix(q_matrix_file)

  nstate <- nrow(params_df)
  nleaf_node <- phy[["Nnode"]]
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- nleaf_node * 2 + 1

  node_depths <- ape::node.depth.edgelength(phy)
  max_depth <- max(node_depths)
  post_order_edges <- ape::reorder.phylo(phy, "postorder")[["edge"]]
  topology_df <- get_topology_df(nnode,
                                 node_depths,
                                 max_depth,
                                 post_order_edges)

  backwards_likelihoods_list <- get_backwards_likelihoods_list(params_df,
                                                               q_matrix,
                                                               nstate,
                                                               nnode,
                                                               post_order_edges,
                                                               topology_df)

  return("Hello world!")
}

#' TODO
#'
#' @param q_matrix_csv TODO
#' @return TODO
#' @noRd
get_q_matrix <- function(q_matrix_file) {
  q_matrix <- as.matrix(read.csv(file.path(q_matrix_file)))
  rownames(q_matrix) <- q_matrix[, 1]
  q_matrix <- q_matrix[, -c(1)]
  diag(q_matrix) <- NA
  class(q_matrix) <- "numeric"
  return(q_matrix)
}

#' TODO
#'
#' @param nnode TODO
#' @param node_depths TODO
#' @param max_depth TODO
#' @param post_order_edges TODO
#' @return TODO
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
#' @param params_df TODO
#' @param q_matrix TODO
#' @param nstate TODO
#' @param nnode TODO
#' @param post_order_edges TODO
#' @param topology_df TODO
#' @return List of state likelihoods used in backwards time equations.
#' list[[x]][[y]] is the likelihood for state y in node x.
#' @noRd
get_backwards_likelihoods_list <- function(params_df,
                                           q_matrix,
                                           nstate,
                                           nnode,
                                           post_order_edges,
                                           topology_df) {
  backwards_likelihoods_list <- rep(list(rep(0, nstate)), nnode)

  # Populate leaf node likelihoods for backwards time equations
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_freq <- params_df$freq[params_df$state == state]
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
                                             params_df, q_matrix)
    backwards_likelihoods_list[[node]] <<- likelihoods
  }))

  return(backwards_likelihoods_list)
}
