source(file.path("ode_solvers.R"))
source(file.path("tree_generator.R"))

# TODO need to decide data input standards and implement validation
sampled_states_params_csv_path <- file.path("three_sampled_states_params.csv")
params_df <- read.csv(sampled_states_params_csv_path)

#Transition matrix
sampled_states_q_csv_path <- file.path("three_sampled_states_q.csv")
q_matrix <- as.matrix(read.csv(sampled_states_q_csv_path))
rownames(q_matrix) <- q_matrix[, 1]
q_matrix <- q_matrix[, -c(1)]
diag(q_matrix) <- NA
class(q_matrix) <- "numeric"

# TODO test input is a binomial tree
phy <- get_rand_phy(3, 20, params_df, q_matrix)
# TODO temporarily for comparison purposes
nstate <- nrow(params_df)
plot(history.from.sim.discrete(phy, 1:nstate), phy)
tiplabels(frame = "circle", cex = 0.5)
nodelabels(frame = "circle", cex = 0.5)

# Total number of nodes == number of non-leaf nodes * 2 + 1
nnode <- phy[["Nnode"]] * 2 + 1
# Df with node ids, with leaf node ids coming first. Will also be populated by
# time distances from present day and parent/children connections.
topology_df <- data.frame(
  id = seq_len(nnode),
  t_root = NA,
  left = NA,
  right = NA,
  parent = NA)

# Populate topology df with time distances from present day
node_depths <- node.depth.edgelength(phy)
max_depth <- max(node_depths)
invisible(lapply(seq_len(length(node_depths)), function(i) {
  topology_df$t_root[topology_df$id == i] <<- max_depth - node_depths[i]
}))

# Populate topology df with parent/children connections
post_order_edges <- reorder.phylo(phy, "postorder")[["edge"]]
invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
  node <- post_order_edges[i, 1]
  left <- post_order_edges[i, 2]
  right <- post_order_edges[i + 1, 2]

  topology_df$parent[topology_df$id == left] <<- node
  topology_df$parent[topology_df$id == right] <<- node
  topology_df$left[topology_df$id == node] <<- left
  topology_df$right[topology_df$id == node] <<- right
}))

# To be populated with state likelihoods used in backwards time equations.
# list[[x]][[y]] is the likelihood for state y in node x. Note: likelihoods for
# leaf nodes are sampling probabilities, but the internal node likelihoods must
# be calculated.
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
  likelihoods <- get_backwards_likelihoods(left_likelihoods, right_likelihoods,
                                           left_t0, right_t0, tf,
                                           params_df, q_matrix)
  backwards_likelihoods_list[[node]] <<- likelihoods
}))

# To be populated with ancestral state reconstruction probabilities.
# list[[x]][[y]] is the probability of state y in node x. Note: this will also
# include leaf nodes, which have probabilities == 1 for the observed state.
state_probabilities_list <- rep(list(rep(0, nstate)), nnode)

# Populate leaf node state probabilities
invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
  state <- phy[["tip.state"]][[i]]
  state_freq <- params_df$freq[params_df$state == state]
  state_probabilities_list[[i]][[state]] <<- 1
}))

# Populate root node state probabilities. Root node ID == number of leaf
# nodes + 1 == number of internal nodes + 2.
root_node <- phy[["Nnode"]] + 2
state_probabilities_list[[root_node]] <- (
  backwards_likelihoods_list[[root_node]] * params_df$freq
  / (sum(backwards_likelihoods_list[[root_node]] * params_df$freq))
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
                                           params_df, q_matrix)
  state_probabilities_list[[node]] <<- (
    backwards_likelihoods_list[[node]] * likelihoods
    / sum(backwards_likelihoods_list[[node]] * likelihoods)
  )
}))

# Generate log_df to compare internal node results with diversitree states
internal_node_ids <- root_node:nnode
log_df <- data.frame(
  id = internal_node_ids,
  diversitree_state = phy[["node.state"]],
  match = (
    phy[["node.state"]]
    == unlist(lapply(state_probabilities_list[internal_node_ids], which.max))
  ),
  state_probabilities = unlist(
    lapply(state_probabilities_list[internal_node_ids], function(e) {
      return(paste(e, collapse = ", "))
    })
  )
)
