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

# TODO Populate root node state probabilities

# TODO Populate internal node state probabilities

# # Forward-time equations
# invisible(lapply(seq(length(post_order_edges[, 1]), 1, -2), function(i) {
#   node <- post_order_edges[[i, 1]]
#   left <- post_order_edges[[i - 1, 2]]
#   right <- post_order_edges[[i, 2]]
# 
#   backwards_state_1 <- df$backwards_state_1[df$id == node]
#   backwards_state_2 <- df$backwards_state_2[df$id == node]
# 
#   # If root else not root
#   if (i == length(post_order_edges[, 1])) {
#     denominator <- (backwards_state_1 * state_1_freq
#                     + backwards_state_2 * state_2_freq)
#     df$ancestral_state_1[df$id == node] <<-
#       backwards_state_1 * state_1_freq / denominator
#     df$ancestral_state_2[df$id == node] <<-
#       backwards_state_2 * state_2_freq / denominator
#   } else {
#     parent <- df$parent[df$id == node]
#     parent_ancestral_state_1 <- df$ancestral_state_1[df$id == parent]
#     parent_ancestral_state_2 <- df$ancestral_state_2[df$id == parent]
# 
#     t0 <- df$t_root[df$id == parent]
#     tf <- df$t_root[df$id == node]
# 
#     forwards_sol <- get_forwards_sol(parent_ancestral_state_1,
#                                      parent_ancestral_state_2,
#                                      t0, tf,
#                                      λ, μ, Ψ, q, node)
#     forwards_state_1 <- forwards_sol[["forwards_state_1"]]
#     forwards_state_2 <- forwards_sol[["forwards_state_2"]]
# 
#     numerator_1 <- backwards_state_1 * forwards_state_1
#     numerator_2 <- backwards_state_2 * forwards_state_2
#     denominator <- numerator_1 + numerator_2
#     df$ancestral_state_1[df$id == node] <<- numerator_1 / denominator
#     df$ancestral_state_2[df$id == node] <<- numerator_2 / denominator
#   }
# }))
# 
# # Output log report
# log_df <- data.frame(
#   ancestor_node_id = ancestor_node_ids,
#   ancestral_state_1 = tail(df$ancestral_state_1, phy[["Nnode"]]),
#   ancestral_state_2 = tail(df$ancestral_state_2, phy[["Nnode"]]),
#   bisse_state = phy[["node.state"]],
#   match = NA
# )
# match <- log_df$ancestral_state_1 > log_df$ancestral_state_2
# log_df$match <- match != log_df$bisse_state
