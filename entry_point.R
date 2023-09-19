library("rjson")

source(file.path("ode_solvers.R"))
source(file.path("tree_generator.R"))

# TODO need to decide data input standards and implement validation
params_json_path <- file.path("params.json")
params <- fromJSON(file = params_json_path)
# Speciation rate
λ <- params$lambda
# Extinction rate
μ <- params$mu
# Sampling rate over time and lineages
Ψ <- params$psi
# Transition matrix
q <- do.call("cbind", params$q)

sampled_states_freqs_csv_path <- file.path("two_sampled_states_freqs.csv")
state_freqs_df <- read.csv(sampled_states_csv_path, header = FALSE)

#Transition matrix
sampled_states_q_csv_path <- file.path("two_sampled_states_q.csv")
q_matrix <- as.matrix(read.csv(sampled_states_q_csv_path))
rownames(q_matrix) <- q_matrix[,1]
q_matrix <- q_matrix[, -c(1)]

# TODO test input is a binomial tree
phy <- get_rand_phy(3, 100, λ, μ, q)
# TODO temporarily for comparison purposes
plot(history.from.sim.discrete(phy, 0:1), phy)
tiplabels(frame = "circle", cex = 0.5)
nodelabels(frame = "circle", cex = 0.5)

leaf_node_ids <- seq_along(phy[["tip.state"]])
leaf_nodes_df <- data.frame(
  id = leaf_node_ids,
  # Time to root
  t_root = NA,
  # Children ids
  left = NA,
  right = NA,
  # Parent id
  parent = NA,
  # State probabilities used in backwards time eqns
  backwards_state_1 = ifelse(phy[["tip.state"]] == 0, state_1_freq, 0),
  backwards_state_2 = ifelse(phy[["tip.state"]] == 1, state_2_freq, 0),
  # Ancestral state reconstructions
  ancestral_state_1 = NA,
  ancestral_state_2 = NA)

ancestor_node_ids <- seq(length(leaf_node_ids) + 1,
                         length(leaf_node_ids) + phy[["Nnode"]])
ancestor_nodes_df <- data.frame(
  id = ancestor_node_ids,
  t_root = NA,
  left = NA,
  right = NA,
  parent = NA,
  backwards_state_1 = NA,
  backwards_state_2 = NA,
  ancestral_state_1 = NA,
  ancestral_state_2 = NA)

df <- rbind(leaf_nodes_df, ancestor_nodes_df)

# Populate df with `t_root` vals
node_depths <- node.depth.edgelength(phy)
max_depth <- max(node_depths)
invisible(lapply(seq_len(length(node_depths)), function(i) {
  df$t_root[df$id == i] <<- max_depth - node_depths[i]
}))

# Populate df with connections
post_order_edges <- reorder.phylo(phy, "postorder")[["edge"]]
invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
  node <- post_order_edges[[i, 1]]
  left <- post_order_edges[[i, 2]]
  right <- post_order_edges[[i + 1, 2]]

  df$parent[df$id == left] <<- node
  df$parent[df$id == right] <<- node
  df$left[df$id == node] <<- left
  df$right[df$id == node] <<- right
}))

# Backward-time equations
invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
  node <- post_order_edges[[i, 1]]
  left <- post_order_edges[[i, 2]]
  right <- post_order_edges[[i + 1, 2]]

  tf <- df$t_root[df$id == node]
  left_t0 <- df$t_root[df$id == left]
  right_t0 <- df$t_root[df$id == right]

  left_backwards_state_1 <- df$backwards_state_1[df$id == left]
  left_backwards_state_2 <- df$backwards_state_2[df$id == left]
  left_backwards_sol <- get_backwards_sol(left_backwards_state_1,
                                          left_backwards_state_2,
                                          left_t0, tf,
                                          λ, μ, Ψ, q, node)

  right_backwards_state_1 <- df$backwards_state_1[df$id == right]
  right_backwards_state_2 <- df$backwards_state_2[df$id == right]
  right_backwards_sol <- get_backwards_sol(right_backwards_state_1,
                                           right_backwards_state_2,
                                           right_t0, tf,
                                           λ, μ, Ψ, q, node)

  lb1 <- left_backwards_sol[["backwards_state_1"]]
  rb1 <- right_backwards_sol[["backwards_state_1"]]
  backwards_state_1 <- λ * lb1 * rb1
  df$backwards_state_1[df$id == node] <<- backwards_state_1

  lb2 <- left_backwards_sol[["backwards_state_2"]]
  rb2 <- right_backwards_sol[["backwards_state_2"]]
  backwards_state_2 <- λ * lb2 * rb2
  df$backwards_state_2[df$id == node] <<- backwards_state_2
}))

# Forward-time equations
invisible(lapply(seq(length(post_order_edges[, 1]), 1, -2), function(i) {
  node <- post_order_edges[[i, 1]]
  left <- post_order_edges[[i - 1, 2]]
  right <- post_order_edges[[i, 2]]

  backwards_state_1 <- df$backwards_state_1[df$id == node]
  backwards_state_2 <- df$backwards_state_2[df$id == node]

  # If root else not root
  if (i == length(post_order_edges[, 1])) {
    denominator <- (backwards_state_1 * state_1_freq
                    + backwards_state_2 * state_2_freq)
    df$ancestral_state_1[df$id == node] <<-
      backwards_state_1 * state_1_freq / denominator
    df$ancestral_state_2[df$id == node] <<-
      backwards_state_2 * state_2_freq / denominator
  } else {
    parent <- df$parent[df$id == node]
    parent_ancestral_state_1 <- df$ancestral_state_1[df$id == parent]
    parent_ancestral_state_2 <- df$ancestral_state_2[df$id == parent]

    t0 <- df$t_root[df$id == parent]
    tf <- df$t_root[df$id == node]

    forwards_sol <- get_forwards_sol(parent_ancestral_state_1,
                                     parent_ancestral_state_2,
                                     t0, tf,
                                     λ, μ, Ψ, q, node)
    forwards_state_1 <- forwards_sol[["forwards_state_1"]]
    forwards_state_2 <- forwards_sol[["forwards_state_2"]]

    numerator_1 <- backwards_state_1 * forwards_state_1
    numerator_2 <- backwards_state_2 * forwards_state_2
    denominator <- numerator_1 + numerator_2
    df$ancestral_state_1[df$id == node] <<- numerator_1 / denominator
    df$ancestral_state_2[df$id == node] <<- numerator_2 / denominator
  }
}))

# Output log report
log_df <- data.frame(
  ancestor_node_id = ancestor_node_ids,
  ancestral_state_1 = tail(df$ancestral_state_1, phy[["Nnode"]]),
  ancestral_state_2 = tail(df$ancestral_state_2, phy[["Nnode"]]),
  bisse_state = phy[["node.state"]],
  match = NA
)
match <- log_df$ancestral_state_1 > log_df$ancestral_state_2
log_df$match <- match != log_df$bisse_state
