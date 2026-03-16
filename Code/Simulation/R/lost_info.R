setwd("/Users/yexuan-magpie/Desktop/pastml_treetime/")

get_trans_hist <- function(true_tree_info,meta_data){
  n_tips <- length(true_tree_info$tip.label)
  n_nodes <- true_tree_info$Nnode
  n_total <- n_tips + n_nodes
  
  states <- rep(NA, n_total)
  
  for (i in 1:n_tips) {
    tip_name <- true_tree_info$tip.label[i]
    states[i] <- meta_data$State[meta_data$label == tip_name]
  }
  
  for (i in 1:n_nodes) {
    node_num <- n_tips + i
    if (!is.null(true_tree_info$node.label) && length(true_tree_info$node.label) >= i) {
      node_name <- true_tree_info$node.label[i]
      states[node_num] <- meta_data$State[meta_data$label == node_name]
    }
  }
  
  transitions <- data.frame(
    parent_node = true_tree_info$edge[, 1],
    child_node = true_tree_info$edge[, 2],
    parent_state = states[true_tree_info$edge[, 1]],
    child_state = states[true_tree_info$edge[, 2]]
  )
  
  
  transition_counts <- list(
    state1_to_state2 = sum(transitions$parent_state == 1 & 
                             transitions$child_state == 2, na.rm = TRUE),
    state2_to_state1 = sum(transitions$parent_state == 2 & 
                             transitions$child_state == 1, na.rm = TRUE),
    state1_to_state1 = sum(transitions$parent_state == 1 & 
                             transitions$child_state == 1, na.rm = TRUE),
    state2_to_state2 = sum(transitions$parent_state == 2 & 
                             transitions$child_state == 2, na.rm = TRUE)
  )
  return(transition_counts$state1_to_state2+transition_counts$state2_to_state1+
           transition_counts$state1_to_state1 + transition_counts$state2_to_state2)
  
}

get_diff_hist <- function(tree1,meta_data1,tree2,meta_data2){
  return(abs(get_trans_hist(tree1,meta_data1) - get_trans_hist(tree2,meta_data2)))
}


# testing

tree1 <- read.tree("result_baseline_small_tree_subsample/sim_1/tree_1.nwk")
meta_data1 <- read.csv("result_baseline_small_tree_subsample/sim_1/true_tree_1.csv")

get_trans_hist(tree1,meta_data1)

tree2 <- read.tree("result_baseline_small_tree_subsample/sim_1/tree_before_ds_1.nwk")
meta_data2 <- read.csv("result_baseline_small_tree_subsample/sim_1/true_tree_before_ds_1.csv")
get_trans_hist(tree2,meta_data2)

get_diff_hist(tree1,meta_data1,tree2,meta_data2)
