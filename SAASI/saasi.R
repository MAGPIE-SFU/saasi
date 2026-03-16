# This is probably going to solve all the problems of having so many codes

saasi <- function(pars,phy,q_matrix){
  k <- (sqrt(4+4*length(pars))-2)/2 
  
  # get the parameters
  param <- cbind(matrix(pars[seq_len(3*k)], k, 3),
                 matrix(pars[(3*k+1):(k*(k+2))], k, k-1, TRUE))
  
  tran <- matrix(pars[seq_len(3*k)], k, 3)
  params_df <- data.frame(state = c(1:k), freq = exp(tran[,3]/sum(tran[,3])),lambda = tran[1:k,1],mu = tran[1:k,2], psi = tran[1:k,3])
  nstate <- nrow(params_df)
  

  nleaf_node <- phy[["Nnode"]]
  nnode <- nleaf_node * 2 + 1
  topology_df <- data.frame(
    id = seq_len(nnode),
    t_root = NA,
    left = NA,
    right = NA,
    parent = NA)
  
  #Populate topology df with time distances from present day
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
  #print(backwards_likelihoods_list)
  backwards_likelihoods_list <- rep(list(rep(0, nstate)), nnode)
  
  # Populate leaf node likelihoods for backwards time equations
  invisible(lapply(seq_along(phy[["tip.state"]]), function(i) {
    state <- phy[["tip.state"]][[i]]
    state_freq <- params_df$freq[params_df$state == state]
    backwards_likelihoods_list[[i]][[state]] <<- log(state_freq)
  }))
  #print("backwards_likelihoods_list")
  #print(backwards_likelihoods_list)
  
  #backwards_likelihoods_list <- backwards_likelihoods_list/sum(backwards_likelihoods_list)
  
  #print("backwards_likelihoods_list after")
  #print(backwards_likelihoods_list)
  # Populate internal node likelihoods for backwards time equations
  invisible(lapply(seq(1, length(post_order_edges[, 1]), 2), function(i) {
    node <- post_order_edges[[i, 1]]
    left <- post_order_edges[[i, 2]]
    right <- post_order_edges[[i + 1, 2]]
    
    # print(node)
    
    tf <- topology_df$t_root[topology_df$id == node]
    left_t0 <- topology_df$t_root[topology_df$id == left]
    right_t0 <- topology_df$t_root[topology_df$id == right]
    #print(tf)
    #print(left_t0)
    #print(right_t0)
    left_likelihoods <- backwards_likelihoods_list[[left]]
    right_likelihoods <- backwards_likelihoods_list[[right]]
    #print("left")
    #print(left_likelihoods)
    #print("right")
    #print(right_likelihoods)
    
    likelihoods <- get_backwards_likelihoods(abs(left_likelihoods), abs(right_likelihoods),
                                             left_t0, right_t0, tf,
                                             params_df, q_matrix)
    
    likelihoods <- abs(likelihoods)
    #print("likelihoods before normalizing")
    #print(likelihoods)
    #print("backwards_likelihoods_list[[node]]")
    #print(backwards_likelihoods_list[[node]])
    #backwards_likelihoods_list[[node]] <<- likelihoods/sum(likelihoods)
    backwards_likelihoods_list[[node]] <<- likelihoods/sum(likelihoods)
    #print("likelihoods")
    #print(likelihoods)
    
    #print("backwards_likelihoods_list[[node]]")
    #print(backwards_likelihoods_list[[node]])
    
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
  root_node <- nleaf_node + 2
  # state_probabilities_list[[root_node]] <- (
  #   backwards_likelihoods_list[[root_node]] * log(rev(params_df$freq))
  #   / (sum(backwards_likelihoods_list[[root_node]] * log(rev(params_df$freq)))
  # ))
  state_probabilities_list[[root_node]] <- (
    backwards_likelihoods_list[[root_node]] * log(rev(params_df$freq))
    / (sum(backwards_likelihoods_list[[root_node]] * log(rev(params_df$freq)))
    ))
  #print(state_probabilities_list[[root_node]])
  #print(backwards_likelihoods_list)
  # Populate internal node state probabilities
  invisible(lapply(seq(length(post_order_edges[, 1]), 1, -2), function(i) {
    node <- post_order_edges[[i, 1]]
    if (node == root_node) {
      return()
    }
    #print(node)
    parent <- topology_df$parent[topology_df$id == node]
    #print(parent)

    if(node == root_node){
      parent_state_probabilities <- state_probabilities_list[[parent]]
    }else{
      parent_state_probabilities <- state_probabilities_list[[parent]]*backwards_likelihoods_list[[node]]
    }
    
    t0 <- topology_df$t_root[topology_df$id == parent]
    tf <- topology_df$t_root[topology_df$id == node]
    parent_state_probabilities <- abs(parent_state_probabilities)/sum(abs(parent_state_probabilities))
    
    likelihoods <- get_forwards_likelihoods(parent_state_probabilities,
                                            t0, tf,
                                            params_df, q_matrix)
    likelihoods <- abs(likelihoods)
    state_probabilities_list[[node]] <<- (
      backwards_likelihoods_list[[node]] * likelihoods
      / sum(backwards_likelihoods_list[[node]] * likelihoods)
    )
  }))
  
  internal_node_ids <- root_node:nnode
  
  node_prob <- matrix(unlist(state_probabilities_list[internal_node_ids]),ncol=k,byrow = TRUE)
  node_lik <- as.data.frame(node_prob,row.names = internal_node_ids)
  colnames(node_lik) = c(1:k)
  return(node_lik)
}